from __future__ import absolute_import, division, print_function

# More Python 3 vs 2 mojo
try:
    import queue
except ImportError:
    import Queue as queue

#import pygraphviz as pgv
import targetc as genc

from constructs import *
from expression import *
from codegen import *
from schedule import *
from poly import *
from bounds import *
from inline import *

# LOG CONFIG #
pipe_logger = logging.getLogger("pipe.py")
pipe_logger.setLevel(logging.INFO)
LOG = pipe_logger.log

def get_parents_from_comp_obj(comp):
    refs = comp.getObjects(Reference)
    # Filter out self and image references 
    refs = [ ref for ref in refs if not ref.objectRef == comp and \
                                 not isinstance(ref.objectRef, Image) ]
    return list(set([ref.objectRef for ref in refs]))

def get_comp_objs_and_dep_maps(outputs):
    """
    Find all the compute objects required for the outputs and
    also builds parent and children maps for the compute objects
    """
    comp_objs = []
    comp_objs_parents = {}
    comp_objs_children = {}
    q = queue.Queue()
    for comp_obj in outputs:
        q.put(comp_obj)
    while not q.empty():
        obj = q.get()
        parent_objs = get_parents_from_comp_obj(obj)
        if obj not in comp_objs:
            comp_objs.append(obj)
            comp_objs_parents[obj] = parent_objs
            for parobj in parent_objs:
                if parobj in comp_objs_children:
                    if obj not in comp_objs_children[parobj]:
                        comp_objs_children[parobj].append(obj)
                else:
                    comp_objs_children[parobj] = [obj]
            if len(parent_objs) != 0:
                for r in parent_objs:
                    q.put(r)

    for comp in comp_objs:
        if comp not in comp_objs_parents:
            comp_objs_parents[comp] = []
        if comp not in comp_objs_children:
            comp_objs_children[comp] = []

    return comp_objs, comp_objs_parents, comp_objs_children

class Group:
    """ 
        Group is a part of the pipeline which realizes a set of computation
        objects. Scheduling and storage allocation is done at the level of a 
        group. A group also maintains a polyhedral representation of the 
        computation objects when possible.
    """
    # Construct a group from a set of language functions / reductions
    def __init__(self, _ctx, _comp_objs, _param_constraints):
        # All the computation constructs in the language derive from the
        # Function class. Input images cannot be part of a group.
        for comp in _comp_objs:
            assert(isinstance(comp, Function))
            #assert(not isinstance(comp, Image))

        self._comp_objs  = _comp_objs
        self._polyrep = None
        refs = []

        for comp in self._comp_objs:
            refs += comp.getObjects(Reference)

        self._inputs = list(set([ref.objectRef for ref in refs \
                            if isinstance(ref.objectRef, Image)]))

        # Create a polyhedral representation if possible.
        # Currently doing extraction only when all the compute_objs
        # domains are affine. This can be revisited later.
        if self.isPolyhedral():
            self._polyrep = PolyRep(_ctx, self, [], _param_constraints)

    @property
    def compute_objs(self):
        return self._comp_objs
    @property
    def polyRep(self):
        return self._polyrep
    @property
    def inputs(self):
        return self._inputs

    def is_fused(self):
        return len(self.compute_objs) > 1

    def getParameters(self):
        params = []
        for comp in self._comp_objs:
            params = params + comp.getObjects(Parameter)
        return list(set(params))

    def isPolyhedral(self):
        polyhedral = True
        for comp in self._comp_objs:
            if (not comp.hasBoundedIntegerDomain()):
                polyhedral = False
        return polyhedral

    def orderComputeObjs(self):
        # Order stores the numbering of each compute object 
        # when topologically sorted.
        order = {}
        # Initialize all the initial numbering to zero for
        # all compute objects in the group
        for comp in self._comp_objs:
            order[comp] = 0
        # Doing a topological sort in an iterative fashion
        change = True
        while(change):
            change = False
            for comp in self._comp_objs:
                parentObjs = get_parents_from_comp_obj(comp)
                for pobj in parentObjs:
                    if (pobj in order  and (order[pobj] >= order[comp])):
                        order[comp] = order[pobj] + 1
                        change = True
        return order

    def __str__(self):
        comp_str  = "\n\n".join([comp.__str__() \
                    for comp in self._comp_objs]) + '\n'
        return comp_str + '\n' + self._polyrep.__str__()

class Pipeline:
    def __init__(self, _ctx, _outputs,
                 _param_estimates, _param_constraints,
                 _grouping, _inline_directives,
                 _options, _name = None):
        # Name of the pipleline is a concatenation of the names of the 
        # pipeline outputs, unless it is explicitly named.
        if _name is None:
            _name = ''.join([out.name for out in _outputs])

        self._name = _name

        self._ctx = _ctx
        self._org_outputs = _outputs
        self._param_estimates = _param_estimates
        self._param_constraints = _param_constraints
        self._grouping = _grouping
        self._inline_directives = _inline_directives
        self._options = _options

        ''' CONSTRUCT DAG '''
        # Maps from a compute object to its parents and children by
        # backtracing starting from given live-out functions.
        # TODO: see if there is a cyclic dependency between the compute
        # objects. Self references are not treated as cycles.
        self._org_comp_objs, \
         _, \
          _ = \
            get_comp_objs_and_dep_maps(self._org_outputs)

        ''' CLONING '''
        # Clone the computation objects i.e. functions and reductions
        self._clone_map = {}
        for comp in self._org_comp_objs:
            self._clone_map[comp] = comp.clone()
        self._outputs = [self._clone_map[obj] for obj in self._org_outputs]
        # Modify the references in the cloned objects (which refer to
        # the original objects)
        for comp in self._org_comp_objs:
            cln = self._clone_map[comp]
            refs = cln.getObjects(Reference)
            for ref in refs:
                if not isinstance(ref.objectRef, Image):
                    ref._replace_ref_object(self._clone_map[ref.objectRef])

        ''' DAG OF CLONES '''
        self._comp_objs, \
         self._comp_objs_parents, \
          self._comp_objs_children = \
            get_comp_objs_and_dep_maps(self._outputs)

        ''' INITIAL GROUPING '''
        # Create a group for each pipeline function / reduction and compute
        # maps for parent and child relations between the groups
        self._groups, \
         self._group_parents, \
          self._group_children = \
            self.build_initial_groups()

        # Store the initial pipeline graph. The compiler can modify 
        # the pipeline by inlining functions.
        #self._initialGraph = self.drawPipelineGraph()

        # Make a list of all the input groups
        inputs = []
        for f in self._groups:
            inputs = inputs + self._groups[f].inputs
        self._inputs = list(set(inputs))

        # Checking bounds
        self.bounds_check_pass()

        # inline pass
        self.inline_pass()

        # make sure the set of functions to be inlined and those to be grouped
        # are disjoint
        if self._inline_directives and self._grouping:
            group_comps = []
            for g in self._grouping:
                group_comps += g
            a = set(self._inline_directives)
            b = set(group_comps)
            assert a.isdisjoint(b)

        ''' GROUPING '''
        # TODO check grouping validity
        if self._grouping:
            # for each group
            for g in self._grouping:
                # get clones of all functions
                merge_group_list = \
                    [self._groups[self._clone_map[f]] for f in g]
                if len(merge_group_list) > 1:
                     merged = merge_group_list[0]
                     for i in range(1, len(merge_group_list)):
                        merged = \
                            self.merge_groups(merged, merge_group_list[i])
                        # to be done after each merge, to know if the merging
                        # was valid.
                        align_and_scale_parts(self, merged)
        else:
            # Run the grouping algorithm
            pass

        # ***
        log_level = logging.INFO
        LOG(log_level, "Grouped compute objects:")
        glist = list(set([g for g in self._groups.values()]))
        for g in glist:
            LOG(log_level, [comp.name for comp in g.compute_objs])
        # ***

        ''' BASE SCHEDULE AND CODEGEN '''
        for g in list(set(self._groups.values())):
            gparts = base_schedule(g)
            for p in gparts:
                p.liveout = True
            #fused_schedule(g)  # (g, param_estimates)

    @property
    def groups(self):
        return self._groups

    @property
    def name(self):
        return self._name

    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs

    @property
    def originalGraph(self):
        return self._initialGraph

    def getParameters(self):
        params=[]
        for group in self._groups.values():
            params = params + group.getParameters()
        return list(set(params))

    def build_initial_groups(self):
        """
        Place each compute object of the pipeline in its own Group, and set the
        dependence relations between the created Group objects.
        """
        comp_objs = self._comp_objs
        comp_parents = self._comp_objs_parents
        comp_children = self._comp_objs_children
        group_map = {}
        group_parents = {}
        group_children = {}
        for comp in comp_objs:
            group_map[comp] = Group(self._ctx, [comp], self._param_constraints)

        for comp in group_map:
            group_parents[group_map[comp]] = \
                [ group_map[p] for p in comp_parents[comp] ]
            group_children[group_map[comp]] = \
                [ group_map[c] for c in comp_children[comp] ]

        return group_map, group_parents, group_children

    '''
    def drawPipelineGraph(self):
        G = pgv.AGraph(strict=False, directed=True)
        groupList = list(set([self._groups[f] for f in self._groups]))

        # TODO add input nodes to the graph
        for i in range(0, len(groupList)):
            subGraphNodes = []
            for obj in self._comp_objs:
                if self._groups[obj] == groupList[i]:
                    subGraphNodes.append(obj.name)
            G.add_nodes_from(subGraphNodes)
            G.add_subgraph(nbunch = subGraphNodes, 
                           name = "cluster_" + str(i))

        for obj in self._comp_objs:
            for pobj in self._comp_objs_parents[obj]:
                G.add_edge(pobj.name, obj.name)

        G.layout(prog='dot')
        return G
    '''

    def generate_code(self, outputs_no_alloc=False,
                            is_extern_c_func=False,
                            are_io_void_ptrs=False):

        """
        Code generation for the entire pipeline starts here.

        Flags:

        1. "outputs_no_alloc"
        (*) True => memory allocation code is not generated
        (*) False => malloc / pool_malloc / multidim array


        2. "is_extern_c_func"
        (*) True => function declaration generated with ' extern "C" ' string
                    (used when dynamic libs are needed for python wrapping)
        (*) False => normal C function declaration


        3. "are_io_void_ptrs"
        (*) True => all inputs and outputs of the pipeline are expected, by the
                    C function declaration, to be passed as 'void *'
                    (used when dynamic libs are needed for python wrapping)
        (*) False => inputs and outputs are to be passed as pointers of their
                     data type. E.g: 'float *'
        """

        return generate_code_for_pipeline(self, outputs_no_alloc,
                                                is_extern_c_func,
                                                are_io_void_ptrs)

    def drop_compute_obj(self, comp_obj):
        # if the compute object is a child of any other
        if comp_obj in self._comp_objs_parents:
            for p_comp in self._comp_objs_parents[comp_obj]:
                self._comp_objs_children[p_comp].remove(comp_obj)
            self._comp_objs_parents.pop(comp_obj)
        # if the compute object is a parent of any other
        if comp_obj in self._comp_objs_children:
            for c_comp in self._comp_objs_children[comp_obj]:
                self._comp_objs_parents[c_comp].remove(comp_obj)
            self._comp_objs_children.pop(comp_obj)
        # remove comp_obj
        self._comp_objs.remove(comp_obj)

        return

    def add_group(self, group):
        """
        add a new group to the pipeline
        """
        if not group._comp_objs:
            return

        parents = []
        children = []
        for comp in group.compute_objs:
            # add group's parents
            parents.extend( \
                [self._groups[p] for p in self._comp_objs_parents[comp] \
                                   if p not in group.compute_objs] )
            # add group's children
            children.extend( \
                [self._groups[c] for c in self._comp_objs_children[comp] \
                                   if c not in group.compute_objs] )
            # add group to comp -> group mapping
            self._groups[comp] = group
        self._group_parents[group] = list(set(parents))
        self._group_children[group] = list(set(children))

        # add group as child for all its parents
        for parent_group in self._group_parents[group]:
            self._group_children[parent_group].append(group)
        # add group as parent for all its children
        for child_group in self._group_children[group]:
            self._group_parents[child_group].append(group)

        return

    def drop_group(self, group):
        """
        drop the group from the pipeline
        """
        # if group is a child of any other group
        if group in self._group_parents:
            for p_group in self._group_parents[group]:
                self._group_children[p_group].remove(group)
            self._group_parents.pop(group)
        # if group is a parent of any other group
        if group in self._group_children:
            for c_group in self._group_children[group]:
                self._group_parents[c_group].remove(group)
            self._group_children.pop(group)
        # remove group from comp -> group mapping
        for comp in group.compute_objs:
            self._groups.pop(comp)

        return

    def merge_groups(self, g1, g2):
        # Get comp objects from both groups
        comp_objs = g1.compute_objs + g2.compute_objs
        comp_objs = list(set(comp_objs))
        # Create a new group
        merged = Group(self._ctx, comp_objs, self._param_constraints)

        self.drop_group(g1)
        self.drop_group(g2)
        self.add_group(merged)

        return merged

    def get_ordered_groups(self):
        # Assign level numbers to each group and sort accourding to the level
        group_order = {}
        groups = set(self._groups.values())
        group_list = [ [g, len(self._group_parents[g])] for g in groups ]

        level = 0
        while group_list:
            # find all the groups whose parents have their levels assigned
            level_assigned = [ t for t in group_list if t[1] == 0 ]
            for assgn in level_assigned:
                group_order[assgn[0]] = level
                group_list.remove(assgn)
                # reduce the unassigned parent count for all the children
                child_groups = self._group_children[assgn[0]]
                for assgn in group_list:
                    if assgn[0] in child_groups:
                        assgn[1] -= 1
            level = level + 1
        
        return sorted(group_order.items(), key=lambda x: x[1])

    def bounds_check_pass(self):
        """ 
        Bounds check pass analyzes if function values used in the compute
        objects are within the domain of the functions. Static analysis is
        only possible when the references to function values are regular
        i.e. they are not data dependent. We restrict ourselves to affine
        references.
        """
        for group in self._groups.values():
            for child in self._group_children[group]:
                check_refs(child, group)
            for inp in group.inputs:
                # Creating a computation group for an input which is given
                # is meaningless. Ideally it should be done in a clean way
                # currently abusing group for construction of a polyhedral
                # representation
                inp_group = Group(self._ctx, [inp], self._param_constraints)

                check_refs(group, inp_group)

        return

    def inline_pass(self):
        """ 
        Inline pass takes all the inlining decisions and inlines functions at
        their use points in other groups.
        """
        # Guidelines for inlining
        # -- Inlining functions which have multiple case constructs can quickly
        #    turn out to be messy code generation nightmare.
        # -- Splitting a use site can occur when the inlined function is defined
        #    by different expression in different parts of the function domain.
        #    Functions which are defined by a single expression over the entire
        #    domain are good candidates for inlining.
        #    Example :
        #    f(x) = g(x-1) + g(x) + g(x+1) when (1 <= x <= N-1)
        #         = 0                      otherwise
        #    h(x) = f(x-1) + f(x) + f(x+1) when (1 <= x <= N-1)
        #         = 0                      otherwise
        #    Inlining f into h will result in splitting as shown below
        #    h(x) = g(x-2) + g(x-1) + g(x) +      when (2 <= x <= N-2)
        #           g(x-1) + g(x) + g(x+1) +
        #           g(x)   + g(x+1) + g(x+2)
        #         = 2*g(x-1) + 3*g(x) + 2*g(x+1)  when (x == 1)
        #           + g(3)
        #         = g(x-2) + 2*g(x-1) +           when (x == N-1)
        #           3*g(x) + 2*g(x+1)
        #         = 0                             otherwise
        #    For multiple dimensions and complex expression it gets ugly fast
        # -- Inlining without splitting is possible even when the function is
        #    defined by mutliple expressions. When the function references at
        #    the use site can all be proved to be generated from a single
        #    expression.
        #    Example :
        #    f(x) = g(x-1) + g(x) + g(x+1) when (1 <= x <= N-1)
        #         = 0                      otherwise
        #    h(x) = f(x-1) + f(x) + f(x+1) when (2 <= x <= N-2)
        #         = 0                      otherwise
        #    Inlining f into h will not result in splitting as shown below
        #    h(x) = g(x-2) + g(x-1) + g(x) +      when (2 <= x <= N-2)
        #           g(x-1) + g(x) + g(x+1) +
        #           g(x)   + g(x+1) + g(x+2)
        #         = 0                             otherwise
        # -- Users should be encouraged to write functions in a form that makes
        #    inlining easy.
        inlined_comp_objs = []
        groups = list(set(self._groups.values()))
        # TODO: _inline_directives flag
        for directive in self._inline_directives:
            # Only function constructs can be inlined for now
            assert isinstance(directive, Function)

            comp = self._clone_map[directive]
            # Does inling into a fused group cause problems?
            group = self._groups[comp]

            # One simply does not walk into Inlining
            assert group._comp_objs[0] not in self._outputs

            for child in self._group_children[group]:
                ref_to_inline_expr_map = \
                    inline_piecewise(child, group, no_split=True)
                if ref_to_inline_expr_map:
                    child._comp_objs[0].replace_refs(ref_to_inline_expr_map)
            # Recompute group graph
            # TODO: tweak the graph edges instead of rebuilding the whole graph
            #self._groups = self.buildGroupGraph()
            self.drop_group(self._groups[comp])
            self.drop_compute_obj(comp)

    # TODO printing the pipeline
    def __str__(self):
        return_str = "Final Group: " + self._name + "\n"
        for s in self._groups:
            return_str = return_str + s.__str__() + "\n"
        return return_str

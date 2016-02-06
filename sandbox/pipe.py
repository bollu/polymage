from __future__ import absolute_import, division, print_function

# More Python 3 vs 2 mojo
try:
    import queue
except ImportError:
    import Queue as queue

import pygraphviz as pgv
import targetc as genc

from grouping import *

from constructs import *
from expression import *
from codegen import *
from schedule import *
from poly_schedule import *
from align_scale import *
from poly import *
from bounds import *
from inline import *
from storage_mapping import *

# LOG CONFIG #
pipe_logger = logging.getLogger("pipe.py")
pipe_logger.setLevel(logging.INFO)
LOG = pipe_logger.log

def get_parents_from_func(comp):
    refs = comp.getObjects(Reference)
    # Filter out self and image references 
    refs = [ ref for ref in refs if not ref.objectRef == comp and \
                                 not isinstance(ref.objectRef, Image) ]
    return list(set([ref.objectRef for ref in refs]))

def get_funcs_and_dep_maps(outputs):
    """
    Find all the compute objects required for the outputs and
    also builds parent and children maps for the compute objects
    """
    funcs = []
    funcs_parents = {}
    funcs_children = {}
    # queue of compute objects
    q = queue.Queue()
    for func in outputs:
        q.put(func)
    while not q.empty():
        obj = q.get()
        parent_objs = get_parents_from_func(obj)
        if obj not in funcs:
            funcs.append(obj)
            funcs_parents[obj] = parent_objs
            for parobj in parent_objs:
                if parobj in funcs_children:
                    if obj not in funcs_children[parobj]:
                        funcs_children[parobj].append(obj)
                else:
                    funcs_children[parobj] = [obj]
            if len(parent_objs) != 0:
                for r in parent_objs:
                    q.put(r)

    for comp in funcs:
        if comp not in funcs_parents:
            funcs_parents[comp] = []
        if comp not in funcs_children:
            funcs_children[comp] = []

    return funcs, funcs_parents, funcs_children

def get_funcs(outputs):
    """
    Find all the compute objects required for the outputs and
    also builds parent and children maps for the compute objects
    """
    funcs = []
    # queue of compute objects
    q = queue.Queue()
    for func in outputs:
        q.put(func)
    while not q.empty():
        obj = q.get()
        parent_objs = get_parents_from_func(obj)
        if obj not in funcs:
            funcs.append(obj)
            if len(parent_objs) != 0:
                for r in parent_objs:
                    q.put(r)

    return funcs

def level_order(objs, parent_map):
    # Order stores the numbering of each object when topologically sorted.
    order = {}
    # Initialize all the initial numbering to zero for all objects
    for obj in objs:
        order[obj] = 0
    # Doing a topological sort in an iterative fashion
    change = True
    while(change):
        change = False
        for obj in objs:
            parent_objs = parent_map[obj]
            if parent_objs is None:
                continue
            for p_obj in parent_objs:
                if (p_obj in order and (order[p_obj] >= order[obj])):
                    order[obj] = order[p_obj] + 1
                    change = True
    return order

class ComputeObject:
    def __init__(self, _func):
        self._func = _func
        self._parents = []
        self._children = []
        self._group = None
        self._is_liveout = True

        self.set_flags()

    @property
    def func(self):
        return self._func

    @property
    def is_parents_set(self):
        return self._is_parents_set
    @property
    def is_children_set(self):
        return self._is_children_set
    @property
    def is_group_set(self):
        return self._is_group_set

    @property
    def parents(self):
        assert self.is_parents_set
        return self._parents
    @property
    def children(self):
        assert self.is_children_set
        return self._children
    @property
    def group(self):
        assert self.is_group_set
        return self._group

    @property
    def is_liveout(self):
        return self._is_liveout

    def set_flags(self):
        self._is_parents_set = False
        self._is_children_set = False
        self._is_group_set = False

    def set_parents(self, parents):
        assert self._parents_set == False
        # assert parents is a list of ComputeObject
        for p in parents:
            assert isinstance(p, ComputeObject)
        self._parents = parents

    def set_children(self, children):
        assert self._children_set == False
        # assert children is a list of ComputeObject
        for p in children:
            assert isinstance(p, ComputeObject)
        self._children = children

    def set_group(self, group):
        assert isinstance(group, Group)
        self._group = group

class Group:
    """ 
        Group is a part of the pipeline which realizes a set of computation
        objects. Scheduling and storage allocation is done at the level of a 
        group. A group also maintains a polyhedral representation of the 
        computation objects when possible.
    """
    # Construct a group from a set of language functions / reductions
    def __init__(self, _ctx, _comp_objs, \
                 _param_constraints):
        # All the computation constructs in the language derive from the
        # Function class. Input images cannot be part of a group.
        for comp in _comp_objs:
            assert(isinstance(comp, Function))
            #assert(not isinstance(comp, Image))

        self._comp_objs  = _comp_objs
        self._level_order_comps = self.order_compute_objs()
        self._root_comps = [comp for comp in self._comp_objs \
                                   if self._level_order_comps[comp] == 0]
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
    def comps(self):
        return self._comp_objs
    @property
    def polyRep(self):
        return self._polyrep
    @property
    def inputs(self):
        return self._inputs
    @property
    def name(self):
        return [comp.func.name for comp in self.comps]

    @property
    def get_ordered_compobjs(self):
        return self._level_order_comps

    # DEAD?
    def is_fused(self):
        return len(self.comps) > 1

    def getParameters(self):
        params = []
        for comp in self.comps:
            params = params + comp.func.getObjects(Parameter)
        return list(set(params))

    def isPolyhedral(self):
        polyhedral = True
        for comp in self.comps:
            if (not comp.func.hasBoundedIntegerDomain()):
                polyhedral = False
        return polyhedral

    def order_compute_objs(self):
        parents = {}
        for comp in self.comps:
            parents[comp] = comp.parents
        order = level_order(self.comps, parents)
        return order

    def get_sorted_compobjs(self):
        sorted_comps = sorted(self._level_order_comps.items(),
                              key=lambda x: x[1])
        sorted_comps = [c[0] for c in sorted_comps]
        return sorted_comps

    def __str__(self):
        comp_str  = '[' + \
                    ', '.join([comp.name \
                        for comp in self._comp_objs]) + \
                    ']'
        return comp_str

class Pipeline:
    def __init__(self, _ctx, _outputs,
                 _param_estimates, _param_constraints,
                 _grouping, _group_size, _inline_directives,
                 _tile_sizes, _size_threshold,
                 _options, _name = None):
        # Name of the pipleline is a concatenation of the names of the 
        # pipeline outputs, unless it is explicitly named.
        if _name is None:
            _name = ''.join([out.name for out in _outputs])

        self._name = _name

        self._ctx = _ctx
        self._orig_outputs = _outputs
        self._param_estimates = _param_estimates
        self._param_constraints = _param_constraints
        self._grouping = _grouping
        self._group_size = _group_size
        self._inline_directives = _inline_directives
        self._options = _options
        self._size_threshold = _size_threshold
        self._tile_sizes = _tile_sizes

        ''' CONSTRUCT DAG '''
        # Maps from a compute object to its parents and children by
        # backtracing starting from given live-out functions.
        # TODO: see if there is a cyclic dependency between the compute
        # objects. Self references are not treated as cycles.
        self._orig_funcs = get_funcs(self._orig_outputs)

        ''' CLONING '''
        # Clone the computation objects i.e. functions and reductions
        self._clone_map = {}
        for func in self._orig_funcs:
            self._clone_map[func] = func.clone()
        self._outputs = [self._clone_map[obj] for obj in self._orig_outputs]
        # Modify the references in the cloned objects (which refer to
        # the original objects)
        for func in self._orig_funcs:
            cln = self._clone_map[func]
            refs = cln.getObjects(Reference)
            for ref in refs:
                if not isinstance(ref.objectRef, Image):
                    ref._replace_ref_object(self._clone_map[ref.objectRef])

        ''' DAG OF CLONES '''
        _funcs, _func_parents, _func_children = \
            get_funcs_and_dep_maps(self._outputs)

        self._comps = \
            self.create_compute_objects(_funcs, _funcs_parents, _func_children)

        self._level_order_comps = self.order_compute_objs()

        ''' INITIAL GROUPING '''
        # Create a group for each pipeline function / reduction and compute
        # maps for parent and child relations between the groups
        self._groups = self.build_initial_groups()

        # Store the initial pipeline graph. The compiler can modify 
        # the pipeline by inlining functions.
        self._initial_graph = self.draw_pipeline_graph()

        # Make a list of all the input groups
        inputs = []
        for f in self._groups:
            inputs = inputs + self._groups[f].inputs
        self._inputs = list(set(inputs))

        # Checking bounds
        bounds_check_pass(self)

        # inline pass
        #inline_pass(self)

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
        if self._grouping and False:
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
        else:
            # Run the grouping algorithm
            auto_group(self)
            pass

        # create list of Groups
        self._group_list = list(set(self._groups.values()))

        # level order traversal of groups
        self._level_order_groups = self.order_group_objs()

        # ***
        log_level = logging.INFO
        LOG(log_level, "Grouped compute objects:")
        glist = list(set([g for g in self._groups.values()]))
        for g in glist:
            LOG(log_level, [comp.name for comp in g.comps])
        # ***

        ''' SCHEDULING '''
        for g in self._group_list:
            # alignment and scaling
            align_and_scale(self, g)
            # base schedule
            gparts = base_schedule(g)
            for p in gparts:
                p.compute_liveness(self)
            # grouping and tiling
            fused_schedule(self, g, self._param_estimates)

        self._grp_schedule = schedule_groups(self)

        # use graphviz to create pipeline graph
        self._pipeline_graph = self.draw_pipeline_graph()

        ''' STORAGE OPTIMIZATION '''

        # storage opt for liveout (full array) allocations
        # 1. derive liveout comps schedule from group schedule
        # 2. classify the storage based on type, dimensionality and size
        self.create_storage_classes()

    @property
    def comps(self):
        return self._comps
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
    def original_graph(self):
        return self._initial_graph
    @property
    def pipeline_graph(self):
        return self._pipeline_graph
    @property
    def get_ordered_compobjs(self):
        return self._level_order_comps
    @property
    def get_ordered_groups(self):
        return self._level_order_groups
    @property
    def group_schedule(self):
        return self._grp_schedule

    def get_parameters(self):
        params=[]
        for group in self._groups.values():
            params = params + group.getParameters()
        return list(set(params))

    def create_compute_objects(self, funcs, parents, children):
        comps = []
        func_map = {}
        for func in funcs:
            comp = ComputeObject(func, parents[func], children[func])
            func_map[func] = comp
            comps.append(comp)

        for func in func_map:
            comp = func_map[func]
            # set parents
            comp_parents = []
            for p_func in parents[func]:
                comp_parents.append(func_map[p_func])
            comp.set_parents(comp_parents)
            # set children
            comp_children = []
            for c_func in children[func]:
                comp_children.append(func_map[c_func])
            comp.set_children(comp_children)

        return comps

    def order_compute_objs(self):
        parents = {}
        for comp in comps:
            parents[comp] = comp.parents
        order = level_order(self._comps, parents)
        return order

    def order_group_objs(self):
        order = level_order(self._group_list, self._group_parents)
        return order

    def get_sorted_compobjs(self):
        sorted_comps = sorted(self._level_order_comps.items(),
                              key=lambda x: x[1])
        sorted_comps = [c[0] for c in sorted_comps]
        return sorted_comps

    def get_sorted_groups(self):
        sorted_groups = sorted(self._level_order_groups.items(),
                               key=lambda x: x[1])
        sorted_groups = [g[0] for g in sorted_groups]
        return sorted_groups

    def build_initial_groups(self):
        """
        Place each compute object of the pipeline in its own Group, and set the
        dependence relations between the created Group objects.
        """
        comps = self.comps
        group_map = {}
        for comp in comps:
            group = Group(self._ctx, [comp], self._param_constraints)
            group_map[comp] = group
            comp.set_group(group)

        for comp in comps:
            group = group_map[comp]
            # set group parents
            g_parents = [ group_map[p_comp] for p_comp in comp.parents ]
            group.set_parents(g_parents)
            # set group children
            g_children = [ group_map[p_comp] for p_comp in comp.children ]
            group.set_children(g_children)

        return group_map

    def draw_pipeline_graph(self):
        gr = pgv.AGraph(strict=False, directed=True)
        group_list = list(set([self.groups[comp] for comp in self.groups]))

        # TODO add input nodes to the graph
        for i in range(0, len(group_list)):
            sub_graph_nodes = []
            for comp in self.comps:
                if self.groups[comp] == group_list[i]:
                    sub_graph_nodes.append(comp.func.name)
            gr.add_nodes_from(sub_graph_nodes)
            gr.add_subgraph(nbunch = sub_graph_nodes,
                           name = "cluster_" + str(i))

        for comp in self.comps:
            for p_comp in comp.parents:
                gr.add_edge(p_comp.func.name, comp.func.name)

        gr.layout(prog='dot')
        return gr

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

    '''
    Pipelne graph operations
    '''

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
        for comp in group.comps:
            # add group's parents
            parents.extend( \
                [self._groups[p] for p in self._comp_objs_parents[comp] \
                                   if p not in group.comps] )
            # add group's children
            children.extend( \
                [self._groups[c] for c in self._comp_objs_children[comp] \
                                   if c not in group.comps] )
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
        for comp in group.comps:
            self._groups.pop(comp)

        return

    def merge_groups(self, g1, g2):
        # Get comp objects from both groups
        comp_objs = g1.comps + g2.comps
        comp_objs = list(set(comp_objs))
        # Create a new group
        merged = Group(self._ctx, comp_objs,
                       self._param_constraints)

        self.drop_group(g1)
        self.drop_group(g2)
        self.add_group(merged)

        return merged

    def replace_group(self, old_group, new_group):
        # if old_group has any child
        if old_group in self._group_children:
            for child in self._group_children[old_group]:
                self._group_parents[child].append(new_group)
                self._group_children[new_group].append(child)
        # if old_group has any parent
        if old_group in self._group_parents:
            for parent in self._group_parents[old_group]:
                self._group_children[parent].append(new_group)
                self._group_parents[new_group].append(parent)

        # replace old_group with new_group in groups list
        comp = old_group.comps[0]
        self.drop_group(old_group)
        self._groups[comp] = new_group

        return

    def make_comp_independent(self, comp_a, comp_b):
        """
        makes comp_b independent of comp_b and updates parent children
        relations in the graph structure
        [ assumes that comp_b is a child of comp_a, and that comp_a is inlined
        into comp_b ]
        """
        group_a = self._groups[comp_a]
        group_b = self._groups[comp_b]
        # if parent_comp has any parent
        if comp_a in self._comp_objs_parents:
            parents_of_a = self._comp_objs_parents[comp_a]
            parents_of_b = self._comp_objs_parents[comp_b]
            parents_of_grp_a = self._group_parents[group_a]
            parents_of_grp_b = self._group_parents[group_b]

            # remove relation between a and b
            self._comp_objs_children[comp_a].remove(comp_b)
            self._group_children[group_a].remove(group_b)
            parents_of_b.remove(comp_a)
            parents_of_grp_b.remove(group_a)

            # new parents list for b
            parents_of_b.extend(parents_of_a)
            parents_of_b = list(set(parents_of_b))
            parents_of_grp_b.extend(parents_of_grp_a)
            parents_of_grp_b = list(set(parents_of_grp_b))

            self._comp_objs_parents[comp_b] = parents_of_b
            self._group_parents[group_b] = parents_of_grp_b

            # new children list for parents_of_b
            for p in parents_of_b:
                children_of_p = self._comp_objs_children[p]
                children_of_p.append(comp_b)
                self._comp_objs_children[p] = list(set(children_of_p))
            for gp in parents_of_grp_b:
                children_of_gp = self._group_children[gp]
                children_of_gp.append(group_b)
                self._group_children[gp] = list(set(children_of_gp))

        return

    def __str__(self):
        return_str = "Final Group: " + self._name + "\n"
        for s in self._groups:
            return_str = return_str + s.__str__() + "\n"
        return return_str

    def create_storage_classes(self):
        return storage_classification(self._comp_objs)

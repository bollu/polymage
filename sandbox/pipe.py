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

def get_parents_from_func(func):
    refs = func.getObjects(Reference)
    # Filter out self and image references 
    refs = [ ref for ref in refs if not ref.objectRef == func and \
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

    for func in funcs:
        if func not in funcs_parents:
            funcs_parents[func] = []
        if func not in funcs_children:
            funcs_children[func] = []

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
    def __init__(self, _func, _is_output=False):
        assert isinstance(_func, Function)
        self._func = _func
        self._parents = []
        self._children = []

        self._group = None

        self._is_output = _is_output
        self._is_liveout = True
        self.set_flags()

        self._level_no = 0

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
    def is_image_typ(self):
        return self._is_image_typ

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
    def level(self):
        return self._level_no
    @property
    def is_output(self):
        return self._is_output
    @property
    def is_liveout(self):
        return self._is_liveout

    def set_flags(self):
        self._is_parents_set = False
        self._is_children_set = False
        self._is_group_set = False
        self._is_image_typ = isinstance(self.func, Image)
        return

    def add_child(self, comp):
        assert isinstance(comp, ComputeObject)
        self._children.append(comp)
        self._children = list(set(self._children))
        self._is_children_set = True
        return
    def add_parent(self, comp):
        assert isinstance(comp, ComputeObject)
        self._parents.append(comp)
        self._parents = list(set(self._parents))
        self._is_parents_set = True
        return

    def remove_child(self, comp):
        if comp in self._children:
            self._children.remove(comp)
        return
    def remove_parent(self, comp):
        if comp in self._parents:
            self._parents.remove(comp)
        return

    def set_parents(self, parents):
        # empty list of parents => root level comp
        if not parents:
            self._is_parents_set = True
            return
        for p in parents:
            assert isinstance(p, ComputeObject)
        self._parents = parents
        self._is_parents_set = True
        return

    def set_children(self, children):
        # empty list of children => leaf level comp
        if not children:
            self._is_children_set = True
            return
        for p in children:
            assert isinstance(p, ComputeObject)
        self._children = children
        self._is_children_set = True
        return

    def set_group(self, group):
        assert isinstance(group, Group)
        self._group = group
        self._is_group_set = True
        return

    def unset_group(self):
        self._group = None
        return

    def compute_liveness(self):
        # if there any children
        if not self.children:
            # no child => live_out
            self._is_liveout = True
            return

        self._is_liveout = False
        for child in self.children:
            # if any child is in another group
            if child.group != self.group:
                self._is_liveout = True
                break

        return

    def set_level(self, _level_no):
        self._level_no = _level_no


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
        self._is_image_typ = False
        for comp in _comp_objs:
            assert(isinstance(comp, ComputeObject))
            if comp.is_image_typ:
                self._is_image_typ = True

        self._comps  = _comp_objs
        self._parents = []
        self._children = []

        self.set_comp_group()

        self._level_order_comps = self.order_compute_objs()
        self._root_comps = self.find_root_comps()

        self._polyrep = None
        refs = []

        for comp in self.comps:
            refs += comp.func.getObjects(Reference)

        self._inputs = list(set([ref.objectRef for ref in refs \
                            if isinstance(ref.objectRef, Image)]))

        # Create a polyhedral representation if possible.
        # Currently doing extraction only when all the compute_objs
        # domains are affine. This can be revisited later.
        if self.isPolyhedral():
            self._polyrep = PolyRep(_ctx, self, [], _param_constraints)

    @property
    def comps(self):
        return self._comps
    @property
    def parents(self):
        return self._parents
    @property
    def children(self):
        return self._children

    @property
    def is_image_typ(self):
        return self._is_image_typ

    @property
    def polyRep(self):
        return self._polyrep
    @property
    def inputs(self):
        return self._inputs
    @property
    def name(self):
        return str([comp.func.name for comp in self.comps])

    @property
    def get_ordered_comps(self):  # <- cant have such a name for property
        return self._level_order_comps
    @property
    def root_comps(self):
        return self._root_comps

    def set_comp_group(self):
        for comp in self.comps:
            comp.set_group(self)
        return

    def find_and_set_parents(self):
        parents = []
        for comp in self.comps:
            comp_parent_groups = [p_comp.group for p_comp in comp.parents]
            parents.extend(comp_parent_groups)
        parents = list(set(parents))
        if self in parents:
            parents.remove(self)
        self.set_parents(parents)
        return parents
    def find_and_set_children(self):
        children = []
        for comp in self.comps:
            comp_children_groups = [c_comp.group for c_comp in comp.children]
            children.extend(comp_children_groups)
        children = list(set(children))
        if self in children:
            children.remove(self)
        self.set_children(children)
        return children

    def add_child(self, group):
        assert isinstance(group, Group)
        self._children.append(group)
        self._children = list(set(self._children))
        return
    def add_parent(self, group):
        assert isinstance(group, Group)
        self._parents.append(group)
        self._parents = list(set(self._parents))
        return

    def remove_child(self, group):
        if group in self._children:
            self._children.remove(group)
        return
    def remove_parent(self, group):
        if group in self._parents:
            self._parents.remove(group)
        return

    def set_parents(self, parents):
        for group in parents:
            assert isinstance(group, Group)
        self._parents = parents
        return
    def set_children(self, children):
        for group in children:
            assert isinstance(group, Group)
        self._children = children
        return

    def compute_liveness(self):
        for comp in self.comps:
            comp.compute_liveness()
            parts = self.polyRep.poly_parts[comp]
            for part in parts:
                part.compute_liveness()
        return

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

    def find_root_comps(self):
        root_comps = [comp for comp in self.comps \
                             if self._level_order_comps[comp] == 0]
        return root_comps

    def get_sorted_comps(self):
        sorted_comps = sorted(self._level_order_comps.items(),
                              key=lambda x: x[1])
        sorted_comps = [c[0] for c in sorted_comps]
        return sorted_comps

    def __str__(self):
        comp_str  = '[' + \
                    ', '.join([comp.func.name \
                        for comp in self.comps]) + \
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
        self._func_map, self._comps = \
            self.create_compute_objects()

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
        for comp in self.comps:
            inputs = inputs + comp.group.inputs
        self._inputs = list(set(inputs))

        # Checking bounds
        #bounds_check_pass(self)  # <-

        # inline pass
        inline_pass(self)  # <-

        # make sure the set of functions to be inlined and those to be grouped
        # are disjoint
        if self._inline_directives and self._grouping:
            group_comps = []
            for g in self._grouping:
                group_funcs += g
            a = set(self._inline_directives)
            b = set(group_funcs)
            assert a.isdisjoint(b)

        ''' GROUPING '''
        # TODO check grouping validity
        if self._grouping and False:
            # for each group
            for g in self._grouping:
                # get clones of all functions
                clones = [self._clone_map[f] for f in g]
                comps = [self.func_map[f] for f in clones]
                # list of group objects to be grouped
                merge_group_list = \
                    [self.groups[comp] for comp in comps]
                if len(merge_group_list) > 1:
                     merged = merge_group_list[0]
                     for i in range(1, len(merge_group_list)):
                        merged = self.merge_groups(merged, merge_group_list[i])
        else:
            # Run the grouping algorithm
            auto_group(self)
            pass

        # level order traversal of groups
        self._level_order_groups = self.order_group_objs()

        # update liveness of compute objects in each new group
        for group in self.groups:
            group.compute_liveness()

        # ***
        log_level = logging.INFO
        LOG(log_level, "Grouped compute objects:")
        for g in self.groups:
            LOG(log_level, g.name)
        # ***

        ''' SCHEDULING '''
        for g in self.groups:
            # alignment and scaling
            align_and_scale(self, g)
            # base schedule
            gparts = base_schedule(g)
            # grouping and tiling
            fused_schedule(self, g, self._param_estimates)

        self._grp_schedule = schedule_groups(self)

        # use graphviz to create pipeline graph
        self._pipeline_graph = self.draw_pipeline_graph()

        ''' STORAGE OPTIMIZATION '''
        '''
        # storage opt for liveout (full array) allocations
        # 1. derive liveout comps schedule from group schedule
        # 2. classify the storage based on type, dimensionality and size
        self.create_storage_classes()
        '''

    @property
    def func_map(self):
        return self._func_map
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
    def get_ordered_groups(self):  # <- naming
        return self._level_order_groups
    @property
    def group_schedule(self):
        return self._grp_schedule

    def get_parameters(self):
        params=[]
        for group in self._groups.values():
            params = params + group.getParameters()
        return list(set(params))

    def create_compute_objects(self):
        funcs, parents, children = \
            get_funcs_and_dep_maps(self._outputs)

        comps = []
        func_map = {}
        for func in funcs:
            output = False
            if func in self._outputs:
                output = True
            comp = ComputeObject(func, output)
            func_map[func] = comp
            comps.append(comp)

        # set parents, children information
        for func in func_map:
            comp = func_map[func]
            # set parents
            comp_parents = [func_map[p_func] for p_func in parents[func]]
            comp.set_parents(comp_parents)
            # set children
            comp_children = [func_map[c_func] for c_func in children[func]]
            comp.set_children(comp_children)

        return func_map, comps

    def order_compute_objs(self):
        parents = {}
        for comp in self.comps:
            parents[comp] = comp.parents
        order = level_order(self.comps, parents)
        for comp in order:
            comp.set_level(order[comp])
        return order

    def order_group_objs(self):
        parents = {}
        for group in self.groups:
            parents[group] = group.parents
        order = level_order(self.groups, parents)
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
        groups = []
        for comp in comps:
            group = Group(self._ctx, [comp], self._param_constraints)
            groups.append(group)

        for group in groups:
            group.find_and_set_parents()
            group.find_and_set_children()

        return groups

    def draw_pipeline_graph(self):
        gr = pgv.AGraph(strict=False, directed=True)

        # TODO add input nodes to the graph
        for i in range(0, len(self.groups)):
            sub_graph_nodes = []
            for comp in self.comps:
                if comp.group == self.groups[i]:
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

    def drop_comp(self, comp):
        # if the compute object is a child of any other
        if comp.parents:
            for p_comp in comp.parents:
                p_comp.remove_child(comp)
        # if the compute object is a parent of any other
        if comp.children:
            for c_comp in comp.children:
                c_comp.remove_parent(comp)
        # remove comp_obj
        self._comps.remove(comp)
        func = comp.func
        self._func_map.pop(func)

        return

    def add_group(self, group):
        """
        add a new group to the pipeline
        """
        if not group.comps:
            return

        self._groups.append(group)

        group.find_and_set_parents()
        group.find_and_set_children()

        # add group as child for all its parents
        for p_group in group.parents:
            p_group.add_child(group)
        # add group as parent for all its children
        for c_group in group.children:
            c_group.add_parent(group)

        return

    def drop_group(self, group):
        """
        drop the group from the pipeline
        """
        # if group is a child of any other group
        if group.parents:
            for p_group in group.parents:
                p_group.remove_child(group)
        # if group is a parent of any other group
        if group.children:
            for c_group in group.children:
                c_group.remove_parent(group)
        self._groups.remove(group)
        for comp in group.comps:
            comp.unset_group()

        return

    def merge_groups(self, g1, g2):
        # Get comp objects from both groups
        comps = g1.comps + g2.comps
        comps = list(set(comps))

        self.drop_group(g1)
        self.drop_group(g2)

        # Create a new group
        merged = Group(self._ctx, comps,
                       self._param_constraints)

        self.add_group(merged)

        return merged

    def replace_group(self, old_group, new_group):
        # if old_group has any child
        if old_group.children:
            for child in old_group.children:
                child.add_parent(new_group)
                new_group.add_child(child)
        # if old_group has any parent
        if old_group.parents:
            for parent in old_group.parents:
                parent.add_child(new_group)
                new_group.add_parent(parent)

        # replace old_group with new_group in groups list
        comp = old_group.comps[0]
        self.drop_group(old_group)
        self._groups[comp] = new_group

        return

    def make_func_independent(self, func_a, func_b):
        """
        makes func_b independent of func_b and updates parent children
        relations in the graph structure
        [ assumes that func_b is a child of func_a, and that func_a is inlined
        into func_b ]
        """
        comp_a = self.func_map[func_a]
        comp_b = self.func_map[func_b]
        group_a = comp_a.group
        group_b = comp_b.group
        # if parent_comp has any parent
        if comp_a.parents:
            parents_of_a = comp_a.parents
            parents_of_b = comp_b.parents
            parents_of_grp_a = group_a.parents
            parents_of_grp_b = group_b.parents

            # remove relation between a and b
            comp_a.remove_child(comp_b)
            group_a.remove_child(group_b)

            parents_of_b.remove(comp_a)
            parents_of_grp_b.remove(group_a)

            # new parents list for b
            # compute object
            parents_of_b.extend(parents_of_a)
            parents_of_b = list(set(parents_of_b))
            # group object
            parents_of_grp_b.extend(parents_of_grp_a)
            parents_of_grp_b = list(set(parents_of_grp_b))

            comp_b.set_parents(parents_of_b)
            group_b.set_parents(parents_of_grp_b)

            # new children list for parents_of_b
            for p_comp in parents_of_b:
                p_comp.add_child(comp_b)
            for p_group in parents_of_grp_b:
                p_group.add_child(group_b)

        return

    def __str__(self):
        return_str = "Final Group: " + self._name + "\n"
        for s in self._groups:
            return_str = return_str + s.__str__() + "\n"
        return return_str

    def create_storage_classes(self):
        return storage_classification(self.comps)

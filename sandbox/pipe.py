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
            assert(not isinstance(comp, Image))

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
    def computeObjs(self):
        return self._comp_objs
    @property
    def polyRep(self):
        return self._polyrep
    @property
    def inputs(self):
        return self._inputs

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
                 _grouping, _options, _name = None):
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
                    ref._replaceRefObject(self._clone_map[ref.objectRef])

        ''' DAG OF CLONES '''
        self._comp_objs, \
         self._comp_objs_parents, \
          self._comp_objs_children = \
            get_comp_objs_and_dep_maps(self._outputs)

        ''' INITIAL GROUPING '''
        # Create a group for each pipeline function / reduction and compute
        # maps for parent and child relations between the groups
        self._groups, \
         self._groupParents, \
          self._groupChildren = \
            self.buildInitialGroups(self._comp_objs,
                                    self._comp_objs_parents,
                                    self._comp_objs_children)

        # Store the initial pipeline graph. The compiler can modify 
        # the pipeline by inlining functions.
        #self._initialGraph = self.drawPipelineGraph()

        # Make a list of all the input groups
        inputs = []
        for f in self._groups:
            inputs = inputs + self._groups[f].inputs
        self._inputs = list(set(inputs))

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
                        # to be done after each merge, to know if the
                        # merging was valid.
                        align_and_scale_parts(self, merged)
        else:
            # Run the grouping algorithm
            pass

        ''' BASE SCHEDULE AND CODEGEN '''
        for g in list(set(self._groups.values())):
            gparts = base_schedule(g)
            for p in gparts:
                p.liveout = True
            #fused_schedule(g)  # (g, param_estimates)
            #g.polyRep.generate_code()

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
    
    def buildInitialGroups(self, compObjs, compObjsParents, compObjsChildren):
        """
        Create a pipeline graph where the nodes are a group and the edges
        represent the dependences between the groups. The initial graph
        has each compute object in a separate group. 
        """
        # TODO correct the comment
        groupMap = {}
        groupParents = {}
        groupChildren = {}
        for comp in compObjs:
            groupMap[comp] = Group(self._ctx, [comp], self._param_constraints)

        for comp in groupMap:
            groupParents[groupMap[comp]] = [ groupMap[p] for p in \
                                                        compObjsParents[comp] ]
            groupChildren[groupMap[comp]] = [ groupMap[c] for c in \
                                                        compObjsChildren[comp] ]

        return groupMap, groupParents, groupChildren
    
    def getParameters(self):
        params=[]
        for group in self._groups.values():
            params = params + group.getParameters()
        return list(set(params))

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

    def generate_code(self):
        return generate_code_for_pipeline(self)

    def merge_groups(self, g1, g2):
        # Get comp objects from both groups 
        comp_objs = g1.computeObjs + g2.computeObjs
        comp_objs = list(set(comp_objs))
        # Create a new group 
        merged = Group(self._ctx, comp_objs, self._param_constraints)
        # Update the group map
        for comp in comp_objs:
            self._groups.pop(comp)
            self._groups[comp] = merged

        # Update the group parent map
        parents = self._groupParents[g1] + self._groupParents[g2]
        parents = [ p for p in parents if p is not g1 and p is not g2 ]
        parents = list(set(parents))

        self._groupParents.pop(g1)
        self._groupParents.pop(g2)

        for g in self._groupParents:
            if g1 in self._groupParents[g] or g2 in self._groupParents[g]:
                self._groupParents[g].append(merged)
            if g1 in self._groupParents[g]:
                self._groupParents[g].remove(g1)
            if g2 in self._groupParents[g]:
                self._groupParents[g].remove(g2)
        
        self._groupParents[merged] = parents

        # update the group child map
        children = self._groupChildren[g1] + self._groupChildren[g2]
        children = [ c for c in children if c is not g1 and c is not g2 ]
        children = list(set(children))

        self._groupChildren.pop(g1)
        self._groupChildren.pop(g2)

        for g in self._groupChildren:
            if g1 in self._groupChildren[g] or g2 in self._groupChildren[g]:
                self._groupChildren[g].append(merged)
            if g1 in self._groupChildren[g]:
                self._groupChildren[g].remove(g1)
            if g2 in self._groupChildren[g]:
                self._groupChildren[g].remove(g2)
        
        self._groupChildren[merged] = children

        return merged

    def get_ordered_groups(self):
        # Assign level numbers to each group and sort accourding to the level
        group_order = {}
        groups = set(self._groups.values())
        group_list = [ [g, len(self._groupParents[g])] for g in groups ]

        level = 0
        while group_list:
            # find all the groups whose parents have their levels assigned
            level_assigned = [ t for t in group_list if t[1] == 0 ]
            for assgn in level_assigned:
                group_order[assgn[0]] = level
                group_list.remove(assgn)
                # reduce the unassigned parent count for all the children
                child_groups = self._groupChildren[assgn[0]]
                for assgn in group_list:
                    if assgn[0] in child_groups:
                        assgn[1] -= 1
            level = level + 1
        
        return sorted(group_order.items(), key=lambda x: x[1])

    def boundsCheckPass(self):
        """ 
            Bounds check pass analyzes if function values used in the compute
            objects are within the domain of the functions. Static analysis is
            only possible when the references to function values are regular
            i.e. they are not data dependent. We restrict ourselves to affine
            references.
        """
        for group in self._groups.values():
            for child in group.childGroups:
                checkRefs(child, group)              
            for inp in group.inputs:
                # Creating a computation group for an input which is given
                # is meaningless. Ideally it should be done in a clean way
                # currently abusing group for construction of a polyhedral
                # representation
                inpGroup = Group([inp], self._ctx, self._param_constraints,
                                 self._param_estimates, self._tileSizes,
                                 self._sizeThreshold, self._groupSize,
                                 self._outputs)
                checkRefs(group, inpGroup)

    def inlinePass(self):
        """ 
            Inline pass takes all the inlining decisions and inlines functions 
            at their use points in other groups.
        """
        # Guidelines for inlining
        # -- Inlining functions which have multiple case constructs can quickly 
        #    turn out to be messy code generation nightmare.
        # -- Splitting an use site can occur when the inlined function is defined 
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
        # -- Users should be enouraged to write functions in a form that makes
        #    inlining easy.
        inlinedCompObjs = []
        for directive in self._inlineDirectives:
            # Only function constructs can be inlined for now
            assert isinstance(directive, Function)
            # Does inling into a fused group cause problems?
            parentGroup = self._groups[directive]
            assert parentGroup.computeObjs[0] not in self._outputs
            for child in parentGroup.childGroups:
                refToInlineExprMap = inline(child, parentGroup, noSplit = True)
                child.computeObjs[0].inlineRefs(refToInlineExprMap)
            # Recompute group graph
            self._groups = self.buildGroupGraph()

    # TODO printing the pipeline 
    def __str__(self):
        return_str = "Final Group: " + self._name + "\n"
        for s in self._groups:
            return_str = return_str + s.__str__() + "\n"
        return return_str

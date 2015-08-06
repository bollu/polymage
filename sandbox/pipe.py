from __future__ import absolute_import, division, print_function

# More Python 3 vs 2 mojo
try:
    import queue
except ImportError:
    import Queue as queue

import pygraphviz as pgv

from constructs import *
from codegen import *
from schedule import *
from poly import *

def getParentsFromCompObj(comp):
    refs = comp.getObjects(Reference)
    # Filter out self and image references 
    refs = [ ref for ref in refs if not ref.objectRef == comp and \
                                 not isinstance(ref.objectRef, Image) ]
    return list(set([ref.objectRef for ref in refs]))

def getCompObjsAndDependenceMaps(outputs):
    """
    Find all the compute objects required for the outputs and
    also builds parent and children maps for the compute objects
    """
    compObjs = []
    compObjsParents = {}
    compObjsChildren = {}
    q = queue.Queue()
    for compObj in outputs:
        q.put(compObj)
    while not q.empty():
        obj = q.get()
        parentObjs = getParentsFromCompObj(obj)
        if obj not in compObjs:
            compObjs.append(obj)
            compObjsParents[obj] = parentObjs
            for parobj in parentObjs:
                if parobj in compObjsChildren:
                    if obj not in compObjsChildren[parobj]:
                        compObjsChildren[parobj].append(obj)
                else:
                    compObjsChildren[parobj] = [obj]                        
            if len(parentObjs) != 0:
                for r in parentObjs:
                    q.put(r)

    for comp in compObjs:
        if comp not in compObjsParents:
            compObjsParents[comp] = []
        if comp not in compObjsChildren:
            compObjsChildren[comp] = []

    return compObjs, compObjsParents, compObjsChildren

class Group:
    """ 
        Group is a part of the pipeline which realizes a set of computation
        objects. Scheduling and storage allocation is done at the level of a 
        group. A group also maintains a polyhedral representation of the 
        computation objects when possible.
    """
    # Construct a group from a set of language functions / reductions
    def __init__(self, _ctx, _compObjs, _paramConstraints):
        # All the computation constructs in the language derive from the
        # Function class. Input images cannot be part of a group.
        for comp in _compObjs:
            assert(isinstance(comp, Function))
            assert(not isinstance(comp, Image))

        self._compObjs  = _compObjs
        self._polyrep = None
        refs = []

        for comp in self._compObjs:
            refs += comp.getObjects(Reference)

        self._inputs = list(set([ref.objectRef for ref in refs \
                            if isinstance(ref.objectRef, Image)]))

        # Create a polyhedral representation if possible
        polyhedral = True

        # Currently doing extraction only when all the computeObjs
        # domains are affine. This can be revisited later. 
        for comp in self._compObjs:
            if (not comp.hasBoundedIntegerDomain()):
                polyhedral = False
        if polyhedral:
            self._polyrep = PolyRep(_ctx, self, _paramConstraints) 
        
    @property
    def computeObjs(self):
        return self._compObjs
    @property
    def polyRep(self):
        return self._polyrep
    @property
    def inputs(self):
        return self._inputs

    def getParameters(self):
        params = []
        for comp in self._compObjs:
            params = params + comp.getObjects(Parameter)
        return list(set(params))

    def isPolyhedral(self):
        polyhedral = True
        for comp in self._compObjs:
            if (not comp.hasBoundedIntegerDomain()):
                polyhedral = False
        return polyhedral

    def orderComputeObjs(self):
        # Order stores the numbering of each compute object 
        # when topologically sorted.
        order = {}
        # Initialize all the initial numbering to zero for
        # all compute objects in the group
        for comp in self._compObjs:
            order[comp] = 0
        # Doing a topological sort in an iterative fashion
        change = True
        while(change):
            change = False
            for comp in self._compObjs:
                parentObjs = getParentsFromCompObj(comp)
                for pobj in parentObjs:
                    if (pobj in order  and (order[pobj] >= order[comp])):
                        order[comp] = order[pobj] + 1
                        change = True
        return order

    def __str__(self):
        comp_str  = "\n\n".join([comp.__str__() \
                    for comp in self._compObjs]) + '\n'
        return comp_str + '\n' + self._polyrep.__str__()

class Pipeline:
    def __init__(self, _ctx, _outputs, \
                 _paramConstraints, _grouping, \
                 _options, _name=None):
        # Name of the pipleline is a concatenation of the names of the 
        # pipeline outputs, unless it is explicitly named.
        if _name is None:
            _name = ""
            for out in _outputs:
                _name = _name + out.name

        self._name   = _name

        self._ctx = _ctx
        self._orgOutputs = _outputs
        self._paramConstraints = _paramConstraints
        self._grouping = _grouping
        self._options = _options

        ''' CONSTRUCT DAG '''
        # Maps from a compute object to its parents and children by
        # backtracing starting from given live-out functions.
        # TODO: see if there is a cyclic dependency between the compute
        # objects. Self references are not treated as cycles.
        self._orgCompObjs, \
         _, \
          _ = \
            getCompObjsAndDependenceMaps(self._orgOutputs)

        ''' CLONING '''
        # Clone the computation objects i.e. functions and reductions
        self._cloneMap = {}
        for comp in self._orgCompObjs:
            self._cloneMap[comp] = comp.clone()
        self._outputs = [self._cloneMap[obj] for obj in self._orgOutputs]
        # Modify the references in the cloned objects (which refer to 
        # the original objects)  
        for comp in self._orgCompObjs:
            cln = self._cloneMap[comp]
            refs = cln.getObjects(Reference)
            for ref in refs:
                if not isinstance(ref.objectRef, Image): 
                    ref._replaceRefObject(self._cloneMap[ref.objectRef])

        ''' DAG OF CLONES '''
        self._compObjs, self._compObjsParents, self._compObjsChildren = \
                                getCompObjsAndDependenceMaps(self._outputs)

        ''' INITIAL GROUPING '''
        # Create a group for each pipeline function / reduction and compute
        # maps for parent and child relations between the groups
        self._groups, \
         self._groupParents, \
          self._groupChildren = \
                              self.buildInitialGroups(self._compObjs,
                                                      self._compObjsParents,
                                                      self._compObjsChildren)

        # Store the initial pipeline graph. The compiler can modify 
        # the pipeline by inlining functions.
        self._initialGraph = self.drawPipelineGraph()

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
                merge_group_list = [self._groups[self._cloneMap[f]] for f in g]
                if len(merge_group_list) > 1:
                     merged = merge_group_list[0]
                     for i in range(1, len(merge_group_list)):
                        merged = self.merge_groups(merged, merge_group_list[i])
                        # to be done after each merging, to know if the
                        # merging was valid
                        align_and_scale_parts(self, merged)
        else:
            # Run the grouping algorithm
            pass

        ''' BASE SCHEDULE AND CODEGEN '''
        for g in list(set(self._groups.values())):
            baseSchedule(g)
            g.polyRep.generateCode()

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
            groupMap[comp] = Group(self._ctx, [comp], self._paramConstraints)

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

    def drawPipelineGraph(self):
        G = pgv.AGraph(strict=False, directed=True)
        groupList = list(set([self._groups[f] for f in self._groups]))

        # TODO add input nodes to the graph
        for i in range(0, len(groupList)):
            subGraphNodes = []
            for obj in self._compObjs:
                if self._groups[obj] == groupList[i]:
                    subGraphNodes.append(obj.name)
            G.add_nodes_from(subGraphNodes)
            G.add_subgraph(nbunch = subGraphNodes, 
                           name = "cluster_" + str(i))

        for obj in self._compObjs:
            for pobj in self._compObjsParents[obj]:
                G.add_edge(pobj.name, obj.name)

        G.layout(prog='dot')
        return G

    def generateCode(self):
        return generate_code_for_pipeline(self)

    def merge_groups(self, g1, g2):
        # Get comp objects from both groups 
        comp_objs = g1.computeObjs + g2.computeObjs
        # Create a new group 
        merged = Group(self._ctx, comp_objs, self._paramConstraints)
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

    def getOrderedGroups(self):
        # Assign level numbers to each group and sort accourding to the 
        # level
        groupOrder = {}
        groups = set(self._groups.values())
        groupList = [ [g, len(self._groupParents[g])] for g in groups ]

        level = 0
        while groupList:
            # find all the groups whose parents have their levels assigned
            levelAssigned = [ t for t in groupList if t[1] == 0 ]
            for assgn in levelAssigned:
                groupOrder[assgn[0]] = level
                groupList.remove(assgn)
                # reduce the unassigned parent count for all the children
                childGroups = self._groupChildren[assgn[0]]
                for assgn in groupList:
                    if assgn[0] in childGroups:
                        assgn[1] = assgn[1] - 1
            level = level + 1
        
        return sorted(groupOrder.items(), key=lambda x: x[1])

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
                inpGroup = Group([inp], self._ctx, self._paramConstraints,
                                 self._paramEstimates, self._tileSizes,
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

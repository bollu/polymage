# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

import sys
from Constructs import *

# More Python 3 vs 2 mojo
try:
    import queue
except ImportError:
    import Queue as queue

import pygraphviz as pgv

# Function to construct a stage from a language function / reduction

class Stage:
    """ 
        Stage is a part of the pipeline which realizes a set of computation
        objects. Scheduling and storage allocation for the computation objects
        is done at the level of a stage. A stage also maintains a polyhedral 
        representation of the computation objects when possible.
    """
    def __init__(self, _computeObjs, _ctx, _paramConstraints, _paramEstimates,
                 _tileSizes, _sizeThreshold, _groupSize, _outputs):
        for comp in _computeObjs:
            assert(isinstance(comp, Function))

        self._computeObjs  = _computeObjs
        self._childStages = []
        self._parentStages = []
        self._outputs = _outputs
        refs = []

        for comp in self._computeObjs:
            # Filter out self references
            currRefs = comp.getObjects(Reference)
            currRefs = [ ref for ref in currRefs if not ref.objectRef == comp ]
            refs = refs + currRefs

        self._parentObjs = list(set([ref.objectRef for ref in refs \
                                if not isinstance(ref.objectRef, Image)]))
        self._inputs = list(set([ref.objectRef for ref in refs \
                            if isinstance(ref.objectRef, Image)]))
        # Create a polyhedral representation
        self._polyrep = opt.PolyRep(_ctx, self, _paramConstraints,
                                    _paramEstimates, _tileSizes,
                                    _sizeThreshold, _groupSize, _outputs)

    @property
    def childStages(self):
        return self._childStages
    @property
    def parentStages(self):
        return self._parentStages
    @property
    def parentObjs(self):
        return self._parentObjs
    @property
    def computeObjs(self):
        return self._computeObjs
    @property
    def inputs(self):
        return self._inputs
    @property
    def polyRep(self):
        return self._polyrep

    def getParameters(self):
        params = []
        for comp in self._computeObjs:
            params = params + comp.getObjects(Parameter)
        return list(set(params))

        #return list(set([comp.getObjects(Parameter) for comp in self._computeObjs]))

    def isFused(self):
        return len(self._computeObjs) > 1

    def isPolyhedral(self):
        polyhedral = True
        for comp in self._computeObjs:
            if (not comp.hasBoundedIntegerDomain()):
                polyhedral = False
        return polyhedral

    def orderComputeObjs(self):
        order = {}
        for comp in self.computeObjs:
            order[comp] = 0
        # Doing a topological sort the easy way. There is a topological
        # sort else where in the code hopefully the can be unifed at
        # some point.
        change = True
        while(change):
            change = False
            for comp in self.computeObjs:
                # filter self references
                refs = comp.getObjects(Reference)
                refs = [ref for ref in refs if not ref.objectRef == comp]
                parentObjs = list(set([ref.objectRef for ref in refs\
                                  if not isinstance(ref.objectRef, Image)]))
                for pobj in parentObjs:
                    if (pobj in order  and (order[pobj] >= order[comp])):
                        order[comp] += 1
                        change = True
        return order

    def createLoopVariables(self, variables):
        cvarMap = {}
        # Create iterator variables and bind them to DSL variables
        for i in xrange(0, len(variables)):
            varType = Cgen.cgenType.get(variables[i].typ)
            var = Cgen.cVariable(varType, Cgen.cNameGen.getIteratorName())
            # Binding variables to C iterators
            cvarMap[variables[i]] = var
        return cvarMap

    def createPerfectNestedLoop(self, body, variables, domains,
                                cfuncMap, cparamMap, cvarMap):
        lbody = body
        for i in xrange(0, len(variables)):    
            var = cvarMap[variables[i]]
            # Convert lb and ub expressions to C expressions
            lb  = generateCExpr(domains[i].lowerBound, 
                                cparamMap, cvarMap, cfuncMap)                        
            ub  = generateCExpr(domains[i].upperBound, 
                                cparamMap, cvarMap, cfuncMap)                        
            step  = generateCExpr(domains[i].step, 
                                  cparamMap, cvarMap, cfuncMap)                        
            varDecl =  Cgen.cDeclaration(var.typ, var, lb)
            comp = '<='
            if domains[i].step.value < 0:
                comp = '>='
            cond = Cgen.cCond(var, comp, ub)
            inc = Cgen.cAssign(var, var+step)
            loop = Cgen.cFor(varDecl, cond, inc)
            lbody.add(loop, False)
            lbody = loop.body
        return lbody     

    def generateFunctionScanLoops(self, obj, body, cparamMap, cfuncMap):
        # Compute function points in lexicographic order of domain
        cvarMap = self.createLoopVariables(obj.variables)

        # Generate loops. lbody is the body of the innermost loop.
        lbody = self.createPerfectNestedLoop(body, obj.variables, 
                                             obj.domain, cfuncMap, 
                                             cparamMap, cvarMap)
        arglist = obj.variables
        # Convert function definition into a C expression and add it to loop body
        #hasExpr = False
        for case in obj.definition:
            if(isinstance(case, AbstractExpression)):
                caseexpr = generateCExpr(case, cparamMap, 
                                         cvarMap, cfuncMap)
                arrayRef = generateCExpr(obj(*arglist), cparamMap,
                                         cvarMap, cfuncMap)
                assign = Cgen.cAssign(arrayRef, caseexpr)
                lbody.add(assign, False)
                #hasExpr = True
            elif(isinstance(case, Case)):
                ccond = generateCCond(case.condition, cparamMap,
                                      cvarMap, cfuncMap)
                condexpr = generateCExpr(case.expression, cparamMap,
                                         cvarMap, cfuncMap)
                cif = Cgen.cIfThen(ccond)

                if (isinstance(case.expression, AbstractExpression)):
                    arrayRef = generateCExpr(obj(*arglist), cparamMap,
                                         cvarMap, cfuncMap)
                    assign = Cgen.cAssign(arrayRef, condexpr)
                    with cif.ifBlock as ifblock:
                        ifblock.add(assign)
                        #ifblock.add(Cgen.cContinue())
                else:
                    assert False
                lbody.add(cif, False)
            else:
                assert False
        #if not hasExpr:
        #    caseexpr = generateCExpr(obj.default, cparamMap, 
        #                             cvarMap, cfuncMap)
        #    assign = Cgen.cAssign(array(*arglist), caseexpr)
        #    lbody.add(assign, False)

    def generateAccumulatorScanLoops(self, obj, body, cparamMap, cfuncMap):
        # Compute accumulator points in lexicographic order of reduction domain
        cvarMap = self.createLoopVariables(obj.reductionVariables)

        # Generate loops. lbody is the body of the innermost loop.
        lbody = self.createPerfectNestedLoop(body, obj.reductionVariables, 
                                             obj.reductionDomain, cfuncMap, 
                                             cparamMap, cvarMap)

        # Convert function definition into a C expression and add it to loop body
        #hasExpr = False
        for case in obj.definition:
            if(isinstance(case, Accumulate)):
                caseexpr = generateCExpr(case.expression, cparamMap, 
                                         cvarMap, cfuncMap)
                refArgs = case.accumulateRef.arguments
                accumRef = generateCExpr(obj(*refArgs), cparamMap,
                                             cvarMap,cfuncMap)
                assign = Cgen.cAssign(accumRef, accumRef + caseexpr)                         
                lbody.add(assign, False)
                #hasExpr = True
            elif(isinstance(case, Case)):
                ccond = generateCCond(case.condition, cparamMap,
                                      cvarMap, cfuncMap)
                condexpr = generateCExpr(case.expression, cparamMap,
                                         cvarMap, cfuncMap)
                cif = Cgen.cIfThen(ccond)

                if(isinstance(case.expression, Accumulate)):
                    refArgs = case.accumulateRef.arguments
                    accumRef = generateCExpr(obj(*refArgs), cparamMap,
                                             cvarMap,cfuncMap)
                    assign = Cgen.cAssign(accumRef, accumRef + condexpr)                         
                    with cif.ifBlock as ifblock:
                        ifblock.add(assign)
                        #ifblock.add(Cgen.cContinue())
                else:
                    assert False
                lbody.add(cif, False)
            else:
                assert False
        #if not hasExpr:
        #    caseexpr = generateCExpr(obj.default, cparamMap, 
        #                             cvarMap, cfuncMap)
        #    assign = Cgen.cAssign(array(*arglist), caseexpr)
        #    lbody.add(assign, False)    

    def generateCNaive(self, body, cfuncMap, cparamMap, outputs, isExternAlloc=False):
        self._polyrep.generateCode()
        compObjs = self.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        freeList = []

        for obj in [item[0] for item in sortedCompObjs]:
            if(isinstance(obj, Image)):
                continue
            # Allocate storage
            arraydim = len(obj.variables)
            # Conservative storage estimation for each dimension
            # - Check if interval step is increasing or decreasing and pick the 
            #   bound expression accordingly
            notLiveOut = not obj in outputs
            if obj in self.polyRep.polyParts:
                for part in self.polyRep.polyParts[obj]:
                   notLiveOut = notLiveOut and not part.liveout
            else:
                notLiveOut = False

            redDims = [ -1 for i in xrange(0, len(obj.domain))]
            scratch = [ False for i in xrange(0, len(obj.domain))]
            if obj in self.polyRep.polyParts and notLiveOut:
                for part in self.polyRep.polyParts[obj]:
                    for i in xrange(0, len(obj.domain)):
                        if i in part.dimScratchSize:
                            redDims[i] = max(redDims[i], part.dimScratchSize[i])
                            scratch[i] = True
            dims = []
            for i in xrange(0, len(obj.domain)):
                interval = obj.domain[i]
                if redDims[i] == -1:
                    if interval.step.value > 0 :
                        dims.append(simplifyExpr(interval.upperBound - 
                                                 interval.lowerBound + 1))
                    else:
                        dims.append(simplifyExpr(interval.lowerBound - 
                                                interval.upperBound + 1))
                else:
                    dims.append(redDims[i])

            # Simplify the expressions for the array dimensions
            arrayType = Cgen.cgenType.get(obj.typ)
            array = Cgen.cArray(arrayType, obj.name, dims)
            arrayPtr = Cgen.cPointer(arrayType, 1)
            cfuncMap[obj] = (array, scratch)

            if obj not in outputs and not notLiveOut:
                arrayDecl = Cgen.cDeclaration(arrayPtr, array)
                body.add(arrayDecl)

            if not notLiveOut:
                # do not allocate for output arrays if they are
                # already allocated
                if obj in outputs and isExternAlloc:
                    # FIXME: accessing a class member directly,
                    # as a makeshift alternative
                    array.layout = 'contigous'
                    pass
                else:
                    array.allocate_contigous(body)

            if obj not in outputs and not notLiveOut:
                freeList.append(array)

            if obj in self.polyRep.polyParts:
                continue
            if type(obj) == Function:
                self.generateFunctionScanLoops(obj, body, cparamMap, cfuncMap)
            elif type(obj) == Accumulator:
                self.generateAccumulatorScanLoops(obj, body, cparamMap, cfuncMap)
            else:
                assert False

        if self.polyRep.polyast != []:
            for ast in self.polyRep.polyast:
                generateCNaiveFromIslAst(ast, body, cfuncMap, cparamMap)

        # Deallocate storage
        for array in freeList:
            array.deallocate(body)

    def __str__(self):
        comp_str  = "\n\n".join([comp.__str__() for comp in self._computeObjs]) + '\n'
        parent_str = "Parent Objects: " + \
                     ", ".join([ item.name for item in self._parentObjs]) + '\n' 
        return comp_str + '\n' + parent_str + '\n' + self._polyrep.__str__()

class PipeLine:
    def __init__(self, _outputs):
        # Name of the pipleline is a concatenation of the names of the pipeline
        # outputs.
        _name = ""
        for out in _outputs:
            _name = _name + out.name
        self._name   = _name
        
        self._outputs = _outputs
        
        # Create a stage for each pipeline input / function / reduction
        self._stages = self.createStages()

        # Determine the parent and child stages for each stage
        # _parents is a map from a stage to a list of parent stages
        # _children is a map from a stage to a list of child stages

        self._parents, self._children = self.buildStageGraph();
        self._initialGraph = self.drawPipelineGraph()

        # Make a list of all the input stages
        inputs = []
        for s in self._stages:
            if s.isInput(): 
                inputs = inputs + self._stages[s]

        self._inputs = inputs

    @property
    def stages(self):
        return self._stages

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
    def islContext(self):
        return self._ctx

    @property
    def graph(self):
        return self._intialGraph

    def getParameters(self):
        params=[]
        for stage in self._stages.values():
            params = params + stage.getParameters()
        return list(set(params))

    def drawPipelineGraph(self):
        G = pgv.AGraph(strict=False, directed=True)
        for f in self._stages:
            s = self._stages[f] 
            assert len(s.computeObjs) == 1
            for obj in s.computeObjs:
                for pobj in s.parentObjs:
                    G.add_edge(pobj.name, obj.name)

        G.layout(prog='dot')
        return G

    def drawPipelineGraphWithGroups(self):
        G = pgv.AGraph(strict=False, directed=True)
        for f in self._stages:
            s = self._stages[f]
            for i in xrange(0, len(s.polyRep.groups)):
                subGraphNodes = []
                for p in s.polyRep.groups[i]:
                    if p.comp.name not in subGraphNodes:
                        subGraphNodes.append(p.comp.name)
                G.add_nodes_from(subGraphNodes)
                G.add_subgraph(nbunch = subGraphNodes, 
                               name = "cluster_" + str(i))

        # Temporary hack getting edges from the orginal graph
        # this may not reflect the current pipeline due to 
        # inlining.
        G.add_edges_from(self.graph.edges())
        G.layout(prog='dot')
        return G

    def buildStageGraph(self):
        # Find all the computation objects that are required for the pipeline,
        # create stages for the computation and the depedency graph. This step 
        # assumes that there are no cycles in the pipeline stage graph this has 
        # assumption has to be revisited when dealing with more complex pipelines.
        stages = {}
        q = queue()
        for compObj in self._outputs:
            q.put(compObj)
        while not q.empty():
            obj = q.get()
            if obj not in stages:
                stages[obj] = Stage([obj], self._ctx, self._paramConstraints,
                                    self._paramEstimates, self._tileSizes,
                                    self._sizeThreshold, self._groupSize, 
                                    self._outputs)                
                if len(stages[obj].parentObjs) != 0:
                    for r in stages[obj].parentObjs:
                        q.put(r)
        
        for obj in stages:
            for pobj in stages[obj].parentObjs:
                stages[pobj].childStages.append(stages[obj])
                stages[obj].parentStages.append(stages[pobj])
        
        return stages

    def boundsCheckPass(self):
        """ Bounds check pass analyzes if function values used in the compute 
            objects are within the domain of the functions. Static analysis 
            is only possible when the references to function values are regular
            i.e. they are not data dependent. We restrict ourselves to affine
            references."""
        for stage in self._stages.values():
            for child in stage.childStages:
                opt.checkRefs(child, stage)              
            for inp in stage.inputs:
                # Creating a computation stage for an input which is given
                # is meaningless. Ideally it should be done in a clean way
                # currently abusing stage for construction of a polyhedral
                # representation
                inpStage = Stage([inp], self._ctx, self._paramConstraints,
                                 self._paramEstimates, self._tileSizes,
                                 self._sizeThreshold, self._groupSize,
                                 self._outputs)
                opt.checkRefs(stage, inpStage)

    def inlinePass(self):
        """ Inline pass takes all the inlining decisions and inlines functions 
            at their use points in other stages."""
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
            # Does inling into a fused stage cause problems?
            parentStage = self._stages[directive]
            assert parentStage.computeObjs[0] not in self._outputs
            for child in parentStage.childStages:
                refToInlineExprMap = opt.inline(child, parentStage, noSplit = True)
                child.computeObjs[0].inlineRefs(refToInlineExprMap)
            # Recompute stage graph
            self._stages = self.buildStageGraph()

    def getOrderedStages(self):
        # Topological sorting of stages.
        # Quick and dirty way of doing things have to revisit
        stageOrder = {}
        stageList = [ (len(s.parentStages), s) for s in self._stages.values()]
        level = 0
        while len(stageList) > 0:
            levelAssigned = [ s for s in stageList if s[0] == 0]
            for assigned in levelAssigned:
                stageOrder[assigned[1]] = level
                stageList.remove(assigned)
                def decChild(s):
                    if s[1] in assigned[1].childStages:
                        return (s[0]-1, s[1])
                    else: 
                        return s
                stageList = [ decChild(s) for s in stageList ] 
            level = level + 1
        return sorted(stageOrder.items(), key=lambda s: s[1])
   
    def generateCNaive(self, isExternFunc=False, areParamsVoidPtrs=False, isExternAlloc=False):
        # Order stages
        stageOrder = self.getOrderedStages()
        stageOrder = [ stage[0] for stage in stageOrder ]
        # Generate header and define block
        m = Cgen.cModule('Pipeline')
        with m.includes as incblock:
            incblock.add(Cgen.cInclude('stdio.h'))
            incblock.add(Cgen.cInclude('stdlib.h'))
            incblock.add(Cgen.cInclude('malloc.h'))
            incblock.add(Cgen.cInclude('cmath'))
            incblock.add(Cgen.cInclude('string.h'))
            incblock.add(Cgen.cMacroDecl(Cgen.cMacroMin))
            incblock.add(Cgen.cMacroDecl(Cgen.cMacroMax))
            incblock.add(Cgen.cMacroDecl(Cgen.cMacroFloord))
        with m.funcs as funcblock:
            # Collect all the Inputs and Parameters of the pipeline. 
            # Add them as parameters to the pipeline function.
            # TODO:
            # Also set these data as pipeline properties, since they are
            # required again later during autotuning.
            params = []
            for s in stageOrder:
                for obj in s.computeObjs:
                    params = params + obj.getObjects(Parameter)

            pipelineParams = OrderedDict()
            params = list(set(params))

            cparamMap = {}
            cfuncMap  = {}

            # 1. Collecting scalar parameters
            params.sort(key=lambda x: x.name)
            for param in params:
                cvar = Cgen.cVariable(Cgen.cgenType.get(param.typ), param.name)
                if param.definition is None:
                    pipelineParams[cvar] = cvar.typ
                else:
                    assert False
                # Binding parameters to C variables
                cparamMap[param] = cvar

            # 2. Collecting image parameters
            self.inputs.sort(key=lambda x: x.name)
            for img in self.inputs:
                cpoint = Cgen.cPointer(Cgen.cgenType.get(img.typ), 1)
                cvar = Cgen.cVariable(cpoint, img.name)
                pipelineParams[cvar] = cvar.typ
                # Binding array parameters to C array
                carray = Cgen.cArray(Cgen.cgenType.get(img.typ),
                                     img.name, img.dimensions)
                cfuncMap[img] = (carray,
                                 [ False for i in xrange(0, len(img.dimensions))])
                carray.layout = 'contigous'

            # 3. Collecting output parameters
            self.outputs.sort(key=lambda x: x.name)

            # areParamsVoidPtrs : if the target is to generate shared library
            # implementation using python ctypes

            # isExternAlloc : is the result array allocated outside the
            # generated implementation

            if not isExternAlloc:
                for func in self.outputs:
                    cref = Cgen.cReference(Cgen.cgenType.get(func.typ), 1)
                    cvar = Cgen.cVariable(cref, func.name)
                    pipelineParams[cvar] = cvar.typ
            else:
                for func in self.outputs:
                    cpoint = Cgen.cPointer(Cgen.cgenType.get(func.typ), 1)
                    cvar = Cgen.cVariable(cpoint, func.name)
                    pipelineParams[cvar] = cvar.typ

            pipeline = Cgen.cFunction(Cgen.cVoid, 'pipeline_' + self.name, pipelineParams)
            pipelineDecl = Cgen.cFunctionDecl(pipeline, isExternFunc, areParamsVoidPtrs)
            pipelineBody = Cgen.cFunctionBody(pipelineDecl)

            funcblock.add(pipelineBody)

            with pipelineBody.body as pbody:
                if areParamsVoidPtrs:
                    # Add assignment to inputs
                    for inp in self.inputs:
                        # actual input to be used
                        varType = Cgen.cgenType.get(inp.typ)
                        varPtr = Cgen.cPointer(varType, 1)
                        var = Cgen.cVariable(varType, inp.name)
                        varDecl = Cgen.cDeclaration(varPtr, var)
                        pbody.add(varDecl)

                        # dummy void * input taken as argument
                        dummyType = Cgen.cgenType.get(inp.typ)
                        dummyPtr = Cgen.cPointer(dummyType, 1)
                        dummyCast = Cgen.cCast(dummyPtr, inp.name+'_void_arg')
                        varAssign = Cgen.cAssign(var, dummyCast)
                        pbody.add(varAssign)

                    # Add assignment to output
                    for out in self.outputs:
                        # actual output to be used
                        varType = Cgen.cgenType.get(out.typ)
                        varPtr = Cgen.cPointer(varType, 1)
                        var = Cgen.cVariable(varType, out.name)
                        varDecl = Cgen.cDeclaration(varPtr, var)
                        pbody.add(varDecl)

                        # dummy void * output taken as argument
                        dummyType = Cgen.cgenType.get(out.typ)
                        dummyPtr = Cgen.cPointer(dummyType, 1)
                        dummyCast = Cgen.cCast(dummyPtr, out.name+'_void_arg')
                        varAssign = Cgen.cAssign(var, dummyCast)
                        pbody.add(varAssign)

                for s in stageOrder:
                    s.generateCNaive(pbody, cfuncMap, cparamMap, self.outputs, isExternAlloc)

        return m

    def __str__(self):
        return_str = "Final Stage: " + self._name + "\n"
        for s in self._stages:
            return_str = return_str + s.__str__() + "\n"
        return return_str

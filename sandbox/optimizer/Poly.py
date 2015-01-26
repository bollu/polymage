import islpy as isl
from constructs import *
import math as math

def lcm(a, b):
    return a*b/(gcd(a, b))

def optimizeSchedule(partScheds, dependencies):
    # The pluto optimizer can be used to optimize the schedule for comparision.
    pass

def stripMineSchedule(sched, dim, size):
    sched = sched.insert_dims(isl._isl.dim_type.out, dim, 1)
    name = sched.get_dim_name(isl._isl.dim_type.out, 1 + dim) 
    sched = sched.set_dim_name(isl._isl.dim_type.out, dim, 'S_' + name)
    ineqs = []
    #  size*(Ti) <= i <= size*(Ti) + size - 1
    coeff = {}
    coeff[('out', dim)] = sizes[dim - startDim] 
    coeff[('constant', 0)] = sizes[dim - startDim] - 1
    coeff[('out', numDims + dim)] = -1
    ineqs.append(coeff)

    coeff = {}
    coeff[('out', dim)] = -sizes[dim - startDim] 
    coeff['out', numDims + dim] = 1
    ineqs.append(coeff) 
    sched = addConstriants(sched, ineqs, [])

    return sched

def tileSchedule(sched, dim, size, overlapOffset = 0):
    # Extend space to accomodate the tiling dimensions
    sched = sched.insert_dims(isl._isl.dim_type.out, dim, 1)
    # Create the tile dimensions and their constraints
    name = sched.get_dim_name(isl._isl.dim_type.out, 1 + dim) 
    sched = sched.set_dim_name(isl._isl.dim_type.out, dim, '_T' + name)

    ineqs = []
    #  size*(Ti) <= i <= size*(Ti) + size - 1
    coeff = {}
    coeff[('out', dim)] = size 
    coeff[('constant', 0)] = size - 1 + overlapOffset
    coeff[('out', 1 + dim)] = -1
    ineqs.append(coeff)

    coeff = {}
    coeff[('out', dim)] = -size 
    coeff['out', 1 + dim] = 1
    ineqs.append(coeff) 
    sched = addConstriants(sched, ineqs, [])
    return (sched, ('rect', name, '_T' + name, size))

def addConstraintsFromList(obj, localSpace, constraintList, constraintAlloc):
    for const in constraintList:
        c = constraintAlloc(localSpace)
        m = 1
        for coeff in const:
            if isinstance(const[coeff], Fraction):
                m = (abs(const[coeff].denominator) * m)/gcd(abs(const[coeff].denominator), m)
        assert m.denominator == 1
        m = m.numerator
        for coeff in const:
            if isinstance(const[coeff], Fraction):
               const[coeff] = m*const[coeff]
               assert const[coeff].denominator == 1
               const[coeff] = const[coeff].numerator
            else:
               const[coeff] = m * const[coeff]
        for coeff in const:
            dim = coeff[1]
            if coeff[0] == 'param':                
                if (type(dim) == str):
                    dim = obj.find_dim_by_name(isl._isl.dim_type.param, dim)
                c = c.set_coefficient_val(isl._isl.dim_type.param, dim, const[coeff])
            elif coeff[0] == 'in':
                if (type(dim) == str):
                    dim = obj.find_dim_by_name(isl._isl.dim_type.in_, dim)
                c = c.set_coefficient_val(isl._isl.dim_type.in_, dim, const[coeff])
            elif coeff[0] == 'out':
                if (type(dim) == str):
                    dim = obj.find_dim_by_name(isl._isl.dim_type.out, dim)
                c = c.set_coefficient_val(isl._isl.dim_type.out, dim, const[coeff])
            elif coeff[0] == 'constant':
                c = c.set_constant_val(const[coeff])
            else:
               assert False
        obj = obj.add_constraint(c)
    return obj

def addConstriants(obj, ineqs, eqs):
    space = obj.get_space()
    if (isinstance(obj, isl.Map)):
        for bmap in obj.get_basic_maps():
            localSpace = bmap.get_local_space()
   
            obj = addConstraintsFromList(obj, localSpace, ineqs, 
                                 isl.Constraint.inequality_alloc) 
            obj = addConstraintsFromList(obj, localSpace, eqs, 
                                 isl.Constraint.equality_alloc) 
    elif (isinstance(obj, isl.Set)):
        for bset in obj.get_basic_sets():
            localSpace = bset.get_local_space()
   
            obj = addConstraintsFromList(obj, localSpace, ineqs, 
                                 isl.Constraint.inequality_alloc) 
            obj = addConstraintsFromList(obj, localSpace, eqs, 
                                 isl.Constraint.equality_alloc)
    elif (isinstance(obj, isl.BasicSet) or 
          isinstance(obj, isl.BasicMap)):
        localSpace = obj.get_local_space()
        obj = addConstraintsFromList(obj, localSpace, ineqs, 
                            isl.Constraint.inequality_alloc) 
        obj = addConstraintsFromList(obj, localSpace, eqs, 
                            isl.Constraint.equality_alloc)
    else:
        assert False
        
    return obj

def extractValueDependence(part, ref, refPolyDom):
    # Dependencies are calculated between values. There is no storage
    # mapping done yet.
    deps = []
    accessRegion = isl.BasicSet.universe(refPolyDom.domSet.get_space())
    partDom = part.sched.domain().align_params(refPolyDom.domSet.get_space())    
    accessRegion = accessRegion.align_params(partDom.get_space())
    
    rel = isl.BasicMap.from_domain_and_range(partDom, accessRegion)
    dimOut = rel.dim(isl._isl.dim_type.out)
    sourceDims = [ ('out', i) for i in xrange(0, dimOut)]
    numArgs = len(ref.arguments)                

    for i in xrange(0, numArgs):
        arg = ref.arguments[i]
        # If the argument is not affine the dependence reflects that
        # the computation may depend on any value of the referenced object
        if (isAffine(arg)):
            coeff = getAffineVarAndParamCoeff(arg)
            coeff = mapCoeffToDim(coeff)

            coeff[('constant', 0)] = getConstantFromExpr(arg, affine = True)
            coeff[sourceDims[i]] = -1
            rel = addConstriants(rel, [], [coeff])
    if not rel.is_empty():
        deps.append(PolyDep(ref.objectRef, part.comp, rel))
    return deps 

class PolyPart(object):
    def __init__(self, _sched, _expr, _pred, _comp, _align, 
                 _scale, _stageNo, _liveout):
        self.sched = _sched
        self.expr = _expr
        self.pred = _pred
        self.comp = _comp
        # Dependencies between values of computation objects
        self.deps = []
        # Should depvecs even be a member of polypart?
        # The should be associated with the group?
        self.depVecs = []
        # Mapping between the input variables to the corresponding 
        # schedule dimension. A full affine schedule will need a 
        # transformation matrix. Currently we only shuffle the 
        # dimension order apart from tiling so a simple dimension
        # alignment vector suffices. This has to be changed to 
        # handle more general cases later.
        self.align = _align
        # Scaling factors for each schedule dimension
        self.scale = _scale
        # Default alignment and scaling factors are set while
        # constructing the polypart. These are changed by the
        # alignment and loop scaling passes. Both these passes 
        # attempt to improve locality and uniformize dependencies.
        self.stageNo = _stageNo
        self.dimTileInfo = {}
        self.dimScratchSize = {}
        self.parallelSchedDims = []
        self.vectorSchedDim = []
        self.liveout = _liveout
    
    def getPartRefs(self):
        refs = self.expr.collect(Reference)
        if (self.pred):
            refs += self.pred.collect(Reference)
        return refs

    def __str__(self):
        partStr = "Schedule: " + self.sched.__str__() + '\n'\
                  "Expression: " + self.expr.__str__() + '\n'\
                  "Predicate: " + self.pred.__str__() + '\n'
        depstr = 'Dependencies:' + '\n'
        for dep in self.deps:
            depstr = depstr + dep.__str__() + '\n'
        return partStr + depstr

class PolyDomain(object):
    def __init__(self, _domSet, _comp):
        self.domSet = _domSet
        self.comp = _comp
        
    def __str__(self):
        return "Domain: " + self.domSet.__str__()

class PolyDep(object):
    def __init__(self, _producerObj, _consumerObj, _rel):
        self.producerObj = _producerObj
        self.consumerObj = _consumerObj
        self.rel         = _rel
    
    def __str__(self):
        return self.rel.__str__()

class PolyRep(object):
    """ The PolyRep class is the polyhedral representation of a 
        stage. It gives piece-wise domain and schedule for a stage.
        Polyhedral transformations modify the piece-wise domains as 
        well as the schedules.
    """
    def __init__(self, _ctx, _stage, _paramConstraints, _paramEstimates, 
                 _tileSizes, _sizeThreshold, _groupSize, _outputs):
        self.stage = _stage
        self.ctx = _ctx
        self.polyParts = {}
        self.polyDoms = {}
        self.polyast = []        
        self.groups = []
        self.outputs = _outputs
        self.tileSizes = _tileSizes
        self.sizeThreshold = _sizeThreshold
        self.groupSize = _groupSize

        self._varCount = 0
        self._funcCount = 0

        if(_stage.isPolyhedral()):
            # Currently doing extraction only when all the computeObjs
            # domains are affine. This can be revisited later. 
            self.extractPolyRepFromStage(_paramConstraints)
            self.fusedSchedule(_paramEstimates)
            #self.simpleSchedule(_paramEstimates)
    
    def generateCode(self):
        self.polyast = []
        if self.polyParts:
            self.buildAst()

    def __str__(self):
        polystr = ""
        for comp in self.polyParts:
            for part in self.polyParts[comp]:
                polystr = polystr + part.__str__() + '\n'

        if (self.polyast != []):
            for ast in self.polyast:
                printer = isl.Printer.to_str(self.ctx)
                printer = printer.set_output_format(isl.format.C)
                printOpts = isl.AstPrintOptions.alloc(self.ctx) 
                printer = ast.print_(printer, printOpts)
                aststr = printer.get_str()
                polystr = polystr + '\n' + aststr
        return polystr
   
    def formatParamConstraints(self, paramConstraints, params):
        contextConds = []
        for paramConst in paramConstraints:
            paramsInConst = paramConst.collect(Parameter)
            contextAdd = True
            for p in paramsInConst:
                if p not in params:
                    contextAdd = False
            if contextAdd and isAffine(paramConst):
                paramConstConjunct = paramConst.splitToConjuncts()
                if len(paramConstConjunct) == 1:
                    contextConds.append(paramConst)
        return contextConds             

    def extractPolyRepFromStage(self, paramConstraints):
        compObjs = self.stage.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: s[1])
        stageDomMin = sortedCompObjs[0][1]
        numObjs = len(sortedCompObjs)
        stageDomMax = sortedCompObjs[numObjs -1][1] + 1

        # Comute the max dimensionality of the space and collect the parameters.
        def maxDim(objs):
            dim = 0
            for comp in objs:
                if type(comp) == Accumulator:
                    dim = max(dim, len(comp.reductionVariables))
                    dim = max(dim, len(comp.variables))
                elif type(comp) == Function or type(comp) == Image:
                    dim = max(dim, len(comp.variables))
            return dim

        params = []
        dim = maxDim(compObjs)
        for comp in compObjs:
            params = params + comp.getObjects(Parameter)
        params = list(set(params))
        paramNames = [param.name for param in params]        

        contextConds = self.formatParamConstraints(paramConstraints, params)
                    
        # Have to revisit this hopefully there is a better way to do this
        # than just concatenating an '
        # The [t] is for the stage dimension
        scheduleNames = ['_t'] + [ self.getVarName()  for i in xrange(0, dim) ]

        for comp in compObjs:
            if (type(comp) == Function or type(comp) == Image):
                self.extractPolyRepFromFunction(comp, scheduleNames, paramNames,
                                                contextConds, compObjs[comp] + 1,
                                                paramConstraints)
            elif (type(comp) == Accumulator):
                self.extractPolyRepFromAccumulator(comp, scheduleNames, paramNames,
                                                   contextConds, compObjs[comp] + 1,
                                                   paramConstraints)
            else:
                assert False
                    
        # Check if the basic sets in the domain are non-empty and non-overlapping
        # This might not always be feasible since we allow for non-affine 
        # approximations. Can runtime checks be generated?

    def createSchedSpace(self, variables, domains, scheduleNames, paramNames,
                         contextConds):
        # Variable names for refrerring to dimensions
        varNames = [ variables[i].name for i in xrange(0, len(variables)) ]        
        space = isl.Space.create_from_names(self.ctx, in_ = varNames,
                                                     out = scheduleNames,
                                                     params = paramNames)

        schedMap = isl.BasicMap.universe(space)
        # Adding the domain constraints
        [ineqs, eqs] = formatDomainConstraints(domains, varNames)
        schedMap = addConstriants(schedMap, ineqs, eqs)

        [paramIneqs, paramEqs] = formatConjunctConstraints(contextConds)
        schedMap = addConstriants(schedMap, paramIneqs, paramEqs)
        return schedMap
   
    def makePolyParts(self, sched, expr, pred, comp, align, scale, 
                      stageNo, liveout):
        # Detect selects with modulo constraints and split into 
        # multiple parts. This technique can also be applied to the
        # predicate but for now we focus on selects.
        polyParts = []
        # This is very very temporary solution there should be a 
        # better way of doing this. Only targetting conditions 
        # of the form (affine)%constant == constant.
        brokenParts = []
        if isinstance(expr, Select):
            conjuncts = expr.condition.splitToConjuncts()
            if len(conjuncts) == 1 and len(conjuncts[0]) == 1:
                cond = conjuncts[0][0]
                leftExpr = cond.lhs
                rightExpr = cond.rhs
                isLeftModulo = isAffine(leftExpr, modulo = True) and \
                               not isAffine(leftExpr)
                isRightConstant = isConstantExpr(rightExpr)
                breakSelect = False
                if isLeftModulo and isRightConstant and \
                    cond.conditional == '==' and \
                    isinstance(leftExpr, AbstractBinaryOpNode)\
                    and leftExpr.op == '%' and isAffine(leftExpr.left)\
                    and isConstantExpr(leftExpr.right):
                    breakSelect = True
                if breakSelect:
                    leftCoeff = getAffineVarAndParamCoeff(leftExpr.left)
                    leftCoeff = mapCoeffToDim(leftCoeff)
                    leftConst = getConstantFromExpr(leftExpr.left, affine = True)
                    rightConst = getConstantFromExpr(rightExpr, affine = True)
                    modConst = getConstantFromExpr(leftExpr.right, affine = True)
                
                    mulName = '_Mul_'
                    remName = '_Rem_'
                    trueSched = sched.copy()
                    dimIn = trueSched.dim(isl._isl.dim_type.in_)
                    trueSched = trueSched.insert_dims(isl._isl.dim_type.in_, 
                                                      dimIn, 1)
                    trueSched = trueSched.set_dim_name(isl._isl.dim_type.in_, 
                                                       dimIn, mulName)
                
                    eqs = []
                    leftCoeff[('constant', 0)] = leftConst - rightConst 
                    leftCoeff[('in', dimIn)] = -modConst
                    eqs.append(leftCoeff)

                    trueSched = addConstriants(trueSched, [], eqs)
                    trueSched = trueSched.project_out(isl._isl.dim_type.in_, 
                                                      dimIn, 1)
                    brokenParts.append((trueSched, expr.trueExpression))

                    falseSched = sched.copy()
                    dimIn = falseSched.dim(isl._isl.dim_type.in_)
                    falseSched = falseSched.insert_dims(isl._isl.dim_type.in_, 
                                                        dimIn, 2)
                    falseSched = falseSched.set_dim_name(isl._isl.dim_type.in_,
                                                         dimIn, mulName)
                    falseSched = falseSched.set_dim_name(isl._isl.dim_type.in_, 
                                                         dimIn + 1, remName)
                
                    eqs = []
                    leftCoeff[('constant', 0)] = leftConst - rightConst 
                    leftCoeff[('in', dimIn)] = -modConst
                    leftCoeff[('in', dimIn+1)] = -1
                    eqs.append(leftCoeff)

                    ineqs = []
                    coeff = {}
                    coeff[('in', dimIn+1)] = 1
                    coeff[('constant', 0)] = -1
                    ineqs.append(coeff)
                
                    coeff = {}
                    coeff[('in', dimIn+1)] = -1
                    coeff[('constant', 0)] = modConst - 1
                    ineqs.append(coeff)

                    falseSched = addConstriants(falseSched, ineqs, eqs)
                    falseSched = falseSched.project_out(isl._isl.dim_type.in_,
                                                        dimIn, 2)                    
                    brokenParts.append((falseSched, expr.falseExpression))

                    #print trueSched
                    #print falseSched
        if not brokenParts:
            polyPart = PolyPart(sched, expr, pred, comp, list(align), 
                                list(scale), stageNo, liveout)
            id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)                            
            polyPart.sched = polyPart.sched.set_tuple_id(
                                          isl._isl.dim_type.in_, id_)
            polyParts.append(polyPart)
        else:
            for bsched, bexpr in brokenParts:
                polyPart = PolyPart(bsched, bexpr, pred, comp, list(align), 
                                    list(scale), stageNo, liveout)
                id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)                            
                polyPart.sched = polyPart.sched.set_tuple_id(
                                          isl._isl.dim_type.in_, id_)
                polyParts.append(polyPart)
        return polyParts       

    def createPolyPartsFromDefinition(self, comp, schedMap, stageNo, 
                                      scheduleNames, domain):
        self.polyParts[comp] = []
        for case in comp.definition:
            sched = schedMap.copy()
            liveout = comp in self.outputs

            # The basic schedule is an identity schedule appended with 
            # a stage dimension. The stage dimension gives the ordering 
            # of the stages.

            align, scale = self.defaultAlignAndScale(sched, domain)
            
            if (isinstance(case, Case)):
                # Dealing with != and ||. != can be replaced with < || >. 
                # and || splits the domain into two.
                splitConjuncts = case.condition.splitToConjuncts()
                for conjunct in splitConjuncts:
                    # If the condition is non-affine it is stored as a 
                    # predicate for the expression. An affine condition 
                    # is added to the domain.
                    affine = True
                    for cond in conjunct:
                        affine = affine and isAffine(cond.lhs) and isAffine(cond.rhs)
                    if(affine):
                        [conjunctIneqs, conjunctEqs] = \
                                formatConjunctConstraints(conjunct)
                        sched = addConstriants(sched, conjunctIneqs, conjunctEqs)
                        # Create a user pointer, tuple name and add it to the map
                        parts = self.makePolyParts(sched, case.expression, None,
                                                   comp, align, scale, stageNo, 
                                                   liveout)
                        for part in parts:
                            self.polyParts[comp].append(part)
                    else:
                        parts = self.makePolyParts(sched, case.expression, 
                                                   case.condition, comp, align, 
                                                   scale, stageNo, liveout)

                        for part in parts:
                            self.polyParts[comp].append(part)
            else:
                assert(isinstance(case, AbstractExpression) or 
                        isinstance(case, Accumulate))
                parts = self.makePolyParts(sched, case, None, comp, 
                                           align, scale, stageNo, liveout)
                for part in parts:
                    self.polyParts[comp].append(part)
        # Subtract all the part domains to find the domain where the
        # default expression has to be applied
        #sched = isl.BasicMap.identity(self.polyspace)
        #sched = addConstriants(sched, ineqs, eqs)
        # Adding stage identity constraint
        #stageCoeff = {}
        #stageCoeff[varDims[0]] = -1
        #stageCoeff[('constant', 0)] = compObjs[comp]
        #sched = addConstriants(sched, [], [stageCoeff])
        #sched = addConstriants(sched, paramIneqs, paramEqs)

        #for part in self.polyParts[comp]:
        #    sched = sched.subtract_range(part.sched.range())
        #    if (sched.is_empty()):
        #        break
        #if(not sched.fast_is_empty()):
        #    bmapList = []
        #    if (isinstance(sched, isl.BasicMap)):
        #        bmapList.append(sched)
        #    else:
        #        sched.foreach_basic_map(bmapList.append)
        #    for bmap in bmapList:    
        #        polyPart = PolyPart(bmap, comp.default, None, comp)
        #        id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)                            
        #        polyPart.sched = polyPart.sched.set_tuple_id(
        #                                   isl._isl.dim_type.in_, id_)
        #        self.polyParts[comp].append(polyPart)
        
    def extractPolyRepFromFunction(self, comp, scheduleNames, paramNames,
                                   contextConds, stageNo, paramConstraints):
        self.polyDoms[comp] = self.extractPolyDomFromComp(comp, paramConstraints)
        schedMap = self.createSchedSpace(comp.variables, comp.domain, scheduleNames,
                                         paramNames, contextConds)
        self.createPolyPartsFromDefinition(comp, schedMap, stageNo, scheduleNames,
                                           comp.domain)

    def extractPolyRepFromAccumulator(self, comp, scheduleNames, paramNames,
                                      contextConds, stageNo, paramConstraints):
        self.polyDoms[comp] = self.extractPolyDomFromComp(comp, paramConstraints)
        schedMap = self.createSchedSpace(comp.reductionVariables, comp.reductionDomain, 
                                         scheduleNames, paramNames, contextConds)
        self.createPolyPartsFromDefinition(comp, schedMap, stageNo, scheduleNames,
                                           comp.reductionDomain)
        domMap = self.createSchedSpace(comp.variables, comp.domain, scheduleNames, 
                                       paramNames, contextConds)
        self.createPolyPartsFromDefault(comp, domMap, stageNo, scheduleNames)

    def defaultAlignAndScale(self, sched, domain):
        dimOut = sched.dim(isl._isl.dim_type.out)
        dimIn = sched.dim(isl._isl.dim_type.in_)
        align = [ i+1 for i in xrange(0, dimIn) ]
        scale = []
        for i in xrange(0, dimIn):
            scale.append(domain[i].step.value)
        return (align, scale)

    def createPolyPartsFromDefault(self, comp, schedMap, stageNo, scheduleNames):
        sched = schedMap.copy()
        align, scale = self.defaultAlignAndScale(sched, comp.domain)

        assert(isinstance(comp.default, AbstractExpression))
        polyPart = PolyPart(sched, comp.default, None, comp,
                            align, scale, stageNo - 1, False)
        id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)
        polyPart.sched = polyPart.sched.set_tuple_id(isl._isl.dim_type.in_, id_)
        self.polyParts[comp].append(polyPart)

    def extractPolyDomFromComp(self, comp, paramConstraints):
        varNames = [ var.name for var in comp.variables ]
        domMapNames = [ name +'\'' for name in varNames]
        params = []
        for interval in comp.domain:
            params = params + interval.collect(Parameter)
        params = list(set(params))
        paramNames = [ param.name for param in params ]

        space = isl.Space.create_from_names(self.ctx, in_ = varNames,
                                                      out = domMapNames,
                                                      params = paramNames)
        domMap = isl.BasicMap.universe(space)
        # Adding the domain constraints
        [ineqs, eqs] = formatDomainConstraints(comp.domain, varNames)
        domMap = addConstriants(domMap, ineqs, eqs)
       
        paramConds = self.formatParamConstraints(paramConstraints, params) 
        [paramIneqs, paramEqs] = formatConjunctConstraints(paramConds)
        domMap = addConstriants(domMap, paramIneqs, paramEqs)

        polyDom = PolyDomain(domMap.domain(), comp)
        id_ = isl.Id.alloc(self.ctx, comp.name, polyDom)
        polyDom.domSet = polyDom.domSet.set_tuple_id(id_)
        return polyDom

    def computeDependencies(self):
        # Extract dependencies. In case the dependencies cannot be exactly 
        # represented they are approximated.
        for comp in self.polyParts:
            for part in self.polyParts[comp]:
                refs = part.getPartRefs()
                # There could be multiple references to the same value. So the 
                # dependencies might be redundant this has to be eliminated 
                # later.
                for ref in refs:
                    # Considering only dependencies within the same stage
                    if ref.objectRef in self.polyParts:
                        part.deps += extractValueDependence(part, ref, 
                                                self.polyDoms[ref.objectRef])

    def getDependenceVector(self, parentPart, childPart, ref, scaleMap = None):
        def getScale(sMap, p, i):
            if sMap is not None:
                return sMap[p][i]
            return p.scale[i]

        numArgs = len(ref.arguments)
        dimOut = parentPart.sched.dim(isl._isl.dim_type.out)
        depVec = [ '-' for i in xrange(0, dimOut) ]

        if isinstance(parentPart.comp, Accumulator):
            for i in xrange(1, dimOut):
                depVec[i] = '*'
            depVec[0] = childPart.stageNo - parentPart.stageNo
            return (depVec, parentPart.stageNo)

        for i in xrange(0, numArgs):
            arg = ref.arguments[i]
            parentVarSchedDim = parentPart.align[i]            
            if (isAffine(arg)):
                domDimCoeff = self.getDomainDimCoeffs(childPart.sched, arg)
                paramCoeff = self.getParamCoeffs(childPart.sched, arg)
                # Parameter coefficents can also be considered to
                # generate parametric shifts. Yet to be seen.
                
                # Indexed with multiple variables.
                if (len(domDimCoeff) > 1 or 
                    (len(domDimCoeff) == 1 and len(paramCoeff) >=1)):
                    # Although there are multiple coefficients. If 
                    # there is only one variable coefficient and other
                    # parametric coefficients. Uniformization can be 
                    # done with parametric shifts. Full affine scheduling 
                    # might be able to find a way to uniformize 
                    # dependencies. This has to be further explored.
                    #print ref
                    #assert False
                    depVec[parentVarSchedDim] = '*'
                # Indexed with a single variable. This can either be an 
                # uniform access or can be uniformized with scaling when 
                # possible.
                elif len(domDimCoeff) == 1 and len(paramCoeff) == 0:                    
                    dim = (domDimCoeff.keys())[0]
                    childVarSchedDim = childPart.align[dim]

                    pscale = getScale(scaleMap, parentPart, i)
                    cscale = getScale(scaleMap, childPart, dim)

                    assert Fraction(pscale).denominator == 1
                    assert Fraction(cscale).denominator == 1

                    if ((childVarSchedDim == parentVarSchedDim) and 
                        (domDimCoeff[dim] * pscale == cscale)):
                        depVec[parentVarSchedDim] = \
                                -getConstantFromExpr(arg, affine=True)
                        accessScale = pscale
                        if depVec[parentVarSchedDim] > 0:
                            depVec[parentVarSchedDim] = \
                                (int(math.ceil(depVec[parentVarSchedDim] *
                                accessScale)))
                        else:         
                            depVec[parentVarSchedDim] = \
                               (int(math.floor(depVec[parentVarSchedDim] *
                                accessScale)))
                            #print parentPart.sched
                            #print childPart.sched
                            #print childPart.expr
                            #print depVec, ref
                    else:
                        depVec[parentVarSchedDim] = '*'
                elif len(domDimCoeff) == 0 and len(paramCoeff) > 0:
                    #print ref
                    #assert False
                    depVec[parentVarSchedDim] = '*'
                # Only parametric or Constant access. The schedule in this 
                # dimension can be shifted to this point to uniformize the 
                # dependence.
                # In case the dimension in the parent has a constant size
                # an upper and lower bound on the dependence vector can 
                # be computed.
                elif len(domDimCoeff) + len(paramCoeff) == 0:
                    # offsets should be set here.
                    accessConstant = getConstantFromExpr(arg, affine = True)
                    parentLowerBound = parentPart.sched.domain().dim_min(i)
                    parentUpperBound = parentPart.sched.domain().dim_max(i)
                    if ((parentLowerBound.is_cst() and 
                        parentLowerBound.n_piece() == 1) and
                        (parentUpperBound.is_cst() and
                        parentUpperBound.n_piece() == 1)):
                          
                        pscale = getScale(scaleMap, parentPart, i)

                        lowVecAff = (parentLowerBound.get_pieces())[0][1]
                        val = lowVecAff.get_constant_val()
                        assert(val.get_den_si() == 1)
                        lowVec = int(math.floor((accessConstant - val.get_num_si()) *
                                             pscale))

                        highVecAff = (parentUpperBound.get_pieces())[0][1]
                        val = highVecAff.get_constant_val()
                        assert(val.get_den_si() == 1)
                        highVec = int(math.ceil((accessConstant - val.get_num_si()) *
                                             pscale))

                        if highVec == lowVec:
                            depVec[parentVarSchedDim] = highVec
                        else:
                            # Unpack dependence vectors when this hits
                            #print ref
                            #assert False
                            #depVec[parentVarSchedDim] = (lowVec, highVec)
                            depVec[parentVarSchedDim] = '*'
                    else:
                        depVec[parentVarSchedDim] = '*'
                else:
                    assert False
            else:
                #print ref
                #assert False
                depVec[parentVarSchedDim] = '*'

        assert depVec[0] == '-'
        depVec[0] = childPart.stageNo - parentPart.stageNo
        for i in xrange(0, dimOut):
            if (depVec[i] == '-'):
                depVec[i] = 0
        #for i in xrange(0, dimOut):
        #    if (depVec[i] == '-'):
        #        depVec[i] = '*'
        #        print parentPart.sched
        #        print childPart.sched
        #        print i, ref
        #        parentLowerBound = parentPart.sched.range().dim_min(i)
        #        parentUpperBound = parentPart.sched.range().dim_max(i)
                    
        #        childLowerBound  = childPart.sched.range().dim_min(i)
        #        childUpperBound = childPart.sched.range().dim_max(i)

        #        if (childLowerBound.is_equal(childUpperBound) and
        #            parentLowerBound.is_equal(parentUpperBound)):
        #            dimDiff = childUpperBound.sub(parentUpperBound)
        #            if (dimDiff.is_cst() and dimDiff.n_piece() == 1):
        #                aff = (dimDiff.get_pieces())[0][1]
        #                val = aff.get_constant_val()
        #                depVec[i] = (val.get_num_si())/(val.get_den_si())
        return (depVec, parentPart.stageNo)

    def getGroupDependenceVectors(self, group, scaleMap = None):
        depVecs = []   
        for part in group:
            refs = part.getPartRefs()
            for ref in refs:
                if ref.objectRef in self.polyParts:
                    for pp in self.polyParts[ref.objectRef]:
                        if pp not in group:
                            continue
                        depVec = self.getDependenceVector(pp, part, ref, 
                                                          scaleMap)
                        depVecs.append(depVec)  
        return depVecs

    def alignWithGroup(self, part, group):

        def getDomainDimsInvolved(sched, arg):
            domDims = []
            if (isAffine(arg)):
                coeff = getAffineVarAndParamCoeff(arg)
                for item in coeff:
                    if type(item) == Variable:
                        dim = sched.find_dim_by_name(isl._isl.dim_type.in_,
                                                     item.name)
                        domDims.append(dim)
            return domDims

        dimIn      = part.sched.dim(isl._isl.dim_type.in_)
        refs       = part.getPartRefs()
        # Avoiding Input references this should be revisited at some point
        parentRefs = [ ref for ref in refs \
                       if ref.objectRef in self.polyParts ]
        dimAlignMap = {}           
        for ref in parentRefs:
            parentParts = self.polyParts[ref.objectRef]
            # Filter out self references
            groupParentParts = [ pp for pp in parentParts \
                                 if pp in group and pp != part]
            if not groupParentParts:
                continue
            numArgs = len(ref.arguments)
            for i in xrange(0, numArgs):
                arg = ref.arguments[i]
                domDims = getDomainDimsInvolved(part.sched, arg)
                # This can get tricky with multiple parts have to revisit
                for repParentPart in groupParentParts:
                    for dim in domDims:
                        if dim not in dimAlignMap:
                            dimAlignMap[dim] = [repParentPart.align[i]]
                        if repParentPart.align[i] not in dimAlignMap[dim]:
                            dimAlignMap[dim].append(repParentPart.align[i])

        newAlign = [ '-' for i in xrange(0, part.sched.dim(isl._isl.dim_type.in_))]
        for i in xrange(0, dimIn):
            alignDim = None
            # Check if the dimension is uniquely mapped
            if i in dimAlignMap and len(dimAlignMap[i]) == 1:
                alignDim = dimAlignMap[i][0]
            # Remove the assigned dimension from the maps    
            for dim in dimAlignMap:
                if alignDim in dimAlignMap[dim]:
                    dimAlignMap[dim].remove(alignDim)
            if alignDim is not None:
                newAlign[i] = alignDim
        return newAlign 

    def alignParts(self):
        """ Embeds parts whose dimension is smaller than the schedule space."""
        # Go through the parts in a sorted order and compute alignments
        compObjs = self.stage.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        # Alignments in certain cases me result in chaning the relative order of
        # dimensions. This is only valid if there are no dependencies between the
        # dimensions being reordered. Dependence vectors or polyhedra can be used
        # to decide if the dimensions can be reordered.

        # For now we only realign parts which do not have self dependencies
        # Find all parts which do not have self dependencies
        noDepParts = []
        for comp in [item[0] for item in sortedCompObjs]:
            for p in self.polyParts[comp]:
                if not self.isPartSelfDependent(p):
                    noDepParts.append(p)

        # Find an alignment map for parts which are not full dimensional
        def completeMap(newAlign, currAlign):
            for i in xrange(0, len(newAlign)):
                if newAlign[i] in currAlign:
                    currAlign.remove(newAlign[i])
            for i in xrange(0, len(newAlign)):         
                if newAlign[i] == '-':
                    newAlign[i] = currAlign.pop(0)
            return newAlign

        def compatibleAlign(align1, align2):
            compatible = True
            if len(align1) == len(align2):
                for i in xrange(0, len(align1)):
                    if not ((align1[i] == '-' or align2[i] == '-')
                            or (align1[i] == align2[i])):
                        compatible = False
            else:
                compatible = False
            return compatible

        def alignGroupToPart(group, part, align):
            for g in group:
                dimIn = g.sched.dim(isl._isl.dim_type.in_)
                dimOut = g.sched.dim(isl._isl.dim_type.out)
                newAlign = [ '-' for i in xrange(0, dimIn)]
                for i in xrange(0, dimIn):
                    for j in xrange(0, len(align)):
                        if g.align[i] == align[j]:
                            assert newAlign[i] == '-'
                            newAlign[i] = part.align[j]
                g.align = completeMap(newAlign, g.align)

        # Progressive alignment algorithm. Need to revisit when we encounter
        # stranger things. alignedPartGroups consists of groups of parts whose
        # alignments are linked. All the alignments in a group should be reordered 
        # in the same way so as to not disturb the previous alignments.
        alignedPartGroups = []
        for p in noDepParts:
            # Find if there is a unique alignment with the aligned parts
            # i.e the current part can be aligned with already aligned 
            # parts or all the aligned parts can be aligned to the current 
            # part.
            parentGroups = self.findParentGroups(p, alignedPartGroups)
            otherGroups = [ g for g in alignedPartGroups \
                            if g not in parentGroups ]
            aligns = {}

            for i in xrange(0, len(parentGroups)):
                aligns[i] = self.alignWithGroup(p, parentGroups[i])
            
            mergeGroups = []
            # If the alignment has alteast one valid reordering. Add the
            # group to the list of groups to be aligned and merged.
            for i in xrange(0, len(parentGroups)):
                addGroup = False
                for dim in aligns[i]:
                    if dim != '-':
                        addGroup = True
                if addGroup:
                    mergeGroups.append(i)

            mergedGroup = []
            for i in mergeGroups:
                alignGroupToPart(parentGroups[i], p, aligns[i])
                mergedGroup = mergedGroup + parentGroups[i]
            parentGroups = [i for j, i in enumerate(parentGroups)\
                            if j not in mergeGroups ]
            
            mergedGroup = mergedGroup + [p]
            parentGroups.append(mergedGroup)
            alignedPartGroups = parentGroups + otherGroups
         
        # The alignment procedure above can move the fast varying
        # dimension outside. This has to be fixed.


    def isCompSelfDependent(self, comp):
        parts = self.polyParts[comp]        
        for p in parts:
            if self.isPartSelfDependent(p):
                return True
        return False

    def isPartSelfDependent(self, part):
        refs    = part.getPartRefs()
        objRefs = [ ref.objectRef for ref in refs\
                         if ref.objectRef == part.comp]
        if len(objRefs) > 0:
            return True
        return False

    def isGroupDependentOnPart(self, group, parentPart):
        for part in group:
            refs = part.getPartRefs()
            # This can be more precise
            objRefs = [ ref.objectRef for ref in refs\
                         if ref.objectRef == parentPart.comp]
            if len(objRefs) > 0:
                return True
        return False

    def isParent(self, part1, part2): 
        refs    = part2.getPartRefs()
        objRefs = [ ref.objectRef for ref in refs\
                         if ref.objectRef == part1.comp]
        if len(objRefs) > 0:
            return True
        return False

    def estimateReuse(self, part1, part2):
        """Gives an estimate of data reuse between the two poly parts."""
        # A very naive cost model
        if self.isParent(part1, part2):
            return 1
        return 0

    def getDomainDimCoeffs(self, sched, arg):
        domDimCoeff = {}
        if (isAffine(arg)):
            coeff = getAffineVarAndParamCoeff(arg)
            for item in coeff:
                if type(item) == Variable:
                    dim = sched.find_dim_by_name(isl._isl.dim_type.in_,
                                                 item.name)
                    domDimCoeff[dim] = coeff[item]
        return domDimCoeff

    def getParamCoeffs(self, sched, arg):
        paramCoeff = {}
        if (isAffine(arg)):
            coeff = getAffineVarAndParamCoeff(arg)
            for item in coeff:
                if type(item) == Parameter:
                    dim = sched.find_dim_by_name(isl._isl.dim_type.param,
                                                 item.name)
                    paramCoeff[dim] == coeff[item]
        return paramCoeff

    def computeRelativeScalingFactors(self, parentPart, childPart):
        """Computes the relative scaling factor in each dimension necessary to 
           uniformize dependencies between two parts. If the dependencies are 
           already uniform the scaling factor for the dimension is set to 1. If
           the dependencies cannot be uniformized along a particular dimension 
           the scaling factor for that dimension is set to *."""
        
        assert self.isParent(parentPart, childPart)

        # Scaling and offset factors structure.
        # ([scale factors], [offsets])
        #
        # [scale factors] determine the amount by which the dimension has to be
        # scaled to uniformize the dependencies. Each reference to the parent is
        # considered while determing the scaling factors. All the references 
        # should have the same scaling factor in a particular dimension otherwise
        # the scaling factor for the dimension cannot be determined uniquely.
        #
        # [offsets] specify the shift in each dimension that is require to 
        # uniformize dependencies. Simliar to dimension scale factors all 
        # offsets for a dimension should agree otherwise the offset for the 
        # dimension cannot be determied uniquely.
        refs       = childPart.getPartRefs()
        # Filtering out self references
        parentRefs = [ ref for ref in refs \
                       if ref.objectRef == parentPart.comp] 

        dimIn = childPart.sched.dim(isl._isl.dim_type.in_)
        scale = [ '-' for i in xrange(0, dimIn) ]
        offset = [ '-' for i in xrange(0, dimIn) ]

        def findDimScheduledTo(part, schedDim):
            for i in xrange(0, len(part.align)):
                if part.align[i] == schedDim:
                    return i
            return -1

        for ref in parentRefs:
            numArgs = len(ref.arguments)
            for i in xrange(0, numArgs):
                arg = ref.arguments[i]
                parentVarSchedDim = parentPart.align[i]            
                if (isAffine(arg)):
                    domDimCoeff = self.getDomainDimCoeffs(childPart.sched, arg)
                    paramCoeff = self.getParamCoeffs(childPart.sched, arg)
                    # Parameter coefficents can also be considered to
                    # generate parametric shifts. Yet to be seen.
                    
                    # Indexed with multiple variables.
                    if (len(domDimCoeff) > 1 or 
                        (len(domDimCoeff) == 1 and len(paramCoeff) >=1)):
                        # Although there are multiple coefficients. If 
                        # there is only one variable coefficient and other
                        # parametric coefficients. Uniformization can be 
                        # done with parametric shifts. Full affine scheduling 
                        # might be able to find a way to uniformize 
                        # dependencies. This has to be further explored.
                        dim = findDimScheduledTo(childPart, parentVarSchedDim)
                        if dim != -1:
                            scale[dim] = '*'
                    # Indexed with a single variable. This can either be an 
                    # uniform access or can be uniformized with scaling when 
                    # possible.
                    elif len(domDimCoeff) == 1 and len(paramCoeff) == 0:
                        dim = (domDimCoeff.keys())[0]
                        childVarSchedDim = childPart.align[dim]
                        # Checking if the schedule dimensions match only 
                        # then can the dependence be uniformized.
                        if childVarSchedDim != parentVarSchedDim:
                            continue
                        if scale[dim] == '-':
                            scale[dim] = domDimCoeff[dim]
                        elif scale[dim] != domDimCoeff[dim]:
                            scale[dim] = '*'
                    elif len(domDimCoeff) == 0 and len(paramCoeff) > 0:
                        continue
                    # Only parametric or Constant access. The schedule in this 
                    # dimension can be shifted to this point to uniformize the 
                    # dependence.
                    # In case the dimension in the parent has a constant size
                    # an upper and lower bound on the dependence vector can 
                    # be computed.
                    elif len(domDimCoeff) + len(paramCoeff) == 0:
                        # offsets should be set here
                        continue
                else:
                    dim = findDimScheduledTo(childPart, parentVarSchedDim)
                    if dim != -1:
                        scale[dim] = '*'

        for i in xrange(0, dimIn):
            if scale[i] == '-':
                scale[i] = 1
            if offset[i] == '-':
                offset[i] = 0
        return (scale, offset)

    def findParentGroups(self, part, groups):
        parentParts = []
        for comp in self.polyParts:
            for p in self.polyParts[comp]:
                if self.isParent(p, part):
                    parentParts.append(p)
        parentGroups = []        
        for p in parentParts:        
            for g in groups:
                if (p in g) and (g not in parentGroups):
                    parentGroups.append(g)
        return parentGroups

    def findChildGroups(self, parentGroup, groups):
        childGroups = []
        for p in parentGroup:
            for child in groups:
                if child != parentGroup and\
                  self.isGroupDependentOnPart(child, p) and\
                  child not in childGroups:
                    childGroups.append(child)
        return childGroups             

    def findLeafGroups(self, groups):
        leafGroups = []
        for i in xrange(0, len(groups)):
            isLeaf = True
            for p in groups[i]:
                for j in xrange(0, len(groups)):
                    if j!=i and self.isGroupDependentOnPart(groups[j], p):
                        isLeaf = False
            if isLeaf:
                leafGroups.append(groups[i])
        return leafGroups

    def simpleGroup(self, paramEstimates, single = True):
        compObjs = self.stage.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        # Create a reuse matrix among the poly parts
        assert(len(sortedCompObjs) >= 1)
        parts = []
        groups = []
        if single:
            group = []
            for comp in [item[0] for item in sortedCompObjs]:
                for part in self.polyParts[comp]:
                    group.append(part)
            if group:
                groups.append(group)
        else:    
            for comp in [item[0] for item in sortedCompObjs]:
                if self.isCompSelfDependent(comp) or True:
                    group = []
                    for part in self.polyParts[comp]:
                        group.append(part)
                    if group:
                        groups.append(group)
                else:
                    for part in self.polyParts[comp]:
                        groups.append([part])
        return groups        
    
    def createGroupScaleMap(self, group):
        groupScaleMap = {}
        for part in group:
            groupScaleMap[part] = list(part.scale)
        return groupScaleMap

    def isStencilGroup(self, group):
        depVecsGroup = self.getGroupDependenceVectors(group)
        for depVec, h in depVecsGroup:
            if '*' in depVec:
                return False
        return True    

    def getPartSize(self, part, paramEstimates):
        size = None
        domain = part.comp.domain
        if isinstance(part.comp, Accumulator):
            domain = part.comp.reductionDomain
        for interval in domain:
            subsSize = self.getDimSize(interval, paramEstimates)
            if isConstantExpr(subsSize):
                if size is None:
                    size = getConstantFromExpr(subsSize)
                else:
                    size = size * getConstantFromExpr(subsSize)
            else:
                size = '*'
                break
        assert size is not None
        return size
    
    def getDimSize(self, interval, paramEstimates):
        paramValMap = {}
        for est in paramEstimates:
            assert isinstance(est[0], Parameter)
            paramValMap[est[0]] = Value.numericToValue(est[1])

        if interval.step.value > 0:
            dimSize = interval.upperBound - interval.lowerBound + 1
        else:
            dimSize = interval.lowerBound - interval.upperBound + 1
        return substituteVars(dimSize, paramValMap)

    def groupStages(self, paramEstimates):
        compObjs = self.stage.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        # Create a reuse matrix among the poly parts
        parts = []
        for comp in [item[0] for item in sortedCompObjs]:
            parts = parts + self.polyParts[comp]
        reuse = {}
        relScaleFactors = {}
        for i in xrange(0, len(parts)):
            for j in xrange(0, len(parts)):
                reuse[(parts[i], parts[j])] = self.estimateReuse(parts[i], parts[j])
                # Skip scaling factor computation for reductions
                isReduction = isinstance(parts[i].comp, Accumulator) or \
                              isinstance(parts[j].comp, Accumulator)
                if (reuse[(parts[i], parts[j])] > 0) and (parts[j] != parts[i]) and \
                        not isReduction:
                    relScaleFactors[(parts[i], parts[j])] = \
                        self.computeRelativeScalingFactors(parts[i], parts[j])

        def scaleToParentGroup(part, group):
            refs       = part.getPartRefs()
            # Avoiding Input references this should be revisited at some point
            parentRefs = [ ref for ref in refs \
                           if ref.objectRef in self.polyParts ]
            dimIn = part.sched.dim(isl._isl.dim_type.in_)
            newScale = [ '-' for i in xrange(0, dimIn)]
            for ref in parentRefs:
                parentParts = self.polyParts[ref.objectRef]
                groupParentParts = [ pp for pp in parentParts \
                                     if pp in group ]
                if not groupParentParts:
                    continue
                for gpp in groupParentParts:
                    relScale, offset = relScaleFactors[(gpp, part)]
                    for i in xrange(0, dimIn):
                        # Check if the dimension can be scaled
                        if relScale[i] == '*':
                            newScale[i] = '*'
                        else:
                            # Get the current scaling of the dimension
                            # aligned to.
                            alignDim = None
                            for j in xrange(0, len(gpp.align)):
                                if gpp.align[j] == part.align[i]:
                                    if alignDim is None:
                                        alignDim = j
                                    else:
                                        # A dimension cannot be aligned to 
                                        # multiple schedule dimensions.
                                        assert False
                            if alignDim is not None:
                                currScale = Fraction(gpp.scale[alignDim],
                                                     part.scale[i])
                                if newScale[i] == '-':
                                    newScale[i] = relScale[i] * currScale
                                elif newScale[i] != currScale * relScale[i]:
                                    newScale[i] = '*'
            for i in xrange(0, dimIn):
                if newScale[i] == '-':
                    newScale[i] = 1
            return newScale                

        def isGroupWithPartStencil(group, part, scale):
            if '*' in scale:
                return False
            sm = self.createGroupScaleMap(group)            
            scaleGroupToPart(group, part, scale, sm)
            sm[part] = list(part.scale)
            assert part not in group           
            normalizeGroupScaling(group + [part], sm)
            depVecs = self.getGroupDependenceVectors(group + [part], sm)

            for vec, h in depVecs:
                if '*' in vec:
                    return False
            return True

        def getScale(sMap, p, i):
            if sMap is not None:
                return sMap[p][i]
            return p.scale[i]

        def setScale(sMap, p, i, val):
            if sMap is not None:
                sMap[p][i] = val
            else:
                p.scale[i] = val

        def scaleGroupToPart(group, part, scale, scaleMap = None):
            for g in group:
                dimIn = g.sched.dim(isl._isl.dim_type.in_)
                for j in xrange(0, len(scale)):
                    if scale[j] != '*':
                        scaled = False
                        for i in xrange(0, dimIn):
                            if g.align[i] == part.align[j]:
                                assert not scaled
                                scaled = True
                                s = Fraction(getScale(scaleMap, g, i), 
                                             scale[j])
                                setScale(scaleMap, g, i, s)
        
        def normalizeGroupScaling(group, scaleMap = None):
            dimOut = group[0].sched.dim(isl._isl.dim_type.out)
            norm = [ 1 for i in xrange(0, dimOut)]
            for g in group:
                dimIn = g.sched.dim(isl._isl.dim_type.in_)
                for j in xrange(0, dimOut):
                    scaled = False
                    for i in xrange(0, dimIn):
                        if g.align[i] == j:
                            assert not scaled
                            scaled = True
                            d = Fraction(getScale(scaleMap, g, i)).denominator
                            norm[j] = lcm(norm[j], d)
            # Normalizing the scales of the group
            for g in group:
                dimIn = g.sched.dim(isl._isl.dim_type.in_)
                dimOut = g.sched.dim(isl._isl.dim_type.out)
                for j in xrange(0, len(norm)):
                    for i in xrange(0, dimIn):
                        if g.align[i] == j:
                            s = getScale(scaleMap, g, i)
                            setScale(scaleMap, g, i, s * norm[j])

        def getGroupCost(group):
            return 1

        def estimateGroupCostWithBundle(group, bundle, scale, partSizeMap):
            return 1

        # Progressive algorithm for grouping stages assigns groups level by 
        # level trying to maximze reuse and size of stencil groups.

        # Filter out small computations. We characterize a computation as 
        # small when the domains of parents and the computation itself is 
        # small to benefit from tiling or parallelization. The parents are
        # included so that we do not miss out on storage optimzations.
        def getSmallComputations(pts, estimates):
            smallParts = []
            # This can be more precise but for now just estimating the 
            # size of the part by the size of the computation object.

            # Currentlty the size is just being estimated by the number
            # of points in the domain. It can be more accurately done 
            # by considering the arithmetic intensity in the expressions.
            partSizeMap = {}
            for p in pts:
                partSizeMap[p] = self.getPartSize(p, estimates)

            for i in xrange(0, len(pts)):
                smallPart = False
                if partSizeMap[pts[i]] != '*':
                    smallPart = partSizeMap[pts[i]] <= self.sizeThreshold
                for j in xrange(i+1, len(pts)):
                    if self.isParent(pts[j], pts[i]):
                        if partSizeMap[pts[j]] != '*':
                            smallPart = smallPart and \
                                        (partSizeMap[pts[j]] > self.sizeThreshold)
                        else:
                            smallPart = False
                if smallPart:
                    smallParts.append(pts[i])
            
            return smallParts, partSizeMap 

        smallParts, partSizeMap = getSmallComputations(parts, paramEstimates)

        # All the parts of a computation which has self dependencies should be
        # in the same group. Bundle such parts together.
        smallGroups = []
        optGroups = self.simpleGroup(paramEstimates, single=False)
        #print "intial groups begin"
        #for g in optGroups:
        #    for p in g:
        #        print p.comp.name,
        #    print
        #print "intial groups end"
        initialParts = 0
        for g in optGroups:
            initialParts += len(g)
        for g in optGroups:
            smallGroup = True
            for p in g:
                if not p in smallParts:
                    smallGroup = False
            if smallGroup:
                smallGroups.append(g)

        opt = True
        while opt:
            children = {}
            opt = False
            for gi in xrange(0, len(optGroups)):
                children[gi] = self.findChildGroups(optGroups[gi], optGroups)
            newGroups = [ group for group in optGroups ]
            for gi in children:
                isSmall = True
                isReduction = False
                for p in optGroups[gi]:
                    if not p in smallParts:
                        isSmall = False
                    if isinstance(p.comp, Accumulator):
                        isReduction = True
                #print len(children[gi]), isSmall
                if not isSmall and not isReduction and len(optGroups[gi]) < self.groupSize:
                    if (len(children[gi]) > 1) and False:
                        newGroup = [ p for p in optGroups[gi] ]
                        merge = True
                        for childGroup in children[gi]:
                            # Check if all the children can be fused
                            parentGroups = []
                            for cp in childGroup:
                                pgs = self.findParentGroups(cp, newGroups)
                                for pg in pgs:
                                    if pg not in parentGroups:
                                        parentGroups.append(pg)
                            for pg in parentGroups:
                                if pg not in children[gi] and pg != optGroups[gi]:
                                    merge = False
                        if merge:
                            print "parent group begin"
                            for p in optGroups[gi]:
                                print p.comp.name,
                            print
                            print "parent groups end"
                            print "adding groups begin"
                            for g in children[gi]:
                                for p in g:
                                    print p.comp.name,
                                print
                            print "adding groups end"
                            for childGroup in children[gi]:
                                for p in childGroup:
                                    scale = scaleToParentGroup(p, newGroup)
                                    scaleGroupToPart(newGroup, p, scale)
                                newGroup = newGroup + childGroup
                                newGroups.remove(childGroup)
                            newGroups.remove(optGroups[gi])
                            newGroups.append(newGroup)
                            opt = True
                            break
                    elif (len(children[gi]) == 1):
                        print "parent group begin"
                        for p in optGroups[gi]:
                            print p.comp.name,
                        print
                        print "parent groups end"
                        print "adding groups begin"
                        for g in children[gi]:
                            for p in g:
                                print p.comp.name,
                            print
                        print "adding groups end"
                        newGroup = [ p for p in optGroups[gi] ]
                        for p in children[gi][0]:
                            scale = scaleToParentGroup(p, newGroup)
                            scaleGroupToPart(newGroup, p, scale)
                        newGroup = newGroup + children[gi][0]
                        newGroups.remove(children[gi][0])
                        newGroups.remove(optGroups[gi])
                        newGroups.append(newGroup)
                        opt = True
                        break
            optGroups = newGroups         
        #print "final groups begin"
        #for g in optGroups:
        #    for p in g:
        #        print p.comp.name,
        #    print
        #print "final groups end"    
        #print optGroups
        #bundles = [ b for b in bundles if b not in smallGroups ]
        """
        for b in bundles:
            parentGroups = []
            # Find the leaf parents in the parent groups. This avoids cycles
            # by construction. The fact that a there are two separate dependent
            # parent groups means that fusing them was deemed suboptimal. They
            # are not considered for fusion again. This would leave only the 
            # leaf group as a choice for fusion.

            #Filter out non leaf groups
            leafGroups = self.findLeafGroups(optGroups + smallGroups)
            for p in b:
                parents = self.findParentGroups(p, leafGroups)
                parentGroups = parentGroups + parents

            # Filter out the small groups from the parent groups. Since they
            # are not considered for fusion.
            parentGroups = [ g for g in parentGroups \
                             if g not in smallGroups] 
            otherGroups = [ g for g in optGroups \
                            if g not in parentGroups ]
            scales = {}

            for i in xrange(0, len(parentGroups)):
                scaleList = []
                for p in b:
                    scaleList.append(scaleToParentGroup(p, parentGroups[i]))
                scales[i] = scaleList    
            
            # Compute the cost of adding the bundle to each of the candidate
            # groups. The approach we take now is to either add the bundle to
            # one of the groups or to merge all of the groups. Parital merging 
            # is not explored. 

            # The characteristics of a group are captured in the following 
            # structure - (tile dims, storage, overlap)
            # (tile dims) indicates the number of dimensions that can be tiled 
            # profitably. Fusion which reduces the number of tile dims is avoided. 
            # (storage) gives an estimate of the storage savings original function
            # sizes vs scratch or modulo buffer sizes.
            # (overlap) gives an estimate of the fraction of redundant computation
            # relative to tile sizes that will be done in the group.
           
            def isProfitable(cCost, nCost):
                return True

            candGroups = []
            for i in xrange(0, len(parentGroups)):
                currCost = getGroupCost(parentGroups[i])
                newCost = estimateGroupCostWithBundle(parentGroups[i], b,
                                                      scales[i], partSizeMap)
                if isProfitable(currCost, newCost):
                    candGroups.append(i)
            
            mergedGroup = []
            for i in candGroups:
                for j in xrange(0, len(b)):
                    scaleGroupToPart(parentGroups[i], b[j], scales[i][j])
                mergedGroup = mergedGroup + parentGroups[i]
            parentGroups = [i for j, i in enumerate(parentGroups)\
                            if j not in candGroups ]
            
            mergedGroup = mergedGroup + b
            parentGroups.append(mergedGroup)
            optGroups = parentGroups + otherGroups
        """ 
        for group in optGroups:
            normalizeGroupScaling(group)
       
        finalParts = 0
        for group in optGroups:
            finalParts += len(group)
        
        assert initialParts == finalParts
        return optGroups

    def baseSchedule(self, paramEstimates):
        self.alignParts()
        stageGroups = self.groupStages(paramEstimates)
        #stageGroups = self.simpleGroup(paramEstimates, single = False)

        for group in stageGroups:
            time = {}

            for part in group:
                time[part] = 0

            change = True
            while change:
                change = False
                for part in group:
                    for parent in group:
                        if self.isParent(parent, part):
                            if self.isParent(part, parent):
                                continue
                            elif time[part] <= time[parent]:
                                time[part] = time[parent] + 1
                                change = True

            for part in group:
                part.stageNo = time[part]
                dimIn = part.sched.dim(isl._isl.dim_type.in_)
                dimOut = part.sched.dim(isl._isl.dim_type.out)
                [ineqs, eqs] = formatScheduleConstraints(dimIn, dimOut, 
                                                         part.align, 
                                                         part.scale,
                                                         part.stageNo)
                part.sched = addConstriants(part.sched, ineqs, eqs)
        self.groups = stageGroups
        return stageGroups

    def simpleSchedule(self, paramEstimates):
        """Generate a simple schedule for the stage."""
        stageGroups = self.baseSchedule(paramEstimates)
        for i in xrange(0, len(stageGroups)):
            for p in stageGroups[i]:
                p.liveout = True           

    def computeTileSlope(self, stageDeps, hmax):
        # Compute slopes
        # -- The first dimension in the domain gives the stage order. The slope of
        #    the tile in each dimension is computed with respect to the stage order.
        #    The min extent and max extent in the each dimension are computed. The
        #    hyperplanes representing the min and max extent give the shape of the
        #    tile in that dimension.
        if len(stageDeps) < 1 :
            return ([], [])

        vecLen = len(stageDeps[0][0])
        slopeMin = [ (0, 1) for i in xrange(0, vecLen - 1) ]
        slopeMax = [ (0, 1) for i in xrange(0, vecLen - 1) ]
        # Find max and min widths of dependencies at the base
        widths = []
        hmin = min([ dep[1] for dep in stageDeps ])
        minWidth = [ 0 for i in xrange(0, vecLen - 1)]
        maxWidth = [ 0 for i in xrange(0, vecLen - 1)]
        depUnknown = [ False for i in xrange(0, vecLen - 1) ] 
        for currh in xrange(hmax - 1, hmin - 1, -1):
            maxW = [ 0 for i in xrange(0, vecLen - 1)]
            minW = [ 0 for i in xrange(0, vecLen - 1)]
            hDepVecs = [ depVec for depVec in stageDeps if \
                         depVec[1] == currh]
            for depVec, h in hDepVecs:             
                for i in xrange(0, len(depVec)-1):
                    if depVec[i+1] == '*':
                        depUnknown[i] = True
                        continue
                    if depVec[i+1] > 0:
                        maxW[i] = max(maxW[i], depVec[i+1])
                    if depVec[i+1] < 0:
                        minW[i] = min(minW[i], depVec[i+1])
            for i in xrange(0, len(depVec)-1):
                minWidth[i] = minWidth[i] + minW[i]
                maxWidth[i] = maxWidth[i] + maxW[i]
            widths.append((list(minWidth), currh))
            widths.append((list(maxWidth), currh))
            #print widths
                    
        for width, h in widths:
            scale = hmax - h 
            for i in xrange(0, vecLen-1):  
                if ((Fraction(width[i], scale) < 
                     Fraction(slopeMin[i][0], slopeMin[i][1])) and width[i] < 0):
                    slopeMin[i] = (width[i], scale)
                if ((Fraction(width[i], scale) >  
                     Fraction(slopeMax[i][0], slopeMax[i][1])) and width[i] > 0):
                    slopeMax[i] = (width[i], scale)

        for i in xrange(0, vecLen-1):             
            if depUnknown[i]:
                slopeMin[i] = '*'
                slopeMax[i] = '*'

        return (slopeMin, slopeMax)           

    def fusedSchedule(self, paramEstimates):
        """Generate an optimized schedule for the stage."""
        # Overall Approach
        # -- Partition the stages into groups
        #    -- Group together stages which have uniform dependencies across 
        #       or dependencies that can be uniformized.
        #    -- Try to group stages which only have inter-stage dependencies.
        #    -- Intra-stage dependencies are dealt separately. Since they 
        #       generally inhibit concurrent start.
        #    -- While grouping the stages use scaled schedules to uniformize
        #       dependencies. Algorithm to determine scaling factors.
        #    -- Align the dimensions of stages based on the parameters defining
        #       the dimension as well as subsequent access of the dimension.
        #       Can this be done while extracting the polyhedral representation?
        #    -- Try to reduce the live-range of the stages while grouping. A
        #       stage which only has consumers within the group can be optimized
        #       for storage.
        #    -- Use the estimates of input sizes and number of threads to formulate        
        #       simple heuristics.
        stageGroups = self.baseSchedule(paramEstimates)
        # -- Compute dependencies
        #    -- Since partitioning might introduce scaling factors. The 
        #       dependencies have to be computed based on the schedule
        #       How to extract dependence vectors from dependence polyhedra?
        #       This might be a better approach than trying to finding the vectors
        #       in an independent step.
        stageDeps = {}
        for i in xrange(0, len(stageGroups)):            
            stageDeps[i] = self.getGroupDependenceVectors(stageGroups[i])
            #for g in stageGroups[i]:
            #    print g.sched
            #    print g.expr
            #print stageDeps[i]

        # -- Generate a tiled schedule for the group
        #    -- Stencil groups are groups which have only uniform inter stage 
        #       dependencies. These stages can be tiled using the overlap or split
        #       tiling approach.
        #    -- Intra tile uniform dependencies should be tiled in a pipeline 
        #       fashion. Can this be folded into the overlap or split tiling strategy
        #       or needs to be dealt separately?
        #       Integral images and time iterated computations are important patterns 
        #       that fall into this category.
        #    -- For general affine dependencies the pluto algorithm should be used.
        #       We currently do not focus on general affine dependencies.
        stencilGroups = []
        for i in xrange(0, len(stageGroups)):
            # No point in tiling a group that has no dependencies
            isStencil = len(stageDeps[i]) > 0 and len(stageGroups[i]) > 1
            for dep, h in stageDeps[i]:
                # Skips groups which have self deps
                if dep[0] == 0:
                    isStencil = False
            if isStencil:
                stencilGroups.append(i)
            else:
                for p in stageGroups[i]:
                    partSize = self.getPartSize(p, paramEstimates)
                    bigPart = partSize != '*' and partSize > self.sizeThreshold/2
                    if not self.isPartSelfDependent(p) and bigPart:
                        # Determine the outer most dim and mark it parallel
                        # the inner most dim and mark it as vector
                        parallelDim = None
                        vecDim = None
                        for dom in xrange(0, len(p.align)):
                            interval = p.comp.domain[dom]
                            if isinstance(p.comp, Accumulator):
                                interval = p.comp.reductionDomain[dom]
                            # Since size could be estimated so can interval
                            # size no need to check.
                            intSize = self.getDimSize(interval, paramEstimates)
                            if(getConstantFromExpr(intSize) >= 32):
                                if parallelDim is not None:
                                    parallelDim = min(p.align[dom], parallelDim)
                                else:
                                    parallelDim = p.align[dom]
                                    
                            if(getConstantFromExpr(intSize) >= 4):
                                if vecDim is not None:
                                    vecDim = max(p.align[dom], vecDim)
                                else:
                                    vecDim = p.align[dom]
                        if parallelDim is not None:
                            pDimName = p.sched.get_dim_name(isl._isl.dim_type.out,
                                                           parallelDim)
                            p.parallelSchedDims.append(pDimName)
                        if vecDim is not None:
                            vDimName = p.sched.get_dim_name(isl._isl.dim_type.out,
                                                           vecDim)
                            p.vectorSchedDim.append(vDimName)

            # Find the stages which are not liveout
            maxStage = max([ p.stageNo for p in stageGroups[i] ])
            for p in stageGroups[i]:
                isLiveOut = not isStencil
                #isLiveOut = True
                for gn in xrange(0, len(stageGroups)):
                    if gn != i:
                        isLiveOut = isLiveOut or self.isGroupDependentOnPart(
                                                            stageGroups[gn], p)                        
                if p.stageNo == maxStage:        
                    p.liveout = True
                p.liveout = p.liveout or isLiveOut     
                    
        for gi in stencilGroups:
            assert(len(stageGroups[gi]) > 1)
            hmax = max( [ s.stageNo for s in stageGroups[gi] ] )
            hmin = min( [ s.stageNo for s in stageGroups[gi] ] )
            slopeMin, slopeMax = self.computeTileSlope(stageDeps[gi], hmax)
            #print slopeMin, slopeMax, hmax - hmin
            
            #self.splitTile(stageGroups[gi], slopeMin, slopeMax)
            self.overlapTile(stageGroups[gi], slopeMin, slopeMax)
            print stageDeps[gi]
            print slopeMin, slopeMax, hmax, len(stageGroups[gi])
            #for p in stageGroups[gi]:
            #    print p.scale, p.comp.name + ' = ' +  p.expr.__str__()
            #for p in stageGroups[gi]:
            #    print p.dimTileInfo

            # Determine the buffer sizes for stages in each dimension
            for p in stageGroups[gi]:
                for dom in p.dimTileInfo:
                    if p.dimTileInfo[dom][0] != 'none': 
                        dimName = p.dimTileInfo[dom][1]
                        tileDimName = p.dimTileInfo[dom][2]
                        extent = p.dimTileInfo[dom][3]
                        if p.dimTileInfo[dom][0] == 'overlap':
                            # Accounting for the overlap region
                            L = p.dimTileInfo[dom][4]
                            R = p.dimTileInfo[dom][5]
                            h = p.dimTileInfo[dom][6]
                            extent += abs(L * h) + abs(R * h)
                            baseWidth = h - p.stageNo
                            #extent += abs(L * h) + abs(R * baseWidth) 
                        p.dimScratchSize[dom] = \
                            int(math.ceil(Fraction(extent, p.scale[dom])))
                        #accName = '_Acc_' + p.sched.get_dim_name(isl._isl.dim_type.in_, dom)
                        #remName = '_Rem_' + p.sched.get_dim_name(isl._isl.dim_type.in_, dom)
                        mulName = '_Mul_' + p.sched.get_dim_name(isl._isl.dim_type.in_, dom)
                        dimIn = p.sched.dim(isl._isl.dim_type.in_)
                        domId =  p.sched.get_tuple_id(isl._isl.dim_type.in_)
                        p.sched = p.sched.insert_dims(isl._isl.dim_type.in_, dimIn, 1)
                        p.sched = p.sched.set_tuple_id(isl._isl.dim_type.in_, domId)
                        #p.sched = p.sched.set_dim_name(isl._isl.dim_type.in_, dimIn, accName)
                        #p.sched = p.sched.set_dim_name(isl._isl.dim_type.in_, dimIn+1, remName)
                        p.sched = p.sched.set_dim_name(isl._isl.dim_type.in_, dimIn, mulName)
                        schedDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, dimName)
                        tileDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, tileDimName)
                        
                        eqs = []
                        coeff = {}
                        coeff[('in', dimIn)] = p.scale[dom]
                        coeff[('out', schedDim)] = -1
                        coeff[('out', tileDim)] = p.dimTileInfo[dom][3]
                        eqs.append(coeff)
                       
                        ineqs = []
                        #coeff = {}
                        #coeff[('in', dimIn+2)] = p.scale[dom]
                        #coeff[('in', dimIn+1)] = 1
                        #coeff[('in', dimIn)] = -1
                        #eqs.append(coeff)

                        #coeff = {}
                        #coeff[('in', dimIn+1)] = 1
                        #coeff[('constant', 0)] = 0
                        #ineqs.append(coeff)

                        #coeff = {}
                        #coeff[('in', dimIn+1)] = -1
                        #coeff[('constant', 0)] = p.scale[dom] - 1
                        #ineqs.append(coeff)

                        p.sched = addConstriants(p.sched, ineqs, eqs)
            
            # Second level storage savings can be achieved by utilizing modulo buffers
            # in the non-vector dimension. The fastest varying dimension is considered
            # the vector dimension and by this point should be the inner-most dimension.

            # Disabling this for two reasons
            # 1) The code generator generates awful code. There is no reason to expect
            #    it to generate anything nice.
            # 2) The dimension that has skewing applied to it need not be tiled. This 
            #    has to be integrated into scheduling itself.
            """
            for p in stageGroups[gi]:
                oneDim = True
                for dom in p.dimTileInfo:
                    if p.dimTileInfo[dom][0] == 'overlap' and oneDim:
                        oneDim = False
                        dimName = p.dimTileInfo[dom][1]
                        
                        # Skewing the dimension
                        schedDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, dimName)
                        p.sched = p.sched.insert_dims(isl._isl.dim_type.out, schedDim  + 1, 1)
                        p.sched = p.sched.set_dim_name(isl._isl.dim_type.out, 
                                                        schedDim + 1, '_shift' + dimName)
                        timeDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')

                        R = p.dimTileInfo[dom][5]
                        eqs = []
                        coeff = {}
                        coeff[('out', schedDim)] = 1
                        coeff[('out', timeDim)] = abs(R)
                        coeff[('out', schedDim + 1)] = -1
                        eqs.append(coeff)
                        p.sched = addConstriants(p.sched, [], eqs)
                        p.sched = p.sched.remove_dims(isl._isl.dim_type.out, schedDim, 1)
                        
                        # Moving time inside
                        timeDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')
                        p.sched = p.sched.insert_dims(isl._isl.dim_type.out, timeDim, 1)
                        p.sched = p.sched.set_dim_name(isl._isl.dim_type.out, 
                                                        timeDim, '_tmp' + dimName)
                        schedDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, '_shift' + dimName)

                        eqs = []
                        coeff = {}
                        coeff[('out', timeDim)] = 1
                        coeff[('out', schedDim)] = -1
                        eqs.append(coeff)
                        p.sched = addConstriants(p.sched, [], eqs)
                        p.sched = p.sched.remove_dims(isl._isl.dim_type.out, schedDim, 1)
                        p.sched = p.sched.set_dim_name(isl._isl.dim_type.out, 
                                                        timeDim, '_shift' + dimName)
            """
            # -- Mark parallel dimensions and vector dimensions in each group
            #    -- Find the outer most parallel dimension which can generate "enough"
            #       tasks for the given number of threads.
            #    -- Partial and full tile separation to enable better vectorization.
            #    -- We currently rely on compiler vectorization. This is quite unreliable.
            #       We need to revisit the vectorization strategy.
            for p in stageGroups[gi]:
                outerParallelDim = None
                innerVecDim = None
                for dom in p.dimTileInfo:
                    if p.dimTileInfo[dom][0] == 'none':
                        # Either the dimension is too small to be parallelized or 
                        # is skewed. In both cases the dimension cannot be parallel.
                        # This can change when we choose to not tile a dimension.
                        continue
                    elif p.dimTileInfo[dom][0] == 'overlap':
                        dimName = p.dimTileInfo[dom][1]
                        tileDimName = p.dimTileInfo[dom][2]
                        schedDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, 
                                                            dimName)
                        tileDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, 
                                                            tileDimName)
                        if outerParallelDim is not None:
                            outerParallelDim = min(tileDim, outerParallelDim)
                        else:
                            outerParallelDim = tileDim
                        if innerVecDim is not None:
                            innerVecDim = max(schedDim, innerVecDim)
                        else:
                            innerVecDim = schedDim

                if outerParallelDim is not None:
                    pDimName = p.sched.get_dim_name(isl._isl.dim_type.out,
                                                    outerParallelDim)
                    p.parallelSchedDims.append(pDimName)
                if innerVecDim is not None:
                    vDimName = p.sched.get_dim_name(isl._isl.dim_type.out,
                                                    innerVecDim)
                    p.vectorSchedDim.append(vDimName)

            # Computations which have different scale but map to the same time
            # generate a lot of conditionals which can hinder performance. This
            # step separates all computations in a time step by adding an additional 
            # dimension.
            compParts = {}
            for p in stageGroups[gi]:
                if p.comp in compParts:
                    compParts[p.comp].append(p)
                else:
                    compParts[p.comp] = [p]

            pi = 0
            for comp in compParts:
                for p in compParts[comp]:
                    timeDim = p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')
                    p.sched = p.sched.insert_dims(isl._isl.dim_type.out, timeDim + 1, 1)
                    p.sched = p.sched.set_dim_name(isl._isl.dim_type.out, 
                                               timeDim + 1, '_o')
                    eqs = []
                    coeff = {}
                    coeff[('constant', 0)] = -pi
                    coeff[('out', timeDim + 1)] = 1
                    eqs.append(coeff)
                    p.sched = addConstriants(p.sched, [], eqs)
                    pi += 1

            #for p in stageGroups[gi]:
            #    print p.sched
            #assert False


    def moveIndependentDim(self, dim, group, stageDim):
        # Move the independent dimensions outward of the stage dimension.
        for part in group:
            part.sched = part.sched.insert_dims(isl._isl.dim_type.out, 
                                                stageDim, 1)
            noDepId = part.sched.get_dim_id(
                            isl._isl.dim_type.out, dim + 1)
            noDepName = part.sched.get_dim_name(
                            isl._isl.dim_type.out, dim + 1)
            eqs = []
            coeff = {}
            coeff[('out', dim+1)] = -1
            coeff[('out', stageDim)] = 1
            eqs.append(coeff)
            part.sched = addConstriants(part.sched, [], eqs)
            part.sched = part.sched.remove_dims(
                                isl._isl.dim_type.out, dim+1, 1)
            part.sched = part.sched.set_dim_name(
                                    isl._isl.dim_type.out, 
                                    stageDim, noDepName)

    def getGroupHeight(self, group):
        minHeight = min( [ part.stageNo for part in group ] )
        maxHeight = max( [ part.stageNo for part in group ] )
        return maxHeight - minHeight

    def overlapTile(self, group, slopeMin, slopeMax):
        stageDim = 0
        tileDims = 0
        noTileDims = 0
        h = self.getGroupHeight(group)
        numTileDims = 0
        for i in xrange(1, len(slopeMin) + 1):                    
            # Check if every stage in the group has enough iteration 
            # points in the dimension to benefit from tiling.
            tile = False
            for part in group:
                currDim = stageDim + noTileDims + 2*tileDims + 1
                lowerBound = part.sched.range().dim_min(currDim)
                upperBound = part.sched.range().dim_max(currDim)
                size = upperBound.sub(lowerBound)
                if (size.is_cst() and size.n_piece() == 1):
                    aff = (size.get_pieces())[0][1]
                    val = aff.get_constant_val()
                    if val > self.tileSizes[numTileDims]:
                        tile = True
                else:
                    tile = True
            if tile and slopeMin[i-1] != '*':        
                # Altering the schedule by constructing overlapped tiles.
                for part in group:
                    # Extend space to accomodate the tiling dimensions
                    part.sched = part.sched.insert_dims(
                                    isl._isl.dim_type.out, 
                                    stageDim + tileDims, 1)
                    name = part.sched.get_dim_name(
                                isl._isl.dim_type.out, 
                                stageDim + noTileDims + 2*tileDims + 2)
                    part.sched = part.sched.set_dim_name(
                                    isl._isl.dim_type.out, 
                                    stageDim + tileDims, 
                                    '_T' + name)
                    R = int(math.floor(Fraction(slopeMin[i-1][0], 
                                                slopeMin[i-1][1])))
                    L = int(math.ceil(Fraction(slopeMax[i-1][0], 
                                               slopeMax[i-1][1])))
                    # L and R are normals to the left and the right 
                    # bounding hyperplanes of the uniform dependencies
                
                    tileSize = self.tileSizes[numTileDims]
                    # Compute the overlap shift
                    #print slopeMax, slopeMin, h, L, R, i-1
                    overlapShift = abs(L * (h)) + abs(R * (h))
                    for j in xrange(0, len(part.align)):
                        if i == part.align[j]:
                            assert j not in part.dimTileInfo
                            if tileSize%part.scale[j] != 0:
                                tileSize = int(math.ceil(part.scale[j]))
                            part.dimTileInfo[j] = ('overlap', name, '_T' + name, 
                                                     tileSize, L, R, h)
                    ineqs = []
                    eqs = []
                    coeff = {}
                    itDim = stageDim + noTileDims + 2*tileDims + 2
                    tileDim = stageDim + tileDims
                    timeDim = stageDim + tileDims + 1
                    coeff[('out', timeDim)] = -L
                    coeff[('out', itDim)] = 1
                    coeff[('out', tileDim)] = -tileSize
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('out', timeDim)] = L
                    coeff[('out', itDim)] = -1
                    coeff[('out', tileDim)] = tileSize
                    coeff[('constant', 0)] = tileSize - 1 + overlapShift
                    ineqs.append(coeff)
                
                    coeff = {}
                    coeff[('out', timeDim)] = -R
                    coeff[('out', itDim)] = 1
                    coeff[('out', tileDim)] = -tileSize
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('out', timeDim)] = R
                    coeff[('out', itDim)] = -1
                    coeff[('out', tileDim)] = tileSize
                    coeff[('constant', 0)] = tileSize + overlapShift - 1
                    ineqs.append(coeff)

                    priorDom = part.sched.domain()
                    part.sched = addConstriants(part.sched, ineqs, eqs)
                    postDom = part.sched.domain()
               
                    assert(part.sched.is_empty() == False)
                    # Tiling should not change the domain that is iterated over               
                    assert(priorDom.is_equal(postDom))
                tileDims += 1
                numTileDims += 1
            else:
                #self.moveIndependentDim(i, group, stageDim)
                name = part.sched.get_dim_name(isl._isl.dim_type.out, stageDim) 
                for part in group:                        
                    for j in xrange(0, len(part.align)):
                        if i == part.align[j]:
                            assert j not in part.dimTileInfo
                            part.dimTileInfo[j] = ('none', name)
                noTileDims += 1

    def splitTile(self, group, slopeMin, slopeMax):
        stageDim = 0
        dtileDims = 0
        numTileDims = 0
        for i in xrange(1, len(slopeMin) + 1):                    
            if ((slopeMin[i-1][0] != 0 or slopeMax[i-1][0] !=0)):
                # Altering the schedule by constructing split tiles.
                for part in group:
                    # Extend space to accomodate the tiling dimensions
                    part.sched = part.sched.insert_dims(
                                    isl._isl.dim_type.out, 
                                    stageDim + 2*dtileDims, 2)
                    # Dimension i is for the orientation of the tiles 
                    # upward or inverted.
                    name = part.sched.get_dim_name(
                                isl._isl.dim_type.out, 
                                stageDim + 3*dtileDims + 3)
                    part.sched = part.sched.set_dim_name(
                                    isl._isl.dim_type.out, 
                                    stageDim + 2*dtileDims + 1, 
                                    '_T' + name)
                    part.sched = part.sched.set_dim_name(
                                    isl._isl.dim_type.out, 
                                    stageDim + 2*dtileDims, 
                                    '_Dir' + name)
                    
                    L = (slopeMin[i-1][0], slopeMin[i-1][1])
                    R = (slopeMax[i-1][0], slopeMax[i-1][1])
                    # L and R are normals to the left and the right 
                    # bounding hyperplanes of the uniform dependencies
                    
        # Tile size
        #   -- Pick tile sizes such that there are only two sets of tiles 
        #      in the time sense .i.e there should be only one fused stage. 
        #      This has to be revisited when time iterated computations are 
        #      incorporated
                    #offset = 3*tileSize/4
                    tileSize = self.tileSizes[numTileDims]
                    offset = tileSize/2
                    ineqs = []
                    eqs = []
                    coeff = {}
                    coeff[('out', stageDim + 2*dtileDims + 2)] = L[0]
                    coeff[('out', stageDim + 3*dtileDims + 3)] = L[1]
                    coeff[('out', stageDim + 2*dtileDims + 1)] = -tileSize
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('out', stageDim + 2*dtileDims + 2)] = -L[0]
                    coeff[('out', stageDim + 3*dtileDims + 3)] = -L[1]
                    coeff[('out', stageDim + 2*dtileDims + 1)] = tileSize
                    coeff[('constant', 0)] = tileSize - 1
                    ineqs.append(coeff)
                    
                    coeff = {}
                    coeff[('out', stageDim + 2*dtileDims + 2)] = R[0]
                    coeff[('out', stageDim + 3*dtileDims + 3)] = R[1]
                    coeff[('out', stageDim + 2*dtileDims + 1)] = -tileSize 
                    coeff[('out', stageDim + 2*dtileDims)] = -tileSize 
                    coeff[('constant', 0)] = -offset
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('out', stageDim + 2*dtileDims + 2)] = -R[0]
                    coeff[('out', stageDim + 3*dtileDims + 3)] = -R[1]
                    coeff[('out', stageDim + 2*dtileDims + 1)] = tileSize
                    coeff[('out', stageDim + 2*dtileDims)] = tileSize 
                    coeff[('constant', 0)] = tileSize + offset - 1
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('out', stageDim + 2*dtileDims)] = 1
                    coeff[('constant', 0)] = 1
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('out', stageDim + 2*dtileDims)] = -1
                    coeff[('constant', 0)] = 0
                    ineqs.append(coeff)

                    #eqsUpward = eqs[:]
                    #eqsDown = eqs[:]
                    #coeff = {}
                    #coeff[('out', stageDim + 2*dtileDims)] = -1
                    #coeff[('constant', 0)] = 0
                    #eqsUpward.append(coeff)

                    #coeff = {}
                    #coeff[('out', stageDim + 2*dtileDims)] = 1
                    #coeff[('constant', 0)] = 1
                    #eqsDown.append(coeff)

                    #schedUp = addConstriants(part.sched, ineqs, eqsUpward)
                    #schedDown = addConstriants(part.sched, ineqs, eqsDown)                                
                    #part.sched = schedUp.union(schedDown)
                    part.sched = addConstriants(part.sched, ineqs, eqs)
                    assert(part.sched.is_empty() == False)
                dtileDims += 1
                numTileDims += 1
            else:
                stageDim = self.moveIndependentDim(i, group, stageDim)

    def buildAst(self):
        sortedGroups = []
        tempGroups = [ g for g in self.groups ]
        while tempGroups:
            leafGroups = self.findLeafGroups(tempGroups)
            for l in leafGroups:
                sortedGroups.insert(0, l)
            for l in leafGroups:
                tempGroups.remove(l)
        for group in sortedGroups:
            #astbld =  isl.AstBuild.from_context(isl.BasicSet("[C, R]->{: R>=1 and C>=1}", self.ctx))
            astbld =  isl.AstBuild.from_context(group[0].sched.params())
            #astbld =  astbld.set_options(isl.UnionMap("{ }"))

            schedMap = None
            optMap = None
            for part in group:
                if schedMap is None:
                    schedMap = isl.UnionMap.from_map(part.sched)
                else:
                    partMap = isl.UnionMap.from_map(part.sched)
                    schedMap = schedMap.union(partMap)
                if optMap is None:
                    #unrollUnionSet = isl.UnionSet.from_set(isl.Set("{unroll[x] : x = 0 or x = 2}", self.ctx))
                    unrollUnionSet = isl.UnionSet.from_set(isl.Set("{:}", self.ctx))
                    domUnionSet = isl.UnionSet.universe(isl.UnionSet.from_set(part.sched.range()))
                    optMap = isl.UnionMap.from_domain_and_range(domUnionSet, unrollUnionSet)
                else:
                    #unrollUnionSet = isl.UnionSet.from_set(isl.Set("{unroll[x] : x = 0 or x = 2}", self.ctx))
                    unrollUnionSet = isl.UnionSet.from_set(isl.Set("{:}", self.ctx))
                    domUnionSet = isl.UnionSet.universe(isl.UnionSet.from_set(part.sched.range()))
                    optMap = optMap.union(isl.UnionMap.from_domain_and_range(domUnionSet, unrollUnionSet))
            astbld = astbld.set_options(optMap)

            numIds = group[0].sched.dim(isl._isl.dim_type.out)
            ids = isl.IdList.alloc(self.ctx, numIds)
            for i in xrange(0, numIds):
                schedName = group[0].sched.get_dim_name(isl._isl.dim_type.out, i)
                ids = ids.add(isl.Id.alloc(self.ctx, schedName, None))
            astbld = astbld.set_iterators(ids)
            
            #assert False
            def printer(arg):
                print arg
            schedMap.foreach_map(printer)
            self.polyast.append(astbld.ast_from_schedule(schedMap))
       
    def getVarName(self):
        name = "_i" + str(self._varCount)
        self._varCount+=1
        return name

    def getFuncName(cls):
        name = "_f" + str(self._funcCount)
        self._funcCount+=1
        return name

def mapCoeffToDim(coeff):
    for var in coeff.keys():
        coeffval = coeff[var]
        coeff.pop(var)
        if (isinstance(var, Parameter)):
            coeff[('param', var.name)] = coeffval
        elif (isinstance(var, Variable)):
            coeff[('in', var.name)] = coeffval
    return coeff

def formatScheduleConstraints(dimIn, dimOut, align, scale, stageNo):
    ineqCoeff = []
    eqCoeff   = []
    dimSet = [ False for i in xrange(0, dimOut) ]
    for i in xrange(0, dimIn):
        coeff = {}
        # Have to convert the iterator to non-zero range. Currently if the 
        # step is negative the iterator will take negative values.
        coeff[('out', align[i])] = 1
        assert scale[i] >= 1
        coeff[('in', i)] = -1 * scale[i]
        eqCoeff.append(coeff)
        dimSet[align[i]] = True

    # Adding stage identity constraint
    stageCoeff = {}
    stageCoeff[('out', 0)] = -1
    stageCoeff[('constant', 0)] = stageNo
    eqCoeff.append(stageCoeff)

    for i in xrange(1, dimOut):
        if not dimSet[i]:
            coeff = {}
            coeff[('out', i)] = 1
            coeff[('constant', 0)] = 0
            eqCoeff.append(coeff)
    return [ineqCoeff, eqCoeff]     

def formatDomainConstraints(domain, varNames):
    ineqCoeff = []
    eqCoeff   = []
    domLen = len(domain)
    for i in xrange(0, domLen):
        coeff = {}
        interval = domain[i]

        lbCoeff = getAffineVarAndParamCoeff(interval.lowerBound)
        # Mapping from variable names to the corresponding dimension
        lbCoeff = mapCoeffToDim(lbCoeff)
        lbConst = getConstantFromExpr(interval.lowerBound, affine = True)

        # Normalizing into >= format
        coeff = dict( (n, -lbCoeff.get(n)) for n in  lbCoeff)
        coeff[('constant', 0)] = -lbConst
        coeff[('in', varNames[i])] = 1
        ineqCoeff.append(coeff)

        ubCoeff = getAffineVarAndParamCoeff(interval.upperBound)
        # M_apping from variable names to the corresponding dimension 
        ubCoeff = mapCoeffToDim(ubCoeff)
        ubConst = getConstantFromExpr(interval.upperBound, affine = True)

        # Normalizing into >= format
        coeff = ubCoeff
        coeff[('constant', 0)] = ubConst
        coeff[('in', varNames[i])] = -1
        ineqCoeff.append(coeff)

    return [ineqCoeff, eqCoeff]

def formatConjunctConstraints(conjunct):
    # Check if the condition is a conjunction
    ineqCoeff = []
    eqCoeff = []
    for cond in conjunct:
        coeff = {}
        leftCoeff = getAffineVarAndParamCoeff(cond.lhs)
        rightCoeff = getAffineVarAndParamCoeff(cond.rhs)
        leftConst = getConstantFromExpr(cond.lhs, affine = True)
        rightConst = getConstantFromExpr(cond.rhs, affine = True)

        # Mapping from variable names to the corresponding dimension 
        leftCoeff = mapCoeffToDim(leftCoeff)
        rightCoeff = mapCoeffToDim(rightCoeff)

        # Normalizing >= format
        if (cond.conditional in ['<=','<']):
            coeff = dict( (n, -leftCoeff.get(n, 0) + rightCoeff.get(n, 0))\
                          for n in set(leftCoeff)| set(rightCoeff) )
            coeff[('constant', 0)] = -leftConst + rightConst - int(cond.conditional == '<')
            ineqCoeff.append(coeff)
        elif(cond.conditional in ['>=','>']):
            coeff = dict( (n, leftCoeff.get(n, 0) - rightCoeff.get(n, 0))\
                           for n in set(leftCoeff)| set(rightCoeff) )
            coeff[('constant', 0)] = leftConst - rightConst + int(cond.conditional == '>')
            ineqCoeff.append(coeff)
        else:
            # Weird
            assert(cond.conditional == '==')
            coeff = dict( (n, leftCoeff.get(n, 0) - rightCoeff.get(n, 0))\
                          for n in set(leftCoeff)| set(rightCoeff) )
            coeff[('constant', 0)] = leftConst - rightConst
            eqCoeff.append(coeff)
    return [ineqCoeff, eqCoeff]

def checkRefs(childStage, parentStage):
    # Check refs works only on non-fused stages. It can be made to
    # work with fused stages as well. However, it might serve very
    # little use.
    assert (not childStage.isFused() and not parentStage.isFused())
    parentFunc = parentStage.computeObjs[0]
    childObj   = childStage.computeObjs[0]

    # Only verifying if both child and  parent stage have a polyhedral 
    # representation
    if childStage.polyRep.polyParts and parentStage.polyRep.polyDoms:
        for childPart in childStage.polyRep.polyParts[childObj]:
            # Compute dependence relations between child and parent
            childRefs = childPart.getPartRefs()
            if childPart.pred:
                childRefs += childPart.pred.collect(Reference)
            # It is not generally feasible to check the validity of
            # and access when the reference is not affine. 
            # Approximations can be done but for now skipping them.
            def affineParentRef(ref, parentFunc):
                affine = True
                for arg in ref.arguments:
                    affine = affine and isAffine(arg) 
                return affine and ref.objectRef == parentFunc    
            childRefs = [ ref for ref in childRefs if \
                            affineParentRef(ref, parentFunc)]

            deps = []
            for ref in childRefs:
                deps += extractValueDependence(childPart, ref, 
                             parentStage.polyRep.polyDoms[parentFunc])
            for dep in deps:
                diff = dep.rel.range().subtract(
                        parentStage.polyRep.polyDoms[parentFunc].domSet)
                if(not diff.is_empty()):
                    raise TypeError("Reference out of domain", childStage, 
                                     parentStage, diff)

def inline(childStage, parentStage, noSplit = False):
    refToInlineExprMap = {}
    # Inling currently only handles non-fused stages
    assert (not childStage.isFused() and not parentStage.isFused())
    # Computation object in the parent stage can only be a function     
    parentFunc = parentStage.computeObjs[0]
    assert isinstance(parentFunc, Function)
    childObj  = childStage.computeObjs[0]

    # Simple scalar functions which are defined on a non-bounded 
    # integer domain can be inlined.
    # TODO

    # Inlining only if both child and parent stage have a polyhedral
    # representation
    if childStage.polyRep.polyParts and parentStage.polyRep.polyParts: 
        for childPart in childStage.polyRep.polyParts[childObj]:
            # Compute dependence relations between child and parent
            childRefs = childPart.getPartRefs()
            if childPart.pred:
                childRefs += childPart.pred.collect(Reference)
            childRefs = [ ref for ref in childRefs if ref.objectRef == parentFunc]

            deps = []
            for ref in childRefs:
                deps += extractValueDependence(childPart, ref, 
                            parentStage.polyRep.polyDoms[parentFunc])
            
            # Check if all the values come from the same parent part
            depToPartMap = {}
            for dep in deps:
                accessRegion = dep.rel.range().copy().reset_tuple_id()
                diff = dep.rel.range().copy().reset_tuple_id()
                for parentPart in parentStage.polyRep.polyParts[parentFunc]:
                    partRegion = parentPart.sched.domain().copy().reset_tuple_id()
                    partdiff = accessRegion.subtract(partRegion)
                    diff = diff.subtract(partRegion)
                    if(partdiff.is_empty()):
                        depToPartMap[dep] = parentPart
                if (not diff.is_empty()):
                    assert False, "Inlining cannot be done."

            parts = list(set(depToPartMap.values()))
            singlePart = (len(parts) == 1)

            if(singlePart):
                parentExpr = parts[0].expr
                if parts[0].pred:
                    inline = False
                else:
                    inlineDeps = []
                    for ref in childRefs:
                        refToInlineExprMap[ref] = parentExpr
            elif(noSplit):
                pass
            else:
                pass
    else:
        pass
    return refToInlineExprMap

def getParamsInvolved(sched, dim):
    numParams = sched.domain().n_param()
    paramNames = [ ]
    constraints = sched.domain().get_constraints()
    dimName = sched.domain().get_dim_name(isl._isl.dim_type.set, dim)
    for p in xrange(0, numParams): 
        for const in constraints:
           pname = sched.domain().get_dim_name(isl._isl.dim_type.param, p)
           dimId  = const.get_space().find_dim_by_name(isl._isl.dim_type.set, dimName)
           pId = const.get_space().find_dim_by_name(isl._isl.dim_type.param, pname)
           dimInvolved = const.involves_dims(isl._isl.dim_type.set, dimId, 1)
           paramInvolved = const.involves_dims(isl._isl.dim_type.param, pId, 1)
           if (dimInvolved and paramInvolved):
                   if pname not in paramNames:
                       paramNames.append(pname)
    return paramNameis

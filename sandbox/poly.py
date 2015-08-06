from __future__ import absolute_import, division, print_function

import math

import islpy as isl

from constructs import *
from expression import *

def lcm(a, b):
    return a*b/(gcd(a, b))

def optimizeSchedule(partScheds, dependencies):
    # The pluto optimizer can be used to optimize the schedule for comparision.
    pass

def addConstraintsFromList(obj, localSpace, constraintList, constraintAlloc):
    for const in constraintList:
        c = constraintAlloc(localSpace)
        m = 1
        for coeff in const:
            if isinstance(const[coeff], Fraction):
                den = int(gcd(abs(const[coeff].denominator), m))
                m = (abs(const[coeff].denominator) * m)//den
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

def addConstraints(obj, ineqs, eqs):
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
    assert(part.sched)
    deps = []
    accessRegion = isl.BasicSet.universe(refPolyDom.domSet.get_space())
    partDom = part.sched.domain().align_params(refPolyDom.domSet.get_space())    
    accessRegion = accessRegion.align_params(partDom.get_space())
    
    rel = isl.BasicMap.from_domain_and_range(partDom, accessRegion)
    dimOut = rel.dim(isl._isl.dim_type.out)
    sourceDims = [ ('out', i) for i in range(0, dimOut)]
    numArgs = len(ref.arguments)                

    for i in range(0, numArgs):
        arg = ref.arguments[i]
        # If the argument is not affine the dependence reflects that
        # the computation may depend on any value of the referenced object
        if (isAffine(arg)):
            coeff = get_affine_var_and_param_coeff(arg)
            coeff = mapCoeffToDim(coeff)

            coeff[('constant', 0)] = get_constant_from_expr(arg, affine = True)
            coeff[sourceDims[i]] = -1
            rel = addConstraints(rel, [], [coeff])
    if not rel.is_empty():
        deps.append(PolyDep(ref.objectRef, part.comp, rel))
    return deps 

class PolyPart(object):
    def __init__(self, _schedMap, _expr, _pred, _comp, 
                 _align, _scale, _levelNo):
        self.schedMap = _schedMap
        self.expr = _expr
        self.pred = _pred
        self.comp = _comp
        self.sched = None
        # Dependencies between values of computation objects
        self.deps = []
        # Mapping between the input variables to the corresponding 
        # schedule dimension. A full affine schedule will need a 
        # transformation matrix. Currently we only shuffle the 
        # dimension order apart from tiling so a simple dimension
        # alignment vector suffices. This has to be changed to 
        # handle more general cases later.
        self._align = _align
        # Scaling factors for each schedule dimension
        self._scale = _scale
        # Default alignment and scaling factors are set while
        # constructing the polypart. These are changed by the
        # alignment and loop scaling passes. Both these passer
        # attempt to improve locality and uniformize dependencies.
        self.levelNo = _levelNo

    @property
    def align(self):
        align_clone = [i for i in self._align]
        return align_clone

    @property
    def scale(self):
        scale_clone = [i for i in self._scale]
        return scale_clone

    def set_align(self, align):
        self._align = [i for i in align]
        return

    def set_scale(self, _scale):
        self._scale = [i for i in _scale]
        return

    def getPartRefs(self):
        refs = self.expr.collect(Reference)
        if (self.pred):
            refs += self.pred.collect(Reference)
        return refs

    def __str__(self):
        partStr = "Schedule: " + self.sched.__str__() + '\n'\
                  "Expression: " + self.expr.__str__() + '\n'\
                  "Predicate: " + self.pred.__str__() + '\n'
        depstr = ""
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
        group. It gives piece-wise domain and schedule for each compute
        object in the group. Polyhedral transformations modify the 
        piece-wise domains as well as the schedules.
    """
    def __init__(self, _ctx, _group, _paramConstraints):
        self.group = _group
        self.ctx = _ctx
        self.polyParts = {}
        self.polyDoms = {}
        self.polyast = []

        self._varCount = 0
        self._funcCount = 0

        self.extractPolyRepFromGroup(_paramConstraints)

        #self.fusedSchedule(_paramEstimates)
        #self.simpleSchedule(_paramEstimates)

    def extractPolyRepFromGroup(self, paramConstraints):
        compObjs = self.group.orderComputeObjs()
        numObjs = len(compObjs.items())

        # Comute the max dimensionality of the compute objects
        def maxDim(objs):
            dim = 0
            for comp in objs:
                if type(comp) == Reduction:
                    dim = max(dim, len(comp.reductionVariables))
                    dim = max(dim, len(comp.variables))
                elif type(comp) == Function or type(comp) == Image:
                    dim = max(dim, len(comp.variables))
            return dim

        dim = maxDim(compObjs)
        # Get all the parameters used in the group compute objects
        params = []
        for comp in compObjs:
            params = params + comp.getObjects(Parameter)
        params = list(set(params))
        paramNames = [param.name for param in params]        

        # Represent all the constraints specified on the parameters relevant
        # to the group.
        contextConds = self.formatParamConstraints(paramConstraints, params)
                    
        # The [t] is for the stage dimension
        scheduleNames = ['_t'] + [ self.getVarName()  for i in range(0, dim) ]

        for comp in compObjs:
            if (type(comp) == Function or type(comp) == Image):
                self.extractPolyRepFromFunction(comp, scheduleNames, paramNames,
                                                contextConds, compObjs[comp] + 1,
                                                paramConstraints)
            elif (type(comp) == Reduction):
                self.extractPolyRepFromReduction(comp, scheduleNames, paramNames,
                                                 contextConds, compObjs[comp] + 1,
                                                 paramConstraints)
            else:
                assert False

    def formatParamConstraints(self, paramConstraints, params):
        contextConds = []
        for paramConst in paramConstraints:
            # Only consider parameter constraints of parameters
            # given in params.
            paramsInConst = paramConst.collect(Parameter)
            contextAdd = True
            for p in paramsInConst:
                if p not in params:
                    contextAdd = False
            # Only add the constraint if it is affine and has no conjunctions. 
            # Handling conjunctions can be done but will require more care.
            if contextAdd and isAffine(paramConst):
                paramConstConjunct = paramConst.splitToConjuncts()
                if len(paramConstConjunct) == 1:
                    contextConds.append(paramConst)
        return contextConds

    def extractPolyRepFromFunction(self, comp, scheduleNames, paramNames,
                                   contextConds, levelNo, paramConstraints):
        schedMap = self.createSchedSpace(comp.variables, comp.domain, 
                                         scheduleNames, paramNames, 
                                         contextConds)

        polyDom = PolyDomain(schedMap.domain(), comp)
        id_ = isl.Id.alloc(self.ctx, comp.name, polyDom)
        polyDom.domSet = polyDom.domSet.set_tuple_id(id_)
        self.polyDoms[comp] = polyDom

        self.createPolyPartsFromDefinition(comp, schedMap, levelNo, 
                                           scheduleNames, comp.domain)

    def extractPolyRepFromReduction(self, comp, scheduleNames, paramNames,
                                    contextConds, levelNo, paramConstraints):
        schedMap = self.createSchedSpace(comp.reductionVariables, 
                                         comp.reductionDomain, 
                                         scheduleNames, paramNames, 
                                         contextConds)

        polyDom = PolyDomain(schedMap.domain(), comp)
        id_ = isl.Id.alloc(self.ctx, comp.name, polyDom)
        polyDom.domSet = polyDom.domSet.set_tuple_id(id_)
        self.polyDoms[comp] = polyDom

        self.createPolyPartsFromDefinition(comp, schedMap, levelNo, 
                                           scheduleNames,
                                           comp.reductionDomain)

        domMap = self.createSchedSpace(comp.variables, comp.domain, 
                                       scheduleNames, paramNames, 
                                       contextConds)

        # Initializing the reduction earlier than any other function
        self.createPolyPartsFromDefault(comp, domMap, -1 , scheduleNames)

    def createSchedSpace(self, variables, domains, scheduleNames, paramNames,
                         contextConds):
        # Variable names for refrerring to dimensions
        varNames = [ variables[i].name for i in range(0, len(variables)) ]
        space = isl.Space.create_from_names(self.ctx, in_ = varNames,
                                                     out = scheduleNames,
                                                     params = paramNames)

        schedMap = isl.BasicMap.universe(space)
        # Adding the domain constraints
        [ineqs, eqs] = formatDomainConstraints(domains, varNames)
        schedMap = addConstraints(schedMap, ineqs, eqs)

        # Adding the parameter constraints
        [paramIneqs, paramEqs] = formatConjunctConstraints(contextConds)
        schedMap = addConstraints(schedMap, paramIneqs, paramEqs)

        return schedMap

    def createPolyPartsFromDefinition(self, comp, schedMap, levelNo, 
                                      scheduleNames, domain):
        self.polyParts[comp] = []
        for case in comp.defn:
            schedM = schedMap.copy()

            # The basic schedule is an identity schedule appended with 
            # a level dimension. The level dimension gives the ordering 
            # of the compute objects within a group.

            align, scale = self.defaultAlignAndScale(schedM)

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
                        schedM = addConstraints(schedM, conjunctIneqs, conjunctEqs)
                        parts = self.makePolyParts(schedM, case.expression, None,
                                                   comp, align, scale, levelNo) 
                        for part in parts:
                            self.polyParts[comp].append(part)
                    else:
                        parts = self.makePolyParts(schedM, case.expression, 
                                                   case.condition, comp, align, 
                                                   scale, levelNo)

                        for part in parts:
                            self.polyParts[comp].append(part)
            else:
                assert(isinstance(case, AbstractExpression) or 
                        isinstance(case, Accumulate))
                parts = self.makePolyParts(schedM, case, None, comp, 
                                           align, scale, levelNo)
                for part in parts:
                    self.polyParts[comp].append(part)

        # TODO adding a boundary padding and default to the function 
        # will help DSL usability. 

        # An attempt to subtract all the part domains to find the domain 
        # where the default expression has to be applied. 

        #schedM = isl.BasicMap.identity(self.polyspace)
        #schedM = addConstraints(sched, ineqs, eqs)
        # Adding stage identity constraint
        #levelCoeff = {}
        #levelCoeff[varDims[0]] = -1
        #levelCoeff[('constant', 0)] = compObjs[comp]
        #schedM = addConstraints(schedM, [], [levelCoeff])
        #schedM = addConstraints(schedM, paramIneqs, paramEqs)

        #for part in self.polyParts[comp]:
        #    schedM = schedM.subtract_range(part.schedMap.range())
        #    if (schedM.is_empty()):
        #        break
        #if(not schedM.fast_is_empty()):
        #    bmapList = []
        #    if (isinstance(schedM, isl.BasicMap)):
        #        bmapList.append(schedM)
        #    else:
        #        schedM.foreach_basic_map(bmapList.append)
        #    for bmap in bmapList:    
        #        polyPart = PolyPart(bmap, comp.default, None, comp)
        #        id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)                            
        #        polyPart.schedMap = polyPart.schedMap.set_tuple_id(
        #                                   isl._isl.dim_type.in_, id_)
        #        self.polyParts[comp].append(polyPart)

    def createPolyPartsFromDefault(self, comp, schedMap, levelNo, 
                                   scheduleNames):
        schedM = schedMap.copy()
        align, scale = self.defaultAlignAndScale(sched)

        assert(isinstance(comp.default, AbstractExpression))
        polyPart = PolyPart(schedM, comp.default, None, comp,
                            align, scale, levelNo, False)

        id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)
        polyPart.schedMap = \
                polyPart.schedMap.set_tuple_id(isl._isl.dim_type.in_, id_)
        self.polyParts[comp].append(polyPart)

    def makePolyParts(self, schedMap, expr, pred, comp, 
                      align, scale, levelNo):
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
                isLeftModulo = isAffine(leftExpr, includeModulo = True) and \
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
                    leftCoeff = get_affine_var_and_param_coeff(leftExpr.left)
                    leftCoeff = mapCoeffToDim(leftCoeff)
                    leftConst = get_constant_from_expr(leftExpr.left, affine = True)
                    rightConst = get_constant_from_expr(rightExpr, affine = True)
                    modConst = get_constant_from_expr(leftExpr.right, affine = True)

                    mulName = '_Mul_'
                    remName = '_Rem_'
                    trueSched = schedMap.copy()
                    dimIn = trueSched.dim(isl._isl.dim_type.in_)
                    trueSched = trueSched.insert_dims(isl._isl.dim_type.in_, 
                                                      dimIn, 1)
                    trueSched = trueSched.set_dim_name(isl._isl.dim_type.in_, 
                                                       dimIn, mulName)
                
                    eqs = []
                    leftCoeff[('constant', 0)] = leftConst - rightConst 
                    leftCoeff[('in', dimIn)] = -modConst
                    eqs.append(leftCoeff)

                    trueSched = addConstraints(trueSched, [], eqs)
                    trueSched = trueSched.project_out(isl._isl.dim_type.in_, 
                                                      dimIn, 1)
                    brokenParts.append((trueSched, expr.trueExpression))

                    falseSched = schedMap.copy()
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

                    falseSched = addConstraints(falseSched, ineqs, eqs)
                    falseSched = falseSched.project_out(isl._isl.dim_type.in_,
                                                        dimIn, 2)                    
                    brokenParts.append((falseSched, expr.falseExpression))

        # Note the align and scale lists are cloned otherwise all the 
        # parts will be sharing the same alignment and scaling
        if not brokenParts:
            polyPart = PolyPart(schedMap, expr, pred, comp, list(align), 
                                list(scale), levelNo)
            # Create a user pointer, tuple name and add it to the map
            id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)                            
            polyPart.schedMap = polyPart.schedMap.set_tuple_id(
                                          isl._isl.dim_type.in_, id_)
            polyParts.append(polyPart)
        else:
            for bschedMap, bexpr in brokenParts:
                polyPart = PolyPart(bschedMap, bexpr, pred, comp, list(align), 
                                    list(scale), levelNo)
                # Create a user pointer, tuple name and add it to the map
                id_ = isl.Id.alloc(self.ctx, comp.name, polyPart)                            
                polyPart.schedMap = polyPart.schedMap.set_tuple_id(
                                          isl._isl.dim_type.in_, id_)
                polyParts.append(polyPart)
        return polyParts

    def defaultAlignAndScale(self, sched):
        dimOut = sched.dim(isl._isl.dim_type.out)
        dimIn = sched.dim(isl._isl.dim_type.in_)
        align = {}
        # align[i] = j means input dimension i is mapped to output 
        # dimension j
        for i in range(0, dimIn):
            align[i] = [i+1]
        # the default scaling in each dimension is set to 1 i.e., the
        # schedule dimension correspoinding to input dimension will be 
        # scaled by 1
        scale = [1 for i in range(0, dimIn)]
        return (align, scale)

    def generateCode(self):
        self.polyast = []
        if self.polyParts:
            self.buildAst()

    def buildAst(self):
        #astbld =  isl.AstBuild.from_context(isl.BasicSet("[C, R]->{: R>=1 and C>=1}", self.ctx))
        parts = []
        for plist in self.polyParts.values():
            parts.extend(plist)
       
        # TODO figure out a way to create the correct parameter context
        # since the parameters for all the parts may not be the same
        astbld =  isl.AstBuild.from_context(parts[0].sched.params())
        #astbld =  astbld.set_options(isl.UnionMap("{ }"))

        schedMap = None
        optMap = None
        for part in parts:
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

        # All parts in the group will have the same schedule dimension
        # using the first part as the canonical one
        numIds = parts[0].sched.dim(isl._isl.dim_type.out)
        ids = isl.IdList.alloc(self.ctx, numIds)
        for i in range(0, numIds):
            schedName = parts[0].sched.get_dim_name(isl._isl.dim_type.out, i)
            ids = ids.add(isl.Id.alloc(self.ctx, schedName, None))
        astbld = astbld.set_iterators(ids)

        def printer(arg):
            #print(arg)
            pass

        schedMap.foreach_map(printer)
        self.polyast.append(astbld.ast_from_schedule(schedMap))

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
        depVec = [ '-' for i in range(0, dimOut) ]

        if isinstance(parentPart.comp, Accumulator):
            for i in range(1, dimOut):
                depVec[i] = '*'
            depVec[0] = childPart.levelNo - parentPart.levelNo
            return (depVec, parentPart.levelNo)

        for i in range(0, numArgs):
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
                    #print(ref)
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
                                -get_constant_from_expr(arg, affine=True)
                        accessScale = pscale
                        if depVec[parentVarSchedDim] > 0:
                            depVec[parentVarSchedDim] = \
                                (int(math.ceil(depVec[parentVarSchedDim] *
                                accessScale)))
                        else:         
                            depVec[parentVarSchedDim] = \
                               (int(math.floor(depVec[parentVarSchedDim] *
                                accessScale)))
                            #print(parentPart.sched)
                            #print(childPart.sched)
                            #print(childPart.expr)
                            #print(depVec, ref)
                    else:
                        depVec[parentVarSchedDim] = '*'
                elif len(domDimCoeff) == 0 and len(paramCoeff) > 0:
                    #print(ref)
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
                    accessConstant = get_constant_from_expr(arg, affine = True)
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
                            #print(ref)
                            #assert False
                            #depVec[parentVarSchedDim] = (lowVec, highVec)
                            depVec[parentVarSchedDim] = '*'
                    else:
                        depVec[parentVarSchedDim] = '*'
                else:
                    assert False
            else:
                #print(ref)
                #assert(False)
                depVec[parentVarSchedDim] = '*'

        assert depVec[0] == '-'
        depVec[0] = childPart.levelNo - parentPart.levelNo
        for i in range(0, dimOut):
            if (depVec[i] == '-'):
                depVec[i] = 0
        #for i in range(0, dimOut):
        #    if (depVec[i] == '-'):
        #        depVec[i] = '*'
        #        print(parentPart.sched)
        #        print(childPart.sched)
        #        print(i, ref)
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
        return (depVec, parentPart.levelNo)

    def alignParts(self):
        """ Embeds parts whose dimension is smaller than the schedule space."""
        # Go through the parts in a sorted order and compute alignments
        compObjs = self.stage.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        # Alignments in certain cases may result in changing the relative order of
        # dimensions. This is valid only if there are no dependencies between the
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
            for i in range(0, len(newAlign)):
                if newAlign[i] in currAlign:
                    currAlign.remove(newAlign[i])
            for i in range(0, len(newAlign)):         
                if newAlign[i] == '-':
                    newAlign[i] = currAlign.pop(0)
            return newAlign

        def compatibleAlign(align1, align2):
            compatible = True
            if len(align1) == len(align2):
                for i in range(0, len(align1)):
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
                newAlign = [ '-' for i in range(0, dimIn)]
                for i in range(0, dimIn):
                    for j in range(0, len(align)):
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

            for i in range(0, len(parentGroups)):
                aligns[i] = self.alignWithGroup(p, parentGroups[i])
            
            mergeGroups = []
            # If the alignment has alteast one valid reordering. Add the
            # group to the list of groups to be aligned and merged.
            for i in range(0, len(parentGroups)):
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

    def isParent(self, part1, part2): 
        refs    = part2.getPartRefs()
        objRefs = [ ref.objectRef for ref in refs\
                         if ref.objectRef == part1.comp]
        if len(objRefs) > 0:
            return True
        return False

    def getDomainDimCoeffs(self, sched, arg):
        domDimCoeff = {}
        if (isAffine(arg)):
            coeff = get_affine_var_and_param_coeff(arg)
            for item in coeff:
                if type(item) == Variable:
                    dim = sched.find_dim_by_name(isl._isl.dim_type.in_,
                                                 item.name)
                    domDimCoeff[dim] = coeff[item]
        return domDimCoeff

    def getParamCoeffs(self, sched, arg):
        paramCoeff = {}
        if (isAffine(arg)):
            coeff = get_affine_var_and_param_coeff(arg)
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
        scale = [ '-' for i in range(0, dimIn) ]
        offset = [ '-' for i in range(0, dimIn) ]

        def findDimScheduledTo(part, schedDim):
            for i in range(0, len(part.align)):
                if part.align[i] == schedDim:
                    return i
            return -1

        for ref in parentRefs:
            numArgs = len(ref.arguments)
            for i in range(0, numArgs):
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

        for i in range(0, dimIn):
            if scale[i] == '-':
                scale[i] = 1
            if offset[i] == '-':
                offset[i] = 0
        return (scale, offset)
    
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
                    size = get_constant_from_expr(subsSize)
                else:
                    size = size * get_constant_from_expr(subsSize)
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

        dimSize = interval.upperBound - interval.lowerBound + 1
        return substituteVars(dimSize, paramValMap)

       
    def getVarName(self):
        name = "_i" + str(self._varCount)
        self._varCount+=1
        return name

    def getFuncName(cls):
        name = "_f" + str(self._funcCount)
        self._funcCount+=1
        return name

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

def mapCoeffToDim(coeff):
    variables = list(coeff.keys())
    for var in variables:
        coeffval = coeff[var]
        coeff.pop(var)
        if (isinstance(var, Parameter)):
            coeff[('param', var.name)] = coeffval
        elif (isinstance(var, Variable)):
            coeff[('in', var.name)] = coeffval
    return coeff

def formatScheduleConstraints(dimIn, dimOut, align, scale, levelNo):
    ineqCoeff = []
    eqCoeff   = []
    dimSet = [ False for i in range(0, dimOut) ]
    # Adding identity constraint for each dimension
    for i in range(0, dimIn):
        coeff = {}
        coeff[('out', align[i])] = 1
        assert scale[i] >= 1
        coeff[('in', i)] = -1 * scale[i]
        eqCoeff.append(coeff)
        dimSet[align[i]] = True

    # Setting the leading schedule dimension to level
    levelCoeff = {}
    levelCoeff[('out', 0)] = -1
    levelCoeff[('constant', 0)] = levelNo
    eqCoeff.append(levelCoeff)

    # Setting the remaining dimensions to zero
    for i in range(1, dimOut):
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
    for i in range(0, domLen):
        coeff = {}
        interval = domain[i]

        lbCoeff = get_affine_var_and_param_coeff(interval.lowerBound)
        # Mapping from variable names to the corresponding dimension
        lbCoeff = mapCoeffToDim(lbCoeff)
        lbConst = get_constant_from_expr(interval.lowerBound, affine = True)

        # Normalizing into >= format
        coeff = dict( (n, -lbCoeff.get(n)) for n in  lbCoeff)
        coeff[('constant', 0)] = -lbConst
        coeff[('in', varNames[i])] = 1
        ineqCoeff.append(coeff)

        ubCoeff = get_affine_var_and_param_coeff(interval.upperBound)
        # M_apping from variable names to the corresponding dimension 
        ubCoeff = mapCoeffToDim(ubCoeff)
        ubConst = get_constant_from_expr(interval.upperBound, affine = True)

        # Normalizing into >= format
        coeff = ubCoeff
        coeff[('constant', 0)] = ubConst
        coeff[('in', varNames[i])] = -1
        ineqCoeff.append(coeff)

    return [ineqCoeff, eqCoeff]

def formatConjunctConstraints(conjunct):
    # TODO check if the condition is a conjunction
    ineqCoeff = []
    eqCoeff = []
    for cond in conjunct:
        coeff = {}
        leftCoeff = get_affine_var_and_param_coeff(cond.lhs)
        rightCoeff = get_affine_var_and_param_coeff(cond.rhs)
        leftConst = get_constant_from_expr(cond.lhs, affine = True)
        rightConst = get_constant_from_expr(cond.rhs, affine = True)

        # Mapping from variable names to the corresponding dimension 
        leftCoeff = mapCoeffToDim(leftCoeff)
        rightCoeff = mapCoeffToDim(rightCoeff)
        
        def constantDivFactor(const):
            m = 1
            for coeff in const:
                if isinstance(const[coeff], Fraction):
                    m = (abs(const[coeff].denominator) * m)//gcd(abs(const[coeff].denominator), m)
            assert m.denominator == 1
            m = m.numerator
            return m

        # Normalizing >= format
        if (cond.conditional in ['<=','<']):
            coeff = dict( (n, -leftCoeff.get(n, 0) + rightCoeff.get(n, 0))\
                          for n in set(leftCoeff)| set(rightCoeff) )
            d = constantDivFactor(coeff)
            coeff[('constant', 0)] = -leftConst + rightConst - int(cond.conditional == '<') - Fraction(d-1, d)
            ineqCoeff.append(coeff)
        elif(cond.conditional in ['>=','>']):
            coeff = dict( (n, leftCoeff.get(n, 0) - rightCoeff.get(n, 0))\
                           for n in set(leftCoeff)| set(rightCoeff) )
            d = constantDivFactor(coeff)
            coeff[('constant', 0)] = leftConst - rightConst - int(cond.conditional == '>') + Fraction(d-1, d)
            ineqCoeff.append(coeff)
        else:
            # Weird
            assert(cond.conditional == '==')
            coeff = dict( (n, leftCoeff.get(n, 0) - rightCoeff.get(n, 0))\
                          for n in set(leftCoeff)| set(rightCoeff) )
            coeff[('constant', 0)] = leftConst - rightConst
            eqCoeff.append(coeff)
    return [ineqCoeff, eqCoeff]

def getParamsInvolved(sched, dim):
    numParams = sched.domain().n_param()
    paramNames = [ ]
    constraints = sched.domain().get_constraints()
    dimName = sched.domain().get_dim_name(isl._isl.dim_type.set, dim)
    for p in range(0, numParams): 
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

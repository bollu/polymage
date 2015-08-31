from __future__ import absolute_import, division, print_function

def isCompSelfDependent(self, comp):
    parts = self.poly_parts[comp]
    for p in parts:
        if self.isPartSelfDependent(p):
            return True
    return False

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

def getPartSize(self, part, param_estimates):
    size = None
    domain = part.comp.domain
    if isinstance(part.comp, Accumulator):
        domain = part.comp.reductionDomain
    for interval in domain:
        subsSize = self.getDimSize(interval, param_estimates)
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

def getDimSize(self, interval, param_estimates):
    paramValMap = {}
    for est in param_estimates:
        assert isinstance(est[0], Parameter)
        paramValMap[est[0]] = Value.numericToValue(est[1])

    dimSize = interval.upperBound - interval.lowerBound + 1
    return substituteVars(dimSize, paramValMap)

def getGroupDependenceVectors(self, group, scaleMap = None):
    depVecs = []   
    for part in group:
        refs = part.getPartRefs()
        for ref in refs:
            if ref.objectRef in self.poly_parts:
                for pp in self.poly_parts[ref.objectRef]:
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
                   if ref.objectRef in self.poly_parts ]
    dimAlignMap = {}           
    for ref in parentRefs:
        parentParts = self.poly_parts[ref.objectRef]
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

def isGroupDependentOnPart(self, group, parentPart):
    for part in group:
        refs = part.getPartRefs()
        # This can be more precise
        objRefs = [ ref.objectRef for ref in refs\
                     if ref.objectRef == parentPart.comp]
        if len(objRefs) > 0:
            return True
    return False

def estimateReuse(self, part1, part2):
        """Gives an estimate of data reuse between the two poly parts."""
        # A very naive cost model
        if self.isParent(part1, part2):
            return 1
        return 0

def findParentGroups(self, part, groups):
    parentParts = []
    for comp in self.poly_parts:
        for p in self.poly_parts[comp]:
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

def groupStages(self, param_estimates):
    compObjs = self.stage.orderComputeObjs()
    sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
    # Create a reuse matrix among the poly parts
    parts = []
    for comp in [item[0] for item in sortedCompObjs]:
        parts = parts + self.poly_parts[comp]
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
                       if ref.objectRef in self.poly_parts ]
        dimIn = part.sched.dim(isl._isl.dim_type.in_)
        newScale = [ '-' for i in xrange(0, dimIn)]
        for ref in parentRefs:
            parentParts = self.poly_parts[ref.objectRef]
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

    smallParts, partSizeMap = getSmallComputations(parts, param_estimates)

    # All the parts of a computation which has self dependencies should be
    # in the same group. Bundle such parts together.
    smallGroups = []
    optGroups = self.simpleGroup(param_estimates, single=False)
    #print("intial groups begin")
    #for g in optGroups:
    #    for p in g:
    #        print(p.comp.name)
    #print("intial groups end")
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
            #print(len(children[gi]), isSmall)
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
                        print("parent group begin")
                        for p in optGroups[gi]:
                            print(p.comp.name)
                        print("parent groups end")
                        print("adding groups begin")
                        for g in children[gi]:
                            for p in g:
                                print(p.comp.name)
                        print("adding groups end")
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
                    print("parent group begin")
                    for p in optGroups[gi]:
                        print(p.comp.name)
                    print("parent groups end")
                    print("adding groups begin")
                    for g in children[gi]:
                        for p in g:
                            print(p.comp.name)
                    print("adding groups end")
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
    #print("final groups begin")
    #for g in optGroups:
    #    for p in g:
    #        print(p.comp.name)
    #print("final groups end")
    #print(optGroups)
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

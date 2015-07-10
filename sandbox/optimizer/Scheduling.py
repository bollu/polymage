from Poly import *
def getParentParts(part, group):
     refs = part.getPartRefs()
     parentParts = []
     for ref in refs:
         if ref.objectRef != part.comp:
             if ref.objectRef in group.polyRep.polyParts:
                 parentParts.extend(group.polyRep.polyParts[ref.objectRef])
     return list(set(parentParts))

def baseSchedule(group):
    """
         Construct the base schedule for a group with a polyhedral 
         representation.
    """
    assert(group.isPolyhedral)

    time = {}

    parts = []
    for sublist in group.polyRep.polyParts.values():
        parts.extend(sublist)
    
    for part in parts:
        time[part] = 0

    change = True
    while change:
        change = False
        for part in parts:
            parentParts = getParentParts(part, group)          
            for p in parentParts:
                if time[part] <= time[p]:
                    time[part] = time[p] + 1
                    change = True

    for part in parts:
        part.levelNo = time[part]
        dimIn = part.sched.dim(isl._isl.dim_type.in_)
        dimOut = part.sched.dim(isl._isl.dim_type.out)
        [ineqs, eqs] = formatScheduleConstraints(dimIn, dimOut, 
                                                 part.align, 
                                                 part.scale,
                                                 part.levelNo)
        part.sched = addConstraints(part.sched, ineqs, eqs)

def alignParts():
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
        #print(widths)
                
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
        #    print(g.sched)
        #    print(g.expr)
        #print(stageDeps[i])

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
        maxStage = max([ p.levelNo for p in stageGroups[i] ])
        for p in stageGroups[i]:
            isLiveOut = not isStencil
            #isLiveOut = True
            for gn in xrange(0, len(stageGroups)):
                if gn != i:
                    isLiveOut = isLiveOut or self.isGroupDependentOnPart(
                                                        stageGroups[gn], p)                        
            if p.levelNo == maxStage:        
                p.liveout = True
            p.liveout = p.liveout or isLiveOut     
                
    for gi in stencilGroups:
        assert(len(stageGroups[gi]) > 1)
        hmax = max( [ s.levelNo for s in stageGroups[gi] ] )
        hmin = min( [ s.levelNo for s in stageGroups[gi] ] )
        slopeMin, slopeMax = self.computeTileSlope(stageDeps[gi], hmax)
        #print(slopeMin, slopeMax, hmax - hmin)
        
        #self.splitTile(stageGroups[gi], slopeMin, slopeMax)
        self.overlapTile(stageGroups[gi], slopeMin, slopeMax)
        print(stageDeps[gi])
        print(slopeMin, slopeMax, hmax, len(stageGroups[gi]))
        #for p in stageGroups[gi]:
        #    print(p.scale, p.comp.name + ' = ' +  p.expr.__str__())
        #for p in stageGroups[gi]:
        #    print(p.dimTileInfo)

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
                        baseWidth = h - p.levelNo
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
        #    print(p.sched)
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
    minHeight = min( [ part.levelNo for part in group ] )
    maxHeight = max( [ part.levelNo for part in group ] )
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
                #print(slopeMax, slopeMin, h, L, R, i-1)
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

def simpleSchedule(self, paramEstimates):
    """Generate a simple schedule for the stage."""
    stageGroups = self.baseSchedule(paramEstimates)
    for i in xrange(0, len(stageGroups)):
        for p in stageGroups[i]:
            p.liveout = True

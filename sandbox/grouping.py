from __future__ import absolute_import, division, print_function

def get_group_dep_vecs(group, parts_list=[], scale_map = None):
    dep_vecs = []
    g_poly_rep = group.polyRep
    g_poly_parts = g_poly_rep.poly_parts
    if parts_list == []:
        for comp in g_poly_parts:
            parts_list.extend(g_poly_parts[comp])
    for part in parts_list:
        for ref in part.refs:
            # if the parent is in the same group
            if ref.objectRef in g_poly_parts:
                for pp in g_poly_parts[ref.objectRef]:
                    if pp not in parts_list:
                        continue
                    dep_vec = \
                        part.compute_dependence_vector(pp, ref, scale_map)
                    dep_vecs.append(dep_vec)
    return dep_vecs

def isGroupDependentOnPart(self, group, parentPart):
    for part in group:
        refs = part.refs
        # This can be more precise
        objRefs = [ ref.objectRef for ref in refs\
                     if ref.objectRef == parentPart.comp]
        if len(objRefs) > 0:
            return True
    return False

def findLeafGroups(self, groups):
    leafGroups = []
    for i in range(0, len(groups)):
        isLeaf = True
        for p in groups[i]:
            for j in range(0, len(groups)):
                if j!=i and self.isGroupDependentOnPart(groups[j], p):
                    isLeaf = False
        if isLeaf:
            leafGroups.append(groups[i])
    return leafGroups

def auto_group(pipeline, param_estimates):
    comps = pipeline._comp_objs
    sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
    # Create a reuse matrix among the poly parts
    parts = []
    for comp in [item[0] for item in sortedCompObjs]:
        parts = parts + self.poly_parts[comp]

    def get_group_cost(group):
        return 1

    # Progressive algorithm for grouping stages assigns groups level by 
    # level trying to maximze reuse and size of stencil groups.

    # Filter out small computations. We characterize a computation as 
    # small when the domains of parents and the computation itself is 
    # small to benefit from tiling or parallelization. The parents are
    # included so that we do not miss out on storage optimzations.
    def get_small_computations(parts, estimates, size_thresh):
        small_parts = []
        # This can be more precise but for now just estimating the 
        # size of the part by the size of the computation object.

        # Currentlty the size is just being estimated by the number
        # of points in the domain. It can be more accurately done 
        # by considering the arithmetic intensity in the expressions.
        part_size_map = {}
        for p in parts:
            part_size_map[p] = p.get_size(estimates)

        for i in range(0, len(parts)):
            is_small_part = False
            if part_size_map[parts[i]] != '*':
                is_small_part = (part_size_map[parts[i]] <= size_threshold)
            for j in range(i+1, len(parts)):
                if parts[j].is_parent_of(parts[i]):
                    if part_size_map[parts[j]] != '*':
                        is_small_part = is_small_part and \
                            (part_size_map[parts[j]] > size_threshold)
                    else:
                        is_small_part = False
            if is_small_part:
                small_parts.append(parts[i])

        return smallParts, part_size_map

    small_parts, part_size_map = \
        get_small_computations(parts, param_estimates, size_threshold)

    # All the parts of a computation which has self dependencies should be
    # in the same group. Bundle such parts together.
    small_groups = []
    optGroups = self.simpleGroup(param_estimates, single=False)

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
        for gi in range(0, len(optGroups)):
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
    #bundles = [ b for b in bundles if b not in smallGroups ]
    """
    def estimateGroupCostWithBundle(group, bundle, scale, partSizeMap):
        return 1

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

        for i in range(0, len(parentGroups)):
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
        for i in range(0, len(parentGroups)):
            currCost = getGroupCost(parentGroups[i])
            newCost = estimateGroupCostWithBundle(parentGroups[i], b,
                                                  scales[i], partSizeMap)
            if isProfitable(currCost, newCost):
                candGroups.append(i)
        
        mergedGroup = []
        for i in candGroups:
            for j in range(0, len(b)):
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

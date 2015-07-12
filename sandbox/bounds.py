from __future__ import absolute_import, division, print_function

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


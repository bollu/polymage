from __future__ import absolute_import, division, print_function

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
            childRefs = childPart.refs
            if childPart.pred:
                childRefs += childPart.pred.collect(Reference)
            childRefs = [ ref for ref in childRefs if ref.objectRef == parentFunc]

            deps = []
            for ref in childRefs:
                deps += extract_value_dependence(childPart, ref,
                            parentStage.polyRep.poly_doms[parentFunc])
            
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


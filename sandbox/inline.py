from __future__ import absolute_import, division, print_function

from constructs import *
from poly import *

def inline_piecewise(child_group, parent_group, no_split = False):
    ref_to_inline_expr_map = {}
    # Inling currently only handles non-fused stages
    assert (not child_group.is_fused() and not parent_group.is_fused())
    # Computation object in the parent group can only be a function
    parent_comp = parent_group.compute_objs[0]
    assert isinstance(parent_comp, Function)
    child_comp = child_group.compute_objs[0]

    # Simple scalar functions which are defined on a non-bounded 
    # integer domain can be inlined.
    # TODO

    # Inlining only if both child and parent stage have a polyhedral
    # representation
    if child_group.polyRep.poly_parts and parent_group.polyRep.poly_parts:
        child_parts = child_group.polyRep.poly_parts
        parent_parts = parent_group.polyRep.poly_parts
        parent_doms = parent_group.polyRep.poly_doms
        for child_part in child_parts[child_comp]:
            # - Collect all the refs and pick out refs to only parent_comp
            child_refs = child_part.refs
            child_refs = [ ref for ref in child_refs \
                                 if ref.objectRef == parent_comp ]
            # Compute dependence relations between child and parent
            deps = []
            for ref in child_refs:
                deps += extract_value_dependence(child_part, ref,
                            parent_doms[parent_comp])
            # Check if all the values come from the same parent part
            dep_to_part_map = {}
            for dep in deps:
                # total accessible region by this dependence
                access_region = dep.rel.range().copy().reset_tuple_id()
                # holds the remaining available region after each piece's
                # region is cut off from the access region
                left_over = dep.rel.range().copy().reset_tuple_id()

                # analyze for each parent piece
                for parent_part in parent_parts[parent_comp]:
                    # access region of this part
                    part_region = \
                        parent_part.sched.domain().copy().reset_tuple_id()
                    # portion of the access_region touched by this part's
                    # region
                    partdiff = access_region.subtract(part_region)
                    # chop off the part's region from the left-over region
                    left_over = left_over.subtract(part_region)
                    if(partdiff.is_empty()):
                        dep_to_part_map[dep] = parent_part
                if (not left_over.is_empty()):
                    assert False, "Inlining cannot be done."

            parts = list(set(dep_to_part_map.values()))
            single_part = (len(parts) == 1)

            if(single_part):
                parent_expr = parts[0].expr
                if parts[0].pred:
                    inline = False
                else:
                    inline_deps = []
                    for ref in child_refs:
                        ref_to_inline_expr_map[ref] = parent_expr
            elif(no_split):
                pass
            else:
                pass
    else:
        pass
    return ref_to_inline_expr_map


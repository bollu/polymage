from __future__ import absolute_import, division, print_function

from expression import isAffine
from poly import extract_value_dependence
import logging

# LOG CONFIG #
bounds_logger = logging.getLogger("bounds.py")
bounds_logger.setLevel(logging.ERROR)

LOG = bounds_logger.log

def check_refs(child_group, parent_group):
    # Check refs works only on non-fused groups. It can be made to
    # work with fused groups as well. However, it might serve very
    # little use.
    # assert (not child_group.isFused() and not parent_group.isFused())

    # get the only comp_obj in child and parent groups
    parent_func = parent_group._comp_objs[0]
    child_obj = child_group._comp_objs[0]

    # Only verifying if both child and  parent group have a polyhedral
    # representation
    if child_group.polyRep.poly_parts and parent_group.polyRep.poly_doms:
        for child_part in child_group.polyRep.poly_parts[child_obj]:
            # Compute dependence relations between child and parent
            child_refs = child_part.refs
            if child_part.pred:
                child_refs += child_part.pred.collect(Reference)

            # It is not generally feasible to check the validity of
            # and access when the reference is not affine. 
            # Approximations can be done but for now skipping them.
            def affine_ref(ref):
                affine = True
                for arg in ref.arguments:
                    affine = affine and isAffine(arg)
                return affine

            # filter out only the affine refs to parent_func
            child_refs = [ ref for ref in child_refs \
                                 if ref.objectRef == parent_func and
                                    affine_ref(ref) ]

            deps = []
            parent_dom = parent_group.polyRep.poly_doms[parent_func]
            for ref in child_refs:
                deps += extract_value_dependence(child_part, ref, parent_dom)
            for dep in deps:
                diff = dep.rel.range().subtract(parent_dom.dom_set)
                # ***
                ref_str = "referenced    = "+str(dep.rel.range())
                dom_str = "parent domain = "+str(parent_dom.dom_set)
                log_level = logging.DEBUG
                LOG(log_level, ref_str)
                LOG(log_level, dom_str)
                # ***
                if(not diff.is_empty()):
                    # ***
                    log_level = logging.ERROR
                    LOG(log_level, "_______________________")
                    LOG(log_level, "Reference out of domain")
                    LOG(log_level, ref_str)
                    LOG(log_level, dom_str)
                    LOG(log_level, "_______________________")
                    # ***
                    raise TypeError("Reference out of domain", child_group,
                                     parent_group, diff)

    return

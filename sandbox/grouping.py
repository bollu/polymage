from __future__ import absolute_import, division, print_function

from constructs import *

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

def auto_group(pipeline):
    param_est = pipeline._param_estimates
    size_thresh = pipeline._size_threshold
    grp_size = pipeline._group_size

    comps = pipeline._comp_objs
    group_map = pipeline._groups
    g_child_par_map = pipeline._group_parents
    g_par_child_map = pipeline._group_children
    groups = list(set([group_map[comp] for comp in group_map]))

    def get_group_cost(group):
        return 1

    # Progressive algorithm for grouping stages assigns groups level by level
    # trying to maximze reuse and size of stencil groups.

    # Filter out small computations. We characterize a computation as small
    # when the domains of parents and the computation itself is too small to
    # benefit from tiling or parallelization. The parents are included so that
    # we do not miss out on storage optimzations.
    def get_small_comps(pipeline, comps):
        comp_parents = pipeline._comp_objs_parents

        small_comps = []
        # Currentlty the size is just being estimated by the number of points
        # in the domain. It can be more accurately done by considering the
        # arithmetic intensity in the expressions.
        comp_size_map = {}
        for comp in comps:
            parts = group_map[comp].polyRep.poly_parts[comp]
            p_sizes = []
            for p in parts:
                p_size = p.get_size(param_est)
                if p_size == '*':
                    p_sizes.append(0)
                else:
                    p_sizes.append(p_size)
            comp_size_map[comp] = sum(p_sizes)

        for comp in comps:
            is_small_comp = False
            if comp_size_map[comp] != '*':
                is_small_comp = (comp_size_map[comp] <= size_thresh)
            # iterate over parents of comp
            for pcomp in comp_parents:
                if comp_size_map[pcomp] != '*':
                    is_small_comp = is_small_comp and \
                        (comp_size_map[pcomp] > size_thresh)
                else:
                    is_small_comp = False
            if is_small_comp:
                small_comps.append(comp)

        return small_comps, comp_size_map


    small_comps, comp_size_map = get_small_comps(pipeline, comps)

    # loop termination boolean
    opt = True
    while opt:
        opt = False

        for group in pipeline._group_children:
            is_small_grp = True
            is_reduction_grp = False
            for comp in group._comp_objs:
                if not comp in small_comps:
                    is_small_grp = False
                if isinstance(comp, Reduction):
                    is_reduction_grp = True
            for g_child in pipeline._group_children[group]:
                for comp in g_child._comp_objs:
                    if isinstance(comp, Reduction):
                        is_reduction_grp = True
            # merge if
            # 1. big enough
            # 2. does not contain reduction
            # 3. number of comps in group < grp_size
            if not is_small_grp and \
               not is_reduction_grp and \
               len(group._comp_objs) < grp_size:
                merge = True

                all_children_of_g = pipeline._group_children[group]

                # check if group and its children are merged the total
                # group size exceeds grp_size
                child_comps_count = 0
                for g_child in all_children_of_g:
                    child_comps_count += len(g_child._comp_objs)

                merge_count = len(group._comp_objs)+child_comps_count
                if merge_count > grp_size:
                    merge = False

                # - if group has many children
                if (len(all_children_of_g) > 1):
                    if merge:
                        # check if its possible to group with all children.
                        # collect all the parents of all children of group
                        all_parents = []
                        for g_child in all_children_of_g:
                            all_parents += pipeline._group_parents[g_child]
                        all_parents = list(set(all_parents))
                        all_parents.remove(group)

                        # if all_parents are children of group => OK to merge
                        if not set(all_parents).issubset(set(all_children_of_g)):
                            merge = False

                    if merge:
                        new_grp = group
                        for g_child in all_children_of_g:
                            new_grp = pipeline.merge_groups(new_grp, g_child)
                        opt = True
                        break
                # - if group has only one child
                elif (len(all_children_of_g) == 1):
                    if merge:
                        pipeline.merge_groups(group, all_children_of_g[0])
                        opt = True
                        break

    return

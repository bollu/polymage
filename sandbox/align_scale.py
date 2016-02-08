from __future__ import absolute_import, division, print_function

import logging
from utils import *
from poly import *

# LOG CONFIG #
align_scale_logger = logging.getLogger("align_scale.py")
align_scale_logger.setLevel(logging.INFO)
LOG = align_scale_logger.log

class ASPacket(object):
    """
    contains alignment and scaling for a part
    """
    def __init__(self, _align, _scale):
        self.align = _align
        self.scale = _scale

    def clone(self):
        align = list(self.align)
        scale = list(self.scale)
        return ASPacket(align, scale)

class ASInfo(object):
    """
    contains the following types of scaling and alignment for a part:
    1. true : contains only those that can be inferred from the parent part
    2. full : [true] + random alignment and scaling for the remaining dims
    """
    def __init__(self, _true, _full, _borrowed=False):
        # _true, _full -> type : ASPacket
        # absolute alignment and scaling w.r.t base part
        assert type(_true) == ASPacket and type(_full) == ASPacket
        self.true = _true
        self.full = _full
        self.borrowed = bool(_borrowed)

    def clone(self):
        true = self.true.clone()
        full = self.full.clone()
        return ASInfo(true, full, self.borrowed)

def default_align_and_scale(sched, max_dim=None, shift=False):
    dim_in = sched.dim(isl._isl.dim_type.in_)

    # align[i] = j means input dimension i is mapped to output
    # dimension j

    # the default scaling in each dimension is set to 1 i.e., the
    # schedule dimension correspoinding to input dimension will be
    # scaled by 1

    #shift=False
    if shift:
        assert max_dim != None
        off = max_dim - dim_in
        align = [ i for i in range(off, max_dim) ]
    else:
        align = [ i for i in range(0, dim_in) ]

    scale = [1 for i in range(0, dim_in)]

    if max_dim:
        if dim_in < max_dim:
            for i in range(dim_in, max_dim):
                align.append(NULL)
                scale.append(NULL)

    return (align, scale)

def align_and_scale(pipeline, group):

    """
    Alignment structure:
    [alignment]

    [alignment] of a dimension 'i', of a polypart is given by the i'th member
    of the list part.align. This value is a dimension in the root part or the
    part that is at the first level of the topologically sorted group. In case
    there are multiple parent parts, the same alignment should hold for all
    those parts. A mapping to NULL instead of an integer dimension indicates
    alignment to none.

    Scaling factor structure:
    [scale factors]

    [scale factors] determine the amount by which the dimension has to be
    scaled to uniformize the dependencies. Each reference to the parent is
    considered while determing the scaling factors. All the references
    should have the same scaling factor in a particular dimension otherwise
    the scaling factor for the dimension cannot be determined uniquely.

    * The part.align does not contain topological level information in it,
    i.e, the dimensions handled here are spatial.
    """

    # ***
    log_level = logging.DEBUG
    LOG(log_level, "___________________________________")
    LOG(log_level, "in align_and_scale_parts()")
    # ***

    def log_a_s(level, al_str, sc_str, a_s_pack):
        log_str1 = al_str+" = "+str([i for i in a_s_pack.align])
        log_str2 = sc_str+" = "+str([i for i in a_s_pack.scale])
        LOG(level, log_str1)
        LOG(level, log_str2)
        return

    def non_null(_list, null=NULL):
        n_list = [x for x in _list if x != null]
        return n_list

    def new_list(fill=NULL):
        ''' return a new empty list of length max_dim '''
        return [fill for i in range(0, max_dim)]

    def extract_arg_vars_coefs(ref_args):
        '''
        Extracts variables and their co-efficients from each argument
        in ref_args
        '''
        ref_arg_vars = []
        ref_arg_coeffs = {}
        for i in range(0, len(ref_args)):
            arg = ref_args[i]
            if isAffine(arg):
                arg_vars = arg.collect(Variable)
                nvars = len(arg_vars)
                if nvars == 0:
                    # example: dim 0 in (2, x+1, y-1)
                    ref_arg_vars.append(NULL)
                    pass
                elif nvars == 1:
                    # example: dims 0, 1 and 2 in (c, x+1, y-1)
                    ref_arg_vars.append(arg_vars[0])
                    arg_coeff = get_affine_var_and_param_coeff(arg)
                    ref_arg_coeffs.update(arg_coeff)
                # Restricting the arg expression to have only one variable
                else:
                    # example: a) dim 1 in (c, x+x, y-1)
                    #          b) dim 2 in (c, x+1, x-y)
                    # asserting the type to be not (b), a rare corner case,
                    # which we dont want to handle now
                    assert(arg_vars[0] == arg_vars[i] \
                           for i in range(1, len(arg_vars)))
                    ref_arg_vars.append(arg_vars[0])
                    arg_coeff = get_affine_var_and_param_coeff(arg)
                    ref_arg_coeffs.update(arg_coeff)

        return ref_arg_vars, ref_arg_coeffs

    def pick_best(aln_scl_list):
        '''
        Picks the align_scale which contains maximum information.
        aln_scl_list is a list of ASInfo objects, with info across multiple
        parents/children.
        '''
        assert aln_scl_list
        dim_max = 0
        best = None
        # info with minimum number of unknowns is better
        for aln_scl in aln_scl_list:
            if not aln_scl:  # empty entry
                continue
            dims = non_null(aln_scl.true.align)
            if dim_max < len(dims):
                dim_max = len(dims)
                best = aln_scl

        if best == None:
            best = aln_scl_list[0]

        return best

    def prune_align_scale(aln_scl, dim_in, info):
        '''
        prune the align_scale information to dim_in from max_dim
        '''
        max_dim = info.max_dim
        null = [NULL for d in range(0, max_dim-dim_in)]
        aln_scl.true.align[dim_in:max_dim] = null
        aln_scl.full.align[dim_in:max_dim] = null
        aln_scl.true.scale[dim_in:max_dim] = null
        aln_scl.full.scale[dim_in:max_dim] = null

        return

    def complete_align_scale(part, align, scale):
        '''
        extend the align_scale to a list of length max_dim using special NULLs
        '''
        dim_in = part.sched.dim(isl._isl.dim_type.in_)
        assert dim_in > 0
        new_align = align
        new_scale = scale
        if not all(align[dim] != NULL for dim in range(0, dim_in)):
            used_dims = [dim for dim in align if dim != NULL]
            universal = set(range(0, len(align)))
            avail_dims = \
                list(universal.difference(set(used_dims)))
            avail_dims.sort(reverse=True)
            for dim in range(0, dim_in):
                if align[dim] == NULL:
                    new_align[dim] = avail_dims.pop()
                    new_scale[dim] = 1
        return (new_align, new_scale)

    def borrow_align_scale(comp_parts, info):
        '''
        borrows align_scale information for a part with NULL true align_scale,
        from another part with maximum true align_scale info
        '''

        if len(comp_parts) > 1:
            # list all parts of the comp that has NULL true alignemnt
            empty = [p for p in comp_parts \
                         if not non_null(info.align_scale[p].true.align)]

            # pick the best align_scale from the parts with non NULL true
            # alignment
            non_empty = list(set(comp_parts).difference(set(empty)))
            aln_scl_list = [info.align_scale[p] for p in non_empty]
            best_aln_scl = pick_best(aln_scl_list)

            # set the 'borrowed' flag, and borrow the alignment and scaling
            # info from the part picked as best. prune the info to the part's
            # actual dimensionality.
            best_aln_scl.borrowed = True
            for p in empty:
                dim_in = p.sched.dim(isl._isl.dim_type.in_)
                best_clone = best_aln_scl.clone()
                prune_align_scale(best_clone, dim_in, info)
                info.align_scale[p] = best_clone

        return

    def align_scale_vars(child_part, parent_part,
                         ref_arg_vars, ref_arg_coeffs,
                         info, reverse=False):
        '''
        Finds an alignment and scaling factor for each dimension associated
        with the reference argument variable

        reverse flag is used to set the alignment direction:
        align parent to child : reverse = True
        align child to parent : reverse = False
        '''

        def get_domain_dims(sched, var_list):
            '''
            Gets the dimension associated with each of the varibale in the
            var_list, assuming that it is present in the schedule
            '''
            # dict: var -> dim
            domain_dims = {}
            for var in var_list:
                dim = sched.find_dim_by_name(isl._isl.dim_type.in_, var.name)
                domain_dims[var] = dim
            return domain_dims

        def get_argvar_order(var_list):
            '''
            Assumes that the var_list was built in the reference order of the
            argument. Each var relatively aligns to the dimension of the
            parent, that is same as its occurence position in the argument.
            '''
            dim_map = {}
            dim = 0
            for var in var_list:
                if var != NULL:
                    dim_map[var] = dim
                dim += 1
            assert(dim <= max_dim)
            return dim_map

        def new_dict():
            ''' return a new empty dictionary of length max_dim '''
            _dict = {}
            for i in range(0, max_dim):
                _dict[i] = NULL
            return _dict

        def compute_abs(dim, dst_pack, rel_align_map, rel_scale_map, \
                        reverse=False, pseudo=False):
            ''' compute the absolute align_scale using relative align_scale '''
            root_dim = rel_align_map[dim]
            dim_align = NULL
            dim_scale = NULL
            if root_dim != NULL:
                dim_align = dst_pack.align[root_dim]
                if dim_align != NULL:
                    dim_scale = 1
                    if not pseudo:
                        if reverse:
                            dim_scale = Fraction(dst_pack.scale[root_dim]) / \
                                        Fraction(rel_scale_map[dim])
                        else:
                            dim_scale = dst_pack.scale[root_dim] * \
                                        rel_scale_map[dim]

            return dim_align, dim_scale

        def unit_test(align, scale):
            # 1. test for unique alignment
            n_align = non_null(align)
            assert (len(n_align) == len(set(n_align)))

            # 2. test for alignment boundary
            out_of_bound = all(max_dim > dim for dim in n_align)\
                           and \
                           all(0 <= dim for dim in n_align)
            assert (out_of_bound)

            # 3. test for scaling
            for dim in range(0, max_dim):
                if align[dim] == NULL:
                    assert scale[dim] == NULL
                else:
                    assert scale[dim] != NULL
                    # precisely, this should be tested for a constant

            # 4. dimensionality test
            assert len(align) == info.max_dim
            assert len(scale) == info.max_dim

            return

        def align_src_to_dst(dst, \
                             rel_align, rel_scale, \
                             src_map, dst_map, \
                             ref_arg_coeffs, reverse=False):
            '''
            Given a target and a relative alignment-scaling info, find the
            true and full alignment-scaling
            '''

            # for each variable in the reference argument, get the alignment of
            # the part relative to the destination
            # dim:dim map from src to dst
            for var in ref_arg_vars:
                if var != NULL:
                    rel_align[src_map[var]] = dst_map[var]
                    rel_scale[src_map[var]] = ref_arg_coeffs[var]

            # only those dims whose align,scale can be infered from the parent
            true_align = new_list()
            true_scale = new_list()
            for dim in range(0, max_dim):
                true_align[dim], true_scale[dim] = \
                    compute_abs(dim, dst.true, rel_align, rel_scale, reverse)

            # dims of the src part which are not related to dst
            rem_src_dims = [dim for dim in src_map.values() \
                                  if rel_align[dim] == NULL]

            # dims of the dst part which are not related to src
            rem_dst_dims = [dim for dim in dst_map.values() \
                                  if dim not in rel_align.values()]
            rem_dst_dims.sort(reverse=True)

            # align each of the rem_src_map to any of rem_dst_map
            for dim in rem_src_dims:
                rel_scale[dim] = 1
                if rem_dst_dims:
                    rel_align[dim] = rem_dst_dims.pop()
                else:
                    rel_align[dim] = '*'  # dangling

            # aligned to one of the dst dims
            aligned_dims = [dim for dim in src_map.values()
                                  if rel_align[dim] != '*']
            # not aligned to any of the dst dims
            dangling_dims = [dim for dim in src_map.values()
                                   if dim not in aligned_dims]

            # normalize to the base align_scale using the relative align_scale
            full_align = new_list()
            full_scale = new_list()
            for dim in aligned_dims:
                if true_align[dim] != NULL:
                    full_align[dim] = true_align[dim]
                    full_scale[dim] = true_scale[dim]
                else:
                    full_align[dim], full_scale[dim] = \
                        compute_abs(dim, dst.full, rel_align, rel_scale, \
                                    reverse, pseudo=True)

            # dangling_dims are assigned any available dim of the base alignment
            avail_dims = [dim for dim in range(0, max_dim) \
                                if dim not in full_align]
            avail_dims.sort(reverse=True)
            for dim in dangling_dims:
                full_align[dim] = avail_dims.pop()
                full_scale[dim] = 1

            # validity tests
            unit_test(true_align, true_scale)
            unit_test(full_align, full_scale)

            true = ASPacket(true_align, true_scale)
            full = ASPacket(full_align, full_scale)
            part_solution = ASInfo(true, full)

            return part_solution

        # BEGIN
        max_dim = info.max_dim

        # relative alignment and scaling
        # dict: dim -> dim
        rel_align = new_dict()
        rel_scale = new_dict()

        # var:dim map of variable domain of the child part
        child_map = get_domain_dims(child_part.sched, \
                                    child_part.comp.func.variableDomain[0])
        # var:dim map of parent part in argument order
        parent_map = get_argvar_order(ref_arg_vars)

        # set the source and the destination for alignment/scaling
        if reverse == True:
            dst_part = child_part
            src_part = parent_part

            dst_map = child_map
            src_map = parent_map
        else:
            dst_part = parent_part
            src_part = child_part

            dst_map = parent_map
            src_map = child_map

        # child info
        dst_pack = info.align_scale[dst_part]
        assert(dst_pack)

        # align and scale source to destination
        src_part_solution = align_src_to_dst(dst_pack,
                                             rel_align, rel_scale,
                                             src_map, dst_map,
                                             ref_arg_coeffs, reverse)

        if src_part in info.align_scale:  # already has a solution
            src_part_solution = \
                pick_best([info.align_scale[src_part], src_part_solution])
        info.align_scale[src_part] = src_part_solution  # add to global set

        return

    def solve_from_ref_part(part, ref, ref_part, info, reverse=False):
        '''
        given the Reference object's poly part, solve for alignment and scaling
        reverse = True  : solve for parent
        reverse = False : solve for child
        '''
        # process the argument list
        ref_args = ref.arguments
        ref_arg_vars, ref_arg_coeffs = extract_arg_vars_coefs(ref_args)

        # match the part variables with the reference variables
        align_scale_vars(part, ref_part, \
                         ref_arg_vars, ref_arg_coeffs, \
                         info, reverse)

        return

    def solve_from_ref(part, ref, info, reverse=False):
        '''
        given the Reference object, solve for alignment and scaling
        reverse = True  : solve for parent
        reverse = False : solve for child
        '''
        max_dim = info.max_dim
        ref_comp = info.pipe.func_map[ref.objectRef]

        ref_poly_parts = group.polyRep.poly_parts[ref_comp]
        for ref_part in ref_poly_parts:
            solve_from_ref_part(part, ref, ref_part, info, reverse)

        return

    ''' bottom-up phase '''

    def solve_comp_from_children(comp, info):
        """
        solve for alignment and scaling of comp by referring to the solved
        child comps.
        """

        # if already visited
        if comp in info.solved:
            return

        part_map = info.group.polyRep.poly_parts
        all_children = info.pipe.comp.children
        solved_children = [child for child in all_children \
                                   if child in info.solved]

        for child in solved_children:
            for child_part in part_map[child]:
                # collect the references made only to comp
                refs = [ref for ref in child_part.refs
                              if ref.objectRef == comp.func]

                # if the poly part makes no reference to any other compute object
                if not refs:
                    continue
                else:
                    # compute the alignment and scaling across references
                    for ref in refs:
                        solve_from_ref(child_part, ref, info, reverse=True)

        # check if all parts of the comp have an alignment and scaling
        for part in part_map[comp]:
            assert part in info.align_scale

        # borrow align_scale for the under previleged!
        borrow_align_scale(part_map[comp], info)

        return

    def solve_comp_parents(parents, info):
        '''
        Solves for alignment and scaling for poly parts of comps in the parents
        list, using the information of only the already solved children.
        Traverses in a breadth-first fashion.
        '''
        if parents == []:
            return

        poly_parts = info.group.polyRep.poly_parts
        parent_order = {}
        for parent in parents:
            parent_order[parent] = poly_parts[parent][0].level
            # index 0 picks the fist poly part

        # parents near leaf level shall be solved at the earliest
        sorted_order = sorted(parent_order.items(), key=lambda x:x[1], \
                              reverse=True)
        sorted_parents = [x[0] for x in sorted_order]

        for parent in sorted_parents:
            solve_comp_from_children(parent, info)
            info.solved.append(parent)

        # collect unsolved parents of parents within the group
        all_grand_parents = []
        for parent in sorted_parents:
            if parent.parents:
                grand_parents = [gp for gp in parent.parents \
                                      if gp not in info.solved and \
                                         gp.group == info.group]
                all_grand_parents += grand_parents
        all_grand_parents = list(set(all_grand_parents))

        solve_comp_parents(all_grand_parents, info)

        return

    ''' top-down phase '''

    def solve_comp_from_parents(comp, info):
        """
        solve for alignment and scaling of comp by referring to the solved
        parent comps.
        """
        # if already visited
        if comp in info.solved:
            return

        func_map = info.pipe.func_map

        comp_parts = info.group.polyRep.poly_parts[comp]
        all_parents = [p for p in comp.parents \
                           if p.group == info.group]
        solved_pars = [par for par in all_parents \
                             if par in info.solved]

        # unsolved parents are the newly discovered parents
        discovered_parents = set(all_parents).difference(set(solved_pars))
        info.discovered = list(set(info.discovered).union(discovered_parents))

        # solve for each part
        for p in comp_parts:
            # collect the references to solved parents
            refs = [ref for ref in p.refs
                          if func_map[ref.objectRef] in solved_pars]

            # if the poly part makes no reference to any other compute object
            if not refs:
                aln_scl = default_align_and_scale(p.sched, info.max_dim, shift=True)
                true_default = ASPacket(new_list(), new_list())
                full_default = ASPacket(aln_scl[0], aln_scl[1])
                info.align_scale[p] = ASInfo(true_default, full_default)
            else:
                # compute the alignment and scaling across references
                for ref in refs:
                    solve_from_ref(p, ref, info)

            assert p in info.align_scale

        # borrow align_scale for the under previleged!
        borrow_align_scale(comp_parts, info)

        return

    def solve_comp_children(children, info):
        '''
        Solves for alignment and scaling for poly parts of comps in the
        children list, using the information of only the already solved
        parents. Traverses in a breadth-first fashion.
        '''
        # leaf comp node
        if children == []:
            return

        group = info.group
        poly_parts = group.polyRep.poly_parts
        sorted_children = group.get_sorted_comps()

        for child in sorted_children:
            solve_comp_from_parents(child, info)
            info.solved.append(child)

        # collect unsolved children of children within the group
        all_grand_children = []
        for child in sorted_children:
            if child.children:
                grand_children = [gc for gc in child.children \
                                       if gc not in info.solved and \
                                          gc in info.group.comps]
                all_grand_children += grand_children
        all_grand_children = list(set(all_grand_children))

        solve_comp_children(all_grand_children, info)

        return

    def plus_one(align):
        plus_align = []
        for dim in align:
            if dim != NULL:
                plus_align.append(dim+1)
            else:
                plus_align.append(NULL)
        return plus_align

    def increment_aligns(comps, part_map):
        for comp in comps:
            for part in part_map[comp]:
                plus_align = plus_one(part.align)
                part.set_align(plus_align)
        return

    def find_scale_norm(info):
        '''
        find the lcm of scales of each dim, of all aligned parts
        '''
        norm = [1 for i in range(0, max_dim)]

        # compute the lcm of the Fraction denominators of scaling factors of
        # all poly parts in the group, for each dimension
        for part in info.align_scale:
            align = part.align
            scale = part.scale
            for dim in range(0, max_dim):
                if align[dim] != NULL:
                    d = Fraction(scale[dim].denominator)
                    align_dim = align[dim]
                    norm[align_dim] = lcm(d, norm[align_dim])

        return norm

    def normalize_scale(norm, info):
        '''
        normalize the scaling factors, so that none of them is lesser than 1
        '''
        for part in info.align_scale:
            align = part.align
            scale = part.scale
            new_scale = [1 for i in range(0, max_dim)]
            for dim in range(0, max_dim):
                if align[dim] != NULL:
                    align_dim = align[dim]
                    new_scale[dim] = int(norm[align_dim] * part.scale[dim])
                else:
                    new_scale[dim] = NULL
            part.set_scale(new_scale)
        return

    ''' main '''
    comps = group.get_sorted_comps()
    group_part_map = group.polyRep.poly_parts

    # list all parts with no self references and find the max dim
    max_dim = 0
    no_self_dep_parts = []
    for comp in comps:
        for p in group_part_map[comp]:
            p_dim_in = p.sched.dim(isl._isl.dim_type.in_)
            if not p.is_self_dependent:
                no_self_dep_parts.append(p)
                if max_dim < p_dim_in:
                    max_dim = p_dim_in

    # begin from the topologically earliest comp parts as the base parts for
    # scaling and alignment reference
    root_comps = group.root_comps

    # prefer the comp closer to the root_comp of the pipeline
    abs_min_level = 1000000
    for comp in root_comps:
        if abs_min_level > comp.level:
            abs_min_level = comp.level

    abs_min_comps = []
    abs_min_parts = []
    for comp in root_comps:
        if comp.level == abs_min_level:
            abs_min_comps.append(comp)
            abs_min_parts += group_part_map[comp]

    # pick the min-level part with the highest dimensionality as base part
    base_part = None
    dim_max = 0
    for p in abs_min_parts:
        p_dim_in = p.sched.dim(isl._isl.dim_type.in_)
        if p_dim_in > dim_max:
            dim_max = p_dim_in
            base_part = p

    assert (base_part != None)
    base_comp = base_part.comp

    # initial alignment and scaling for the base comp parts
    # ***
    log_level = logging.DEBUG
    LOG(log_level, "____")
    LOG(log_level, str(base_comp.func.name)+\
                   " (level : "+str(base_part.level)+")")

    class Info(object):
        def __init__(self, _pipe, _group, _max_dim):
            self.pipe = _pipe
            self.group = _group
            self.max_dim = _max_dim  # max dimensionality of the group
            self.solved = []
            # list of comps encountered while solving for base comp family
            self.discovered = []
            # mapping from poly-part to ASPacket
            self.align_scale = {}

    info = Info(pipeline, group, max_dim)
    info.solved.append(base_comp)
    for part in group_part_map[base_comp]:
        # set default values for base parts
        align, scale = default_align_and_scale(part.sched, max_dim, shift=True)
        # update to the temporary info
        true_pack = ASPacket(align, scale)
        full_pack = ASPacket(align, scale)
        info.align_scale[part] = ASInfo(true_pack, full_pack)

    # recursively compute alignment and scaling for the family of base_comp
    if base_comp.children:
        solve_comp_children(base_comp.children, info)  # <-

    # compute newly discovered parents iteratively until no new parent is
    # discovered.
    solve_comp_parents(info.discovered, info)  # <-

    # set all the solutions into polypart object members
    for comp in comps:
        for part in group_part_map[comp]:
            # after solving, no poly part cannot not have a solution
            assert part in info.align_scale
            # short hand
            align = info.align_scale[part].full.align
            scale = info.align_scale[part].full.scale
            # if the align_scale is not set for range(0, dim_in):
            align_scale = complete_align_scale(part, align, scale)
            align = align_scale[0]
            scale = align_scale[1]

            # set the final values into the poly part object
            part.set_align(align)
            part.set_scale(scale)

    ''' normalizing the scaling factors '''
    norm = find_scale_norm(info)
    normalize_scale(norm, info)

    ''' increment alignments '''
    increment_aligns(comps, group_part_map)

    # ***
    log_level = logging.DEBUG
    LOG(log_level, "")
    LOG(log_level, "Final alignment and scaling")
    for comp in comps:
        for part in group_part_map[comp]:
            LOG(log_level, part.comp.func.name)
            log_str1 = "part.align = "+str([i for i in part.align])
            log_str2 = "part.scale = "+str([i for i in part.scale])
            LOG(log_level, log_str1)
            LOG(log_level, log_str2)
    # ***

    # ***
    log_level = logging.DEBUG
    LOG(log_level, "")
    LOG(log_level, "done ... align_scale_parts()")
    LOG(log_level, "___________________________________")
    # ***

    return True


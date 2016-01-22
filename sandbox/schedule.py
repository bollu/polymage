from __future__ import absolute_import, division, print_function

import logging
from grouping import get_group_dep_vecs
import expression
from utils import *
import poly as poly
from poly import *

# LOG CONFIG #
schedule_logger = logging.getLogger("schedule.py")
schedule_logger.setLevel(logging.DEBUG)
LOG = schedule_logger.log

class ASPacket(object):
    """
    contains alignment and scaling for a part
    """
    def __init__(self, _align, _scale):
        self.align = _align
        self.scale = _scale

class ASInfo(object):
    """
    contains the following types of scaling and alignment for a part:
    1. true : contains only those that can be inferred from the parent part
    2. full : [true] + random alignment and scaling for the remaining dims
    """
    def __init__(self, _true, _full):
        # _true, _full -> type : ASPacket
        # absolute alignment and scaling w.r.t base part
        self.true = _true
        self.full = _full

def format_schedule_constraints(dim_in, dim_out, align, scale, level_no):
    ineq_coeff = []
    eq_coeff   = []
    dim_set = [ False for i in range(0, dim_out) ]
    # Adding identity constraint for each dimension
    for i in range(0, dim_in):
        coeff = {}
        coeff[('out', align[i])] = 1
        if scale[i] != '-':
            assert scale[i] >= 1
        coeff[('in', i)] = -1 * scale[i]
        eq_coeff.append(coeff)
        dim_set[align[i]] = True

    # Setting the leading schedule dimension to level
    level_coeff = {}
    level_coeff[('out', 0)] = -1
    level_coeff[('constant', 0)] = level_no-1
    eq_coeff.append(level_coeff)

    # Setting the remaining dimensions to zero
    for i in range(1, dim_out):
        if not dim_set[i]:
            coeff = {}
            coeff[('out', i)] = 1
            coeff[('constant', 0)] = 0
            eq_coeff.append(coeff)
    return [ineq_coeff, eq_coeff]

def base_schedule(group):
    """
    Construct the base schedule for a group with a polyhedral representation.
    """

    assert(group.isPolyhedral)

    parts = []
    for sublist in group.polyRep.poly_parts.values():
        parts.extend(sublist)

    for part in parts:
        dim_in = part.sched.dim(isl._isl.dim_type.in_)
        dim_out = part.sched.dim(isl._isl.dim_type.out)
        [ineqs, eqs] = format_schedule_constraints(dim_in, dim_out,
                                                   part.align,
                                                   part.scale,
                                                   part._level_no)
        part.sched = poly.add_constraints(part.sched, ineqs, eqs)

    return parts

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
                align.append('-')
                scale.append('-')

    return (align, scale)

def align_and_scale(pipeline, group):

    """
    Alignment structure:
    [alignment]

    [alignment] of a dimension 'i', of a polypart is given by the i'th member
    of the list part.align. This value is a dimension in the root part or the
    part that is at the first level of the topologically sorted group. In case
    there are multiple parent parts, the same alignment should hold for all
    those parts. A mapping to '-' instead of an integer dimension indicates
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
                    ref_arg_vars.append('-')
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
            dims = [dim for dim in aln_scl.true.align
                          if dim != '-']
            if dim_max < len(dims):
                dim_max = len(dims)
                best = aln_scl

        if best == None:
            best = aln_scl_list[0]

        return best

    def complete_align_scale(part, align, scale):
        dim_in = part.sched.dim(isl._isl.dim_type.in_)
        assert dim_in > 0
        new_align = align
        new_scale = scale
        if not all(align[dim] != '-' for dim in range(0, dim_in)):
            used_dims = [dim for dim in align if dim != '-']
            universal = set(range(0, len(align)))
            avail_dims = \
                list(universal.difference(set(used_dims)))
            avail_dims.sort(reverse=True)
            for dim in range(0, dim_in):
                if align[dim] == '-':
                    new_align[dim] = avail_dims.pop()
                    new_scale[dim] = 1
        return (new_align, new_scale)

    def non_null(_list, null='-'):
        n_list = [x for x in _list if x != null]
        return n_list

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
                if var != '-':
                    dim_map[var] = dim
                dim += 1
            assert(dim <= max_dim)
            return dim_map

        def new_list():
            ''' return a new empty list of length max_dim '''
            return ['-' for i in range(0, max_dim)]

        def new_dict():
            ''' return a new empty dictionary of length max_dim '''
            _dict = {}
            for i in range(0, max_dim):
                _dict[i] = '-'
            return _dict

        def compute_abs(dim, dst_pack, rel_align_map, rel_scale_map, \
                        reverse=False, pseudo=False):
            root_dim = rel_align_map[dim]
            dim_align = '-'
            dim_scale = '-'
            if root_dim != '-':
                dim_align = dst_pack.align[root_dim]
                if dim_align != '-':
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
                if align[dim] == '-':
                    assert scale[dim] == '-'
                else:
                    assert scale[dim] != '-'
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
                if var != '-':
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
                                  if rel_align[dim] == '-']

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
                if true_align[dim] != '-':
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
                                     child_part.comp.variableDomain[0])
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
        ref_comp = ref.objectRef

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

        def collect_child_refs(children):
            refs = []
            for child in children:
                refs[child_part] = [ref for ref in child.refs]

        # if already visited
        if comp in info.solved:
            return

        part_map = info.group.polyRep.poly_parts
        all_children = info.pipe._comp_objs_children[comp]
        solved_children = [child for child in all_children \
                                   if child in info.solved]

        for child in solved_children:
            for child_part in part_map[child]:
                # collect the references made only to comp
                refs = [ref for ref in child_part.refs
                              if ref.objectRef == comp]

                # if the poly part makes no reference to any other compute object
                if not refs:
                    continue
                else:
                    # compute the alignment and scaling across references
                    for ref in refs:
                        solve_from_ref(child_part, ref, info, reverse=True)

        for part in part_map[comp]:
            assert part in info.align_scale

        return

    def solve_comp_parents(parents, info):
        '''
        Solves for alignment and scaling for poly parts of comps in the parents
        list, using the information of only the already solved children.
        Traverses in a breadth-first fashion.
        '''
        if parents == []:
            return

        parent_map = info.pipe._comp_objs_parents
        poly_parts = info.group.polyRep.poly_parts
        parent_order = {}
        for parent in parents:
            parent_order[parent] = poly_parts[parent][0]._level_no
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
            if parent in parent_map:
                grand_parents = [gp for gp in parent_map[parent] \
                                      if gp not in info.solved and \
                                         gp in info.group._comp_objs]
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

        comp_parts = info.group.polyRep.poly_parts[comp]
        all_parents = [p for p in info.pipe._comp_objs_parents[comp] \
                           if p in info.group._comp_objs]
        solved_pars = [par for par in all_parents \
                             if par in info.solved]

        # unsolved parents are the newly discovered parents
        discovered_parents = set(all_parents).difference(set(solved_pars))
        info.discovered = list(set(info.discovered).union(discovered_parents))

        # solve for each part
        for p in comp_parts:
            # collect the references to solved parents
            refs = [ref for ref in p.refs
                          if ref.objectRef in solved_pars]

            # if the poly part makes no reference to any other compute object
            if not refs:
                aln_scl = default_align_and_scale(p.sched, info.max_dim, shift=True)
                true_default = ASPacket(aln_scl[0], aln_scl[1])
                full_default = ASPacket(aln_scl[0], aln_scl[1])
                info.align_scale[p] = ASInfo(true_default, full_default)
            else:
                # compute the alignment and scaling across references
                for ref in refs:
                    solve_from_ref(p, ref, info)

            assert p in info.align_scale

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

        children_map = info.pipe._comp_objs_children
        poly_parts = info.group.polyRep.poly_parts
        child_order = {}
        for child in children:
            if child in info.group._comp_objs:
                child_order[child] = poly_parts[child][0]._level_no
                # index 0 picks the fist poly part
        sorted_order = sorted(child_order.items(), key=lambda x: x[1])
        sorted_children = [x[0] for x in sorted_order]

        for child in sorted_children:
            solve_comp_from_parents(child, info)
            info.solved.append(child)

        # collect unsolved children of children within the group
        all_grand_children = []
        for child in sorted_children:
            if child in children_map:
                grand_children = [gc for gc in children_map[child] \
                                       if gc not in info.solved and \
                                          gc in info.group._comp_objs]
                all_grand_children += grand_children
        all_grand_children = list(set(all_grand_children))

        solve_comp_children(all_grand_children, info)

        return

    ''' main '''
    comp_objs = group._comp_objs

    # list all parts with no self references and find the max dim
    max_dim = 0
    min_level = 1000000
    no_self_dep_parts = []
    for comp in comp_objs:
        for p in group.polyRep.poly_parts[comp]:
            p_dim_in = p.sched.dim(isl._isl.dim_type.in_)
            if not p.is_self_dependent:
                no_self_dep_parts.append(p)
                if max_dim < p_dim_in:
                    max_dim = p_dim_in
                if min_level > p._level_no:
                    min_level = p._level_no

    # begin from the topologically earliest part as the base for
    # alignment reference
    min_level_parts = [part for part in no_self_dep_parts \
                              if part._level_no == min_level]
    min_level_comps = list(set([p.comp for p in min_level_parts]))

    # prefer the comp appearing at a higher level in the 'pipeline'
    abs_min_level = 1000000
    for comp in min_level_comps:
        if abs_min_level > pipeline._level_order[comp]:
            abs_min_level = pipeline._level_order[comp]

    abs_min_comps = []
    abs_min_parts = []
    for comp in min_level_comps:
        if pipeline._level_order[comp] == abs_min_level:
            abs_min_comps.append(comp)
            abs_min_parts += group.polyRep.poly_parts[comp]

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
    log_level = logging.DEBUG-1
    LOG(log_level, "____")
    LOG(log_level, str(base_part.comp.name)+\
                   " (level : "+str(base_part._level_no)+")")

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
    for part in group.polyRep.poly_parts[base_comp]:
        # set default values for base parts
        align, scale = default_align_and_scale(part.sched, max_dim, shift=True)
        # update to the temporary info
        true_pack = ASPacket(align, scale)
        full_pack = ASPacket(align, scale)
        info.align_scale[part] = ASInfo(true_pack, full_pack)

    # recursively compute alignment and scaling for the family of base_comp
    if base_comp in pipeline._comp_objs_children:
        children = [c for c in pipeline._comp_objs_children[base_comp] \
                        if c in group._comp_objs]
        if children:
            solve_comp_children(pipeline._comp_objs_children[base_comp], info)

    # compute newly discovered parents iteratively until no new parent is
    # discovered.
    solve_comp_parents(info.discovered, info)

    # increment alignments
    def plus_one(align):
        plus_align = []
        for dim in align:
            if dim != '-':
                plus_align.append(dim+1)
            else:
                plus_align.append('-')
        return plus_align

    # set all the solutions into polypart object members
    for comp in comp_objs:
        for part in group.polyRep.poly_parts[comp]:
            # after solving, no poly part cannot not have a solution
            assert part in info.align_scale
            # short hand
            align = info.align_scale[part].full.align
            scale = info.align_scale[part].full.scale
            # if the align_scale is not set for range(0, dim_in):
            align_scale = complete_align_scale(part, align, scale)
            align = align_scale[0]
            scale = align_scale[1]

            plus_align = plus_one(align)

            # set the final values into the poly part object
            part.set_align(plus_align)
            part.set_scale(scale)

    ''' normalizing the scaling factors '''
    # normalize the scaling factors, so that none of them is lesser than 1
    norm = [1 for i in range(0, max_dim)]

    # compute the lcm of the Fraction denominators of scaling factors of all
    # poly parts in the group, for each dimension
    for part in info.align_scale:
        scale = part.scale
        for dim in range(0, max_dim):
            if scale[dim] != '-':
                d = Fraction(scale[dim].denominator)
                norm[dim] = lcm(d, norm[dim])

    LOG(logging.DEBUG, "")
    LOG(logging.DEBUG, "Final alignment and scaling")

    for part in info.align_scale:
        scale = part.scale
        new_scale = [1 for i in range(0, max_dim)]
        for dim in range(0, max_dim):
            if scale[dim] != '-':
                new_scale[dim] = int(norm[dim] * part.scale[dim])
            else:
                new_scale[dim] = '-'
        part.set_scale(new_scale)

        # ***
        log_level = logging.DEBUG
        LOG(log_level, part.comp.name)

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

def compute_tile_slope(dep_vecs, hmax):
    # Compute slopes
    # -- The first dimension in the domain gives the comp obj order. The slope
    #    of the tile in each dimension is computed w.r.t the comp obj order.
    #    The min extent and max extent in each dimension are computed. The
    #    hyperplanes representing the min and max extent give the shape of the
    #    tile in that dimension.
    #
    # =================================================================
    # '*' = live-out / live-in
    # '+' = intermediate
    # '-' = spurious compute due to over-approximation
    # (edges bw '*' and '+' are the actual dependences)
    #
    #  ^  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #  |               /|\               X                |\
    #  |              + + + + + + + + + + + + + + + + + + + +
    #  |             /| |              /   \              |  \
    #  |            - + + + + + + + + + + + + + + + + + + + + -
    #  |           / /|\|            /       \            |    \
    #  h          - + + + + + + + + + + + + + + + + + + + + + + -
    #  |         / /|\  |          /           \          |      \
    #  |        - + + + + + + + + + + + + + + + + + + + + + + + + -
    #  |       / /|\    |        /               \        |        \
    #  v  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #         <--- a --->       <------- o ------->       <--- b --->
    #
    # a = max_width
    # b = min_width
    # o = overlap
    # =================================================================

    if len(dep_vecs) < 1 :
        return ([], [])

    vec_len = len(dep_vecs[0][0])
    slope_min = [ (0, 1) for i in range(0, vec_len - 1) ]
    slope_max = [ (0, 1) for i in range(0, vec_len - 1) ]
    # Find max and min widths of dependencies at the base
    widths = []
    hmin = min([ dep[1] for dep in dep_vecs ])
    # extent cumulated at each level
    min_width = [ 0 for i in range(0, vec_len - 1) ]
    max_width = [ 0 for i in range(0, vec_len - 1) ]
    dep_unknown = [ False for i in range(0, vec_len - 1) ]
    # for each level
    for currh in range(hmax - 1, hmin - 1, -1):
        local_max_w = [ 0 for i in range(0, vec_len - 1) ]
        local_min_w = [ 0 for i in range(0, vec_len - 1) ]
        # dep vecs for this level
        h_dep_vecs = [ dep_vec for dep_vec in dep_vecs
                                 if dep_vec[1] == currh ]
        # h corresponds to topo-order height in the pipeline
        for dep_vec, h in h_dep_vecs:
            for i in range(0, len(dep_vec)-1):
                if dep_vec[i+1] == '*':
                    dep_unknown[i] = True
                    continue
                if dep_vec[i+1] > 0:
                    local_max_w[i] = max(local_max_w[i], dep_vec[i+1])
                if dep_vec[i+1] < 0:
                    local_min_w[i] = min(local_min_w[i], dep_vec[i+1])
        for i in range(0, len(dep_vec)-1):
            max_width[i] = max_width[i] + local_max_w[i]
            min_width[i] = min_width[i] + local_min_w[i]
        # remember the cumukated extent at each level
        widths.append((list(min_width), currh))
        widths.append((list(max_width), currh))

    # use the level-specific extents to determine the slope
    # (h will be descending, as 'widths' list was populated in that order)
    for width, h in widths:
        scale = hmax - h
        for i in range(0, vec_len-1):
            # be careful while constructing Fraction object
            # min slope
            if ((Fraction(width[i], scale) < Fraction(slope_min[i][0],
                                                      slope_min[i][1])) and \
                width[i] < 0):
                slope_min[i] = (width[i], scale)
            # max slope
            if ((Fraction(width[i], scale) > Fraction(slope_max[i][0],
                                                      slope_max[i][1])) and \
                width[i] > 0):
                slope_max[i] = (width[i], scale)

    for i in range(0, vec_len-1):
        if dep_unknown[i]:
            slope_min[i] = '*'
            slope_max[i] = '*'

    return (slope_min, slope_max)

def mark_par_and_vec_for_tile(poly_part):
    p = poly_part
    # -- Mark parallel dimensions and vector dimensions for tiles
    #    -- Find the outer most parallel dimension which can generate "enough"
    #       tasks for the given number of threads.
    #    -- Partial and full tile separation to enable better vectorization.
    outer_parallel_dim = None
    inner_vec_dim = None
    for dim in p.dim_tile_info:
        if p.dim_tile_info[dim][0] == 'none':
            # Either the dimension is too small to be parallelized or
            # is skewed. In both cases the dimension cannot be parallel.
            # This can change when we choose to not tile a dimension.
            continue
        elif p.dim_tile_info[dim][0] == 'overlap':
            dim_name = p.dim_tile_info[dim][1]
            tile_dim_name = p.dim_tile_info[dim][2]
            sched_dim = p.sched.find_dim_by_name(isl._isl.dim_type.out,
                                                 dim_name)
            tile_dim = p.sched.find_dim_by_name(isl._isl.dim_type.out,
                                                tile_dim_name)
            # update outermost parallel dim
            if outer_parallel_dim is not None:
                outer_parallel_dim = min(tile_dim, outer_parallel_dim)
            else:
                outer_parallel_dim = tile_dim
            # update outermost vector dim
            if inner_vec_dim is not None:
                inner_vec_dim = max(sched_dim, inner_vec_dim)
            else:
                inner_vec_dim = sched_dim

    # mark parallel
    if outer_parallel_dim is not None:
        p_dim_name = p.sched.get_dim_name(isl._isl.dim_type.out,
                                          outer_parallel_dim)
        p.parallel_sched_dims.append(p_dim_name)
    # mark vector
    if inner_vec_dim is not None:
        v_dim_name = p.sched.get_dim_name(isl._isl.dim_type.out,
                                        inner_vec_dim)
        p.vector_sched_dim.append(v_dim_name)

    return

def mark_par_and_vec(poly_part, param_estimates):
    p = poly_part
    # Determine the outer most dim and mark it parallel,
    # the inner most dim and mark it as vector
    parallel_dim = None
    vec_dim = None
    for dim in range(0, len(p.align)):
        interval = p.comp.domain[dim]
        if isinstance(p.comp, Reduction):
            interval = p.comp.reductionDomain[dim]
        # Since size could be estimated so can interval size be
        intr_size = \
            poly.get_dim_size(interval, param_estimates)

        # outer parallel dim
        if(get_constant_from_expr(intr_size) >= 32):
            if parallel_dim is not None:
                parallel_dim = min(p.align[dim], parallel_dim)
            else:
                parallel_dim = p.align[dim]

        # inner vector dim
        if(get_constant_from_expr(intr_size) >= 8):
            if vec_dim is not None:
                vec_dim = max(p.align[dim], vec_dim)
            else:
                vec_dim = p.align[dim]

    # mark parallel
    if parallel_dim is not None:
        p_dim_name = p.sched.get_dim_name(isl._isl.dim_type.out,
                                          parallel_dim)
        p.parallel_sched_dims.append(p_dim_name)
    # mark vector
    if vec_dim is not None:
        v_dim_name = p.sched.get_dim_name(isl._isl.dim_type.out,
                                          vec_dim)
        p.vector_sched_dim.append(v_dim_name)

    return

def enable_tile_scratchpad(group_parts):
    # Determine the buffer sizes for stages in each dimension
    for p in group_parts:
        for dim in p.dim_tile_info:
            if p.dim_tile_info[dim][0] == 'none':
                continue
            dim_name = p.dim_tile_info[dim][1]
            tile_dim_name = p.dim_tile_info[dim][2]
            extent = p.dim_tile_info[dim][3]
            if p.dim_tile_info[dim][0] == 'overlap':
                # Accounting for the overlap region
                left = p.dim_tile_info[dim][4]
                right = p.dim_tile_info[dim][5]
                h = p.dim_tile_info[dim][6]
                extent += abs(left * h) + abs(right * h)
            p.dim_scratch_size[dim] = \
                int(math.ceil(Fraction(extent, p.scale[dim])))
            mul_name = \
              '_Mul_'+p.sched.get_dim_name(isl._isl.dim_type.in_, dim)
            dim_in = p.sched.dim(isl._isl.dim_type.in_)
            dim_id =  p.sched.get_tuple_id(isl._isl.dim_type.in_)
            p.sched = p.sched.insert_dims(isl._isl.dim_type.in_, dim_in, 1)
            p.sched = p.sched.set_tuple_id(isl._isl.dim_type.in_, dim_id)
            p.sched = \
              p.sched.set_dim_name(isl._isl.dim_type.in_, dim_in, mul_name)
            sched_dim = \
              p.sched.find_dim_by_name(isl._isl.dim_type.out, dim_name)
            tile_dim = \
              p.sched.find_dim_by_name(isl._isl.dim_type.out, tile_dim_name)

            eqs = []
            coeff = {}
            coeff[('in', dim_in)] = p.scale[dim]
            coeff[('out', sched_dim)] = -1
            coeff[('out', tile_dim)] = p.dim_tile_info[dim][3]
            eqs.append(coeff)

            ineqs = []

            p.sched = poly.add_constraints(p.sched, ineqs, eqs)

    return

def schedule_time_within_group(part_comp_map):
    # Computations which have different scale but map to the same time
    # generate a lot of conditionals which can hinder performance. This
    # step separates all computations in a time step by adding an additional
    # dimension.
    pi = 0
    for comp in part_comp_map:
        for p in part_comp_map[comp]:
            time_dim = \
              p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')
            p.sched = \
              p.sched.insert_dims(isl._isl.dim_type.out, time_dim + 1, 1)
            p.sched = p.sched.set_dim_name(isl._isl.dim_type.out,
                                           time_dim + 1, '_o')

            eqs = []
            coeff = {}
            coeff[('constant', 0)] = -pi
            coeff[('out', time_dim + 1)] = 1
            eqs.append(coeff)
            p.sched = poly.add_constraints(p.sched, [], eqs)
            pi += 1

    return

def fused_schedule(pipeline, group, param_estimates):
    """Generate an optimized schedule for the stage."""

    g_poly_rep = group.polyRep
    g_poly_parts = g_poly_rep.poly_parts
    g_all_parts = []
    for comp in g_poly_parts:
        g_all_parts.extend(g_poly_parts[comp])

    # get dependence vectors between each part of the group and each of its
    # parents' part
    comp_deps = get_group_dep_vecs(group, g_all_parts)

    # No point in tiling a group that has no dependencies
    is_stencil = len(comp_deps) > 0 and len(g_all_parts) > 1
    for dep, h in comp_deps:
        # Skips groups which have self deps
        if dep[0] == 0:
            is_stencil = False

    # threshold for parallelism
    if not is_stencil:
        for p in g_all_parts:
            part_size = p.get_size(param_estimates)
            big_part = (part_size != '*' and \
                        part_size > pipeline._size_threshold)
            if not p.is_self_dependent and big_part:
                mark_par_and_vec(p, pipeline._param_estimates)

    # Find the parts which are not liveout
    for p in g_all_parts:
        is_liveout = not is_stencil
        #is_liveout = True
        p.is_liveout = p.is_liveout or is_liveout

    if is_stencil:
        assert(len(g_all_parts) > 1)
        hmax = max( [ p._level_no for p in g_all_parts ] )
        hmin = min( [ p._level_no for p in g_all_parts ] )
        slope_min, slope_max = compute_tile_slope(comp_deps, hmax)

        #splitTile(stageGroups[gi], slopeMin, slopeMax)
        overlap_tile(pipeline, g_all_parts, slope_min, slope_max)

        enable_tile_scratchpad(g_all_parts)

        for p in g_all_parts:
            mark_par_and_vec_for_tile(p)

        '''
        for p in g_all_parts:
            skewed_schedule(p)
        '''
        schedule_time_within_group(g_poly_parts)

    return

def move_independent_dim(dim, group_parts, stageDim):
    # Move the independent dimensions outward of the stage dimension.
    for part in group_parts:
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

def get_group_height(group_parts):
    min_height = min( [ part._level_no for part in group_parts ] )
    max_height = max( [ part._level_no for part in group_parts ] )
    return max_height - min_height

def overlap_tile(pipeline, group_parts, slope_min, slope_max):
    comp_dim = 0
    tile_dims = 0
    no_tile_dims = 0
    h = get_group_height(group_parts)
    num_tile_dims = 0
    for i in range(1, len(slope_min) + 1):
        # Check if every part in the group has enough iteration
        # points in the dimension to benefit from tiling.
        tile = False
        for part in group_parts:
            curr_dim = comp_dim + no_tile_dims + 2*tile_dims + 1
            lower_bound = part.sched.range().dim_min(curr_dim)
            upper_bound = part.sched.range().dim_max(curr_dim)
            size = upper_bound.sub(lower_bound)
            if (size.is_cst() and size.n_piece() == 1):
                aff = (size.get_pieces())[0][1]
                val = aff.get_constant_val()
                if val > pipeline._tile_sizes[num_tile_dims]:
                    tile = True
            else:
                tile = True
        if tile and slope_min[i-1] != '*':
            # Altering the schedule by constructing overlapped tiles
            for part in group_parts:
                # Extend space to accomodate the tiling dimensions
                part.sched = part.sched.insert_dims(
                                isl._isl.dim_type.out,
                                comp_dim + tile_dims, 1)
                # get the name of the untiled dim to name its corresponding
                # tiled dimension
                name = part.sched.get_dim_name(
                            isl._isl.dim_type.out,
                            comp_dim + no_tile_dims + 2*tile_dims + 2)
                part.sched = part.sched.set_dim_name(
                                isl._isl.dim_type.out,
                                comp_dim + tile_dims,
                                '_T' + name)
                right = int(math.floor(Fraction(slope_min[i-1][0],
                                                slope_min[i-1][1])))
                left = int(math.ceil(Fraction(slope_max[i-1][0],
                                              slope_max[i-1][1])))
                # L and R are normals to the left and the right
                # bounding hyperplanes of the uniform dependencies

                tile_size = pipeline._tile_sizes[num_tile_dims]
                # Compute the overlap shift
                #print(slope_max, slope_min, h, L, R, i-1)
                overlap_shift = abs(left * (h)) + abs(right * (h))
                for j in range(0, len(part.align)):
                    if i == part.align[j]:
                        assert j not in part.dim_tile_info
                        if tile_size%part.scale[j] != 0:
                            tile_size = int(math.ceil(part.scale[j]))
                        part.dim_tile_info[j] = ('overlap', name, '_T' + name,
                                                 tile_size, left, right, h)
                ineqs = []
                eqs = []
                coeff = {}
                it_dim = comp_dim + no_tile_dims + 2*tile_dims + 2
                tile_dim = comp_dim + tile_dims
                time_dim = comp_dim + tile_dims + 1

                coeff[('out', time_dim)] = -left
                coeff[('out', it_dim)] = 1
                coeff[('out', tile_dim)] = -tile_size
                ineqs.append(coeff)

                coeff = {}
                coeff[('out', time_dim)] = left
                coeff[('out', it_dim)] = -1
                coeff[('out', tile_dim)] = tile_size
                coeff[('constant', 0)] = tile_size - 1 + overlap_shift
                ineqs.append(coeff)
            
                coeff = {}
                coeff[('out', time_dim)] = -right
                coeff[('out', it_dim)] = 1
                coeff[('out', tile_dim)] = -tile_size
                ineqs.append(coeff)

                coeff = {}
                coeff[('out', time_dim)] = right
                coeff[('out', it_dim)] = -1
                coeff[('out', tile_dim)] = tile_size
                coeff[('constant', 0)] = tile_size + overlap_shift - 1
                ineqs.append(coeff)

                prior_dom = part.sched.domain()
                part.sched = poly.add_constraints(part.sched, ineqs, eqs)
                post_dom = part.sched.domain()

                assert(part.sched.is_empty() == False)
                # Tiling should not change the domain that is iterated over
                assert(prior_dom.is_equal(post_dom))
            tile_dims += 1
            num_tile_dims += 1
        else:
            #self.move_independent_dim(i, group_parts, comp_dim)
            name = part.sched.get_dim_name(isl._isl.dim_type.out, comp_dim)
            for part in group_parts:
                for j in range(0, len(part.align)):
                    if i == part.align[j]:
                        assert j not in part.dim_tile_info
                        part.dim_tile_info[j] = ('none', name)
            no_tile_dims += 1

    return

def splitTile(self, group, slopeMin, slopeMax):
    stageDim = 0
    dtileDims = 0
    numTileDims = 0
    for i in range(1, len(slopeMin) + 1):
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

def skewed_schedule(poly_part):
    # Second level storage savings can be achieved by utilizing modulo buffers
    # in the non-vector dimension. The fastest varying dimension is considered
    # the vector dimension and by this point should be the inner-most dimension.

    # Disabling this for two reasons
    # 1) The code generator generates awful code. There is no reason to expect
    #    it to generate anything nice.
    # 2) The dimension that has skewing applied to it need not be tiled. This
    #    has to be integrated into scheduling itself.
    p = poly_part
    one_dim = True
    for dim in p.dim_tile_info:
        if p.dim_tile_info[dim][0] == 'overlap' and one_dim:
            one_dim = False
            dim_name = p.dim_tile_info[dim][1]

            # Skewing the dimension
            sched_dim = \
              p.sched.find_dim_by_name(isl._isl.dim_type.out, dim_name)
            p.sched = \
              p.sched.insert_dims(isl._isl.dim_type.out, sched_dim  + 1, 1)
            p.sched = p.sched.set_dim_name(isl._isl.dim_type.out,
                                           sched_dim + 1, '_shift' + dim_name)
            time_dim = p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')
            right = p.dim_tile_info[dim][5]

            eqs = []
            coeff = {}
            coeff[('out', sched_dim)] = 1
            coeff[('out', time_dim)] = abs(right)
            coeff[('out', sched_dim + 1)] = -1
            eqs.append(coeff)
            p.sched = add_constraints(p.sched, [], eqs)

            p.sched = p.sched.remove_dims(isl._isl.dim_type.out, sched_dim, 1)

            # Moving time inside
            time_dim = p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')
            p.sched = p.sched.insert_dims(isl._isl.dim_type.out, time_dim, 1)
            p.sched = p.sched.set_dim_name(isl._isl.dim_type.out,
                                           time_dim, '_tmp' + dim_name)
            sched_dim = p.sched.find_dim_by_name(isl._isl.dim_type.out,
                                                 '_shift' + dim_name)

            eqs = []
            coeff = {}
            coeff[('out', time_dim)] = 1
            coeff[('out', sched_dim)] = -1
            eqs.append(coeff)
            p.sched = add_constraints(p.sched, [], eqs)
            p.sched = p.sched.remove_dims(isl._isl.dim_type.out, sched_dim, 1)
            p.sched = p.sched.set_dim_name(isl._isl.dim_type.out,
                                           time_dim, '_shift' + dim_name)

    return

def getDomainDimCoeffs(self, sched, arg):
    domDimCoeff = {}
    if (isAffine(arg)):
        coeff = getAffineVarAndParamCoeff(arg)
        for item in coeff:
            if type(item) == Variable:
                dim = sched.find_dim_by_name(isl._isl.dim_type.in_,
                                             item.name)
                domDimCoeff[dim] = coeff[item]
    return domDimCoeff

def getParamCoeffs(self, sched, arg):
    paramCoeff = {}
    if (isAffine(arg)):
        coeff = getAffineVarAndParamCoeff(arg)
        for item in coeff:
            if type(item) == Parameter:
                dim = sched.find_dim_by_name(isl._isl.dim_type.param,
                                             item.name)
                paramCoeff[dim] == coeff[item]
    return paramCoeff

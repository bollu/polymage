from __future__ import absolute_import, division, print_function

from poly import *
import logging
from grouping import get_group_dep_vecs
import expression

# LOG CONFIG #
schedule_logger = logging.getLogger("schedule.py")
schedule_logger.setLevel(logging.INFO)
LOG = schedule_logger.log

def get_parent_parts(part, group):
     refs = part.refs
     parent_parts = []
     for ref in refs:
         if ref.objectRef != part.comp:
             if ref.objectRef in group.polyRep.poly_parts:
                 parent_parts.extend(group.polyRep.poly_parts[ref.objectRef])
     return list(set(parent_parts))

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
        part.sched = add_constraints(part.sched.copy(), ineqs, eqs)

    return parts

def align_and_scale_parts(pipeline, group):

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

    def compatible_align(align1, align2):
        '''
        Treats alignment vectors of different length as incompatible
        '''
        compatible = True

        if not align1 and not align2:
            return compatible
        elif not align1 or not align2:
            return not compatible
        elif len(align1) == len(align2):
            for i in range(0, len(align1)):
                if not ((align1[i] == '-' or align2[i] == '-')
                        or (align1[i] == align2[i])):
                    compatible = False
                    break
        else:
            compatible = False
        return compatible

    def compatible_scale(scale1, scale2):
        '''
        Treats scaling vectors of different length as incompatible
        '''
        compatible = True

        if not scale1 and not scale2:
            return compatible
        elif not scale1 or not scale2:
            return not compatible
        elif len(scale1) == len(scale2):
            for i in range(0, len(scale1)):
                # fix this for scale[1]
                if not ((scale1[i] == '-' or scale2[i] == '-')
                        or (scale1[i] == scale2[i])):
                    compatible = False
                    break
        else:
            compatible = False
        return compatible

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

    def align_scale_vars(part, parent_part, ref_arg_vars, ref_arg_coeffs):
        '''
        Finds an alignment and scaling factor for each dimension associated
        with the reference argument variable
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

        # parent info
        parent_align = parent_part.align
        parent_scale = parent_part.scale
        max_dim = len(parent_align)

        # relative alignment and scaling
        # dict: dim -> dim
        rel_align = {}
        rel_scale = {}

        # start from null alignment
        part_align = ['-' for i in range(0, max_dim)]
        part_scale = ['-' for i in range(0, max_dim)]

        # dimensionality of the part
        part_dim_in = part.sched.dim(isl._isl.dim_type.in_)

        # dimensions of variable domain of the part and its parent
        part_dims = get_domain_dims(part.sched, part.comp.variableDomain[0])
        parent_dims = get_argvar_order(ref_arg_vars)

        # ***
        log_level = logging.DEBUG-2
        log_str = "aligning and scaling with all ref-part variables..."
        LOG(log_level, "")
        LOG(log_level, log_str)
        # ***

        # for each variable in the reference argument, get the alignment of
        # the part relative to the parent
        for var in ref_arg_vars:
            if var != '-':
                rel_align[part_dims[var]] = parent_dims[var]
                rel_scale[part_dims[var]] = ref_arg_coeffs[var]

        # dims of the part which didn't get an alignment
        rem_part_dims = [dim for dim in part_dims.values() \
                                 if dim not in rel_align]
        # dims of the parent part to which no part dim was aligned
        rem_parent_dims = [dim for dim in parent_dims.values() \
                                   if dim not in rel_align.values()]

        # align each of the remaining part dims to any remaining dim of
        # the parent part
        for dim in rem_part_dims:
            rel_scale[dim] = 1
            if rem_parent_dims:
                rel_align[dim] = rem_parent_dims.pop()
            else:
                rel_align[dim] = '*'

        aligned_dims = [dim for dim in part_dims.values()
                                if rel_align[dim] != '*']
        dangling_dims = [dim for dim in part_dims.values()
                                 if dim not in aligned_dims]

        # normalize to the base alignment, scaling using the relative
        # alignment, scaling
        for dim in aligned_dims:
            root_dim = rel_align[dim]
            part_align[dim] = parent_align[root_dim]
            part_scale[dim] = parent_scale[root_dim] * rel_scale[dim]

        # dangling_dims are assigned any available dim of the base alignment
        avail_dims = [dim for dim in range(0, max_dim) \
                            if dim not in part_align]
        for dim in dangling_dims:
            part_align[dim] = avail_dims.pop()
            part_scale[dim] = 1

        # test for unique alignment
        assert (len(part_align) == len(set(part_align)))

        # test for alignment boundary
        out_of_bound = all(max_dim >= dim for dim in part_align if dim != '-') \
                       and \
                       all(0 < dim for dim in part_align if dim != '-')
        assert (out_of_bound)

        return part_align, part_scale

    def align_scale_with_ref(part, ref, max_dim):
        ref_comp = ref.objectRef

        # initialize new alignment
        part_align = ['-' for i in range(0, max_dim)]
        part.set_align(part_align)

        old_align = part.align

        # initialize new scaling
        part_scale = ['-' for i in range(0, max_dim)]
        part.set_scale(part_scale)

        old_scale = part.scale

        # ***
        log_level = logging.DEBUG-2
        log_str = "aligning and scaling with all ref-parts..."
        LOG(log_level, "")
        LOG(log_level, log_str)
        # ***

        ref_poly_parts = group.polyRep.poly_parts[ref_comp]
        for ref_part in ref_poly_parts:
            part_align = part.align
            part_scale = part.scale

            # if old_align is not empty(init) and the alignment has changed
            no_conflict = compatible_align(old_align, part_align)
            if old_align and not no_conflict:
                return False, True  # or False, False

            # if old_scale is not empty(init) and the alignment has changed
            no_conflict = compatible_align(old_scale, part_scale)
            if old_scale and not no_conflict:
                return True, False

            old_align = part.align
            old_scale = part.scale

            ref_part_align = ref_part.align
            ref_part_scale = ref_part.scale

            # ***
            log_level = logging.DEBUG-2
            LOG(log_level, "")
            log_str1 = "ref_part_align = "+str([i for i in ref_part_align])
            log_str2 = "ref_part_scale = "+str([i for i in ref_part_scale])
            LOG(log_level, "ref = %s", str(ref_part.comp.name))
            LOG(log_level, log_str1)
            LOG(log_level, log_str2)
            # ***

            # process the argument list
            ref_args = ref.arguments
            ref_arg_vars, ref_arg_coeffs = extract_arg_vars_coefs(ref_args)

            # ***
            log_level = logging.DEBUG-2
            log_str1 = "ref_arg_vars  = "+ \
                      str([i.__str__() for i in ref_arg_vars])
            log_str2 = "ref_arg_coeffs = "+ \
                      str([(i.__str__(), ref_arg_coeffs[i]) \
                            for i in ref_arg_coeffs])
            LOG(log_level, log_str1)
            LOG(log_level, log_str2)
            # ***

            # match the part variables with the reference variables
            part_align, part_scale = \
                align_scale_vars(part, ref_part, ref_arg_vars, ref_arg_coeffs)

            # ***
            log_level = logging.DEBUG-2
            log_str1 = "old_align = "+str([i for i in old_align])
            log_str2 = "old_scale = "+str([i for i in old_scale])
            LOG(log_level, log_str1)
            LOG(log_level, log_str2)

            log_level = logging.DEBUG-1
            log_str1 = "part_align = "+str([i for i in part_align])
            log_str2 = "part_scale = "+str([i for i in part_scale])
            LOG(log_level, log_str1)
            LOG(log_level, log_str2)
            # ***

            part.set_align(part_align)
            part.set_scale(part_scale)

        no_conflict = compatible_align(old_align, part_align)
        if old_align and not no_conflict:
            return False, True  # or False, False

        no_conflict = compatible_align(old_scale, part_scale)
        if old_scale and not no_conflict:
            return True, False

        return True, True

    # BEGIN
    comp_objs = group._comp_objs

    # list all parts with no self references and find the max dim
    max_dim = 0
    no_self_dep_parts = []
    for comp in comp_objs:
        for p in group.polyRep.poly_parts[comp]:
            p_align = p.align
            if not p.is_self_dependent():
                no_self_dep_parts.append(p)
                # update size of align vector to max dim
                # assuming that 'align' has only spatial dims
                if max_dim < len(p_align):
                    max_dim = len(p_align)

    sorted_parts = sorted(no_self_dep_parts, \
                          key = lambda part:part._level_no)

    # begin from the topologically earliest part as the base for
    # alignment reference
    base_parts = [part for part in sorted_parts \
                       if part._level_no == sorted_parts[0]._level_no]

    # the alignment positions and scaling factors for variables follows
    # domain order of base parts
    base_align = [i+1 for i in range(0, max_dim)]
    base_scale = [1 for i in range(0, max_dim)]

    # initial alignment and scaling for all the base parts
    for part in base_parts:
        # ***
        log_level = logging.DEBUG-1
        LOG(log_level, "____")
        LOG(log_level, str(part.comp.name)+\
                       " (level : "+str(part._level_no)+")")
        # ***

        part.set_align(base_align)
        part.set_scale(base_scale)

    other_parts = [part for part in sorted_parts \
                        if part not in base_parts]
    for part in other_parts:
        # ***
        log_level = logging.DEBUG-1
        LOG(log_level, "____")
        LOG(log_level, str(part.comp.name)+\
                       " (level : "+str(part._level_no)+")")
        # ***

        old_align = []
        old_scale = []
        part.set_align(old_align)
        part.set_scale(old_scale)

        part_align = part.align
        part_scale = part.scale

        # ***
        log_level = logging.DEBUG-2
        LOG(log_level, "")
        LOG(log_level, "aligning and scaling with all refs...")
        # ***

        refs = part.refs
        if not refs:
            part.set_align(base_align)
            part.set_scale(base_scale)
        for ref in refs:
            no_conflict = compatible_align(part_align, old_align)
            if old_align and not no_conflict:
                LOG(logging.ERROR, "Conflict in alignment across refs")
                return False

            no_conflict = compatible_scale(part_scale, old_scale)
            if old_scale and not no_conflict:
                LOG(logging.ERROR, "Conflict in scaling across refs")
                return False

            old_align = part.align
            old_scale = part.scale

            # Alignment and scaling with references
            no_align_conflict, \
              no_scale_conflict = \
                align_scale_with_ref(part, ref, max_dim)

            if old_align and not no_align_conflict:
                LOG(logging.ERROR, \
                    "Conflict in alignment across ref parts")
                return False
            if old_scale and not no_scale_conflict:
                LOG(logging.ERROR, \
                    "Conflict in scaling across ref parts")
                return False

            part_align = part.align
            part_scale = part.scale

        no_conflict = compatible_align(part_align, old_align)
        if old_align and not no_conflict:
            LOG(logging.ERROR, "Conflict in alignment across refs")
            return False

        no_conflict = compatible_scale(part_scale, old_scale)
        if old_scale and not no_conflict:
            LOG(logging.ERROR, "Conflict in scaling across refs")
            return False

    # normalize the scaling factors, so that none of them is lesser than 1
    norm = [1 for i in range(0, max_dim)]

    # compute the lcm of the Fraction denominators of all scaling factors
    # for each dimension
    for part in other_parts:
        scale = part.scale
        for dim in range(0, max_dim):
            if scale[dim] != '-':
                d = Fraction(scale[dim].denominator)
                norm[dim] = lcm(d, norm[dim])

    LOG(logging.DEBUG, "")
    LOG(logging.DEBUG, "Final alignment and scaling")

    for part in sorted_parts:
        scale = part.scale
        new_scale = [1 for i in range(0, max_dim)]
        for dim in range(0, max_dim):
            if scale[dim] != '-':
                new_scale[dim] = norm[dim] * part.scale[dim]
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
    LOG(log_level, "done ... align_parts()")
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
            get_dim_size(interval, param_estimates)

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

            p.sched = add_constraints(p.sched, ineqs, eqs)

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
            p.sched = add_constraints(p.sched, [], eqs)
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
            if not p.is_self_dependent() and big_part:
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
                part.sched = add_constraints(part.sched, ineqs, eqs)
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

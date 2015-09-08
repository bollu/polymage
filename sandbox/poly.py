from __future__ import absolute_import, division, print_function

import math

import islpy as isl

from constructs import *
from expression import *

def lcm(a, b):
    return a*b/(gcd(a, b))

def optimizeSchedule(partScheds, dependencies):
    # The pluto optimizer can be used to optimize the schedule for
    # comparision.
    pass

def add_constraints_from_list(obj, local_space, constraint_list,
                              constraint_alloc):
    for constr in constraint_list:
        c = constraint_alloc(local_space)

        # find the normalization factor
        m = 1
        for coeff in constr:
            if isinstance(constr[coeff], Fraction):
                den = int(gcd(abs(constr[coeff].denominator), m))
                m = (abs(constr[coeff].denominator) * m)//den
        assert m.denominator == 1
        m = m.numerator

        # normalize
        for coeff in constr:
            if isinstance(constr[coeff], Fraction):
               constr[coeff] = m * constr[coeff]
               assert constr[coeff].denominator == 1
               constr[coeff] = constr[coeff].numerator
            else:
               constr[coeff] = m * constr[coeff]

        for coeff in constr:
            dim = coeff[1]
            try:
                if coeff[0] == 'param':
                    if (type(dim) == str):
                        dim = \
                            obj.find_dim_by_name(isl._isl.dim_type.param, dim)
                    c = c.set_coefficient_val(isl._isl.dim_type.param,
                                              dim, constr[coeff])
                elif coeff[0] == 'in':
                    if (type(dim) == str):
                        dim = obj.find_dim_by_name(isl._isl.dim_type.in_, dim)
                    c = c.set_coefficient_val(isl._isl.dim_type.in_,
                                              dim, constr[coeff])
                elif coeff[0] == 'out':
                    if (type(dim) == str):
                        dim = obj.find_dim_by_name(isl._isl.dim_type.out, dim)
                    c = c.set_coefficient_val(isl._isl.dim_type.out,
                                              dim, constr[coeff])
                elif coeff[0] == 'constant':
                    c = c.set_constant_val(constr[coeff])
                else:
                   assert False
            except isl.Error:
                # Ignore this constraint conjunct since the referred dimension
                # is not scheduled in the obj. This happens when we try to add
                # constraint for a dimension that is not at all used by a part.
                # FIXME: isl's find_dim_by_name throws exception on not finding
                # any scheduled dimension. It's better to replace the exception
                # handling with an isl function, if any, to test for the
                # existence of a dimension in that part.
                pass

        obj = obj.add_constraint(c)
    return obj

def add_constraints(obj, ineqs, eqs):

    def add_constraints_for_element(obj, local_space, ineqs, eqs):
        obj = add_constraints_from_list(obj, local_space, ineqs,
                                        isl.Constraint.inequality_alloc)
        obj = add_constraints_from_list(obj, local_space, eqs,
                                        isl.Constraint.equality_alloc)
        return obj

    space = obj.get_space()
    if (isinstance(obj, isl.Map)):
        for bmap in obj.get_basic_maps():
            local_space = bmap.get_local_space()
            obj = add_constraints_for_element(obj, local_space, ineqs, eqs)
    elif (isinstance(obj, isl.Set)):
        for bset in obj.get_basic_sets():
            local_space = bset.get_local_space()
            obj = add_constraints_for_element(obj, local_space, ineqs, eqs)
    elif (isinstance(obj, isl.BasicSet) or
          isinstance(obj, isl.BasicMap)):
        local_space = obj.get_local_space()
        obj = add_constraints_for_element(obj, local_space, ineqs, eqs)
    else:
        assert False

    return obj

def extract_value_dependence(part, ref, refPolyDom):
    # Dependencies are calculated between values. There is no storage
    # mapping done yet.
    assert(part.sched)
    deps = []
    accessRegion = isl.BasicSet.universe(refPolyDom.dom_set.get_space())
    partDom = part.sched.domain().align_params(refPolyDom.dom_set.get_space())
    accessRegion = accessRegion.align_params(partDom.get_space())

    rel = isl.BasicMap.from_domain_and_range(partDom, accessRegion)
    dimOut = rel.dim(isl._isl.dim_type.out)
    sourceDims = [ ('out', i) for i in range(0, dimOut)]
    numArgs = len(ref.arguments)                

    for i in range(0, numArgs):
        arg = ref.arguments[i]
        # If the argument is not affine the dependence reflects that
        # the computation may depend on any value of the referenced object
        if (isAffine(arg)):
            coeff = get_affine_var_and_param_coeff(arg)
            coeff = map_coeff_to_dim(coeff)

            coeff[('constant', 0)] = get_constant_from_expr(arg, affine = True)
            coeff[sourceDims[i]] = -1
            rel = add_constraints(rel, [], [coeff])
    if not rel.is_empty():
        deps.append(PolyDep(ref.objectRef, part.comp, rel))
    return deps 

class PolyPart(object):
    def __init__(self, _sched, _expr, _pred, _comp,
                 _align, _scale, _level_no):
        self.sched = _sched
        self.expr = _expr
        self.pred = _pred
        self.comp = _comp
        # Dependencies between values of computation objects
        self.deps = []
        # Mapping between the input variables to the corresponding 
        # schedule dimension. A full affine schedule will need a 
        # transformation matrix. Currently we only shuffle the 
        # dimension order apart from tiling so a simple dimension
        # alignment vector suffices. This has to be changed to 
        # handle more general cases later.
        self._align = _align
        # Scaling factors for each schedule dimension
        self._scale = _scale
        # Default alignment and scaling factors are set while
        # constructing the polypart. These are changed by the
        # alignment and loop scaling passes. Both these passer
        # attempt to improve locality and uniformize dependencies.
        self._level_no = _level_no

        # maps tiled dimensions to their respective scratchpad sizes
        self.dim_scratch_size = {}

    @property
    def align(self):
        align_clone = [i for i in self._align]
        return align_clone

    @property
    def scale(self):
        scale_clone = [i for i in self._scale]
        return scale_clone

    def set_align(self, align):
        self._align = [i for i in align]
        return

    def set_scale(self, _scale):
        self._scale = [i for i in _scale]
        return

    def getPartRefs(self):
        refs = self.expr.collect(Reference)
        if (self.pred):
            refs += self.pred.collect(Reference)
        return refs

    def __str__(self):
        partStr = "Schedule: " + self.sched.__str__() + '\n'\
                  "Expression: " + self.expr.__str__() + '\n'\
                  "Predicate: " + self.pred.__str__() + '\n'
        depstr = ""
        for dep in self.deps:
            depstr = depstr + dep.__str__() + '\n'
        return partStr + depstr

class PolyDomain(object):
    def __init__(self, _dom_set, _comp):
        self.dom_set = _dom_set
        self.comp = _comp

    def __str__(self):
        return "Domain: " + self.dom_set.__str__()

# NOTE: Dead Code?
class PolyDep(object):
    def __init__(self, _producer_obj, _consumer_obj, _rel):
        self.producer_obj = _producer_obj
        self.consumer_obj = _consumer_obj
        self.rel         = _rel

    def __str__(self):
        return self.rel.__str__()

class PolyRep(object):
    """ The PolyRep class is the polyhedral representation of a 
        group. It gives piece-wise domain and schedule for each compute
        object in the group. Polyhedral transformations modify the 
        piece-wise domains as well as the schedules.
    """
    def __init__(self, _ctx, _group, _outputs,
                 _param_constraints):

        self.group = _group
        self.outputs = _outputs
        self.param_constraints = _param_constraints
        self.ctx = _ctx

        self.poly_parts = {}
        self.poly_doms = {}
        self.polyast = []

        self._var_count = 0
        self._func_count = 0

        # TODO: move the following outside __init__()
        # For now, let this be. Compilation optimizations can come later.

        self.extract_polyrep_from_group(_param_constraints)

        #self.fusedSchedule(_param_estimates)
        #self.simpleSchedule(_param_estimates)

    def extract_polyrep_from_group(self, param_constraints):
        # dict: comp_obj -> level_no
        comp_objs = self.group.orderComputeObjs()
        num_objs = len(comp_objs.items())

        # Comute the max dimensionality of the compute objects
        def max_dim(objs):
            dim = 0
            for comp in objs:
                if type(comp) == Reduction:
                    dim = max(dim, len(comp.reductionVariables))
                    dim = max(dim, len(comp.variables))
                elif type(comp) == Function or type(comp) == Image:
                    dim = max(dim, len(comp.variables))
            return dim

        dim = max_dim(comp_objs)

        # Get all the parameters used in the group compute objects
        grp_params = []
        for comp in comp_objs:
            grp_params = grp_params + comp.getObjects(Parameter)
        grp_params = list(set(grp_params))

        param_names = [param.name for param in grp_params]

        # Represent all the constraints specified on the parameters relevant
        # to the group.
        context_conds = \
            self.format_param_constraints(param_constraints, grp_params)

        # The [t] is for the stage dimension
        schedule_names = ['_t'] + \
                         [ self.getVarName()  for i in range(0, dim) ]

        for comp in comp_objs:
            if (type(comp) == Function or type(comp) == Image):
                self.extract_polyrep_from_function(comp, schedule_names,
                                                   param_names, context_conds,
                                                   comp_objs[comp]+1,
                                                   param_constraints)
            elif (type(comp) == Reduction):
                self.extract_polyrep_from_reduction(comp, schedule_names,
                                                    param_names, context_conds,
                                                    comp_objs[comp]+1,
                                                    param_constraints)
            else:
                assert False

    def format_param_constraints(self, param_constraints, grp_params):
        context_conds = []
        grp_params_set = set(grp_params)
        for param_constr in param_constraints:
            # Only consider parameter constraints of parameters
            # given in params.
            params_in_constr = param_constr.collect(Parameter)
            context_add = set(params_in_constr).issubset(grp_params_set)

            # Only add the constraint if it is affine and has no conjunctions.
            # Handling conjunctions can be done but will require more care.
            if context_add and isAffine(param_constr):
                param_constr_conjunct = param_constr.splitToConjuncts()
                if len(param_constr_conjunct) == 1:
                    context_conds.append(param_constr)
        return context_conds

    def extract_polyrep_from_function(self, comp,
                                      schedule_names, param_names,
                                      context_conds, level_no,
                                      param_constraints):

        sched_map = self.create_sched_space(comp.variables, comp.domain,
                                            schedule_names, param_names,
                                            context_conds)

        poly_dom = PolyDomain(sched_map.domain(), comp)
        id_ = isl.Id.alloc(self.ctx, comp.name, poly_dom)
        poly_dom.dom_set = poly_dom.dom_set.set_tuple_id(id_)
        self.poly_doms[comp] = poly_dom

        self.create_poly_parts_from_definition(comp, sched_map, level_no, 
                                               schedule_names, comp.domain)

    def extract_polyrep_from_reduction(self, comp,
                                       schedule_names, param_names,
                                       context_conds, level_no,
                                       param_constraints):

        sched_map = self.create_sched_space(comp.reductionVariables,
                                            comp.reductionDomain,
                                            schedule_names, param_names,
                                            context_conds)

        poly_dom = PolyDomain(sched_map.domain(), comp)
        id_ = isl.Id.alloc(self.ctx, comp.name, poly_dom)
        poly_dom.dom_set = poly_dom.dom_set.set_tuple_id(id_)
        self.poly_doms[comp] = poly_dom

        self.create_poly_parts_from_definition(comp, sched_map, level_no,
                                               schedule_names,
                                               comp.reductionDomain)

        dom_map = self.create_sched_space(comp.variables, comp.domain,
                                          schedule_names, param_names,
                                          context_conds)

        # Initializing the reduction earlier than any other function
        self.create_poly_parts_from_default(comp, dom_map, -1, schedule_names)


    def create_sched_space(self, variables, domains,
                           schedule_names, param_names, context_conds):
        # Variable names for referring to dimensions
        var_names = [ var.name for var in variables ]
        space = isl.Space.create_from_names(self.ctx, in_ = var_names,
                                                      out = schedule_names,
                                                      params = param_names)

        sched_map = isl.BasicMap.universe(space)
        # Adding the domain constraints
        [ineqs, eqs] = format_domain_constraints(domains, var_names)
        sched_map = add_constraints(sched_map, ineqs, eqs)

        # Adding the parameter constraints
        [param_ineqs, param_eqs] = format_conjunct_constraints(context_conds)
        sched_map = add_constraints(sched_map, param_ineqs, param_eqs)

        return sched_map

    def create_poly_parts_from_definition(self, comp,
                                          sched_map, level_no,
                                          schedule_names, domain):
        self.poly_parts[comp] = []
        for case in comp.defn:
            sched_m = sched_map.copy()

            # The basic schedule is an identity schedule appended with
            # a level dimension. The level dimension gives the ordering
            # of the compute objects within a group.

            align, scale = self.default_align_and_scale(sched_m)

            if (isinstance(case, Case)):
                # Dealing with != and ||. != can be replaced with < || >.
                # and || splits the domain into two.
                split_conjuncts = case.condition.splitToConjuncts()
                for conjunct in split_conjuncts:
                    # If the condition is non-affine it is stored as a
                    # predicate for the expression. An affine condition
                    # is added to the domain.
                    affine = True
                    for cond in conjunct:
                        affine = affine and \
                                 isAffine(cond.lhs) and isAffine(cond.rhs)
                    if(affine):
                        [conjunct_ineqs, conjunct_eqs] = \
                            format_conjunct_constraints(conjunct)
                        sched_m = add_constraints(sched_m,
                                                  conjunct_ineqs,
                                                  conjunct_eqs)
                        parts = self.make_poly_parts(sched_m, case.expression,
                                                     None, comp,
                                                     align, scale, level_no)
                        # FIXME: Is a loop required here? make_poly_part
                        # seems to return a list of one part
                        for part in parts:
                            self.poly_parts[comp].append(part)
                    else:
                        parts = self.make_poly_parts(sched_m, case.expression,
                                                     case.condition, comp,
                                                     align, scale, level_no)

                        # FIXME: Is a loop required here? make_poly_part
                        # seems to return a list of one part
                        for part in parts:
                            self.poly_parts[comp].append(part)
            else:
                assert(isinstance(case, AbstractExpression) or
                       isinstance(case, Accumulate))
                parts = self.make_poly_parts(sched_m, case,
                                             None, comp,
                                             align, scale, level_no)
                # FIXME: Is a loop required here? make_poly_part
                # seems to return a list of one part
                for part in parts:
                    self.poly_parts[comp].append(part)

        # TODO adding a boundary padding and default to the function 
        # will help DSL usability. 

        # An attempt to subtract all the part domains to find the domain
        # where the default expression has to be applied. 

        #sched_m = isl.BasicMap.identity(self.polyspace)
        #sched_m = add_constraints(sched, ineqs, eqs)
        # Adding stage identity constraint
        #level_coeff = {}
        #level_coeff[varDims[0]] = -1
        #level_coeff[('constant', 0)] = compObjs[comp]
        #sched_m = add_constraints(sched_m, [], [level_coeff])
        #sched_m = add_constraints(sched_m, param_ineqs, param_eqs)

        #for part in self.poly_parts[comp]:
        #    sched_m = sched_m.subtract_range(part.sched.range())
        #    if (sched_m.is_empty()):
        #        break
        #if(not sched_m.fast_is_empty()):
        #    bmap_list = []
        #    if (isinstance(sched_m, isl.BasicMap)):
        #        bmap_list.append(sched_m)
        #    else:
        #        sched_m.foreach_basic_map(bmap_list.append)
        #    for bmap in bmap_list:
        #        poly_part = PolyPart(bmap, comp.default, None, comp)
        #        id_ = isl.Id.alloc(self.ctx, comp.name, poly_part)
        #        poly_part.sched = poly_part.sched.set_tuple_id(
        #                                   isl._isl.dim_type.in_, id_)
        #        self.poly_parts[comp].append(poly_part)

    def create_poly_parts_from_default(self, comp, sched_map,
                                       level_no, schedule_names):
        sched_m = sched_map.copy()
        align, scale = self.default_align_and_scale(sched)

        assert(isinstance(comp.default, AbstractExpression))
        poly_part = PolyPart(sched_m, comp.default,
                             None, comp,
                             align, scale, level_no)

        id_ = isl.Id.alloc(self.ctx, comp.name, poly_part)
        poly_part.sched = \
                poly_part.sched.set_tuple_id(isl._isl.dim_type.in_, id_)
        self.poly_parts[comp].append(poly_part)

    def make_poly_parts(self, sched_map, expr, pred, comp,
                        align, scale, level_no):
        # Detect selects with modulo constraints and split into 
        # multiple parts. This technique can also be applied to the
        # predicate but for now we focus on selects.
        poly_parts = []
        # This is very very temporary solution there should be a 
        # better way of doing this. Only targetting conditions 
        # of the form (affine)%constant == constant.
        broken_parts = []
        if isinstance(expr, Select):
            conjuncts = expr.condition.splitToConjuncts()
            print("conjuncts =", conjuncts)
            if len(conjuncts) == 1 and len(conjuncts[0]) == 1:
                cond = conjuncts[0][0]
                left_expr = cond.lhs
                right_expr = cond.rhs
                is_left_modulo = isAffine(left_expr, includeModulo=True) and \
                                 not isAffine(left_expr)
                is_right_constant = is_constant_expr(right_expr)
                break_select = False
                # check for 'affine % constant == constant'
                if is_left_modulo and is_right_constant and \
                    cond.conditional == '==' and \
                    isinstance(left_expr, AbstractBinaryOpNode)\
                    and left_expr.op == '%' and isAffine(left_expr.left)\
                    and is_constant_expr(left_expr.right):
                    break_select = True
                if break_select:
                    left_const = get_constant_from_expr(left_expr.left,
                                                        affine = True)
                    right_const = get_constant_from_expr(right_expr,
                                                         affine = True)
                    mod_const = get_constant_from_expr(left_expr.right,
                                                       affine = True)

                    left_coeff = get_affine_var_and_param_coeff(left_expr.left)
                    left_coeff = map_coeff_to_dim(left_coeff)

                    mul_name = '_Mul_'
                    rem_name = '_Rem_'

                    # true branch schedule
                    true_sched = sched_map.copy()
                    dim_in = true_sched.dim(isl._isl.dim_type.in_)
                    true_sched = \
                        true_sched.insert_dims(isl._isl.dim_type.in_,
                                               dim_in, 1)
                    true_sched = \
                        true_sched.set_dim_name(isl._isl.dim_type.in_,
                                                dim_in, mul_name)

                    eqs = []
                    left_coeff[('constant', 0)] = left_const - right_const
                    left_coeff[('in', dim_in)] = -mod_const
                    eqs.append(left_coeff)

                    true_sched = add_constraints(true_sched, [], eqs)
                    true_sched = true_sched.project_out(isl._isl.dim_type.in_,
                                                        dim_in, 1)
                    broken_parts.append((true_sched, expr.trueExpression))

                    # false branch schedule
                    false_sched = sched_map.copy()
                    dim_in = false_sched.dim(isl._isl.dim_type.in_)
                    false_sched = \
                        false_sched.insert_dims(isl._isl.dim_type.in_,
                                                dim_in, 2)
                    false_sched = \
                        false_sched.set_dim_name(isl._isl.dim_type.in_,
                                                 dim_in, mul_name)
                    false_sched = \
                        false_sched.set_dim_name(isl._isl.dim_type.in_,
                                                 dim_in+1, rem_name)

                    eqs = []
                    left_coeff[('constant', 0)] = left_const - right_const
                    left_coeff[('in', dim_in)] = -mod_const
                    left_coeff[('in', dim_in+1)] = -1
                    eqs.append(left_coeff)

                    ineqs = []
                    coeff = {}
                    coeff[('in', dim_in+1)] = 1
                    coeff[('constant', 0)] = -1
                    ineqs.append(coeff)

                    coeff = {}
                    coeff[('in', dim_in+1)] = -1
                    coeff[('constant', 0)] = mod_const-1
                    ineqs.append(coeff)

                    false_sched = add_constraints(false_sched, ineqs, eqs)
                    false_sched = \
                        false_sched.project_out(isl._isl.dim_type.in_,
                                                dim_in, 2)
                    broken_parts.append((false_sched, expr.falseExpression))

        # Note the align and scale lists are cloned otherwise all the parts
        # will be sharing the same alignment and scaling
        if not broken_parts:
            poly_part = PolyPart(sched_map, expr, pred, comp,
                                 list(align), list(scale), level_no)
            # Create a user pointer, tuple name and add it to the map
            id_ = isl.Id.alloc(self.ctx, comp.name, poly_part)
            poly_part.sched = poly_part.sched.set_tuple_id(
                                          isl._isl.dim_type.in_, id_)
            poly_parts.append(poly_part)
        else:
            for bsched_map, bexpr in broken_parts:
                poly_part = PolyPart(bsched_map, bexpr, pred, comp,
                                     list(align), list(scale), level_no)
                # Create a user pointer, tuple name and add it to the map
                id_ = isl.Id.alloc(self.ctx, comp.name, poly_part)
                poly_part.sched = poly_part.sched.set_tuple_id( \
                                        isl._isl.dim_type.in_, id_)
                poly_parts.append(poly_part)
        return poly_parts

    def default_align_and_scale(self, sched):
        dim_in = sched.dim(isl._isl.dim_type.in_)
        align = {}
        # align[i] = j means input dimension i is mapped to output 
        # dimension j
        for i in range(0, dim_in):
            align[i] = [i+1]
        # the default scaling in each dimension is set to 1 i.e., the
        # schedule dimension correspoinding to input dimension will be 
        # scaled by 1
        scale = [1 for i in range(0, dim_in)]
        return (align, scale)

    def generateCode(self):
        self.polyast = []
        if self.poly_parts:
            self.buildAst()

    def buildAst(self):
        #astbld =  isl.AstBuild.from_context(isl.BasicSet("[C, R]->{: R>=1 and C>=1}", self.ctx))
        parts = []
        for plist in self.poly_parts.values():
            parts.extend(plist)
       
        # TODO figure out a way to create the correct parameter context
        # since the parameters for all the parts may not be the same
        astbld =  isl.AstBuild.from_context(parts[0].sched.params())
        #astbld =  astbld.set_options(isl.UnionMap("{ }"))

        sched_map = None
        optMap = None
        for part in parts:
            if sched_map is None:
                sched_map = isl.UnionMap.from_map(part.sched)
            else:
                partMap = isl.UnionMap.from_map(part.sched)
                sched_map = sched_map.union(partMap)
            if optMap is None:
                #unrollUnionSet = isl.UnionSet.from_set(isl.Set("{unroll[x] : x = 0 or x = 2}", self.ctx))
                unrollUnionSet = isl.UnionSet.from_set(isl.Set("{:}", self.ctx))
                domUnionSet = isl.UnionSet.universe(isl.UnionSet.from_set(part.sched.range()))
                optMap = isl.UnionMap.from_domain_and_range(domUnionSet, unrollUnionSet)
            else:
                #unrollUnionSet = isl.UnionSet.from_set(isl.Set("{unroll[x] : x = 0 or x = 2}", self.ctx))
                unrollUnionSet = isl.UnionSet.from_set(isl.Set("{:}", self.ctx))
                domUnionSet = isl.UnionSet.universe(isl.UnionSet.from_set(part.sched.range()))
                optMap = optMap.union(isl.UnionMap.from_domain_and_range(domUnionSet, unrollUnionSet))
        astbld = astbld.set_options(optMap)

        # All parts in the group will have the same schedule dimension
        # using the first part as the canonical one
        numIds = parts[0].sched.dim(isl._isl.dim_type.out)
        ids = isl.IdList.alloc(self.ctx, numIds)
        for i in range(0, numIds):
            schedName = parts[0].sched.get_dim_name(isl._isl.dim_type.out, i)
            ids = ids.add(isl.Id.alloc(self.ctx, schedName, None))
        astbld = astbld.set_iterators(ids)

        def printer(arg):
            #print(arg)
            pass

        sched_map.foreach_map(printer)
        self.polyast.append(astbld.ast_from_schedule(sched_map))

    # NOTE: Dead code?
    def computeDependencies(self):
        # Extract dependencies. In case the dependencies cannot be exactly 
        # represented they are approximated.
        for comp in self.poly_parts:
            for part in self.poly_parts[comp]:
                refs = part.getPartRefs()
                # There could be multiple references to the same value. So the 
                # dependencies might be redundant this has to be eliminated 
                # later.
                for ref in refs:
                    # Considering only dependencies within the same stage
                    if ref.objectRef in self.poly_parts:
                        part.deps += extract_value_dependence(part, ref,
                                                self.poly_doms[ref.objectRef])

    def getDependenceVector(self, parentPart, childPart, ref, scaleMap = None):
        def getScale(sMap, p, i):
            if sMap is not None:
                return sMap[p][i]
            return p.scale[i]

        numArgs = len(ref.arguments)
        dimOut = parentPart.sched.dim(isl._isl.dim_type.out)
        depVec = [ '-' for i in range(0, dimOut) ]

        if isinstance(parentPart.comp, Accumulator):
            for i in range(1, dimOut):
                depVec[i] = '*'
            depVec[0] = childPart.level_no - parentPart.level_no
            return (depVec, parentPart.level_no)

        for i in range(0, numArgs):
            arg = ref.arguments[i]
            parentVarSchedDim = parentPart.align[i]            
            if (isAffine(arg)):
                domDimCoeff = self.getDomainDimCoeffs(childPart.sched, arg)
                paramCoeff = self.getParamCoeffs(childPart.sched, arg)
                # Parameter coefficents can also be considered to
                # generate parametric shifts. Yet to be seen.
                
                # Indexed with multiple variables.
                if (len(domDimCoeff) > 1 or 
                    (len(domDimCoeff) == 1 and len(paramCoeff) >=1)):
                    # Although there are multiple coefficients. If 
                    # there is only one variable coefficient and other
                    # parametric coefficients. Uniformization can be 
                    # done with parametric shifts. Full affine scheduling 
                    # might be able to find a way to uniformize 
                    # dependencies. This has to be further explored.
                    #print(ref)
                    #assert False
                    depVec[parentVarSchedDim] = '*'
                # Indexed with a single variable. This can either be an 
                # uniform access or can be uniformized with scaling when 
                # possible.
                elif len(domDimCoeff) == 1 and len(paramCoeff) == 0:                    
                    dim = (domDimCoeff.keys())[0]
                    childVarSchedDim = childPart.align[dim]

                    pscale = getScale(scaleMap, parentPart, i)
                    cscale = getScale(scaleMap, childPart, dim)

                    assert Fraction(pscale).denominator == 1
                    assert Fraction(cscale).denominator == 1

                    if ((childVarSchedDim == parentVarSchedDim) and 
                        (domDimCoeff[dim] * pscale == cscale)):
                        depVec[parentVarSchedDim] = \
                                -get_constant_from_expr(arg, affine=True)
                        accessScale = pscale
                        if depVec[parentVarSchedDim] > 0:
                            depVec[parentVarSchedDim] = \
                                (int(math.ceil(depVec[parentVarSchedDim] *
                                accessScale)))
                        else:         
                            depVec[parentVarSchedDim] = \
                               (int(math.floor(depVec[parentVarSchedDim] *
                                accessScale)))
                            #print(parentPart.sched)
                            #print(childPart.sched)
                            #print(childPart.expr)
                            #print(depVec, ref)
                    else:
                        depVec[parentVarSchedDim] = '*'
                elif len(domDimCoeff) == 0 and len(paramCoeff) > 0:
                    #print(ref)
                    #assert False
                    depVec[parentVarSchedDim] = '*'
                # Only parametric or Constant access. The schedule in this 
                # dimension can be shifted to this point to uniformize the 
                # dependence.
                # In case the dimension in the parent has a constant size
                # an upper and lower bound on the dependence vector can 
                # be computed.
                elif len(domDimCoeff) + len(paramCoeff) == 0:
                    # offsets should be set here.
                    accessConstant = get_constant_from_expr(arg, affine = True)
                    parentLowerBound = parentPart.sched.domain().dim_min(i)
                    parentUpperBound = parentPart.sched.domain().dim_max(i)
                    if ((parentLowerBound.is_cst() and 
                        parentLowerBound.n_piece() == 1) and
                        (parentUpperBound.is_cst() and
                        parentUpperBound.n_piece() == 1)):
                          
                        pscale = getScale(scaleMap, parentPart, i)

                        lowVecAff = (parentLowerBound.get_pieces())[0][1]
                        val = lowVecAff.get_constant_val()
                        assert(val.get_den_si() == 1)
                        lowVec = int(math.floor((accessConstant - val.get_num_si()) *
                                             pscale))

                        highVecAff = (parentUpperBound.get_pieces())[0][1]
                        val = highVecAff.get_constant_val()
                        assert(val.get_den_si() == 1)
                        highVec = int(math.ceil((accessConstant - val.get_num_si()) *
                                             pscale))

                        if highVec == lowVec:
                            depVec[parentVarSchedDim] = highVec
                        else:
                            # Unpack dependence vectors when this hits
                            #print(ref)
                            #assert False
                            #depVec[parentVarSchedDim] = (lowVec, highVec)
                            depVec[parentVarSchedDim] = '*'
                    else:
                        depVec[parentVarSchedDim] = '*'
                else:
                    assert False
            else:
                #print(ref)
                #assert(False)
                depVec[parentVarSchedDim] = '*'

        assert depVec[0] == '-'
        depVec[0] = childPart.level_no - parentPart.level_no
        for i in range(0, dimOut):
            if (depVec[i] == '-'):
                depVec[i] = 0
        #for i in range(0, dimOut):
        #    if (depVec[i] == '-'):
        #        depVec[i] = '*'
        #        print(parentPart.sched)
        #        print(childPart.sched)
        #        print(i, ref)
        #        parentLowerBound = parentPart.sched.range().dim_min(i)
        #        parentUpperBound = parentPart.sched.range().dim_max(i)
                    
        #        childLowerBound  = childPart.sched.range().dim_min(i)
        #        childUpperBound = childPart.sched.range().dim_max(i)

        #        if (childLowerBound.is_equal(childUpperBound) and
        #            parentLowerBound.is_equal(parentUpperBound)):
        #            dimDiff = childUpperBound.sub(parentUpperBound)
        #            if (dimDiff.is_cst() and dimDiff.n_piece() == 1):
        #                aff = (dimDiff.get_pieces())[0][1]
        #                val = aff.get_constant_val()
        #                depVec[i] = (val.get_num_si())/(val.get_den_si())
        return (depVec, parentPart.level_no)

    def alignParts(self):
        """ Embeds parts whose dimension is smaller than the schedule space."""
        # Go through the parts in a sorted order and compute alignments
        compObjs = self.stage.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        # Alignments in certain cases may result in changing the relative order of
        # dimensions. This is valid only if there are no dependencies between the
        # dimensions being reordered. Dependence vectors or polyhedra can be used
        # to decide if the dimensions can be reordered.

        # For now we only realign parts which do not have self dependencies
        # Find all parts which do not have self dependencies
        noDepParts = []
        for comp in [item[0] for item in sortedCompObjs]:
            for p in self.poly_parts[comp]:
                if not self.isPartSelfDependent(p):
                    noDepParts.append(p)

        # Find an alignment map for parts which are not full dimensional
        def completeMap(newAlign, currAlign):
            for i in range(0, len(newAlign)):
                if newAlign[i] in currAlign:
                    currAlign.remove(newAlign[i])
            for i in range(0, len(newAlign)):         
                if newAlign[i] == '-':
                    newAlign[i] = currAlign.pop(0)
            return newAlign

        def compatibleAlign(align1, align2):
            compatible = True
            if len(align1) == len(align2):
                for i in range(0, len(align1)):
                    if not ((align1[i] == '-' or align2[i] == '-')
                            or (align1[i] == align2[i])):
                        compatible = False
            else:
                compatible = False
            return compatible

        def alignGroupToPart(group, part, align):
            for g in group:
                dimIn = g.sched.dim(isl._isl.dim_type.in_)
                dimOut = g.sched.dim(isl._isl.dim_type.out)
                newAlign = [ '-' for i in range(0, dimIn)]
                for i in range(0, dimIn):
                    for j in range(0, len(align)):
                        if g.align[i] == align[j]:
                            assert newAlign[i] == '-'
                            newAlign[i] = part.align[j]
                g.align = completeMap(newAlign, g.align)

        # Progressive alignment algorithm. Need to revisit when we encounter
        # stranger things. alignedPartGroups consists of groups of parts whose
        # alignments are linked. All the alignments in a group should be reordered 
        # in the same way so as to not disturb the previous alignments.
        alignedPartGroups = []
        for p in noDepParts:
            # Find if there is a unique alignment with the aligned parts
            # i.e the current part can be aligned with already aligned 
            # parts or all the aligned parts can be aligned to the current 
            # part.
            parentGroups = self.findParentGroups(p, alignedPartGroups)
            otherGroups = [ g for g in alignedPartGroups \
                            if g not in parentGroups ]
            aligns = {}

            for i in range(0, len(parentGroups)):
                aligns[i] = self.alignWithGroup(p, parentGroups[i])
            
            mergeGroups = []
            # If the alignment has alteast one valid reordering. Add the
            # group to the list of groups to be aligned and merged.
            for i in range(0, len(parentGroups)):
                addGroup = False
                for dim in aligns[i]:
                    if dim != '-':
                        addGroup = True
                if addGroup:
                    mergeGroups.append(i)

            mergedGroup = []
            for i in mergeGroups:
                alignGroupToPart(parentGroups[i], p, aligns[i])
                mergedGroup = mergedGroup + parentGroups[i]
            parentGroups = [i for j, i in enumerate(parentGroups)\
                            if j not in mergeGroups ]
            
            mergedGroup = mergedGroup + [p]
            parentGroups.append(mergedGroup)
            alignedPartGroups = parentGroups + otherGroups
         
        # The alignment procedure above can move the fast varying
        # dimension outside. This has to be fixed.


    def isCompSelfDependent(self, comp):
        parts = self.poly_parts[comp]
        for p in parts:
            if self.isPartSelfDependent(p):
                return True
        return False

    def isPartSelfDependent(self, part):
        refs    = part.getPartRefs()
        objRefs = [ ref.objectRef for ref in refs\
                         if ref.objectRef == part.comp]
        if len(objRefs) > 0:
            return True
        return False

    def isParent(self, part1, part2): 
        refs    = part2.getPartRefs()
        objRefs = [ ref.objectRef for ref in refs\
                         if ref.objectRef == part1.comp]
        if len(objRefs) > 0:
            return True
        return False

    def getDomainDimCoeffs(self, sched, arg):
        domDimCoeff = {}
        if (isAffine(arg)):
            coeff = get_affine_var_and_param_coeff(arg)
            for item in coeff:
                if type(item) == Variable:
                    dim = sched.find_dim_by_name(isl._isl.dim_type.in_,
                                                 item.name)
                    domDimCoeff[dim] = coeff[item]
        return domDimCoeff

    def getParamCoeffs(self, sched, arg):
        paramCoeff = {}
        if (isAffine(arg)):
            coeff = get_affine_var_and_param_coeff(arg)
            for item in coeff:
                if type(item) == Parameter:
                    dim = sched.find_dim_by_name(isl._isl.dim_type.param,
                                                 item.name)
                    paramCoeff[dim] == coeff[item]
        return paramCoeff

    def computeRelativeScalingFactors(self, parentPart, childPart):
        """Computes the relative scaling factor in each dimension necessary to 
           uniformize dependencies between two parts. If the dependencies are 
           already uniform the scaling factor for the dimension is set to 1. If
           the dependencies cannot be uniformized along a particular dimension 
           the scaling factor for that dimension is set to *."""
        
        assert self.isParent(parentPart, childPart)

        # Scaling and offset factors structure.
        # ([scale factors], [offsets])
        #
        # [scale factors] determine the amount by which the dimension has to be
        # scaled to uniformize the dependencies. Each reference to the parent is
        # considered while determing the scaling factors. All the references 
        # should have the same scaling factor in a particular dimension otherwise
        # the scaling factor for the dimension cannot be determined uniquely.
        #
        # [offsets] specify the shift in each dimension that is require to 
        # uniformize dependencies. Simliar to dimension scale factors all 
        # offsets for a dimension should agree otherwise the offset for the 
        # dimension cannot be determied uniquely.
        refs       = childPart.getPartRefs()
        # Filtering out self references
        parentRefs = [ ref for ref in refs \
                       if ref.objectRef == parentPart.comp] 

        dimIn = childPart.sched.dim(isl._isl.dim_type.in_)
        scale = [ '-' for i in range(0, dimIn) ]
        offset = [ '-' for i in range(0, dimIn) ]

        def findDimScheduledTo(part, schedDim):
            for i in range(0, len(part.align)):
                if part.align[i] == schedDim:
                    return i
            return -1

        for ref in parentRefs:
            numArgs = len(ref.arguments)
            for i in range(0, numArgs):
                arg = ref.arguments[i]
                parentVarSchedDim = parentPart.align[i]            
                if (isAffine(arg)):
                    domDimCoeff = self.getDomainDimCoeffs(childPart.sched, arg)
                    paramCoeff = self.getParamCoeffs(childPart.sched, arg)
                    # Parameter coefficents can also be considered to
                    # generate parametric shifts. Yet to be seen.
                    
                    # Indexed with multiple variables.
                    if (len(domDimCoeff) > 1 or 
                        (len(domDimCoeff) == 1 and len(paramCoeff) >=1)):
                        # Although there are multiple coefficients. If 
                        # there is only one variable coefficient and other
                        # parametric coefficients. Uniformization can be 
                        # done with parametric shifts. Full affine scheduling 
                        # might be able to find a way to uniformize 
                        # dependencies. This has to be further explored.
                        dim = findDimScheduledTo(childPart, parentVarSchedDim)
                        if dim != -1:
                            scale[dim] = '*'
                    # Indexed with a single variable. This can either be an 
                    # uniform access or can be uniformized with scaling when 
                    # possible.
                    elif len(domDimCoeff) == 1 and len(paramCoeff) == 0:
                        dim = (domDimCoeff.keys())[0]
                        childVarSchedDim = childPart.align[dim]
                        # Checking if the schedule dimensions match only 
                        # then can the dependence be uniformized.
                        if childVarSchedDim != parentVarSchedDim:
                            continue
                        if scale[dim] == '-':
                            scale[dim] = domDimCoeff[dim]
                        elif scale[dim] != domDimCoeff[dim]:
                            scale[dim] = '*'
                    elif len(domDimCoeff) == 0 and len(paramCoeff) > 0:
                        continue
                    # Only parametric or Constant access. The schedule in this 
                    # dimension can be shifted to this point to uniformize the 
                    # dependence.
                    # In case the dimension in the parent has a constant size
                    # an upper and lower bound on the dependence vector can 
                    # be computed.
                    elif len(domDimCoeff) + len(paramCoeff) == 0:
                        # offsets should be set here
                        continue
                else:
                    dim = findDimScheduledTo(childPart, parentVarSchedDim)
                    if dim != -1:
                        scale[dim] = '*'

        for i in range(0, dimIn):
            if scale[i] == '-':
                scale[i] = 1
            if offset[i] == '-':
                offset[i] = 0
        return (scale, offset)

    def createGroupScaleMap(self, group):
        groupScaleMap = {}
        for part in group:
            groupScaleMap[part] = list(part.scale)
        return groupScaleMap

    def isStencilGroup(self, group):
        depVecsGroup = self.getGroupDependenceVectors(group)
        for depVec, h in depVecsGroup:
            if '*' in depVec:
                return False
        return True    

    def getPartSize(self, part, param_estimates):
        size = None
        domain = part.comp.domain
        if isinstance(part.comp, Accumulator):
            domain = part.comp.reductionDomain
        for interval in domain:
            subsSize = self.getDimSize(interval, param_estimates)
            if isConstantExpr(subsSize):
                if size is None:
                    size = get_constant_from_expr(subsSize)
                else:
                    size = size * get_constant_from_expr(subsSize)
            else:
                size = '*'
                break
        assert size is not None
        return size
    
    def getDimSize(self, interval, param_estimates):
        paramValMap = {}
        for est in param_estimates:
            assert isinstance(est[0], Parameter)
            paramValMap[est[0]] = Value.numericToValue(est[1])

        dimSize = interval.upperBound - interval.lowerBound + 1
        return substituteVars(dimSize, paramValMap)

       
    def getVarName(self):
        name = "_i" + str(self._var_count)
        self._var_count+=1
        return name

    def getFuncName(cls):
        name = "_f" + str(self._func_count)
        self._func_count+=1
        return name

    def __str__(self):
        polystr = ""
        for comp in self.poly_parts:
            for part in self.poly_parts[comp]:
                polystr = polystr + part.__str__() + '\n'

        if (self.polyast != []):
            for ast in self.polyast:
                printer = isl.Printer.to_str(self.ctx)
                printer = printer.set_output_format(isl.format.C)
                printOpts = isl.AstPrintOptions.alloc(self.ctx) 
                printer = ast.print_(printer, printOpts)
                aststr = printer.get_str()
                polystr = polystr + '\n' + aststr
        return polystr


def map_coeff_to_dim(coeff):
    variables = list(coeff.keys())
    for var in variables:
        coeffval = coeff[var]
        coeff.pop(var)
        if (isinstance(var, Parameter)):
            coeff[('param', var.name)] = coeffval
        elif (isinstance(var, Variable)):
            coeff[('in', var.name)] = coeffval
    return coeff

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
    level_coeff[('constant', 0)] = level_no
    eq_coeff.append(level_coeff)

    # Setting the remaining dimensions to zero
    for i in range(1, dim_out):
        if not dim_set[i]:
            coeff = {}
            coeff[('out', i)] = 1
            coeff[('constant', 0)] = 0
            eq_coeff.append(coeff)
    return [ineq_coeff, eq_coeff]

def format_domain_constraints(domain, var_names):
    ineq_coeff = []
    eq_coeff   = []
    dom_len = len(domain)
    for i in range(0, dom_len):
        coeff = {}
        interval = domain[i]

        lb_coeff = get_affine_var_and_param_coeff(interval.lowerBound)
        # Mapping from variable names to the corresponding dimension
        lb_coeff = map_coeff_to_dim(lb_coeff)
        lb_const = get_constant_from_expr(interval.lowerBound, affine = True)

        # Normalizing into >= format
        coeff = dict( (n, -lb_coeff.get(n)) for n in lb_coeff )
        coeff[('constant', 0)] = -lb_const
        coeff[('in', var_names[i])] = 1
        ineq_coeff.append(coeff)

        ub_coeff = get_affine_var_and_param_coeff(interval.upperBound)
        # Mapping from variable names to the corresponding dimension
        ub_coeff = map_coeff_to_dim(ub_coeff)
        ub_const = get_constant_from_expr(interval.upperBound, affine = True)

        # Normalizing into >= format
        coeff = ub_coeff
        coeff[('constant', 0)] = ub_const
        coeff[('in', var_names[i])] = -1
        ineq_coeff.append(coeff)

    return [ineq_coeff, eq_coeff]

def format_conjunct_constraints(conjunct):
    # TODO check if the condition is a conjunction
    # print([ cond.__str__() for cond in conjunct ])
    ineq_coeff = []
    eq_coeff = []
    for cond in conjunct:
        coeff = {}
        left_coeff = get_affine_var_and_param_coeff(cond.lhs)
        right_coeff = get_affine_var_and_param_coeff(cond.rhs)
        left_const = get_constant_from_expr(cond.lhs, affine = True)
        right_const = get_constant_from_expr(cond.rhs, affine = True)

        # Mapping from variable names to the corresponding dimension
        left_coeff = map_coeff_to_dim(left_coeff)
        right_coeff = map_coeff_to_dim(right_coeff)
        
        def constant_div_factor(const):
            m = 1
            for coeff in const:
                if isinstance(const[coeff], Fraction):
                    m = (abs(const[coeff].denominator) * m) // \
                        gcd(abs(const[coeff].denominator), m)
            assert m.denominator == 1
            m = m.numerator
            return m

        # Normalizing >= format
        if (cond.conditional in ['<=','<']):
            coeff = dict( (n, -left_coeff.get(n, 0) + right_coeff.get(n, 0)) \
                          for n in set(left_coeff) | set(right_coeff) )
            d = constant_div_factor(coeff)
            coeff[('constant', 0)] = -left_const + right_const - \
                                     int(cond.conditional == '<') - \
                                     Fraction(d-1, d)
            ineq_coeff.append(coeff)
        elif(cond.conditional in ['>=','>']):
            coeff = dict( (n, left_coeff.get(n, 0) - right_coeff.get(n, 0)) \
                           for n in set(left_coeff) | set(right_coeff) )
            d = constant_div_factor(coeff)
            coeff[('constant', 0)] = left_const - right_const - \
                                     int(cond.conditional == '>') + \
                                     Fraction(d-1, d)
            ineq_coeff.append(coeff)
        else:
            # Weird
            assert(cond.conditional == '==')
            coeff = dict( (n, left_coeff.get(n, 0) - right_coeff.get(n, 0)) \
                          for n in set(left_coeff) | set(right_coeff) )
            coeff[('constant', 0)] = left_const - right_const
            eq_coeff.append(coeff)

    return [ineq_coeff, eq_coeff]

# NOTE: Dead code?
def getParamsInvolved(sched, dim):
    numParams = sched.domain().n_param()
    paramNames = [ ]
    constraints = sched.domain().get_constraints()
    dimName = sched.domain().get_dim_name(isl._isl.dim_type.set, dim)
    for p in range(0, numParams): 
        for const in constraints:
           pname = sched.domain().get_dim_name(isl._isl.dim_type.param, p)
           dimId  = const.get_space().find_dim_by_name(isl._isl.dim_type.set, dimName)
           pId = const.get_space().find_dim_by_name(isl._isl.dim_type.param, pname)
           dimInvolved = const.involves_dims(isl._isl.dim_type.set, dimId, 1)
           paramInvolved = const.involves_dims(isl._isl.dim_type.param, pId, 1)
           if (dimInvolved and paramInvolved):
                   if pname not in paramNames:
                       paramNames.append(pname)
    return paramNameis

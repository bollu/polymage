from __future__ import absolute_import, division, print_function

from collections import OrderedDict

import time
import pipe
from constructs import *
from expression import *
from schedule import *
from poly import *
import islpy as isl
import expr_ast as expr
import targetc as genc
import logging

# LOG CONFIG #
codegen_logger = logging.getLogger("codegen.py")
codegen_logger.setLevel(logging.INFO)
LOG = codegen_logger.log

# short hands for methods
new_temp = genc.CNameGen.get_temp_var_name
new_iter = genc.CNameGen.get_iterator_name

def isl_expr_to_cgen(expr, prologue_stmts = None):

    prolog = prologue_stmts
    if expr.get_type() == isl._isl.ast_expr_type.op:
        # short hand
        op_typ = expr.get_op_type()

        if (op_typ == isl._isl.ast_op_type.access or  # 23
            op_typ == isl._isl.ast_op_type.call or  # 22
            op_typ == isl._isl.ast_op_type.cond or  # 15
            op_typ == isl._isl.ast_op_type.member or  # 24
            op_typ == isl._isl.ast_op_type.select):  # 16
            assert ("Unexpected isl ast_expr op_type:"+str(op_typ) \
                    and False)
        if op_typ == isl._isl.ast_op_type.min:
            num_args = expr.get_op_n_arg()
            cmin = genc.CMin(isl_expr_to_cgen(expr.get_op_arg(0), prolog),
                             isl_expr_to_cgen(expr.get_op_arg(1), prolog))
            for i in range(2, num_args):
                cmin = genc.CMin(cmin,
                                 isl_expr_to_cgen(expr.get_op_arg(i), prolog))
            if prolog is not None:
                var_type = genc.TypeMap.convert(Int)
                min_cvar = genc.CVariable(var_type, new_temp())
                decl = genc.CDeclaration(var_type, min_cvar, cmin)
                prolog.append(decl)
                cmin = min_cvar
            return cmin
        if op_typ == isl._isl.ast_op_type.max:
            num_args = expr.get_op_n_arg()
            cmax = genc.CMax(isl_expr_to_cgen(expr.get_op_arg(0), prolog),
                             isl_expr_to_cgen(expr.get_op_arg(1), prolog))
            for i in range(2, num_args):
                cmax = genc.CMax(cmax,
                                 isl_expr_to_cgen(expr.get_op_arg(i), prolog))
            if prolog is not None:
                var_type = genc.TypeMap.convert(Int)
                max_cvar = genc.CVariable(var_type, new_temp())
                decl = genc.CDeclaration(var_type, max_cvar, cmax)
                prolog.append(decl)
                cmax = max_cvar
            return cmax
        if op_typ == isl._isl.ast_op_type.fdiv_q:
            assert("Division must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return \
                genc.c_macro_floord(isl_expr_to_cgen(expr.get_op_arg(0), prolog),
                                  isl_expr_to_cgen(expr.get_op_arg(1), prolog))
        if op_typ == isl._isl.ast_op_type.add:
            assert("Addition must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) + \
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)
        if op_typ == isl._isl.ast_op_type.mul:
            assert("Multiplication must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) * \
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)
        if (op_typ == isl._isl.ast_op_type.div or
            op_typ == isl._isl.ast_op_type.pdiv_q):
            assert("Division must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) / \
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)
        if op_typ == isl._isl.ast_op_type.pdiv_r:
            assert("Division must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) % \
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)
        if op_typ == isl._isl.ast_op_type.sub:
            assert("Subtraction must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) - \
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)
        if op_typ == isl._isl.ast_op_type.minus:
            assert("Unary minus must have exactly 1 argument!" and \
                   expr.get_op_n_arg() == 1)
            return -isl_expr_to_cgen(expr.get_op_arg(0), prolog)
        if (op_typ == isl._isl.ast_op_type.and_ or
            op_typ == isl._isl.ast_op_type.and_then):
            assert("Logical And must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) & \
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)
        if (op_typ == isl._isl.ast_op_type.or_ or
            op_typ == isl._isl.ast_op_type.or_else):
            assert("Logical Or must have exactly 2 arguments!" and \
                   expr.get_op_n_arg() == 2)
            return isl_expr_to_cgen(expr.get_op_arg(0), prolog) |\
                   isl_expr_to_cgen(expr.get_op_arg(1), prolog)

    if expr.get_type() == isl._isl.ast_expr_type.int:
        return expr.get_val().to_python()
    if expr.get_type() == isl._isl.ast_expr_type.id:
        return genc.CVariable(genc.c_int, expr.get_id().get_name())

def isl_cond_to_cgen(cond, prologue_stmts = None):
    comp_dict = { isl._isl.ast_op_type.eq: '==',
                  isl._isl.ast_op_type.ge: '>=',
                  isl._isl.ast_op_type.gt: '>',
                  isl._isl.ast_op_type.le: '<=',
                  isl._isl.ast_op_type.lt: '<' ,
                  isl._isl.ast_op_type.and_: '&&',
                  isl._isl.ast_op_type.and_then: '&&',
                  isl._isl.ast_op_type.or_: '||',
                  isl._isl.ast_op_type.or_else: '||' }

    assert("Unsupported isl conditional type:"+str(cond.get_op_type())+\
           " of condition:"+str(cond) and \
           cond.get_op_type() in comp_dict)
    comp_op = comp_dict[cond.get_op_type()]
    assert("Conditional expressions must have exactly 2 arguments!" and \
           cond.get_op_n_arg() == 2)
    if (comp_op == '&&' or comp_op == '||'):
        left = isl_cond_to_cgen(cond.get_op_arg(0), prologue_stmts)
        right = isl_cond_to_cgen(cond.get_op_arg(1), prologue_stmts)
    else:
        left = isl_expr_to_cgen(cond.get_op_arg(0), prologue_stmts)
        right = isl_expr_to_cgen(cond.get_op_arg(1), prologue_stmts)

    return genc.CCond(left, comp_op, right)

def is_inner_most_parallel(node):
    # Check if there are no parallel loops within the node
    no_inner_parallel = True
    if node.get_type() == isl._isl.ast_node_type.block:
        num_nodes = (node.block_get_children().n_ast_node())
        for i in range(0, num_nodes):
            child =  node.block_get_children().get_ast_node(i)
            no_inner_parallel = no_inner_parallel and \
                                is_inner_most_parallel(child)
    else:
        if node.get_type() == isl._isl.ast_node_type.for_:
            user_nodes = get_user_nodes_in_body(node.for_get_body())
            var = isl_expr_to_cgen(node.for_get_iterator())
            if is_sched_dim_parallel(user_nodes, var.name):
                no_inner_parallel = False
        elif node.get_type() == isl._isl.ast_node_type.if_:
            no_inner_parallel = no_inner_parallel and \
                              is_inner_most_parallel(node.if_get_then())
            has_else = node.if_has_else()
            if has_else:
              no_inner_parallel = no_inner_parallel and \
                                is_inner_most_parallel(node.if_get_else())
        elif node.get_type() == isl._isl.ast_node_type.user:
            pass
        else:
            assert("Unexpected isl ast node type:"+str(node.get_type()) and \
                   False)
    return no_inner_parallel

def get_user_nodes_in_body(body):
    user_nodes = []
    if body.get_type() == isl._isl.ast_node_type.block:
        num_nodes = (body.block_get_children().n_ast_node())
        for i in range(0, num_nodes):
            child = body.block_get_children().get_ast_node(i)
            user_nodes += get_user_nodes_in_body(child)
    else:
        if body.get_type() == isl._isl.ast_node_type.for_:
            user_nodes += get_user_nodes_in_body(body.for_get_body())
        elif body.get_type() == isl._isl.ast_node_type.if_:
            user_nodes += get_user_nodes_in_body(body.if_get_then())
            if body.if_has_else():
              user_nodes += get_user_nodes_in_body(body.if_get_else())
        elif body.get_type() == isl._isl.ast_node_type.user:
            user_nodes += [body]
        else:
            assert("Unexpected isl node type:"+str(node.get_type()) and \
                   False)
    return user_nodes

def is_sched_dim_parallel(polyrep, user_nodes, sched_dim_name):
    is_parallel = True
    for node in user_nodes:
        part_id = node.user_get_expr().get_op_arg(0).get_id()
        part = isl_get_id_user(part_id)
        if sched_dim_name not in part.parallel_sched_dims:
            is_parallel = False
    return is_parallel

def is_sched_dim_vector(polyrep, user_nodes, sched_dim_name):
    is_vector = True
    for node in user_nodes:
        part_id = node.user_get_expr().get_op_arg(0).get_id()
        part = isl_get_id_user(part_id)
        if sched_dim_name not in part.vector_sched_dim:
            is_vector = False
    return is_vector

def get_arrays_for_user_nodes(pipe, polyrep, user_nodes, cfunc_map):
    arrays = []
    for node in user_nodes:
        part_id = node.user_get_expr().get_op_arg(0).get_id()
        part = isl_get_id_user(part_id)
        array, scratch = cfunc_map[part.comp.func]
        if (array not in arrays) and (True in scratch):
            arrays.append(array)
    return arrays

def cvariables_from_variables_and_sched(node, variables, sched):
    cvar_map = {}
    for i in range(0, len(variables)):
        var_name = variables[i].name
        dim = sched.find_dim_by_name(isl._isl.dim_type.in_,
                                     var_name)
        cvar_map[variables[i]] = \
            isl_expr_to_cgen(node.user_get_expr().get_op_arg(dim+1))
    return cvar_map

# TESTME
def generate_c_naive_from_accumlate_node(pipe, polyrep, node, body,
                                         cparam_map, cfunc_map):
    part_id = node.user_get_expr().get_op_arg(0).get_id()
    poly_part = isl_get_id_user(part_id)
    dom_len = len(poly_part.comp.func.reductionVariables)
    # FIXME: this has to be changed similar to Expression Node
    cvar_map = cvariables_from_variables_and_sched(node,
                 poly_part.comp.func.reductionVariables, poly_part.sched)
    expr = generate_c_expr(pipe, poly_part.expr.expression,
                           cparam_map, cvar_map, cfunc_map)
    prologue = []
    array_ref = generate_c_expr(pipe, poly_part.expr.accumulate_ref,
                                cparam_map, cvar_map, cfunc_map,
                                prologue_stmts = prologue)
    assign = genc.CAssign(array_ref, array_ref + expr)

    if prologue is not None:
        for s in prologue:
            body.add(s)

    if poly_part.pred:
        ccond = generate_c_cond(pipe, poly_part.pred,
                                cparam_map, cvar_map, cfunc_map)
        cif = genc.CIfThen(ccond)
        with cif.if_block as ifblock:
            ifblock.add(assign)
        body.add(cif)
    else:
        body.add(assign)

def generate_c_naive_from_expression_node(pipe, polyrep, node, body,
                                          cparam_map, cfunc_map):
    part_id = node.user_get_expr().get_op_arg(0).get_id()
    poly_part = isl_get_id_user(part_id)
    variables = poly_part.comp.func.variables
    dom_len = len(variables)

    # Get the mapping to the array
    array, scratch = cfunc_map[poly_part.comp.func]

    acc_scratch = [ False for i in range(0, dom_len) ]
    for i in range(0, dom_len):
        if i in poly_part.dim_tile_info:
            if (poly_part.dim_tile_info[i][0] != 'none'):
                acc_scratch[i] = True

    cvar_map = \
        cvariables_from_variables_and_sched(node, variables,
                                            poly_part.sched)
    arglist = []
    scratch_map = {}
    for i in range(0, dom_len):
        acc_expr = variables[i] - \
                  poly_part.comp.func.domain[i].lowerBound
        if acc_scratch[i]:
            var_name = variables[i].name
            #dim = \
            #    poly_part.sched.find_dim_by_name(isl._isl.dim_type.in_,
            #                                     '_Acc_' + var_name)
            #dim_rem = \
            #    poly_part.sched.find_dim_by_name(isl._isl.dim_type.in_,
            #                                     '_Rem_' + var_name)
            mul_rem = \
                poly_part.sched.find_dim_by_name(isl._isl.dim_type.in_,
                                                 '_Mul_' + var_name)
            #org_var = Variable(Int, '_Acc_' + var_name)
            #rem_var = Variable(Int, '_Rem_' + var_name)
            #cvar_map[org_var] = \
            #    isl_expr_to_cgen(node.user_get_expr().get_op_arg(dim+1))
            #cvar_map[rem_var] = \
            #   isl_expr_to_cgen(node.user_get_expr().get_op_arg(dim_rem+1))
            mul_var = Variable(Int, '_Mul_' + var_name)
            cvar_map[mul_var] = \
                isl_expr_to_cgen(node.user_get_expr().get_op_arg(mul_rem+1))
            scratch_map[variables[i]] = (mul_var)
            if scratch[i]:
                acc_expr = (mul_var)
        arglist.append(generate_c_expr(pipe, acc_expr,
                                       cparam_map, cvar_map, cfunc_map,
                                       scratch_map))
    prologue = []
    expr = generate_c_expr(pipe, poly_part.expr, cparam_map, cvar_map, cfunc_map,
                           scratch_map, prologue_stmts = prologue)
    assign = genc.CAssign(array(*arglist), expr)

    if prologue is not None:
        for s in prologue:
            body.add(s)

    if poly_part.pred:
        ccond = generate_c_cond(pipe, poly_part.pred,
                                cparam_map, cvar_map, cfunc_map, scratch_map)
        cif = genc.CIfThen(ccond)
        with cif.if_block as ifblock:
            ifblock.add(assign)
        body.add(cif)
        #var = genc.CVariable(genc.c_int, "_c_" + poly_part.comp.func.name)
        #incr = genc.CAssign(var, var + 1)
        #body.add(incr)
    else:
        body.add(assign)
        #var = genc.CVariable(genc.c_int, "_c_" + poly_part.comp.func.name)
        #incr = genc.CAssign(var, var + 1)
        #body.add(inc)

def generate_c_naive_from_isl_ast(pipe, polyrep, node, body,
                                  cparam_map, cfunc_map, pooled):
    if node.get_type() == isl._isl.ast_node_type.block:
        num_nodes = (node.block_get_children().n_ast_node())
        for i in range(0, num_nodes):
            child = node.block_get_children().get_ast_node(i)
            generate_c_naive_from_isl_ast(pipe, polyrep, child, body,
                                          cparam_map, cfunc_map, pooled)
    else:
        if node.get_type() == isl._isl.ast_node_type.for_:
            # Convert lb and ub expressions to C expressions
            prologue = []
            cond = isl_cond_to_cgen(node.for_get_cond(), prologue)
            var = isl_expr_to_cgen(node.for_get_iterator())
            var_inc = isl_expr_to_cgen(node.for_get_inc())
            incr = genc.CAssign(var, var+var_inc)
            if prologue is not None:
                for s in prologue:
                    body.add(s)

            prologue = []
            init = isl_expr_to_cgen(node.for_get_init(), prologue)
            if prologue is not None:
                for s in prologue:
                    body.add(s)
            var_decl =  genc.CDeclaration(var.typ, var, init)
            loop = genc.CFor(var_decl, cond, incr)

            # Check if the loop is a parallel or a vector dimension by
            # examining the loop body.
            user_nodes = get_user_nodes_in_body(node.for_get_body())

            dim_parallel = is_sched_dim_parallel(polyrep, user_nodes, var.name)
            dim_vector = is_sched_dim_vector(polyrep, user_nodes, var.name)
            arrays = get_arrays_for_user_nodes(pipe, polyrep, user_nodes, cfunc_map)

            if dim_parallel:
                omp_pragma = genc.CPragma("omp parallel for schedule(static)")
                body.add(omp_pragma)
            if dim_vector:
                vec_pragma = genc.CPragma("ivdep")
                body.add(vec_pragma)

            #flat_scratch = 'flatten_scratchpad' in pipeline._options
            flat_scratch = True

            body.add(loop)
            # Assuming only one parallel dimension and a whole lot
            # of other things.
            freelist = []
            if dim_parallel:
                with loop.body as lbody:
                    for array in arrays:
                        #if array.is_constant_size() and False:
                        #if array.is_constant_size():
                        if array.is_constant_size() or True:
                            array_decl = genc.CArrayDecl(array, flat_scratch)
                            lbody.add(array_decl)
                        else:
                            array_ptr = genc.CPointer(array.typ, 1)
                            array_decl = genc.CDeclaration(array_ptr, array)
                            lbody.add(array_decl)
                            array.allocate_contiguous(lbody)
                            freelist.append(array)
            with loop.body as lbody:
                generate_c_naive_from_isl_ast(pipe, polyrep,
                                              node.for_get_body(),
                                              lbody, cparam_map, cfunc_map,
                                              pooled)
                # Deallocate storage
                for array in freelist:
                    array.deallocate(lbody, pooled)

        if node.get_type() == isl._isl.ast_node_type.if_:
            if_cond = isl_cond_to_cgen(node.if_get_cond())
            if node.if_has_else():
                cif_else = genc.CIfThenElse(if_cond)
                with cif_else.if_block as ifblock:
                    generate_c_naive_from_isl_ast(pipe, polyrep,
                                                  node.if_get_then(), ifblock,
                                                  cparam_map, cfunc_map, pooled)
                with cif_else.else_block as elseblock:
                    generate_c_naive_from_isl_ast(pipe, polyrep,
                                                  node.if_get_else(),
                                                  elseblock,
                                                  cparam_map, cfunc_map, pooled)
                body.add(cif_else)
            else:
                cif = genc.CIfThen(if_cond)
                with cif.if_block as ifblock:
                    generate_c_naive_from_isl_ast(pipe,
                                                  polyrep, node.if_get_then(),
                                                  ifblock,
                                                  cparam_map, cfunc_map, pooled)
                body.add(cif)

        if node.get_type() == isl._isl.ast_node_type.user:
            # The first argument is the computation object.
            # Retrieving the polyPart.
            part_id = node.user_get_expr().get_op_arg(0).get_id()
            poly_part = isl_get_id_user(part_id)
            if isinstance(poly_part.expr, Reduce):
                generate_c_naive_from_accumlate_node(pipe, polyrep, node, body,
                                                     cparam_map, cfunc_map)
            elif isinstance(poly_part.expr, AbstractExpression):
                generate_c_naive_from_expression_node(pipe, polyrep, node, body,
                                                      cparam_map, cfunc_map)
            else:
                assert ("Invalid pipeline stage type:"+str(poly_part.expr) \
                        and False)

def generate_c_cond(pipe, cond,
                    cparam_map, cvar_map, cfunc_map, 
                    scratch_map = {}, prologue_stmts = None):
    if (cond.conditional in ['&&', '||']):
        left_cond = generate_c_cond(pipe, cond.lhs,
                                    cparam_map, cvar_map, cfunc_map,
                                    scratch_map)
        right_cond = generate_c_cond(pipe, cond.rhs,
                                     cparam_map, cvar_map, cfunc_map,
                                     scratch_map, prologue_stmts)
        return genc.CCond(left_cond, cond.conditional, right_cond)
    elif (cond.conditional in ['<', '<=', '>', '>=', '==', '!=']):
        left_cond = generate_c_expr(pipe, cond.lhs,
                                   cparam_map, cvar_map, cfunc_map,
                                   scratch_map, prologue_stmts)
        right_cond = generate_c_expr(pipe, cond.rhs,
                                    cparam_map, cvar_map, cfunc_map,
                                    scratch_map, prologue_stmts)
        return genc.CCond(left_cond, cond.conditional, right_cond)
    assert("Unsupported or invalid conditional:"+str(cond.conditional) and \
           False)

def generate_c_expr(pipe, exp, cparam_map, cvar_map, cfunc_map,
                    scratch_map = {}, prologue_stmts = None):
    if isinstance(exp, AbstractBinaryOpNode):
        left = generate_c_expr(pipe, exp.left, cparam_map, cvar_map,
                               cfunc_map, scratch_map, prologue_stmts)
        right = generate_c_expr(pipe, exp.right, cparam_map, cvar_map,
                                cfunc_map, scratch_map, prologue_stmts)
        return AbstractBinaryOpNode(left, right, exp.op)
    if isinstance(exp, AbstractUnaryOpNode):
        child = generate_c_expr(pipe, exp.child, cparam_map, cvar_map,
                                cfunc_map, scratch_map, prologue_stmts)
        return AbstractUnaryOpNode(child, exp.op)
    if isinstance(exp, Value):
        return genc.CValue(exp.value, exp.typ)
    if isinstance(exp, Variable):
        if isinstance(exp, Parameter):
            return cparam_map[exp]
        return cvar_map[exp]
    if isinstance(exp, Reference):
        array, scratch = cfunc_map[exp.objectRef]
        num_args = len(exp.objectRef.domain)
        shifted_args = []
        for i in range(num_args):
            scratch_arg = (exp.arguments[i] -
                           exp.objectRef.domain[i].lowerBound)
            if scratch and scratch[i]:
                scratch_arg = substitute_vars(exp.arguments[i], scratch_map)
            shifted_args.append(simplify_expr(scratch_arg))
        args = [ generate_c_expr(pipe, arg, cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
                 for arg in shifted_args ]
        return array(*args)
    if isinstance(exp, Select):
        c_cond = generate_c_cond(pipe, exp.condition,
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        true_expr = generate_c_expr(pipe, exp.trueExpression,
                                    cparam_map, cvar_map, cfunc_map,
                                    scratch_map, prologue_stmts)
        false_expr = generate_c_expr(pipe, exp.falseExpression,
                                     cparam_map, cvar_map, cfunc_map,
                                     scratch_map, prologue_stmts)
        if prologue_stmts is not None:
            # selects are executed with both true and false paths computed
            # before the condition check

            # true path
            var_type = genc.TypeMap.convert(getType(exp.trueExpression))
            true_c_var = genc.CVariable(var_type, new_temp())
            decl = genc.CDeclaration(var_type, true_c_var, true_expr)
            prologue_stmts.append(decl)

            # false path
            var_type = genc.TypeMap.convert(getType(exp.falseExpression))
            false_c_var = genc.CVariable(var_type, new_temp())
            decl = genc.CDeclaration(var_type, false_c_var, false_expr)
            prologue_stmts.append(decl)

            var_type = genc.TypeMap.convert(getType(exp))
            sel_c_var = genc.CVariable(var_type, new_temp())
            decl = genc.CDeclaration(var_type, sel_c_var,
                                     genc.CSelect(c_cond,
                                                  true_c_var,
                                                  false_c_var))
            prologue_stmts.append(decl)

            return sel_c_var

        return genc.CSelect(c_cond, true_expr, false_expr)
    if isinstance(exp, Max):
        cexpr1 = generate_c_expr(pipe, exp.arguments[0],
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        cexpr2 = generate_c_expr(pipe, exp.arguments[1],
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        return genc.CMax(cexpr1, cexpr2)
    if isinstance(exp, Min):
        cexpr1 = generate_c_expr(pipe, exp.arguments[0],
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        cexpr2 = generate_c_expr(pipe, exp.arguments[1],
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        return genc.CMin(cexpr1, cexpr2)
    if isinstance(exp, Pow):
        cexpr1 = generate_c_expr(pipe, exp.arguments[0],
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        cexpr2 = generate_c_expr(pipe, exp.arguments[1],
                                 cparam_map, cvar_map, cfunc_map,
                                 scratch_map, prologue_stmts)
        return genc.CPow(cexpr1, cexpr2)
    if isinstance(exp, Powf):
        cexpr1 = generate_c_expr(pipe, exp.arguments[0],
                                 cparam_map, cvar_map, cfunc_map)
        cexpr2 = generate_c_expr(pipe, exp.arguments[1],
                                 cparam_map, cvar_map, cfunc_map)
        return genc.CPowf(cexpr1, cexpr2)
    if isinstance(exp, Exp):
        cexpr = generate_c_expr(pipe, exp.arguments[0],
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CExp(cexpr)
    if isinstance(exp, Sqrt):
        cexpr = generate_c_expr(pipe, exp.arguments[0],
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CSqrt(cexpr)
    if isinstance(exp, Sqrtf):
        cexpr = generate_c_expr(pipe, exp.arguments[0],
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CSqrtf(cexpr)
    if isinstance(exp, Sin):
        cexpr = generate_c_expr(pipe, exp.arguments[0],
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CSin(cexpr)
    if isinstance(exp, Cos):
        cexpr = generate_c_expr(pipe, exp.arguments[0],
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CCos(cexpr)
    if isinstance(exp, Abs):
        cexpr = generate_c_expr(pipe, exp.arguments[0],
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CAbs(cexpr)
    if isinstance(exp, Cast):
        cexpr = generate_c_expr(pipe, exp.expression,
                                cparam_map, cvar_map, cfunc_map,
                                scratch_map, prologue_stmts)
        return genc.CCast(genc.TypeMap.convert(exp.typ), cexpr)
    raise TypeError(type(exp))

def create_loop_variables(group, variables):
    cvar_map = {}
    # Create iterator variables and bind them to DSL variables
    for i in range(0, len(variables)):
        var_type = genc.TypeMap.convert(variables[i].typ)
        var = genc.CVariable(var_type, new_iter())
        # Bind variables to C iterators
        cvar_map[variables[i]] = var

    return cvar_map

def create_perfect_nested_loop(pipe, group, pipe_body, variables, domains,
                               cparam_map, cvar_map, cfunc_map):
    lbody = pipe_body
    for i in range(0, len(variables)):
        var = cvar_map[variables[i]]
        # Convert lb and ub expressions to C expressions
        lb = generate_c_expr(pipe, domains[i].lowerBound,
                             cparam_map, cvar_map, cfunc_map)
        ub = generate_c_expr(pipe, domains[i].upperBound,
                             cparam_map, cvar_map, cfunc_map)

        var_decl = genc.CDeclaration(var.typ, var, lb)
        comp_op = '<='
        cond = genc.CCond(var, comp_op, ub)
        incr = genc.CAssign(var, var+1)

        loop = genc.CFor(var_decl, cond, incr)
        lbody.add(loop, False)
        lbody = loop.body

    return lbody

def generate_function_scan_loops(pipe, group, comp, pipe_body,
                                 cparam_map, cfunc_map):
    """
    generates code for Function class
    """
    func = comp.func
    # Compute function points in lexicographic order of domain
    cvar_map = create_loop_variables(group, func.variables)

    # Generate loops. lbody is the body of the innermost loop.
    lbody = \
        create_perfect_nested_loop(pipe, group, pipe_body,
                                   func.variables, func.domain,
                                   cparam_map, cvar_map, cfunc_map)

    arglist = func.variables
    # Convert function definition into a C expression and add it to
    # loop body
    for case in func.defn:
        if(isinstance(case, AbstractExpression)):
            case_expr = generate_c_expr(pipe, case,
                                        cparam_map, cvar_map, cfunc_map)
            array_ref = generate_c_expr(pipe, func(*arglist),
                                        cparam_map, cvar_map, cfunc_map)
            assign = genc.CAssign(array_ref, case_expr)
            lbody.add(assign, False)
        elif(isinstance(case, Case)):
            c_cond = generate_c_cond(pipe, case.condition,
                                     cparam_map, cvar_map, cfunc_map)
            case_expr = generate_c_expr(pipe, case.expression,
                                        cparam_map, cvar_map, cfunc_map)
            cif = genc.CIfThen(c_cond)

            if (isinstance(case.expression, AbstractExpression)):
                array_ref = generate_c_expr(pipe, func(*arglist),
                                            cparam_map, cvar_map, cfunc_map)
                assign = genc.CAssign(array_ref, case_expr)
                # FIXME: aliased referencing works, but direct call to
                # add method fails with assertion on block._is_open()
                with cif.if_block as ifblock:
                    ifblock.add(assign)
            else:
                assert("Invalid expression instance:"+str(expr) and \
                       False)
            lbody.add(cif, False)
        else:
            assert("Invalid case instance:"+str(case) and \
                   False)

# TESTME
def generate_reduction_scan_loops(pipe, group, comp, pipe_body,
                                  cparam_map, cfunc_map):
    """
    generates code for Reduction class
    """
    func = comp.func
    # Compute Reduction points in lexicographic order of reduction domain
    cvar_map = create_loop_variables(group, func.reductionVariables)

    # Generate loops. lbody is the body of the innermost loop.
    lbody = \
        create_perfect_nested_loop(pipe, group, pipe_body,
                                   func.reductionVariables,
                                   func.reductionDomain,
                                   cparam_map, cvar_map, cfunc_map)

    # Convert function definition into a C expression and add it to loop body
    for case in func.defn:
        if(isinstance(case, Reduce)):
            case_expr = generate_c_expr(pipe, case.expression,
                                        cparam_map, cvar_map, cfunc_map)
            ref_args = case.accumulate_ref.arguments
            accum_ref = generate_c_expr(pipe, obj(*refArgs),
                                        cparam_map, cvar_map, cfunc_map)
            assign = genc.CAssign(accum_ref, accum_ref + case_expr)
            lbody.add(assign, False)
        elif(isinstance(case, Case)):
            c_cond = generate_c_cond(pipe, case.condition,
                                     cparam_map, cvar_map, cfunc_map)
            cond_expr = generate_c_expr(pipe, case.expression,
                                        cparam_map, cvar_map, cfunc_map)
            cif = genc.CIfThen(c_cond)

            if(isinstance(case.expression, Reduce)):
                ref_args = case.accumulate_ref.arguments
                accum_ref = generate_c_expr(pipe, func(*refArgs),
                                            cparam_map, cvar_map, cfunc_map)
                assign = genc.CAssign(accum_ref, accum_ref + cond_expr)
                with cif.if_block as ifblock:
                    ifblock.add(assign)
            else:
                assert("Invalid expression instance:"+str(case.expression) and \
                       False)

            lbody.add(cif, False)
        else:
            assert("Invalid case instance:"+str(case) and \
                   False)

def generate_code_for_group(pipeline, g, body, options,
                            cparam_map, cfunc_map,
                            outputs, outputs_no_alloc):

    g.polyRep.generate_code()

    # NOTE uses the level_no of the first polypart of each compute object of
    # the group as the key for sorting compare operator. *Idea is that all
    # parts of a compute object bears the same level_no*, thus repeated calling
    # of 'order_compute_objs' can be avoided.
    group_part_map = g.polyRep.poly_parts
    sorted_comps = g.get_sorted_comps()

    # ***
    log_str = str([comp.func.name for comp in sorted_comps])
    LOG(logging.DEBUG, log_str)
    # ***

    # list of arrays to be freed
    group_freelist = []

    pooled = 'pool_alloc' in pipeline._options
    flatten_scratchpad = 'flatten_scratchpad' in pipeline._options

    for comp in sorted_comps:
        func = comp.func

        # ***
        LOG(logging.DEBUG, func.name)
        # ***

        # nothing to do!
        if comp.is_image_typ:
            continue

        # 1. Allocate storage
        array_dims = len(func.variables)

        # 1.1. scratchpad allocation, wherever applicable
        reduced_dims = [ -1 for i in range(0, len(func.domain))]
        scratch = [ False for i in range(0, len(func.domain))]
        if comp in group_part_map and not comp.is_liveout:
            for part in group_part_map[comp]:
                for i in range(0, len(func.domain)):
                    if i in part.dim_scratch_size:  # as a key
                        reduced_dims[i] = max(reduced_dims[i],
                                              part.dim_scratch_size[i])
                        scratch[i] = True

        # 1.2. prepare the sizes of each dimension for allocation
        dims = []
        for i in range(0, len(func.domain)):
            interval = func.domain[i]
            if reduced_dims[i] == -1:
                # NOTE interval step is always +1
                dim_expr = simplify_expr(interval.upperBound -
                                         interval.lowerBound + 1)

                # FIXME Creating both a cVariable (for C declaration) and a
                # Variable (to append to dims) with same namestring.
                # Link both these properly.
                dim_var_name = new_temp()

                dim_var = Variable(Int, dim_var_name)
                dim_c_var = genc.CVariable(genc.c_int, dim_var_name)

                dim_var_decl = \
                            genc.CDeclaration(genc.c_int, dim_c_var, dim_expr)

                body.add(dim_var_decl)
                dims.append(dim_var)
            else:
                dims.append(reduced_dims[i])

        # 1.3. declare and allocate arrays
        array_type = genc.TypeMap.convert(func.typ)
        array = genc.CArray(array_type, func.name, dims)

        # full array
        if comp.is_liveout:
            array.layout = 'contiguous'
            # do not allocate output arrays if they are already allocated
            if not comp.is_output or not outputs_no_alloc:
                array_ptr = genc.CPointer(array_type, 1)
                array_decl = genc.CDeclaration(array_ptr, array)
                body.add(array_decl)
                array.allocate_contiguous(body, pooled)
        else:
            if flatten_scratchpad:
                array.layout = 'contiguous_static'

        cfunc_map[func] = (array, scratch)

        # array is freed, if comp is a group liveout and not an output
        if not comp.is_output and comp.is_liveout:
            group_freelist.append(array)

        if comp in g.polyRep.poly_parts:
            continue

        # 1.4. generate scan loops
        if type(func) == Function:
            generate_function_scan_loops(pipeline, g, comp, body,
                                         cparam_map, cfunc_map)
        elif type(func) == Reduction:
            generate_reduction_scan_loops(pipeline, g, comp, body,
                                          cparam_map, cfunc_map)
        else:  # less likely
            assert("Invalid compute object type:"+str(type(func))+\
                   " of object:"+str(func.name) and \
                   False)

    # 2. generate code for built isl ast
    polyrep = g.polyRep
    if polyrep.polyast != []:
        for ast in polyrep.polyast:
            generate_c_naive_from_isl_ast(pipe, polyrep, ast, body,
                                          cparam_map, cfunc_map, pooled)
            pass

    return group_freelist

def generate_code_for_pipeline(pipeline,
                               outputs_no_alloc=False,
                               is_extern_c_func=False,
                               are_io_void_ptrs=False):

    sorted_groups = sort_scheduled_objs(pipeline.group_schedule)

    # Create a top level module for the pipeline
    m = genc.CModule('Pipeline')

    # 1. Add header files which are requried by the pipeline
    with m.includes as inc_block:
        inc_block.add(genc.CInclude('stdio.h'))
        inc_block.add(genc.CInclude('stdlib.h'))
        inc_block.add(genc.CInclude('malloc.h'))
        inc_block.add(genc.CInclude('cmath'))
        inc_block.add(genc.CInclude('string.h'))
        inc_block.add(genc.CMacroDecl(genc.c_macro_min))
        inc_block.add(genc.CMacroDecl(genc.c_macro_max))
        inc_block.add(genc.CMacroDecl(genc.c_macro_floord))
        # TODO add the pool allocation header and option

    # 2. Add function blocks
    with m.funcs as func_block:

        # Maps from pipeline parameters and functions to c variables and
        # arrays. These maps are sent to each group for code generation.
        # They are updated during the code generation of each group to
        # include liveout functions of the group.
        cparam_map = {}
        cfunc_map = {}

        # Dictonary with all the pipeline arguments
        pipeline_args = OrderedDict()

        # 2.1. Collect all the inputs and parameters of the pipeline and
        # add them as pipeline function arguments.
        params = []
        for g in sorted_groups:
            params = params + g.getParameters()
        # Remove duplicates and sort by name
        params = list(set(params))
        params.sort(key=lambda x: x.name)

        # 2.1. collect pipeline parameters
        for param in params:
            cvar_type = genc.TypeMap.convert(param.typ)
            cvar = genc.CVariable(cvar_type, param.name)
            # Bind parameters to C variables
            cparam_map[param] = cvar
            pipeline_args[cvar] = cvar.typ

        # 2.2. collect inputs
        inputs = sorted(pipeline.inputs, key=lambda x: x.name)
        for img in inputs:
            if are_io_void_ptrs:
                img_type = genc.c_void
                cptr = genc.CPointer(img_type, 1)
                cvar = genc.CVariable(cptr, img.name+'_void_arg')
            else:
                img_type = genc.TypeMap.convert(img.typ)
                cptr = genc.CPointer(img_type, 1)
                cvar = genc.CVariable(cptr, img.name)
            pipeline_args[cvar] = cvar.typ
            # Bind input functions to C arrays
            carr = genc.CArray(img_type, img.name, img.dimensions,
                               'contiguous')
            cfunc_map[img] = (carr, [])

        # 2.3. collect outputs
        outputs = sorted(pipeline.outputs, key=lambda x: x.name)
        if outputs_no_alloc:
            pass_by_type = genc.CReference
        else:
            pass_by_type = genc.CPointer

        for out in outputs:
            if are_io_void_ptrs:
                out_typ = genc.c_void
                cptr = pass_by_type(out_typ, 1)
                cvar = genc.CVariable(cptr, out.name+'_void_arg')
            else:
                out_typ = genc.TypeMap.convert(out.typ)
                cptr = pass_by_type(out_typ, 1)
                cvar = genc.CVariable(cptr, out.name)

            pipeline_args[cvar] = cvar.typ

        # 2.4. function name and declaration
        cpipe_name = 'pipeline_' + pipeline.name
        cpipe = genc.CFunction(genc.c_void, cpipe_name, pipeline_args)
        cpipe_decl = genc.CFunctionDecl(cpipe,
                                        is_extern_c_func,
                                        are_io_void_ptrs)
        cpipe_body = genc.CFunctionBody(cpipe_decl)

        func_block.add(cpipe_body)

        # 2.5. function body
        with cpipe_body.body as pbody:

            # If the code being generated is going to be compiled as shared
            # library, and used by python (through ctypes), the i/o data array
            # pointers will be given as void pointers. These should be casted
            # to their respective types first.

            # 2.5.1. typecast the i/o array ptrs
            if are_io_void_ptrs:
                inouts = list(set(inputs) | set(outputs))
                for inout in inouts:
                    # actual input to be used
                    var_type = genc.TypeMap.convert(inout.typ)
                    var_ptr = genc.CPointer(var_type, 1)
                    var = genc.CVariable(var_type, inout.name)
                    var_decl = genc.CDeclaration(var_ptr, var)
                    pbody.add(var_decl)

                    # dummy void * input/output taken as argument
                    dummy_type = genc.TypeMap.convert(inout.typ)
                    dummy_ptr = genc.CPointer(dummy_type, 1)
                    dummy_cast = genc.CCast(dummy_ptr, inout.name+'_void_arg')
                    var_assign = genc.CAssign(var, dummy_cast)
                    pbody.add(var_assign)

            # Boolean to check if the memory allocations should be done using
            # malloc or the custom pool allocator
            pooled = 'pool_alloc' in pipeline._options

            pipe_freelist = []
            for g in sorted_groups:
                group_freelist = \
                    generate_code_for_group(pipeline, g, pbody,
                                            pipeline._options,
                                            cparam_map, cfunc_map,
                                            outputs, outputs_no_alloc)
                pipe_freelist.extend(group_freelist)

            # TODO free the arrays ASAP (compaction)
            # 3. Deallocate storage
            for array in pipe_freelist:
                array.deallocate(pbody, pooled)

    return m

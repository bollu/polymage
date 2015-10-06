from __future__ import absolute_import, division, print_function

from constructs import *

def get_affine_var_and_param_coeff(expr):
    expr = Value.numericToValue(expr)

    if (not isAffine(expr) or \
            isinstance(expr, Value) or \
            is_constant_expr(expr)):
        return {}
    elif (isinstance(expr, Variable)):
        return {expr: 1}
    elif (isinstance(expr, AbstractBinaryOpNode)):
        coeff = {}
        left_coeff = get_affine_var_and_param_coeff(expr.left)
        right_coeff = get_affine_var_and_param_coeff(expr.right)
        if (expr.op == '+'):
            coeff = dict((n, left_coeff.get(n, 0) + right_coeff.get(n, 0)) \
                             for n in set(left_coeff) | set(right_coeff))
        elif (expr.op == '-'):
            coeff = dict((n, left_coeff.get(n, 0) - right_coeff.get(n, 0)) \
                             for n in set(left_coeff) | set(right_coeff))
        elif (expr.op == '*'):
            left_is_constant = is_constant_expr(expr.left, affine=True)
            right_is_constant = is_constant_expr(expr.right, affine=True)
            #sanity check should be true if the expression is affine
            assert(not (left_is_constant and right_is_constant))

            if (left_is_constant and not right_is_constant):
                coeff = dict((n, get_constant_from_expr(expr.left, \
                                                        affine=True) \
                                 * right_coeff.get(n, 0)) \
                                 for n in set(right_coeff))
            elif(right_is_constant and not left_is_constant):
                coeff = dict((n, get_constant_from_expr(expr.right, \
                                                        affine=True) \
                                 * left_coeff.get(n, 0))
                                 for n in set(left_coeff))
        elif (expr.op == '/'):
            right_is_constant = is_constant_expr(expr.right, affine=True)
            #sanity check should be true if the expression is affine
            assert(right_is_constant)

            coeff = dict((n, Fraction(1, get_constant_from_expr(expr.right, \
                                                             affine=True)) \
                                         * left_coeff.get(n, 0)) \
                                         for n in set(left_coeff))
        return coeff
    elif (isinstance(expr, AbstractUnaryOpNode)):
        child_coeff = get_affine_var_and_param_coeff(expr.child)
        if (expr.op == '-'):
            child_coeff = dict((n, -child_coeff.get(n)) for n in child_coeff)
        return child_coeff
    raise TypeError(type(expr))

def evaluateBinaryOp(leftVal, rightVal, op, affine):
    if op == '/':
        if affine:
            return Fraction(leftVal, rightVal)
        else:
            return leftVal/rightVal
    elif op == '+':
        return leftVal + rightVal
    elif op == '*':
        return leftVal * rightVal
    elif op == '-':
        return leftVal - rightVal
    elif op == '%':
        return leftVal % rightVal
    elif op == '>>':
        return leftVal >> rightVal
    elif op == '<<':
        return leftVal << rightVal
    elif op == '|':
        return leftVal | rightVal
    elif op == '^':
        return leftVal ^ rightVal
    elif op == '&':
        return leftVal & rightVal
    raise TypeError(type(op))

def evaluateUnaryOp(val, op):
    if op == '-':
        return -val
    elif op == '+':
        return val
    raise TypeError(type(val))

def get_constant_from_expr(expr, affine=False):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        if affine:
            assert (expr.typ is Int) or (expr.typ is Rational)
        return expr.value 
    elif (isinstance(expr, Variable)):
        return 0
    elif (isinstance(expr, Reference)):
        return 0    
    elif (isinstance(expr, AbstractBinaryOpNode)):
        leftConst = get_constant_from_expr(expr.left, affine)
        rightConst = get_constant_from_expr(expr.right, affine)
        return evaluateBinaryOp(leftConst, rightConst, expr.op, affine)
    elif (isinstance(expr, AbstractUnaryOpNode)):
        childConst = get_constant_from_expr(expr.child, affine)
        return evaluateUnaryOp(childConst, expr.op)
    elif (isinstance(expr, (Select, Cast, InbuiltFunction))):
        return 0
    raise TypeError(type(expr))

def is_constant_expr(expr, affine = False):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        if affine:
            return (expr.typ is Int) or (expr.typ is Rational)
        return True 
    elif (isinstance(expr, Variable)):
        return False
    elif (isinstance(expr, Reference)):
        return False    
    elif (isinstance(expr, AbstractBinaryOpNode)):
        return (is_constant_expr(expr.left, affine) and 
                is_constant_expr(expr.right, affine))
    elif (isinstance(expr, AbstractUnaryOpNode)):
        return is_constant_expr(expr.child, affine)
    elif (isinstance(expr, (Select, InbuiltFunction, Cast))):
        return False
    raise TypeError(type(expr))

def simplify_expr(expr):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.clone()
    elif (isinstance(expr, Variable)):
        return expr.clone()
    elif (isinstance(expr, Reference)):
        simple_args = []
        for arg in expr.arguments:
            simple_args.append(simplify_expr(arg))
        # Equivalent to cloning
        return expr.objectRef(*simple_args)
    elif (isinstance(expr, AbstractBinaryOpNode) or
          isinstance(expr, AbstractUnaryOpNode)):
        if (isAffine(expr, includeDiv=False)):
            coeff = get_affine_var_and_param_coeff(expr)
            variables = list(set(expr.collect(Variable)))
            params = list(set(expr.collect(Parameter)))
            simple_expr = get_constant_from_expr(expr)
            for var in variables:
                if (coeff[var] == 1):
                    simple_expr = simple_expr + var
                elif (coeff[var] == -1):
                    simple_expr = simple_expr - var
                else:
                    simple_expr = simple_expr + coeff[var] * var
            for param in params:
                if (coeff[param] == 1):
                    simple_expr = simple_expr + param
                elif (coeff[param] == -1):
                    simple_expr = simple_expr - param
                else:
                    simple_expr = simple_expr + coeff[param] * param
            return Value.numericToValue(simple_expr)
        else:
            return expr.clone()
    elif (isinstance(expr, (Select, Cast, InbuiltFunction))):
        # Some simplification can be done but ignoring for now
        return expr.clone()
    raise TypeError(type(expr))

def substitute_refs(expr, ref_to_expr_map):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.clone()
    elif (isinstance(expr, Variable)):
        return expr.clone()
    elif (isinstance(expr, Reference)):
        if (expr in ref_to_expr_map):
            ref_vars = expr.objectRef.variables
            var_to_expr_map = {}
            assert len(expr.arguments) == len(ref_vars)
            for i in range(0, len(ref_vars)):
                var_to_expr_map[ref_vars[i]] = expr.arguments[i]
            # Equivalent to cloning
            return substitute_vars(ref_to_expr_map[expr], var_to_expr_map)
        else:
            return expr.clone()
    elif (isinstance(expr, AbstractBinaryOpNode)):
        left = substitute_refs(expr.left, ref_to_expr_map)
        right = substitute_refs(expr.right, ref_to_expr_map)
        op = expr.op
        new_expr = AbstractBinaryOpNode(left, right, op)
        return simplify_expr(new_expr)
    elif (isinstance(expr, AbstractUnaryOpNode)):
        child = substitute_refs(expr.child, ref_to_expr_map)
        op = expr.op
        new_expr = AbstractUnaryOpNode(child, op)
        return simplify_expr(new_expr)
    elif (isinstance(expr, InbuiltFunction)):
        expr = expr.clone()
        expr.inline_refs(ref_to_expr_map)
        return expr
    elif (isinstance(expr, (Select, Cast))):
        expr = expr.clone()
        expr.inline_refs(ref_to_expr_map)
        return expr
    raise TypeError(type(expr))

def substitute_vars(expr, var_to_expr_map):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.clone()
    elif (isinstance(expr, Variable)):
        if expr in var_to_expr_map:
            return var_to_expr_map[expr].clone()
        return expr.clone()
    elif (isinstance(expr, Reference)):
        num_args = len(expr.arguments)
        args = []
        for i in range(0, num_args):
            args.append(substitute_vars(expr.arguments[i], var_to_expr_map))
        # Equivalent to cloning
        return expr.objectRef(*args)
    elif (isinstance(expr, AbstractBinaryOpNode)):
        left = substitute_vars(expr.left, var_to_expr_map)
        right = substitute_vars(expr.right, var_to_expr_map)
        op = expr.op
        new_expr = AbstractBinaryOpNode(left, right, op)
        return simplify_expr(new_expr)
    elif (isinstance(expr, AbstractUnaryOpNode)):
        child = substitute_vars(expr.child, var_to_expr_map)
        op = expr.op
        new_expr = AbstractUnaryOpNode(child, op)
        return simplify_expr(new_expr)
    elif (isinstance(expr, Cast)):
        typ = expr.typ
        return Cast(typ, substitute_vars(expr.expression, var_to_expr_map))
    elif (isinstance(expr, Select)):
        new_cond = substitute_vars(expr.condition, var_to_expr_map)
        new_true = substitute_vars(expr.trueExpression, var_to_expr_map)
        new_false = substitute_vars(expr.falseExpression, var_to_expr_map)
        return Select(new_cond, new_true, new_false)
    elif (isinstance(expr, InbuiltFunction)):
        expr = expr.clone()
        expr.substitute_vars(var_to_expr_map)
        return expr
    raise TypeError(type(expr))


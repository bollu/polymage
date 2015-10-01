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

def simplifyExpr(expr):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.clone()
    elif (isinstance(expr, Variable)):
        return expr.clone()
    elif (isinstance(expr, Reference)):
        simpleArgs = []
        for arg in expr.arguments:
            simpleArgs.append(simplifyExpr(arg))
        # Equivalent to cloning
        return expr.objectRef(*simpleArgs)
    elif (isinstance(expr, AbstractBinaryOpNode) or
          isinstance(expr, AbstractUnaryOpNode)):
        if (isAffine(expr, includeDiv=False)):
            coeff = get_affine_var_and_param_coeff(expr)
            variables = list(set(expr.collect(Variable)))
            params = list(set(expr.collect(Parameter)))
            simpleExpr = get_constant_from_expr(expr)
            for var in variables:
                if (coeff[var] == 1):
                    simpleExpr = simpleExpr + var
                elif (coeff[var] == -1):
                    simpleExpr = simpleExpr - var
                else:
                    simpleExpr = simpleExpr + coeff[var] * var
            for param in params:
                if (coeff[param] == 1):
                    simpleExpr = simpleExpr + param
                elif (coeff[param] == -1):
                    simpleExpr = simpleExpr - param
                else:
                    simpleExpr = simpleExpr + coeff[param] * param
            return Value.numericToValue(simpleExpr)
        else:
            return expr.clone()
    elif (isinstance(expr, (Select, Cast, InbuiltFunction))):
        # Some simplification can be done but ignoring for now
        return expr.clone()
    raise TypeError(type(expr))

def substituteRefs(expr, refToExprMap):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.clone()
    elif (isinstance(expr, Variable)):
        return expr.clone()
    elif (isinstance(expr, Reference)):
        if (expr in refToExprMap):
            refVars = expr.objectRef.variables
            varToExprMap = {}
            assert len(expr.arguments) == len(refVars)
            for i in range(0, len(refVars)):
                varToExprMap[refVars[i]] = expr.arguments[i]
            # Equivalent to cloning
            return substituteVars(refToExprMap[expr], varToExprMap)
        else:
            return expr.clone()
    elif (isinstance(expr, AbstractBinaryOpNode)):
        left = substituteRefs(expr.left, refToExprMap)
        right = substituteRefs(expr.right, refToExprMap)
        op = expr.op
        newExpr = AbstractBinaryOpNode(left, right, op)
        return simplifyExpr(newExpr)
    elif (isinstance(expr, AbstractUnaryOpNode)):
        child = substituteRefs(expr.child, refToExprMap)
        op = expr.op
        newExpr = AbstractUnaryOpNode(child, op)
        return simplifyExpr(newExpr)
    elif (isinstance(expr, InbuiltFunction)):
        expr = expr.clone()
        expr.inlineRefs(refToExprMap)
        return expr
    elif (isinstance(expr, (Select, Cast))):
        expr = expr.clone()
        expr.inlineRefs(refToExprMap)
        return expr
    raise TypeError(type(expr))

def substituteVars(expr, varToExprMap):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.clone()
    elif (isinstance(expr, Variable)):
        if expr in varToExprMap:
            return varToExprMap[expr].clone()
        return expr.clone()
    elif (isinstance(expr, Reference)):
        numArgs = len(expr.arguments)
        args = []
        for i in range(0, numArgs):
            args.append(substituteVars(expr.arguments[i], varToExprMap))
        # Equivalent to cloning
        return expr.objectRef(*args)
    elif (isinstance(expr, AbstractBinaryOpNode)):
        left = substituteVars(expr.left, varToExprMap)
        right = substituteVars(expr.right, varToExprMap)
        op = expr.op
        newExpr = AbstractBinaryOpNode(left, right, op)
        return simplifyExpr(newExpr)
    elif (isinstance(expr, AbstractUnaryOpNode)):
        child = substituteVars(expr.child, varToExprMap)
        op = expr.op
        newExpr = AbstractUnaryOpNode(child, op)
        return simplifyExpr(newExpr)
    elif (isinstance(expr, Cast)):
        typ = expr.typ
        return Cast(typ, substituteVars(expr.expression, varToExprMap))
    elif (isinstance(expr, Select)):
        newCond = substituteVars(expr.condition, varToExprMap)
        newTrue = substituteVars(expr.trueExpression, varToExprMap)
        newFalse = substituteVars(expr.falseExpression, varToExprMap)
        return Select(newCond, newTrue, newFalse)
    elif (isinstance(expr, InbuiltFunction)):
        expr = expr.clone()
        expr.substituteVars(varToExprMap)
        return expr
    raise TypeError(type(expr))


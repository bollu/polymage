# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

from Constructs import *

def getAffineVarAndParamCoeff(expr):
    expr = Value.numericToValue(expr)
    if (not isAffine(expr) or isinstance(expr, Value)):
        return {}
    elif (isinstance(expr, Variable)):
        return {expr: 1}
    elif (isinstance(expr, AbstractBinaryOpNode)):
        coeff = {}
        leftCoeff = getAffineVarAndParamCoeff(expr.left)
        rightCoeff = getAffineVarAndParamCoeff(expr.right)
        if (expr.op == '+'):
            coeff = dict( (n, leftCoeff.get(n, 0) + rightCoeff.get(n, 0))\
                          for n in set(leftCoeff) | set(rightCoeff) )
        elif (expr.op == '-'):
            coeff = dict( (n, leftCoeff.get(n, 0) - rightCoeff.get(n, 0))\
                          for n in set(leftCoeff) | set(rightCoeff) )
        elif (expr.op == '*'):
            leftIsConstant = isConstantExpr(expr.left, affine = True)
            rightIsConstant = isConstantExpr(expr.right, affine = True)
            #sanity check should be true if the expression is affine
            assert(not (leftIsConstant and rightIsConstant))
            if (leftIsConstant and not rightIsConstant):
                coeff = dict( (n, getConstantFromExpr(expr.left, affine = True) *
                                  rightCoeff.get(n, 0))
                              for n in set(rightCoeff) )
            elif(rightIsConstant and not leftIsConstant):
                coeff = dict( (n, getConstantFromExpr(expr.right, affine = True) *
                                  leftCoeff.get(n, 0))
                              for n in set(leftCoeff) )
        elif (expr.op == '/'):
            rightIsConstant = isConstantExpr(expr.right, affine = True)
            #sanity check should be true if the expression is affine
            assert(rightIsConstant)
            coeff = dict( (n, Fraction(1, getConstantFromExpr(expr.right, affine = True)) *
                                       leftCoeff.get(n, 0))
                           for n in set(leftCoeff) )
        return coeff
    elif (isinstance(expr, AbstractUnaryOpNode)):
        childCoeff = getAffineVarAndParamCoeff(expr.child)
        if (expr.op == '-'):
            childCoeff = dict( (n, -childCoeff.get(n)) for n in childCoeff )
        return childCoeff
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

def getConstantFromExpr(expr, affine=False):
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
        leftConst = getConstantFromExpr(expr.left, affine)
        rightConst = getConstantFromExpr(expr.right, affine)
        return evaluateBinaryOp(leftConst, rightConst, expr.op, affine)
    elif (isinstance(expr, AbstractUnaryOpNode)):
        childConst = getConstantFromExpr(expr.child, affine)
        return evaluateUnaryOp(childConst, expr.op)
    elif (isinstance(expr, (Select, Cast, InbuiltFunction))):
        return 0
    raise TypeError(type(expr))

def isConstantExpr(expr, affine = False):
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
        return (isConstantExpr(expr.left, affine) and 
                isConstantExpr(expr.right, affine))
    elif (isinstance(expr, AbstractUnaryOpNode)):
        return isConstantExpr(expr.child, affine)
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
            coeff = getAffineVarAndParamCoeff(expr)
            variables = list(set(expr.collect(Variable)))
            params = list(set(expr.collect(Parameter)))
            simpleExpr = getConstantFromExpr(expr)
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

def getType(expr):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression))
    if (isinstance(expr, Value)):
        return expr.typ
    elif (isinstance(expr, Variable)):
        return expr.typ
    elif (isinstance(expr, Reference)):
        return expr.objectRef.typ
    elif (isinstance(expr, AbstractBinaryOpNode)):
        leftType = getType(expr.left)
        rightType = getType(expr.right)
        return typeImplictConversion(leftType, rightType) 
    elif (isinstance(expr, AbstractUnaryOpNode)):
        return getType(expr.child)
    elif (isinstance(expr, Cast)):
        return expr.typ
    elif (isinstance(expr, Select)):
        trueType = getType(expr.trueExpression)
        falseType = getType(expr.falseExpression)
        assert trueType == falseType
        return trueType
    elif (isinstance(expr, InbuiltFunction)):
        return expr.getType()
    raise TypeError(type(expr))

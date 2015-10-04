from __future__ import absolute_import, division, print_function

# TODO remove this at some point
from expr_ast import *
from expr_types import *
import logging

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s")

"""
class Max(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        if typeCheck:
            _leftTyp = getType(_leftExpr)
            _rightTyp = getType(_rightExpr)
            assert _leftTyp == _rightTyp
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
   
    def getType(self):
        return getType(self._args[0])

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Max(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "isl_max(" + leftStr + ", " + rightStr + ")"  

class Min(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        if typeCheck:
            _leftTyp = getType(_leftExpr)
            _rightTyp = getType(_rightExpr)
            assert _leftTyp == _rightTyp
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
    
    def getType(self):
        return getType(self._args[0])

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Min(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "isl_min(" + leftStr + ", " + rightStr + ")"
"""

class Pow(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr):
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
   
    def getType(self):
        return Double

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Pow(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "pow(" + leftStr + ", " + rightStr + ")"

class Powf(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr):
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
    
    def getType(self):
        return Float

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Powf(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "powf(" + leftStr + ", " + rightStr + ")"

class Exp(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return Double

    def clone(self):
        return Exp(self._args[0].clone())

    def __str__(self):
        return "std::exp(" +  self._args[0].__str__() +  ")"

class Sin(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return Double

    def clone(self):
        return Sin(self._args[0].clone())

    def __str__(self):
        return "std::sin(" +  self._args[0].__str__() +  ")"

class Cos(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Double

    def clone(self):
        return Cos(self._args[0].clone())

    def __str__(self):
        return "std::cos(" +  self._args[0].__str__() +  ")"

class Sqrt(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Double

    def clone(self):
        return Sqrt(self._args[0].clone())

    def __str__(self):
        return "std::sqrt(" +  self._args[0].__str__() +  ")"

class Sqrtf(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Float

    def clone(self):
        return Sqrtf(self._args[0].clone())

    def __str__(self):
        return "std::sqrtf(" +  self._args[0].__str__() +  ")"

class Abs(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return getType(self._args[0])

    def clone(self):
        return Abs(self._args[0].clone())

    def __str__(self):
        return "std::abs(" +  self._args[0].__str__() +  ")"

class Cast(AbstractExpression):
    def __init__(self, _typ, _expr):
        _expr = Value.numericToValue(_expr)
        assert _typ in [Float, Double, 
                        UChar, Char, 
                        UShort, Short, 
                        UInt, Int, 
                        ULong, Long]
        assert(isinstance(_expr, AbstractExpression))
        self._typ  = _typ
        self._expr = _expr

    @property
    def typ(self):
        return self._typ

    @property 
    def expression(self):
        return self._expr

    def collect(self, objType):
        objs = []
        objs += self._expr.collect(objType)
        if (type(self) == objType):
            objs += [self]
        return list(set(objs))

    def clone(self):
        return Cast(self._typ, self._expr.clone())
    
    def replaceReferences(self, refToExprMap):
        self._expr = substituteRefs(self._expr, refToExprMap)

    def __str__(self):
        exprStr = self._expr.__str__()
        return "(" + self._typ.cTypeName() + ") " + "(" + exprStr + ")"

class Select(AbstractExpression):
    def __init__(self, _cond, _trueExpr, _falseExpr, typeCheck = True):
        assert(isinstance(_cond, Condition))
        _trueExpr = Value.numericToValue(_trueExpr)
        _falseExpr = Value.numericToValue(_falseExpr)
        assert(isinstance(_trueExpr, AbstractExpression))
        assert(isinstance(_falseExpr, AbstractExpression))
        if typeCheck:
            trueType = getType(_trueExpr)
            falseType = getType(_falseExpr)
            assert trueType == falseType
        self._trueExpr = _trueExpr
        self._falseExpr = _falseExpr
        self._cond = _cond

    @property 
    def condition(self):
        return self._cond
    @property
    def trueExpression(self):
        return self._trueExpr     
    @property
    def falseExpression(self):
        return self._falseExpr

    def collect(self, objType):
        objs = []
        objs += self._cond.collect(objType)
        objs += self._trueExpr.collect(objType)
        objs += self._falseExpr.collect(objType)
        if (type(self) == objType):
            objs += [self]
        return list(set(objs))

    def replaceReferences(self, refToExprMap):
        self._cond.replaceReferences(refToExprMap)
        self._trueExpr = substituteRefs(self._trueExpr, refToExprMap)
        self._falseExpr = substituteRefs(self._falseExpr, refToExprMap)
    
    def clone(self):
        return Select(self._cond.clone(),
                      self._trueExpr.clone(),
                      self._falseExpr.clone())

    def __str__(self):
        condStr = self._cond.__str__()
        trueStr = self._trueExpr.__str__()
        falseStr = self._falseExpr.__str__()
        return "(" + condStr + "? " + trueStr + ": " + falseStr + ")"

class Max(Select):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        Select.__init__(self, Condition(_leftExpr, '>', _rightExpr),
                                        _leftExpr, _rightExpr, typeCheck)
class Min(Select):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        Select.__init__(self, Condition(_leftExpr, '<', _rightExpr), 
                                        _leftExpr, _rightExpr, typeCheck)

class Variable(AbstractExpression):
    def __init__(self, _typ, _name):
        self._name    = _name
        self._typ     = _typ
    
    @property 
    def name(self):
        return self._name

    @property 
    def typ(self):
        return self._typ

    def clone(self):
        return self

    def __str__(self):
        return self._name.__str__()

class Parameter(Variable):
    def __init__(self, _typ, _name):
        Variable.__init__(self, _typ, _name)
    
class Interval(object):
    def __init__(self, _typ,  _lb, _ub):
        _lb   = Value.numericToValue(_lb)
        _ub   = Value.numericToValue(_ub)
        assert(isinstance(_lb, AbstractExpression))
        assert(isinstance(_ub, AbstractExpression))
        self._lb   = _lb
        self._ub   = _ub
        self._typ  = _typ
 
    @property 
    def lowerBound(self):
        return self._lb
    @property 
    def upperBound(self):
        return self._ub
    @property 
    def typ(self):
        return self._typ

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._lb.collect(objType)
        objs += self._ub.collect(objType)
        return list(set(objs))

    def clone(self):
        return Interval(self._typ, 
                        self._lb.clone(), 
                        self._ub.clone())

    def __str__(self):
        return '(' + self._lb.__str__() + ', ' +\
               self._ub.__str__() + ')'

class Reference(AbstractExpression):
    def __init__(self, _obj, _args):
        _args = [ Value.numericToValue(arg) for arg in _args]
        for arg in _args:
            assert(isinstance(arg, AbstractExpression))
        self._obj  = _obj
        self._args = _args

    @property
    def objectRef(self):
        return self._obj

    def _replaceRefObject(self, cloneObj):
        self._obj = cloneObj

    @property
    def arguments(self):
        return self._args

    def clone(self):
        cloneArgs = [arg.clone() for arg in self._args]
        return Reference(self._obj, cloneArgs)

    def collect(self, objType):
        objs = []
        for arg in self._args:
            objs += arg.collect(objType)
        if (type(self) is objType):
            objs += [self]
        return list(set(objs))

    def __str__(self):
        arg_str = ", ".join([arg.__str__() for arg in self._args])
        return self._obj.name + "(" + arg_str + ")"

class Condition(object):
    def __init__(self, _left, _cond, _right):
        _left  = Value.numericToValue(_left)
        _right = Value.numericToValue(_right)
        assert(_cond in ['<', '<=', '>', '>=', '==', '!=', '&&', '||'])
        if _cond in ['<', '<=', '>', '>=', '==', '!=']:
            assert(isinstance(_left, AbstractExpression))
            assert(isinstance(_right, AbstractExpression))
        if _cond in ['&&', '||']:
            assert(isinstance(_left, Condition))
            assert(isinstance(_right, Condition))
        self._left  = _left
        self._right = _right
        self._cond  = _cond

    @property
    def lhs(self):
        return self._left
    @property
    def rhs(self):
        return self._right
    @property
    def conditional(self):
        return self._cond
   
    def clone(self):
        return Condition(self._left.clone(), self._cond, 
                         self._right.clone())

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._left.collect(objType) + self._right.collect(objType)
        return list(set(objs))

    def replaceReferences(self, refToExprMap):
        if(isinstance(self._left, Condition)):
            self._left.replaceReferences(refToExprMap)
        else:
            self._left = substituteRefs(self._left, refToExprMap)
        if(isinstance(self._right, Condition)):
            self._right.replaceReferences(refToExprMap)
        else:
            self._right = substituteRefs(self._right, refToExprMap)

    def splitToConjuncts(self):
        conjuncts = []
        if self._cond in ['<', '<=', '>', '>=', '==']:
            conjuncts.append([self])
        elif (self._cond  == '!='):
            lessThan = Condition(self._left, '<', self._right)
            conjuncts.append([lessThan])
            greaterThan = Condition(self._left, '>', self._right)
            conjuncts.append([greaterThan])
        elif (self._cond == '||'):
            conjuncts = self._left.splitToConjuncts() + \
                        self._right.splitToConjuncts()
        elif (self._cond == '&&'):
            leftConjuncts = self._left.splitToConjuncts()
            rightConjuncts = self._right.splitToConjuncts()
            for lconjunct in leftConjuncts:
                for rconjunct in rightConjuncts:
                    conjuncts.append(lconjunct + rconjunct)
        else:
            assert False
        return conjuncts

    def __and__(self, other):
        assert(isinstance(other, Condition))
        return Condition(self, '&&', other)

    def __or__(self, other):
        assert(isinstance(other, Condition))
        return Condition(self, '||', other)
    
    def __str__(self):        
        if (self._cond is None):
            assert self._left is None and self._right is None
            return ""
        left_str = self._left.__str__()
        right_str = self._right.__str__()
        return "(" + left_str + " " + self._cond + " " + right_str + ")"

class Case(object):
    def __init__(self, _cond, _expr):
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_cond, Condition))
        assert(isinstance(_expr, (AbstractExpression, Reduce)))
        self._cond  = _cond
        self._expr  = _expr
    @property
    def condition(self):
        return self._cond
    @property
    def expression(self):
        return self._expr

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._cond.collect(objType) + self._expr.collect(objType)
        return list(set(objs))

    def replaceReferences(self, refToExprMap):
        self._cond.replaceReferences(refToExprMap)
        self._expr = substituteRefs(self._expr, refToExprMap)

    def clone(self):
        return Case(self._cond.clone(), self._expr.clone())

    def __str__(self):
        return 'Case(' + self._cond.__str__() + ')' +\
                '{ ' + self._expr.__str__() + ' }'

class Op(object):
    Sum = 0
    Mul = 1
    Min = 2
    Max = 3

class Reduce(object):
    def __init__(self, _red_ref, _expr, _op_typ):
        assert isinstance(_red_ref, Reference)
        assert isinstance(_red_ref.objectRef, Reduction)
        _expr = Value.numericToValue(_expr)
        assert isinstance(_expr, AbstractExpression)
        assert _op_typ in [Op.Sum, Op.Mul, Op.Min, Op.Max]

        self._red_ref = _red_ref
        self._expr = _expr
        self._op_typ = _op_typ

    @property 
    def accumulate_ref(self):
        return self._red_ref
    @property 
    def expression(self):
        return self._expr
    @property 
    def op_type(self):
        return self._expr

    def replaceReferences(self, refToExprMap):
        self._expr = substituteRefs(self._expr, refToExprMap)
        self._red_ref = substituteRefs(self._red_ref, refToExprMap)
   
    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        
        objs = self._red_ref.collect(objType) + self._expr.collect(objType)
        return list(set(objs))  

    def clone(self):
        return Reduce(self._red_ref.clone(), self._expr.clone(), self._op_typ)

    def __str__(self):
        op_str = None
        op_sep = None
        if (self._op_typ == Op.Sum):
            op_str = ''
            op_sep = ' + '
        elif (self._op_typ == Op.Mul):
            op_str = ''
            op_sep = ' * '
        elif (self._op_typ == Op.Min):
            op_str = 'Min'
            op_sep = ', '
        elif (self._op_typ == Op.Max):
            op_str = 'Max'
            op_sep = ', '
        else:
            assert False

        ret_str = 'Reduce [ ' + self._red_ref.__str__() + ' = ' + \
                  op_str + '(' + \
                  self._red_ref.__str__() + op_sep + self._expr.__str__() + ') ]'
        return ret_str

class Function(object):
    def __init__(self, _varDom, _typ, _name):
        self._name      = _name
        # Type of the scalar range of the function
        self._typ       = _typ
        # Variables of the function
        self._variables = None

        # Gives the domain of each variable. Domain of each variable is
        # expected to be over integers. Function evaluation in the
        # lexicographic order of the domain is assumed to be valid.

        assert(len(_varDom[0]) == len(_varDom[1]))
        for i in range(0, len(_varDom[0])):
            assert(isinstance(_varDom[0][i], Variable))
            assert(isinstance(_varDom[1][i], Interval))
            assert(_varDom[0][i].typ ==  _varDom[1][i].typ)
        # add check to ensure that upper bound and lower bound
        # expressions for each variable are only defined in 
        # terms of variables of the function and global parameters

        # Should bounds be restricted only to parameters or function
        # variables be allowed? No for now

        # Should the domain be restricted to the positive quadrant?
        # Can this be done automatically
        self._variables = _varDom[0]
        self._varDomain = _varDom[1]

        # * Body of a function is composed of Case and Expression constructs.
        # * The Case constructs are expected to be non-overlapping. Therefore,
        #   value at each point in the function domain is uniquely defined.
        self._body      = []

    @property
    def name(self):
        return self._name
    @property
    def typ(self):
        return self._typ

    @property
    def variableDomain(self):
        return (self._variables, self._varDomain)
           
    @property
    def domain(self):
        return self._varDomain

    @property
    def variables(self):
        return self._variables

    @property
    def defn(self):
        return self._body
    @defn.setter
    def defn(self, _def):
        assert(self._body == [])
        for case in _def:
            case = Value.numericToValue(case)
            assert(isinstance(case, (Case, AbstractExpression)))
            # check if the Case and Expression constructs only use
            # function variables and global parameters

            # MOD -> if _def is not a Case, shouldnt it be disallowed after
            # the first definition?
            self._body.append(case)

    def __call__(self, *args):
        assert(len(args) == len(self._variables))
        for arg in args:
            arg = Value.numericToValue(arg)
            assert(isinstance(arg, AbstractExpression))
        return Reference(self, args)

    def replaceReferences(self, refToExprMap):
        numCases = len(self._body)
        for i in range(0, numCases):
            if isinstance(self._body[i], Case):
                self._body[i].replaceReferences(refToExprMap)
            else:
                self._body[i] = substituteRefs(self._body[i], refToExprMap)

    def getObjects(self, objType):
        objs = []
        for case in self._body:
            objs += case.collect(objType)
        for interval in self._varDomain:
            objs += interval.collect(objType)
        return list(set(objs))

    def hasBoundedIntegerDomain(self):
        boundedIntegerDomain = True
        for varDom in self._varDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
                    break
            else:
                boundedIntegerDomain = False
                break

        return boundedIntegerDomain

    def clone(self):
        newBody = [ c.clone() for c in self._body ]
        varDom = ( [ v.clone() for v in self._variables], 
                   [ d.clone() for d in self._varDomain] )
        newFunc = Function(varDom, self._typ, self._name)
        newFunc.defn = newBody
        return newFunc
    
    def __str__(self):
        if (self._body):
            var_str = ", ".join([var.__str__() for var in self._variables])
            dom_str = ', '.join([self._variables[i].__str__() + self._varDomain[i].__str__()\
                                for i in range(len(self._varDomain))])
            case_str = "{ " + "\n ".join([case.__str__() for case in self._body]) + " }"
            return "Domain: " + dom_str + '\n' + self._name + "(" + var_str + ") = " +\
                    case_str + '\n'
        else:
            return self._name

class Image(Function):
    def __init__(self, _typ, _name, _dims):
        _dims = [ Value.numericToValue(dim) for dim in _dims ]
        # Have to evaluate if a  stronger constraint 
        # can be imposed. Only AbstractExpression in parameters?
        for dim in _dims:
            assert(isinstance(dim, AbstractExpression))
        self._dims = _dims
        intervals = []
        variables = []
        i = 0
        for dim in self._dims:
            # Just assuming it will not be more that UInt
            intervals.append(Interval(UInt, 0, dim-1))
            variables.append(Variable(UInt, "_" + _name + str(i)))
            i = i + 1
        Function.__init__(self, (variables, intervals), _typ, _name)

    @property
    def dimensions(self):
        return tuple(self._dims)

    def __str__(self):
        dim_str = ", ".join([dim.__str__() for dim in self._dims])
        return self._name.__str__() + "(" + dim_str + ")"

class Reduction(Function):
    def __init__(self, _varDom, _redDom, _typ, _name):
        Function.__init__(self, _varDom, _typ, _name)
        # Gives the domain of the reduction. Reduction domain of each variable is 
        # expected to be over integers. Reduction evaluation in the lexicographic 
        # order of the domain is assumed to be valid.
        assert(len(_redDom[0]) == len(_redDom[1]))
        for i in range(0, len(_redDom[0])):
            assert(isinstance(_redDom[0][i], Variable))
            assert(isinstance(_redDom[1][i], Interval))
            assert(_redDom[0][i].typ ==  _redDom[1][i].typ)
        # add check to ensure that upper bound and lower bound
        # expressions for each variable are only defined in
        # terms of variables of the function and global parameters

        # Should bounds be restricted only to parameters or function
        # variables be allowed? No for now

        # Should the domain be restricted to the positive quadrant?
        # Can this be done automatically
        self._redVariables = _redDom[0]
        self._redDomain = _redDom[1]

        # Intial value of each accumulator cell. Default is set to zero of the 
        # given type
        self._default   = Value(0, _typ)

    @property
    def default(self):
        return self._default
    @default.setter
    def default(self, _expr):
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_expr, AbstractExpression))
        self._default = _expr

    @property
    def reductionDomain(self):
        return self._redDomain

    @property
    def reductionVariables(self):
        return self._redVariables
        
    @property
    def defn(self):
        return self._body

    @defn.setter
    def defn(self, _def):
        assert(self._body == [])
        for case in _def:
            case = Value.numericToValue(case)
            assert(isinstance(case, (Case, Reduce)))
            # check if the Case and Expression constructs only use
            # function variables and global parameters
    
            # Which way is better Case inside accumulate or accumulate inside Case
            
            # MOD -> if _def is not a Case, shouldnt it be disallowed after
            # the first definition?
            self._body.append(case)

    def hasBoundedIntegerDomain(self):
        boundedIntegerDomain = True
        for varDom in self._varDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False

        for redDom in self._redDomain:
            if isinstance(redDom, Interval):
                if(not isAffine(redDom.lowerBound) or
                   not isAffine(redDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False

        return boundedIntegerDomain

    def getObjects(self, objType):
        objs = []
        for case in self._body:
            objs += case.collect(objType)
        for interval in self._varDomain:
            objs += interval.collect(objType)
        for interval in self._redDomain:
            objs += interval.collect(objType)
        objs += self._default.collect(objType)
        return list(set(objs))
   
    def clone(self):
        newBody = [ r.clone() for r in self._body ]
        varDom = ( [ v.clone() for v in self._variables], 
                   [ d.clone() for d in self._varDomain] )
        redDom = ( [ r.clone() for r in self._redVariables],
                   [ d.clone() for d in self._redDomain] )
        newRed = Reduction(varDom, redDom, self._typ, self._name)
        newRed.defn = newBody
        newRed.default = self._default.clone()
        return newRed    

    def __str__(self):
        if (self._body):
            varStr = ", ".join([var.__str__() for var in self._variables])
            domStr = ', '.join([self._variables[i].__str__() + self._varDomain[i].__str__()\
                                for i in range(len(self._varDomain))])
            redDomStr = ', '.join([self._redVariables[i].__str__() + self._redDomain[i].__str__()\
                                for i in range(len(self._redDomain))])
            caseStr = "{ " + "\n ".join([case.__str__() for case in self._body]) + " }"
            return "Domain: " + domStr + '\n' + "Reduction Domain: " + redDomStr + '\n' +\
                    self._name + "(" + varStr + ") = " +\
                    caseStr + '\n' + "Default: " + self._default.__str__()
        else:
            return self._name
        
def isAffine(expr, includeDiv = True, includeModulo = False):
    """
        Function to determine if an expression is affine or not. The input
        is a binary expression tree. It recursively checks if the left and right 
        sub expressions are affine. Determines if the entire expression is 
        affine using the following rules:

        affine +,-,* constant = affine 
        constant +,-,* affine = affine 
        affine +,- affine  = affine
        affine *,/ affine = non-affine
        non-affine operand always results in a non-affine expression

        Divisions and modulo operators are considered affine if the appropriate 
        option is specified. 

        This function is meant to work on straight forward expressions it can
        be easily tricked into conservatively saying that the expression is 
        not affine. For example -x^3 + x^3 will be non affine. 

        Making the expression anaylsis more robust will require integration 
        with a symbolic math package. 
    """
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression)
           or isinstance(expr, Condition))
    if (isinstance(expr, Value)):
        return (expr.typ is Int) or (includeDiv and (expr.typ is Rational))
    elif (isinstance(expr, Variable)):
        return True
    elif (isinstance(expr, Reference)):
        return False
    elif (isinstance(expr, AbstractBinaryOpNode)):
        leftCheck = isAffine(expr.left, includeDiv, includeModulo)
        rightCheck = isAffine(expr.right, includeDiv, includeModulo)
        # Bad coding style I suppose. Have to fix this
        if (leftCheck and rightCheck):
            if (expr.op in ['+','-']):
                return True
            elif(expr.op in ['*']):
                if(not (expr.left.has(Variable) or expr.left.has(Parameter)) or 
                    not (expr.right.has(Variable) or expr.right.has(Parameter))):
                    return True
                else:
                    return False
            elif(includeDiv and expr.op in ['/']):
                if (not (expr.right.has(Variable)) and not(expr.right.has(Parameter))):
                    return True
                else:
                    return False
            elif(includeModulo and expr.op in ['%']):
                if (not (expr.right.has(Variable)) and not(expr.right.has(Parameter))):
                    return isAffine(expr.left, includeDiv, False)
                else:
                    return False
            else:
                return False
        else:
            return False
    elif (isinstance(expr, AbstractUnaryOpNode)):
        return isAffine(expr.child, includeDiv, includeModulo)
    elif (isinstance(expr, Condition)):
        return isAffine(expr.lhs, includeDiv, includeModulo) and \
               isAffine(expr.rhs, includeDiv, includeModulo)
    elif (isinstance(expr, (Select, Cast, InbuiltFunction))):
        return False
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
        return result_type(leftType, rightType)
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

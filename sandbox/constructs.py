from fractions import Fraction
from fractions import gcd
class AbstractExpression(object):
    """ AbstractExpression class is a tree representation for expressions
        involving binary arithmetic comparison and unary operations over
        function, accumulator and image references, variables, parameters and
        constants. Users are not expected to directly create expressions, 
        the arithmetic operators are overloaded. Therefore expression in
        parameters, variables and references will automatically be converted 
        in to an expression tree.
    """
    def typeCheck(func):
        def checked(self, other):
            other = Value.numericToValue(other)
            assert isinstance(other, AbstractExpression)
            return func(self, other)
        return checked

    @typeCheck
    def __add__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return Add(self, other)
    @typeCheck
    def __radd__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return Add(other, self)

    @typeCheck
    def __sub__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return Sub(self, other)
    @typeCheck
    def __rsub__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return UnaryMinus(self)
        return Sub(other, self)


    @typeCheck
    def __mul__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
                return self
            if other.value == 0:
                return Value(0, Int)
        return Mul(self, other)
    @typeCheck
    def __rmul__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
                return self
            if other.value == 0:
                return Value(0, Int)
        return Mul(other, self)

    @typeCheck
    def __div__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
                return self
        return Div(self, other)
    @typeCheck
    def __rdiv__(self, other):
        return Div(other, self)

    @typeCheck
    def __lshift__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return LShift(self, other)
    def __rlshift__(self, other):
        return LShift(other, self)

    @typeCheck
    def __rshift__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return RShift(self, other)
    @typeCheck
    def __rrshift__(self, other):
        return RShift(other, self)

    @typeCheck
    def __mod__(self, other):
        return Mod(self, other)
    @typeCheck
    def __rmod__(self, other):
        return Mod(other, self)
    
    @typeCheck
    def __and__(self, other):
        return LAnd(self, other)
    @typeCheck
    def __rand__(self, other):
        return LAnd(other, self)
    
    @typeCheck
    def __xor__(self, other):
        return Xor(self, other)
    @typeCheck
    def __rxor__(self, other):
        return Xor(other, self)

    @typeCheck
    def __or__(self, other):
        return LOr(self, other)
    @typeCheck
    def __ror__(self, other):
        return LOr(other, self)

    def __neg__(self):
        return UnaryMinus(self)
    def __pos__(self):
        return UnaryPlus(self)

    def has(self, objType):
        return len(self.collect(objType)) > 0

    def collect(self, objType):
        objs = []
        if (type(self) is objType):
            objs = [self]
        return objs

class AbstractBinaryOpNode(AbstractExpression):
    def __init__(self, _left, _right, _op=None): 
        self._left  = _left
        self._right = _right
        self._op    = _op
    
    @property 
    def left(self):
        return self._left
    @property
    def right(self):
        return self._right
    @property
    def op(self):
        return self._op

    def collect(self, objType):
        objs = []
        if (self._left is not None):
            objs = objs + self._left.collect(objType)
        if (self._right is not None):
            objs = objs + self._right.collect(objType)
        if (type(self) == objType):
            objs = objs = [self]
        return list(set(objs))  

    def clone(self):
        return AbstractBinaryOpNode(self._left.clone(), self._right.clone(), 
                                    self._op)

    def __str__(self):        
        if (self._op is None):
            assert self._left is None and self._right is None
            return ""
        left_str = ""
        right_str = ""
        if (self._left is not None):
            left_str = self._left.__str__()
        if (self._right is not None):
            right_str = self._right.__str__()
        return "(" + left_str + " " + self._op + " " + right_str + ")"


class Add(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '+')        
class Mul(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '*')
class Div(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '/')
class Sub(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '-')
class LShift(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '<<')
class RShift(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '>>')
class Mod(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '%')
class LAnd(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '&')
class LOr(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '|')
class Xor(AbstractBinaryOpNode):
    def __init__(self, _left, _right):
        AbstractBinaryOpNode.__init__(self, _left, _right, '^')

class InbuiltFunction(AbstractExpression):
    def __init__(self, *_args):
        _args = [ Value.numericToValue(arg) for arg in _args] 
        for arg in _args:
            assert(isinstance(arg, AbstractExpression))
        self._args = _args    

    @property
    def arguments(self):
        return self._args 
   
    def getType(self):
        raise TypeError(self)

    def collect(self, objType):
        objs = []
        for arg in self._args:
            objs = objs + arg.collect(objType)
        if (type(self) == objType):
            objs = objs + [self]
        return list(set(objs))

    def substituteVars(self, varToExprMap):
        for i in xrange(0, len(self._args)):
            self._args[i] = substituteVars(self._args[i], varToExprMap)

    def inlineRefs(self, refToExprMap):
        self._args = [ substituteRefs(arg, refToExprMap) for arg in self._args]

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
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Exp(cloneArgs[0])

    def __str__(self):
        return "std::exp(" +  self._args[0].__str__() +  ")"

class Sin(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return Double

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Sin(cloneArgs[0])

    def __str__(self):
        return "std::sin(" +  self._args[0].__str__() +  ")"

class Cos(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Double

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Cos(cloneArgs[0])

    def __str__(self):
        return "std::cos(" +  self._args[0].__str__() +  ")"

class Sqrt(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Double

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Sqrt(cloneArgs[0])

    def __str__(self):
        return "std::sqrt(" +  self._args[0].__str__() +  ")"

class Sqrtf(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Float

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Sqrtf(cloneArgs[0])

    def __str__(self):
        return "std::sqrtf(" +  self._args[0].__str__() +  ")"

class Abs(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return getType(self._args[0])

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Abs(cloneArgs[0])

    def __str__(self):
        return "std::abs(" +  self._args[0].__str__() +  ")"

class Cast(AbstractExpression):
    def __init__(self, _typ, _expr):
        _expr = Value.numericToValue(_expr)
        assert _typ in [Float, Int, UInt, UChar, Char, Double, Long, UInt, Short, UShort] 
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
        if (self._expr is not None):
            objs = objs + self._expr.collect(objType)
        if (type(self) == objType):
            objs = objs + [self]
        return list(set(objs))

    def clone(self):
        return Cast(self._typ, self._expr.clone())
    
    def inlineRefs(self, refToExprMap):
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
        if (self._cond is not None):
            objs = objs + self._cond.collect(objType)
        if (self._trueExpr is not None):
            objs = objs + self._trueExpr.collect(objType)
        if (self._falseExpr is not None):
            objs = objs + self._falseExpr.collect(objType)
        if (type(self) == objType):
            objs = objs + [self]
        return list(set(objs))

    def inlineRefs(self, refToExprMap):
        self._cond.inlineRefs(refToExprMap)
        self._trueExpr = substituteRefs(self._trueExpr, refToExprMap)
        self._falseExpr = substituteRefs(self._falseExpr, refToExprMap)
    
    def clone(self):
        return Select(self._cond.clone(), self._trueExpr.clone(), 
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

class AbstractUnaryOpNode(AbstractExpression):
    def __init__(self, _child, _op=None): 
        self._child  = _child
        self._op    = _op
    
    @property 
    def child(self):
        return self._child
    @property
    def op(self):
        return self._op
        
    def collect(self, objType):
        objs = []
        if (self._child is not None):
            objs = objs + self._child.collect(objType)
        if (type(self) == objType):
            objs = objs + [self]
        return list(set(objs))

    def clone(self):
        return AbstractUnaryOpNode(self._child.clone(), self._op)

    def __str__(self):        
        if (self._op is None):
            assert self._child is None
            return ""
        child_str = ""
        if (self._child is not None):
            child_str = self._child.__str__()
        return self._op + "(" + child_str + ")"

class UnaryPlus(AbstractUnaryOpNode):
    def __init__(self, _child):
        AbstractUnaryOpNode.__init__(self, _child, '+')

class UnaryMinus(AbstractUnaryOpNode):
    def __init__(self, _child):
        AbstractUnaryOpNode.__init__(self, _child, '-')

class Value(AbstractExpression):

    @classmethod
    def numericToValue(cls, _value):
        if type(_value) is int:
            _value = Value(_value, Int)
        elif type(_value) is float:
            _value = Value(_value, Float)
        elif type(_value) is long:
            _value = Value(_value, Long)
        elif type(_value) is Fraction:
            _value = Value(_value, Rational)
        return _value

    def __init__(self, _value, _typ):
        self._value = _value
        self._typ = _typ

    def clone(self):
        return Value(self._value, self._typ)

    @property
    def value(self):
        return self._value
    @property
    def typ(self):
        return self._typ
    def __eq__(self, _value):
        return self._value == _value
    def __str__(self):
        return self._value.__str__()

class Float(object):
    @staticmethod
    def cTypeName():
        return 'float'
    pass
class Double(object):
    @staticmethod
    def cTypeName():
        return 'double'
    pass
class Int(object):
    @staticmethod
    def cTypeName():
        return 'int32'
    pass
class UInt(object):
    @staticmethod
    def cTypeName():
        return 'uint32'
    pass
class Short(object):
    @staticmethod
    def cTypeName():
        return 'int16'
    pass
class UShort(object):
    @staticmethod
    def cTypeName():
        return 'uint64'
    pass
class Char(object):
    @staticmethod
    def cTypeName():
        return 'int8'
    pass
class UChar(object):
    @staticmethod
    def cTypeName():
        return 'uint8'
    pass
class Long(object):
    @staticmethod
    def cTypeName():
        return 'int64'
    pass
class ULong(object):
    @staticmethod
    def cTypeName():
        return 'float'
    pass
class Rational(object):
    pass

def typeImplictConversion(a, b):
    if a is Double or b is Double:
        return Double
    elif a is Float or b is Float:
        return Float
    elif a is ULong or b is ULong:
        return ULong
    elif a is Long or b is Long:
        return Long
    elif a is UInt or b is UInt:
        return UInt
    elif a is Int or b is Int:
        return Int
    elif a is UShort or b is UShort:
        return UShort
    elif a is Short or b is Short:
        return Short
    elif a is UChar or b is UChar:
        return UChar
    elif a is Char and b is Char:
        return Char
    raise TypeError((a, b))

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
        self._definition = None

    @property
    def definition(self):
        return self._definition

    @definition.setter
    def definition(self, _def):
        assert(self._definition is None) 
        _def = Value.numericToValue(_def)
        assert(isinstance(_def, AbstractExpression))
        self._definition = _def
    
class Interval(object):
    def __init__(self, _typ,  _lb, _ub, _step=1):
        _lb   = Value.numericToValue(_lb)
        _ub   = Value.numericToValue(_ub)
        _step = Value.numericToValue(_step)
        assert(isinstance(_lb, AbstractExpression))
        assert(isinstance(_ub, AbstractExpression))        
        assert(isinstance(_step, AbstractExpression))
        assert(type(_step.value) == int)
        self._lb   = _lb
        self._ub   = _ub
        self._step = _step
        self._typ  = _typ
 
    @property 
    def lowerBound(self):
        return self._lb
    @property 
    def upperBound(self):
        return self._ub
    @property 
    def step(self):
        return self._step
    @property 
    def typ(self):
        return self._typ

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._lb.collect(objType)        
        objs = objs + self._ub.collect(objType)
        objs = objs + self._step.collect(objType)
        return list(set(objs))    
    
    def __str__(self):
        return '(' + self._lb.__str__() + ', ' +\
               self._ub.__str__() + ', ' + self._step.__str__() + ')'

class Domain(object):
    def __init__(self, _intervals):
        self._intervals = _intervals

    @property
    def intervals(self):
        return tuple(self._intervals)

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = []
        for interval in self._intervals:
            objs = objs +  interval.collect(objType)
        return list(set(objs))    

    def __str__(self):
        interval_str = ", ".join([interval.__str__() for interval in self._intervals])
        return "(" + interval_str + ")"

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
    @property
    def arguments(self):
        return self._args

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args]
        return Reference(self._obj, cloneArgs)

    def collect(self, objType):
        objs = []
        for arg in self._args:     
            objs = objs + arg.collect(objType)
        if (type(self) is objType):
            objs = objs + [self]
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

    def inlineRefs(self, refToExprMap):
        if(isinstance(self._left, Condition)):
            self._left.inlineRefs(refToExprMap) 
        else:
            self._left = substituteRefs(self._left, refToExprMap) 
        if(isinstance(self._right, Condition)):
            self._right.inlineRefs(refToExprMap)
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
        left_str = ""
        right_str = ""
        if (self._left is not None):
            left_str = self._left.__str__()
        if (self._right is not None):
            right_str = self._right.__str__()
        return "(" + left_str + " " + self._cond + " " + right_str + ")"

class Case(object):
    def __init__(self, _cond, _expr):
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_cond, Condition))
        assert(isinstance(_expr, (AbstractExpression, Accumulate)))
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

    def inlineRefs(self, refToExprMap):
        self._cond.inlineRefs(refToExprMap)
        self._expr = substituteRefs(self._expr, refToExprMap)
    
    def __str__(self):
        return 'Case(' + self._cond.__str__() + ')' +\
                '{ ' + self._expr.__str__() + ' }'

class Op(object):
    Sum = 0
    Min = 1
    Max = 2    
    
class Accumulate(object):
    def __init__(self, _accRef, _expr, _opTyp):
        assert isinstance(_accRef, Reference)
        assert isinstance(_accRef.objectRef, Accumulator)
        _expr = Value.numericToValue(_expr)
        assert isinstance(_expr, AbstractExpression)
        assert _opTyp in [Op.Sum, Op.Min, Op.Max]

        self._accRef = _accRef
        self._expr = _expr        
        self._opTyp = _opTyp

    @property 
    def accumulateRef(self):
        return self._accRef
    @property 
    def expression(self):
        return self._expr

    def inlineRefs(self, refToExprMap):
        self._expr = substituteRefs(self._expr, refToExprMap)     
        self._accRef = substituteRefs(self._accRef, refToExprMap)     
    
    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        
        objs = self._accRef.collect(objType) + self._expr.collect(objType)
        return list(set(objs))  

    def __str__(self):
        opStr = None
        if (self._opTyp == Op.Sum):
            opStr = '+'
        elif (self._opTyp == Op.Min):
            opStr = 'Min'
        elif (self._opTyp == Op.Max):
            opStr = 'Max'
        else:
            assert False
        return 'Accumulate(' + self._accRef.__str__() + ' = ' + \
                 self._accRef.__str__() + ',' + self._expr.__str__() + \
                 ' ' + opStr + ')' 

class Function(object):
    def __init__(self, _typ, _name):
        self._name      = _name
        # Type of the scalar range of the function
        self._typ       = _typ
        # Variables of the function 
        self._variables = None

        # Gives the domain of each variable. Domain of each variable is expected
        # to be over integers. Function evaluation in the lexicographic order of 
        # the domain is assumed to be valid.
        self._varDomain    = None 

        # * Body of a function is composed of Case and Expression constructs. 
        # * The Case constructs are expected to be non-overlapping.
        # * If multiple Case constructs are satisfied at a point the value 
        #   can be defined by any one of them
        # * All the Case constructs followed by an Expression construct are 
        #   ignored.
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
    @variableDomain.setter
    def variableDomain(self, _varDom):
        assert(self._variables is None)
        assert(self._varDomain is None)
        assert(len(_varDom[0]) == len(_varDom[1]))
        for i in xrange(0, len(_varDom[0])):
            assert(isinstance(_varDom[0][i], Variable))
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
    
    @property
    def domain(self):
        return self._varDomain

    @property
    def variables(self):
        return self._variables

    @property
    def definition(self):
        return self._body
    @definition.setter
    def definition(self, _def):
        _def = Value.numericToValue(_def)
        assert(isinstance(_def, (Case, AbstractExpression)))
        # check if the Case and Expression constructs only use
        # function variables and global parameters
        if(isinstance(_def, Case)):
            for part in self._body:
               assert(isinstance(part, Case))
        self._body.append(_def)

    def __call__(self, *args):
        assert(len(args) == len(self._variables))
        for arg in args:
            arg = Value.numericToValue(arg)
            assert(isinstance(arg, AbstractExpression))
        return Reference(self, args)

    def inlineRefs(self, refToExprMap):
        numCases = len(self._body)
        for i in xrange(0, numCases):
            if isinstance(self._body[i], Case):
                self._body[i].inlineRefs(refToExprMap)
            else:
                self._body[i] = substituteRefs(self._body[i], refToExprMap)

    def getObjects(self, objType):
        objs = []
        for case in self._body:
            objs = objs + case.collect(objType)
        for interval in self._varDomain:
            objs = objs + interval.collect(objType)
        return list(set(objs))

    def hasBoundedIntegerDomain(self):
        boundedIntegerDomain = True
        for varDom in self._varDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False
        return boundedIntegerDomain

    def __str__(self):
        if (self._body):
            var_str = ", ".join([var.__str__() for var in self._variables])
            dom_str = ', '.join([self._variables[i].__str__() + self._varDomain[i].__str__()\
                                for i in xrange(len(self._varDomain))])
            case_str = "{ " + "\n ".join([case.__str__() for case in self._body]) + " }"
            return "Domain: " + dom_str + '\n' + self._name + "(" + var_str + ") = " +\
                    case_str + '\n'
        else:
            return self._name

class Image(Function):
    def __init__(self, _typ, _name, _dims):
        Function.__init__(self, _typ, _name)
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
            intervals.append(Interval(UInt, 0, dim-1, 1))
            variables.append(Variable(UInt, "_" + _name + str(i)))
            i = i + 1
        self.variableDomain = (variables, intervals)

    @property
    def dimensions(self):
        return tuple(self._dims)

    def __str__(self):
        dim_str = ", ".join([dim.__str__() for dim in self._dims])
        return self._name.__str__() + "(" + dim_str + ")"

class Accumulator(Function):
    def __init__(self, _typ, _name):
        Function.__init__(self, _typ, _name)  
        # Gives the domain of the reduction. Reduction domain of each variable is 
        # expected to be over integers. Reduction evaluation in the lexicographic 
        # order of the domain is assumed to be valid.
        self._redDomain = None
        self._redVariables = None
        # Intial value of each accumulator cell. Default is set to zero of the given 
        # type
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

    @reductionDomain.setter
    def reductionDomain(self, _redDom):
        assert(self._redVariables is None)
        assert(self._redDomain is None)
        assert(len(_redDom[0]) == len(_redDom[1]))
        for i in xrange(0, len(_redDom[0])):
            assert(isinstance(_redDom[0][i], Variable))
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

    @property
    def defintion(self):
        return self._body

    @defintion.setter
    def definition(self, _def):
        _def = Value.numericToValue(_def)
        assert(isinstance(_def, (Case, Accumulate)))
        # check if the Case and Expression constructs only use
        # function variables and global parameters

        # Which way is better Case inside accumulate or accumulate inside Case
        if(isinstance(_def, Case)):
            for part in self._body:
               assert(isinstance(_def.expression, Accumulate))
               assert(isinstance(part, Case))
        self._body.append(_def)

    def hasBoundedIntegerDomain(self):
        boundedIntegerDomain = True
        for varDom in self._varDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False

        for varDom in self._redDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False
        
        return boundedIntegerDomain
    
    def getObjects(self, objType):
        objs = []
        for case in self._body:
            objs = objs + case.collect(objType)
        for interval in self._varDomain:
            objs = objs + interval.collect(objType)
        for interval in self._redDomain:
            objs = objs + interval.collect(objType)
        objs = objs + self._default.collect(objType)
        return list(set(objs))
    
    def __str__(self):
        if (self._body):
            varStr = ", ".join([var.__str__() for var in self._variables])
            domStr = ', '.join([self._variables[i].__str__() + self._varDomain[i].__str__()\
                                for i in xrange(len(self._varDomain))])
            redDomStr = ', '.join([self._redVariables[i].__str__() + self._redDomain[i].__str__()\
                                for i in xrange(len(self._redDomain))])
            caseStr = "{ " + "\n ".join([case.__str__() for case in self._body]) + " }"
            return "Domain: " + domStr + '\n' + "Reduction Domain: " + redDomStr + '\n' +\
                    self._name + "(" + varStr + ") = " +\
                    caseStr + '\n' + "Default: " + self._default.__str__()
        else:
            return self._name

def isAffine(expr, div = True, modulo = False):
    expr = Value.numericToValue(expr)
    assert(isinstance(expr, AbstractExpression)
           or isinstance(expr, Condition))
    if (isinstance(expr, Value)):
        return (expr.typ == Int) or ((expr.typ == Rational) and div)    
    elif (isinstance(expr, Variable)):
        return True
    elif (isinstance(expr, Reference)):
        return False    
    elif (isinstance(expr, AbstractBinaryOpNode)):
        leftCheck = isAffine(expr.left, div, modulo)
        rightCheck = isAffine(expr.right, div, modulo)
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
            elif(expr.op in ['/'] and div):
                if (not (expr.right.has(Variable) or expr.right.has(Parameter))):
                    return True
                else:
                    return False
            elif(expr.op in ['%'] and modulo):
                if (not (expr.right.has(Variable) or expr.right.has(Parameter))):
                    return isAffine(expr.left, div, False)
                else:
                    return False
            else:
                return False
        else:
            return False
    elif (isinstance(expr, AbstractUnaryOpNode)):
        return isAffine(expr.child, div, modulo)
    elif (isinstance(expr, Condition)):
        return isAffine(expr.lhs, div, modulo) and \
               isAffine(expr.rhs, div, modulo)
    elif (isinstance(expr, (Select, Cast, InbuiltFunction))):
        return False
    raise TypeError(type(expr))

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
            leftCheck = isConstantExpr(expr.left, affine = True)
            rightCheck = isConstantExpr(expr.right, affine = True)
            #sanity check should be true if the expression is affine
            assert(not (leftCheck and rightCheck))           
            if (leftCheck and not rightCheck):
                coeff = dict( (n, getConstantFromExpr(expr.left, affine = True) * 
                                  rightCoeff.get(n, 0))
                              for n in set(rightCoeff) )
            elif(rightCheck and not leftCheck):
                coeff = dict( (n, getConstantFromExpr(expr.right, affine = True) * 
                                  leftCoeff.get(n, 0))
                              for n in set(leftCoeff) )
        elif (expr.op == '/'):
            rightCheck = isConstantExpr(expr.right, affine = True)
            #sanity check should be true if the expression is affine
            assert(rightCheck)           
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
            assert (expr.typ == Int) or (expr.typ == Rational)
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
            return (expr.typ == Int) or (expr.typ == Rational)
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
        if (isAffine(expr, div=False)):
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
            for i in xrange(0, len(refVars)):
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
        for i in xrange(0, numArgs):
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

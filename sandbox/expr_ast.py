#
# Copyright 2014-2016 Vinay Vasista, Ravi Teja Mullapudi, Uday Bondhugula,
# and others from Multicore Computing Lab, Department of Computer Science
# and Automation, Indian Institute of Science
#

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# expr_ast.py : IR for PolyMage expressions.
#

from __future__ import absolute_import, division, print_function

from fractions import Fraction
from fractions import gcd

# TODO remove this at some point
from expr_types import *

class AbstractExpression(object):
    """ AbstractExpression class is a tree representation for expressions
    involving binary arithmetic comparison and unary operations over function,
    reduction and image references, variables, parameters and constants. Users
    are not expected to directly create expressions.  The arithmetic operators
    are overloaded therefore expressions in parameters, variables and
    references will automatically be converted in to an expression tree.
    """
    def type_check(func):
        def checked(self, other):
            other = Value.numericToValue(other)
            assert isinstance(other, AbstractExpression)
            return func(self, other)
        return checked

    @type_check
    def __add__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return Add(self, other)
    @type_check
    def __radd__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return Add(other, self)

    @type_check
    def __sub__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return Sub(self, other)
    @type_check
    def __rsub__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return UnaryMinus(self)
        return Sub(other, self)


    @type_check
    def __mul__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
                return self
            if other.value == 0:
                return Value(0, Int)
        return Mul(self, other)
    @type_check
    def __rmul__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
               return self
            if other.value == 0:
                return Value(0, Int)
        return Mul(other, self)

    @type_check
    def __truediv__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
                return self
        return Div(self, other)

    @type_check
    def __rtruediv__(self, other):
        return Div(other, self)

    @type_check
    def __floordiv__(self, other):
        if isinstance(other, Value):
            if other.value == 1:
                return self
        return Div(self, other)

    @type_check
    def __rfloordiv__(self, other):
        return Div(other, self)

    @type_check
    def __lshift__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return LShift(self, other)
    def __rlshift__(self, other):
        return LShift(other, self)

    @type_check
    def __rshift__(self, other):
        if isinstance(other, Value):
            if other.value == 0:
                return self
        return RShift(self, other)
    @type_check
    def __rrshift__(self, other):
        return RShift(other, self)

    @type_check
    def __mod__(self, other):
        return Mod(self, other)
    @type_check
    def __rmod__(self, other):
        return Mod(other, self)
    
    @type_check
    def __and__(self, other):
        return LAnd(self, other)
    @type_check
    def __rand__(self, other):
        return LAnd(other, self)
    
    @type_check
    def __xor__(self, other):
        return Xor(self, other)
    @type_check
    def __rxor__(self, other):
        return Xor(other, self)

    @type_check
    def __or__(self, other):
        return LOr(self, other)
    @type_check
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

class Value(AbstractExpression):

    @classmethod
    def numericToValue(cls, _value):
        if type(_value) is int:
            _value = Value(_value, Int)
        elif type(_value) is float:
            _value = Value(_value, Float)
        # Python 3 has no long
        #elif type(_value) is long:
        #    _value = Value(_value, Long)
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
        objs += self._left.collect(objType)
        objs += self._right.collect(objType)
        if (type(self) == objType):
            objs += [self]
        return list(set(objs))  

    def clone(self):
        return AbstractBinaryOpNode(self._left.clone(), 
                                    self._right.clone(), 
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
   
    def get_type(self):
        raise TypeError(self)

    def collect(self, objType):
        objs = []
        for arg in self._args:
            objs = objs + arg.collect(objType)
        if (type(self) == objType):
            objs = objs + [self]
        return list(set(objs))

    def sub_vars(self, varToExprMap):
        for i in range(0, len(self._args)):
            self._args[i] = sub_vars(self._args[i], varToExprMap)

    def inline_refs(self, refToExprMap):
        self._args = [ substitute_refs(arg, refToExprMap) for arg in self._args]

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
        objs += self._child.collect(objType)
        if (type(self) == objType):
            objs += [self]
        return list(set(objs))

    def clone(self):
        return AbstractUnaryOpNode(self._child.clone(), self._op)

    def __str__(self):
        child_str = self._child.__str__()
        return self._op + "(" + child_str + ")"

class UnaryPlus(AbstractUnaryOpNode):
    def __init__(self, _child):
        AbstractUnaryOpNode.__init__(self, _child, '+')

class UnaryMinus(AbstractUnaryOpNode):
    def __init__(self, _child):
        AbstractUnaryOpNode.__init__(self, _child, '-')

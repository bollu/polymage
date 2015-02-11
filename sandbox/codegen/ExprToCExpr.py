# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

import sys
sys.path.insert(0, '../frontend')

from Constructs import *

class cValue(Value):
    def __init__(self, _value, _typ):
        Value.__init__(self, _value, _typ)
    def __str__(self):
        if (self._typ == Float):
            return self._value.__str__() + 'f'
        return self._value.__str__()

class cSelect(Select):
    def __init__(self, _cCond, _trueCExpr, _falseCExpr):
        Select.__init__(self, _cCond, _trueCExpr, _falseCExpr, False)

class cMax(Max):
    def __init__(self, _expr1, _expr2):
        Max.__init__(self, _expr1, _expr2, False)

class cMin(Min):
    def __init__(self, _expr1, _expr2):
        Min.__init__(self, _expr1, _expr2, False)

class cPow(Pow):
    def __init__(self, _expr1, _expr2):
        Pow.__init__(self, _expr1, _expr2)

class cPowf(Powf):
    def __init__(self, _expr1, _expr2):
        Powf.__init__(self, _expr1, _expr2)

class cExp(Exp):
    def __init__(self, _expr):
        Exp.__init__(self, _expr)

class cSin(Sin):
    def __init__(self, _expr):
        Sin.__init__(self, _expr)

class cCos(Cos):
    def __init__(self, _expr):
        Cos.__init__(self, _expr)

class cSqrt(Sqrt):
    def __init__(self, _expr):
        Sqrt.__init__(self, _expr)

class cSqrtf(Sqrtf):
    def __init__(self, _expr):
        Sqrtf.__init__(self, _expr)

class cAbs(Abs):
    def __init__(self, _expr):
        Abs.__init__(self, _expr)


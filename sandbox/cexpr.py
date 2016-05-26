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
# cexpr.py : C Expression classes
#

from __future__ import absolute_import, division, print_function

from constructs import *


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


class cRandomFloat(RandomFloat):
    def __init__(self):
        RandomFloat.__init__(self)


class cLog(Log):
    def __init__(self, _expr):
        Log.__init__(self, _expr)


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

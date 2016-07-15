from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../')

from .constructs import *
from .expression import *

def test_affine():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    assert(isAffine(x + y) == True)
    assert(isAffine(3) == True)
    assert(isAffine(x*y) == False)
    assert(isAffine(-x + N + 3*y) == True)
    assert(isAffine(2*x + N/2 + 3*y) == True)
    c1 = Condition(x, '<', 2*y)
    c2 = Condition(x, '>', 2-y)
    c3 = Condition(x, '>=', x*y)
    c4 = Condition(x + 2*N, '<=', y + N)
    c5 = Condition(x*N, '!=', y)
    assert(isAffine(c1) == True)
    assert(isAffine(c2) == True)
    assert(isAffine(c3) == False)
    assert(isAffine(c4) == True)
    assert(isAffine(c5) == False)

def test_coeff():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    coeff = get_affine_var_and_param_coeff(1+x)
    assert(coeff[x] == 1)
    coeff = get_affine_var_and_param_coeff(1+x +y)
    assert(coeff[x] == 1 and coeff[y] == 1)
    coeff = get_affine_var_and_param_coeff(3)
    assert(coeff == {})
    coeff = get_affine_var_and_param_coeff(N*x + y)
    assert(coeff == {})
    coeff = get_affine_var_and_param_coeff(x*y)
    assert(coeff == {})
    coeff = get_affine_var_and_param_coeff(2*(x*3+y +N +x + y -5) 
                                      + 3*(-x) + 4*(-y) + N)
    assert(coeff[x] == 5 and coeff[y] == 0 and coeff[N] == 3)


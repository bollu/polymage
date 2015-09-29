from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../../frontend')

from constructs import *
from expression import *

def test_simplify():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")

    expr = (2*x + 3*x + 3 + 4 + 1 + 5*y - 4*y + N + 2 * 3 + 4//2)//2 
    expr = simplifyExpr(expr)
    coeff = get_affine_var_and_param_coeff(expr)
    
    assert coeff[x] == Fraction(5, 2)
    assert coeff[y] == Fraction(1, 2)
    assert coeff[N] == Fraction(1, 2)
    assert get_constant_from_expr(expr, affine = True) == 8


from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def interpolate(g, u, l, impipeDict, name):
    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    innerBox = interior[l]['innerBox']

    uu = Function(([z, y, x], \
                   [extent[l], extent[l], extent[l]]), \
                   Double, \
                   str(name))

    zz = z/2
    yy = y/2
    xx = x/2

    def z1(xx):
        return g(zz  , yy+1, xx) + g(zz  , yy  , xx)
    def z2(xx):
        return g(zz+1, yy  , xx) + g(zz  , yy  , xx)
    def z3(xx):
        return g(zz+1, yy+1, xx) + g(zz+1, yy, xx) \
             + g(zz  , yy+1, xx) + g(zz  , yy, xx)

    expr_000 = g(zz, yy, xx)
    expr_001 = 0.500 * (g(zz, yy, xx) + g(zz, yy, xx+1))
    expr_010 = 0.500 * z1(xx)
    expr_011 = 0.250 * (z1(xx) + z1(xx+1))
    expr_100 = 0.500 * z2(xx)
    expr_101 = 0.250 * (z2(xx) + z2(xx+1))
    expr_110 = 0.250 * z3(xx)
    expr_111 = 0.125 * (z3(xx) + z3(xx+1))

    even_x = Condition(x%2, '==', 0)
    even_y = Condition(y%2, '==', 0)
    even_z = Condition(z%2, '==', 0)

    odd_x  = Condition(x%2, '==', 1)
    odd_x  = Condition(y%2, '==', 1)
    odd_x  = Condition(z%2, '==', 1)

    if u == None:
        correct = 0.0
    else:
        correct = u(z, y, x)

    uu.defn = [ correct + \
                Select(even_z,
                  Select(even_y,
                    Select(even_x,
                      expr_000,
                      expr_001),
                    Select(even_x,
                      expr_010,
                      expr_011)),
                  Select(even_y,
                    Select(even_x,
                      expr_100,
                      expr_101),
                    Select(even_x,
                      expr_110,
                      expr_111))) ]

    return uu

from __init__ import *

import sys
sys.path.insert(0, ROOT)

from fractions  import Fraction
from compiler import *
from constructs import *

def interpolate(G, U, l, pipe_data, name):
    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l]['inner_box']

    UU = Function(([z, y, x], [extent[l], extent[l], extent[l]]),
                   Double, str(name))

    zz = z/2
    yy = y/2
    xx = x/2

    def z1(xx):
        return G(zz  , yy+1, xx) + G(zz  , yy  , xx)
    def z2(xx):
        return G(zz+1, yy  , xx) + G(zz  , yy  , xx)
    def z3(xx):
        return G(zz+1, yy+1, xx) + G(zz+1, yy, xx) \
             + G(zz  , yy+1, xx) + G(zz  , yy, xx)

    expr_000 = G(zz, yy, xx)
    expr_001 = 0.500 * (G(zz, yy, xx) + G(zz, yy, xx+1))
    expr_010 = 0.500 * z1(xx)
    expr_011 = 0.250 * (z1(xx) + z1(xx+1))
    expr_100 = 0.500 * z2(xx)
    expr_101 = 0.250 * (z2(xx) + z2(xx+1))
    expr_110 = 0.250 * z3(xx)
    expr_111 = 0.125 * (z3(xx) + z3(xx+1))

    even_x = Condition(x%2, '==', 0)
    even_y = Condition(y%2, '==', 0)
    even_z = Condition(z%2, '==', 0)

    if U == None:
        correct = 0.0
    else:
        correct = U(z, y, x)

    UU.defn = [ correct + \
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

    return UU

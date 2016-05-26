from __init__ import *
import sys
from polymage_common import set_zero_ghosts

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def restrict(R, l, pipe_data, name):
    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l-1]['inner_box']

    S = Function(([z, y, x], [extent[l-1], extent[l-1], extent[l-1]]),
                 Double, str(name))

    zz = 2*z
    yy = 2*y
    xx = 2*x

    def x1(xx):
        return (R(zz  , yy-1, xx) + R(zz  , yy+1, xx)
              + R(zz-1, yy  , xx) + R(zz+1, yy  , xx))

    def y1(xx):
        return (R(zz-1, yy-1, xx) + R(zz-1, yy+1, xx)
              + R(zz+1, yy-1, xx) + R(zz+1, yy+1, xx))

    S.defn = [ Case(inner_box,
                    0.5000 * R(zz, yy, xx)
                  + 0.2500 * (R(zz  , yy  , xx-1)
                            + R(zz  , yy  , xx+1)
                            + x1(xx  ))
                  + 0.1250 * (x1(xx-1)
                            + x1(xx+1)
                            + y1(xx  ))
                  + 0.0625 * (y1(xx-1)
                            + y1(xx+1))
                  ) ]

    set_zero_ghosts(S, ghosts[l-1])

    return S



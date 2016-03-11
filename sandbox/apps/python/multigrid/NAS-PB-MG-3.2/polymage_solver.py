import sys
from polymage_common import set_zero_ghosts

sys.path.insert(0, '../../../../')
sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def psinv(R, U, l, app_data, name):
    pipe_data = app_data['pipe_data']

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    c = app_data['c']

    inner_box = interior[l]['inner_box']

    W = Function(([z, y, x], [extent[l], extent[l], extent[l]]),
                 Double, str(name))

    def r1(x):
        return (R(z  , y-1, x) + R(z  , y+1, x)
              + R(z-1, y  , x) + R(z+1, y  , x))

    def r2(x):
        return (R(z-1, y-1, x) + R(z-1, y+1, x)
              + R(z+1, y-1, x) + R(z+1, y+1, x))

    if U == None:
        u = 0.0
    else:
        u = U(z, y, x)

    W.defn = [ Case(inner_box,
                   u \
                 + c[0] * R(z, y, x)
                 + c[1] * (R(z, y, x-1) + R(z, y, x+1) + r1(x))
                 + c[2] * (r2(x) + r1(x-1) + r1(x+1))
#                + c[3] * (r2(x-1) + r2(x+1))
                 ) ]

    set_zero_ghosts(W, ghosts[l])

    return W

import sys
from polymage_common import set_zero_ghosts

sys.path.insert(0, '../../../../')
sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def residual(U, V, l, app_data, name):
    pipe_data = app_data['pipe_data']

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    a = app_data['a']

    inner_box = interior[l]['inner_box']

    R = Function(([z, y, x], [extent[l], extent[l], extent[l]]),
                 Double, str(name))

    def u1(x):
        return (U(z  , y-1, x) + U(z  , y+1, x)
              + U(z-1, y  , x) + U(z+1, y  , x))

    def u2(x):
        return (U(z-1, y-1, x) + U(z-1, y+1, x)
              + U(z+1, y-1, x) + U(z+1, y+1, x))

    R.defn = [ Case(inner_box,
                   V(z, y, x) \
                 - a[0] * U(z, y, x)
#                - a[1] * (U(z, y, x-1) + U(z, y, x+1) + u1(x))
                 - a[2] * (u2(x) + u1(x-1) + u1(x+1))
                 - a[3] * (u2(x-1) + u2(x+1))
                 ) ]

    set_zero_ghosts(r, ghosts[l])

    return r

def resid_pipe(app_data):
    pipe_data = app_data['pipe_data']

    n = pipe_data['n']

    u = Image(Double, "u", [n, n, n])
    v = Image(Double, "v", [n, n, n])

    lt = app_data['lt']

    r = residual(u, v, lt, pipe_data, app_data, "resid")

    return r

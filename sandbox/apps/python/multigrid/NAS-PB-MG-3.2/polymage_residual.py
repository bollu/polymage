import sys
from polymage_common import setZeroGhosts

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def residual(u, v, l, impipeDict, dataDict, name):
    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    a = dataDict['a']

    innerBox = interior[l]['innerBox']

    r = Function(([z, y, x], \
                   [extent[l], extent[l], extent[l]]), \
                   Double, \
                   str(name))

    def u1(x):
        return (u(z  , y-1, x) + u(z  , y+1, x)
              + u(z-1, y  , x) + u(z+1, y  , x))

    def u2(x):
        return (u(z-1, y-1, x) + u(z-1, y+1, x)
              + u(z+1, y-1, x) + u(z+1, y+1, x))

    r.defn = [ Case(innerBox,
                   v(z, y, x)
                 - a[0] * u(z, y, x)
#                - a[1] * (u(z, y, x-1) + u(z, y, x+1) + u1(x))
                 - a[2] * (u2(x) + u1(x-1) + u1(x+1))
                 - a[3] * (u2(x-1) + u2(x+1))
                 ) ]

    setZeroGhosts(r, ghosts[l])

    return r

def residPipe(impipeDict, dataDict):
    n = impipeDict['n']

    u = Image(Double, "u", [n, n, n])
    v = Image(Double, "v", [n, n, n])

    lt = dataDict['lt']

    r = residual(u, v, lt, impipeDict, dataDict, "resid")

    return r

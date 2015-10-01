import sys
from polymage_common import setZeroGhosts

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def psinv(r, u_, l, impipeDict, dataDict, name):
    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    c = dataDict['c']

    innerBox = interior[l]['innerBox']

    w = Function(([z, y, x], \
                   [extent[l], extent[l], extent[l]]), \
                   Double, \
                   str(name))

    def r1(x):
        return (r(z  , y-1, x) + r(z  , y+1, x)
              + r(z-1, y  , x) + r(z+1, y  , x))

    def r2(x):
        return (r(z-1, y-1, x) + r(z-1, y+1, x)
              + r(z+1, y-1, x) + r(z+1, y+1, x))

    if u_ == None:
        u = 0.0
    else:
        u = u_(z, y, x)

    w.defn = [ Case(innerBox,
                   u
                 + c[0] * r(z, y, x)
                 + c[1] * (r(z, y, x-1) + r(z, y, x+1) + r1(x))
                 + c[2] * (r2(x) + r1(x-1) + r1(x+1))
#                + c[3] * (u2(x-1) + u2(x+1))
                 ) ]

    setZeroGhosts(w, ghosts[l])

    return w

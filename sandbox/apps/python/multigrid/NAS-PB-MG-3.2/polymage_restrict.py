import sys
from polymage_common import setZeroGhosts

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def restrict(r, l, impipeDict, name):
    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    innerBox = interior[l-1]['innerBox']

    s = Function(([z, y, x], \
                   [extent[l-1], extent[l-1], extent[l-1]]), \
                   Double, \
                   str(name))

    zz = 2*z
    yy = 2*y
    xx = 2*x

    def x1(xx):
        return (r(zz  , yy-1, xx) + r(zz  , yy+1, xx)
              + r(zz-1, yy  , xx) + r(zz+1, yy  , xx))

    def y1(xx):
        return (r(zz-1, yy-1, xx) + r(zz-1, yy+1, xx)
              + r(zz+1, yy-1, xx) + r(zz+1, yy+1, xx))

    s.defn = [ Case(innerBox,
                   0.5000 * r(zz, yy, xx)
                 + 0.2500 * (r(zz  , yy  , xx-1)
                           + r(zz  , yy  , xx+1)
                           + x1(xx  ))
                 + 0.1250 * (x1(xx-1)
                           + x1(xx+1)
                           + y1(xx  ))
                 + 0.0625 * (y1(xx-1)
                           + y1(xx+1))
                 )]

    setZeroGhosts(s, ghosts[l-1])

    return s



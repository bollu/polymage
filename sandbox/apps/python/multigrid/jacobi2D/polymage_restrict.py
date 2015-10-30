from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')
sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def restrict(U_, l, name, impipeDict):
    y = impipeDict['y']
    x = impipeDict['x']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    innerBox = interior[l-1]['innerBox']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))
    W_.defn = [ Case(innerBox,
# corners
                     (U_(2*y-1, 2*x-1)           \
                    + U_(2*y-1, 2*x+1)           \
                    + U_(2*y+1, 2*x-1)           \
                    + U_(2*y+1, 2*x+1)) * 0.0625 \
# edge centers
                    +(U_(2*y-1, 2*x  )           \
                    + U_(2*y+1, 2*x  )           \
                    + U_(2*y  , 2*x-1)           \
                    + U_(2*y  , 2*x+1)) * 0.125  \
# center
                    + U_(2*y  , 2*x  )  * 0.25) ]

    setGhosts(W_, ghosts[l-1], 0.0)

    return W_

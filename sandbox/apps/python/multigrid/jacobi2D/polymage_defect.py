import sys
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def defect(U_, F_, l, name, impipeData):
    if U_ == None:
        return F_

    y = impipeData['y']
    x = impipeData['x']

    invhh = impipeData['invhh']

    extent = impipeData['extent']
    interior = impipeData['interior']
    ghosts = impipeData['ghosts']

    innerBox = interior[l]['innerBox']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))
    W_.defn = [ Case(innerBox,
                      F_(y, x)
                   - (U_(y  , x  ) * 4.0 \
                   -  U_(y-1, x  )       \
                   -  U_(y+1, x  )       \
                   -  U_(y  , x-1)       \
                   -  U_(y  , x+1)       \
                     ) * invhh[l]) ]

    setGhosts(W_, ghosts[l], 0.0)

    return W_

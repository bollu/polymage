import sys
from polymage_common import setGhosts

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def defect(U_, F_, l, name, impipeDict):
    if U_ == None:
        return F_

    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    invhh = impipeDict['invhh']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    innerBox = interior[l]['innerBox']

    W_ = Function(([z, y, x], \
                   [extent[l], extent[l], extent[l]]), \
                   Double, \
                   str(name))
    W_.defn = [ Case(innerBox,
                      F_(z  , y  , x  )
                   - (U_(z  , y  , x  ) * 6.0 \
                   -  U_(z-1, y  , x  )       \
                   -  U_(z+1, y  , x  )       \
                   -  U_(z  , y-1, x  )       \
                   -  U_(z  , y+1, x  )       \
                   -  U_(z  , y  , x-1)       \
                   -  U_(z  , y  , x+1)       \
                     ) * invhh[l]) ]

    setGhosts(W_, ghosts[l], 0.0)

    return W_

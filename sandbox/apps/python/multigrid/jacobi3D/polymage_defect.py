import sys
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def defect(U_, F_, l, name, pipeData):
    if U_ == None:
        return F_

    z = pipeData['z']
    y = pipeData['y']
    x = pipeData['x']

    invhh = pipeData['invhh']

    extent = pipeData['extent']
    interior = pipeData['interior']
    ghosts = pipeData['ghosts']

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

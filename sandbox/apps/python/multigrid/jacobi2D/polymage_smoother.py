import sys
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')
sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def wJacobi(U_, F_, l, name, impipeDict, dataDict):
    y = impipeDict['y']
    x = impipeDict['x']

    L = dataDict['L']

    invhh = impipeDict['invhh']

    jacobi_c = impipeDict['jacobi_c']
    c = jacobi_c[l]

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    innerBox = interior[l]['innerBox']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))
    if U_ != None:
        W_.defn = [ Case(innerBox,
                         U_(y, x) - c * (( \
                         U_(y  , x  ) * 4.0 \
                       - U_(y-1, x  )       \
                       - U_(y+1, x  )       \
                       - U_(y  , x-1)       \
                       - U_(y  , x+1)       \
                       ) * invhh[l]         \
                       - F_(y, x))) ]
    else:
        W_.defn = [ Case(innerBox, c * F_(y, x)) ]

    if l == L:
        setGhosts(W_, ghosts[l], U_(y, x))
    else:
        setGhosts(W_, ghosts[l], 0.0)

    return W_

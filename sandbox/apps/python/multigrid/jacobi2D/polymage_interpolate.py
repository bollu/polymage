from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def interpolate(U_, correction, l, name, impipeData):
    if U_ == None:
        return correction

    y = impipeData['y']
    x = impipeData['x']

    if correction == None:
        correct = 0.0
    else:
        correct = correction(y, x)

    extent = impipeData['extent']
    interior = impipeData['interior']
    ghosts = impipeData['ghosts']

    innerBox = interior[l]['innerBox']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))

    even_y = Condition(y%2, '==', 0)
    even_x = Condition(x%2, '==', 0)

    expr_00 =  U_(y/2  , x/2  )
    expr_01 = (U_(y/2  , x/2  )
             + U_(y/2  , x/2+1)) * 0.5
    expr_10 = (U_(y/2  , x/2  )
             + U_(y/2+1, x/2  )) * 0.5
    expr_11 = (U_(y/2  , x/2  )
             + U_(y/2+1, x/2  )
             + U_(y/2  , x/2+1)
             + U_(y/2+1, x/2+1)) * 0.25

    W_.defn = [ Case(innerBox,
                    correct + \
                      Select(even_y,
                        Select(even_x,
                          expr_00,
                          expr_01),
                        Select(even_x,
                          expr_10,
                          expr_11))) ]

    setGhosts(W_, ghosts[l], correct)

    return W_

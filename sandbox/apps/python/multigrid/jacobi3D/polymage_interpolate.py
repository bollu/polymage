from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def interpolate(U_, correction, l, name, pipeData):
    if U_ == None:
        return correction

    z = pipeData['z']
    y = pipeData['y']
    x = pipeData['x']

    if correction == None:
        correct = 0.0
    else:
        correct = correction(z, y, x)

    extent = pipeData['extent']
    interior = pipeData['interior']
    ghosts = pipeData['ghosts']

    innerBox = interior[l]['innerBox']

    W_ = Function(([z, y, x], \
                   [extent[l], extent[l], extent[l]]), \
                   Double, \
                   str(name))

    even_z = Condition(z%2, '==', 0)
    even_y = Condition(y%2, '==', 0)
    even_x = Condition(x%2, '==', 0)

    expr_000 =  U_(z/2  , y/2  , x/2  )
    expr_001 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2  , x/2+1)) / 2.0
    expr_010 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )) / 2.0
    expr_100 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2+1, y/2  , x/2  )) / 2.0
    expr_011 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )
              + U_(z/2  , y/2  , x/2+1)
              + U_(z/2  , y/2+1, x/2+1)) / 4.0
    expr_101 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2+1, y/2  , x/2  )
              + U_(z/2  , y/2  , x/2+1)
              + U_(z/2+1, y/2  , x/2+1)) / 4.0
    expr_110 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )
              + U_(z/2+1, y/2  , x/2  )
              + U_(z/2+1, y/2+1, x/2  )) / 4.0
    expr_111 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )
              + U_(z/2  , y/2  , x/2+1)
              + U_(z/2  , y/2+1, x/2+1)
              + U_(z/2+1, y/2  , x/2  )
              + U_(z/2+1, y/2+1, x/2  )
              + U_(z/2+1, y/2  , x/2+1)
              + U_(z/2+1, y/2+1, x/2+1)) / 8.0

    W_.defn = [ Case(innerBox,
                    correct + \
                    Select(even_z,
                      Select(even_y,
                        Select(even_x,
                          expr_000,
                          expr_001),
                        Select(even_x,
                          expr_010,
                          expr_011)),
                      Select(even_y,
                        Select(even_x,
                          expr_100,
                          expr_101),
                        Select(even_x,
                          expr_110,
                          expr_111)))) ]

    setGhosts(W_, ghosts[l], correct)

    return W_

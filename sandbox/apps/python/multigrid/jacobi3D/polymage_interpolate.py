from __init__ import *

import sys
from fractions  import Fraction
from polymage_common import set_ghosts

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def interpolate(U_, correction, l, name, pipe_data):
    if U_ == None:
        return correction

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    if correction == None:
        correct = 0.0
    else:
        correct = correction(z, y, x)

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l]['inner_box']

    W_ = Function(([z, y, x], [extent[l], extent[l], extent[l]]),
                  Double, str(name))

    even_z = Condition(z%2, '==', 0)
    even_y = Condition(y%2, '==', 0)
    even_x = Condition(x%2, '==', 0)

    expr_000 =  U_(z/2  , y/2  , x/2  )
    expr_001 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2  , x/2+1)) * 0.5
    expr_010 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )) * 0.5
    expr_100 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2+1, y/2  , x/2  )) * 0.5
    expr_011 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )
              + U_(z/2  , y/2  , x/2+1)
              + U_(z/2  , y/2+1, x/2+1)) * 0.25
    expr_101 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2+1, y/2  , x/2  )
              + U_(z/2  , y/2  , x/2+1)
              + U_(z/2+1, y/2  , x/2+1)) * 0.25
    expr_110 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )
              + U_(z/2+1, y/2  , x/2  )
              + U_(z/2+1, y/2+1, x/2  )) * 0.25
    expr_111 = (U_(z/2  , y/2  , x/2  )
              + U_(z/2  , y/2+1, x/2  )
              + U_(z/2  , y/2  , x/2+1)
              + U_(z/2  , y/2+1, x/2+1)
              + U_(z/2+1, y/2  , x/2  )
              + U_(z/2+1, y/2+1, x/2  )
              + U_(z/2+1, y/2  , x/2+1)
              + U_(z/2+1, y/2+1, x/2+1)) * 0.125

    W_.defn = [ Case(inner_box,
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

    set_ghosts(W_, ghosts[l], correct)

    return W_

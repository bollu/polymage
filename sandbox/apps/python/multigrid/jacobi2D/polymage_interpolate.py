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

    y = pipe_data['y']
    x = pipe_data['x']

    if correction == None:
        correct = 0.0
    else:
        correct = correction(y, x)

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l]['inner_box']

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

    W_.defn = [ Case(inner_box,
                     correct + \
                       Select(even_y,
                         Select(even_x,
                           expr_00,
                           expr_01),
                         Select(even_x,
                           expr_10,
                           expr_11))) ]

    set_ghosts(W_, ghosts[l], correct)

    return W_

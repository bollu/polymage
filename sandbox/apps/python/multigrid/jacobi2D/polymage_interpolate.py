from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction
from polymage_common import set_ghosts

sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def interpolate(U_, correction, l, name, impipe_data):
    if U_ == None:
        return correction

    y = impipe_data['y']
    x = impipe_data['x']

    if correction == None:
        correct = 0.0
    else:
        correct = correction(y, x)

    extent = impipe_data['extent']
    interior = impipe_data['interior']
    ghosts = impipe_data['ghosts']

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

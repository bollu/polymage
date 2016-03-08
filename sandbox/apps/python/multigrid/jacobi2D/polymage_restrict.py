from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction
from polymage_common import set_ghosts

sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def restrict(U_, l, name, pipe_data):
    y = pipe_data['y']
    x = pipe_data['x']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l-1]['inner_box']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))
    W_.defn = [ Case(inner_box,
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

    set_ghosts(W_, ghosts[l-1], 0.0)

    return W_

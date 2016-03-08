import sys
from polymage_common import set_ghosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def defect(U_, F_, l, name, pipe_data):
    if U_ == None:
        return F_

    y = pipe_data['y']
    x = pipe_data['x']

    invhh = pipe_data['invhh']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l]['inner_box']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))
    W_.defn = [ Case(inner_box,
                     F_(y, x) \
                   - (U_(y  , x  ) * 4.0 \
                   -  U_(y-1, x  )       \
                   -  U_(y+1, x  )       \
                   -  U_(y  , x-1)       \
                   -  U_(y  , x+1)       \
                     ) * invhh[l]) ]

    set_ghosts(W_, ghosts[l], 0.0)

    return W_

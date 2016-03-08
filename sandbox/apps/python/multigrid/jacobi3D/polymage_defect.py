import sys
from polymage_common import set_ghosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def defect(U_, F_, l, name, pipe_data):
    if U_ == None:
        return F_

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    invhh = pipe_data['invhh']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l]['inner_box']

    W_ = Function(([z, y, x], [extent[l], extent[l], extent[l]]),
                   Double, str(name))
    W_.defn = [ Case(inner_box,
                      F_(z  , y  , x  ) \
                   - (U_(z  , y  , x  ) * 6.0 \
                   -  U_(z-1, y  , x  )       \
                   -  U_(z+1, y  , x  )       \
                   -  U_(z  , y-1, x  )       \
                   -  U_(z  , y+1, x  )       \
                   -  U_(z  , y  , x-1)       \
                   -  U_(z  , y  , x+1)       \
                     ) * invhh[l]) ]

    set_ghosts(W_, ghosts[l], 0.0)

    return W_

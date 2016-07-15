from __init__ import *
import sys
from polymage_common import set_ghosts

sys.path.insert(0, ROOT)

from compiler   import *
from constructs import *

def w_jacobi(U_, F_, l, name, app_data):
    pipe_data = app_data['pipe_data']

    y = pipe_data['y']
    x = pipe_data['x']

    L = app_data['L']

    invhh = pipe_data['invhh']

    jacobi_c = pipe_data['jacobi_c']
    c = jacobi_c[l]

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l]['inner_box']

    W_ = Function(([y, x], [extent[l], extent[l]]), Double, str(name))
    if U_ != None:
        W_.defn = [ Case(inner_box,
                         U_(y, x) - c * (( \
                         U_(y  , x  ) * 4.0 \
                       - U_(y-1, x  )       \
                       - U_(y+1, x  )       \
                       - U_(y  , x-1)       \
                       - U_(y  , x+1)       \
                       ) * invhh[l]         \
                       - F_(y, x))) ]
    else:
        W_.defn = [ Case(inner_box, c * F_(y, x)) ]

    if l == L:
        set_ghosts(W_, ghosts[l], U_(y, x))
    else:
        set_ghosts(W_, ghosts[l], 0.0)

    return W_

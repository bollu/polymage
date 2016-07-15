from __init__ import *
import sys
from fractions  import Fraction
from polymage_common import set_ghosts

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def restrict(U_, l, name, pipe_data):
    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    inner_box = interior[l-1]['inner_box']

    W_ = Function(([z, y, x], [extent[l], extent[l], extent[l]]),
                  Double, str(name))

    W_.defn = [ Case(inner_box,
# corners
                     (U_(2*z-1, 2*y-1, 2*x-1)             \
                    + U_(2*z-1, 2*y-1, 2*x+1)             \
                    + U_(2*z-1, 2*y+1, 2*x-1)             \
                    + U_(2*z-1, 2*y+1, 2*x+1)             \
                    + U_(2*z+1, 2*y-1, 2*x-1)             \
                    + U_(2*z+1, 2*y-1, 2*x+1)             \
                    + U_(2*z+1, 2*y+1, 2*x-1)             \
                    + U_(2*z+1, 2*y+1, 2*x+1)) * 1.0/64.0 \
# edge centers
                    +(U_(2*z-1, 2*y-1, 2*x  )             \
                    + U_(2*z-1, 2*y  , 2*x-1)             \
                    + U_(2*z-1, 2*y  , 2*x+1)             \
                    + U_(2*z-1, 2*y+1, 2*x  )             \
                    + U_(2*z  , 2*y-1, 2*x-1)             \
                    + U_(2*z  , 2*y-1, 2*x+1)             \
                    + U_(2*z  , 2*y+1, 2*x-1)             \
                    + U_(2*z  , 2*y+1, 2*x+1)             \
                    + U_(2*z+1, 2*y-1, 2*x  )             \
                    + U_(2*z+1, 2*y  , 2*x-1)             \
                    + U_(2*z+1, 2*y  , 2*x+1)             \
                    + U_(2*z+1, 2*y+1, 2*x  )) * 1.0/32.0 \
# face centers
                    +(U_(2*z-1, 2*y  , 2*x  )             \
                    + U_(2*z+1, 2*y  , 2*x  )             \
                    + U_(2*z  , 2*y-1, 2*x  )             \
                    + U_(2*z  , 2*y+1, 2*x  )             \
                    + U_(2*z  , 2*y  , 2*x-1)             \
                    + U_(2*z  , 2*y  , 2*x+1)) * 1.0/16.0 \
# cube center
                    + U_(2*z  , 2*y  , 2*x  )  * 1.0/8.0) ]

    set_ghosts(W_, ghosts[l-1], 0.0)

    return W_

from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction
from polymage_common import setGhosts

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def restrict(U_, l, name, pipeData):
    z = pipeData['z']
    y = pipeData['y']
    x = pipeData['x']

    extent = pipeData['extent']
    interior = pipeData['interior']
    ghosts = pipeData['ghosts']

    innerBox = interior[l-1]['innerBox']

    W_ = Function(([z, y, x], \
                   [extent[l], extent[l], extent[l]]), \
                   Double, \
                   str(name))

    W_.defn = [ Case(innerBox,
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

    setGhosts(W_, ghosts[l-1], 0.0)

    return W_

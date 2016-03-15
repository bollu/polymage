from __future__ import absolute_import, division, print_function

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def interpolate(pipe_data):
    L = 10 # number of pyramid levels

    # Params
    R = Parameter(Int, "R") # image rows
    C = Parameter(Int, "C") # image cols

    pipe_data['R'] = R
    pipe_data['C'] = C

    # Vars
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")

    # Input Image
    img = Image(Float, "img", [4, R+2, C+2])

    # Intervals
    rgba = Interval(Int, 0, 3) # colour with alpha channel
    rgb = Interval(Int, 0, 2) # colour channels
    row = Interval(Int, 0, R+1)
    col = Interval(Int, 0, C+1)
    row2 = Interval(Int, 0, R-1)
    col2 = Interval(Int, 0, C-1)

    #####################################################################################

    # DOWNSAMPLE
    def pyrDown(f, l, name):
        # decrement factor
        decFactor = 1 << l
        # original factor
        orgFactor = 1 << (l-1)

        # domain (downsampled)
        decRowr = Interval(Int, 0, (R/decFactor)+1)
        decColr = Interval(Int, 0, (C/decFactor)+1)
        # domain (original)
        colr = Interval(Int, 0, (C/orgFactor)+1)

        # downsample in 'x' dimension (using [1 2 1]' filter)

        # body case
        condx = Condition(x, '>=', 1) & \
                Condition(x, '<=', (R/decFactor)) & \
                Condition(y, '>=', 1) & \
                Condition(y, '<=', (C/orgFactor))
        # boundary cases
        condxLeft   = Condition(y, '<=', 0)
        condxRight  = Condition(y, '>=', (C/orgFactor)+1)
        condxTop    = Condition(x, '<=', 0)
        condxBottom = Condition(x, '>=', (R/decFactor)+1)

        downx = Function(([c, x, y], [rgba, decRowr, colr]), Float, "Dx_" + str(l) + "_" + name)
        downx.defn = [Case(condx, (f(c, 2*x-2, y) + \
                                               f(c, 2*x-1, y) * 2.0 + \
                                               f(c, 2*x  , y) \
                                              ) * 0.25),
                     Case(condxLeft, 0),
                     Case(condxRight, 0),
                     Case(condxBottom, 0),
                     Case(condxTop, 0)]
    
        # Please refer to 'local_laplacian' app to check how to set
        # boundary conditions if necessary

        # downsample in 'y' dimension (using [1 2 1] filter)

        condy = Condition(x, '>=', 1) & \
                Condition(x, '<=', (R/decFactor)) & \
                Condition(y, '>=', 1) & \
                Condition(y, '<=', (C/decFactor))
        condyLeft   = Condition(y, '<=', 0)
        condyRight  = Condition(y, '>=', (C/decFactor)+1)
        condyTop    = Condition(x, '<=', 0)
        condyBottom = Condition(x, '>=', (R/decFactor)+1)

        downy = Function(([c, x, y], [rgba, decRowr, decColr]), Float, "D_" + str(l) + "_" + name)
        downy.defn = [Case(condy, (downx(c, x, 2*y-2) + \
                                               downx(c, x, 2*y-1) * 2.0 + \
                                               downx(c, x, 2*y  )
                                              ) * 0.25),
                     Case(condyLeft, 0),
                     Case(condyRight, 0),
                     Case(condyBottom, 0),
                     Case(condyTop, 0)]


        return downy

    # UPSAMPLE & INTERPOLATE
    def pyrUp_Xpolate(f, d, l, name):
        decFactor = 1 << l+1
        orgFactor = 1 << l

        # domain (original)
        decRowr = Interval(Int, 0, (R/decFactor)+1)
        decColr = Interval(Int, 0, (C/decFactor)+1)

        # domain (upsampled)
        colr = Interval(Int, 0, (C/orgFactor)+1)
        rowr = Interval(Int, 0, (R/orgFactor)+1)

        # upsample in 'x' dimension
        condx = Condition(x, '>=', 1) & \
                Condition(x, '<=', (R/orgFactor)) & \
                Condition(y, '>=', 1) & \
                Condition(y, '<=', (C/decFactor))
        condxLeft   = Condition(y, '<=', 0)
        condxRight  = Condition(y, '>=', (C/decFactor)+1)
        condxTop    = Condition(x, '<=', 0)
        condxBottom = Condition(x, '>=', (R/orgFactor)+1)

        upx = Function(([c, x, y], [rgba, rowr, decColr]), Float, "Ux_" + name)
        upx.defn = [Case(condx, (f(c, (x+1)/2, y) + \
                                             f(c, (x+2)/2, y) \
                                            ) / 2.0),
                   Case(condxLeft, 0),
                   Case(condxRight, 0),
                   Case(condxBottom, 0),
                   Case(condxTop, 0)]

        # upsample in 'y' dimension and interpolate
        condy = Condition(x, '>=', 1) & \
                Condition(x, '<=', (R/orgFactor)) & \
                Condition(y, '>=', 1) & \
                Condition(y, '<=', (C/orgFactor))
        condyLeft   = Condition(y, '<=', 0)
        condyRight  = Condition(y, '>=', (C/orgFactor)+1)
        condyTop    = Condition(x, '<=', 0)
        condyBottom = Condition(x, '>=', (R/orgFactor)+1)

        upy = (upx(c, x, (y+1)/2) + \
               upx(c, x, (y+2)/2)) / 2.0

        interpolate = Function(([c, x, y], [rgba, rowr, colr]), Float, name)
        interpolate.defn = [Case( condy, d(c, x, y) + \
                                                     (1.0 - d(3, x, y)) * upy ),
                           Case(condyLeft, 0),
                           Case(condyRight, 0),
                           Case(condyBottom, 0),
                           Case(condyTop, 0)]

        return interpolate


    #####################################################################################
    # MULTISCALE INTERPOLATION

    # 1.
    downsampled = {}
    # level : 0
    downsampled[0] = Function(([c, x, y], [rgba, row, col]), Float, "downsampled_L0")
    downsampled[0].defn = [img(c, x, y) * img(3, x, y)]
    # level : everything else
    for l in range(1, L):
        downsampled[l] = pyrDown(downsampled[l-1], l, "downsampled")


    # 2.
    interpolated = {}
    # level : L-1
    interpolated[L-1] = downsampled[L-1]
    # level : everything else
    for l in range(L-2, -1, -1):
        interpolated[l] = pyrUp_Xpolate(interpolated[l+1], downsampled[l], l, "interpolated_L" + str(l))


    # 3. Normalize
    normalized = Function(([c, x, y], [rgb, row2, col2]), Float, "interpolate")
    normalized.defn = [interpolated[0](c, x+1, y+1) / interpolated[0](3, x+1, y+1)]

    return normalized

    #####################################################################################
# END

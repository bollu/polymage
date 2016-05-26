from __init__ import *

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def pyramid_blending(pipe_data):

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")

    pipe_data['R'] = R
    pipe_data['C'] = C

    L = 4

    # Input Images
    img1 = Image(Float, "img1", [3, R, C])
    img2 = Image(Float, "img2", [3, R, C])
    mask = Image(Float, "mask", [R, C])

    # Intervals
    cr = Interval(Int, 0, 2)

    # Gaussian pyramid construction
    Gauss_img1 = {}
    Gauss_img2 = {}
    Gauss_mask = {}
    Gauss_img1[0] = img1
    Gauss_img2[0] = img2
    Gauss_mask[0] = mask

    def pyr_down(f, l, color, name):
        dec_factor = 1 << l
        org_factor = 1 << (l-1)
        dec_rowr = Interval(Int, l, (R/dec_factor)-2)
        dec_colr = Interval(Int, l, (C/dec_factor)-2)

        colr = Interval(Int, l-1, (C/org_factor)-2)

        if color:
            downx = Function(([c, x, y], [cr, dec_rowr, colr]),
                             Float, "Dx_" + str(l) + "_" + name)
            downx.defn = [ (1 * f(c, 2*x-2, y) + \
                            4 * f(c, 2*x-1, y) + \
                            6 * f(c, 2*x  , y) + \
                            4 * f(c, 2*x+1, y) + \
                            1 * f(c, 2*x+2, y)) * 0.0625 ]
        else:
            downx = Function(([x, y], [dec_rowr, colr]),
                             Float, "Dx_" + str(l) + "_" + name)
            downx.defn = [ (1 * f(2*x-2, y) + \
                            4 * f(2*x-1, y) + \
                            6 * f(2*x  , y) + \
                            4 * f(2*x+1, y) + \
                            1 * f(2*x+2, y)) * 0.0625 ]

        if color:
            downy = Function(([c, x, y], [cr, dec_rowr, dec_colr]),
                             Float, "Dy_" + str(l) + "_" + name)
            downy.defn = [ (1 * downx(c, x, 2*y-2) + \
                            4 * downx(c, x, 2*y-1) + \
                            6 * downx(c, x, 2*y  ) + \
                            4 * downx(c, x, 2*y+1) + \
                            1 * downx(c, x, 2*y+2)) * 0.0625 ]
        else:
            downy = Function(([x, y], [dec_rowr, dec_colr]),
                             Float, "Dy_" + str(l) + "_" + name)
            downy.defn = [ (1 * downx(x, 2*y-2) + \
                            4 * downx(x, 2*y-1) + \
                            6 * downx(x, 2*y  ) + \
                            4 * downx(x, 2*y+1) + \
                            1 * downx(x, 2*y+2)) * 0.0625 ]

        return downy

    for l in range(1, L):
        Gauss_img1[l] = pyr_down(Gauss_img1[l-1], l, True, "img1")
        Gauss_img2[l] = pyr_down(Gauss_img2[l-1], l, True, "img2")
        Gauss_mask[l] = pyr_down(Gauss_mask[l-1], l, False, "mask")

    def lapl(fd, f, l, name):

        # 'l' is the output level

        # input interval
        dec_colr = Interval(Int, (1<<(L-l))-1, C/(1<<(l+1)) - (1<<(L-l)) + 2)

        # output interval
        colr = Interval(Int, (1<<(L-l+1))-1, C/(1<<l) - (1<<(L-l+1)) + 2)
        rowr = Interval(Int, (1<<(L-l+1))-1, R/(1<<l) - (1<<(L-l+1)) + 2)

        upx = Function(([c, x, y], [cr, rowr, dec_colr]),
                       Float, "Ux_" + str(l) + "_" + name)

        evenxexpr = (1 * fd(c, (x/2)-1, y) + \
                     6 * fd(c, (x/2)  , y) + \
                     1 * fd(c, (x/2)+1, y)) * 0.125
        oddxexpr = (4 * fd(c, (x-1)/2, y) + \
                    4 * fd(c, (x+1)/2, y)) * 0.125

        upxexpr = Select(Condition(x%2, "==", 0), evenxexpr, oddxexpr)
        upx.defn = [ upxexpr ]

        upy = Function(([c, x, y], [cr, rowr, colr]),
                       Float, "Uy_" + str(l) + "_" + name)

        evenyexpr =  f(c, x, y) - ((1 * upx(c, x, (y/2)-1) + \
                                    6 * upx(c, x, (y/2)  ) + \
                                    1 * upx(c, x, (y/2)+1)) * 0.125)
        oddyexpr  =  f(c, x, y) - ((4 * upx(c, x, (y-1)/2) + \
                                    4 * upx(c, x, (y+1)/2)) * 0.125)

        upyexpr = Select(Condition(y%2, "==", 0), evenyexpr, oddyexpr)
        upy.defn = [ upyexpr ]

        return upy

    # Laplacian pyramid construction
    Lapl_img1 = {}
    Lapl_img2 = {}
    for l in range(0, L-1):
        Lapl_img1[l] = lapl(Gauss_img1[l+1], Gauss_img1[l], l, "img1")
        Lapl_img2[l] = lapl(Gauss_img2[l+1], Gauss_img2[l], l, "img2")

    Lapl_img1[L-1] = Gauss_img1[L-1]
    Lapl_img2[L-1] = Gauss_img2[L-1]

    def blend(lap1, lap2, m, l):
        colr = Interval(Int, (1<<(L-l+1))-1, C/(1<<l) - (1<<(L-l+1)) + 2)
        rowr = Interval(Int, (1<<(L-l+1))-1, R/(1<<l) - (1<<(L-l+1)) + 2)

        reslap = Function(([c, x, y], [cr, rowr, colr]),
                          Float, "Res_" + str(l))
        reslap.defn = [ lap1(c, x, y) * m(x, y) + lap2(c, x, y) * (1-m(x, y)) ]
        return reslap

    # Pyramid blending
    res_lapl = {}
    for l in range(0, L):
        res_lapl[l] = blend(Lapl_img1[l], Lapl_img2[l], Gauss_mask[l], l)
 
    def collapse(fd, f, l, name):
        # 'l' is the output level

        # input interval
        dec_colr = Interval(Int, (1<<(L-l))-1, C/(1<<(l+1)) - (1<<(L-l)) + 2)

        # output interval
        colr = Interval(Int, (1<<(L-l+1))-1, C/(1<<l) - (1<<(L-l+1)) + 2)
        rowr = Interval(Int, (1<<(L-l+1))-1, R/(1<<l) - (1<<(L-l+1)) + 2)

        upx = Function(([c, x, y], [cr, rowr, dec_colr]),
                       Float, str(name)+'_x')

        evenxexpr =  (1 * fd(c, (x/2)-1, y) + \
                      6 * fd(c, (x/2)  , y) + \
                      1 * fd(c, (x/2)+1, y)) * 0.125
        oddxexpr  =  (4 * fd(c, (x-1)/2, y) + \
                      4 * fd(c, (x+1)/2, y)) * 0.125

        upxexpr = Select(Condition(x%2, "==", 0), evenxexpr, oddxexpr)
        upx.defn = [ upxexpr ]

        evenyexpr =  f(c, x, y) + ((1 * upx(c, x, (y/2)-1) + \
                                    6 * upx(c, x, (y/2)  ) + \
                                    1 * upx(c, x, (y/2)+1)) * 0.125)
        oddyexpr  =  f(c, x, y) + ((4 * upx(c, x, (y-1)/2) + \
                                    4 * upx(c, x, (y+1)/2)) * 0.125)

        upyexpr = Select(Condition(y%2, "==", 0), evenyexpr, oddyexpr)

        upy = Function(([c, x, y], [cr, rowr, colr]),
                       Float, str(name))
        upy.defn = [ upyexpr ]

        return upy

    # Collapsing blended pyramid
    out_pyr = {}
    out_pyr[L-1] = res_lapl[L-1]
    for l in range(L-2, 0, -1):
        out_pyr[l] = collapse(out_pyr[l+1], res_lapl[l], l, "Col_"+str(l))

    out_pyr[0] = collapse(out_pyr[1], res_lapl[0], 0, "blend")

    return out_pyr[0]

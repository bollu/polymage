from __future__ import absolute_import, division, print_function

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

# PolyMage Specification
# ======================

def pyramid_blending(pipe_data):
    L = 4 # number of pyramid levels

    # Params
    R = Parameter(Int, "R") # image rows
    C = Parameter(Int, "C") # image cols

    # Vars
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")

    # Input Images
    img1 = Image(Float, "img1", [3, R+4, C+4])
    img2 = Image(Float, "img2", [3, R+4, C+4])
    mask = Image(Float, "mask", [R+4, C+4])

    # register in the dictionary
    pipe_data['R'] = R
    pipe_data['C'] = C

    # Intervals
    cr = Interval(Int, 0, 2) # colour channels

    #####################################################################################

    # DOWNSAMPLE
    def pyrDown(f, l, colour, name):
        # decrement factor
        decFactor = 1 << l
        # original factor
        orgFactor = 1 << (l-1)

        # domain (downsampled)
        decRowr = Interval(Int, 0, (R/decFactor) + 3)
        decColr = Interval(Int, 0, (C/decFactor) + 3)
        # domain (original)
        colr = Interval(Int, 0, (C/orgFactor)+3)

        # downsample in 'x' dimension (using [1 4 6 4 1]' filter)
        # body case
        condx = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/decFactor) + 1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/orgFactor) + 1)
     
        # if the image is coloured
        if colour:
            downx = Function(([c, x, y], [cr, decRowr, colr]), Float, "Dx_" + str(l) + "_" + name)
            downx.defn = [Case(condx, (f(c, 2*(x-1)-2, y) + \
                                                   f(c, 2*(x-1)-1, y) * 4.0 + \
                                                   f(c, 2*(x-1)  , y) * 6.0 + \
                                                   f(c, 2*(x-1)+1, y) * 4.0 + \
                                                   f(c, 2*(x-1)+2, y) \
                                                  ) * 0.0625)]
        # if grayscale image
        else:            
            downx = Function(([x, y], [decRowr, colr]), Float, "Dx_" + str(l) + "_" + name)
            downx.defn = [Case(condx, (f(2*(x-1)-2, y) + \
                                                   f(2*(x-1)-1, y) * 4.0 + \
                                                   f(2*(x-1)  , y) * 6.0 + \
                                                   f(2*(x-1)+1, y) * 4.0 + \
                                                   f(2*(x-1)+2, y) \
                                                  ) * 0.0625)]

        #fi

        # Please refer to 'local_laplacian' app to check how to set
        # boundary conditions if necessary

        # downsample in 'y' dimension (using [1 4 6 4 1] filter)

        condy = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/decFactor) + 1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/decFactor) + 1)

        if colour:
            downy = Function(([c, x, y], [cr, decRowr, decColr]), Float, "Dy_" + str(l) + "_" + name)
            downy.defn = [Case(condy, (downx(c, x, 2*(y-1)-2) + \
                                                   downx(c, x, 2*(y-1)-1) * 4.0 + \
                                                   downx(c, x, 2*(y-1)  ) * 6.0 + \
                                                   downx(c, x, 2*(y-1)+1) * 4.0 + \
                                                   downx(c, x, 2*(y-1)+2) \
                                                  ) * 0.0625)]
        else:
            downy = Function(([x, y], [decRowr, decColr]), Float, "Dy_" + str(l) + "_" + name)
            downy.defn = [Case(condy, (downx(x, 2*(y-1)-2) + \
                                                   downx(x, 2*(y-1)-1) * 4.0 + \
                                                   downx(x, 2*(y-1)  ) * 6.0 + \
                                                   downx(x, 2*(y-1)+1) * 4.0 + \
                                                   downx(x, 2*(y-1)+2) \
                                                  ) * 0.0625)]

        #fi

        return downy


    # UPSAMPLE
    def pyrUp(f, l, name):
        decFactor = 1 << l+1
        orgFactor = 1 << l

        # domain (original)
        decColr = Interval(Int, 0, (C/decFactor) + 3)

        # domain (upsampled)
        colr = Interval(Int, 0, (C/orgFactor)+3)
        rowr = Interval(Int, 0, (R/orgFactor)+3)

        # upsample in 'x' dimension
        condx = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/orgFactor) + 1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/decFactor) + 1)

        evenxexpr =  (f(c, x/2  , y) + \
                      f(c, x/2+1, y) * 6.0 + \
                      f(c, x/2+2, y) \
                     ) * 0.125 
        oddxexpr  =  (f(c, (x-1)/2+1, y) * 4.0 + \
                      f(c, (x+1)/2+1, y) * 4.0 \
                     ) * 0.125
        upxexpr = Select(Condition(x%2, "==", 0), \
                                evenxexpr, \
                                oddxexpr)

        upx = Function(([c, x, y], [cr, rowr, decColr]), Float, "Ux_" + str(l) + "_" + name)
        upx.defn = [Case(condx, upxexpr)]

        # upsample in 'y' dimension
        condy = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/orgFactor) + 1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/orgFactor) + 1)

        evenyexpr =  (upx(c, x, y/2  ) + \
                      upx(c, x, y/2+1) * 6.0 + \
                      upx(c, x, y/2+2) \
                     ) * 0.125
        oddyexpr  =  (upx(c, x, (y-1)/2+1) * 4.0 + \
                      upx(c, x, (y+1)/2+1) * 4.0 \
                     ) * 0.125
        upyexpr = Select(Condition(y%2, "==", 0), \
                                evenyexpr, \
                                oddyexpr)

        upy = Function(([c, x, y], [cr, rowr, colr]), Float, "Uy_" + str(l) + "_" + name)
        upy.defn = [Case(condy, upyexpr)]
        return upy


    def laplace(f1, f2, l, name):
        colr = Interval(Int, 0, (C/(1<<l))+3)
        rowr = Interval(Int, 0, (R/(1<<l))+3)

        lapl = Function(([c, x, y], [cr, rowr, colr]), Float, name+"_"+str(l))
        lapl.defn = [f2(c, x, y) - f1(c, x, y)]

        return lapl

    def blend(lap1, lap2, m, l):
        rowr = Interval(Int, 0, (R/(1<<l)) + 3)
        colr = Interval(Int, 0, (C/(1<<l)) + 3)

        reslap = Function(([c, x, y], [cr, rowr, colr]), Float, "Res_" + str(l))
        reslap.defn = [lap1(c, x, y) * (    m(x, y)) + \
                            lap2(c, x, y) * (1 - m(x, y))]

        return reslap
 
    def collapse(f1, f2, l, name):
        colr = Interval(Int, 0, (C/(1<<l))+3)
        rowr = Interval(Int, 0, (R/(1<<l))+3)

        coll = Function(([c, x, y], [cr, rowr, colr]), Float, name+"_"+str(l))
        coll.defn = [f1(c, x, y) + f2(c, x, y)]

        return coll


    #####################################################################################
    # PYRAMID BLENDING
 
    # 1. Gaussian pyramid construction
    GaussImg1 = {}
    GaussImg2 = {}
    GaussMask = {}
    GaussImg1[0] = img1
    GaussImg2[0] = img2
    GaussMask[0] = mask
    for l in range(1, L):
        GaussImg1[l] = pyrDown(GaussImg1[l-1], l, 1, "img1")
        GaussImg2[l] = pyrDown(GaussImg2[l-1], l, 1, "img2")
        GaussMask[l] = pyrDown(GaussMask[l-1], l, 0, "mask")


    # 2. Laplacian pyramid construction
    LaplImg1 = {}
    LaplImg2 = {}
    for l in range(0, L-1):
        LaplImg1[l] = laplace(pyrUp(GaussImg1[l+1], l, "img1"), GaussImg1[l], l, "lapl1")
        LaplImg2[l] = laplace(pyrUp(GaussImg2[l+1], l, "img2"), GaussImg2[l], l, "lapl2")
    LaplImg1[L-1] = GaussImg1[L-1]
    LaplImg2[L-1] = GaussImg2[L-1]


    # 3. Blend the images using mask
    resLapl = {}
    for l in range(0, L):
        resLapl[l] = blend(LaplImg1[l], LaplImg2[l], GaussMask[l], l)


    # 4. Collapsing pyramid of blended image samples
    outPyr = {}
    outPyr[L-1] = resLapl[L-1]
    for l in range(L-2, 0, -1):
        outPyr[l] = collapse(pyrUp(outPyr[l+1], l, "coll"), resLapl[l], l, "coll")
    outPyr[0] = collapse(pyrUp(outPyr[1], 0, "blend"), resLapl[0], 0, "blend")

    return outPyr[0]

    #####################################################################################
# END

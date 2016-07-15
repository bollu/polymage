from __init__ import *

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

# PolyMage Specification
# ======================

def unsharp_mask(pipe_data):

    # Params
    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    threshold = Parameter(Float, "threshold")
    weight = Parameter(Float, "weight")
 
    pipe_data['R'] = R
    pipe_data['C'] = C
    pipe_data['threshold'] = threshold
    pipe_data['weight'] = weight

    # Vars
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")

    # Input Image
    img = Image(Float, "img", [3, R+4, C+4])

    # Intervals
    cr = Interval(Int, 0, 2)
    xrow = Interval(Int, 2, R+1)
    xcol = Interval(Int, 0, C+3)
    yrow = Interval(Int, 2, R+1)
    ycol = Interval(Int, 2, C+1)

    cond = Condition(x, '>=', 1) & Condition(x, '<=', R) & \
           Condition(y, '<=', C) & Condition(y, '>=', 1)

    #####################################################################################
    # UNSHARP MASK

    blurx = Function(([c, x, y], [cr, xrow, xcol]), Float, "blurx")
    blurx.defn = [ Case(cond,( img(c, x-2, y) \
                       + img(c, x-1, y) * 4.0 \
                       + img(c, x  , y) * 6.0 \
                       + img(c, x+1, y) * 4.0 \
                       + img(c, x+2, y) ) * 0.0625)]

    blury = Function(([c, x, y], [cr, yrow, ycol]), Float, "blury")
    blury.defn = [ Case(cond,( blurx(c, x, y-2) \
                       + blurx(c, x, y-1) * 4.0 \
                       + blurx(c, x, y  ) * 6.0 \
                       + blurx(c, x, y+1) * 4.0 \
                       + blurx(c, x, y+2) ) * 0.0625)]

    sharpen = Function(([c, x, y], [cr, yrow, ycol]), Float, "sharpen")
    sharpen.defn = [ Case(cond,img(c, x, y)   * ( 1 + weight ) \
                       - blury(c, x, y) * (     weight ))]

    masked = Function(([c, x, y], [cr, yrow, ycol]), Float, "mask")
    masked.defn = [ Case(cond,Select( Condition( \
                                         Abs(img(c, x, y) - blury(c, x, y)), \
                                         '<', \
                                         threshold), \
                                       img(c, x, y), \
                                       sharpen(c, x, y) ))]

    #####################################################################################
    return masked
# END

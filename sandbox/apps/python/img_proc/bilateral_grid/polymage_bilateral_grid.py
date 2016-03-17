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

def bilateral_grid(pipe_data):
    # Params
    rows = Parameter(Int, "rows")
    cols = Parameter(Int, "cols")

    pipe_data['R'] = rows
    pipe_data['C'] = cols

    # Vars
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    z = Variable(Int, "z")
    c = Variable(Int, "c")

    # Intervals
    rowr = Interval(Int, 0, rows-1)
    colr = Interval(Int, 0, cols-1)

    # Input Images
    img = Image(Float, "input", [rows, cols])


    #####################################################################################
    # BILATERAL GRID

    # sigma_s is hardcoded to 8 and sigma_r is hardcoded to 0.1 
    # Bilateral grid
    invSigma_r = 10
    sigma_s = 8

    gridrowr = Interval(Int, 0, rows/sigma_s + 3)
    gridcolr = Interval(Int, 0, cols/sigma_s + 3)
    gridintr = Interval(Int, 0, invSigma_r + 3)

    cr = Interval(Int, 0, 1)
    cd = Interval(Int, 0, 0)

    cond0 = Condition(c, '==', 0)
    cond1 = Condition(c, '==', 1)

    # Accumulation
    grid = Reduction(([c, x, y, z], [cr, gridrowr, gridcolr, gridintr]), ([x, y], [rowr, colr]), Float, "grid")
    intaccess = Cast(Int, (img(x, y) * invSigma_r))
    grid.defn = [ Case(cond0, Reduce(grid(0, 2 + (x/sigma_s), 2 + (y/sigma_s), 2 + intaccess), \
                                        img(x, y), Op.Sum)),
                  Case(cond1, Reduce(grid(1, 2 + (x/sigma_s), 2 + (y/sigma_s), 2 + intaccess), \
                                        1, Op.Sum)) ]
    grid.default = 0

    cond = Condition(x, '>=', 2) & \
           Condition(x, '<=', rows/sigma_s + 1) & \
           Condition(y, '>=', 2) & \
           Condition(y, '<=', cols/sigma_s + 1) & \
           Condition(z, '>=', 2) & \
           Condition(z, '<=', invSigma_r + 1)

    blurz = Function(([c, x, y, z], [cr, gridrowr, gridcolr, gridintr]), Float, "blurz")
    blurz.defn = [Case(cond, (grid(c, x, y, z-2) + \
                                          grid(c, x, y, z-1) * 4 + \
                                          grid(c, x, y, z  ) * 6 + \
                                          grid(c, x, y, z+1) * 4 + \
                                          grid(c, x, y, z+2))) ]

    blurx = Function(([c, x, y, z], [cr, gridrowr, gridcolr, gridintr]), Float, "blurx")
    blurx.defn = [ Case(cond, (blurz(c, x-2, y, z) + \
                                          blurz(c, x-1, y, z) * 4 + \
                                          blurz(c, x  , y, z) * 6 + \
                                          blurz(c, x+1, y, z) * 4 + \
                                          blurz(c, x+2, y, z))) ]

    blury = Function(([c, x, y, z], [cr, gridrowr, gridcolr, gridintr]), Float, "blury")
    blury.defn = [ Case(cond, (blurx(c, x, y-2, z) + \
                                          blurx(c, x, y-1, z) * 4 + \
                                          blurx(c, x, y  , z) * 6 + \
                                          blurx(c, x, y+1, z) * 4 + \
                                          blurx(c, x, y+2, z))) ]

    zv = img(x, y) * invSigma_r
    zi = Cast(Int, zv)
    zf = zv - zi
    xf = Cast(Float, x % sigma_s) / sigma_s
    yf = Cast(Float, y % sigma_s) / sigma_s
    xi = x/sigma_s
    yi = y/sigma_s

   
    def lerp(a, b, w):
        return a + w*(b-a)

    xintery0z0 = lerp(blury(c, 2+xi  , 2+yi  , 2+zi  ), \
                      blury(c, 2+xi+1, 2+yi  , 2+zi  ), xf)
    xintery1z0 = lerp(blury(c, 2+xi  , 2+yi+1, 2+zi  ), \
                      blury(c, 2+xi+1, 2+yi+1, 2+zi  ), xf) 
    xinteryinterz0 = lerp(xintery0z0, xintery1z0, yf)

    xintery0z1 = lerp(blury(c, 2+xi  , 2+yi  , 2+zi+1), \
                      blury(c, 2+xi+1, 2+yi  , 2+zi+1), xf)
    xintery1z1 = lerp(blury(c, 2+xi  , 2+yi+1, 2+zi+1), \
                      blury(c, 2+xi+1, 2+yi+1, 2+zi+1), xf)
    xinteryinterz1 = lerp(xintery0z1, xintery1z1, yf)

    #xinteryinterzinter = lerp(xinteryinterz0, xinteryinterz1, zf)

    
    interpolated = Function(([c, x, y], [cr, rowr, colr]), Float, "interpolated")
    interpolated.defn = [lerp(xinteryinterz0, xinteryinterz1, zf)]


    filterexpr = Select(Condition(interpolated(1, x, y), '>', 0), \
                               interpolated(0, x, y) / interpolated(1, x, y), \
                               img(x, y)) 
    filtered = Function(([c, x, y], [cd, rowr, colr]), Float, "filtered")
    filtered.defn = [filterexpr]

    return filtered
    #####################################################################################
# END

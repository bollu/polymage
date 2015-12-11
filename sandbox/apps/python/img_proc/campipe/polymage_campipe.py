from __future__ import absolute_import, division, print_function

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def camera_pipe(pipe_data):

    # Pipeline Variables
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")
    v = Variable(Int, "v")

    # Pipeline Parameters
    R = Parameter(Int, "R")  # image rows
    C = Parameter(Int, "C")  # image cols
    kelvin = Parameter(Float, "colour_temp")  # temperature
    gamma = Parameter(Float, "gamma")  # gamma value
    contrast = Parameter(Float, "contrast")  # colour contrast

    # register in the dictionary
    pipe_data['R'] = R
    pipe_data['C'] = C

    # Intervals:

    # shifts in each dimension
    ghost_x = 12
    ghost_y = 16

    # camera_pipe intervals :
    # bounds for input image
    row = Interval(Int, 0, R-1)
    col = Interval(Int, 0, C-1)
    # bounds for denoise function
    ghost_zone_2x = Interval(Int, (ghost_x-4), (R-24-ghost_x+6) - 1)
    ghost_zone_2y = Interval(Int, (ghost_y-4), (C-32-ghost_y+6) - 1)
    # ghost zone without any offset(ghost }:))
    ghost_zone_0x = Interval(Int, (ghost_x-0), (R-24-ghost_x+0) - 1)
    ghost_zone_0y = Interval(Int, (ghost_y-0), (C-32-ghost_y+0) - 1)
    # bounds for g_gr, r_r, b_b, g_gb
    half_ghost_zone_2x = Interval(Int, ((ghost_x)//2) - 2, ((R-24-ghost_x)//2 + 2))
    half_ghost_zone_2y = Interval(Int, ((ghost_y)//2) - 2, ((C-32-ghost_y)//2 + 2))
    # bounds for g_r, g_b
    half_ghost_zone_1x = Interval(Int, ((ghost_x)//2) - 1, ((R-24-ghost_x)//2 + 1))
    half_ghost_zone_1y = Interval(Int, ((ghost_y)//2) - 1, ((C-32-ghost_y)//2 + 1))
    # bounds for b_r, b_gr, r_gr, b_gb, r_gb, r_b
    half_ghost_zone_0x = Interval(Int, ((ghost_x)//2) - 0, ((R-24-ghost_x)//2 + 0))
    half_ghost_zone_0y = Interval(Int, ((ghost_y)//2) - 0, ((C-32-ghost_y)//2 + 0))
    # bounds for colour channels
    rgb = Interval(Int, 0, 2)
    grbg = Interval(Int, 0, 3)
    # bound for LUT
    lut_range = Interval(Int, -32768, 32767)

    # Image Inputs
    matrix_3200 = Image(Float, "matrix_3200", [3, 4])
    matrix_7000 = Image(Float, "matrix_7000", [3, 4])
    img = Image(Short, "img", [R, C])

    # Pipeline
    # ========

    # 1. Hot Pixel Suppression / Denoising
    x_max = Max(img(x-2, y), img(x+2, y))
    y_max = Max(img(x, y-2), img(x, y+2))
    x_min = Min(img(x-2, y), img(x+2, y))
    y_min = Min(img(x, y-2), img(x, y+2))
    max_ = Max(x_max, y_max)
    min_ = Min(x_min, y_min)
    clamp = Min(max_, img(x, y))
    clamp = Max(min_, clamp)

    denoised = Function(([x, y], [ghost_zone_2x, ghost_zone_2y]), \
                        Short, "denoised")
    denoised.defn = [clamp]

    # 2. Deinterleave the Bayer Array
    deinterleaved = \
        Function(([c, x, y], [grbg, half_ghost_zone_2x, half_ghost_zone_2y]), \
                 Short, "deinterleaved")
    deinterleaved.defn = [Case(Condition(c, '==', 0), denoised(2*x, 2*y)),
                          Case(Condition(c, '==', 1), denoised(2*x, 2*y+1)),
                          Case(Condition(c, '==', 2), denoised(2*x+1, 2*y)),
                          Case(Condition(c, '==', 3), denoised(2*x+1, 2*y+1))]

    # 3. Perform Demosaicing on the Deinterleaved array
    #
    # Halide :
    # "These are the values we already know from the input
    # x_y = the value of channel x at a site in the input of channel y
    # gb refers to green sites in the blue rows
    # gr refers to green sites in the red rows.
    # We'll give more convenient names to the four channels we know"

    g_gr = Function(([x, y], [half_ghost_zone_2x, half_ghost_zone_2y]), \
                    Short, "g_gr")
    g_gr.defn = [deinterleaved(0, x, y)]

    r_r = Function(([x, y], [half_ghost_zone_2x, half_ghost_zone_2y]), \
                   Short, "r_r")
    r_r.defn = [deinterleaved(1, x, y)]

    b_b = Function(([x, y], [half_ghost_zone_2x, half_ghost_zone_2y]), \
                   Short, "b_b")
    b_b.defn = [deinterleaved(2, x, y)]

    g_gb = Function(([x, y], [half_ghost_zone_2x, half_ghost_zone_2y]), \
                    Short, "g_gb")
    g_gb.defn = [deinterleaved(3, x, y)]

    # Halide :
    # "These are the ones we need to interpolate
    # b_r, g_r, b_gr, r_gr, b_gb, r_gb, r_b, g_b
    #
    # First calculate green at the red and blue sites
    #
    # Try interpolating vertically and horizontally. Also compute
    # differences vertically and horizontally. Use interpolation in
    # whichever direction had the smallest difference."

    def absd(a, b):
        return Select(Condition(a, '>', b), a - b, b - a)

    #
    gv_r  =    (g_gb(x-1, y) + g_gb(x, y))/2
    gvd_r = absd(g_gb(x-1, y), g_gb(x, y))
    gh_r  =    (g_gr(x, y+1) + g_gr(x, y))/2
    ghd_r = absd(g_gr(x, y+1), g_gr(x, y))

    g_r = Function(([x, y], [half_ghost_zone_1x, half_ghost_zone_1y]), \
                   Short, "g_r")
    g_r.defn = [Select(Condition(ghd_r, '<', gvd_r), gh_r, gv_r)]

    #
    gv_b  =    (g_gr(x+1, y) + g_gr(x, y))/2
    gvd_b = absd(g_gr(x+1, y), g_gr(x, y))
    gh_b  =    (g_gb(x, y-1) + g_gb(x, y))/2
    ghd_b = absd(g_gb(x, y-1), g_gb(x, y))

    g_b = Function(([x, y], [half_ghost_zone_1x, half_ghost_zone_1y]), \
                   Short, "g_b")
    g_b.defn = [Select(Condition(ghd_b, '<', gvd_b), gh_b, gv_b)]

    # Halide :
    # "Next interpolate red at gr by first interpolating, then
    # correcting using the error green would have had if we had
    # interpolated it in the same way (i.e. add the second derivative
    # of the green channel at the same place)."

    correction = g_gr(x, y) - (g_r(x, y) + g_r(x, y-1))/2
    r_gr = Function(([x, y], [half_ghost_zone_0x, half_ghost_zone_0y]), \
                    Short, "r_gr")
    r_gr.defn = [correction + (r_r(x, y-1) + r_r(x, y))/2]

    # Halide : "Do the same for other reds and blues at green sites"
    correction = g_gr(x, y) - (g_b(x, y) + g_b(x-1, y))/2
    b_gr = Function(([x, y], [half_ghost_zone_0x, half_ghost_zone_0y]), \
                    Short, "b_gr")
    b_gr.defn = [correction + (b_b(x, y) + b_b(x-1, y))/2]

    correction = g_gb(x, y) - (g_r(x, y) + g_r(x+1, y))/2
    r_gb = Function(([x, y], [half_ghost_zone_0x, half_ghost_zone_0y]), \
                    Short, "r_gb")
    r_gb.defn = [correction + (r_r(x, y) + r_r(x+1, y))/2]

    correction = g_gb(x, y) - (g_b(x, y) + g_b(x, y+1))/2
    b_gb = Function(([x, y], [half_ghost_zone_0x, half_ghost_zone_0y]), \
                    Short, "b_gb")
    b_gb.defn = [correction + (b_b(x, y) + b_b(x, y+1))/2]

    # Halide:
    # "Now interpolate diagonally to get red at blue and blue at
    # red. Hold onto your hats; this gets really fancy. We do the
    # same thing as for interpolating green where we try both
    # directions (in this case the positive and negative diagonals),
    # and use the one with the lowest absolute difference. But we
    # also use the same trick as interpolating red and blue at green
    # sites - we correct our interpolations using the second
    # derivative of green at the same sites."

    correction = g_b(x, y)  - (g_r(x, y) + g_r(x+1, y-1))/2
    rp_b  = correction + (r_r(x, y) + r_r(x+1, y-1))/2
    rpd_b = absd(r_r(x, y), r_r(x+1, y-1))

    correction = g_b(x, y)  - (g_r(x, y-1) + g_r(x+1, y))/2
    rn_b  = correction + (r_r(x, y-1) + r_r(x+1, y))/2
    rnd_b = absd(r_r(x, y-1), r_r(x+1, y))

    r_b = Function(([x, y], [half_ghost_zone_0x, half_ghost_zone_0y]), \
                   Short, "r_b")
    r_b.defn = [Select(Condition(rpd_b, '<', rnd_b), rp_b, rn_b)]

    # Halide : "Same thing for blue at red"
    correction = g_r(x, y)  - (g_b(x, y) + g_b(x-1, y+1))/2;
    bp_r  = correction + (b_b(x, y) + b_b(x-1, y+1))/2;
    bpd_r = absd(b_b(x, y), b_b(x-1, y+1));

    correction = g_r(x, y)  - (g_b(x, y+1) + g_b(x-1, y))/2;
    bn_r  = correction + (b_b(x, y+1) + b_b(x-1, y))/2;
    bnd_r = absd(b_b(x, y+1), b_b(x-1, y));

    b_r = Function(([x, y], [half_ghost_zone_0x, half_ghost_zone_0y]), \
                   Short, "b_r")
    b_r.defn = [Select(Condition(bpd_r, '<', bnd_r), bp_r, bn_r)]

    # 4. Interleave the resulting channels
    def interleave_x(a, b, name):
        out = Function(([x, y], [ghost_zone_0x, ghost_zone_0y]), \
                       Short, name)
        out.defn = [Select(Condition((x%2), '==', 0), a(x/2, y), b(x/2, y))]
        return out

    def interleave_y(a, b, name):
        out = Function(([x, y], [half_ghost_zone_0x, ghost_zone_0y]), \
                       Short, name)
        out.defn = [Select(Condition((y%2), '==', 0), a(x, y/2), b(x, y/2))]
        return out

    red = interleave_x(interleave_y(r_gr, r_r, "red_x1"), \
                       interleave_y(r_b, r_gb, "red_x2"), \
                       "red")(x, y)
    green = interleave_x(interleave_y(g_gr, g_r, "green_x1"), \
                         interleave_y(g_b, g_gb, "green_x2"), \
                         "green")(x, y)
    blue = interleave_x(interleave_y(b_gr, b_r, "blue_x1"), \
                        interleave_y(b_b, b_gb, "blue_x2"), \
                        "blue")(x, y)

    # 5. Colour Correction
    #
    # Halide :
    # "Get a color matrix by linearly interpolating between two
    # calibrated matrices using inverse kelvin."

    alpha = (1.0/kelvin - 1.0/3200) / (1.0/7000 - 1.0/3200)

    #def matrix(i, j):
    #    val = (matrix_3200(i, j) * alpha + matrix_7000(i, j) * (1 - alpha))
    #    val = Cast(Int, (val * 256.0)) # Halide : "Q8.8 fixed point"
    #    return val

    mat = Function(([x, y], [rgb, grbg]), Int, "matrix")
    val = (matrix_3200(x, y) * alpha + matrix_7000(x, y) * (1 - alpha))
    mat.defn = [Cast(Int, val * 256.0)]

    r = mat(0, 3) + mat(0, 0) * red + mat(0, 1) * green + mat(0, 2) * blue
    g = mat(1, 3) + mat(1, 0) * red + mat(1, 1) * green + mat(1, 2) * blue
    b = mat(2, 3) + mat(2, 0) * red + mat(2, 1) * green + mat(2, 2) * blue

    corrected = Function(([c, x, y], [rgb, ghost_zone_0x, ghost_zone_0y]), \
                         Short, "corrected")
    corrected.defn = [ Case(Condition(c, '==', 2), r/256),
                       Case(Condition(c, '==', 1), g/256),
                       Case(Condition(c, '==', 0), b/256) ]

    # 6. Apply Curve
    #
    # Halide :
    # "copied from FCam"

    def lut_value(anything):
        xf = Cast(Float, anything/1024.0)

        # clamp
        xf = Min(1.0, xf)
        xf = Max(0.0, xf)

        g = Powf(xf, 1.0/gamma)
        b = 2.0 - Powf(2.0, contrast/100.0)
        a = 2.0 - 2.0*b
        z = Select(Condition(g, '>', 0.5), \
                            1.0 - (a*(1.0-g)*(1.0-g) + b*(1.0-g)), \
                            a*g*g + b*g)
    
        val = Cast(UChar, Min(Max(z*256.0, 0.0), 255.0))

        return val

    # compute this beforehand (its a LUT)
    curve = Function(([v], [lut_range]), UChar, "curveLUT")
    curve.defn = [lut_value(v)]

    # pick from LUT map
    curved = Function(([c, x, y], [rgb, ghost_zone_0x, ghost_zone_0y]), \
                      UChar, "process")
    # (1) compute only those out of inLUT on-the-fly
    # inLUT = Condition(corrected(c, x, y), '<=', 65535) & \
    #         Condition(corrected(c, x, y), '>=', 0)
    # curved.defn = [Select(inLUT, curve(corrected(c, x, y)), \
    #                              lut_value(corrected(c, x, y)))]
    # (2) with correct range [ -2^15 : 2^15 ]
    curved.defn = [curve(corrected(c, x, y))]
    # (3) compute everything on-the-fly
    # curved.defn = [lut_value(corrected(c, x, y))]

    return curved

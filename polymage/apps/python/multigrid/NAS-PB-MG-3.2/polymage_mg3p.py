from __init__ import *

from polymage_common import set_zero_ghosts
from polymage_residual import residual
from polymage_restrict import restrict
from polymage_solver import psinv
from polymage_interp import interpolate

import sys
sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def mg3p(app_data):
    pipe_data = app_data['pipe_data']

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    lt = app_data['lt']
    lb = app_data['lb']
    n  = pipe_data['n']

    U = Image(Double, "U", [n+2, n+2, n+2])
    V = Image(Double, "V", [n+2, n+2, n+2])
    R = Image(Double, "R", [n+2, n+2, n+2])

    a = app_data['a']
    c = app_data['c']

    extent = pipe_data['extent']
    interior = pipe_data['interior']
    ghosts = pipe_data['ghosts']

    # from finest to the coarsest level
    # restrict the residual to the next coarser level
    restr = {}
    restr[lt] = R
    for l in range(lt, lb, -1):
        restr[l-1] = restrict(restr[l], l, pipe_data, "restrict_"+str(l-1))

    # compute an approximate solution on the coarsest grid
    smooth = {}
    smooth[lb] = psinv(restr[lb], None, lb, app_data, "psinv_"+str(lb))

    residl = {}

    intrp = {}
    for l in range(lb+1, lt):
        # prolongate from level l-1 to l
        intrp[l] = \
          interpolate(smooth[l-1], None, l, pipe_data, "interp_"+str(l))

        # compute residual for level l
        residl[l] = \
          residual(intrp[l], restr[l], l, app_data, "resid_"+str(l))

        # apply smoother
        smooth[l] = \
          psinv(residl[l], intrp[l], l, app_data, "psinv_"+str(l))
    # endfor

    # prolongate from level lt-1 to lt
    intrp[lt] = interpolate(smooth[lt-1], U, lt, pipe_data, "interp_"+str(lt))

    # compute residual for level lt
    residl[lt] = residual(intrp[lt], V, lt, app_data, "MG_R")

    # apply smoother
    smooth[lt] = psinv(residl[lt], intrp[lt], lt, app_data, "MG_U")

    return smooth[lt], residl[lt]

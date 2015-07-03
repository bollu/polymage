import sys
from polymage_common   import setZeroGhosts
from polymage_residual import residual
from polymage_restrict import restrict
from polymage_solver   import psinv
from polymage_interp   import interpolate

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def mg3P(u, v, r, impipeDict, dataDict):
    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    lt = dataDict['lt']
    lb = dataDict['lb']
    N  = dataDict['N']

    a = dataDict['a']
    c = dataDict['c']

    extent = impipeDict['extent']
    interior = impipeDict['interior']
    ghosts = impipeDict['ghosts']

    # from finest to the coarsest level
    # restrict the residual to the next coarser level
    restr = {}
    restr[lt] = r
    for l in range(lt, lb, -1):
        restr[l-1] = restrict(restr[l], l, impipeDict, "restrict_"+str(l-1))

    # compute an approximate solution on the coarsest grid
    smooth = {}
    smooth[lb] = psinv(restr[lb], None, lb, impipeDict, dataDict, "psinv_"+str(lb))

    residl = {}

    intrp = {}
    for l in range(lb+1, lt):
        # prolongate from level l-1 to l
        intrp[l] = interpolate(smooth[l-1], None, l, impipeDict, "interp_"+str(l))

        # compute residual for level l
        residl[l] = residual(intrp[l], restr[l], l, impipeDict, dataDict, "resid_"+str(l))

        # apply smoother
        smooth[l] = psinv(residl[l], intrp[l], l, impipeDict, dataDict, "psinv_"+str(l))
    # endfor

    # prolongate from level lt-1 to lt
    intrp[lt] = interpolate(smooth[lt-1], u, lt, impipeDict, "interp_"+str(lt))

    # compute residual for level lt
    residl[lt] = residual(intrp[lt], v, lt, impipeDict, dataDict, "mgR_")

    # apply smoother
    smooth[lt] = psinv(residl[lt], intrp[lt], lt, impipeDict, dataDict, "mgU_")

    return smooth[lt], residl[lt]

def mg3pPipe(impipeDict, dataDict):
    n = impipeDict['n']

    u = Image(Double, "u", [n, n, n])
    v = Image(Double, "v", [n, n, n])
    r = Image(Double, "r", [n, n, n])

    mg_u, mg_r = mg3P(u, v, r, impipeDict, dataDict)

    return mg_u, mg_r

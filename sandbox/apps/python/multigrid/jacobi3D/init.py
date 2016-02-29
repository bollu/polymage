import sys
import os.path
import numpy as np
from arg_parser import parse_args
from printer         import printHeader, printUsage, \
                            printLine
from polymage_common import setVars, setCases
from execMG          import calcNorm

def initNorm(appData):
    grid_data = appData['grid_data']
    U_ = grid_data['U_']

    appData['resid'] = 0.0
    appData['err']   = 0.0

    # calculate the initial residual norm and error
    print("[init]: calculating the initial norm and error ...")
    calcNorm(U_, appData)
    print("[init]: ... DONE")

    return

def initBorder(grid, borderWidth, borderValues):
    # size of the grid - assumed to be same in
    # all the dimesnions
    n = grid.shape[0]

    w = borderWidth
    v = borderValues

    # z-planes
    grid[  0:w] = v
    grid[n-w:n] = v

    # y-planes
    grid[0:n,   0:w] = v
    grid[0:n, n-w:n] = v

    # x-planes
    grid[0:n, 0:n,   0:w] = v
    grid[0:n, 0:n, n-w:n] = v

    return

def initGrids(appData):
    print("[init_mg.py] : grids")

    N = appData['N']

    # working grid (even step)
    U_ = np.ones((N+2, N+2, N+2), np.float64)
    initBorder(U_, borderWidth=1, borderValues=0.0)

    # working grid (odd step)
    W_ = np.zeros((N+2, N+2, N+2), np.float64)
    initBorder(W_, borderWidth=1, borderValues=0.0)

    # RHS
    F_ = np.zeros((N+2, N+2, N+2), np.float64)

    # exact solution
    U_EXACT_ = np.zeros((N+2, N+2, N+2), np.float64)

    grid_data = {}
    grid_data['U_']       = U_
    grid_data['W_']       = W_
    grid_data['F_']       = F_
    grid_data['U_EXACT_'] = U_EXACT_

    appData['grid_data'] = grid_data

    return 

def initParams(app_args,appData):
    print("[init_mg.py] : parameters")

    # size of each dimension of the coarsest grid
    n = int(app_args.n)
    # number of multigrid levels
    L = int(app_args.L)

    N = n
    # compute the size of the finest grid
    for l in range(0,L):
        N = 2*N+1

    appData['n'] = n
    appData['N'] = N
    appData['L'] = L

    # pre-smoother, post-smoother and
    # coarse-grid relaxation steps
    appData['nu1'] = int(app_args.nu1)
    appData['nuc'] = int(app_args.nuc)
    appData['nu2'] = int(app_args.nu2)

    # pool allocate option
    appData['pool_alloc'] = False

    assert not (appData['nu1'] == 0 and \
                appData['nu2'] == 0 and
                appData['nuc'] == 0)

    return appData

def getInput(app_args,appData):
    appData['app_args'] = app_args
    appData['mode'] = app_args.mode
    appData['cycle'] = app_args.cycle
    appData['nit']  = int(app_args.nit)
    cycle_name = appData['cycle']+"cycle"
    appData['cycle_name'] = app_args.cycle_name

    appData['timer'] = app_args.timer

    return

def initAll(pipeData, appData):
    # TODO init cycle type {V, W}
    app_args = parse_args()
    getInput(app_args,appData)

    initParams(app_args,appData)

    initGrids(appData)

    setVars(pipeData, appData)

    setCases(pipeData, appData)

    return

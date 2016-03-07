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
    # all the dimensions
    n = grid.shape[0]

    w = borderWidth
    v = borderValues

    # y-planes
    grid[  0:w] = v
    grid[n-w:n] = v

    # x-planes
    grid[0:n,   0:w] = v
    grid[0:n, n-w:n] = v

    return

def initBorderPiecewise(grid, borderWidth, borderValues):
    w = borderWidth
    v = borderValues

    N = grid.shape[0]
    n = N-2*w+1

    # copy all the nine pieces
    grid[0:w, 0:w] = v[(0, 0)][0:w, 0:w]
    grid[0:w, w:n] = v[(0, 1)][0:w, w:n]
    grid[0:w, n:N] = v[(0, 2)][0:w, n:N]

    grid[w:n, 0:w] = v[(1, 0)][w:n, 0:w]
    #grid[w:n, w:n] = v[(1, 1)][w:n, w:n]
    grid[w:n, n:N] = v[(1, 2)][w:n, n:N]

    grid[n:N, 0:w] = v[(2, 0)][n:N, 0:w]
    grid[n:N, w:n] = v[(2, 1)][n:N, w:n]
    grid[n:N, n:N] = v[(2, 2)][n:N, n:N]

    return

def initGrids(app_args,appData):
    print("[init_mg.py] : grids")

    N = appData['N']
    problem = appData['problem']

    # working grid (even step)
    U_ = np.ones((N+2, N+2), np.float64)
    # working grid (odd step)
    W_ = np.zeros((N+2, N+2), np.float64)

    if problem == 1:
        initBorder(U_, borderWidth=1, borderValues=0.0)
        initBorder(W_, borderWidth=1, borderValues=0.0)
        # RHS
        F_ = np.zeros((N+2, N+2), np.float64)
        # exact solution
        U_EXACT_ = np.zeros((N+2, N+2), np.float64)
    else:
        #'''
        N1 = {}
        L = int(app_args.L) 
        n = int(app_args.n)
        for l in range(0,L+1):
            if l == 0:
                N1[0] = n
            else:
                N1[l] = 2*N1[l-1]+1

            h = 1.0/(N1[l]+1)
        #'''
        indices = np.indices((N+2, N+2))
        x = indices[0] * h
        y = indices[1] * h

        x_1 = np.array(x)
        x_1 = x_1 + 1.0
        y_1 = np.array(y)
        y_1 = y_1 + 1.0

        x_y = x + y
        if problem == 2:
            b = {(0, 0):x, (0, 1):y  , (0, 2):x_1, \
                 (1, 0):x, (1, 1):U_ , (1, 2):x_1, \
                 (2, 0):x, (2, 1):y_1, (2, 2):x_1}
            initBorderPiecewise(U_, borderWidth=1, borderValues=b)
            initBorderPiecewise(W_, borderWidth=1, borderValues=b)
            # initialize F_
            F_ = np.zeros((N+2, N+2), np.float64) # f(i, j) <- 0
            # initialize U_
            U_EXACT_ = x_y
        elif problem == 3:
            sx = np.sin(x)
            sy = np.sin(y)

            sx_1 = np.sin(x_1)
            sy_1 = np.sin(y_1)

            sx_y = np.sin(x_y)

            b = {(0, 0):sx, (0, 1):sy  , (0, 2):sx_1, \
                 (1, 0):sx, (1, 1):U_  , (1, 2):sx_1, \
                 (2, 0):sx, (2, 1):sy_1, (2, 2):sx_1}
            initBorderPiecewise(U_, borderWidth=1, borderValues=b)
            initBorderPiecewise(W_, borderWidth=1, borderValues=b)
            # initialize F_
            F_ = -2.0 * sx_y
            # initialize U_EXACT_
            U_EXACT_ = sx_y
        elif problem == 4:
            xx = x*x
            yy = y*y

            xx_1 = xx+1.0
            yy_1 = yy+1.0

            xx_yy = xx + yy

            b = {(0, 0):xx, (0, 1):yy  , (0, 2):xx_1, \
                 (1, 0):xx, (1, 1):U_  , (1, 2):xx_1, \
                 (2, 0):xx, (2, 1):yy_1, (2, 2):xx_1}
            initBorderPiecewise(U_, borderWidth=1, borderValues=b)
            initBorderPiecewise(W_, borderWidth=1, borderValues=b)
            # initialize f
            F_ = np.zeros((N+2, N+2), np.float64) # f(i, j) <- 0
            # initialize u
            U_EXACT_ = xx_yy

    grid_data = {}
    grid_data['U_']       = U_
    grid_data['W_']       = W_
    grid_data['F_']       = F_
    grid_data['U_EXACT_'] = U_EXACT_

    appData['grid_data'] = grid_data

    return 

def initParams(app_args, appData):
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
    appData['runs'] = int(app_args.runs)

    # pre-smoother, post-smoother and
    # coarse-grid relaxation steps
    appData['nu1'] = int(app_args.nu1)
    appData['nuc'] = int(app_args.nuc)
    appData['nu2'] = int(app_args.nu2)

    # problem type
    appData['problem'] = int(app_args.problem)

    # pool allocate option
    appData['pool_alloc'] = app_args.pool_alloc

    assert not (appData['nu1'] == 0 and \
                appData['nu2'] == 0 and
                appData['nuc'] == 0)

    return appData

def getInput(app_args, appData):
    appData['app_args'] = app_args
    appData['mode'] = app_args.mode
    appData['cycle'] = app_args.cycle
    appData['nit'] = int(app_args.nit)

    cycle_name = appData['cycle']+"cycle"
    appData['cycle_name'] = cycle_name
    appData['timer'] = app_args.timer
  
    return

def initAll(impipeData, appData):
    # TODO init cycle type {V, W}
    app_args = parse_args()
    getInput(app_args,appData)

    initParams(app_args,appData)

    initGrids(app_args,appData)

    setVars(impipeData, appData)

    setCases(impipeData, appData)

    return

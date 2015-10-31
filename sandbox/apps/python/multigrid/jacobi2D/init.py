import sys
import os.path
import numpy as np

from printer         import printHeader, printUsage, \
                            printLine
from polymage_common import setVars, setCases
from execMG          import calcNorm

def initNorm(dataDict):
    gridDict = dataDict['gridDict']
    U_ = gridDict['U_']

    dataDict['resid'] = 0.0
    dataDict['err']   = 0.0

    # calculate the initial residual norm and error
    print("[init]: calculating the initial norm and error ...")
    calcNorm(U_, dataDict)
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

def initGrids(dataDict):
    print("[init_mg.py] : grids")

    N = dataDict['N']
    problem = dataDict['problem']

    # working grid (even step)
    U_ = np.ones((N+2, N+2), np.float64)
    # working grid (odd step)
    W_ = np.ones((N+2, N+2), np.float64)

    if problem == 1:
        initBorder(U_, borderWidth=1, borderValues=0.0)
        initBorder(W_, borderWidth=1, borderValues=0.0)
        # RHS
        F_ = np.zeros((N+2, N+2), np.float64)
        # exact solution
        U_EXACT_ = np.zeros((N+2, N+2), np.float64)
    else:
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
            U_ = xx_yy

    gridDict = {}
    gridDict['U_']       = U_
    gridDict['W_']       = W_
    gridDict['F_']       = F_
    gridDict['U_EXACT_'] = U_EXACT_

    dataDict['gridDict'] = gridDict

    return 

def initParams(dataDict):
    print("[init_mg.py] : parameters")

    # size of each dimension of the coarsest grid
    n = 255
    # number of multigrid levels
    L = 2

    N = n
    # compute the size of the finest grid
    for l in range(0,L):
        N = 2*N+1

    dataDict['n'] = n
    dataDict['N'] = N
    dataDict['L'] = L

    # pre-smoother, post-smoother and
    # coarse-grid relaxation steps
    dataDict['nu1'] = 10
    dataDict['nuc'] = 0
    dataDict['nu2'] = 0

    # problem type
    dataDict['problem'] = 1

    # pool allocate option
    dataDict['pool_alloc'] = False

    assert not (dataDict['nu1'] == 0 and \
                dataDict['nu2'] == 0 and
                dataDict['nuc'] == 0)

    return dataDict

def getInput(dataDict):
    if len(sys.argv) > 3:
        dataDict['mode'] = sys.argv[1]
        dataDict['cycle'] = sys.argv[2]
        dataDict['nit']  = int(sys.argv[3])
    else:
        printUsage()
        sys.exit(1)

    cycle_name = dataDict['cycle']+"cycle"
    dataDict['cycle_name'] = cycle_name

    dataDict['timer'] = os.path.isfile("timer.flag")

    return

def initAll(impipeDict, dataDict):
    # TODO init cycle type {V, W, S}

    getInput(dataDict)

    initParams(dataDict)

    initGrids(dataDict)

    setVars(impipeDict, dataDict)

    setCases(impipeDict, dataDict)

    return

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

def initGrids(dataDict):
    print("[init_mg.py] : grids")

    N = dataDict['N']

    # working grid (even step)
    U_ = np.ones((N+2, N+2, N+2), np.float64)
    initBorder(U_, borderWidth=1, borderValues=0.0)

    # working grid (odd step)
    W_ = np.zeros((N+2, N+2, N+2), np.float64)

    # RHS
    F_ = np.zeros((N+2, N+2, N+2), np.float64)

    # exact solution
    U_EXACT_ = np.zeros((N+2, N+2, N+2), np.float64)

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
    n = 15
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

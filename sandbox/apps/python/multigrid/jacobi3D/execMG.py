import sys
import os
import ctypes
import numpy as np
import time

from printer import printLine, printLayout, printErrors

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def calcNorm(U_, dataDict):
    N = dataDict['N']

    gridDict = dataDict['gridDict']
    F_       = gridDict['F_']
    U_EXACT_ = gridDict['U_EXACT_']

    # lib function name
    norm = dataDict['pipeline_norm']

    resid = np.zeros((1), np.float64)
    err   = np.zeros((1), np.float64)

    # lib function args
    normArgs = []
    normArgs += [ctypes.c_int(N)]
    normArgs += [ctypes.c_void_p(      F_.ctypes.data)]
    normArgs += [ctypes.c_void_p(U_EXACT_.ctypes.data)]
    normArgs += [ctypes.c_void_p(      U_.ctypes.data)]
    normArgs += [ctypes.c_void_p(     err.ctypes.data)]
    normArgs += [ctypes.c_void_p(   resid.ctypes.data)]

    # call lib function
    norm(*normArgs)

    # save the old norm values
    dataDict['oldResidual'] = dataDict['resid']
    dataDict['oldErr']      = dataDict['err']

    # register the norm values in the data dictionary
    dataDict['resid'] = resid[0]
    dataDict['err']   = err[0]

    return

def callVCycle(U_, W_, dataDict):
    n = dataDict['n']

    gridDict = dataDict['gridDict']
    F_       = gridDict['F_']

    # lib function name
    vCycleFunc = dataDict['pipeline_vcycle']

    # lib function args
    vCycleArgs = []
    vCycleArgs += [ctypes.c_int(n)]
    vCycleArgs += [ctypes.c_void_p(F_.ctypes.data)]
    vCycleArgs += [ctypes.c_void_p(U_.ctypes.data)]
    vCycleArgs += [ctypes.c_void_p(W_.ctypes.data)]

    # call lib function
    vCycleFunc(*vCycleArgs)

    return

def multigrid(dataDict):
    gridDict = dataDict['gridDict']
    U_ = gridDict['U_']
    W_ = gridDict['W_']

    nit = dataDict['nit']
    it  = 0

    printLayout(dataDict)
    printErrors(it, dataDict)

    timer = dataDict['timer']
    if timer == True:
        t1 = time.clock()

    while it < nit :
        it += 1
        if it%2 == 1:
            callVCycle(U_, W_, dataDict)
            if timer == False:
                calcNorm(W_, dataDict)
        else:
            callVCycle(W_, U_, dataDict)
            if timer == False:
                calcNorm(U_, dataDict)

        if timer == False:
            printErrors(it, dataDict)

    if timer == True:
        t2 = time.clock()

        timeTaken = float(t2) - float(t1)
        print
        print '[execMG] : time taken to execute = ', timeTaken*1000, ' ms'

    return

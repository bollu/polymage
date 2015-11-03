import sys
import os
import ctypes
import numpy as np
import time

from printer import printLine, printLayout, printErrors

from compiler   import *
from constructs import *

def minimal_exec_mg(pipe_lib, pipe_lib_func, func_params,
                    func_args, tuner_data, app_data=None):
    pipe_arg_data = tuner_dict['_tuner_pipe_arg_data']

    it = 0
    it_max = app_data['nit']

    pool_alloc = dataDict['pool_alloc']  # bool

    # build function argument list based on the iteration,
    # even : in = U_ : out = W_
    # odd  : in = W_ : out = U_
    func_args = []

    arg_data = {}
    arg_data['n'] = app_data['n']
    arg_data['U_'] = app_data['U_']

    arg_data['W_'] = app_data['W_']
    arg_data['F_'] = app_data['F_']
    func_args.append(map_c_func_args(func_params, arg_data))

    arg_data['U_'] = app_data['W__']
    arg_data['W_'] = app_data['U_']
    func_args.append(map_c_func_args(func_params, arg_data))

    if pool_alloc:
        pipe_lib.pool_init()

    while it < itMax:
        pipe_lib_func(*(func_args[it%2]))
        it += 1

    if pool_alloc:
        pipe_lib.pool_destroy()

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

def callMGCycle(U_, W_, dataDict):
    n = dataDict['n']

    gridDict = dataDict['gridDict']
    F_       = gridDict['F_']

    # lib function name
    func_name = 'pipeline_'+dataDict['cycle_name']
    mgCycleFunc = dataDict[func_name]

    # lib function args
    mgCycleArgs = []
    mgCycleArgs += [ctypes.c_int(n)]
    mgCycleArgs += [ctypes.c_void_p(F_.ctypes.data)]
    mgCycleArgs += [ctypes.c_void_p(U_.ctypes.data)]
    mgCycleArgs += [ctypes.c_void_p(W_.ctypes.data)]

    # call lib function
    mgCycleFunc(*mgCycleArgs)

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
        t1 = time.time()

    while it < nit :
        it += 1
        if it%2 == 1:
            callMGCycle(U_, W_, dataDict)
            if timer == False:
                calcNorm(W_, dataDict)
        else:
            callMGCycle(W_, U_, dataDict)
            if timer == False:
                calcNorm(U_, dataDict)

        if timer == False:
            printErrors(it, dataDict)

    if timer == True:
        t2 = time.time()

        timeTaken = float(t2) - float(t1)
        print("")
        print("[execMG] : time taken to execute = ", timeTaken, " ms")

    return

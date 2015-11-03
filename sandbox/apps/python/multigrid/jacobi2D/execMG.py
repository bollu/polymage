import sys
import os
import ctypes
import numpy as np
import time

from printer import printLine, printLayout, printErrors

from compiler   import *
from constructs import *
from utils import *

def minimal_exec_mg(pipe_lib, pipe_lib_func, func_params,
                    func_args, tuner_data, app_data):
    it = 0
    it_max = app_data['nit']

    pool_alloc = app_data['pool_alloc']  # bool
    grid_data = app_data['grid_data']

    # build function argument list based on the iteration,
    # even : in = U_ : out = W_
    # odd  : in = W_ : out = U_
    func_args = []

    arg_data = {}
    arg_data['n'] = app_data['n']
    arg_data['F_'] = grid_data['F_']

    arg_data['V_'] = grid_data['U_']
    arg_data['Vcycle'] = grid_data['W_']
    func_args.append(map_cfunc_args(func_params, arg_data))

    arg_data['V_'] = grid_data['W_']
    arg_data['Vcycle'] = grid_data['U_']
    func_args.append(map_cfunc_args(func_params, arg_data))

    '''
    if pool_alloc:
        pipe_lib.pool_init()
    '''

    while it < it_max:
        pipe_lib_func(*(func_args[it%2]))
        it += 1

    '''
    if pool_alloc:
        pipe_lib.pool_destroy()
    '''

    return

def calcNorm(U_, app_data):
    N = app_data['N']

    grid_data = app_data['grid_data']
    F_       = grid_data['F_']
    U_EXACT_ = grid_data['U_EXACT_']

    # lib function name
    norm = app_data['pipeline_norm']

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
    app_data['oldResidual'] = app_data['resid']
    app_data['oldErr']      = app_data['err']

    # register the norm values in the data dictionary
    app_data['resid'] = resid[0]
    app_data['err']   = err[0]

    return

def callMGCycle(U_, W_, app_data):
    n = app_data['n']

    grid_data = app_data['grid_data']
    F_       = grid_data['F_']

    # lib function name
    func_name = 'pipeline_'+app_data['cycle_name']
    mgCycleFunc = app_data[func_name]

    # lib function args
    mgCycleArgs = []
    mgCycleArgs += [ctypes.c_int(n)]
    mgCycleArgs += [ctypes.c_void_p(F_.ctypes.data)]
    mgCycleArgs += [ctypes.c_void_p(U_.ctypes.data)]
    mgCycleArgs += [ctypes.c_void_p(W_.ctypes.data)]

    # call lib function
    mgCycleFunc(*mgCycleArgs)

    return

def multigrid(app_data):
    grid_data = app_data['grid_data']
    U_ = grid_data['U_']
    W_ = grid_data['W_']

    nit = app_data['nit']
    it  = 0

    printLayout(app_data)
    printErrors(it, app_data)

    timer = app_data['timer']
    if timer == True:
        t1 = time.time()

    while it < nit :
        it += 1
        if it%2 == 1:
            callMGCycle(U_, W_, app_data)
            if timer == False:
                calcNorm(W_, app_data)
        else:
            callMGCycle(W_, U_, app_data)
            if timer == False:
                calcNorm(U_, app_data)

        if timer == False:
            printErrors(it, app_data)

    if timer == True:
        t2 = time.time()

        timeTaken = float(t2) - float(t1)
        print("")
        print("[execMG] : time taken to execute = ", timeTaken, " ms")

    return

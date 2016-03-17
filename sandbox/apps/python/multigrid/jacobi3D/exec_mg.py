from __init__ import *

import sys
import os
import ctypes
import numpy as np
import time

from printer import print_line, print_layout, print_errors

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
    n = app_data['n']
    arg_data['n'] = n
    arg_data['F_'] = grid_data['F_']

    arg_data['V_'] = grid_data['U_']
    arg_data['Vcycle'] = grid_data['W_']
    func_args.append(map_cfunc_args(func_params, arg_data))

    arg_data['V_'] = grid_data['W_']
    arg_data['Vcycle'] = grid_data['U_']
    func_args.append(map_cfunc_args(func_params, arg_data))

    if pool_alloc:
        pipe_lib.pool_init()

    while it < it_max:
        pipe_lib_func(*(func_args[it%2]))
        it += 1

    if pool_alloc:
        pipe_lib.pool_destroy()

    return

def calc_norm(U_, app_data):
    N = app_data['N']

    grid_data = app_data['grid_data']
    F_ = grid_data['F_']
    U_EXACT_ = grid_data['U_EXACT_']

    # lib function name
    norm = app_data['pipeline_norm']

    resid = np.zeros((1), np.float64)
    err   = np.zeros((1), np.float64)

    # lib function args
    norm_args = []
    norm_args += [ctypes.c_int(N)]
    norm_args += [ctypes.c_void_p(F_.ctypes.data)]
    norm_args += [ctypes.c_void_p(U_EXACT_.ctypes.data)]
    norm_args += [ctypes.c_void_p(U_.ctypes.data)]
    norm_args += [ctypes.c_void_p(err.ctypes.data)]
    norm_args += [ctypes.c_void_p(resid.ctypes.data)]

    # call lib function
    norm(*norm_args)

    # save the old norm values
    app_data['old_residual'] = app_data['resid']
    app_data['old_err'] = app_data['err']

    # register the norm values in the data dictionary
    app_data['resid'] = resid[0]
    app_data['err'] = err[0]

    return

def call_mg_cycle(U_, W_, app_data):
    n = app_data['n']

    grid_data = app_data['grid_data']
    F_ = grid_data['F_']

    # lib function name
    func_name = 'pipeline_'+app_data['cycle_name']
    mg_cycle_func = app_data[func_name]

    # lib function args
    mg_cycle_args = []
    mg_cycle_args += [ctypes.c_int(n)]
    mg_cycle_args += [ctypes.c_void_p(F_.ctypes.data)]
    mg_cycle_args += [ctypes.c_void_p(U_.ctypes.data)]
    mg_cycle_args += [ctypes.c_void_p(W_.ctypes.data)]

    # call lib function
    mg_cycle_func(*mg_cycle_args)

    return

def multigrid(app_data):
    grid_data = app_data['grid_data']
    U_ = grid_data['U_']
    W_ = grid_data['W_']
    UW = [U_, W_]

    nit = app_data['nit']
    time_store = {}
    time_taken = 0

    print_layout(app_data)
    print_errors(0, app_data)

    timer = app_data['timer']
    nruns = int(app_data['runs'])

    run = 0
    while run < nruns:
        it = 0
        if timer:
            t1 = time.time()
        while it < nit:
            it += 1
            call_mg_cycle(UW[(it-1)%2], UW[it%2], app_data)
            if not timer:
                calc_norm(UW[it%2], app_data)
                print_errors(it, app_data)
        if timer:
            t2 = time.time()
            time_store[run] = float(t2) - float(t1)
            #print("Time taken for iter ", run," = ",time_store[run]*1000, "ms")
            time_taken += time_store[run]
        run += 1

    if timer:
        time_taken = time_taken / nruns
        print("")
        print("[exec_mg] : Average time taken to execute = ",
              time_taken*1000, " ms")

    return

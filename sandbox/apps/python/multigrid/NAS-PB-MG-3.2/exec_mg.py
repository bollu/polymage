import sys
import os
import ctypes
import numpy as np
import time

from verify import verify_norm

def calc_norm(V_, app_data):
    N = V_.shape[0]

    # lib function name
    norm2u3 = app_data['norm2u3']

    rnm2 = np.zeros((1), np.float64)
    rnmu = np.zeros((1), np.float64)

    # lib function args
    norm_args = []
    norm_args += [ctypes.c_int(N)]
    norm_args += [ctypes.c_void_p(rnm2.ctypes.data)]
    norm_args += [ctypes.c_void_p(rnmu.ctypes.data)]
    norm_args += [ctypes.c_void_p(V_.ctypes.data)]

    # call lib function
    norm2u3(*norm_args)

    # register the norm values in the data dictionary
    app_data['rnm2'] = rnm2[0]
    app_data['rnmu'] = rnmu[0]

    return

def calc_resid(U_, V_, R_):
    N = r.shape[0]

    # lib function name
    resid = app_data['pipeline_resid']

    # lib function args
    resid_args = []
    resid_args += [ctypes.c_int(N)]
    resid_args += [ctypes.c_void_p(U_.ctypes.data)]
    resid_args += [ctypes.c_void_p(V_.ctypes.data)]
    resid_args += [ctypes.c_void_p(R_.ctypes.data)]

    # call the function
    resid(*resid_args)

    return

def call_mg3P(U_IN_, R_IN_, V_, U_OUT_, R_OUT_):
    N = r.shape[0]

    # lib function name
    mg = app_data['pipeline_nas_mg']

    # lib function args
    mg_args = []
    mg_args += [ctypes.c_int(N)]
    mg_args += [ctypes.c_void_p(R_IN_.ctypes.data)]
    mg_args += [ctypes.c_void_p(U_IN_.ctypes.data)]
    mg_args += [ctypes.c_void_p(V_.ctypes.data)]
    mg_args += [ctypes.c_void_p(R_OUT_.ctypes.data)]
    mg_args += [ctypes.c_void_p(U_OUT_.ctypes.data)]

    # call the function
    mg(*mg_args)

    return

def multigrid(app_data):
    grid_data = app_data['grid_data']

    U0_ = grid_data['U0_']
    R0_ = grid_data['R0_']
    V_ = grid_data['V_']
    U1_ = grid_data['U1_']
    R1_ = grid_data['R1_']

    UU = [U0_, U1_]
    RR = [R0_, R1_]

    lib_mg = app_data['nas_mg.so']

    # compute the initial residual
    print("[exec]: computing the initial residual ...")
    calc_resid(U0_, V_, R0_, app_data)
    print("[exec]: ... DONE")

    # calculate the initial residual norm
    print("[exec]: calculating the initial norm ...")
    calc_norm(R0_, app_data)
    print("[exec]: ... DONE")

    print()
    print("[exec]: norm =", app_data['rnm2'])
    print("        err  =", app_data['rnmu'])

    app_args = app_data['app_args']
    pool_alloc = bool(app_args.pool_alloc)

    if pool_alloc:
        lib_mg.pool_init()

    # mg3p time
    t_total = 0.0

    nit = app_data['nit']
    nit = 10
    # call 'nit' v-cycles
    for it in range(1, nit+1):
        print()
        print("[exec]: iter", it)
        in_ = (it-1)%2
        out_ = it%2

        # timer ON
        if timer:
            t1 = time.time()
        
        # call the v-cycle
        call_mg3P(UU[in_], RR[in_], V_, UU[out_], RR[out_])

        # timer OFF
        if timer:
            t2 = time.time()
            # mg3P time +=
            t_total += (float(t2) - float(t1))

        # update the residual
        calc_resid(UU[out_], V_, RR[out_], app_data)

    #endfor

    if pool_alloc:
        lib_mg.pool_destroy()

    print()
    print("[exec]: calculating the final norm ...")
    # calculate the final residual norm
    calc_norm(RR[nit%2], app_data)
    print("[exec]: ... DONE")

    print()
    print("[exec]: time taken by multigrid pipeline =", t_total, "s")

    return

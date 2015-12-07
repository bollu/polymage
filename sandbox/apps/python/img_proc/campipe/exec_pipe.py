import sys
import os
import ctypes
import numpy as np
import time

from printer import print_line

from compiler   import *
from constructs import *
from utils import *

def call_pipe(U_, W_, app_data):
    rows = app_data['rows']
    cols = app_data['cols']

    img_data = app_data['img_data']
    IN = img_data['IN']
    M3200 = img_data['M3200']
    M7000 = img_data['M7000']
    OUT = img_data['OUT']

    # lib function name
    func_name = 'pipeline_'+app_data['app_name']
    pipe_func = app_data[func_name]

    # lib function args
    pipe_args = []
    pipe_args += [ctypes.c_int(rows)]
    pipe_args += [ctypes.c_int(cols)]
    pipe_args += [ctypes.c_void_p(IN.ctypes.data)]
    pipe_args += [ctypes.c_void_p(M3200.ctypes.data)]
    pipe_args += [ctypes.c_void_p(M7000.ctypes.data)]
    pipe_args += [ctypes.c_void_p(OUT.ctypes.data)]

    # call lib function
    pipe_func(*pipe_args)

    return

def campipe(app_data):
    grid_data = app_data['grid_data']
    U_ = grid_data['U_']
    W_ = grid_data['W_']

    runs = app_data['runs']
    it  = 0

    timer = app_data['timer']
    if timer == True:
        t1 = time.time()

    while it < runs :
        call_pipe(U_, W_, app_data)
        if timer == False:
            calcNorm(W_, app_data)

    if timer == True:
        t2 = time.time()

        time_taken = float(t2) - float(t1)
        print("")
        print("[exec_pipe] : time taken to execute = ", time_taken, " ms")

    return

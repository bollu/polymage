import sys
import os
import ctypes
import numpy as np
import time

from printer import print_line

from compiler   import *
from constructs import *
from utils import *

def call_pipe(app_data):
    rows = app_data['rows']
    cols = app_data['cols']

    app_args = app_data['app_args']
    colour_temp = float(app_args.colour_temp)
    contrast = float(app_args.contrast)
    gamma = float(app_args.gamma)

    img_data = app_data['img_data']
    IN = img_data['IN']
    M3200 = img_data['M3200']
    M7000 = img_data['M7000']
    OUT = img_data['OUT']

    # lib function name
    func_name = 'pipeline_'+app_data['app']
    pipe_func = app_data[func_name]

    # lib function args
    pipe_args = []
    pipe_args += [ctypes.c_int(cols)]
    pipe_args += [ctypes.c_int(rows)]
    pipe_args += [ctypes.c_float(colour_temp)]
    pipe_args += [ctypes.c_float(contrast)]
    pipe_args += [ctypes.c_float(gamma)]
    pipe_args += [ctypes.c_void_p(IN.ctypes.data)]
    pipe_args += [ctypes.c_void_p(M3200.ctypes.data)]
    pipe_args += [ctypes.c_void_p(M7000.ctypes.data)]
    pipe_args += [ctypes.c_void_p(OUT.ctypes.data)]

    # call lib function
    pipe_func(*pipe_args)

    return

def campipe(app_data):
    app_args = app_data['app_args']
    runs = int(app_args.runs)
    it  = 0

    timer = bool(app_args.timer)
    if timer == True:
        t1 = time.time()

    while it < runs :
        call_pipe(app_data)
        it += 1

    if timer == True:
        t2 = time.time()

        time_taken = float(t2) - float(t1)
        print("")
        print("[exec_pipe] : time taken to execute = ", time_taken*1000, " ms")

    return

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
    rows = app_data['R']
    cols = app_data['C']

    app_args = app_data['app_args']
    patch_size = int(app_args.patch_size)
    search_area = int(app_args.search_area)

    img_data = app_data['img_data']
    INL = img_data['INL']
    INR = img_data['INR']
    OUT = img_data['OUT']

    # lib function name
    func_name = 'pipeline_'+app_data['app']
    pipe_func = app_data[func_name]

    # lib function args
    pipe_args = []
    pipe_args += [ctypes.c_int(cols)]
    pipe_args += [ctypes.c_int(rows)]
    pipe_args += [ctypes.c_void_p(INL.ctypes.data)]
    pipe_args += [ctypes.c_void_p(INR.ctypes.data)]
    pipe_args += [ctypes.c_void_p(OUT.ctypes.data)]

    # call lib function
    pipe_func(*pipe_args)
    
    return

def lensblur(app_data):

    it  = 0
    app_args = app_data['app_args']
   
    runs = int(app_args.runs)
    timer = app_args.timer
    if timer == True:
        t1 = time.time()

    while it < runs :
        call_pipe(app_data)
        it += 1

    if timer == True:
        t2 = time.time()

        time_taken = float(t2) - float(t1)
        print("")
        print("[exec_pipe] : time taken to execute = ",
              (time_taken * 1000) / runs, " ms")

    return


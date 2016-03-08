import numpy as np
import time
import sys

from init import init_all, init_norm
from printer import print_header, print_config, print_line
from builder import create_lib, build_mGCycle
from exec_mg import multigrid
from app_tuner import auto_tune

app = "jacobi-2d"

def main():
    #-------------------------------------------------------------------
    print
    print_line()
    print_header()
    print_line()
    #-------------------------------------------------------------------
    app_data = {}
    pipe_data = {}

    print("[main]: initializing...")
    print("")

    # init all the required data
    init_all(pipe_data, app_data)

    print_config(app_data)
    cycle_name = app_data['cycle']+"cycle"
    if app_data['mode'] == 'tune':
        auto_tune(pipe_data, app_data)
    else:
        #-------------------------------------------------------------------
        create_lib(        None,    "norm", pipe_data, app_data, app_data['mode'])
        create_lib(build_mGCycle, cycle_name, pipe_data, app_data, app_data['mode'])
        #-------------------------------------------------------------------
        init_norm(app_data)
        multigrid(app_data)
        #-------------------------------------------------------------------

    return

main()

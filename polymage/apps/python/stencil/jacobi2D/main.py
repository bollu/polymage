import numpy as np
import time
import sys

from __init__ import *

from init import init_all, init_norm
from printer import print_header, print_config, print_line
from builder import create_lib, build_jacobi
from exec_pipe import exec_jacobi
from app_tuner import auto_tune

app = "jacobi"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    # init all the required data
    init_all(app_data)
    print_config(app_data)

    if app_data['mode'] == 'tune':
        print("Tuning...")
        auto_tune(app_data)
    else:
        #-------------------------------------------------------------------
        create_lib(None, "norm", app_data)
        create_lib(build_jacobi, app, app_data)
        #-------------------------------------------------------------------
        init_norm(app_data)
        exec_jacobi(app_data)
        #-------------------------------------------------------------------

    return

main()

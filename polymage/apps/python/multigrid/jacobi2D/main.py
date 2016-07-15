import numpy as np
import time
import sys

from __init__ import *

from init import init_all, init_norm
from printer import print_header, print_config, print_line
from builder import create_lib, build_mg_cycle
from exec_mg import multigrid
from app_tuner import auto_tune

app = "jacobi-2d"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    # init all the required data
    init_all(app_data)

    print_config(app_data)
    if app_data['mode'] == 'tune':
        auto_tune(app_data)
    else:
        #-------------------------------------------------------------------
        create_lib(None, "norm", app_data)
        create_lib(build_mg_cycle, app_data['cycle_name'], app_data)
        #-------------------------------------------------------------------
        init_norm(app_data)
        multigrid(app_data)
        #-------------------------------------------------------------------

    return

main()

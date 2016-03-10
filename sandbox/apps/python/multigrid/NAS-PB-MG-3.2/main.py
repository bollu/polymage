import numpy as np
import time
import sys

from init import init_all
from verify import verify_norm
from builder import create_lib, build_resid, build_mg_3p
from exec_mg import multigrid
from printer import print_line, print_header, print_config

app = 'nas-pb-mg-3.2'
prob_classes = ['S', 'W', 'A', 'B', 'C', 'D']

def main():
    print_header()

    app_data = {}

    # init all the required data
    init_all(app_data)

    print_config(app_data)

    app_name = "nas_mg_class_"+app_data['class']
    app_data['app'] = app_name

    if mode == 'tune':
        #app_tune(app_data)
        pass
    else:
        #-------------------------------------------------------------------
        # setting up multigrid v-cycle computation
        create_pipe_lib(build_mg3p, pipe_name, app_data)
        # setting up standalone version of residual computation
        create_pipe_lib(build_resid, "resid", app_data)
        #-------------------------------------------------------------------
        multigrid(app_data)
        verify_norm(app_data)
        #-------------------------------------------------------------------

    return

main()

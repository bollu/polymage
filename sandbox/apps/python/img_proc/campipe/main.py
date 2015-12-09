import numpy as np
import time
import sys

sys.path.insert(0, '../../../../')

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib, build_campipe
from exec_pipe import campipe
#from app_tuner import auto_tune

app = "camera_pipeline"

def main():
    print_header()

    print("[main]: initializing...")
    print("")

    app_data = {}
    pipe_data = {}

    app_data['app'] = app

    # init all the required data
    init_all(sys.argv, pipe_data, app_data)

    print_config(app_data)
    if app_data['mode'] == 'tune':
        auto_tune(pipe_data, app_data)
    else:
        # create shared lib
        create_lib(build_campipe, app, pipe_data, app_data, app_data['mode'])
        # execute the compiled pipeline
        campipe(app_data)

    return

main()

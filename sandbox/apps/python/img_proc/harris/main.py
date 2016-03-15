import numpy as np
import time
import sys

sys.path.insert(0, '../../../../')

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib,build_harris
from exec_pipe import harrispipe
from app_tuner import auto_tune

app = "harris"

def main():
    print_header()
    
    print("[main]: initializing...")
    print("")

    app_data = {}
    pipe_data = {}

    app_data['app'] = app
    app_data['app_name'] = app

    init_all(sys.argv, pipe_data, app_data)
    print_config(app_data)
    if app_data['mode'] == 'tune':
        print("Tuning")
        auto_tune(pipe_data,app_data)
    else:
        create_lib(build_harris, app, app_data)
        harrispipe(app_data)

    return

main()

import numpy as np
import time
import sys

from __init__ import *

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib, build_campipe
from exec_pipe import campipe
from app_tuner import auto_tune

app = "camera_pipeline"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    init_all(app_data)
    print_config(app_data)

    if app_data['mode'] == 'tune':
        auto_tune(app_data)
    else:
        # create shared lib
        create_lib(build_campipe, app, app_data)
        # execute the compiled pipeline
        campipe(app_data)

    return

main()

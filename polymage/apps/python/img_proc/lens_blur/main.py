import numpy as np
import time
import sys

from __init__ import *

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib,build_lensblur
from exec_pipe import lensblur
#from app_tuner import auto_tune

app = "lens_blur"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    init_all(app_data)
    print_config(app_data)

    if app_data['mode'] == 'tune':
        #auto_tune(pipe_data,app_data)
        pass
    else:
        create_lib(build_lensblur, app, app_data)
        lensblur(app_data)

    return

main()

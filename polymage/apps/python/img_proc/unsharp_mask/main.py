import numpy as np
import time
import sys

from __init__ import *

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib,build_unsharp
from exec_pipe import unsharp_mask
from app_tuner import auto_tune

app = "unsharp_mask"

def main():
    print_header()
    
    print("[main]: initializing...")
    print("")

    app_data = {}

    app_data['app'] = app
    app_data['app_name'] = app
    app_data['ROOT'] = ROOT

    init_all(app_data)
    print_config(app_data)
    if app_data['mode'] == 'tune':
        print("Tuning")
        auto_tune(app_data)
    else:
        create_lib(build_unsharp, app, app_data)
        unsharp_mask(app_data)

    return

main()

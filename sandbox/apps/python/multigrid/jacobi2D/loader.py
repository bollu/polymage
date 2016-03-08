import ctypes
import numpy as np

def load_lib(lib_file, lib_func_name, app_data):
    print("")
    print("[loader]: loading the shared library", lib_file, "...")

    # assuming that it is present,
    # load the shared library
    lib = ctypes.cdll.LoadLibrary('./'+lib_file)

    print("[loader]: ... DONE")

    # name of the lib function
    lib_func = lib[lib_func_name]

    # register the library and the function in the data dictionary
    app_data[str(lib_file)] = lib
    app_data[str(lib_func_name)] = lib_func

    return

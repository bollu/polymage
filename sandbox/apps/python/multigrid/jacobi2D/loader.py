import ctypes
import numpy as np

def loadLib(libFile, libFuncName, appData):
    print("")
    print("[loader]: loading the shared library", libFile, "...")

    # assuming that it is present,
    # load the shared library
    lib = ctypes.cdll.LoadLibrary('./'+libFile)

    print("[loader]: ... DONE")

    # name of the lib function
    libFunc = lib[libFuncName]

    # register the library and the function in the data dictionary
    appData[str(libFile)] = lib
    appData[str(libFuncName)] = libFunc

    return

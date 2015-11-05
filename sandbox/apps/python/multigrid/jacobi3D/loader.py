import ctypes
import numpy as np

def loadLib(libFile, libFuncName, appData):
    print("")
    print("[misc]: loading the shared library", libFile, "...")

    # assuming that it is present,
    # load the shared library
    lib = ctypes.cdll.LoadLibrary('./'+libFile)

    print("[misc]: ... DONE")

    # name of the lib function
    libFunc = lib[libFuncName]

    # register the library and the function in the data dictionary
    appData[str(libFile)] = lib
    appData[str(libFuncName)] = libFunc

    return

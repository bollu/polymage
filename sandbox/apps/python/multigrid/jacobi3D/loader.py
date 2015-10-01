import ctypes
import numpy as np

def loadLib(libFile, libFuncName, dataDict):
    print("")
    print("[misc]: loading the shared library", libFile, "...")

    # assuming that it is present,
    # load the shared library
    lib = ctypes.cdll.LoadLibrary('./'+libFile)

    print("[misc]: ... DONE")

    # name of the lib function
    libFunc = lib[libFuncName]

    # register the library and the function in the data dictionary
    dataDict[str(libFile)] = lib
    dataDict[str(libFuncName)] = libFunc

    return

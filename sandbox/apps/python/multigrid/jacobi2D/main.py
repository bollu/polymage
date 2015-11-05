import numpy as np
import time
import sys

from init      import initAll, initNorm
from printer   import printHeader, printConfig, printLine
from builder   import createLib, buildMGCycle
from execMG    import multigrid
from app_tuner import auto_tune

app = "jacobi-2d"

#-------------------------------------------------------------------
# initialize parameters

#-------------------------------------------------------------------
# main

def main():
    #-------------------------------------------------------------------
    print
    printLine()
    printHeader()
    printLine()
    #-------------------------------------------------------------------
    appData = {}
    pipeData = {}

    print("[main]: initializing...")
    print("")

    # init all the required data
    initAll(pipeData, appData)

    printConfig(appData)
    cycleName = appData['cycle']+"cycle"
    if appData['mode'] == 'tune':
        auto_tune(pipeData, appData)
    else:
        #-------------------------------------------------------------------
        createLib(        None,    "norm", pipeData, appData, appData['mode'])
        createLib(buildMGCycle, cycleName, pipeData, appData, appData['mode'])
        #-------------------------------------------------------------------
        initNorm(appData)
        multigrid(appData)
        #-------------------------------------------------------------------

    return

main()

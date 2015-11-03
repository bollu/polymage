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
    impipeData = {}

    print("[main]: initializing...")
    print("")

    # init all the required data
    initAll(impipeData, appData)

    printConfig(appData)
    cycleName = appData['cycle']+"cycle"
    if appData['mode'] == 'tune':
        auto_tune(impipeData, appData)
    else:
        #-------------------------------------------------------------------
        createLib(        None,    "norm", impipeData, appData, appData['mode'])
        createLib(buildMGCycle, cycleName, impipeData, appData, appData['mode'])
        #-------------------------------------------------------------------
        initNorm(appData)
        multigrid(appData)
        #-------------------------------------------------------------------

    return

main()

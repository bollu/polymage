import numpy as np
import time
import sys

from init    import initAll, initNorm
from printer import printHeader, printConfig, printLine
from builder import createLib, buildMGCycle
from execMG  import multigrid

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
    dataDict = {}
    impipeDict = {}

    print("[main]: initializing...")
    print("")

    # init all the required data
    initAll(impipeDict, dataDict)

    printConfig(dataDict)
    cycleName = dataDict['cycle']+"cycle"
    if dataDict['mode'] == tune:
        app_tuner(impipeDict, dataDict)
    else:
        #-------------------------------------------------------------------
        createLib(        None,    "norm", impipeDict, dataDict, dataDict['mode'])
        createLib(buildMGCycle, cycleName, impipeDict, dataDict, dataDict['mode'])
        #-------------------------------------------------------------------
        initNorm(dataDict)
        multigrid(dataDict)
        #-------------------------------------------------------------------

    return

main()

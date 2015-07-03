import numpy as np
import time
import sys

from init    import initAll, initNorm
from printer import printHeader, printConfig, printLine
from builder import createLib, buildVCycle
from execMG  import multigrid

app = "jacobi-3d"

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

    print "[main]: initializing..."
    print

    # init all the required data
    initAll(impipeDict, dataDict)

    printConfig(dataDict)
    #-------------------------------------------------------------------
    createLib(       None,   "norm", impipeDict, dataDict, dataDict['mode'])
    createLib(buildVCycle, "vcycle", impipeDict, dataDict, dataDict['mode'])
    #-------------------------------------------------------------------
    initNorm(dataDict)
    multigrid(dataDict)

    #-------------------------------------------------------------------

    return

main()

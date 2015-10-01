import sys

def printLine(toFile=None):
    if toFile:
        print("--------------------------------------------------", file=to_file)
    else:
        print("--------------------------------------------------")

def printHeader():
    print("              Multigrid Jacobi-3D")

def printUsage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <#iters>")
    print("[main]: 'mode'  :: {'new', 'existing', 'tune'}")

def printConfig(dataDict):
    nx = ny = nz = dataDict['N']
    printLine()
    print("# Problem Settings #")
    print("")
    print("[main]: multigrid levels =", dataDict['L'])
    print("[main]: grid size        =", nx, "x", ny, "x", nz)
    print("[main]: n-iterations     =", dataDict['nit'])

    print("[main]: nu1              =", dataDict['nu1'])
    print("[main]: nu2              =", dataDict['nu2'])
    print("[main]: nuc              =", dataDict['nuc'])
    printLine()

def printLayout(dataDict):
    L = dataDict['L']
    n = dataDict['n']
    nu1 = dataDict['nu1']
    nu2 = dataDict['nu2']
    nuc = dataDict['nuc']

    print(": levels: 0.."+str(L)+" , coarse grid: "+str(n)+"x"+str(n)+"x"+str(n))
    print("# pre-smoothing: ",nu1,", post-smoothing: ",nu2,", coarse relaxation: ",nuc)
    print("")
    print("# discr                 iter       error         residual       rho-error    rho-residual")
    print("# ---------------       ----       -----         --------       ---------    ------------")


def printErrors(it, dataDict):
    N = dataDict['N']
    oldResidual = dataDict['oldResidual']
    oldErr      = dataDict['oldErr']
    resid       = dataDict['resid']
    err         = dataDict['err']

    if it == 0:
        print("")
        print("# "+"%4d"%N+"x"+"%d"%N+"x"+"%-4d"%N+"       :  "+ \
              "%3d" % it+"  :  "+ \
              "%0.6e" % err+"   "+"%0.6e" % resid)
    else:
        if oldResidual == 0.0:
            rhoResidual = 1.0
        else:
            rhoResidual =  (resid/oldResidual)

        if oldErr == 0.0:
            rhoErr   = 1.0
        else:
            rhoErr = (err/oldErr)

        print("# "+"%4d"%N+"x"+"%d"%N+"x"+"%-4d"%N+"       :  "+ \
              "%3d" % it+"  :  "+ \
              "%0.6e" % err+"   "+"%0.6e" % resid+"  :  "+ \
              "%0.6f" % rhoErr+"     "+"%0.6f" % rhoResidual)

    return


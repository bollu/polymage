import sys

def printLine(toFile=None):
    if toFile:
        print("--------------------------------------------------", file=to_file)
    else:
        print("--------------------------------------------------")

def printHeader():
    print("              Multigrid Jacobi-2D")

def printUsage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <#iters>")
    print("[main]: 'mode'  :: {'new', 'existing', 'tune'}")

def printConfig(appData):
    nx = ny = appData['N']
    printLine()
    print("# Problem Settings #")
    print("")
    print("[main]: multigrid levels =", appData['L'])
    print("[main]: grid size        =", nx, "x", ny)
    print("[main]: n-iterations     =", appData['nit'])

    print("[main]: nu1              =", appData['nu1'])
    print("[main]: nu2              =", appData['nu2'])
    print("[main]: nuc              =", appData['nuc'])
    printLine()

def printLayout(appData):
    L = appData['L']
    n = appData['n']
    nu1 = appData['nu1']
    nu2 = appData['nu2']
    nuc = appData['nuc']

    print(": levels: 0.."+str(L)+" , coarse grid: "+str(n)+"x"+str(n))
    print("# pre-smoothing: ",nu1,", post-smoothing: ",nu2,", coarse relaxation: ",nuc)
    print("")
    print("# discr            iter       error         residual       rho-error    rho-residual")
    print("# ----------       ----       -----         --------       ---------    ------------")


def printErrors(it, appData):
    N = appData['N']
    oldResidual = appData['oldResidual']
    oldErr      = appData['oldErr']
    resid       = appData['resid']
    err         = appData['err']

    if it == 0:
        print("")
        print("# "+"%4d"%N+"x"+"%d"%N+"       :  "+ \
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

        print("# "+"%4d"%N+"x"+"%d"%N+"       :  "+ \
              "%3d" % it+"  :  "+ \
              "%0.6e" % err+"   "+"%0.6e" % resid+"  :  "+ \
              "%0.6f" % rhoErr+"     "+"%0.6f" % rhoResidual)

    return


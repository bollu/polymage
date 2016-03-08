import sys

def print_line(to_file=None):
    if to_file:
        print("--------------------------------------------------",
              file=to_file)
    else:
        print("--------------------------------------------------")

def print_header():
    print("              Multigrid Jacobi-2D")

def print_usage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <#iters>")
    print("[main]: 'mode'  :: {'new', 'existing', 'ready', 'tune'}")

def print_config(app_data):
    nx = ny = app_data['N']
    print_line()
    print("# Problem Settings #")
    print("")
    print("[main]: multigrid levels =", app_data['L'])
    print("[main]: grid size        =", nx, "x", ny)
    print("[main]: n-iterations     =", app_data['nit'])

    print("[main]: nu1              =", app_data['nu1'])
    print("[main]: nu2              =", app_data['nu2'])
    print("[main]: nuc              =", app_data['nuc'])
    print_line()

def print_layout(app_data):
    L = app_data['L']
    n = app_data['n']
    nu1 = app_data['nu1']
    nu2 = app_data['nu2']
    nuc = app_data['nuc']

    print(": levels: 0.."+str(L)+" , coarse grid: "+str(n)+"x"+str(n))
    print("# pre-smoothing: ", nu1,", post-smoothing: ", nu2,
          ", coarse relaxation: ", nuc)
    print("")
    print("# discr            iter       error         residual       rho-error    rho-residual")
    print("# ----------       ----       -----         --------       ---------    ------------")


def print_errors(it, app_data):
    N = app_data['N']
    old_residual = app_data['old_residual']
    old_err = app_data['old_err']
    resid = app_data['resid']
    err = app_data['err']

    if it == 0:
        print("")
        print("# "+"%4d"%N+"x"+"%d"%N+"       :  "+ \
              "%3d" % it+"  :  "+ \
              "%0.6e" % err+"   "+"%0.6e" % resid)
    else:
        if old_residual == 0.0:
            rho_residual = 1.0
        else:
            rho_residual =  (resid/old_residual)

        if old_err == 0.0:
            rho_err   = 1.0
        else:
            rho_err = (err/old_err)

        print("# "+"%4d"%N+"x"+"%d"%N+"       :  "+ \
              "%3d" % it+"  :  "+ \
              "%0.6e" % err+"   "+"%0.6e" % resid+"  :  "+ \
              "%0.6f" % rho_err+"     "+"%0.6f" % rho_residual)

    return


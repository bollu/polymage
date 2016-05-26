import sys

def print_line(to_file=None):
    if to_file:
        print("--------------------------------------------------",
              file=to_file)
    else:
        print("--------------------------------------------------")
    return

def print_header():
    print_line()
    print("              Jacobi-2D")
    print_line()
    print("[main]: initializing...")
    print("")
    return

def print_usage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <flags>")
    print("[main]: 'mode'  :: {'new', 'existing', 'ready', 'tune'}")
    return

def print_config(app_data):
    nx = ny = nz = app_data['N']
    print_line()
    print("# Problem Settings #")
    print("")
    print("[main]: grid size          =", ny, "x", nx)
    print("[main]: stencil-iterations =", app_data['T'])

    print_line()
    return

def print_layout(app_data):
    N = app_data['N']
    T = app_data['T']

    print(": levels: 1.."+str(1)+" , coarse grid: "+str(n)+"x"+str(n))
    print("# pre-smoothing: ", T,", post-smoothing: ", 0,
          ", coarse relaxation: ", 0)
    print("")
    print("# discr                 iter       residual norm")
    print("# ---------------       ----       -------------")
    return

def print_norm(it, app_data):
    N = app_data['N']
    resid = app_data['resid']
    err = app_data['err']

    print("")
    print("# "+"%4d"%N+"x"+"%d"%N+"       :  "+ \
          "%3d" % it+"  :  "+ \
          "%0.6e" % err+"   "+"%0.6e" % resid)

    return
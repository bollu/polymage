import sys

def print_line(to_file=None):
    if to_file:
        print("--------------------------------------------------",
              file=to_file)
    else:
        print("--------------------------------------------------")

def print_header():
    print_line()
    print("NAS Parallel Benchmark v3.2")
    print("            MG")
    print("[main]: initializing...")
    print_line()
    print()

def print_usage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <class> <mode>")
    print("[main]: 'class' :: {'S', 'W', 'A', 'B', 'C', 'D'}")
    print("[main]: 'mode'  :: {'new', 'existing', 'ready', 'tune'}")

def print_config(app_data):
    nx = ny = nz = app_data['prob_size']
    print_line()

    print("# Problem Settings #")
    print("[main]: CLASS        = \""+app_data['prob_class']+"\"")
    print("[main]: top level    =", app_data['lt'])
    print("[main]: bottom level =", app_data['lb'])
    print("[main]: grid size    =", nx, "x", ny, "x", nz)
    print("[main]: n-iterations =", app_data['nit'])

    print()
    print("# Stencil Co-efficients #")
    print("[main]: a =", app_data['a'])
    print("[main]: c =", app_data['c'])

    verify_data = app_data['verify_data']
    print()
    print("# Verification Values #")
    print("[main]: threshold         =", verify_data['epsilon'])
    print("[main]: Class \""+app_data['prob_class']+"\" " \
                 + "L2 norm =", verify_data['verify_value'])

    print()
    print("# Initial Norms #")
    print("[main]: initial norm =", app_data['rnm2'])
    print("[main]: initial err  =", app_data['rnmu'])

    print_line()

    return

import sys

def print_line(to_file=None):
    if to_file:
        print("--------------------------------------------------", file=to_file)
    else:
        print("--------------------------------------------------")

def print_header():
    print("Max Filter")

def print_usage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <image> ", end=" ")
    print("<rows> <cols> <patch_size> <search_area> <nruns>")
    print("[main]: <mode>  :: {'new', 'existing', 'tune'}")

def print_config(app_data):
    app_args = app_data['app_args']
    rows = app_data['R']
    cols = app_data['C']
    print_line()
    print("# Problem Settings #")
    print("")
    print("[main]: mode        =", app_args.mode)
    print("[main]: image       =", app_args.img_file)
    print("[main]: image size  =", rows, "x", cols)
    print("[main]: patch size  =", app_args.patch_size)
    print("[main]: search area =", app_args.search_area)
    print("[main]: nruns       =", app_args.runs)
    print_line()

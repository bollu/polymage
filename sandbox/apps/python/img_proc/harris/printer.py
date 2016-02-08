import sys

def print_line(to_file=None):
    if to_file:
        print("--------------------------------------------------", file=to_file)
    else:
        print("--------------------------------------------------")

def print_header():
    print("Harris Corner")

def print_usage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <image> ", end=" ")
    print("<rows> <cols> <nruns>")
    print("[main]: <mode>  :: {'new', 'existing', 'tune'}")

def print_config(app_data):
    app_args = app_data['app_args']
    rows = app_data['rows']
    cols = app_data['cols']
    print_line()
    print("# Problem Settings #")
    print("")
    print("[main]: mode        =", app_args.mode)
    print("[main]: image       =", app_args.img_file)
#    print("[main]: image size  =", rows, "x", cols)
#    print("[main]: colour_temp =", app_args.colour_temp)
#    print("[main]: contrast    =", app_args.contrast)
#    print("[main]: gamma       =", app_args.gamma)
    print("[main]: nruns       =", app_args.runs)
    print_line()

import sys

def print_line(to_file=None):
    if to_file:
        print("--------------------------------------------------", file=to_file)
    else:
        print("--------------------------------------------------")

def print_header():
    print("               Camera Pipeline")

def print_usage():
    print("[main]: Usage: ")
    print("[main]: "+sys.argv[0]+" <mode> <image> <colour_temp>", end=" ")
    print("<contrast> <gamma> <rows> <cols> <nruns>")
    print("[main]: <mode>  :: {'new', 'existing', 'tune'}")

def print_config(app_data):
    rows = app_data['rows']
    cols = app_data['cols']
    print_line()
    print("# Problem Settings #")
    print("")
    print("[main]: image       =", app_data['image_path'])
    print("[main]: image size  =", nx, "x", ny)
    print("[main]: colour_temp =", app_data['colour_temp'])
    print("[main]: contrast    =", app_data['contrast'])
    print("[main]: gamma       =", app_data['gamma'])
    print("[main]: nruns       =", app_data['nruns'])
    print_line()

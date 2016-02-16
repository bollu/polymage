import sys
import subprocess

from compiler import *
from constructs import *

def c_compile(in_file, out_file, c_compiler=None):
    if c_compiler == None or c_compiler == "gnu":
        cxx = "g++"
        opt = "-fopenmp -march=native -O3 -ftree-vectorize"
    elif c_compiler == "intel":
        cxx = "icpc"
        opt = "-openmp -xhost -O3 -ipo -ansi-alias"

    #fi

    #inc = "-I../../../memory_allocation/ "+\
    #      "../../../memory_allocation/simple_pool_allocator.cpp"
    inc = ""

    shared = "-fPIC -shared"
    out = "-o "+out_file

    compile_str = cxx + " " \
                + opt + " " \
                + inc + " " \
                + shared + " " \
                + in_file + " " \
                + out

    print("")
    print("[cpp_compiler]: compiling", in_file, "to", out_file, "...")
    print(">", compile_str)
    subprocess.check_output(compile_str, shell=True)
    print("[cpp_compiler]: ... DONE")

    return

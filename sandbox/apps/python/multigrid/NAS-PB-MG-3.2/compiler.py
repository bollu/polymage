import sys
import subprocess

from compiler   import *
from constructs import *

def c_compile(in_file, out_file, arg_data):
    # CXX compiler and flags :
    cxx = arg_data.cxx
    cxx_flags = arg_data.cxx_flags
    #fi

    # Include Flags :
    #include = "-I../../../../memory_allocation/ "+\
    #          "../../../../memory_allocation/simple_pool_allocator.cpp"
    include = ""

    # Shared library Flags
    shared = "-fPIC -shared"
    out = "-o "+out_file

    compile_str = cxx + " " \
                + cxx_flags + " " \
                + include + " " \
                + shared + " " \
                + in_file + " " \
                + out

    print("")
    print("[cpp_compiler]: compiling", in_file, "to", out_file, "...")
    print(">", compile_str)
    subprocess.check_output(compile_str, shell=True)
    print("[cpp_compiler]: ... DONE")

    return

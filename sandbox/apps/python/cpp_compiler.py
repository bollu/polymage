import sys
import subprocess

from compiler import *
from constructs import *

def gen_compile_string(app_data):
    ROOT = app_data['ROOT']
    arg_data = app_data['app_args']
    # CXX compiler and flags :
    cxx = arg_data.cxx
    cxx_flags = arg_data.cxx_flags
    #fi

    # Include Flags :
    if bool(arg_data.pool_alloc):
        include = "-I"+ROOT+"/memory_allocation/ "+\
                  ROOT+"/memory_allocation/simple_pool_allocator.cpp"
    else:
        include = ""

    # Shared library Flags
    shared = "-fPIC -shared"

    compile_str = cxx + " " \
                + cxx_flags + " " \
                + include + " " \
                + shared + " " \

    app_data['cxx_string'] = compile_str

    return

def c_compile(in_file, out_file, app_data):
    gen_compile_string(app_data)

    compile_str = app_data['cxx_string']
    compile_str += " " + in_file + " -o " + out_file

    print("")
    print("[cpp_compiler]: compiling", in_file, "to", out_file, "...")
    print(">", compile_str)
    subprocess.check_output(compile_str, shell=True)
    print("[cpp_compiler]: ... DONE")

    return

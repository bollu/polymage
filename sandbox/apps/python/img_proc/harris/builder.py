import sys
import subprocess

from cpp_compiler import c_compile
from loader import load_lib
from polymage_harris import harris_pipe

from compiler import *
from constructs import *

def codegen(pipe, file_name, app_data):
    print("")
    print("[builder]: writing the code to", file_name, "...")

    code = pipe.generate_code(is_extern_c_func=True,
                              outputs_no_alloc=True,
                              are_io_void_ptrs=True)

    f = open(file_name, 'w')
    f.write(code.__str__())
    f.close()

    return

def build_harris(pipe_data, app_data):
    print("Inside build_harris function")
    
    out_harrispipe = harris_pipe(pipe_data)
    
    R = pipe_data['R']
    C = pipe_data['C']

    live_outs = [out_harrispipe]
    pipe_name = app_data['app']
    p_estimates = [(R, app_data['rows']), (C, app_data['cols'])]
    p_constraints = [ Condition(R, "==", app_data['rows']), \
                      Condition(C, "==", app_data['cols']) ]
    t_size = [16, 16, 16]
    g_size = 10
    opts = []
    if app_data['pool_alloc'] == True:
        opts += ['pool_alloc']

    pipe = buildPipeline(live_outs,
                         param_estimates=p_estimates,
                         param_constraints=p_constraints,
                         #tile_sizes = t_size,
                         group_size = g_size,
                         options = opts,
                         pipe_name = pipe_name)

    return pipe



def create_lib(build_func, pipe_name, impipe_data, app_data, mode):
    pipe_src  = pipe_name+".cpp"
    pipe_so   = pipe_name+".so"

    if build_func != None:
        if mode == 'new':
            # build the polymage pipeline
            pipe = build_func(impipe_data, app_data)

            # draw the pipeline graph to a png file
            #graph_gen(pipe, pipe_name, app_data)

            # generate pipeline cpp source
            codegen(pipe, pipe_src, app_data)

    if mode != 'ready':
        # compile the cpp code
        c_compile(pipe_src, pipe_so, c_compiler="gnu")

    # load the shared library
    pipe_func_name = "pipeline_"+pipe_name
    load_lib(pipe_so, pipe_func_name, app_data)

    return

import sys
import subprocess

from cpp_compiler import c_compile
from loader import load_lib
from polymage_campipe import cam_pipe

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

'''
def graph_gen(pipe, file_name, app_data):
    graph_file = file_name+".dot"
    png_graph = file_name+".png"

    print("")
    print("[builder]: writing the graph dot file to", graph_file, "...")

    #graph = pipe.drawPipelineGraphWithGroups()
    graph = pipe.originalGraph
    graph.write(graph_file)
    print("[builder]: ... DONE")

    dotty_str = "dot -Tpng "+graph_file+" -o "+png_graph

    print("")
    print("[builder]: drawing the graph using dotty to", png_graph)
    print(">", dotty_str)
    subprocess.check_output(dotty_str, shell=True)
    print("[builder]: ... DONE")

    return
'''

def build_campipe(pipe_data, app_data):
    # construct the camera pipeline
    out_campipe = cam_pipe(pipe_data, app_data)

    n = pipe_data['n']

    live_outs = [out_campipe]
    pipe_name = app_data['cycle_name']
    p_estimates = [(n, app_data['n'])]
    p_constraints = [ Condition(n, "==", app_data['n']) ]
    t_size = [16, 16]
    g_size = 10
    opts = []
    if app_data['pool_alloc'] == True:
        opts += ['pool_alloc']

    pipe = build_pipeline(live_outs,
                          param_estimates=p_estimates,
                          param_constraints=p_constraints,
                          tile_sizes = t_size,
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
            code_gen(pipe, pipe_src, app_data)

    if mode != 'ready':
        # compile the cpp code
        c_compile(pipe_src, pipe_so, c_compiler="gnu")

    # load the shared library
    pipe_func_name = "pipeline_"+pipe_name
    load_lib(pipe_so, pipe_func_name, app_data)

    return

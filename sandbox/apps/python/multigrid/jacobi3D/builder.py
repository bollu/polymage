from __init__ import *

import sys
import subprocess

sys.path.insert(0, ROOT+'apps/python/')

from cpp_compiler import c_compile
from loader import load_lib
from polymage_vcycle import v_cycle
from polymage_wcycle import w_cycle

from compiler import *
from constructs import *

def code_gen(pipe, file_name, app_data):
    print("")
    print("[builder]: writing the code to", file_name, "...")

    code = pipe.generate_code(is_extern_c_func=True,
                              are_io_void_ptrs=True)

    f = open(file_name, 'w')
    f.write(code.__str__())
    f.close()

    return

def generate_graph(pipe, file_name, app_data):
    graph_file = file_name+".dot"
    png_graph = file_name+".png"

    print("")
    print("[builder]: writing the graph dot file to", graph_file, "...")

    graph = pipe.pipeline_graph
    graph.write(graph_file)
    print("[builder]: ... DONE")

    dotty_str = "dot -Tpng "+graph_file+" -o "+png_graph

    print("")
    print("[builder]: drawing the graph using dotty to", png_graph)
    print(">", dotty_str)
    subprocess.check_output(dotty_str, shell=True)
    print("[builder]: ... DONE")

    return

def build_mg_cycle(app_data):
    pipe_data = app_data['pipe_data']
    cycle_type = app_data['cycle']

    if cycle_type == 'V':
        # construct the multigrid v-cycle pipeline
        mg = v_cycle(app_data)
    elif cycle_type == 'W':
        # construct the multigrid w-cycle pipeline
        mg = w_cycle(app_data)

    n = pipe_data['n']

    live_outs = [mg]
    pipe_name = app_data['cycle_name']
    p_estimates = [(n, app_data['n'])]
    p_constraints = [ Condition(n, "==", app_data['n']) ]
    t_size = [8, 8, 32]
    g_size = 6
    opts = []
    if app_data['early_free']:
        opts += ['early_free']
    if app_data['optimize_storage']:
        opts += ['optimize_storage']
    if app_data['pool_alloc']:
        opts += ['pool_alloc']

    mg_pipe = buildPipeline(live_outs,
                            param_estimates=p_estimates,
                            param_constraints=p_constraints,
                            tile_sizes = t_size,
                            group_size = g_size,
                            options = opts,
                            pipe_name = pipe_name)

    return mg_pipe

def create_lib(build_func, pipe_name, app_data):
    mode = app_data['mode']
    app_args = app_data['app_args']
    pipe_src = pipe_name+".cpp"
    pipe_so = pipe_name+".so"
    graph_gen = app_data['graph_gen']

    if build_func != None:
        if mode == 'new':
            # build the polymage pipeline
            pipe = build_func(app_data)

            # draw the pipeline graph to a png file
            if graph_gen:
                generate_graph(pipe, pipe_name, app_data)

            # generate pipeline cpp source
            code_gen(pipe, pipe_src, app_data)
        #fi
    #fi

    if mode != 'ready':
        # compile the cpp code
        c_compile(pipe_src, pipe_so, app_data)
    #fi

    # load the shared library
    lib_func_name = "pipeline_"+pipe_name
    load_lib(pipe_so, lib_func_name, app_data)

    return

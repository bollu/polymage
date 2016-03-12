import sys
import subprocess

sys.path.insert(0, '../../')

from cpp_compiler import c_compile
from loader import load_lib
from polymage_residual import resid_pipe
from polymage_mg3p import mg3p

from compiler   import *
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
    #graph = pipe.original_graph
    graph.write(graph_file)
    print("[builder]: ... DONE")

    dotty_str = "dot -Tpng "+graph_file+" -o "+png_graph

    print("")
    print("[builder]: drawing the graph using dotty to", png_graph)
    print(">", dotty_str)
    subprocess.check_output(dotty_str, shell=True)
    print("[builder]: ... DONE")

    return

def build_resid(app_data):
    pipe_data = app_data['pipe_data']

    # construct the residual pipeline on the finest grid
    r = resid_pipe(app_data)

    n = pipe_data['n']
    N = app_data['N']
    lt = app_data['lt']

    live_outs = [r]
    pipe_name = "resid"
    p_estimates = [(n, N[lt])]
    p_constraints = [ Condition(n, "==", N[lt]) ]
    t_size = [16, 16, 16]
    g_size = 1
    opts = []
    if app_data['pool_alloc']:
        opts += ['pool_alloc']


    r_pipe = buildPipeline(live_outs,
                           param_estimates=p_estimates,
                           param_constraints=p_constraints,
                           tile_sizes = t_size,
                           group_size = g_size,
                           options = opts,
                           pipe_name = pipe_name)

    return r_pipe

def build_mg3p(app_data):
    pipe_data = app_data['pipe_data']
    
    # construct the multigrid v-cycle pipeline
    mg_u, mg_r = mg3p(app_data)

    n = pipe_data['n']
    N = app_data['N']
    lt = app_data['lt']

    live_outs = [mg_u, mg_r]

    pipe_name = app_data['app']
    p_estimates = [(n, N[lt])]
    p_constraints = [ Condition(n, "==", N[lt]) ]
    t_size = [16, 16, 16]
    g_size = 1
    opts = []
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
        c_compile(pipe_src, pipe_so, app_args)
    #fi

    # load the shared library
    lib_func_name = "pipeline_"+pipe_name
    load_lib(pipe_so, lib_func_name, app_data)

    return

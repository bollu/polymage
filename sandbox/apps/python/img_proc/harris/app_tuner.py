from __init__ import *

import sys
sys.path.insert(0, ROOT+'apps/python/')

from cpp_compiler import *
from constructs import * #
from compiler import * #
import tuner
from polymage_harris import harris_pipe

def auto_tune(app_data):
    pipe_data = app_data['pipe_data']

    app_name = app_data['app']
    pipe_name = app_data['app']

    out_harrispipe = harris_pipe(pipe_data)
    live_outs = [out_harrispipe]
    
    R = pipe_data['R']
    C = pipe_data['C']

    rows = app_data['rows']-2
    cols = app_data['cols']-2

    param_estimates = [(R, rows), (C, cols)]
    param_constraints = [ Condition(R, "==", rows), \
                          Condition(C, "==", cols) ]
    
    dst_path = "/tmp"

    group_size_configs = [3, 5, 7, 9, 10]

    tile_size_configs = []

    tile_size_configs.append([64, 256])
    tile_size_configs.append([64, 128])

    tile_size_configs.append([32, 512])
    tile_size_configs.append([32, 256])
    tile_size_configs.append([32, 128])
    tile_size_configs.append([32, 64])

    tile_size_configs.append([16, 512])
    tile_size_configs.append([16, 256])
    tile_size_configs.append([16, 128])
    tile_size_configs.append([16, 64])

    tile_size_configs.append([8, 512])
    tile_size_configs.append([8, 256])
    tile_size_configs.append([8, 128])
    tile_size_configs.append([8, 64])
    tile_size_configs.append([8, 32])

    # relative path to root directory from app dir
    ROOT = app_data['ROOT']
    opts = []
    if app_data['early_free']:
        opts += ['early_free']
    if app_data['optimize_storage']:
        opts += ['optimize_storage']
    if app_data['pool_alloc']:
        opts += ['pool_alloc']

    gen_compile_string(app_data)
    cxx_string = app_data['cxx_string']

    # Generate Variants for Tuning
    # ============================

    gen_config = {"_tuner_app_name": app_name,
                  "_tuner_live_outs": live_outs,
                  "_tuner_param_constraints": param_constraints, #optional
                  "_tuner_param_estimates": param_estimates, #optional
                  "_tuner_tile_size_configs": tile_size_configs, #optional
                  "_tuner_group_size_configs": group_size_configs, #optional
                  "_tuner_opts": opts, #optional
                  "_tuner_dst_path" : dst_path, # optional
                  "_tuner_cxx_string" : cxx_string, # optional
                  "_tuner_root_path" : ROOT, # needed if pool_alloc is set
                  "_tuner_debug_flag": True, # optional
                  "_tuner_opt_datadict": app_data
                 }

    _tuner_src_path, _tuner_configs_count, _tuner_pipe = \
        tuner.generate(gen_config)

    '''
    _tuner_src_path = '/tmp/PolycNUUEYoLt2Mage'
    _tuner_configs_count = 75
    _tuner_pipe = buildPipeline(live_outs)
    '''

    img_data = app_data['img_data']
    IN = img_data['IN']
    OUT = img_data['OUT']

    pipe_args = {}
    pipe_args['R'] = rows
    pipe_args['C'] = cols
    pipe_args['img'] = IN 
    pipe_args['harris'] = OUT 


    # Execute the generated variants
    # ==============================

    exec_config = {"_tuner_app_name": app_name,
                   "_tuner_pipe": _tuner_pipe,
                   "_tuner_pipe_arg_data": pipe_args,
                   "_tuner_src_path": _tuner_src_path, # optional
                   "_tuner_configs_count": _tuner_configs_count, # optional
                   "_tuner_omp_threads": 4, # optional
                   "_tuner_nruns": 1, # optional
                   "_tuner_debug_flag": True, # optional
                   #"_tuner_custom_executor": minimal_exec_mg,
                   "_tuner_app_data": app_data
                  }

    tuner.execute(exec_config)

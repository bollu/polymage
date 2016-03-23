from __init__ import *

import sys
sys.path.insert(0, ROOT+"/apps/python/")

from cpp_compiler import *
from polymage_vcycle import v_cycle
from polymage_wcycle import w_cycle
from exec_mg import minimal_exec_mg
from constructs import *

from compiler import *
import tuner

def auto_tune(app_data):
    pipe_data = app_data['pipe_data']

    cycle_type = app_data['cycle']
    if cycle_type == 'V':
        mg = v_cycle(app_data)
    elif cycle_type == 'W':
        mg = w_cycle(app_data)

    app_name = app_data['cycle_name']
    live_outs = [mg]
    n = pipe_data['n']
    param_estimates = [(n, app_data['n'])]
    param_constraints = [ Condition(n, '==', app_data['n']) ]
    dst_path = "/tmp"

    group_size_configs = [2, 4, 6, 8]

    tile_size_configs = []
    tile_size_configs.append([8, 32])
    tile_size_configs.append([8, 64])

    tile_size_configs.append([8, 128])
    tile_size_configs.append([8, 256])
    tile_size_configs.append([8, 512])

    tile_size_configs.append([16, 64])
    tile_size_configs.append([16, 128])
    tile_size_configs.append([16, 256])
    tile_size_configs.append([16, 512])

    tile_size_configs.append([32, 64])
    tile_size_configs.append([32, 128])
    tile_size_configs.append([32, 256])
    tile_size_configs.append([32, 512])

    tile_size_configs.append([64, 128])
    tile_size_configs.append([64, 256])

    # relative path to root directory from app dir
    ROOT = app_data['ROOT']
    opts = []
    if app_data['early_free']:
        opts += ['early_free']
    if app_data['optimize_storage']:
        opts += ['optimize_storage']
    if app_data['pool_alloc']:
        opts += ['pool_alloc']
    if app_data['multipar']:
        opts += ['multipar']

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


    # Execute the generated variants
    # ==============================

    exec_config = {"_tuner_app_name": app_name,
                   "_tuner_pipe": _tuner_pipe,
                   "_tuner_src_path": _tuner_src_path, # optional
                   "_tuner_configs_count": _tuner_configs_count, # optional
                   "_tuner_omp_threads": 48, # optional
                   "_tuner_nruns": 1, # optional
                   "_tuner_debug_flag": True, # optional
                   "_tuner_custom_executor": minimal_exec_mg,
                   "_tuner_app_data": app_data
                  }

    tuner.execute(exec_config)

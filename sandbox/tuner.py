import ctypes
import _ctypes
import pipe
import constructs
import subprocess
import time
import string
import random
import os

# TODO:
# 1. Design in a robust way to handle encompassing of various search    ( )
#    spaces.
# 2. Introduce parallelism in code generation and compilation           ( )
# 3. Make the search configurations in each space a Set before 
#    enumerating                                                        ( )
#
def generate(_tuner_arg_dict):

    # unpack the arguments from the arg dictionary
    try:
        _tuner_app_name = _tuner_arg_dict['_tuner_app_name']
    except KeyError:
        print 'tuner : generator : \'_tuner_app_name\' - \
                not an optional parameter'

    try:
        _tuner_live_outs = _tuner_arg_dict['_tuner_live_outs']
    except KeyError:
        print 'tuner : generator : \'_tuner_live_outs\' - \
                not an optional parameter'

    try:
        _tuner_param_constraints = _tuner_arg_dict['_tuner_param_constraints']
    except KeyError:
        _tuner_param_constraints = None

    try:
        _tuner_param_estimates = _tuner_arg_dict['_tuner_param_estimates']
    except KeyError:
        _tuner_param_estimates = None

    try:
        _tuner_tile_size_configs = _tuner_arg_dict['_tuner_tile_size_configs']
    except KeyError:
        _tuner_tile_size_configs = None

    try:
        _tuner_group_size_configs = _tuner_arg_dict['_tuner_group_size_configs']
    except KeyError:
        _tuner_group_size_configs = None

    try:
        _tuner_opts = _tuner_arg_dict['_tuner_opts']
    except KeyError:
        _tuner_opts = []

    # for now, assume the following option to be deprecated
    '''
    try:
        _tuner_inline_directives = _tuner_arg_dict['_tuner_inline_directives']
    except KeyError:
        _tuner_inline_directives = None
    '''

    try:
        _tuner_dst_path = _tuner_arg_dict['_tuner_dst_path']
    except KeyError:
        _tuner_dst_path = '/tmp'

    try:
        _tuner_should_debug = _tuner_arg_dict['_tuner_should_debug']
    except KeyError:
        _tuner_should_debug = False
   


    # TODO:
    # 1. write a neat printer                   ( )
    #
    def print_params(to_file):
        if to_file:
            print_line(to_file)
            print >>to_file, "App name                : \""+_tuner_app_name+"\""
            print >>to_file, "Code Path               : \""+_tuner_dst_path+"\""
            print >>to_file, "Live-out Functions      :",
            for func in _tuner_live_outs:
                print >>to_file, "\""+func.name+"\",",
            print >>to_file, "\nParameter Constraints   :"
            for constraint in _tuner_param_constraints:
                print >>to_file, constraint
            print >>to_file, "\nParameter Estimates     :"
            for estimate in _tuner_param_estimates:
                print >>to_file, (estimate[0].name, estimate[1])
            print >>to_file, "\nTile Sizes              :"
            for tile_sizes in _tuner_tile_size_configs:
                print >>to_file, tile_sizes
            print >>to_file, "\nGroup Sizes             :"
            for group_sizes in _tuner_group_size_configs:
                print >>to_file, group_sizes
        if _tuner_should_debug:
            print_line()
            print "App name                : \""+_tuner_app_name+"\""
            print "Code Path               : \""+_tuner_dst_path+"\""
            print "Live-out Functions      :",
            for func in _tuner_live_outs:
                print "\""+func.name+"\",",
            print
            print "\nParameter Constraints :"
            for constraint in _tuner_param_constraints:
                print constraint
            print "\nParameter Estimates   :"
            for estimate in _tuner_param_estimates:
                print (estimate[0].name, estimate[1])
            print "\nTile Sizes            :"
            for tile_sizes in _tuner_tile_size_configs:
                print tile_sizes
            print "\nGroup Sizes           :"
            for group_sizes in _tuner_group_size_configs:
                print group_sizes
            '''
            print "\nfunctions to be inlined:"
            for func in _tuner_inline_directives:
                print "\""+func.name+"\",",
            '''

    app_name = _tuner_app_name+'_polymage_'

    def random_string():
        # for now, using a random string of length 10 
        # as subdir name
        return ''.join(random.SystemRandom().choice(
                    string.lowercase + \
                    string.uppercase + \
                    string.digits) \
                    for _ in xrange(10))
 
    # ensure that a subdirectory is created, with a name not conflicting with
    # existing ones
    dst_sub_dir = str(_tuner_dst_path)+'/'+'Poly'+random_string()+'Mage/'
    while os.path.exists(dst_sub_dir):
        dst_sub_dir =str(_tuner_dst_path)+'/'+'Poly'+random_string()+'Mage/'

    # subdirectories and files
    os.makedirs(dst_sub_dir)
    prog_prefix = str(dst_sub_dir)+str(app_name)
    config_file_name = str(dst_sub_dir)+'configurations.txt'
    config_file = open(config_file_name, 'a')

    # Compile String parts
    cxx='icpc'
    #cxx='g++-4.8'
    opt_flags='-openmp -xhost -O3 -ansi-alias'
    #opt_flags='-fopenmp -march=native -O3'
    shared_lib_flags='-fPIC -shared'
    include_flags='-I ../../memory_allocation/'
    other_CXX='../../memory_allocation/simple_pool_allocator.cpp'

    # Generate group sizes automatically, if none is specified by the user.
    # TODO:
    # 1. Limit the configuration space to the total number of
    #    functions in the pipeline, since any group size greater
    #    than that will lead to redundancy.                             ( )
    #
    if _tuner_group_size_configs == None:
        for i in range(1, 4):
            _tuner_group_size_configs.append(i)
        # Assuming that the number of functions does not exceed
        # this maximum, include a "group-all"
        _tuner_group_size_configs.append(200)

    print_params(config_file)
    config_file.close()

    # TODO:
    # 1. modify this whole looping, when 'thresholding' is implemented
    #    in the DSL compiler                                            ( )
    #
    _tuner_config = 0

    total_t1 = clock()
    # iterate over tile_sizes
    for _tuner_tile_size in _tuner_tile_size_configs:
        # iterate over group_sizes
        for _tuner_group_size in _tuner_group_size_configs:
            _tuner_config += 1

            # Update configs file:
            config_file = open(config_file_name, 'a')
            print_line(config_file)
            print >>config_file, "Config     : #"+str(_tuner_config)
            print >>config_file, "Tile Sizes : "+str(_tuner_tile_size)
            print >>config_file, "Group Size : "+str(_tuner_group_size)

            # .cpp and .so files
            c_file_name = str(prog_prefix)+str(_tuner_config)+'.cpp'
            so_file_name = str(prog_prefix)+str(_tuner_config)+'.so'
            #dot_file_name = str(prog_prefix)+str(_tuner_config)+'.dot'

            t1 = clock()

            # building the pipeline :
            _tuner_build_error = False
            try:
                _tuner_pipe = impipe.buildPipeline(
                            _tuner_live_outs,
                            param_constraints=_tuner_param_constraints,
                            param_estimates=_tuner_param_estimates,
                            tile_sizes=_tuner_tile_size,
                            group_size=_tuner_group_size,
                            options=_tuner_opts)
            except:
                _tuner_build_error = True
            finally:
                pass

            t2 = clock()

            # code generation :
            if _tuner_build_error is True:
                print >>config_file, "[ERROR] Build fail ..."
                print "[ERROR] Build fail ..."
            else:
                c_file = open(c_file_name, 'a')
                c_file.write(_tuner_pipe.generate_c_naive(is_extern_alloc=True, \
                                                    is_extern_func=True, \
                                                    are_params_void_ptrs=True).__str__())
                c_file.close()

                codegen_time = float(t2) - float(t1)
                print >>config_file, "Code Generation Time   : "+str(codegen_time*1000)

                '''
                g = _tuner_pipe.draw_pipeline_graph_with_groups()
                g.write(dot_file_name)
                '''

            # TODO:
            # 1. Currently inline_directives changes the whole DAG
            #    of the pipeline without cloning. Add this parameter
            #    after the compiler facilitates modification.           ( )
            #
  
            # TODO:
            # 1. Handle other CXX flags                                 ( )
            # 2. Include g++ support                                    ( )
            # 3. Should this be using just a 'make'?                    ( )
            #
            compile_string = cxx+' '+ \
                            opt_flags+' '+ \
                            shared_lib_flags+' '+ \
                            include_flags+' '+ \
                            other_CXX+' '+ \
                            c_file_name+' -o'+ \
                            so_file_name\

            t1 = clock()

            # compilation :
            _tuner_compile_error = False
            try:
                subprocess.check_output(compile_string, shell=True)
                pass
            except:
                _tuner_compile_error = True
            finally:
                pass

            t2 = clock()

            if _tuner_compile_error is True:
                print >>config_file, "[ERROR] Compilation aborted ..."
                print "[ERROR] Compilation aborted ..."
            else:
                compile_time = float(t2) - float(t1)
                print >>config_file, "Compilation Time       : "+str(compile_time*1000)
                # total time for this variant:
                print >>config_file, "Total                  : "+str((codegen_time+compile_time)*1000)

            config_file.close()

    total_t2 = clock()
    total_time = float(total_t2) - float(total_t1)

    config_file = open(config_file_name, 'a')
    print_line(config_file)
    print >>config_file, "Time taken to generate all variants : "
    print >>config_file, str(total_time*1000)

    print_line(config_file)
    config_file.close()

    return dst_sub_dir, _tuner_config, _tuner_pipe

def execute(_tuner_arg_dict):

    # unpack the arguments from the arg dictionary
    try:
        _tuner_pipe_arg_dict = _tuner_arg_dict['_tuner_pipe_arg_dict']
    except KeyError:
        print 'tuner : executer : \'_tuner_pipe_arg_dict\' - \
                not an optional parameter'

    try:
        _tuner_app_name = _tuner_arg_dict['_tuner_app_name']
    except KeyError:
        print 'tuner : executer : \'_tuner_app_name\' - \
                not an optional parameter'

    try:
        _tuner_pipe = _tuner_arg_dict['_tuner_pipe']
    except KeyError:
        print 'tuner : executer : \'_tuner_pipe\' - \
                not an optional parameter'

    try:
        _tuner_src_path = _tuner_arg_dict['_tuner_src_path']
    except KeyError:
        _tuner_src_path = "/tmp"

    try:
        _tuner_configs_count = _tuner_arg_dict['_tuner_configs_count']
    except KeyError:
        _tuner_configs_count = 0

    try:
        _tuner_omp_threads = _tuner_arg_dict['_tuner_omp_threads']
    except KeyError:
        _tuner_omp_threads = 1

    try:
        _tuner_nruns = _tuner_arg_dict['_tuner_nruns']
    except KeyError:
        _tuner_nruns = 1

    try:
        _tuner_should_debug = _tuner_arg_dict['_tuner_should_debug']
    except KeyError:
        _tuner_should_debug = True

    try:
        _tuner_custom_executor = _tuner_arg_dict['_tuner_custom_executor']
    except KeyError:
        _tuner_custom_executor = None

    def print_params(to_file=None):
        if to_file:
            print_line(to_file)
            print >>to_file, "App Name              : \""+_tuner_app_name+"\""
            print >>to_file, "Total Configurations  :", _tuner_configs_count
            print >>to_file, "OMP_NUM_THREADS       :",_tuner_omp_threads
            print >>to_file, "Number of Tuning Runs :",_tuner_nruns
        if _tuner_should_debug:
            print_line()
            print "app name              : \""+_tuner_app_name+"\""
            print "Total Configurations  :", _tuner_configs_count
            print "OMP_NUM_THREADS       :",_tuner_omp_threads
            print "Number of tuning runs :",_tuner_nruns

    # Parameters
    _tuner_pipe_params = _tuner_pipe.get_parameters()
    _tuner_pipe_params.sort(key=lambda x: x.name)

    # Inputs (Images)
    _tuner_pipe_inputs = _tuner_pipe.inputs
    _tuner_pipe_inputs.sort(key=lambda x: x.name)

    # Outputs
    _tuner_pipe_outputs = _tuner_pipe.outputs
    _tuner_pipe_outputs.sort(key=lambda x: x.name)

    # construct shared library function properties

    def convert_to_c_type(inp_type, inp_value):
        if inp_type == 'void':
            return ctypes.c_void(inp_value)

        if inp_type == 'char':
            return ctypes.c_char(inp_value)
        if inp_type == 'unsigned char':
            return ctypes.c_ubyte(inp_value)

        if inp_type == 'short' or \
            inp_type == 'short int':
            return ctypes.c_short(inp_value)
        if inp_type == 'unsigned short' or \
            inp_type == 'unsigned short int':
            return ctypes.c_ushort(inp_value)

        if inp_type == 'int':
            return ctypes.c_int(inp_value)
        if inp_type == 'unsigned' or \
            inp_type == 'unsigned int':
            return ctypes.c_uint(inp_value)

        if inp_type == 'long' or \
            inp_type == 'long int':
            return ctypes.c_long(inp_value)
        if inp_type == 'unsigned long' or \
            inp_type == 'unsigned long int':
            return ctypes.c_ulong(inp_value)

        if inp_type == 'long long' or \
            inp_type == 'long long int':
            return ctypes.c_longlong(inp_value)
        if inp_type == 'unsigned long long' or \
            inp_type == 'unsigned long long int':
            return ctypes.c_ulonglong(inp_value)

        if inp_type == 'float':
            return ctypes.c_float(inp_value)
        if inp_type == 'double':
            return ctypes.c_double(inp_value)
        if inp_type == 'long double':
            return ctypes.c_double(inp_value)

    lib_function_name = 'pipeline_'+_tuner_pipe.name

    pipe_func_args = []

    for param in _tuner_pipe_params:
        pipe_func_args += [convert_to_c_type(param.typ().c_name(), _tuner_pipe_arg_dict[param.name])]

    for inp in _tuner_pipe_inputs:
        pipe_func_args += [ctypes.c_void_p(_tuner_pipe_arg_dict[inp.name].ctypes.data)]

    for out in _tuner_pipe_outputs:
        pipe_func_args += [ctypes.c_void_p(_tuner_pipe_arg_dict[out.name].ctypes.data)]


    # set other variables
    app_name = _tuner_app_name+'_polymage_'
    prog_prefix = str(_tuner_src_path)+'/'+str(app_name)

    date_time_now = time.strftime("%d-%m-%Y_%H.%M.%S")
    tuning_report_file_name = str(_tuner_src_path)+'tuning_report'+'_'+str(date_time_now)+'.txt'

    tuning_report_file = open(tuning_report_file_name, 'a')
    print_params(tuning_report_file)
    tuning_report_file.close()

    _tuner_max_time = 1000000


    # TODO: iterate over thread count for tuning                ( )

    # set the thread-count
    os.environ["OMP_NUM_THREADS"] = str(_tuner_omp_threads)

    tuning_time_t1 = clock()

    global_min_config = 0
    global_min_time = _tuner_max_time
    for _tuner_config in range(1, _tuner_configs_count+1):
        print_line()
        print "config #"+str(_tuner_config),": ",

        # log to file
        tuning_report_file = open(tuning_report_file_name, 'a')
        print_line(tuning_report_file)
        print >>tuning_report_file, "Config #"+str(_tuner_config),": ",

        # load shared library, name the function
        so_file_name = str(prog_prefix)+str(_tuner_config)+'.so'

        # load shared library
        _tuner_load_error = False
        try:
            lib_pipeline = ctypes.cdll.Load_library(so_file_name)
        except:
            _tuner_load_error = True

        if _tuner_load_error :
            print >>tuning_report_file, "[ERROR] in loading shared library ..."
            print "[ERROR] in loading shared library ..."
        else:
            pipeline_func = lib_pipeline[lib_function_name]

            # TODO:
            # 1. Switch between report types - whether
            #    to print min / average / gmean etc.                ( )
            # 2. Graphics plots                                     ( )
            #

            # TODO:
            # 1. Use compiler info to get the correct prototype
            #    of the pipeline function to be called.             (X)
            # 2. Handle errors and dump them to log file            ( )

            # start timer
            local_min_time = _tuner_max_time

            # run
            for run in range(1, _tuner_nruns+1):
                t1 = clock()
                _tuner_runtime_error = False
                try:
                    if _tuner_custom_executor == None:
                        pipeline_func(*pipe_func_args)
                    else:
                        _tuner_custom_executor(lib_pipeline, \
                                           pipeline_func, \
                                           pipe_func_args, \
                                           _tuner_arg_dict)
                except:
                    _tuner_runtime_error = True
                t2 = clock()
                t_total = float(t2) - float(t1)
                #print "time taken: ", t2*1000 - t1*1000, "ms"

                # local minima
                if float(t_total) < float(local_min_time):
                    local_min_time = t_total

            if _tuner_runtime_error:
                print >>tuning_report_file, "[ERROR] Execution Aborted ..."
                print "[ERROR] Execution Aborted ..."
            else:
                # global minima
                if float(local_min_time) < float(global_min_time):
                    global_min_time = local_min_time
                    global_min_config = _tuner_config

                print local_min_time*1000, \
                      "(", global_min_time*1000, ")"
                # write the same to log file
                # FIXME: fuse this with stdout dump, to work like 'tee'         ( )
                print >>tuning_report_file, local_min_time*1000

                # FIXME:
                # Unloading the shared lib is almost impossible. Nevertheless, try to
                # dispose off the loaded lib object
                del pipeline_func
                _ctypes.dlclose(lib_pipeline._handle)
                del lib_pipeline

        tuning_report_file.close()

    tuning_report_file = open(tuning_report_file_name, 'a')

    tuning_time_t2 = clock()
    tuning_time = float(tuning_time_t2) - float(tuning_time_t1)

    print_line()
    print "Best Config :"
    print "Config #"+str(global_min_config)," -- ",
    print global_min_time*1000
    print "Src Path    : \""+_tuner_src_path+"\""

    print_line()
    print "Tuning Time :"
    print tuning_time*1000

    print_line()

    # write the same to log file
    # FIXME: fuse this with stdout dump, to work like 'tee'             ( )
    print_line(tuning_report_file)
    print >>tuning_report_file, "Best Config :"
    print >>tuning_report_file, "Config #"+str(global_min_config)," -- ",
    print >>tuning_report_file, global_min_time*1000
    print >>tuning_report_file, "Src Path    : \""+_tuner_src_path+"\""

    print_line(tuning_report_file)
    print >>tuning_report_file, "Tuning Time :"
    print >>tuning_report_file, tuning_time*1000

    print_line(tuning_report_file)
    tuning_report_file.close()

    return

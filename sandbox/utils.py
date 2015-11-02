import ctypes
import _ctypes

def convert_to_ctype(inp_type, inp_value):
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

def get_ordered_cfunc_params(pipe_object):
    # Parameters
    params = pipe_object.get_parameters()
    params.sort(key=lambda x: x.name)

    # Inputs (Images)
    inputs = pipe_object.inputs
    inputs.sort(key=lambda x: x.name)

    # Outputs
    outputs = pipe_object.outputs
    outputs.sort(key=lambda x: x.name)

    return params, inputs, outputs

def map_cfunc_args(params, inputs, outputs, arg_dict):
    func_args = []

    for param in params:
        func_args += [convert_to_ctype(param.typ().c_name(), arg_dict[param.name])]

    for inp in inputs:
        func_args += [ctypes.c_void_p(arg_dict[inp.name].ctypes.data)]

    for out in outputs:
        func_args += [ctypes.c_void_p(arg_dict[out.name].ctypes.data)]

    return func_args


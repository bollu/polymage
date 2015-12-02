import ctypes
import _ctypes

from fractions import gcd

def lcm(a, b):
    return a*b/(gcd(a, b))

def convert_to_ctype(inp_type, inp_value):
    if inp_type == 'void':
        return ctypes.c_void(inp_value)

    if inp_type == 'int8':
        return ctypes.c_char(inp_value)
    if inp_type == 'uint8':
        return ctypes.c_ubyte(inp_value)

    if inp_type == 'int16':
        return ctypes.c_short(inp_value)
    if inp_type == 'uint16':
        return ctypes.c_ushort(inp_value)

    if inp_type == 'int32':
        return ctypes.c_int(inp_value)
    if inp_type == 'uint32':
        return ctypes.c_uint(inp_value)

    if inp_type == 'int64':
        return ctypes.c_longlong(inp_value)
    if inp_type == 'uint64':
        return ctypes.c_ulonglong(inp_value)

    if inp_type == 'float':
        return ctypes.c_float(inp_value)
    if inp_type == 'double':
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

def map_cfunc_args(func_params, arg_data):
    func_args = []
    params = func_params[0]
    inputs = func_params[1]
    outputs = func_params[2]

    for param in params:
        func_args += [convert_to_ctype(param.typ().c_type_name(), arg_data[param.name])]

    for inp in inputs:
        func_args += [ctypes.c_void_p(arg_data[inp.name].ctypes.data)]

    for out in outputs:
        func_args += [ctypes.c_void_p(arg_data[out.name].ctypes.data)]

    return func_args


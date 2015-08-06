from __future__ import absolute_import, division, print_function

from collections import OrderedDict

import targetc as genc

def generate_code_for_group(g, body, cfunc_map, cparam_map):
    pass

def generate_code_for_pipeline(pipe, outsExternAlloc=True, is_io_void_ptr=True):
    sorted_groups = pipe.getOrderedGroups()
    # Discard the level order information
    sorted_groups = [ g[0] for g in sorted_groups ]
    # Create a top level module for the pipeline
    m = genc.CModule('Pipeline')

    # 1. Add header files which are requried by the pipeline
    with m.includes as inc_block:
        inc_block.add(genc.CInclude('stdio.h'))
        inc_block.add(genc.CInclude('stdlib.h'))
        inc_block.add(genc.CInclude('malloc.h'))
        inc_block.add(genc.CInclude('cmath'))
        inc_block.add(genc.CInclude('string.h'))
        inc_block.add(genc.CMacroDecl(genc.c_macro_min))
        inc_block.add(genc.CMacroDecl(genc.c_macro_max))
        inc_block.add(genc.CMacroDecl(genc.c_macro_floord))
        # TODO add the pool allocation header and option

    # 2. Add function blocks
    with m.funcs as func_block:

        # Maps from pipeline parameters and functions to c variables and
        # arrays. These maps are sent to each group for code generation.
        # They are updated during the code generation of each group to
        # include liveout functions of the group.
        cparam_map = {}
        cfunc_map = {}

        # Dictonary with all the pipeline arguments
        pipeline_args = OrderedDict()

        # 2.1. Collect all the inputs and parameters of the pipeline and
        # add them as pipeline function arguments.
        params = []
        for g in sorted_groups:
            params = params + g.getParameters()
        # Remove duplicates and sort by name
        params = list(set(params))
        params.sort(key=lambda x: x.name)

        # 2.1. collect pipeline parameters
        for param in params:
            cvar_type = genc.TypeMap.convert(param.typ)
            cvar = genc.CVariable(cvar_type, param.name)
            # Bind parameters to C variables
            cparam_map[param] = cvar
            pipeline_args[cvar] = cvar.typ

        # 2.2. collect inputs
        inputs = sorted(pipe.inputs, key=lambda x: x.name)
        for img in inputs:
            if is_io_void_ptr:
                img_type = genc.c_void
                cptr = genc.CPointer(img_type, 1)
                cvar = genc.CVariable(cptr, img.name+'_void_arg')
            else:
                img_type = genc.TypeMap.convert(img.typ)
                cptr = genc.CPointer(img_type, 1)
                cvar = genc.CVariable(cptr, img.name)
            pipeline_args[cvar] = cvar.typ
            # Bind input functions to C arrays
            carr = genc.CArray(img_type, img.name, img.dimensions,
                              'contigous')
            cfunc_map[img] = carr

        # 2.3. collect outputs
        outputs = sorted(pipe.outputs, key=lambda x: x.name)
        if outsExternAlloc:
            pass_by_type = genc.CReference
        else:
            pass_by_type = genc.CPointer

        for out in outputs:
            if is_io_void_ptr:
                out_typ = genc.c_void
                cptr = pass_by_type(out_typ, 1)
                cvar = genc.CVariable(cptr, out.name+'_void_arg')
            else:
                out_typ = genc.TypeMap.convert(out.typ)
                cptr = pass_by_type(out_typ, 1)
                cvar = genc.CVariable(cptr, out.name)

            pipeline_args[cvar] = cvar.typ

        # 2.4. function name and declaration
        cpipe_name = 'pipeline_' + pipe.name
        cpipe = genc.CFunction(genc.c_void, cpipe_name, pipeline_args)
        cpipe_decl = genc.CFunctionDecl(cpipe)
        cpipe_body = genc.CFunctionBody(cpipe_decl)

        func_block.add(cpipe_body)

        # 2.5. function body
        with cpipe_body.body as pbody:

            # If the code being generated is going to be compiled as shared
            # library, and used by python (through ctypes), the i/o data array
            # pointers will be given as void pointers. These should be casted
            # to their respective types first.

            # 2.5.1. typecast the i/o array ptrs
            if is_io_void_ptr:
                inouts = list(set(inputs) | set(outputs))
                for inout in inouts:
                    # actual input to be used
                    var_type = genc.TypeMap.convert(inout.typ)
                    var_ptr = genc.CPointer(var_type, 1)
                    var = genc.CVariable(var_type, inout.name)
                    var_decl = genc.CDeclaration(var_ptr, var)
                    pbody.add(var_decl)

                    # dummy void * input/output taken as argument
                    dummy_type = genc.TypeMap.convert(inout.typ)
                    dummy_ptr = genc.CPointer(dummy_type, 1)
                    dummy_cast = genc.CCast(dummy_ptr, inout.name+'_void_arg')
                    var_assign = genc.CAssign(var, dummy_cast)
                    pbody.add(var_assign)

            for g in sorted_groups:
                generate_code_for_group(g, pbody, cparam_map, cfunc_map)

    return m

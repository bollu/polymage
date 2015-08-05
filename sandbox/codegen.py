from __future__ import absolute_import, division, print_function

from collections import OrderedDict

import targetc as genc

def generate_code_for_group(g, body, cfunc_map, cparam_map):
    pass

def generate_code_for_pipeline(p, outsExternAlloc=False):
    sorted_groups = p.getOrderedGroups()
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
        inputs = sorted(p.inputs, key=lambda x: x.name)
        for img in inputs:
            img_type = genc.TypeMap.convert(img.typ)
            cptr = genc.CPointer(img_type, 1)
            cvar = genc.CVariable(cptr, img.name)
            pipeline_args[cvar] = cvar.typ
            # Bind input functions to C arrays
            carr = genc.CArray(img_type, img.name, img.dimensions,
                              'contigous')
            cfunc_map[img] = carr

        # 2.3. collect outputs
        outputs = sorted(p.outputs, key=lambda x: x.name)
        if outsExternAlloc:
            for out in outputs:
                out_typ = genc.TypeMap.convert(out.typ)
                cref = genc.CReference(outTyp, 1)
                cvar = genc.CVariable(cref, out.name)
                pipeline_args[cvar] = cvar.typ
        else:
            for out in outputs:
                outTyp = genc.TypeMap.convert(out.typ)
                cptr = genc.CPointer(outTyp, 1)
                cvar = genc.CVariable(cptr, out.name)
                pipeline_args[cvar] = cvar.typ

        # 2.4. function name and declaration
        cpipe_name = 'pipeline_' + p.name
        cpipe = genc.CFunction(genc.c_void, cpipe_name, pipeline_args)
        cpipe_decl = genc.CFunctionDecl(cpipe)
        cpipe_body = genc.CFunctionBody(cpipe_decl)

        func_block.add(cpipe_body)
        with cpipe_body.body as pbody:
            for g in sorted_groups:
                generate_code_for_group(g, pbody, cparam_map, cfunc_map)
    return m

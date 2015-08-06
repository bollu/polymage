from __future__ import absolute_import, division, print_function

from collections import OrderedDict

import pipe
from constructs import *
import expression as expr
import targetc as genc
import logging

# LOG CONFIG #
codegen_logger = logging.getLogger("codegen.py")
codegen_logger.setLevel(logging.DEBUG)
LOG = codegen_logger.log

def generate_code_for_group(pipeline, g, body, options, \
                            cfunc_map, cparam_map, \
                            outputs, outsExternAlloc):

    g.polyRep.generateCode()

    # NOTE uses the levelNo of the first polypart of each compute object of
    # the group as the key for sorting compare operator. *Idea is that all
    # parts of a compute object bears the same levelNo*, thus repeated calling
    # of 'orderComputeObjs' can be avoided.
    group_parts = g.polyRep.polyParts
    sorted_comp_objs = sorted(g._compObjs, \
                              key = lambda \
                              comp : group_parts[comp][0].levelNo)

    # NOTE the last comp obj in the sorted list has the max level number
    # that is shared by all the liveouts of the group
    last_comp = sorted_comp_objs[len(sorted_comp_objs)-1]
    max_level = group_parts[last_comp][0].levelNo
    is_comp_liveout = {}
    is_comp_output = {}
    for comp in sorted_comp_objs:
        is_comp_liveout[comp] = (group_parts[comp][0].levelNo == max_level)
        is_comp_output[comp] = comp in outputs

    # ***
    log_str = str([comp.name for comp in sorted_comp_objs])
    LOG(logging.DEBUG, log_str)
    # ***

    # TODO free the arrays ASAP (compaction)
    # list of arrays to be freed
    freelist = []

    # Boolean to check if the memory allocations should be done using malloc
    # or the custom pool allocator
    pooled = 'pool_alloc' in pipeline._options

    for comp in sorted_comp_objs:
        # ***
        LOG(logging.DEBUG, comp.name)
        # ***

        # nothing to do!
        if isinstance(comp, pipe.Image):
            continue

        # 1. Allocate storage
        array_dims = len(comp.variables)
        is_output = is_comp_output[comp]
        is_liveout = is_comp_liveout[comp]

        # 1.1. scratchpad allocation, wherever applicable
        reduced_dims = [ -1 for i in xrange(0, len(comp.domain))]
        scratch = [ False for i in xrange(0, len(comp.domain))]
        if not is_liveout and not is_output:
            for part in group_parts[comp]:
                for i in xrange(0, len(comp.domain)):
                    if i in part.dim_scratch_size:  # as a key
                        reduced_dims[i] = max(reduced_dims[i], \
                                              part.dim_scratch_size[i])
                        scratch[i] = True

        # 1.2. prepare the sizes of each dimension for allocation
        dims = []
        for i in xrange(0, len(comp.domain)):
            interval = comp.domain[i]
            if reduced_dims[i] == -1:
                # NOTE interval step is always +1
                dim_expr = expr.simplifyExpr(interval.upperBound -
                                             interval.lowerBound + 1)

                # FIXME Creating both a cVariable (for C declaration) and a
                # Variable (to append to dims) with same namestring.
                # Link both these properly.
                dim_var_name = genc.CNameGen.get_temp_var_name()

                dim_var = Variable(Int, dim_var_name)
                dim_c_var = genc.CVariable(genc.c_int, dim_var_name)

                dim_var_decl = \
                            genc.CDeclaration(genc.c_int, dim_c_var, dim_expr)

                body.add(dim_var_decl)
                dims.append(dim_var)
            else:
                dims.append(reduced_dims[i])

        # 1.3. declare and allocate arrays
        array_type = genc.TypeMap.convert(comp.typ)
        array = genc.CArray(array_type, comp.name, dims)
        array_ptr = genc.CPointer(array_type, 1)
        cfunc_map[comp] = (array, scratch)

        if is_liveout:
            array.layout = 'contigous'
            # do not allocate for output arrays if they are already allocated
            if not is_output or not outsExternAlloc:
                array_decl = genc.CDeclaration(array_ptr, array)
                body.add(array_decl)
                array.allocate_contigous(body, pooled)

        # array is freed, if comp is a group liveout and not an output
        if not is_output and is_liveout:
            freelist.append(array)

        # 1.4. generate scan loops

    # 2. generate code for built isl ast
    if g.polyRep.polyast != []:
        for ast in g.polyRep.polyast:
            #generateCNaiveFromIslAst(ast, body, cfuncMap, cparamMap)
            pass

    # 3. Deallocate storage
    for array in freelist:
        array.deallocate(body, pooled)
        pass

    return

def generate_code_for_pipeline(pipeline, outsExternAlloc=True, is_io_void_ptr=True):
    sorted_groups = pipeline.getOrderedGroups()
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
        inputs = sorted(pipeline.inputs, key=lambda x: x.name)
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
                              'contiguous')
            cfunc_map[img] = carr

        # 2.3. collect outputs
        outputs = sorted(pipeline.outputs, key=lambda x: x.name)
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
        cpipe_name = 'pipeline_' + pipeline.name
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
                generate_code_for_group(pipeline, g, pbody, \
                                        pipeline._options, \
                                        cfunc_map, cparam_map, \
                                        outputs, outsExternAlloc)

    return m

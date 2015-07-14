from __future__ import absolute_import, division, print_function

from collections import OrderedDict

import targetc as genc

def generateCodeForGroup(g, body, cfuncMap, cparamMap):
    pass

def generateCodeForPipeline(p, outsExternAlloc = False):
    sortedGroups = p.getOrderedGroups()
    # Discard the level order information 
    sortedGroups = [ g[0] for g in sortedGroups ]
    # Create a top level module for the pipeline
    m = genc.cModule('Pipeline')
    # Add header files which are requried by the pipeline   
    with m.includes as incblock:
        incblock.add(genc.cInclude('stdio.h'))
        incblock.add(genc.cInclude('stdlib.h'))
        incblock.add(genc.cInclude('malloc.h'))
        incblock.add(genc.cInclude('cmath'))
        incblock.add(genc.cInclude('string.h'))
        incblock.add(genc.cMacroDecl(genc.cMacroMin))
        incblock.add(genc.cMacroDecl(genc.cMacroMax))
        incblock.add(genc.cMacroDecl(genc.cMacroFloord))
        # TODO add the pool allocation header and option
    with m.funcs as funcblock:
        # Collect all the inputs and parameters of the pipeline
        # and add them as pipeline function arguments.
        params = []
        for g in sortedGroups:
            params = params + g.getParameters()
       
        # Remove duplicates and sort by name
        params = list(set(params))
        params.sort(key=lambda x: x.name)
        
        # Maps from pipeline parameters and functions to c 
        # variables and arrays. These maps are sent to each 
        # group for code generation. They are updated during
        # the code generation of each group to include liveout
        # functions of the group.

        cparamMap = {}
        cfuncMap = {}

        # Dictonary with all the pipeline arguments
        pipelineArgs = OrderedDict()

        for param in params:
            cvarType = genc.TypeMap.convert(param.typ)
            cvar = genc.cVariable(cvarType, param.name)
            # Binding parameters to C variables
            cparamMap[param] = cvar
            pipelineArgs[cvar] = cvar.typ
       
        inputs = sorted(p.inputs, key=lambda x: x.name)
        for img in inputs:
            imgType = genc.TypeMap.convert(img.typ)
            cptr = genc.cPointer(imgType, 1)
            cvar = genc.cVariable(cptr, img.name)
            pipelineArgs[cvar] = cvar.typ
            # Binding input functions to C arrays
            carr = genc.cArray(imgType, img.name, img.dimensions,
                              'contigous')
            cfuncMap[img] = carr

        outputs = sorted(p.outputs, key=lambda x: x.name) 

        if outsExternAlloc:
            for out in outputs:
                outTyp = genc.TypeMap.convert(out.typ)
                cref = genc.cReference(outTyp, 1)
                cvar = genc.cVariable(cref, out.name)
                pipelineArgs[cvar] = cvar.typ
        else:
            for out in outputs:
                outTyp = genc.TypeMap.convert(out.typ)
                cptr = genc.cPointer(outTyp, 1)
                cvar = genc.cVariable(cptr, out.name)
                pipelineArgs[cvar] = cvar.typ
        
        cpipeName = 'pipeline_' + p.name
        cpipe = genc.cFunction(genc.cVoid, cpipeName, pipelineArgs)
        cpipeDecl = genc.cFunctionDecl(cpipe)
        cpipeBody = genc.cFunctionBody(cpipeDecl)
      
        funcblock.add(cpipeBody)

        with cpipeBody.body as pbody:
            for g in sortedGroups:
                generateCodeForGroup(g, pbody, cparamMap, cfuncMap)
    return m

from __future__ import absolute_import, division, print_function

import targetc as genc

def generateCodeForGroup(g):
    pass

def generateCodeForPipeline(p):
   sortedGroups = p.getOrderedGroups()
   # Discard the level order information 
   sortedGroups = [ g[0] for g in sortedGroups ]
   # Create a top level module for the pipeline
   m = genc.cModule('PipeLine')
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

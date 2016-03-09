import numpy as np
import time
import sys

from init import init_all
from verify import verify_norm
from builder import create_lib, build_resid, build_mg_3p
from execMG import multigrid

sys.path.insert(0, '../../../')
from misc import printLine

app = 'nas-pb-mg-3.2'
probClasses = ['S', 'W', 'A', 'B', 'C', 'D']

#-------------------------------------------------------------------
# initialize parameters

def header():
    print "NAS Parallel Benchmark v3.2"
    print "            MG"

def usage():
    print "[main]: Usage: "
    print "[main]: "+sys.argv[0]+" <class> <mode>"
    print "[main]: 'class' :: {'S', 'W', 'A', 'B', 'C', 'D'}"
    print "[main]: 'mode'  :: {'new', 'existing', 'tune'}"

#-------------------------------------------------------------------
# main

print
printLine()
header()

if len(sys.argv) > 2:
    probClass = sys.argv[1]
    mode = sys.argv[2]
    if probClass not in probClasses:
        print '[main]: Invalid problem Class'
        usage()
        sys.exit(1)
else:
    usage()
    sys.exit(1)

#-------------------------------------------------------------------
dataDict = {}
impipeDict = {}

dataDict['probClass'] = probClass

# init all the required data
initAll(impipeDict, dataDict)
#-------------------------------------------------------------------

if mode != 'tune':
    # setting up multigrid v-cycle computation
    createPipeLib(buildMg3P, "mgU_mgR_", impipeDict, dataDict, mode)

    # setting up standalone version of residual computation
    createPipeLib(buildResid, "resid", impipeDict, dataDict, mode)
#fi
#-------------------------------------------------------------------

nx = ny = nz = dataDict['probSize']

# Setup report
printLine()
print "# Problem Settings #"
print "[main]: CLASS        = \""+dataDict['probClass']+"\""
print "[main]: top level    =", dataDict['lt']
print "[main]: bottom level =", dataDict['lb']
print "[main]: grid size    =", nx, "x", ny, "x", nz
print "[main]: n-iterations =", dataDict['nit']

print
print "# Stencil Co-efficients #"
print "[main]: a =", dataDict['a']
print "[main]: c =", dataDict['c']

verifyDict = dataDict['verifyDict']
print
print "# Verification Values #"
print "[main]: threshold         =", verifyDict['epsilon']
print "[main]: Class \""+dataDict['probClass']+"\" " \
            + "L2 norm =", verifyDict['verifyValue']

print
print "# Initial Norms #"
print "[main]: initial norm =", dataDict['rnm2']
print "[main]: initial err  =", dataDict['rnmu']
#-------------------------------------------------------------------
printLine()
print
print "MULTIGRID EXECUTION STARTS"
print
printLine()
#-------------------------------------------------------------------
multigrid(impipeDict, dataDict)

#-------------------------------------------------------------------
printLine()
print "[main]: Verifying the results ..."
print
verifyNorm(dataDict)
#-------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import sys
import subprocess
from fractions  import Fraction

from polymage_residual import residPipe
from polymage_mg3P     import mg3pPipe
from compiler          import cCompile
from misc              import loadLib

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def codeGen(pipe, fileName, dataDict):
    print
    print("[builder]: writing the code to", fileName, "...")
     
    code = pipe.generateCNaive(isExternAlloc=True,
                               isExternFunc=True,
                               areParamsVoidPtrs=True)

    f = open(fileName, 'w')
    f.write(code.__str__())
    f.close()

    return

def graphGen(pipe, fileName, dataDict):
    graphFile = fileName+".dot"
    pngGraph = fileName+".png"

    print()
    print("[builder]: writing the graph dot file to", graphFile, "...")

    graph = pipe.originalGraph
    graph.write(graphFile)
    print("[builder]: ... DONE")

    dottyStr = "dot -Tpng "+graphFile+" -o "+pngGraph

    print
    print("[builder]: drawing the graph using dotty to", pngGraph)
    print(">", dottyStr)
    subprocess.check_output(dottyStr, shell=True)
    print("[builder]: ... DONE")

    return

def buildResid(impipeDict, dataDict):

    # construct the residual pipeline on the finest grid
    r = residPipe(impipeDict, dataDict)

    n = impipeDict['n']
    N = dataDict['N']
    lt = dataDict['lt']

    liveOuts = [r]

    pEstimates = [(n, N[lt])]
    pConstraints = [ Condition(n, "==", N[lt]) ]
    tSize = [16, 16, 16]
    gSize = 1
    opts = ["pool_alloc"]

    '''
    rPipe = buildPipeline(liveOuts,
                                 paramEstimates=pEstimates,
                                 paramConstraints=pConstraints,
                                 tileSizes=tSize,
                                 groupSize=gSize,
                                 inlineDirectives=[],
                                 options=opts)
    '''
    rPipe = buildPipeline(liveOuts)

    return rPipe

def buildMg3P(impipeDict, dataDict):
    
    # construct the multigrid v-cycle pipeline
    mg_u, mg_r = mg3pPipe(impipeDict, dataDict)

    n = impipeDict['n']
    N = dataDict['N']
    lt = dataDict['lt']

    gridDict = dataDict['gridDict']
    u = gridDict['u']
    v = gridDict['v']
    r = gridDict['r']

    liveOuts = [mg_u, mg_r]

    pEstimates = [(n, N[lt])]
    pConstraints = [ Condition(n, "==", N[lt]) ]
    tSize = [16, 128, 128]
    gSize = 10
    opts = ["pool_alloc"]

    '''
    mgPipe = buildPipeline(liveOuts,
                                  paramEstimates=pEstimates,
                                  paramConstraints=pConstraints,
                                  tileSizes=tSize,
                                  groupSize=gSize,
                                  inlineDirectives=[],
                                  options=opts)
    '''
    mgPipe = buildPipeline(liveOuts)

    return mgPipe

def createPipeLib(buildFunc, pipeName, impipeDict, dataDict, mode):
    pipeSrc  = pipeName+".cpp"
    pipeSo   = pipeName+".so"

    if mode == 'new':
        # build the polymage pipeline
        pipe = buildFunc(impipeDict, dataDict)

        # draw the pipeline graph to a png file
        graphGen(pipe, pipeName, dataDict)

        print("Stop here")
        assert(False)

        # generate pipeline cpp source
        codeGen(pipe, pipeSrc, dataDict)
    #fi

    if mode != 'ready':
        # compile the cpp code
        cCompile(pipeSrc, pipeSo, cCompiler="intel")
    #fi

    # load the shared library
    pipeFuncName = "pipeline_"+pipeName
    loadLib(pipeSo, pipeFuncName, dataDict)

    return

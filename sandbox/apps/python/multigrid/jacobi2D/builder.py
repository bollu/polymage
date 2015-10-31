import sys
import subprocess

from cpp_compiler    import cCompile
from loader          import loadLib
from polymage_vcycle import vCycle
from polymage_wcycle import wCycle

from compiler   import *
from constructs import *

def codeGen(pipe, fileName, dataDict):
    print("")
    print("[builder]: writing the code to", fileName, "...")

    code = pipe.generate_code(is_extern_c_func=True,
                              outputs_no_alloc=True,
                              are_io_void_ptrs=True)

    f = open(fileName, 'w')
    f.write(code.__str__())
    f.close()

    return

'''
def graphGen(pipe, fileName, dataDict):
    graphFile = fileName+".dot"
    pngGraph = fileName+".png"

    print("")
    print("[builder]: writing the graph dot file to", graphFile, "...")

    #graph = pipe.drawPipelineGraphWithGroups()
    graph = pipe.originalGraph
    graph.write(graphFile)
    print("[builder]: ... DONE")

    dottyStr = "dot -Tpng "+graphFile+" -o "+pngGraph

    print("")
    print("[builder]: drawing the graph using dotty to", pngGraph)
    print(">", dottyStr)
    subprocess.check_output(dottyStr, shell=True)
    print("[builder]: ... DONE")

    return
'''

def buildMGCycle(impipeDict, dataDict):
    cycleType = dataDict['cycle']

    if cycleType == 'V':
        # construct the multigrid v-cycle pipeline
        mg = vCycle(impipeDict, dataDict)
    elif cycleType == 'W':
        # construct the multigrid w-cycle pipeline
        mg = wCycle(impipeDict, dataDict)

    n = impipeDict['n']
    N = dataDict['N']
    L = dataDict['L']

    liveOuts = [mg]
    pipeName = dataDict['cycle_name']
    pEstimates = [(n, dataDict['n'])]
    pConstraints = [ Condition(n, "==", dataDict['n']) ]
    tSize = [16, 16, 16]
    gSize = 10
    opts = []
    if dataDict['pool_alloc'] == True:
        opts += ['pool_alloc']

    '''
    mgPipe = buildPipeLine(liveOuts,
                           paramEstimates=pEstimates,
                           paramConstraints=pConstraints,
                           tileSizes=tSize,
                           groupSize=gSize,
                           inlineDirectives=[],
                           options=opts)
    '''

    mgPipe = buildPipeline(liveOuts,
                           param_estimates=pEstimates,
                           param_constraints=pConstraints,
                           tile_sizes = tSize,
                           group_size = gSize,
                           options = opts,
                           pipe_name = pipeName)

    return mgPipe

def createLib(buildFunc, pipeName, impipeDict, dataDict, mode):
    pipeSrc  = pipeName+".cpp"
    pipeSo   = pipeName+".so"

    if buildFunc != None:
        if mode == 'new':
            # build the polymage pipeline
            pipe = buildFunc(impipeDict, dataDict)

            # draw the pipeline graph to a png file
            #graphGen(pipe, pipeName, dataDict)

            # generate pipeline cpp source
            codeGen(pipe, pipeSrc, dataDict)
        #fi
    #fi

    if mode != 'ready':
        # compile the cpp code
        cCompile(pipeSrc, pipeSo, cCompiler="gnu")
    #fi

    # load the shared library
    pipeFuncName = "pipeline_"+pipeName
    loadLib(pipeSo, pipeFuncName, dataDict)

    return

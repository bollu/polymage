import sys
import subprocess

from cpp_compiler    import cCompile
from loader          import loadLib
from polymage_vcycle import vCycle
from polymage_wcycle import wCycle

from compiler   import *
from constructs import *

def codeGen(pipe, fileName, appData):
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
def graphGen(pipe, fileName, appData):
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

def buildMGCycle(pipeData, appData):
    cycleType = appData['cycle']

    if cycleType == 'V':
        # construct the multigrid v-cycle pipeline
        mg = vCycle(pipeData, appData)
    elif cycleType == 'W':
        # construct the multigrid w-cycle pipeline
        mg = wCycle(pipeData, appData)

    n = pipeData['n']

    liveOuts = [mg]
    pipeName = appData['cycle_name']
    pEstimates = [(n, appData['n'])]
    pConstraints = [ Condition(n, "==", appData['n']) ]
    tSize = [16, 16, 16]
    gSize = 5
    opts = []
    if appData['pool_alloc'] == True:
        opts += ['pool_alloc']

    mgPipe = buildPipeline(liveOuts,
                           param_estimates=pEstimates,
                           param_constraints=pConstraints,
                           tile_sizes = tSize,
                           group_size = gSize,
                           options = opts,
                           pipe_name = pipeName)

    return mgPipe

def createLib(buildFunc, pipeName, pipeData, appData, mode):
    pipeSrc  = pipeName+".cpp"
    pipeSo   = pipeName+".so"

    if buildFunc != None:
        if mode == 'new':
            # build the polymage pipeline
            pipe = buildFunc(pipeData, appData)

            # draw the pipeline graph to a png file
            #graphGen(pipe, pipeName, appData)

            # generate pipeline cpp source
            codeGen(pipe, pipeSrc, appData)
        #fi
    #fi

    if mode != 'ready':
        # compile the cpp code
        cCompile(pipeSrc, pipeSo, cCompiler="gnu")
    #fi

    # load the shared library
    pipeFuncName = "pipeline_"+pipeName
    loadLib(pipeSo, pipeFuncName, appData)

    return

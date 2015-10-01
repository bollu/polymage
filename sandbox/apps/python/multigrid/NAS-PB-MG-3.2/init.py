import numpy as np

from verify          import setVerification
from misc            import loadLib, unpackInput, \
                            ilog2, makePeriodic
from execMG          import calcNorm
from polymage_common import setCases, setVars

def initSizes(dataDict):
    # init the problem size and other parameters
    # and the respective solutions based on the
    # input class (from NAS PB)
    # TODO: generalize for custom specs

    Class = dataDict['probClass']

    if Class == 'S':
        dataDict['probSize'] = 32
        dataDict['nit'] = 4
    elif Class == 'W':
        dataDict['probSize'] = 128
        dataDict['nit'] = 4
    elif Class == 'A':
        dataDict['probSize'] = 256
        dataDict['nit'] = 4
    elif Class == 'B':
        dataDict['probSize'] = 256
        dataDict['nit'] = 20
    elif Class == 'C':
        dataDict['probSize'] = 512
        dataDict['nit'] = 20
    elif Class == 'D':
        dataDict['probSize'] = 1024
        dataDict['nit'] = 50
    else:
        dataDict['probSize'] = 1
        dataDict['nit'] = 1

    lt = ilog2(dataDict['probSize'])
    lb = 1

    # top and bottom levels
    # register in the data dictionary
    dataDict['lt'] = lt
    dataDict['lb'] = lb

    # grid sizes (with ghost) at each level
    N = {}
    N[lt] = dataDict['probSize']+2
    for l in range(lt-1, 0, -1):
        N[l] = (N[l+1]-2)/2+2

    # register in the data dictionary
    dataDict['N'] = N

    return

def initCoefs(dataDict):

    a = np.zeros((4), np.float64)
    c = np.zeros((4), np.float64)

    Class = dataDict['probClass']

    a[0] = -8.0/3.0
    a[1] =  0.0
    a[2] =  1.0/6.0
    a[3] =  1.0/12.0

    if Class == 'A' or Class == 'S' or Class == 'W':
        c[0] = -3.0/8.0
        c[1] =  1.0/32.0
        c[2] = -1.0/64.0
        c[3] =  0.0
    else:
        c[0] = -3.0/17.0
        c[1] =  1.0/33.0
        c[2] = -1.0/61.0
        c[3] =  0.0

    # register in the data dictionary
    dataDict['a'] = a
    dataDict['c'] = c

    return

def initGrids(dataDict):
    withGhost = dataDict['probSize'] + 2
    Class = dataDict['probClass']

    # grids with ghost zone
    n1 = n2 = n3 = withGhost
    v = np.zeros((n1, n2, n3), np.float64)
    u = np.zeros((n1, n2, n3), np.float64)
    r = np.zeros((n1, n2, n3), np.float64)

    u1 = np.zeros((n1, n2, n3), np.float64)
    r1 = np.zeros((n1, n2, n3), np.float64)

    # read the input file and unpack the co-ordinates
    index = np.zeros((20, 3), np.int16)
    unpackInput(index, "inputs/input."+Class)

    # set the initial values at these co-ord.s
    for i in range(0, 10):
        v[index[i,0], index[i,1], index[i,2]] = -1.0
    for i in range(10, 20):
        v[index[i,0], index[i,1], index[i,2]] = 1.0

    #makePeriodic(v, dataDict)

    gridDict = {}

    # register in grid dictionary
    gridDict['v'] = v
    gridDict['u'] = u
    gridDict['r'] = r

    gridDict['u1'] = u1
    gridDict['r1'] = r1

    dataDict['gridDict'] = gridDict

    return

def initNorm(dataDict):
    rnm2 = 0.0
    rnmu = 0.0

    gridDict = dataDict['gridDict']
    v = gridDict['v']

    # load the norm computation library
    libFile = "./norm.so"
    libFuncName = "norm2u3"

    loadLib(libFile, libFuncName, dataDict)

    # calculate norm on the initial grid
    calcNorm(v, dataDict)

    return

def initPolyMageData(impipeDict, dataDict):
    # set variables, intervals etc
    setVars(impipeDict, dataDict)

    # set boundary case conditions
    setCases(impipeDict, dataDict)

    return

def initAll(impipeDict, dataDict):
    # initialize problem sizes and relevant params
    initSizes(dataDict)
    
    # initialize the stencil co-efficients
    initCoefs(dataDict)
    
    # set the verification values
    setVerification(dataDict, isPeriodic=False)
    
    # initialize the grid contents
    initGrids(dataDict)
    
    # initialize the norms
    initNorm(dataDict)
    
    # initialize polymage specific parameters
    initPolyMageData(impipeDict, dataDict)

    return

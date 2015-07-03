import ctypes
import numpy as np

def printLine(toFile=None):
    if toFile:
        print >>toFile, "--------------------------------------------------"
    else:
        print "--------------------------------------------------"

def loadLib(libFile, libFuncName, dataDict):
    print
    print "[misc]: loading the shared library", libFile, "..."

    # assuming that it is present,
    # load the shared library
    lib = ctypes.cdll.LoadLibrary('./'+libFile)

    print "[misc]: ... DONE"

    # name of the lib function
    libFunc = lib[libFuncName]

    # register the library and the function in the data dictionary
    dataDict[str(libFile)] = lib
    dataDict[str(libFuncName)] = libFunc

    return

def makePeriodic(v, dataDict):
    N = dataDict['N']
    lt = dataDict['lt']

    n = N[lt]

    # for each ghost plane, copy the
    # first non-ghost plane at the other end

    # communicate left <-> right planes
    v[1:n-2, 1:n-2, 0  ] = v[1:n-2, 1:n-2, n-2]
    v[1:n-2, 1:n-2, n-1] = v[1:n-2, 1:n-2, 1  ]

    # communicate top <-> bottom planes
    v[1:n-2, 0  , 0:n-1] = v[1:n-2, n-2, 0:n-1]
    v[1:n-2, n-1, 0:n-1] = v[1:n-2, 1  , 0:n-1]

    # communicate front <-> back planes
    v[0  , 0:n-1, 0:n-1] = v[n-2, 0:n-1, 0:n-1]
    v[n-1, 0:n-1, 0:n-1] = v[1  , 0:n-1, 0:n-1]

    return

# takes in an integer and returns its
# logarithm to the base 2
def ilog2(n):
    l = -1
    while n >= 1:
        l += 1
        n /= 2

    return l

def unpackInput(index, fileName):
    f = open(fileName, 'r')

    i = 0
    for line in f:
        lineInt = []
        lineInt.append([int(x) for x in line.split()])
        for j in range(0, 3):
            index[i, j] = lineInt[0][j]-1
        i += 1

    f.close()

    return


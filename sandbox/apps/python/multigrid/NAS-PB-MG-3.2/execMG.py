import sys
import os
import ctypes
import numpy as np
import time

from verify import verifyNorm

sys.path.insert(0, '../../../')

def calcNorm(v, dataDict):
    N = v.shape[0]

    # lib function name
    norm2u3 = dataDict['norm2u3']

    rnm2 = np.zeros((1), np.float64)
    rnmu = np.zeros((1), np.float64)

    # lib function args
    normArgs = []
    normArgs += [ctypes.c_int(N)]
    normArgs += [ctypes.c_void_p(rnm2.ctypes.data)]
    normArgs += [ctypes.c_void_p(rnmu.ctypes.data)]
    normArgs += [ctypes.c_void_p(v.ctypes.data)]

    # call lib function
    norm2u3(*normArgs)

    # register the norm values in the data dictionary
    dataDict['rnm2'] = rnm2[0]
    dataDict['rnmu'] = rnmu[0]

    return

def calcResid(u, v, r, dataDict):
    gridDict = dataDict['gridDict']

    N = r.shape[0]

    # lib function name
    resid = dataDict['pipeline_resid']

    # lib function args
    residArgs = []
    residArgs += [ctypes.c_int(N)]
    residArgs += [ctypes.c_void_p(u.ctypes.data)]
    residArgs += [ctypes.c_void_p(v.ctypes.data)]
    residArgs += [ctypes.c_void_p(r.ctypes.data)]

    # call the function
    resid(*residArgs)

    return

def callMg3P(u, v, r, uOut, rOut, dataDict, it):
    gridDict = dataDict['gridDict']

    N = r.shape[0]

    # lib function name
    mg = dataDict['pipeline_mgU_mgR_']

    # lib function args
    mgArgs = []
    mgArgs += [ctypes.c_int(N)]
    mgArgs += [ctypes.c_void_p(r.ctypes.data)]
    mgArgs += [ctypes.c_void_p(u.ctypes.data)]
    mgArgs += [ctypes.c_void_p(v.ctypes.data)]
    mgArgs += [ctypes.c_void_p(rOut.ctypes.data)]
    mgArgs += [ctypes.c_void_p(uOut.ctypes.data)]

    # call the function
    mg(*mgArgs)

    return

def multigrid(impipeDict, dataDict):
    gridDict = dataDict['gridDict']

    u = gridDict['u']
    v = gridDict['v']
    r = gridDict['r']
    r1 = gridDict['r1']
    u1 = gridDict['u1']

    libMg = dataDict['mgU_mgR_.so']

    # compute the initial residual
    print "[exec]: computing the initial residual ..."
    calcResid(u, v, r, dataDict)
    print "[exec]: ... DONE"

    # calculate the initial residual norm
    print "[exec]: calculating the initial norm ..."
    calcNorm(r, dataDict)
    print "[exec]: ... DONE"

    print
    print "[exec]: norm =", dataDict['rnm2']
    print "        err  =", dataDict['rnmu']

    libMg.pool_init()

    # mg3p time
    tTotal = 0.0

    nit = dataDict['nit']
    nit = 10
    # call 'nit' v-cycles
    for it in range(1, nit+1):
        print
        print "[exec]: iter", it

        # timer ON
        t1 = time.clock()
        
        # call the v-cycle
        if it % 2 == 1:
            callMg3P(u, v, r, u1, r1, dataDict, it)
        else:
            callMg3P(u1, v, r1, u, r, dataDict, it)
        #fi

        # timer OFF
        t2 = time.clock()
        # mg3P time +=
        tTotal += (float(t2) - float(t1))

        # update the residual
        if it % 2 == 1:
            calcResid(u1, v, r1, dataDict)
        else:
            calcResid(u, v, r, dataDict)
        #fi

    #endfor
    libMg.pool_destroy()

    print
    print "[exec]: calculating the final norm ..."
    # calculate the final residual norm
    if nit % 2 == 1:
        calcNorm(r1, dataDict)
    else:
        calcNorm(r, dataDict)
    #fi
    print "[exec]: ... DONE"

    print
    print "[exec]: time taken by multigrid pipeline =", tTotal, "s"

    return

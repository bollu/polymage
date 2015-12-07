from __future__ import absolute_import, division, print_function

import sys
sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

from fractions  import Fraction

def setGhosts(r, ghosts, value):
    for ghost in ghosts:
        r.defn.append(Case(ghosts[ghost], value))

    return

def setVars(impipeData, appData):
    L = appData['L']

    y = Variable(Int, "y")
    x = Variable(Int, "x")

    impipeData['y'] = y
    impipeData['x'] = x

    n = Parameter(Int, "n")
    impipeData['n'] = n

    # grid size at each level
    N = {}
    # jacobi weight (2D)
    omega = 4.0/5.0
    # mulitplier in the jacobi computation
    jacobi_c = {}
    # 1.0/(h*h)
    invhh = {}

    for l in range(0,L+1):
        if l == 0:
            N[0] = n
        else:
            N[l] = 2*N[l-1]+1

        h = 1.0/(N[l]+1)
        invhh[l] = 1.0/(h*h)

        # omega.D^-1 = omega.(d^-1.I) = omega * h^2/6.0
        # d^-1 depends on diagonal elements of A^h
        dinv = (h*h)/4.0
        jacobi_c[l] = omega * dinv
    #endfor

    impipeData['N'] = N

    impipeData['invhh']    = invhh
    impipeData['jacobi_c'] = jacobi_c

    # extent in each dimension
    extent = {}
    for l in range(0, L+1):
        extent[l] = Interval(Int, 0, N[l]+1)

    impipeData['extent'] = extent

    return

def setCases(impipeData, appData):
    y = impipeData['y']
    x = impipeData['x']

    N = impipeData['N']
    L = appData['L']

    interior = {}
    ghosts = {}
    for l in range(0, L+1):
        # grid interior
        interior[l] = {}

        interior[l]['inY'] = Condition(y, ">=", 1  ) \
                           & Condition(y, "<=", N[l])
        interior[l]['inX'] = Condition(x, ">=", 1  ) \
                           & Condition(x, "<=", N[l])
 
        interior[l]['innerBox'] = interior[l]['inY'] \
                                & interior[l]['inX']
 
        # grid ghosts
        ghosts[l] = {}
 
        # top and bottom planes
        ghosts[l]['ghostTop']    = Condition(y, "==", 0)
        ghosts[l]['ghostBottom'] = Condition(y, "==", N[l]+1)
 
        # left and right planes
        ghosts[l]['ghostLeft']   = Condition(x, "==", 0) \
                                 & interior[l]['inY']
        ghosts[l]['ghostRight']  = Condition(x, "==", N[l]+1) \
                                 & interior[l]['inY']

    impipeData['interior'] = interior
    impipeData['ghosts'] = ghosts

    return

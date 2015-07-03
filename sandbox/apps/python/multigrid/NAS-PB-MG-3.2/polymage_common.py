from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction

sys.path.insert(0, '../../../../optimizer')
sys.path.insert(0, '../../../../frontend')

from Compiler   import *
from Constructs import *

def setZeroGhosts(r, ghosts):
    for key in ghosts:
        r.defn.append(Case(ghosts[key], 0.0))

    return

def setVars(impipeDict, dataDict):
    lt = dataDict['lt']

    z = Variable(Int, "z")
    y = Variable(Int, "y")
    x = Variable(Int, "x")

    impipeDict['z'] = z
    impipeDict['y'] = y
    impipeDict['x'] = x

    n = Parameter(Int, "n")
    impipeDict['n'] = n

    N = {}
    N[lt] = n
    for l in range(lt-1, 0, -1):
        N[l] = (N[l+1]-2)/2 + 2
    impipeDict['N'] = N

    # extent in each dimension
    extent = {}
    for l in range(lt, 0, -1):
        extent[l] = Interval(Int, 0, N[l]-1)

    impipeDict['extent'] = extent

    return

def setCases(impipeDict, dataDict):
    z = impipeDict['z']
    y = impipeDict['y']
    x = impipeDict['x']

    N = impipeDict['N']
    lt = dataDict['lt']

    interior = {}
    ghosts = {}
    for l in range(lt, 0, -1):
        # grid interior
        interior[l] = {}

        interior[l]['inZ'] = Condition(z, ">=", 1  ) \
                           & Condition(z, "<=", N[l]-2)
        interior[l]['inY'] = Condition(y, ">=", 1  ) \
                           & Condition(y, "<=", N[l]-2)
        interior[l]['inX'] = Condition(x, ">=", 1  ) \
                           & Condition(x, "<=", N[l]-2)
 
        interior[l]['innerBox'] = interior[l]['inZ'] \
                                & interior[l]['inY'] \
                                & interior[l]['inX']
 
        # grid ghosts
        ghosts[l] = {}
 
        # front and back planes
        ghosts[l]['ghostFront']  = Condition(z, "==", 0)
        ghosts[l]['ghostBack']   = Condition(z, "==", N[l]-1)
 
        # top and bottom planes
        ghosts[l]['ghostTop']    = Condition(y, "==", 0) \
                                 & interior[l]['inZ']
        ghosts[l]['ghostBottom'] = Condition(y, "==", N[l]-1) \
                                 & interior[l]['inZ']
 
        # left and right planes
        ghosts[l]['ghostLeft']   = Condition(x, "==", 0) \
                                 & interior[l]['inY'] \
                                 & interior[l]['inZ']
        ghosts[l]['ghostRight']  = Condition(x, "==", N[l]-1) \
                                 & interior[l]['inY'] \
                                 & interior[l]['inZ']

    impipeDict['interior'] = interior
    impipeDict['ghosts'] = ghosts

    return

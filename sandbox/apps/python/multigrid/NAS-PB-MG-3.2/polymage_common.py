from __future__ import absolute_import, division, print_function

import sys
from fractions  import Fraction

sys.path.insert(0, '../../../../')
sys.path.insert(0, '../../../../')

from compiler import *
from constructs import *

def set_zero_ghosts(r, ghosts):
    for key in ghosts:
        r.defn.append(Case(ghosts[key], 0.0))

    return

def set_vars(app_data):
    pipe_data = app_data['pipe_data']

    lt = app_data['lt']

    z = Variable(Int, "z")
    y = Variable(Int, "y")
    x = Variable(Int, "x")

    pipe_data['z'] = z
    pipe_data['y'] = y
    pipe_data['x'] = x

    n = Parameter(Int, "n")
    pipe_data['n'] = n

    N = {}
    N[lt] = n
    for l in range(lt-1, 0, -1):
        N[l] = (N[l+1]-2)/2 + 2
    pipe_data['N'] = N

    # extent in each dimension
    extent = {}
    for l in range(lt, 0, -1):
        extent[l] = Interval(Int, 0, N[l]-1)

    pipe_data['extent'] = extent

    return

def set_cases(app_data):
    pipe_data = app_data['pipe_data']

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    N = pipe_data['N']
    lt = app_data['lt']

    interior = {}
    ghosts = {}
    for l in range(lt, 0, -1):
        # grid interior
        interior[l] = {}

        interior[l]['in_z'] = Condition(z, ">=", 1  ) \
                            & Condition(z, "<=", N[l]-2)
        interior[l]['in_y'] = Condition(y, ">=", 1  ) \
                            & Condition(y, "<=", N[l]-2)
        interior[l]['in_x'] = Condition(x, ">=", 1  ) \
                            & Condition(x, "<=", N[l]-2)

        interior[l]['inner_box'] = interior[l]['in_z'] \
                                 & interior[l]['in_y'] \
                                 & interior[l]['in_x']

        # grid ghosts
        ghosts[l] = {}

        # front and back planes
        ghosts[l]['ghost_front'] = Condition(z, "==", 0)
        ghosts[l]['ghost_back'] = Condition(z, "==", N[l]-1)

        # top and bottom planes
        ghosts[l]['ghost_top'] = Condition(y, "==", 0) \
                               & interior[l]['in_z']
        ghosts[l]['ghost_bottom'] = Condition(y, "==", N[l]-1) \
                                  & interior[l]['in_z']

        # left and right planes
        ghosts[l]['ghost_left'] = Condition(x, "==", 0) \
                                & interior[l]['in_y'] \
                                & interior[l]['in_z']
        ghosts[l]['ghost_right'] = Condition(x, "==", N[l]-1) \
                                 & interior[l]['in_y'] \
                                 & interior[l]['in_z']

    pipe_data['interior'] = interior
    pipe_data['ghosts'] = ghosts

    return

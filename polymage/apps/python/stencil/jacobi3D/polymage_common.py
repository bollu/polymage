from __init__ import *

import sys
sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

from fractions  import Fraction

def set_ghosts(r, ghosts, value):
    for ghost in ghosts:
        r.defn.append(Case(ghosts[ghost], value))

    return

def set_vars(app_data):
    pipe_data = app_data['pipe_data']

    z = Variable(Int, "z")
    y = Variable(Int, "y")
    x = Variable(Int, "x")

    pipe_data['z'] = z
    pipe_data['y'] = y
    pipe_data['x'] = x

    N = Parameter(Int, "N")
    pipe_data['N'] = N

    # jacobi weight (3D)
    omega = 8.0/9.0
    # omega.D^-1 = omega.(d^-1.I) = omega * h^2/6.0
    # d^-1 depends on diagonal elements of A^h
    h = 1.0/(N+1)
    invhh = 1.0/(h*h)
    dinv = (h*h)/6.0

    # mulitplier in the jacobi computation
    jacobi_c = omega * dinv

    pipe_data['invhh']    = invhh
    pipe_data['jacobi_c'] = jacobi_c

    # extent in each dimension
    extent = Interval(Int, 0, N+1)

    pipe_data['extent'] = extent

    return

def set_cases(app_data):
    pipe_data = app_data['pipe_data']
    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    N = pipe_data['N']

    # grid interior
    interior = {}

    interior['in_z'] = Condition(z, ">=", 1) \
                        & Condition(z, "<=", N)
    interior['in_y'] = Condition(y, ">=", 1) \
                        & Condition(y, "<=", N)
    interior['in_x'] = Condition(x, ">=", 1) \
                        & Condition(x, "<=", N)

    interior['inner_box'] = interior['in_z'] \
                             & interior['in_y'] \
                             & interior['in_x']

    # grid ghosts
    ghosts = {}

    # front and back planes
    ghosts['ghost_front']  = Condition(z, "==", 0)
    ghosts['ghost_back']   = Condition(z, "==", N+1)
 
    # top and bottom planes
    ghosts['ghost_top'] = Condition(y, "==", 0) \
                           & interior['in_z']
    ghosts['ghost_bottom'] = Condition(y, "==", N+1) \
                              & interior['in_z']

    # left and right planes
    ghosts['ghost_left'] = Condition(x, "==", 0) \
                            & interior['in_y'] \
                            & interior['in_z']
    ghosts['ghost_right'] = Condition(x, "==", N+1) \
                             & interior['in_y'] \
                             & interior['in_z']

    pipe_data['interior'] = interior
    pipe_data['ghosts'] = ghosts

    return

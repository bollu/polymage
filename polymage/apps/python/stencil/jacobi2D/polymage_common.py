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

    y = Variable(Int, "y")
    x = Variable(Int, "x")

    pipe_data['y'] = y
    pipe_data['x'] = x

    N = Parameter(Int, "N")
    pipe_data['N'] = N

    # jacobi weight (2D)
    omega = 4.0/5.0
    # omega.D^-1 = omega.(d^-1.I) = omega * h^2/4.0
    # d^-1 depends on diagonal elements of A^h
    h = 1.0/(N+1)
    invhh = 1.0/(h*h)
    dinv = (h*h)/4.0

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
    y = pipe_data['y']
    x = pipe_data['x']

    N = pipe_data['N']

    # grid interior
    interior = {}

    interior['in_y'] = Condition(y, ">=", 1) \
                        & Condition(y, "<=", N)
    interior['in_x'] = Condition(x, ">=", 1) \
                        & Condition(x, "<=", N)

    interior['inner_box'] = interior['in_y'] \
                          & interior['in_x']

    # grid ghosts
    ghosts = {}

    # top and bottom planes
    ghosts['ghost_top'] = Condition(y, "==", 0)
    ghosts['ghost_bottom'] = Condition(y, "==", N+1)
 
    # left and right planes
    ghosts['ghost_left'] = Condition(x, "==", 0) \
                         & interior['in_y']
    ghosts['ghost_right'] = Condition(x, "==", N+1) \
                          & interior['in_y']

    pipe_data['interior'] = interior
    pipe_data['ghosts'] = ghosts

    return

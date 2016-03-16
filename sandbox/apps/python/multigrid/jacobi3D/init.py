import sys
import os.path
import numpy as np
from arg_parser import parse_args
from printer import print_header, print_usage, print_line
from polymage_common import set_vars, set_cases
from exec_mg import calc_norm

def init_norm(app_data):
    grid_data = app_data['grid_data']
    U_ = grid_data['U_']

    app_data['resid'] = 0.0
    app_data['err']   = 0.0

    # calculate the initial residual norm and error
    print("[init]: calculating the initial norm and error ...")
    calc_norm(U_, app_data)
    print("[init]: ... DONE")

    return

def init_border(grid, border_width, border_values):
    # size of the grid - assumed to be same in all the dimensions
    n = grid.shape[0]

    w = border_width
    v = border_values

    # z-planes
    grid[  0:w] = v
    grid[n-w:n] = v

    # y-planes
    grid[0:n,   0:w] = v
    grid[0:n, n-w:n] = v

    # x-planes
    grid[0:n, 0:n,   0:w] = v
    grid[0:n, 0:n, n-w:n] = v

    return

def init_grids(app_data):
    print("[init_mg.py] : grids")

    N = app_data['N']

    # working grid (even step)
    U_ = np.ones((N+2, N+2, N+2), np.float64)
    # working grid (odd step)
    W_ = np.zeros((N+2, N+2, N+2), np.float64)

    init_border(U_, border_width=1, border_values=0.0)
    init_border(W_, border_width=1, border_values=0.0)
    # RHS
    F_ = np.zeros((N+2, N+2, N+2), np.float64)
    # exact solution
    U_EXACT_ = np.zeros((N+2, N+2, N+2), np.float64)

    grid_data = {}
    grid_data['U_'] = U_
    grid_data['W_'] = W_
    grid_data['F_'] = F_
    grid_data['U_EXACT_'] = U_EXACT_

    app_data['grid_data'] = grid_data

    return 

def init_params(app_data):
    app_args = app_data['app_args']

    print("[init_mg.py] : parameters")

    # size of each dimension of the coarsest grid
    n = int(app_args.n)
    # number of multigrid levels
    L = int(app_args.L)

    N = n
    # compute the size of the finest grid
    for l in range(0,L):
        N = 2*N+1

    app_data['n'] = n
    app_data['N'] = N
    app_data['L'] = L

    # pre-smoother, post-smoother and
    # coarse-grid relaxation steps
    app_data['nu1'] = int(app_args.nu1)
    app_data['nuc'] = int(app_args.nuc)
    app_data['nu2'] = int(app_args.nu2)

    assert not (app_data['nu1'] == 0 and \
                app_data['nu2'] == 0 and
                app_data['nuc'] == 0)

    return

def get_input(app_data):
    app_args = parse_args()
    app_data['app_args'] = app_args

    app_data['mode'] = app_args.mode
    app_data['cycle'] = app_args.cycle
    app_data['nit'] = int(app_args.nit)

    app_data['runs'] = int(app_args.runs)
    app_data['graph_gen'] = bool(app_args.graph_gen)

    cycle_name = app_data['cycle']+"cycle"
    app_data['cycle_name'] = cycle_name
    app_data['timer'] = app_args.timer

    # storage optimization
    app_data['optimize_storage'] = bool(app_args.optimize_storage)
    # early freeing of allocated arrays
    app_data['early_free'] = bool(app_args.early_free)
    # pool allocate option
    app_data['pool_alloc'] = bool(app_args.pool_alloc)

    return

def init_all(app_data):
    get_input(app_data)
    init_params(app_data)
    init_grids(app_data)

    pipe_data = {}
    app_data['pipe_data'] = pipe_data

    set_vars(app_data)
    set_cases(app_data)

    return

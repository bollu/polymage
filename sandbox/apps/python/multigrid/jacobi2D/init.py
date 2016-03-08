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

    # y-planes
    grid[  0:w] = v
    grid[n-w:n] = v

    # x-planes
    grid[0:n,   0:w] = v
    grid[0:n, n-w:n] = v

    return

def init_border_piecewise(grid, border_width, border_values):
    w = border_width
    v = border_values

    N = grid.shape[0]
    n = N-2*w+1

    # copy all the nine pieces
    grid[0:w, 0:w] = v[(0, 0)][0:w, 0:w]
    grid[0:w, w:n] = v[(0, 1)][0:w, w:n]
    grid[0:w, n:N] = v[(0, 2)][0:w, n:N]

    grid[w:n, 0:w] = v[(1, 0)][w:n, 0:w]
    #grid[w:n, w:n] = v[(1, 1)][w:n, w:n]
    grid[w:n, n:N] = v[(1, 2)][w:n, n:N]

    grid[n:N, 0:w] = v[(2, 0)][n:N, 0:w]
    grid[n:N, w:n] = v[(2, 1)][n:N, w:n]
    grid[n:N, n:N] = v[(2, 2)][n:N, n:N]

    return

def init_grids(app_args,app_data):
    print("[init_mg.py] : grids")

    N = app_data['N']
    problem = app_data['problem']

    # working grid (even step)
    U_ = np.ones((N+2, N+2), np.float64)
    # working grid (odd step)
    W_ = np.zeros((N+2, N+2), np.float64)

    if problem == 1:
        init_border(U_, border_width=1, border_values=0.0)
        init_border(W_, border_width=1, border_values=0.0)
        # RHS
        F_ = np.zeros((N+2, N+2), np.float64)
        # exact solution
        U_EXACT_ = np.zeros((N+2, N+2), np.float64)
    else:
        N1 = {}
        L = int(app_args.L) 
        n = int(app_args.n)
        for l in range(0,L+1):
            if l == 0:
                N1[0] = n
            else:
                N1[l] = 2*N1[l-1]+1

            h = 1.0/(N1[l]+1)

        indices = np.indices((N+2, N+2))
        x = indices[0] * h
        y = indices[1] * h

        x_1 = np.array(x)
        x_1 = x_1 + 1.0
        y_1 = np.array(y)
        y_1 = y_1 + 1.0

        x_y = x + y
        if problem == 2:
            b = {(0, 0):x, (0, 1):y  , (0, 2):x_1,
                 (1, 0):x, (1, 1):U_ , (1, 2):x_1,
                 (2, 0):x, (2, 1):y_1, (2, 2):x_1}
            init_border_piecewise(U_, border_width=1, border_values=b)
            init_border_piecewise(W_, border_width=1, border_values=b)
            # initialize F_
            F_ = np.zeros((N+2, N+2), np.float64) # f(i, j) <- 0
            # initialize U_
            U_EXACT_ = x_y
        elif problem == 3:
            sx = np.sin(x)
            sy = np.sin(y)

            sx_1 = np.sin(x_1)
            sy_1 = np.sin(y_1)

            sx_y = np.sin(x_y)

            b = {(0, 0):sx, (0, 1):sy  , (0, 2):sx_1,
                 (1, 0):sx, (1, 1):U_  , (1, 2):sx_1,
                 (2, 0):sx, (2, 1):sy_1, (2, 2):sx_1}
            init_border_piecewise(U_, border_width=1, border_values=b)
            init_border_piecewise(W_, border_width=1, border_values=b)
            # initialize F_
            F_ = -2.0 * sx_y
            # initialize U_EXACT_
            U_EXACT_ = sx_y
        elif problem == 4:
            xx = x*x
            yy = y*y

            xx_1 = xx+1.0
            yy_1 = yy+1.0

            xx_yy = xx + yy

            b = {(0, 0):xx, (0, 1):yy  , (0, 2):xx_1, \
                 (1, 0):xx, (1, 1):U_  , (1, 2):xx_1, \
                 (2, 0):xx, (2, 1):yy_1, (2, 2):xx_1}
            init_border_piecewise(U_, border_width=1, border_values=b)
            init_border_piecewise(W_, border_width=1, border_values=b)
            # initialize f
            F_ = np.zeros((N+2, N+2), np.float64) # f(i, j) <- 0
            # initialize u
            U_EXACT_ = xx_yy

    grid_data = {}
    grid_data['U_'] = U_
    grid_data['W_'] = W_
    grid_data['F_'] = F_
    grid_data['U_EXACT_'] = U_EXACT_

    app_data['grid_data'] = grid_data

    return 

def init_params(app_args, app_data):
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
    app_data['runs'] = int(app_args.runs)

    # pre-smoother, post-smoother and
    # coarse-grid relaxation steps
    app_data['nu1'] = int(app_args.nu1)
    app_data['nuc'] = int(app_args.nuc)
    app_data['nu2'] = int(app_args.nu2)

    # problem type
    app_data['problem'] = int(app_args.problem)

    # pool allocate option
    app_data['pool_alloc'] = app_args.pool_alloc

    assert not (app_data['nu1'] == 0 and \
                app_data['nu2'] == 0 and
                app_data['nuc'] == 0)

    return app_data

def get_input(app_args, app_data):
    app_data['app_args'] = app_args
    app_data['mode'] = app_args.mode
    app_data['cycle'] = app_args.cycle
    app_data['nit'] = int(app_args.nit)

    cycle_name = app_data['cycle']+"cycle"
    app_data['cycle_name'] = cycle_name
    app_data['timer'] = app_args.timer
  
    return

def init_all(impipe_data, app_data):
    # TODO init cycle type {V, W}
    app_args = parse_args()
    get_input(app_args,app_data)

    init_params(app_args,app_data)

    init_grids(app_args,app_data)

    set_vars(impipe_data, app_data)

    set_cases(impipe_data, app_data)

    return

import numpy as np

from verify  import set_verification
from misc import unpack_input, ilog2
from exec_mg import calc_norm
from polymage_common import set_cases, set_vars

def init_sizes(app_data):
    # init the problem size and other parameters and the respective solutions
    # based on the input class (from NAS PB)
    # TODO: generalize for custom specs

    # takes in an integer and returns its
    # logarithm to the base 2
    def ilog2(n):
        l = -1
        while n >= 1:
            l += 1
            n /= 2

        return l

    _class = app_data['prob_class']

    if _class == 'S':
        app_data['prob_size'] = 32
        app_data['nit'] = 4
    elif _class == 'W':
        app_data['prob_size'] = 128
        app_data['nit'] = 4
    elif _class == 'A':
        app_data['prob_size'] = 256
        app_data['nit'] = 4
    elif _class == 'B':
        app_data['prob_size'] = 256
        app_data['nit'] = 20
    elif _class == 'C':
        app_data['prob_size'] = 512
        app_data['nit'] = 20
    elif _class == 'D':
        app_data['prob_size'] = 1024
        app_data['nit'] = 50
    else:
        app_data['prob_size'] = 1
        app_data['nit'] = 1

    lt = ilog2(app_data['prob_size'])
    lb = 1

    # top and bottom levels
    # register in the data dictionary
    app_data['lt'] = lt
    app_data['lb'] = lb

    # grid sizes (with ghost) at each level
    N = {}
    N[lt] = app_data['prob_size']+2
    for l in range(lt-1, 0, -1):
        N[l] = (N[l+1]-2)/2+2

    # register in the data dictionary
    app_data['N'] = N

    return

def init_coefs(app_data):
    a = np.zeros((4), np.float64)
    c = np.zeros((4), np.float64)

    _class = app_data['prob_class']

    a[0] = -8.0/3.0
    a[1] =  0.0
    a[2] =  1.0/6.0
    a[3] =  1.0/12.0

    if _class == 'A' or _class == 'S' or _class == 'W':
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
    app_data['a'] = a
    app_data['c'] = c

    return

def make_periodic(V, app_data):
    N = app_data['N']
    lt = app_data['lt']

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

def unpack_input(index, file_name):
    f = open(file_name, 'r')

    i = 0
    for line in f:
        line_int = []
        line_int.append([int(x) for x in line.split()])
        for j in range(0, 3):
            index[i, j] = line_int[0][j]-1
        i += 1

    f.close()

    return

def init_grids(app_data):
    with_ghost = app_data['prob_size'] + 2
    _class = app_data['prob_class']

    # grids with ghost zone
    n1 = n2 = n3 = with_ghost
    V_ = np.zeros((n1, n2, n3), np.float64)
    U0_ = np.zeros((n1, n2, n3), np.float64)
    R0_ = np.zeros((n1, n2, n3), np.float64)

    U1_ = np.zeros((n1, n2, n3), np.float64)
    R1_ = np.zeros((n1, n2, n3), np.float64)

    # read the input file and unpack the co-ordinates
    index = np.zeros((20, 3), np.int16)
    unpack_input(index, "inputs/"+_class+".txt")

    # set the initial values at these co-ord.s
    for i in range(0, 10):
        V_[index[i,0], index[i,1], index[i,2]] = -1.0
    for i in range(10, 20):
        V_[index[i,0], index[i,1], index[i,2]] = 1.0

    #make_periodic(v, app_data)

    grid_data = {}

    # register in grid dictionary
    grid_data['V_'] = V_
    grid_data['U0_'] = U0_
    grid_data['R0_'] = R0_

    grid_data['U1_'] = U1_
    grid_data['R1_'] = R1_

    app_data['grid_data'] = grid_data

    return

def init_norm(app_data):
    rnm2 = 0.0
    rnmu = 0.0

    grid_data = app_data['grid_data']
    V_ = grid_data['V_']

    # load the norm computation library
    lib_file = "./norm.so"
    lib_func_name = "norm2u3"

    load_lib(lib_file, lib_func_name, app_data)

    # calculate norm on the initial grid
    calc_norm(V_, app_data)

    return

def init_polymage_data(app_data):
    # set variables, intervals etc
    set_vars(app_data)

    # set boundary case conditions
    set_cases(app_data)

    return

def get_input(app_data):
    app_args = parse_args()
    app_data['app_args'] = app_args

    app_data['mode'] = app_args.mode
    app_data['prob_class'] = app_args.prob_class

    return

def init_all(app_data):
    # collect the input data
    get_input(app_data)

    # initialize problem sizes and relevant params
    init_sizes(app_data)

    # initialize the stencil co-efficients
    init_coefs(app_data)

    # set the verification values
    set_verification(app_data, is_periodic=False)

    # initialize the grid contents
    init_grids(app_data)

    # initialize the norms
    init_norm(app_data)

    # initialize polymage specific parameters
    pipe_data = {}
    app_data['pipe_data'] = pipe_data
    init_polymage_data(app_data)

    return

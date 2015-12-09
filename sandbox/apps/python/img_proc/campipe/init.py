import sys
import os.path
from PIL import Image
import numpy as np
from arg_parser import parse_args

from printer import print_header, print_usage, print_line
#from polymage_common import set_vars, set_cases

def init_images(app_data):
    print("[init.py] : initializing images...")

    app_args = app_data['app_args']

    # input image: 
    img_path = app_args.img_file
    img = np.array(Image.open(img_path))
    rows, cols = img.shape
    # convert to short int image
    IN = np.array(img)
    IN = IN.astype(np.int16).ravel()

    # final output image
    OUT = np.zeros((3, rows-48, cols-64), np.uint8).ravel()

    # campipe parameter matrices
    matrix_3200 = \
        np.float32(
            [[ 1.6697, -0.2693, -0.4004, -42.4346],
             [-0.3576,  1.0615,  1.5949, -37.1158],
             [-0.2175, -1.8751,  6.9640, -26.6970]]
            )
    matrix_3200 = matrix_3200.ravel()

    matrix_7000 = \
        np.float32(
            [[ 2.2997, -0.4478,  0.1706, -39.0923],
             [-0.3826,  1.5906, -0.2080, -25.4311],
             [-0.0888, -0.7344,  2.2832, -20.0826]]
            )
    matrix_7000 = matrix_7000.ravel()

    img_data = {}
    img_data['IN'] = IN
    img_data['OUT'] = OUT
    img_data['M3200'] = matrix_3200
    img_data['M7000'] = matrix_7000

    app_data['img_data'] = img_data
    app_data['rows'] = rows
    app_data['cols'] = cols

    return

def get_input(sys_args, app_data):
    # parse the command-line arguments
    app_args = parse_args()
    app_data['app_args'] = app_args

    app_data['mode'] = app_args.mode
    app_data['pool_alloc'] = app_args.pool_alloc

    # flag to note execution time
    app_data['timer'] = os.path.isfile("timer.flag")

    return

def init_all(args, pipe_data, app_data):
    get_input(args, app_data)

    init_images(app_data)

    return

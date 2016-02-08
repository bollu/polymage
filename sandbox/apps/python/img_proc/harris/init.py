import sys
import os.path
from PIL import Image
import numpy as np
from arg_parser import parse_args

from printer import print_header, print_usage, print_line


def init_images(app_data):
    print("[init.py] : initializing images...")

    app_args = app_data['app_args']

    # input image: 
    img_path = app_args.img_file
    img = np.array(Image.open(img_path).convert('1'))
    rows, cols = img.shape
    # convert to short int image
    IN = np.array(img)
    IN = IN.astype(np.float32).ravel()

    # final output image
    OUT = np.zeros((3, rows, cols), ).ravel()

    img_data = {}
    img_data['IN'] = IN
    img_data['OUT'] = OUT

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

    return

def init_all(args, pipe_data, app_data):
    get_input(args, app_data)

    init_images(app_data)

    return


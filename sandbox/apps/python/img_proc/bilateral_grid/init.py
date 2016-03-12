import sys
import os.path
from PIL import Image
import numpy as np
from arg_parser import parse_args

from printer import print_header, print_usage, print_line
from utils import *

def init_images(app_data):
    print("[init.py] : initializing images...")

    app_args = app_data['app_args']

    off_left = int(app_args.off_left)
    total_pad = int(app_args.total_pad)

    # input image: 
    img_path = app_args.img_file
    img1 = Image.open(img_path)  #.convert('1'))
    img = np.array(img1)
    rows, cols, c = img.shape

    row_base = (rows-int(app_args.rowdiff))/2
    col_base = (cols-int(app_args.coldiff))/2
    image_region = img[row_base:row_base+rows, \
                     col_base:col_base+cols]
    
    # create ghost zone and copy image roi
    image_ghost = np.empty((rows+total_pad, cols+total_pad, 3), image_region.dtype)
    
    image_ghost[off_left:rows+off_left, off_left:cols+off_left, 0:3] = \
        np.array(image_region[0:rows, 0:cols, 0:3], image_region.dtype)
    
    # clamp the boundary portions
    image_clamp(image_region, image_ghost, \
            rows, cols, 3, \
            image_region.dtype, 1, \
            off_left, total_pad)

    image_rgb = Image.fromarray(image_ghost)
    image1 = image_rgb.convert('1')

    pixels = image1.getdata()
    data = np.reshape(pixels,image1.size)
    OUT = np.zeros((rows, cols), np.float32).ravel()
    img_data = {}
    img_data['IN'] = data
    img_data['OUT'] = OUT

    app_data['img_data'] = img_data
    app_data['rows'] = rows
    app_data['cols'] = cols

    app_data['rowdiff'] = int(app_args.rowdiff)
    app_data['coldiff'] = int(app_args.coldiff)

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


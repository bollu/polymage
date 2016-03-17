from __init__ import *

import sys
import os.path
from PIL import Image
import numpy as np
from arg_parser import parse_args

from printer import print_header, print_usage, print_line

sys.path.insert(0, ROOT)
from utils import *

def init_images(app_data):
    print("[init.py] : initializing images...")

    app_args = app_data['app_args']

    off_left = 24
    total_pad = 56

    # input image: 
    img_path = app_args.img_file
    img1 = Image.open(img_path)  #.convert('1'))
    img = np.array(img1)
    rows, cols, c = img.shape
  
    rowdiff = 2832
    coldiff = 4256

    row_base = (rows-rowdiff)/2
    col_base = (cols-coldiff)/2
    image_region = img[row_base:row_base+rows, \
                     col_base:col_base+cols]
    
    # create ghost zone and copy image roi
    image_ghost = \
        np.empty((rows+total_pad, cols+total_pad, 3), image_region.dtype)
    
    image_ghost[off_left:rows+off_left, off_left:cols+off_left, 0:3] = \
        np.array(image_region[0:rows, 0:cols, 0:3], image_region.dtype)
    
    # clamp the boundary portions
    image_clamp(image_region, image_ghost, rows, cols, 3, \
                image_region.dtype, 1, off_left, total_pad)

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

    app_data['rowdiff'] = rowdiff
    app_data['coldiff'] = coldiff

    return

def get_input(app_data):
    # parse the command-line arguments
    app_args = parse_args()
    app_data['app_args'] = app_args

    app_data['mode'] = app_args.mode
    app_data['runs'] = int(app_args.runs)
    app_data['graph_gen'] = bool(app_args.graph_gen)
    app_data['timer'] = app_args.timer

    # storage optimization
    app_data['optimize_storage'] = bool(app_args.optimize_storage)
    # early freeing of allocated arrays
    app_data['early_free'] = bool(app_args.early_free)
    # pool allocate option
    app_data['pool_alloc'] = bool(app_args.pool_alloc)

    return

def init_all(app_data):
    pipe_data = {}
    app_data['pipe_data'] = pipe_data

    get_input(app_data)

    init_images(app_data)

    return


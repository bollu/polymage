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
    rows = 2048
    cols = 2048

    # input image:
    img_path1 = app_args.img_file1
    img1 = np.array(Image.open(img_path1))

    img_path2 = app_args.img_file2
    img2 = np.array(Image.open(img_path2))

    if img1.shape != img2.shape:
        app_data['error'] = 1
        return 

    R = img1.shape[0]
    C = img1.shape[1]

    off_left = 31
    total_pad = 60

    # convert input image to floating point
    image1_f = np.float32(img1) / 255.0
    image2_f = np.float32(img2) / 255.0

    # create ghost zone and copy image roi

    # Image 1:
    image1_ghost = np.empty((rows+total_pad, cols+total_pad, 3), np.float32)
    image1_ghost[off_left:rows+off_left, off_left:cols+off_left, 0:3] = \
        np.array(img1[0:rows, 0:cols, 0:3], np.float32)
    # clamp the boundary portions
    image_clamp(img1, image1_ghost, rows, cols, 3,
                np.float32, 1, off_left, total_pad)

    # Image 2
    image2_ghost = np.empty((rows+total_pad, cols+total_pad, 3), np.float32)
    image2_ghost[off_left:rows+off_left, off_left:cols+off_left, 0:3] = \
        np.array(img2[0:rows, 0:cols, 0:3], np.float32)
    # clamp the boundary portions
    image_clamp(img2, image2_ghost, rows, cols, 3,
                np.float32, 1, off_left, total_pad)

    # create a simple mask of size (rows+total_pad) x (cols+total_pad)
    maskRow = 820
    half1 = np.zeros((maskRow, cols+total_pad), np.float32)
    half2 = np.ones((rows+total_pad-maskRow, cols+total_pad), np.float32)

    mask_ghost = np.vstack((half1, half2))

    # move colour dimension outside
    image1_f_flip = np.rollaxis(image1_ghost, 2).ravel()
    image2_f_flip = np.rollaxis(image2_ghost, 2).ravel()

    # result array
    OUT = np.empty((3, rows, cols), np.float32)            

    img_data = {}
    img_data['IN1'] = image1_f_flip
    img_data['IN2'] = image2_f_flip
    img_data['OUT'] = OUT
    img_data['mask_ghost'] = mask_ghost

    app_data['img_data'] = img_data
    app_data['rows'] = rows
    app_data['cols'] = cols
    app_data['total_pad'] = 60
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


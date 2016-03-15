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
    rows = 2560
    cols = 1536

    # input image: 
    img_path = app_args.img_file
    image = np.array(Image.open(img_path)) #.convert('1'))

    img_path2 = app_args.alpha_file
    alpha = np.array(Image.open(img_path2)) #.convert('1'))

    if image.shape[0] != alpha.shape[0] or \
               image.shape[1] != alpha.shape[1]:
        print("Please use image with the same shape")
        sys.exit(0)

    #rows, cols, c = img1.shape

    R = image.shape[0]
    C = image.shape[1]
    '''
    if R != 2560 or C != 1536:
        print("Please use 1536x2560 image size")
        sys.exit(0)
    '''
    image_flip = np.rollaxis(image, 2)

    # add alpha channel to image along with 
    # other colour channels
    imgalpha = np.append(image_flip, alpha)
    imgalpha = imgalpha.reshape(4, R, C)

    # get image roi
    row_base = (R - rows)/2
    col_base = (C - cols)/2

    imgalpha_region = imgalpha[0:4, \
                           row_base:row_base+rows+2, \
                           col_base:col_base+cols+2]

    # add ghost region
    imgalpha_ghost = \
        np.empty((4, rows+2, cols+2), np.float32)
    imgalpha_ghost[0:4, 1:rows+1, 1:cols+1] = \
        imgalpha_region

    # convert input image to floating point
    imgalpha_f = np.float32(imgalpha_ghost) / 255.0

    # result array
    res = np.empty((3, rows, cols), np.float32)
    

    img_data = {}
    img_data['IN'] = imgalpha_f
    img_data['OUT'] = res

    app_data['img_data'] = img_data
    app_data['R'] = R
    app_data['C'] = C
    app_data['rows'] = rows
    app_data['cols'] = cols
    return

def get_input(sys_args, app_data):
    # parse the command-line arguments
    app_args = parse_args()
    app_data['error'] = 0
    app_data['app_args'] = app_args

    app_data['mode'] = app_args.mode
    app_data['pool_alloc'] = app_args.pool_alloc

    return

def init_all(args, pipe_data, app_data):
    get_input(args, app_data)

    init_images(app_data)

    return


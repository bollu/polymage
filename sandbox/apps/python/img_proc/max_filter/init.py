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
    image = np.array(Image.open(img_path))

    #img_path2 = app_args.alpha_file
    #alpha = np.array(Image.open(img_path2))

    #if image.shape[0] != alpha.shape[0] or image.shape[1] != alpha.shape[1]:
      #  print("Please use alpha image with the same shape as the image")
      #  sys.exit(0)

    R = image.shape[0]
    C = image.shape[1]
    image_flip = np.rollaxis(image, 2)

    # add alpha channel to image along with other colour channels
    #imgalpha = np.append(image_flip, alpha)
    #imgalpha = imgalpha.reshape(4, R, C)

    #imgalpha_region = imgalpha[0:4, 0:R, 0:C]

    # add ghost region
    #imgalpha_ghost = np.empty((4, R+2, C+2), np.float32)
    #imgalpha_ghost[0:4, 1:R+1, 1:C+1] = imgalpha_region

    # convert input image to floating point
    #imgalpha_f = np.float32(imgalpha_ghost) / 255.0

    # result array
    res = np.empty((3, R, C), np.float32)

    img_data = {}
    img_data['IN'] = image
    #img_data['IN'] = imgalpha_f
    img_data['OUT'] = res

    app_data['img_data'] = img_data
    app_data['R'] = R
    app_data['C'] = C
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


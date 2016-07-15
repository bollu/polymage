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
	img_path_left = app_args.img_file_left
	img_left = np.array(Image.open(img_path_left))#.convert('1'))
	rows, cols, c = img_left.shape

	img_path_right = app_args.img_file_right
	img_right = np.array(Image.open(img_path_right))

	if img_right.shape != img_left.shape:
		app_data['error'] = 1
		print('Error! Image shapes are different.')
		return

	#row_base = int(app_args.rowbase)
	#col_base = int(app_args.colbase)

	#image_region = img[row_base:row_base+rows,col_base:col_base+cols]
	image_region_left = img_left[0:rows, 0:cols, 0:3]
	image_region_right = img_right[0:rows, 0:cols, 0:3]

	# create ghost zones
	image_ghost_left = np.zeros((rows+4, cols+4, 4), image_region_left.dtype)
	image_ghost_left[0:rows, 0:cols, 0:3] = image_region_left[0:rows, 0:cols, 0:3]
	image_ghost_right = np.zeros((rows+4, cols+4, 4), image_region_right.dtype)
	image_ghost_right[0:rows, 0:cols, 0:3] = image_region_right[0:rows, 0:cols, 0:3]

	# convert input image to floating point
	#image_f = np.float32(image_ghost) / 255.0

	# move colour dimension outside
	#image_f_flip = np.rollaxis(image_f, 2).ravel()
	#image_f_flip = image_f.ravel()
	image_f_left = np.rollaxis(image_ghost_left, 2).ravel()
	image_f_right = np.rollaxis(image_ghost_right, 2).ravel()

	OUT = np.zeros((4,rows, cols), np.uint8).ravel()

	#image_f = image_region.astype(np.float32).ravel()

	img_data = {}
	img_data['INL'] = image_f_left
	img_data['INR'] = image_f_right
	img_data['OUT'] = OUT

	app_data['img_data'] = img_data
	app_data['R'] = rows
	app_data['C'] = cols

	#app_data['rowbase'] = int(app_args.rowbase)
	#app_data['colbase'] = int(app_args.colbase)

	return

def get_input(app_data):
	app_args = parse_args()
	app_data['app_args'] = app_args

	app_data['mode'] = app_args.mode
	#app_data['threshold'] = float(app_args.threshold)
	#app_data['weight'] = int(app_args.weight)

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


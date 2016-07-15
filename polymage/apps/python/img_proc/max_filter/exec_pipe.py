import sys
import os
import ctypes
import numpy as np
import time

from printer import print_line

from compiler   import *
from constructs import *
from utils import *

from PIL import Image

def call_pipe(app_data):
	rows = app_data['R']
	cols = app_data['C']

	img_data = app_data['img_data']
	IN = img_data['IN']
	OUT = img_data['OUT']

	# lib function name
	func_name = 'pipeline_'+app_data['app']
	pipe_func = app_data[func_name]

	# lib function args
	pipe_args = []
	pipe_args += [ctypes.c_int(cols)]
	pipe_args += [ctypes.c_int(rows)]
	pipe_args += [ctypes.c_void_p(IN.ctypes.data)]
	pipe_args += [ctypes.c_void_p(OUT.ctypes.data)]

	# call lib function
	pipe_func(*pipe_args)
	
	return

def maxfilter(app_data):

	it  = 0
	app_args = app_data['app_args']

	img_path = app_args.img_file
	img = Image.open(img_path)
	#img.show()
	arr1 = np.array(img)
	print("Image dtype: "+str(arr1.dtype))
	mode = img.mode
	size = img.size

	rows = app_data['R']
	cols = app_data['C']
	in1 = app_data['img_data']['IN']
	in1 = np.reshape(in1, (3, rows, cols))
	in1 = np.rollaxis(in1, 2)
	in1 = np.rollaxis(in1, 2)
	print("Input Image shape: "+str(in1.shape))
	in1 = np.uint8(in1).ravel()
	img_in = Image.frombuffer(mode,size,in1)
	#img_in.show()
	img_in.save("in.jpeg","JPEG")
   
	runs = int(app_args.runs)
	timer = app_args.timer
	if timer == True:
		t1 = time.time()

	while it < runs :
		call_pipe(app_data)
		it += 1

		OUT = app_data['img_data']['OUT']
		OUT = OUT.reshape(4, rows, cols)
		OUT= np.rollaxis(OUT, 2)
		OUT= np.rollaxis(OUT, 2)
		print("Output Image shape: "+str(OUT.shape))
		out1 = np.zeros((rows,cols,3),np.uint8)
		out1[0:rows,0:cols,0:3] = OUT[0:rows,0:cols,0:3]
		out1.ravel()
		img_out = Image.frombuffer(mode,size,out1)
		#img_out.show()
		img_out.save('out.jpeg',"JPEG")

	if timer == True:
		t2 = time.time()

		time_taken = float(t2) - float(t1)
		print("")
		print("[exec_pipe] : time taken to execute = ",
			  (time_taken * 1000) / runs, " ms")

	return


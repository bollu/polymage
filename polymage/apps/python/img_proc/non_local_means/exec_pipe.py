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

	app_args = app_data['app_args']
	patch_size = int(app_args.patch_size)
	search_area = int(app_args.search_area)

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
	pipe_args += [ctypes.c_int(patch_size)]
	pipe_args += [ctypes.c_int(search_area)]
	pipe_args += [ctypes.c_void_p(IN.ctypes.data)]
	pipe_args += [ctypes.c_void_p(OUT.ctypes.data)]

	# call lib function
	pipe_func(*pipe_args)
	
	return

def nlmeans(app_data):

	it  = 0
	app_args = app_data['app_args']

	img_path = app_args.img_file
	img = Image.open(img_path)
	mode = img.mode
	print("mode = "+str(mode))
	size = img.size
	#img.show()
	arr1 = np.array(img)
	#imgIN = Image.frombuffer(mode,size,arr1.ravel())
	#imgIN.show()
	print ("size = "+str(size))
	print("dtype = "+str(arr1.dtype))

	rows = app_data['R']
	cols = app_data['C']
	print("rows = "+str(rows))
	print("cols = "+str(cols))
	in1 = app_data['img_data']['IN']
	#print("len in1 = "+str(len(in1)))

	#in1 = np.reshape(in1, (4,rows+4,cols+4))
	in1 = np.reshape(in1, (4,rows,cols))
	print("in1.shape = "+str(in1.shape))
	in1 = np.rollaxis(in1, 2)
	in1 = np.rollaxis(in1, 2)
	#in1 = np.rollaxis(in1, 1, 2)
	print("in1.shape = "+str(in1.shape))
	in_ghost = np.zeros((rows,cols,3),in1.dtype)
	in_ghost[0:rows,0:cols,0:3] = in1[0:rows,0:cols,0:3]
	in_ghost = in_ghost.ravel()
	in_ghost = np.uint8(in_ghost * 255.0)
	#in1 = in1.ravel()
	#in1 = np.int32(in1 * 255.0)

	img2 = Image.frombuffer(mode,size,in_ghost)
	#img2 = Image.frombuffer(mode,size,in1)
	#np.savetxt('validate/poly_in.txt', in_ghost)
	img2.show()
   
	runs = int(app_args.runs)
	timer = app_args.timer
	if timer == True:
		t1 = time.time()

	b4 = app_data['img_data']['OUT']
	#np.savetxt('validate/poly_out.tmp.txt', b4)
	while it < runs :
		call_pipe(app_data)
		it += 1
		OUT = app_data['img_data']['OUT']
		#np.savetxt('validate/poly_out.txt', OUT)
		OUT = OUT.reshape(3,rows,cols)
		OUT = np.rollaxis(OUT, 2)
		OUT = np.rollaxis(OUT, 2)
		OUT = np.uint8(OUT * 255.0)
		print("OUT SHAPE: "+str(OUT.shape))
		out2 = np.zeros((rows,cols,3),np.uint8)
		out2[0:rows,0:cols,0:3] = OUT[0:rows,0:cols,0:3]
		out2 = out2.ravel()
		OUT=OUT.ravel()
		#with open("out"+str(it)+".dat",'w') as fil:
			#fil.write(str(OUT))
		img1 = Image.frombuffer(mode,size,out2)
		#img1 = Image.frombuffer('RGBX',size,OUT)
		img1.show()

	if timer == True:
		t2 = time.time()

		time_taken = float(t2) - float(t1)
		print("")
		print("[exec_pipe] : time taken to execute = ",
			  (time_taken * 1000) / runs, " ms")

	return


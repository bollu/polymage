from __init__ import *

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def nlmeans(pipe_data):
	# Params
	R = Parameter(Int, "R") # image rows
	C = Parameter(Int, "C") # image cols
	patch_size = Parameter(Int, "patch_size")
	search_area = Parameter(Int, "search_area")

	sigma = 0.12

	pipe_data['R'] = R
	pipe_data['C'] = C

	# Vars
	x = Variable(Int, "x")
	y = Variable(Int, "y")
	z = Variable(Int, "z")
	c = Variable(Int, "c")

	# Input Image
	img = Image(Float, "img", [3, R+2, C+2])

	# clamped image - boundary conditions: repeat edge
	rows = Interval(Int, 0, R-1)
	cols = Interval(Int, 0, C-1)
	colours = Interval(Int, 0, 3)
	s_dom_x = Interval(Int, -search_area/2, search_area + 1)
	s_dom_y = Interval(Int, -search_area/2, search_area + 1)

	#clamped = Function(([c, x, y],[colours, rows, cols]), Float, "clamped")
	#xminR = Min(x, R-1)
	#xclamp = Max(xminR, 0)
	#yminC = Min(y, C-1)
	#yclamp = Max(0, yminC)
	#clamped.defn = [ img(c,) ]
	#clamped = img(c, xclamp, yclamp)

	def clamped(k, i, j):
		xminR = Min(i, R-1)
		xclamp = Max(xminR, 0)
		yminC = Min(j, C-1)
		yclamp = Max(yminC, 0)
		return img(k, xclamp, yclamp)

	inv_sigma_sq = -1.0/(sigma*sigma*patch_size*patch_size)

	# Halide: Define the difference Images (Func dc)
	dx = Variable(Int, "dx")
	dy = Variable(Int, "dy")
	
	## Body case
	#case_dc = Condition(x, ">=", 1) & \
			#Condition(x + dx, ">=", 1) & \
			#Condition(x, "<=", R) & \
			#Condition(x + dx, "<=", R) & \
			#Condition(y, ">=", 1) & \
			#Condition(y+dy, ">=", 1) & \
			#Condition(y, "<=", C) & \
			#Condition(y+dy, "<=", C)

	dc = Function(([c, x, y, dx, dy], [colours, rows, cols, s_dom_x, s_dom_y]), Float, "dc")
	dc.defn = [ Powf((clamped(c,x,y) - clamped(c,x+dx,y+dy)),2) ]

	# Halide: Sum across colour channels (Func d)
#	channels = Interval(Int, 0, 3)
#	channels_red = Reduction(([c, x, y, dx, dy], [colours, rows, cols, s_dom_x, s_dom_y]), \
#			([c], [channels]), \
#			Float, \
#			"channels")
#	d_sum = Reduce(channels_red(c, x, y, dx, dy), \
#			dc(c, x, y, dx, dy), \
#			Op.Sum)
	#d = Function(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]), Float, "d")
	#d.defn = [ d_sum ]
	d = Reduction(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]), \
			([c, x, y, dx, dy], [colours, rows, cols, s_dom_x, s_dom_y]), \
			Float, \
			"d")
	d.defn = [ Reduce(d(x, y, dx, dy), dc(c, x, y, dx, dy), Op.Sum) ]

	# Halide: Find the patch differences by blurring the difference image
	# Function blur_d
	patch_var = Variable(Int, "patch_var")
	patch_interval = Interval(Int, -patch_size/2, patch_size + 1)
	y_patch_int = Interval(Int, patch_size/2, C - patch_size - 2)
	x_patch_int = Interval(Int, patch_size/2, R - patch_size - 2)
	blur_d_y = Reduction(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]), \
			([x, y, dx, dy, patch_var], [rows, y_patch_int, s_dom_x, s_dom_y, patch_interval]), \
			Float, \
			"blur_d_y")
	blur_d_y.defn = [ Reduce(blur_d_y(x, y, dx, dy), \
			d(x, y + patch_var, dx, dy), \
			Op.Sum) ]
	#blur_d_y = Function(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]), Float, "blur_d_y")
	#blur_d_y.defn = [ blur_d_y_sum ]

	#blur_d_sum = Reduce(patch_dom(x, y, dx, dy), \
			#		blur_d_y(x + patch_var, y, dx, dy), \
			#Op.Sum)
	#blur_d = Function(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]) , Float, "blur_d")
	#blur_d.defn = [ blur_d_sum ]
	blur_d = Reduction(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]), \
			([x, y, dx, dy, patch_var], [x_patch_int, cols, s_dom_x, s_dom_y, patch_interval]), \
			Float, \
			"blur_d")
	blur_d.defn = [ Reduce(blur_d(x, y, dx, dy), \
			blur_d_y(x+ patch_var, y, dx, dy), \
			Op.Sum) ]

	# Halide: Compute the weights from the patch differences.
	w = Function(([x, y, dx, dy], [rows, cols, s_dom_x, s_dom_y]), \
			Float, \
			"w")
	w.defn = [ Exp(blur_d(x, y, dx, dy) * inv_sigma_sq) ]
	#w.defn = Cast(Float, blur_d(x, y, dx, dy) * inv_sigma_sq)

	# Halide: Add an alpha channel. Func clamped_with_alpha
	#cond_alpha = Condition(x, ">=", 1) & Condition(x, "<=", R) & \
			#		Condition(y, ">=", 1) & Condition(y, "<=", C) & \
			#Condition(c, ">=", 0) & Condition(c, "<=", 2)
	#ca = Variable(Int, "ca")
	cond_alpha = Condition(c, ">=", 0) & Condition(c, "<=", 2)

	clamped_with_alpha = Function(([c, x, y], [colours, rows, cols]), \
			Float, \
			"clamped_with_alpha")
	clamped_with_alpha.defn = [Select(cond_alpha, clamped(c, x, y), 1.0)]

	# Halide: Define a reduction domain for the search area.
	x_search_int = Interval(Int, search_area/2, R - search_area -2)
	y_search_int = Interval(Int, search_area/2, C - search_area -2)
	non_local_means_sum = Reduction(([c, x, y], [colours, rows,cols]), \
			([c, x, y, dx, dy], [colours, x_search_int, y_search_int, s_dom_x, s_dom_y]), \
			Float,
			"non_local_means_sum")
	# Halide: Compute the sum of the pixels in the search area.
	# Func non_local_means_sum
	non_local_means_sum.defn = [ Reduce(non_local_means_sum(c, x, y), \
			w(x, y, dx, dy) * clamped_with_alpha(c, x + dx, y + dy), \
			Op.Sum) ]
	#non_local_means_sum = Function(([c, x, y], [colours, rows, cols]), Float, "non_local_means_sum")
	#non_local_means_sum.defn = [ s_dom_red ]

	# Final function: non_local_means
	nlm_sum = non_local_means_sum(c, x, y) / non_local_means_sum(3, x, y)
	#condf = Condition(nlm_sum, ">=", 0.0) & Condition(nlm_sum, "<=", 1.0)
	non_local_means = Function(([c, x, y], [colours, rows, cols]), Float, "non_local_means")
	#non_local_means.defn = [ Select(condf, nlm_sum, 0.0) ]
	nlm_min = Min(nlm_sum, 1.0)
	nlm_clamp = Max(nlm_min, 0.0)
	non_local_means.defn = [ nlm_clamp ]

	return non_local_means


from __init__ import *

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def lensblur(pipe_data):
	# Params
	R = Parameter(Int, "R") # image rows
	C = Parameter(Int, "C") # image cols

	pipe_data['R'] = R
	pipe_data['C'] = C

	# Vars
	x = Variable(Int, "x")
	y = Variable(Int, "y")
	z = Variable(Int, "z")
	c = Variable(Int, "c")

	# Input Image
	img_left = Image(Int, "imgl", [3, R+2, C+2])
	img_right = Image(Int, "imgr", [3, R+2, C+2])

	# clamped image - boundary conditions: repeat edge
	rows = Interval(Int, 0, R-1)
	cols = Interval(Int, 0, C-1)
	channels = Interval(Int, 0, 2)

	def clamp(e, mn, mx):
		min1 = Min(e, mx)
		max1 = Max(min1, mn)
		return max1

	def right(k, i, j):
		return img_right(k, clamp(i, 0, R-1), clamp(j, 0, C-1))

	def left(k, i, j):
		return img_left(k, clamp(i, 0, R-1), clamp(j, 0, C-1))


	#Halide: The number of displacements to consider
	slices = 32

	#Halide: The depth to focus on
	focus_depth = 13

	#Halide: The increase in blur radius with misfocus depth
	blur_radius_scale = float(0.5)

	#Halide: The number of samples of the aperture to use
	aperture_samples = 32

	z_int = Interval(Int, 0, aperture_samples)

	maximum_blur_radius = int(max(slices - focus_depth, focus_depth) * blur_radius_scale)

	def absd(a, b):
		mn = Min(a, b)
		mx = Max(a, b)
		return mx - mn

	diff = Function(([c,x,y,z],[channels,rows,cols,z_int]), Int, "diff")
	diff.defn = [ Min(absd(left(c, x, y), right(c, x + 2*z, y)),
		absd(left(c, x, y), right(c, x + 2*z + 1, y))) ]

	cost = Function(([x, y, z], [rows, cols, z_int]), Float, "cost")
	cost.defn = [ Powf(Cast(Float, diff(0,x,y,z)), 2) + \
		Powf(Cast(Float, diff(0,x,y,z)), 2) + \
		Powf(Cast(Float, diff(0,x,y,z)), 2) ]
	
	#Halide: Compute the confidence of cost estimate at each pixel by
	# taking the variance across the stack
	r = Variable(Int, "r")
	r_int = Interval(Int, 0, slices)
	#cost_confidence = Reduction(([x, y], [rows, cols]), 
			#([x, y, r], [rows, cols, r_int]), 
			#Float, "cost_confidence")
	#a = Reduce(cost_confidence(x,y), Powf(cost(x, y, r), 2), Op.Sum) #TODO: Fix this...
	#a1 = Function(([x,y],[rows,cols]) , Float, "expr_a1")
	a1 = Reduction(([x, y], [rows, cols]), 
			([x, y, r], [rows, cols, r_int]), 
			Float, "a1")
	a1.defn = [ Reduce(a1(x,y), Powf(cost(x, y, r), 2), Op.Sum) ] 
	a = Function(([x,y],[rows,cols]) , Float, "expr_a")
	a.defn = [ a1(x,y)/slices ]
	#b1 = Function(([x,y],[rows,cols]) , Float, "expr_b1")
	b1 = Reduction(([x, y], [rows, cols]), 
			([x, y, r], [rows, cols, r_int]), 
			Float, "b1")
	b1.defn = [ Reduce(b1(x, y), cost(x, y, r) / slices, Op.Sum) ]
	b = Function(([x,y],[rows,cols]) , Float, "expr_b")
	b.defn = [ Powf(b1(x,y), 2) ]
	cost_confidence = Function(([x,y],[rows,cols]), Float, "cost_confidence")
	cost_confidence.defn = [ a(x,y) - b(x,y) ]

	#Halide: Do a push-pull thing to blur the cost volume with
	# an exponential decay type thing to inpaint over regions with
	# low confidence
	cost_pyramid_push = []
	cost_pyramid_push.append(Function(([c, x, y, z], [channels, rows, cols, z_int]), Float, "cost_pyramid_push"+str(0)))
	cost_pyramid_push[0].defn = [ Select(Condition(c, "==", 0), cost(x,y,z) * cost_confidence(x,y), cost_confidence(x, y)) ]
	
	#Halide: commented - w = left_im.width(), h = left_im.height();
	def downsample_and_edge(f, s, i, w, h):
		nf = Function(([c, x, y, z], [channels, rows, cols, z_int]), Float, s+"_"+str(i))
		#x_downx_int = Interval(Int, 1, R/2)
		downx = Function(([c, x, y, z], [channels, rows, cols, z_int]), Float, "downx"+str(i))
		downx.defn = [ (f(c, 2*x - 1, y, z) + \
				float(3.0) * (f(c, 2*x, y, z) + \
				f(c, 2*x+1, y, z)) + \
				f(c, 2*x+2, y, z)) / float(8.0) ]
		downy = Function(([c,x,y,z], [channels,rows,cols,z_int]), Float, "downy"+str(i))
		downy.defn = [ (downx(c, x, 2*y-1, z) + \
				float(3.0) * (downx(c, x, 2*y, z) + \
				downx(c, x, 2*y+1, z)) + \
				downx(c, x, 2*y+2, z)) / float(8.0)  ]
		cx = clamp(x,0,int(w))
		cy = clamp(y,0,int(h)) #TODO: check if x,y and h,w need to be reversed
		nf.defn = [ downy(c, cx, cy, z) ]
		#cost_pyramid_push[i] = nf
		return nf

	w = 992
	h = 1024
	for i in range(1,8):
		cost_pyramid_push.append(downsample_and_edge(cost_pyramid_push[i-1],"cost_pyramid_push",i,w,h))
		w /= 2
		h /= 2


	def upsample(f,s,i):
		nf = Function(([c, x, y, z], [channels, rows, cols, z_int]), Float, s+"_"+str(i))

		upx = Function(([c, x, y, z], [channels, rows, cols, z_int]), Float, "upx"+str(i))
		upx.defn = [ float(0.25) * f(c, (x/2) - 1 + 2*(x % 2), y, z) + \
				float(0.75) * f(c, x/2, y, z) ]
		upy = Function(([c, x, y, z], [channels, rows, cols, z_int]), Float, "upy"+str(i))
		upy.defn = [ float(0.25) * upx(c, x, (y/2) - 1 + 2*(y % 2), z) + \
				float(0.75) * upx(c, x, y/2, z) ]

		nf.defn = [ upy(c,x,y,z) ]

		return nf

	def lerp(a,b,w):
		return a + w * (b - a)

	cost_pyramid_pull = [ None for _ in range(8)]
	cost_pyramid_pull[7] = Function(([c,x,y,z], [channels, rows, cols, z_int]), Float, "cost_pyramid_pull"+str(7))
	cost_pyramid_pull[7].defn = [ cost_pyramid_push[7](c,x,y,z) ]
	for i in reversed(range(7)):
		cost_pyramid_pull[i] = Function(([c,x,y,z], [channels, rows, cols, z_int]), Float, "cost_pyramid_pull"+str(i))
		#cost_pyramid_pull[i].defn = lerp(upsample())
		us = upsample(cost_pyramid_pull[i+1], "cost_pyramid_pull", i)
		cost_pyramid_pull[i].defn = [ lerp(us(c,x,y,z), cost_pyramid_push[i](c,x,y,z), float(0.5)) ]

	filtered_cost = Function(([x, y, z], [rows, cols, z_int]), Float, "filtered_cost")
	filtered_cost.defn = [ cost_pyramid_pull[0](0,x,y,z) / cost_pyramid_pull[0](1,x,y,z) ]

	#Halide: Assume minimum cost slice is the correct depth
	depth = Function(([x,y], [rows,cols]), Int, "depth")
	depth.defn = [ cost_confidence(x, y) ] #TODO: place holder. fix this
	#depth.defn = #TODO: Implement this => depth(x, y) = argmin(filtered_cost(x, y, r))[0]

	bokeh_radius = Function(([x, y], [rows, cols]), Float, "bokeh_radius")
	bokeh_radius.defn = [ Cast(Float, Abs(depth(x, y) - focus_depth)) * blur_radius_scale ]

	bokeh_radius_squared = Function(([x, y], [rows, cols]), Float, "bokeh_radius_squared")
	bokeh_radius_squared.defn = [ Powf(bokeh_radius(x, y), 2) ]

	#Halide: Take a max filter of the bokeh radius to determine the worst case bokeh radius
	# to consider at each pixel. Makes the sampling more efficient below.
	rw = Variable(Int, "rw")
	rw_int = Interval(Int, -maximum_blur_radius, 2 * maximum_blur_radius + 1)
	y_rw_int = Interval(Int, maximum_blur_radius, C - 2*maximum_blur_radius -2)
	x_rw_int = Interval(Int, maximum_blur_radius, R - 2*maximum_blur_radius -2)
	worst_case_bokeh_radius_y = Reduction(([x, y], [rows, cols]), 
			([x, y, rw], [rows, y_rw_int, rw_int]),
			Int, "worst_case_bokeh_radius_y")
	worst_case_bokeh_radius_y.defn = [ Reduce(worst_case_bokeh_radius_y(x,y), bokeh_radius(x, y+rw), Op.Max) ]
	worst_case_bokeh_radius = Reduction(([x, y], [rows, cols]), 
			([x, y, rw], [x_rw_int, cols, rw_int]),
			Int, "worst_case_bokeh_radius")
	worst_case_bokeh_radius.defn = [ Reduce(worst_case_bokeh_radius(x,y), bokeh_radius(x+rw, y), Op.Max) ]

	cond_alpha = Condition(c, ">=", 0) & Condition(c, "<=", 2)
	input_with_alpha = Function(([c, x, y], [Interval(Int, 0, 3), rows, cols]), Float, "input_with_alpha")
	input_with_alpha.defn = [ Select(cond_alpha, Cast(Float, left(c,x,y)), float(255.0)) ]

	#Halide: Render a blurred image
	output = Function(([c,x,y], [channels,rows,cols]), Float, "output")
	output.defn = [ input_with_alpha(c,x,y) ]

	def random_float():
		return RandomFloat()
	
	#Halide: The sample locations are a random function of x, y, and
	# sample number (not c).
	worst_radius = worst_case_bokeh_radius(x, y)
	sample_u = (random_float() - float(0.5)) * 2 * worst_radius
	sample_v = (random_float() - float(0.5)) * 2 * worst_radius
	sample_u = clamp(Cast(Int, sample_u), -maximum_blur_radius, maximum_blur_radius)
	sample_v = clamp(Cast(Int, sample_v), -maximum_blur_radius, maximum_blur_radius)

	sample_locations_u = Function(([x,y,z], [rows, cols, z_int]), Int, "sample_locations_u")
	sample_locations_u.defn = [ sample_u ]
	sample_locations_v = Function(([x,y,z], [rows, cols, z_int]), Int, "sample_locations_v")
	sample_locations_v.defn = [ sample_v ]

	s = Variable(Int,"s")
	s_int = Interval(Int, 0, aperture_samples)
	sample_u = sample_locations_u(x,y,z)
	sample_v = sample_locations_v(x,y,z)
	sample_x = x + sample_u
	sample_y = y + sample_v
	r_squared = sample_u * sample_u + sample_v * sample_v


	#Halide: We use this sample if it's from a pixel whose bokeh
	# influences this outpus pixel. Here's a crude approx. that
	# ignores some subtleties of occlusion edges and inpaints
	# behind objects.
	sample_is_within_bokeh_of_this_pixel = Condition(r_squared, 
			"<", bokeh_radius_squared(x,y))

	this_pixel_is_within_bokeh_of_sample = Condition(r_squared, 
			"<", bokeh_radius_squared(sample_x, sample_y))

	sample_is_in_front_of_this_pixel = Condition(depth(sample_x, sample_y), 
			"<", depth(x,y))

	sample_weight = Function(([x,y,z],[rows,cols,z_int]) , Float, "sample_weight")
	sample_weight.defn = [ Select((sample_is_within_bokeh_of_this_pixel | sample_is_in_front_of_this_pixel) & this_pixel_is_within_bokeh_of_sample, float(1.0), float(0.0)) ]

	sample_x = sample_locations_u(x,y,s)
	sample_y = sample_locations_v(x,y,s)
	output = Reduction(([c,x,y],[Interval(Int,0,3),rows,cols]),([c,x,y,s],[channels,rows,cols,s_int]), Float, "output")
	output.defn = [ Reduce(output(c,x,y),
		sample_weight(x,y,s) * input_with_alpha(c,sample_x,sample_y),
		Op.Sum) ]

	final = Function(([c,x,y],[channels,rows,cols]), Float, "final")
	final.defn = [ output(c,x,y) / output(3,x,y) ]

	return final



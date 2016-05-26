from __init__ import *

import sys
import subprocess
import numpy as np
from fractions import Fraction
import math

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def maxfilter(pipe_data):
	# Pipeline Variables
	x = Variable(Int, "x")
	y = Variable(Int, "y")
	c = Variable(Int, "c")
	t = Variable(Int, "t")

	# Pipeline Parameters
	R = Parameter(Int, "R") # image rows
	C = Parameter(Int, "C") # image cols

	# Register in the dictionary
	pipe_data['R'] = R
	pipe_data['C'] = C

	# Input Image
	img = Image(Float, "img", [3, R, C])

	radius = 26

	slices = int(math.ceil(math.log(radius,2))) + 1

	rows = Interval(Int, 0, R-1)
	cols = Interval(Int, 0, C-1)
	c_int = Interval(Int, 0, 2)
	#channels = Interval(Int, 0, 2)

	t_int = Interval(Int, 0, slices)

	def clamped(i,j,k):
		xminR = Min(i, R-1)
		xclamp = Max(0, xminR)
		yminC = Min(j, C-1)
		yclamp = Max(0, yminC)
		return img(k,xclamp,yclamp)

	def clamp(e, mn, mx):
		minm = Min(e, mx)
		maxm = Max(minm, mn)
		return maxm

	rx = Variable(Int, "rx")
	ry = Variable(Int, "ry")
	rx_int = Interval(Int, -radius, R + radius) #TODO: re-check this
	ry_int = Interval(Int, 1, slices - 1)

	y_ry_int = Interval(Int, radius, C - R -radius*3 - 1)

	vert_log = Reduction(([x, y, c, t], [rows, cols, c_int, t_int]), 
			([x, rx, c, ry], [rows, y_ry_int, c_int, ry_int]), Int, "vert_log")
	vert_log.defn = [ Reduce(vert_log(x,y,c,t), 
		Max(clamped(x, rx, c), 
			#clamped(x, rx + clamp(((ry-1)),0,radius*2), c) ),	   #TODO: fix this
			#clamped(x, rx + clamp((1<<(ry-1)),0,radius*2), c) ),   #TODO: to this
			clamped(x, rx + clamp(Cast(Int,Pow(2,ry-1)), 0, radius*2),c)), #TODO: This is a last resort
		Op.Max) ]

	slice_for_radius = Function(([t],[t_int]), Int, "slice_for_radius")
	slice_for_radius.defn = [ Cast(Int, Log(2*t + 1)/float(0.693147)) ]

	y_vert_int = Interval(Int, 0, C - slices - 1 - radius*2)
	vert = Function(([x, y, c, t], [rows, cols, c_int, t_int]), Int, "vert")
	eslice = clamp(slice_for_radius(t), 0, slices)
	first_sample = vert_log(x, y-t, c, eslice)
	#second_sample = vert_log(x, y + t + 1 - clamp( eslice, 0, 2*radius), c, eslice)		 #TODO: Fix this
	#second_sample = vert_log(x, y + t + 1 - clamp(1 << eslice, 0, 2*radius), c, eslice)	#TODO: to this
	second_sample = vert_log(x, y + t + 1 - clamp( Cast(Int,Pow(2,eslice)), 0, 2*radius), c, eslice)		 #TODO: Last Resort
	vert.defn = [ Max(first_sample, second_sample) ]

	dx = Variable(Int, "dx") # for final
	dx_int = Interval(Int, -radius, 2*radius+1)

	dy = Variable(Int, "dy")
	dy_int = Interval(Int, 0, radius + 1)
	x_radius_int = Interval(Int, 0, radius)
	lhs = x*x + dy*dy
	rhs = ((radius + 0.25)*(radius + 0.25))
	cond_t = Condition(lhs, "<", rhs)
	cond_f = Condition(lhs, ">=", rhs)
	filter_height = Reduction(([x], [dx_int]), ([x, dy], [x_radius_int, dy_int]), Int, "filter_height")
	#filter_height.defn = [ Reduce(filter_height(x), Select(cond_t, 1, 0), Op.Sum) ]
	filter_height.defn = [ Case(cond_t, Reduce(filter_height(x),1,Op.Sum)), 
			Case(cond_f, Reduce(filter_height(x), 0, Op.Sum)) ]
	filter_height.default = 0

	x_dx_int = Interval(Int, radius, R - 2*radius - 2)
	final = Reduction(([x, y, c], [rows, cols, c_int]), ([x, y, c, dx], [x_dx_int, cols, c_int, dx_int]), Int, "final")
	final.defn = [ Reduce(final(x, y, c), vert(x+ dx, y, c, clamp(filter_height(dx), 0, radius + 1)), Op.Max) ]

	return final


# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../')

from compiler import *
from constructs import *

def test_blur():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")

    cr = Interval(Int, 0, 2)

    xrow = Interval(Int, 1, R)
    xcol = Interval(Int, 1, C)

    yrow = Interval(Int, 2, R-1)
    ycol = Interval(Int, 2, C-1)

    xcond = Condition(x, '>=', 1) & Condition(x, '<=', R) & \
            Condition(y, '>=', 1) & Condition(y, '<=', C)

    ycond = Condition(x, '>=', 2) & Condition(x, '<=', R-1) & \
            Condition(y, '>=', 2) & Condition(y, '<=', C-1)

    img = Image(Float, "input", [3, R+2, C+2])

    blurx = Function(([c, x, y], [cr, xrow, xcol]), Float, "blurx")
    blurx.defn = [ Case(xcond, (img(c, x-1, y) + \
                                img(c, x  , y) + \
                                img(c, x+1, y) / 3.0)) ]

    blury = Function(([c, x, y], [cr, yrow, ycol]), Float, "blury")
    blury.defn = [ Case(ycond, (blurx(c, x, y-1) + \
                                blurx(c, x, y  ) + \
                                blurx(c, x, y+1) / 3.0)) ]

    groups = [[blurx, blury]]

    p_est = [ (R, 1024), (C, 1024) ]

    # build the pipeline
    pipeline = buildPipeline([blury],
                             grouping = groups,
                             param_estimates = p_est,
                             pipe_name = "blur")

    '''
    filename = 'blur_graph.dot'
    pipeline.originalGraph.write(filename)
    
    filename = 'blur_graph_grouped.dot'
    g = pipeline.drawPipelineGraph()
    g.write(filename)
    '''

    filename = 'blur_naive.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()



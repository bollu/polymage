# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../')

from compiler import *
from constructs import *

def test_harris_corner():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")

    row = Interval(Int, 0, R+1)
    col = Interval(Int, 0, C+1)

    cond = Condition(x, '>=', 1) & Condition(x, '<=', R) & \
           Condition(y, '<=', C) & Condition(y, '>=', 1)

    condRed = Condition(x, '>=', 2) & Condition(x, '<=', R-1) & \
              Condition(y, '<=', C-1) & Condition(y, '>=', 2)

    zeroCond = Condition(x, '==', 0) & Condition(y, '==', 0)
    endCond = Condition(x, '==', R) & Condition(y, '==', C)

    img = Image(Float, "img", [R+2, C+2])

    Iy = Function(([x, y], [row, col]), Float, "Iy")
    Iy.defn = [ Case(cond, img(x-1, y-1)*(-1.0/12) + img(x-1, y+1)*(1.0/12) + \
                           img(x, y-1)*(-2.0/12) + img(x, y+1)*(2.0/12) + \
                           img(x+1, y-1)*(-1.0/12) + img(x+1, y+1)*(1.0/12)) ]

    Ix = Function(([x, y], [row, col]), Float, "Ix")
    Ix.defn = [ Case(cond, img(x-1, y-1)*(-1.0/12) + img(x+1, y-1)*(1.0/12) + \
                           img(x-1, y)*(-2.0/12) + img(x+1, y)*(2.0/12) + \
                           img(x-1, y+1)*(-1.0/12) + img(x+1, y+1)*(1.0/12)) ]

    Ixx = Function(([x, y], [row, col]), Float, "Ixx")
    Ixx.defn = [ Case(cond, Ix(x, y) * Ix(x, y)) ]

    Iyy = Function(([x, y], [row, col]), Float, "Iyy")
    Iyy.defn = [ Case(cond, Iy(x, y) * Iy(x, y)) ]

    Ixy = Function(([x, y], [row, col]), Float, "Ixy")
    Ixy.defn = [ Case(cond, Ix(x, y) * Iy(x, y)) ]

    Sxx = Function(([x, y], [row, col]), Float, "Sxx")
    Sxx.defn = [ Case(condRed, Ixx(x-1, y-1) + Ixx(x-1, y) + Ixx(x-1, y+1) +\
                               Ixx(x, y-1)   + Ixx(x, y)   + Ixx(x, y+1) +\
                               Ixx(x+1, y-1) + Ixx(x+1, y) + Ixx(x+1, y+1)) ]

    Syy = Function(([x, y], [row, col]), Float, "Syy")
    Syy.defn = [ Case(condRed, Iyy(x-1, y-1) + Iyy(x-1, y) + Iyy(x-1, y+1) +\
                               Iyy(x, y-1)   + Iyy(x, y)   + Iyy(x, y+1) +\
                               Iyy(x+1, y-1) + Iyy(x+1, y) + Iyy(x+1, y+1)) ]

    Sxy = Function(([x, y], [row, col]), Float, "Sxy")
    Sxy.defn = [ Case(condRed, Ixy(x-1, y-1) + Ixy(x-1, y) + Ixy(x-1, y+1) +\
                               Ixy(x, y-1)   + Ixy(x, y)   + Ixy(x, y+1) +\
                               Ixy(x+1, y-1) + Ixy(x+1, y) + Ixy(x+1, y+1)) ]

    det = Function(([x, y], [row, col]), Float, "det")
    det.defn = [ Case(condRed, Sxx(x, y) * Syy(x, y) - Sxy(x, y) * Sxy(x, y)) ]

    trace = Function(([x, y], [row, col]), Float, "trace")
    trace.defn = [ Case(condRed, Sxx(x, y) + Syy(x, y)) ]

    harris = Function(([x, y], [row, col]), Float, "harris")
    harris.defn = [ Case(condRed, 
                         det(x, y) - 0.04 * trace(x, y) * trace(x, y)) ]

    # test boundary case
    harris.defn.append(Case(zeroCond, 0.0))

    # test self reference
    harris.defn.append(Case(endCond , harris(0, 0)))

    pipeline = buildPipeline([harris], grouping = [[Ix, Iy, Ixx, Iyy, Ixy], [Sxx, Sxy, Syy], [det, trace, harris]], pipe_name="harris")

    filename = 'harris_graph.dot'
    pipeline.originalGraph.write(filename)
    
    filename = 'harris_graph_grouped.dot'
    g = pipeline.drawPipelineGraph()
    g.write(filename)

    filename = 'harris_naive.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generateCode().__str__())
    c_file.close()



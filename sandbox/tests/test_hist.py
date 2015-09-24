# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../')

from compiler import *
from constructs import *

def test_histogram():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")
    c = Variable(Int, "c")

    row = Interval(Int, 0, R+1)
    col = Interval(Int, 0, C+1)
    cr = Interval(Int, 0, 2)

    cond = Condition(x, '>=', 0) & Condition(x, '<=', R-1) & \
           Condition(y, '<=', C-1) & Condition(y, '>=', 0)

    img = Image(Float, "img", [R+2, C+2])


    cond0 = Condition(c, '==', 0)
    cond1 = Condition(c, '==', 1)
    cond2 = Condition(c, '==', 2)

    # Iterates over [c, x, y] and reduces img(x, y) to hist(c, x).
    # Over each row of img(x, y), op(c) is applied and the result is stored
    # at hist(c, x). Here we use op(0) = Min, op(1) = Max, op(2) = Sum.
    hist = Reduction(([c, x], [cr, row]), \
                     ([c, x, y], [cr, row, col]), \
                     Float, "hist")
    hist.defn = [ Case(cond0, Reduce(hist(c, x),
                                     img(x, y),
                                     Op.Min)),
                  Case(cond1, Reduce(hist(c, x),
                                     img(x, y),
                                     Op.Max)),
                  Case(cond2, Reduce(hist(c, x),
                                     img(x, y),
                                     Op.Sum)) ]

    pipeline = buildPipeline([hist], \
                             pipe_name="hist")

    '''
    filename = 'hist_graph.dot'
    pipeline.originalGraph.write(filename)
    
    filename = 'hist_graph_grouped.dot'
    g = pipeline.drawPipelineGraph()
    g.write(filename)
    '''

    filename = 'hist_naive.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()



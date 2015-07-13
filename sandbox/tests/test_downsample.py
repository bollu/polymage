from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../')

from compiler import *
from constructs import *

def test_downsample():

    print()
    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")

    L = 10

    row = {}
    col = {}
    cond = {}
    RL = {}
    CL = {}
    for l in range(0, L):
        RL[l] = (R/(1<<l))
        CL[l] = (C/(1<<l))

        row[l] = Interval(Int, 0, RL[l]+1)
        col[l] = Interval(Int, 0, CL[l]+1)

        cond[l] = Condition(x, '>=', 1) & Condition(x, '<=', RL[l]) & \
                  Condition(y, '>=', 1) & Condition(y, '<=', CL[l])

    img = Image(Float, "img", [R+2, C+2])

    down = {}

    for l in range(0, L):
        if l == 0:
            up = img
        else:
            up = down[l-1]

        down[l] = Function(([x, y], [row[l], col[l]]), 
                           Float, "down_"+str(l))
        down[l].defn = [ Case(cond[l], 
                             (up(2*x-1, 2*y-1) \
                           +  up(2*x+1, 2*y-1) \
                           +  up(2*x-1, 2*y+1) \
                           +  up(2*x+1, 2*y+1)) * 0.0625 \
                           + (up(2*x-1, 2*y  ) \
                           +  up(2*x+1, 2*y  ) \
                           +  up(2*x  , 2*y-1) \
                           +  up(2*x  , 2*y+1)) * 0.125 \
                           +  up(2*x  , 2*y  )  * 0.25) ]

    # manual grouping
    group1 = []
    for l in range(0, 3):
        group1.append(down[l])

    group2 = []
    for l in range(3, L):
        group2.append(down[l])

    groups = [group1, group2]

    # build the pipeline
    pipeline = buildPipeline([down[L-1]], grouping = groups)

    filename = 'down_graph.dot'
    pipeline.originalGraph.write(filename)
 
    filename = 'down_graph_grouped.dot'
    g = pipeline.drawPipelineGraph()
    g.write(filename)

from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../')

from compiler import *
from constructs import *

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

def test_sampling():
    img = Image(Float, "img", [R+2, C+2])

    down = {}
    for l in range(1, L):
        if l == 1:
            upper = img
        else:
            upper = down[l-1]

        down[l] = Function(([x, y], [row[l], col[l]]),
                           Float, "down_"+str(l))
        down[l].defn = [ Case(cond[l],
                             (upper(2*x-1, 2*y-1) \
                           +  upper(2*x+1, 2*y-1) \
                           +  upper(2*x-1, 2*y+1) \
                           +  upper(2*x+1, 2*y+1)) * 0.0625 \
                           + (upper(2*x-1, 2*y  ) \
                           +  upper(2*x+1, 2*y  ) \
                           +  upper(2*x  , 2*y-1) \
                           +  upper(2*x  , 2*y+1)) * 0.125 \
                           +  upper(2*x  , 2*y  )  * 0.25) ]

    # not a standard upsampling -
    up = {}
    for l in range(L-1, 0, -1):
        if l == L-1:
            lower = down[l]
        else:
            lower = up[l]

        up[l-1] = Function(([x, y], [row[l-1], col[l-1]]),
                           Float, "up_"+str(l-1))
        up[l-1].defn = [ Case(cond[l-1],
                             (lower(x/2-1, y/2-1) \
                           +  lower(x/2+1, y/2-1) \
                           +  lower(x/2-1, y/2+1) \
                           +  lower(x/2+1, y/2+1)) * 0.0625 \
                           + (lower(x/2-1, y/2  ) \
                           +  lower(x/2+1, y/2  ) \
                           +  lower(x/2  , y/2-1) \
                           +  lower(x/2  , y/2+1)) * 0.125 \
                           +  lower(x/2  , y/2  )  * 0.25) ]

    # manual grouping
    group1 = []
    for l in range(1, 5):
        group1.append(down[l])

    group2 = []
    for l in range(5, L):
        group2.append(down[l])
    group2.append(up[L-2])
    group2.append(up[L-3])

    group3 = []
    for l in range(3, L-3):
        group3.append(up[l])

    group4 = []
    for l in range(0, 3):
        group4.append(up[l])

    groups = [group1, group2, group3, group4]

    # build the pipeline
    pipeline = buildPipeline([up[0]], grouping = groups)

    filename = 'down_graph.dot'
    pipeline.originalGraph.write(filename)
 
    filename = 'down_graph_grouped.dot'
    g = pipeline.drawPipelineGraph()
    g.write(filename)

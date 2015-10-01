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

    def down_sampler(img_in, l):
        img_out = Function(([x, y], [row[l], col[l]]),
                           Float, "down_"+str(l))
        img_out.defn = [ Case(cond[l],
                             (img_in(2*x-1, 2*y-1) \
                           +  img_in(2*x+1, 2*y-1) \
                           +  img_in(2*x-1, 2*y+1) \
                           +  img_in(2*x+1, 2*y+1)) * 0.0625 \
                           + (img_in(2*x-1, 2*y  ) \
                           +  img_in(2*x+1, 2*y  ) \
                           +  img_in(2*x  , 2*y-1) \
                           +  img_in(2*x  , 2*y+1)) * 0.125 \
                           +  img_in(2*x  , 2*y  )  * 0.25) ]

        return img_out

    # not a standard upsampling -
    def up_sampler(img_in, l):
        img_out = Function(([x, y], [row[l-1], col[l-1]]),
                           Float, "up_"+str(l-1))
        img_out.defn = [ Case(cond[l-1],
                             (img_in(x/2-1, y/2-1) \
                           +  img_in(x/2+1, y/2-1) \
                           +  img_in(x/2-1, y/2+1) \
                           +  img_in(x/2+1, y/2+1)) * 0.0625 \
                           + (img_in(x/2-1, y/2  ) \
                           +  img_in(x/2+1, y/2  ) \
                           +  img_in(x/2  , y/2-1) \
                           +  img_in(x/2  , y/2+1)) * 0.125 \
                           +  img_in(x/2  , y/2  )  * 0.25) ]

        return img_out

    # Pipeline
    img = Image(Float, "img", [R+2, C+2])

    down = {}
    down[1] = down_sampler(img, 1)
    for l in range(2, L):
        down[l] = down_sampler(down[l-1], l)

    up = {}
    up[L-2] = up_sampler(down[L-1], l)
    for l in range(L-2, 0, -1):
        up[l-1] = up_sampler(up[l], l)

    live_outs = [up[0]]

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

    p_est = [ (R, 1024), (C, 1024) ]

    # build the pipeline
    pipeline = buildPipeline(live_outs,
                             #grouping = groups,
                             param_estimates = p_est,
                             pipe_name="up")

    '''
    filename = 'down_graph.dot'
    pipeline.originalGraph.write(filename)
 
    filename = 'down_graph_grouped.dot'
    g = pipeline.drawPipelineGraph()
    g.write(filename)
    '''

    filename = 'sampling_naive.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()


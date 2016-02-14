from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess

sys.path.insert(0, '../')

from compiler import *
from constructs import *

def test_harris_corner():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")

    row = Interval(Int, 0, R-1)
    col = Interval(Int, 0, C-1)

    img = Image(Float, "img", [R, C])

    alpha = 0.6

    F = {}
    L = 3
    for l in range(0, (2**L)-1):
        F[l] = Function(([x, y], [row, col]), Float, "F"+str(l))

    for l in range(0, (2**(L-1))-1):
        F[l].defn = [ (alpha) * F[2*l+1](x, y) + (1-alpha) * F[2*l+2](x, y) ]
    for l in range((2**(L-1))-1, (2**L)-1):
        F[l].defn = [ l * img(x, y) ]

    p_est = [ (R, 1024), (C, 1024) ]

    # build the pipeline
    pipeline = buildPipeline([F[0]],
                             param_estimates = p_est,
                             group_size = 100,
                             pipe_name = "tree")

    filename = "tree_graph"
    dot_file = filename+".dot"
    png_file = filename+".png"
    g = pipeline.pipeline_graph
    g.write(filename+".dot")
    dotty_str = "dot -Tpng "+dot_file+" -o "+png_file
    subprocess.check_output(dotty_str, shell=True)

    filename = 'tree.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()

    return

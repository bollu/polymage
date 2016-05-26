from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess
sys.path.insert(0, '../')

from compiler import *
from constructs import *

def _test_stencil():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")

    xrow = Interval(Int, 1, R)
    xcol = Interval(Int, 1, C)

    yrow = Interval(Int, 2, R-1)
    ycol = Interval(Int, 2, C-1)

    bounds = Condition(x, '>=', 1) & Condition(x, '<=', R) & \
            Condition(y, '>=', 1) & Condition(y, '<=', C)

    img = Image(Float, "input", [R+1, C+1])

    kernel = [[1, 1, 1], [1, 0, 1], [1, 1, 1]]
    stencil = TStencil([x, y], kernel, "stencil")
    print(stencil)


    groups = [stencil]

    p_est = [ (R, 1024), (C, 1024) ]

    # build the pipeline
    pipeline = buildPipeline([stencil],
                             grouping = groups,
                             param_estimates = p_est,
                             pipe_name = "blur")

    filename = "blur_graph"
    dot_file = filename+".dot"
    png_file = filename+".png"
    g = pipeline.pipeline_graph
    g.write(filename+".dot")
    dotty_str = "dot -Tpng "+dot_file+" -o "+png_file
    subprocess.check_output(dotty_str, shell=True)

    filename = 'blur_naive.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()



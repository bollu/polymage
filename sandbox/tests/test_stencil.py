from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess
sys.path.insert(0, '../')

from compiler import *
from constructs import *


def test_stencil():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")

    xrow = Interval(Int, 1, R)
    xcol = Interval(Int, 1, C)

    bounds = Condition(x, '>=', 2) & Condition(x, '<=', R - 2) & \
        Condition(y, '>=', 2) & Condition(y, '<=', C - 2)

    img = Image(Float, "input", [R + 1, C + 1])

    kernel = [[1, 0, 0], [0, 2, 0], [0, 0, 1]]
    stencil = Stencil(img, [x, y], kernel)
    print(stencil)

    f = Function(([x, y], [xrow, xcol]), Float, "stencilfn")
    f.defn = [Case(bounds, stencil)]

    groups = [f]

    p_est = [(R, 1024), (C, 1024)]

    # build the pipeline
    pipeline = buildPipeline([f],
                             grouping=groups,
                             param_estimates=p_est,
                             pipe_name="blur")

    filename = "blur_graph"
    dot_file = filename + ".dot"
    png_file = filename + ".png"
    g = pipeline.pipeline_graph
    g.write(filename+".dot")
    dotty_str = "dot -Tpng "+dot_file+" -o "+png_file
    subprocess.check_output(dotty_str, shell=True)

    filename = 'test_stencil.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()


if __name__ == "__main__":
    print("running test...")
    test_stencil()

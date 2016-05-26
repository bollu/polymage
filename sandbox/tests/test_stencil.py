from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess
sys.path.insert(0, '../')

from compiler import *
from constructs import *

D = Parameter(Int, "D")
R = Parameter(Int, "R")
C = Parameter(Int, "C")
z = Variable(Int, "z")
y = Variable(Int, "y")
x = Variable(Int, "x")

def test_stencil1():

    xrow = Interval(Int, 1, R)
    xcol = Interval(Int, 1, C)

    bounds = Condition(x, '>=', 2) & Condition(x, '<=', R - 2) & \
             Condition(y, '>=', 2) & Condition(y, '<=', C - 2)

    img = Image(Float, "input", [R + 2, C + 2])

    kernel1 = [[1, 0, 0], [0, 2, 0], [0, 0, 1]]
    stencil_f1 = Stencil(img, [x, y], kernel1)
    f1 = Function(([x, y], [xrow, xcol]), Float, "f1")
    f1.defn = [Case(bounds, stencil_f1)]

    kernel2 = [[-1, 2, -3.8], [0, -1, 0]]
    stencil_f2 = Stencil(f1, [x, y], kernel2)
    f2 = Function(([x, y], [xrow, xcol]), Float, "f2")
    f2.defn = [Case(bounds, stencil_f2)]

    assert str(stencil_f1.input_func.name) == 'input'
    assert stencil_f1.iter_vars == [x, y]
    assert stencil_f1.sizes == [3, 3]
    assert stencil_f1.origin == [1, 1]
    assert stencil_f1.kernel == kernel1
    assert str(Stencil.macro_expand(stencil_f1)) == \
        '((input((x + -1), (y + -1)) + ' + \
        '(input(x, y) * 2)) + ' + \
        'input((x + 1), (y + 1)))'

    assert str(stencil_f2.input_func.name) == 'f1'
    assert stencil_f2.iter_vars == [x, y]
    assert stencil_f2.sizes == [2, 3]
    assert stencil_f2.origin == [0, 1]
    assert stencil_f2.kernel == kernel2
    assert str(Stencil.macro_expand(stencil_f2)) == \
        '((((f1(x, (y + -1)) * -1) + ' + \
        '(f1(x, y) * 2)) + ' + \
        '(f1(x, (y + 1)) * -3.8)) + ' + \
        '(f1((x + 1), y) * -1))'

    kernel3 = [[[ 1,  2, -1], [-2,  1, -2], [-1,  2,  1]],
               [[ 2,  1,  2], [ 1, -2,  1], [ 2,  1,  2]],
               [[-1,  2,  1], [-2,  1, -2], [ 1,  2, -1]]]

    groups = [f1, f2]

    p_est = [(R, 1024), (C, 1024)]

    # build the pipeline
    pipeline = buildPipeline([f2],
                             grouping=groups,
                             param_estimates=p_est,
                             pipe_name="stencil")

    filename = "stencil_graph"
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

    return

def test_stencil2():

    grid = Image(Float, "grid", [D + 2, R + 2, C + 2])

    kernel1 = [[[ 1,  2, -1], [-2,  1, -2], [-1,  2,  1]],
               [[ 2,  1,  2], [ 1, -2,  1], [ 2,  1,  2]],
               [[-1,  2,  1], [-2,  1, -2], [ 1,  2, -1]]]
    stencil1 = Stencil(grid, [x, y, z], kernel1)

    assert stencil1.iter_vars == [x, y, z]
    assert stencil1.sizes == [3, 3, 3]
    assert stencil1.origin == [1, 1, 1]
    assert stencil1.kernel == kernel1
    assert str(Stencil.macro_expand(stencil1)) == \
        '(((((((((((((((((((((((((' + \
        '(grid((x + -1), (y + -1), (z + -1)) + ' + \
        '(grid((x + -1), (y + -1), z) * 2)) + ' + \
        '(grid((x + -1), (y + -1), (z + 1)) * -1)) + ' + \
        '(grid((x + -1), y, (z + -1)) * -2)) + ' + \
        'grid((x + -1), y, z)) + ' + \
        '(grid((x + -1), y, (z + 1)) * -2)) + ' + \
        '(grid((x + -1), (y + 1), (z + -1)) * -1)) + ' + \
        '(grid((x + -1), (y + 1), z) * 2)) + ' + \
        'grid((x + -1), (y + 1), (z + 1))) + ' + \
        '(grid(x, (y + -1), (z + -1)) * 2)) + ' + \
        'grid(x, (y + -1), z)) + ' + \
        '(grid(x, (y + -1), (z + 1)) * 2)) + ' + \
        'grid(x, y, (z + -1))) + ' + \
        '(grid(x, y, z) * -2)) + ' + \
        'grid(x, y, (z + 1))) + ' + \
        '(grid(x, (y + 1), (z + -1)) * 2)) + ' + \
        'grid(x, (y + 1), z)) + ' + \
        '(grid(x, (y + 1), (z + 1)) * 2)) + ' + \
        '(grid((x + 1), (y + -1), (z + -1)) * -1)) + ' + \
        '(grid((x + 1), (y + -1), z) * 2)) + ' + \
        'grid((x + 1), (y + -1), (z + 1))) + ' + \
        '(grid((x + 1), y, (z + -1)) * -2)) + ' + \
        'grid((x + 1), y, z)) + ' + \
        '(grid((x + 1), y, (z + 1)) * -2)) + ' + \
        'grid((x + 1), (y + 1), (z + -1))) + ' + \
        '(grid((x + 1), (y + 1), z) * 2)) + ' + \
        '(grid((x + 1), (y + 1), (z + 1)) * -1))'


    kernel2 = [[[ 1,  2], [-2,  1], [-1,  2]],
               [[ 2,  1], [ 1, -2], [ 2,  1]]]
    stencil2 = Stencil(grid, [x, y, z], kernel2)

    assert stencil2.iter_vars == [x, y, z]
    assert stencil2.sizes == [2, 3, 2]
    assert stencil2.origin == [0, 1, 0]
    assert stencil2.kernel == kernel2
    assert str(Stencil.macro_expand(stencil2)) == \
        '((((((((((' + \
        '(grid(x, (y + -1), z) + ' + \
        '(grid(x, (y + -1), (z + 1)) * 2)) + ' + \
        '(grid(x, y, z) * -2)) + ' + \
        'grid(x, y, (z + 1))) + ' + \
        '(grid(x, (y + 1), z) * -1)) + ' + \
        '(grid(x, (y + 1), (z + 1)) * 2)) + ' + \
        '(grid((x + 1), (y + -1), z) * 2)) + ' + \
        'grid((x + 1), (y + -1), (z + 1))) + ' + \
        'grid((x + 1), y, z)) + ' + \
        '(grid((x + 1), y, (z + 1)) * -2)) + ' + \
        '(grid((x + 1), (y + 1), z) * 2)) + ' + \
        'grid((x + 1), (y + 1), (z + 1)))'

    return

def _test_stencil3():
    kernel = [[[ 1,  2], [-2,  1, -2], [-1,  2,  1]],
              [[ 2,  1], [ 1, -2,  1], [ 2,  1,  2]]]
    assert is_valid_kernel(kernel, 3) == False

    return

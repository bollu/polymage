from __init__ import *
import sys

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

def w_jacobi(U, F, name, app_data):
    pipe_data = app_data['pipe_data']

    z = pipe_data['z']
    y = pipe_data['y']
    x = pipe_data['x']

    interior = pipe_data['interior']
    inner_box = interior['inner_box']
    invhh = pipe_data['invhh']
    c = pipe_data['jacobi_c']

    extent = pipe_data['extent']

    W = Function(([z, y, x], [extent, extent, extent]),
                  Double, str(name))

    W.defn = [ Case(inner_box,
                    (U(z, y, x) \
                   - c * (( \
                     U(z  , y  , x  ) * 6.0 \
                   - U(z-1, y  , x  )       \
                   - U(z+1, y  , x  )       \
                   - U(z  , y-1, x  )       \
                   - U(z  , y+1, x  )       \
                   - U(z  , y  , x-1)       \
                   - U(z  , y  , x+1)       \
                   ) * invhh                \
                   - F(z  , y  , x  ) ))) ]

    return W


def stencil_jacobi(app_data):

    pipe_data = app_data['pipe_data']

    T = app_data['T']
    N = pipe_data['N']

    V = Image(Double, "V_", [N+2, N+2, N+2])
    F = Image(Double, "F_", [N+2, N+2, N+2])
    jacobi = {}
    jacobi[0] = V

    fname = 'jacobi'
    for t in range(1, T):
        jacobi[t] = w_jacobi(jacobi[t-1], F, fname+"_"+str(t), app_data)
    jacobi[T] = w_jacobi(jacobi[T-1], F, fname, app_data)

    return jacobi[T]

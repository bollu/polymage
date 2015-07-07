# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../../frontend')
from Constructs import *
from Expression import *

# Have to change the string comparision tests to something
# better. Maybe eval the rexpression at the point.
def test_param():
    N = Parameter(Float, "N")
    assert N.has(Variable) == False
    assert N.has(Parameter) == True
    assert N.name == "N"
    assert N.typ == Float

def test_variable():
    x = Variable(UInt, "x")
    assert x in x.collect(Variable)
    assert x not in x.collect(Parameter)
    assert x.name == "x"
    assert x.typ == UInt

def test_types():
    assert Value(3.0, Float) == 3.0 
    assert Value(3.0, Double) == 3.0
    assert Value(3, Int) == 3 
    assert Value(3, UInt) == 3 
    assert Value(3, Char) == 3 
    assert Value(3, UChar) == 3 
    assert Value(3, Short) == 3 
    assert Value(3, UShort) == 3
    assert Value(3, Long) == 3 
    assert Value(3, ULong) == 3

def test_interval():
    I = Interval(UInt, Value(3, UInt), Value(5, UInt))
    assert I.typ == UInt
    assert I.lowerBound == 3
    assert I.upperBound == 5

def test_image():
    N = Parameter(UInt, "N")
    Img = Image(UChar, "Input", [N, N, 3])
    assert Img.typ  == UChar
    assert Img.name == "Input"
    assert Img.dimensions[0] == N
    assert Img.dimensions[1] == N
    assert Img.dimensions[2] == 3

def test_expr():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    z = Variable(UInt, "z")
    w = x * (y * z) - N
    assert N in w.collect(Parameter)
    assert w.__str__().replace(' ', '') == "((x*(y*z))-N)"
    w = (+(((-(((x - (y + z)) / z) % 2) ^ 3) >> z) << 5) | 1) & y
    for var in [x, y, z]:
        assert var in w.collect(Variable)
    assert w.__str__().replace(' ', '') == "((+((((-((((x-(y+z))/z)%2))^3)>>z)<<5))|1)&y)"
    i = 3
    w = x * i
    assert w.__str__().replace(' ', '') == "(x*3)"

def test_condition():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    c1 = Condition(x+N, '>', y)
    assert x in c1.collect(Variable)
    assert y in c1.collect(Variable)
    assert N in c1.collect(Parameter)
    assert c1.__str__().replace(' ', '') == "((x+N)>y)"
    c2 = Condition(x , '<', 5)
    assert c2.__str__().replace(' ', '') == "(x<5)"
    c3 = c1 & c2
    assert c3.__str__().replace(' ', '') == "(((x+N)>y)&&(x<5))"

def test_case():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    c1 = Case(Condition(x+1, '>', y), x-y)
    assert c1.condition.__str__().replace(' ', '')  == "((x+1)>y)"
    assert c1.expression.__str__().replace(' ', '') == "(x-y)"

def test_function():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    r = Interval(UInt, 0, N-1)
    func1 = Function(([x, y], [r, r]), UInt, "add")
    func1.defn = [ x + y ]
    assert r in  func1.getObjects(Interval)
    assert func1.defn[0].__str__().replace(' ','') == "(x+y)" 
    func2 = Function(([x, y], [r, r]), UInt, "max")
    func2.defn = [ Case(Condition(x, '>', y), x),
                   Case(Condition(x, '<=', y), y) ]
    assert func2.defn[0].__str__().replace(' ','') == "Case((x>y)){x}" 
    assert func2.defn[1].__str__().replace(' ','') == "Case((x<=y)){y}" 

def test_conjuncts():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")
    c1 = Condition(x, '<', y)
    c2 = Condition(x, '>', N)
    c3 = Condition(x, '<', 3)
    c4 = Condition(x + y, '<=',N)
    c5 = Condition(x, '!=', 2*N)
    cond  = (c1 & c2) | (c3 | c4) & c5
    conjuncts = cond.splitToConjuncts()
    conjuncts_str = []
    for conjunct in conjuncts:
        conjunct_str = ""
        for cond in conjunct:
            conjunct_str = conjunct_str + cond.__str__()
        conjuncts_str.append(conjunct_str)

    assert len(conjuncts_str) == 5
    assert conjuncts_str[0].replace(' ','') == "(x<y)(x>N)" 
    assert conjuncts_str[1].replace(' ','') == "(x<3)(x<(2*N))" 
    assert conjuncts_str[2].replace(' ','') == "(x<3)(x>(2*N))" 
    assert conjuncts_str[3].replace(' ','') == "((x+y)<=N)(x<(2*N))" 
    assert conjuncts_str[4].replace(' ','') == "((x+y)<=N)(x>(2*N))" 

def test_divide():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")

    expr = (2*x + N//2 + 3*y)*3//2
    assert isAffine(expr) == True
    coeff = getAffineVarAndParamCoeff(expr)
    assert coeff[x] == Fraction(3, 1) 
    assert coeff[N] == Fraction(3, 4) 
    assert coeff[y] == Fraction(9, 2) 
    expr = (N-1)//2
    assert getConstantFromExpr(expr, affine = True) == Fraction(-1, 2)

def test_overload():
    N = Parameter(UInt, "N")
    x = Variable(UInt, "x")
    y = Variable(UInt, "y")

    expr = 1//N
    assert expr.__str__().replace(' ', '') == "(1/N)"
    expr = N//1
    assert expr.__str__().replace(' ', '') == "N"
    expr = 1*N
    assert expr.__str__().replace(' ', '') == "N"
    expr = N*1
    assert expr.__str__().replace(' ', '') == "N"
    expr = N + 0
    assert expr.__str__().replace(' ', '') == "N"
    expr = N + 0
    assert expr.__str__().replace(' ', '') == "N"
    expr = N >> 0
    assert expr.__str__().replace(' ', '') == "N"
    expr = N << 0
    assert expr.__str__().replace(' ', '') == "N"
    expr = N >> 0
    assert expr.__str__().replace(' ', '') == "N"
    expr = N * 0
    assert expr.__str__().replace(' ', '') == "0"
    expr = 0 * N
    assert expr.__str__().replace(' ', '') == "0"

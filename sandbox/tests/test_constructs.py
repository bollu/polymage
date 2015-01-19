# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
sys.path.insert(0, '../frontend')
import constructs as c

# Have to change the string comparision tests to something
# better. Maybe eval the rexpression at the point.
def test_param():
    N = c.Parameter(c.Float, "N")
    assert N.has(c.Variable) == False
    assert N.has(c.Parameter) == True
    #assert N.isAffine() == True
    assert N.name == "N"
    assert N.typ == c.Float

def test_variable():
    x = c.Variable(c.UInt, "x")
    assert x in x.collect(c.Variable)
    assert x not in x.collect(c.Parameter)
    #assert x.isAffine() == True
    assert x.name == "x"
    assert x.typ == c.UInt

def test_types():
    assert c.Value(3.0, c.Float) == 3.0 
    assert c.Value(3.0, c.Double) == 3.0
    assert c.Value(3, c.Int) == 3 
    assert c.Value(3, c.UInt) == 3 
    assert c.Value(3, c.Char) == 3 
    assert c.Value(3, c.UChar) == 3 
    assert c.Value(3, c.Short) == 3 
    assert c.Value(3, c.UShort) == 3
    assert c.Value(3, c.Long) == 3 
    assert c.Value(3, c.ULong) == 3

def test_interval():
    I = c.Interval(c.UInt, c.Value(3, c.UInt), c.Value(5, c.UInt), c.Value(1, c.UInt))
    assert I.typ == c.UInt
    assert I.lowerBound == 3
    assert I.upperBound == 5
    assert I.step == 1

def test_image():
    N = c.Parameter(c.UInt, "N")
    Img = c.Image(c.UChar, "Input", [N, N, 3])
    assert Img.typ  == c.UChar
    assert Img.name == "Input"
    assert Img.dimensions[0] == N
    assert Img.dimensions[1] == N
    assert Img.dimensions[2] == 3

def test_expr():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    z = c.Variable(c.UInt, "z")
    w = x * (y * z) - N
    assert N in w.collect(c.Parameter)
    assert w.__str__().replace(' ', '') == "((x*(y*z))-N)"
    w = (+(((-(((x - (y + z)) / z) % 2) ^ 3) >> z) << 5) | 1) & y
    for var in [x, y, z]:
        assert var in w.collect(c.Variable)
    assert w.__str__().replace(' ', '') == "((+((((-((((x-(y+z))/z)%2))^3)>>z)<<5))|1)&y)"
    i = 3
    w = x * i
    assert w.__str__().replace(' ', '') == "(x*3)"

def test_condition():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    c1 = c.Condition(x+N, '>', y)
    assert x in c1.collect(c.Variable)
    assert y in c1.collect(c.Variable)
    assert N in c1.collect(c.Parameter)
    assert c1.__str__().replace(' ', '') == "((x+N)>y)"
    c2 = c.Condition(x , '<', 5)
    assert c2.__str__().replace(' ', '') == "(x<5)"
    c3 = c1 & c2
    assert c3.__str__().replace(' ', '') == "(((x+N)>y)&&(x<5))"

def test_case():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    c1 = c.Case(c.Condition(x+1, '>', y), x-y)
    assert c1.condition.__str__().replace(' ', '')  == "((x+1)>y)"
    assert c1.expression.__str__().replace(' ', '') == "(x-y)"

def test_function():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    func1 = c.Function(c.UInt, "add")
    r = c.Interval(c.UInt, 0, N-1, 1)
    func1.variableDomain = ([x, y], [r, r])
    func1.definition = x + y
    assert r in  func1.getObjects(c.Interval)
    assert func1.definition[0].__str__().replace(' ','') == "(x+y)" 
    func2 = c.Function(c.UInt, "max")
    func2.variableDomain = ([x, y], [r, r])
    func2.definition = c.Case(c.Condition(x, '>', y), x)
    func2.definition = c.Case(c.Condition(x, '<=', y), y)
    assert func2.definition[0].__str__().replace(' ','') == "Case((x>y)){x}" 
    assert func2.definition[1].__str__().replace(' ','') == "Case((x<=y)){y}" 

def test_affine():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    assert(c.isAffine(x + y) == True)
    assert(c.isAffine(3) == True)
    assert(c.isAffine(x*y) == False)
    assert(c.isAffine(-x + N + 3*y) == True)
    assert(c.isAffine(2*x + N/2 + 3*y) == True)
    c1 = c.Condition(x, '<', 2*y)
    c2 = c.Condition(x, '>', 2-y)
    c3 = c.Condition(x, '>=', x*y)
    c4 = c.Condition(x + 2*N, '<=', y + N)
    c5 = c.Condition(x*N, '!=', y)
    assert(c.isAffine(c1) == True)
    assert(c.isAffine(c2) == True)
    assert(c.isAffine(c3) == False)
    assert(c.isAffine(c4) == True)
    assert(c.isAffine(c5) == False)

def test_coeff():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    coeff = c.getAffineVarAndParamCoeff(1+x)
    assert(coeff[x] == 1)
    coeff = c.getAffineVarAndParamCoeff(1+x +y)
    assert(coeff[x] == 1 and coeff[y] == 1)
    coeff = c.getAffineVarAndParamCoeff(3)
    assert(coeff == {})
    coeff = c.getAffineVarAndParamCoeff(N*x + y)
    assert(coeff == {})
    coeff = c.getAffineVarAndParamCoeff(x*y)
    assert(coeff == {})
    coeff = c.getAffineVarAndParamCoeff(2*(x*3+y +N +x + y -5) 
                                      + 3*(-x) + 4*(-y) + N)
    assert(coeff[x] == 5 and coeff[y] == 0 and coeff[N] == 3)

def test_conjuncts():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")
    c1 = c.Condition(x, '<', y)
    c2 = c.Condition(x, '>', N)
    c3 = c.Condition(x, '<', 3)
    c4 = c.Condition(x + y, '<=',N)
    c5 = c.Condition(x, '!=', 2*N)
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
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")

    expr = (2*x + N//2 + 3*y)*3//2
    assert c.isAffine(expr) == True
    coeff = c.getAffineVarAndParamCoeff(expr)
    assert coeff[x] == Fraction(3, 1) 
    assert coeff[N] == Fraction(3, 4) 
    assert coeff[y] == Fraction(9, 2) 
    expr = (N-1)//2
    assert c.getConstantFromExpr(expr, affine = True) == Fraction(-1, 2)

def test_simplify():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")

    expr = (2*x + 3*x + 3 + 4 + 1 + 5*y - 4*y + N + 2 * 3 + 4//2)//2 
    expr = c.simplifyExpr(expr)
    coeff = c.getAffineVarAndParamCoeff(expr)
    print(coeff)
    print(expr)
    assert coeff[x] == Fraction(5, 2)
    assert coeff[y] == Fraction(1, 2)
    assert coeff[N] == Fraction(1, 2)
    assert c.getConstantFromExpr(expr, affine = True) == 8

def test_overload():
    N = c.Parameter(c.UInt, "N")
    x = c.Variable(c.UInt, "x")
    y = c.Variable(c.UInt, "y")

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

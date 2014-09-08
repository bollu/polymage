import islpy as isl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
plt.switch_backend("QT4Agg")


# Context
# Spaces
# Creating New Sets and Relations
# Properties
# Operations
# Lexicographic Optimization
# Vertex Enumeration
# Cardinality

# Dependence Analysis
# Scheduling 
# AST Generation

# Resources 
# http://documen.tician.de/islpy/reference.html
# http://isl.gforge.inria.fr/user.html

def Spaces(ctx):
    # Whenever a new set, relation or similiar object is created from scratch, 
    # the space in which it lives needs to be specified using an Space
    
    # It can be allocated by using
    # alloc(ctx, nparam, n_in, n_out)
    # Parameters: ctx - Context, param - unsigned, n_in - unsigned, 
    # n_out - unsigned

    # or it can be created from names
    # create_from_names(ctx, set=None, in_=None, out=None, params=[])
    # Parameters: 
    # set - names of set-type variables, in - names of in-type variables,
    # out - names of out-type variables, params - names of parameter-type variables.

    n1 = ['x', 'y']
    n2 = ['p', 'q']
    p = ['N', 'M']
    # Set
    space = isl.Space.create_from_names(ctx, set = n1, params = p)
    print space
    # Relation
    space =  isl.Space.create_from_names(ctx, in_ = n1, out = n2, params = p)
    print space

def SetsAndRelations(ctx):
    n1 = ['x', 'y']
    p = ['N', 'M']
    # Set
    space = isl.Space.create_from_names(ctx, set = n1, params = p)
    s = isl.BasicSet.empty(space)
    print s
    s = isl.BasicSet.universe(space)
    print s

    # Constraints
    # equality_alloc(), inequality_aloc(), eq_from_names(), ineq_from_names()

    c1 = isl.Constraint.ineq_from_names(space, {'N':-1, 'x':1})
    c2 = isl.Constraint.ineq_from_names(space, {1:0, 'x':-1})
    c3 = isl.Constraint.ineq_from_names(space, {'M':-1, 'y':1})
    c4 = isl.Constraint.ineq_from_names(space, {1:0, 'y':-1})
    print c1, c2, c3, c4

    s = s.add_constraints([c1, c2, c3, c4])
    print s

    # Relation
    n2 = ['p', 'q']    
    space = isl.Space.create_from_names(ctx, in_ = n1, out = n2, params = p)
    r = isl.BasicMap.empty(space)
    print r
    r = isl.BasicMap.universe(space)
    print r
    r = isl.BasicMap.identity(space)
    print r
    print r.domain()
    print r.range()

    r = isl.BasicMap.universe(space)
    c1 = isl.Constraint.eq_from_names(space, {1:0, 'p':1})
    c2 = isl.Constraint.eq_from_names(space, {'x':1, 'y':1, 'q':-1})
    c3 = isl.Constraint.ineq_from_names(space, {'N':-1, 'x':1})
    c4 = isl.Constraint.ineq_from_names(space, {1:0, 'x':-1})
    c5 = isl.Constraint.ineq_from_names(space, {'M':-1, 'y':1})
    c6 = isl.Constraint.ineq_from_names(space, {1:0, 'y':-1})

    r = r.add_constraints([c1, c2, c3, c4, c5, c6])
    print r
   
    c7 = isl.Constraint.ineq_from_names(space, {1:-1, 'N':1})
    r = isl.BasicMap.identity(space)
    r = r.add_constraints([c1, c2, c3, c4, c5, c6, c7])
    # What happens now?
    print r


def SetsAsIterationSpaces(ctx):
    pass

def LocalSpaces(ctx):
    # A local space is essentially a space with zero or more existentially 
    # quantified variables.
    pass 

ctx = isl.Context()
#Spaces(ctx)
#SetsAndRelations(ctx)
#LocalSpaces(ctx)
SetsAsIterationSpaces(ctx)

def islApiExamples():
    ctx = isl.Context()

    ## Spaces
    names = ['s\'', 'x1', 'x2']
    space1 = isl.Space.create_from_names(ctx, set = names)
    names = ['s', 'y1', 'y2']
    space2 = isl.Space.create_from_names(ctx, set = names)
    #Spaces are equal when their dimensions are equal
    print space1.is_equal(space2)
    #isl dim type can be set/in/out/param/all   
    print space1.find_dim_by_name(isl._isl.dim_type.all,'x2')
    #print space1.get_id_dict(isl._isl.dim_type.set)
    newid  = isl.Id(context = ctx, name = 'x0')
    space1 = space1.set_dim_id(isl._isl.dim_type.set, 1, newid)
    print space1.get_dim_id(isl._isl.dim_type.set, 1)
    #Looks like get_id_dict does not work as expected
    #print space1.get_id_dict(isl._isl.dim_type.set)
    print space1.get_var_dict()

    ## Sets
    space = isl.Space.create_from_names(ctx, set=["i", "j"])

    bset1 = (isl.BasicSet.universe(space)
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: -1, "i": 1}))
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: 5, "i": -1}))
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: -1, "j": 1}))
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: 5, "j": -1})))
    print bset1

    bset2 = isl.BasicSet("[N]->{[x, y] : x >= 0 and x < 5 and y >= 0 and y < N+4 and N >= 0 and N < 10}")    
    print bset2
    print bset2.polyhedral_hull()
    print bset2.lexmax()
    print bset2.lexmin()
    print bset2.lexmin().is_singleton()   

    ## Points
    points = []
    bset1.foreach_point(points.append)
    print points
    point = points[0].get_coordinate_val(isl._isl.dim_type.all, 0)
    pointx = point.get_num_si()/point.get_num_si()
    # The same thing can be achieved with to_python() provided by islpy
    print point.to_python()
    point = points[0].get_coordinate_val(isl._isl.dim_type.all, 1)
    pointy = point.get_num_si()/point.get_num_si()
    print (pointx, pointy)

    ## Dependence computation

    mapspace = isl.Space.create_from_names(ctx, in_=["i", "j"], out = ["i1", "j1"], params = ["N"])
    ## Schedule and AST
    
    # Creating an identity schedule
    bmap = (isl.BasicMap.identity(mapspace)
           .add_constraint(isl.Constraint.ineq_from_names(mapspace, {1: -1, "i": 1}))
           .add_constraint(isl.Constraint.ineq_from_names(mapspace, {"N": 1, "i": -1}))
           .add_constraint(isl.Constraint.ineq_from_names(mapspace, {1: -1, "j": 1}))
           .add_constraint(isl.Constraint.ineq_from_names(mapspace, {"N": 1, "j": -1})))
    #bmap = bmap.insert_dims(isl._isl.dim_type.out, 0, 1)
    #name = bmap.get_dim_name(isl._isl.dim_type.out, 1)
    #bmap = bmap.set_dim_name(isl._isl.dim_type.out, dim, 'S_' + name)

    print bmap
    astbld = isl.AstBuild.from_context(isl.BasicSet("[N] -> { : }"))
    astbld = astbld.set_options(isl.UnionMap("{}"))
    # Printing is strange
    printer = isl.Printer.to_str(ctx)
    printer = printer.set_output_format(isl.format.C)
    printer = (isl.AstBuild.ast_from_schedule(astbld, isl.UnionMap.from_map(bmap))).print_(printer, isl.AstPrintOptions.alloc(ctx))
    print printer.get_str()
    
    extSet = isl.BasicSet("[n] -> { [i] : exists (a = [i/10] : 0 <= i and i <= n and i - 10 a <= 6) }")
    print extSet
    print extSet.get_local_space()
    print extSet.get_local_space().get_div(0)
    extMap = isl.BasicMap.from_domain_and_range(extSet, extSet)
    
    astbld = isl.AstBuild.from_context(isl.BasicSet("[n] -> { : }"))
    astbld = astbld.set_options(isl.UnionMap("{}"))
    # Printing is strange
    printer = isl.Printer.to_str(ctx)
    printer = printer.set_output_format(isl.format.C)
    printer = (isl.AstBuild.ast_from_schedule(astbld, isl.UnionMap.from_map(extMap))).print_(printer, isl.AstPrintOptions.alloc(ctx))
    print printer.get_str()

    gen = isl.BasicMap("[C, R] -> { Ixx[x, y] -> [Dir_i0', T_i0', Dir_i1', T_i1', t, i0', i1'] : t = 0 and i0' = x and i1' = y and x >= 1 and x <= R and y >= 1 and y <= C and R >= 1 and C >= 1 and 64T_i0' <= x and 64T_i0' >= -63 + x and 64T_i0' <= -32 + x - 64Dir_i0' and 64T_i0' >= -95 + x - 64Dir_i0' and Dir_i0' >= -1 and Dir_i0' <= 0 and 64T_i1' <= y and 64T_i1' >= -63 + y and 64T_i1' <= -32 + y - 64Dir_i1' and 64T_i1' >= -95 + y - 64Dir_i1' and Dir_i1' >= -1 and Dir_i1' <= 0 }")
    #gen = isl.UnionMap("[R, C] -> { harris[x, y] -> [1, x, y] : C = 10 and R = 10 and x >= 2 and x <= 9 and y >= 2 and y <= 9; Iyy[x, y] -> [0, x, y] : C = 10 and R = 10 and x >= 1 and x <= 10 and y >= 1 and y <= 10}")
    genbld = isl.AstBuild.from_context(isl.BasicSet("[C, R]->{: R > 1 and C > 1}"))
    #genbld = astbld.set_options(isl.UnionMap("{[i,j] -> unroll[0] : i < 4 or i > 99996}"))
    id_ = isl.Id.alloc(ctx, "Test1", "FakeObj")
    gen = gen.set_tuple_id(isl._isl.dim_type.in_, id_)
    
    #id_ = isl.Id.alloc(ctx, "Test2", "FakeObj")
    #print gen.get_tuple_name(isl._isl.dim_type.out)

    genbld = genbld.set_options(isl.UnionMap("[C, R] -> { [Dir_i0', T_i0', Dir_i1', T_i1', t, i0', i1'] -> unroll[x] : x = 0 or x = 2}"))
    idl =  isl.IdList.alloc(ctx, 3)
    print idl
    idl = idl.add(isl.Id.alloc(ctx, 'Dir_i0', None))
    idl = idl.add(isl.Id.alloc(ctx, 'T_i0', None))
    idl = idl.add(isl.Id.alloc(ctx, 'Dir_i1', None))
    idl = idl.add(isl.Id.alloc(ctx, 'T_i1', None))
    idl = idl.add(isl.Id.alloc(ctx, 't', None))
    idl = idl.add(isl.Id.alloc(ctx, 'i0', None))
    idl = idl.add(isl.Id.alloc(ctx, 'i1', None))
    genbld = genbld.set_iterators(idl)
    printer = isl.Printer.to_str(ctx)
    printer = printer.set_output_format(isl.format.C)
    #astRoot = isl.AstBuild.ast_from_schedule(genbld, isl.UnionMap.from_map(gen))
    print genbld.get_schedule_space()
    astRoot = genbld.ast_from_schedule(isl.UnionMap.from_map(gen))
    print genbld.get_schedule_space()
    printer = astRoot.print_(printer, isl.AstPrintOptions.alloc(ctx))
    print printer.get_str()

    applyBase = isl.BasicMap("[rows, cols] -> { Dx_2_img2[c, x, y] -> [3, c, x, y] : c >= 0 and c <= 2 and x >= 2 and 4x <= 4 + rows and y >= 2 and 2y <= 2 + cols and rows >= 1 and cols >= 1 }")
    align = isl.BasicMap("[rows, cols]->{[a, b, c, d] -> [a, d, b, c]}")
    print applyBase.apply_range(align)

    mem = isl.BasicMap("[cols, rows] -> { Uy_2_img2[c, x, y] -> [c, T_i1', T_i2', 7, 4x, 4y] : c >= 0 and c <= 2 and x >= 2 and 4x <= 4 + rows and y >= 2 and 4y <= 4 + cols and rows >= 1 and cols >= 1 and 64T_i1' <= -2 + x and 64T_i1' >= -70 + x and 64T_i2' <= -2 + y and 64T_i2' >= -73 + y }")
    acc = isl.BasicMap("[cols, rows] -> { Uy_2_img2[c, x, y] -> [c, T_i1', T_i2', 7, x-64T_i1', y-64T_i2'] : c >= 0 and c <= 2 and x >= 2 and 4x <= 4 + rows and y >= 2 and 4y <= 4 + cols and rows >= 1 and cols >= 1 and 64T_i1' <= -2 + x and 64T_i1' >= -70 + x and 64T_i2' <= -2 + y and 64T_i2' >= -73 + y and rows >= 256 and cols >= 256}")
    print acc.range().dim_min(4), acc.range().dim_max(4)

    split = isl.BasicMap("[cols, rows] -> { Uy_0_img1[c, x, y, _Mul_x, _Mul_y] -> [_T_i1, _T_i2, 1, 5, c, x, y] : exists (e0 = floor((-1 + y)/2): rows = 768 and 64_T_i1 = x - _Mul_x and 64_T_i2 = y - _Mul_y and cols = 1024 and 2e0 = -1 + y and c >= 0 and c <= 2 and x >= 2 and x <= 769 and y >= 2 and y <= 1025) }")
    print split

def matplotlibTest():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.set_autoscale_on(True)
    ax.locator_params(integer=True)
    ax.grid(True)
    
    ctx = isl.Context()
    space = isl.Space.create_from_names(ctx, set=["x", "y"])
    bset = (isl.BasicSet.universe(space)
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: -1, "x": 1}))
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: 5, "x": -1, "y":-1}))
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: -1, "y": 1}))
            .add_constraint(isl.Constraint.ineq_from_names(space, {1: 5, "y": -1})))
    drawBasicSet(bset, { 0:'x', 1:'y'}, ax)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_autoscale_on(True)
    ax.locator_params(integer=True)
    ax.grid(True)

    drawBasicSet(bset, { 0:'x', 1:'y', 'color':(0, ('red', 'blue'))}, ax)
    plt.show()

def showTiling2d():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_autoscale_on(True)
    ax.locator_params(integer=True)
    ax.set_aspect('equal')
    ax.grid(True)

    for comp in fusedstage.polyRep.polyParts:
        for polyComp in fusedstage.polyRep.polyParts[comp]:
            drawBasicSet(polyComp.sched.range(), {2:'y', 3:'x', 'color':(0, ('red', 'blue'))}, ax)
    plt.show()

def showTiling3d():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection= '3d')
    
    ax.set_autoscale_on(True)
    ax.locator_params(integer=True)
    ax.grid(True)

    for comp in fusedstage.polyRep.polyParts:
        for polyComp in fusedstage.polyRep.polyParts[comp]:
            drawBasicSet(polyComp.sched.range(), {4:'z', 5:'x', 6:'y', 'color':(0, ('red', 'blue'))}, ax)
    plt.show()


def drawBasicSet(bset, drawSpec, axes):
    """ Visualize a basic set in two or three dimensions. The basic set should be enumerable
        i.e all the bounds must be constants. Only 3 dimensions or less can be marked as the axis 
        dimensions. Non-axis dimensions can be color coded or unrolled. Each dimension must be
        marked as axis = {x, y, z}, unroll, color. """
    # Find the dimension of the plot
    dimToAxis = { 'x' : None, 'y' : None, 'z' : None}
    colorDim = None
    for dim in drawSpec:
        if drawSpec[dim] in dimToAxis:
            if dimToAxis[drawSpec[dim]] is None:
                dimToAxis[drawSpec[dim]] = dim
            else:
                # Multiple dimensions cannot map to the same axis
                assert False       
        if dim == 'color': 
            colorDim = drawSpec[dim]
    
    if axes.name == 'rectilinear':
        assert (dimToAxis['z'] is None)

    # Check if all the dimensions are enumerable
    dimLen = bset.dim(isl._isl.dim_type.set)
    for i in xrange(0, dimLen):
       dimLowerBound = bset.dim_min(i)
       dimUpperBound = bset.dim_max(i)
       assert dimLowerBound.is_cst() and dimUpperBound.is_cst()

    points = []
    bset.foreach_point(points.append)
    #print points
    X = 0 
    Y = 0 
    Z = 0
    C = 'blue'
    if dimToAxis['x'] is not None:
        X = [ point.get_coordinate_val(isl._isl.dim_type.all, dimToAxis['x']).to_python()
              for point in points ]
    if dimToAxis['y'] is not None:    
        Y = [ point.get_coordinate_val(isl._isl.dim_type.all, dimToAxis['y']).to_python()
              for point in points ]
    if dimToAxis['z'] is not None:    
        Z = [ point.get_coordinate_val(isl._isl.dim_type.all, dimToAxis['z']).to_python() 
              for point in points ]
    if colorDim is not None:
        C = []
        for point in points:
            p = point.get_coordinate_val(isl._isl.dim_type.all, colorDim[0]).to_python()
            ringLen = len(colorDim[1])
            C.append(colorDim[1][p % ringLen])

    if axes.name == '3d':
        axes.scatter(X, Y, Z, s=50.0, c = C)
        
    if axes.name == 'rectilinear':
        axes.scatter(X, Y, s=50.0, c = C)


def drawBasicMap(bmap, drawSpec):
    pass

def drawTiles():
    pass

#islApiExamples()
#matplotlibTest()
#showTiling2d()
#showTiling3d()

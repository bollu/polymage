from __future__ import absolute_import, division, print_function

import cgen
from cexpr import *

class cNameGen(object):
    _iteratorPrefix = "_ci"
    _temporaryPrefix = "_ct"

    _iteratorCount = 0
    _tempCount = 0

    @classmethod
    def getIteratorName(cls):
        name = cls._iteratorPrefix + str(cls._iteratorCount)
        cls._iteratorCount+=1
        return name

    @classmethod
    def getTempVarName(cls):
        name = cls._temporaryPrefix + str(cls._tempCount)
        cls._tempCount+=1
        return name

class cExpression(AbstractExpression):
    pass

class cBinaryOp(AbstractBinaryOpNode):
    pass

class cUnaryOp(AbstractUnaryOpNode):
    pass

class cCond(Condition):
    pass

class cCast(cExpression):
    def __init__(self, _typ, _expr):
        assert isinstance(_typ, (cType, cPointer))
        self._typ  = _typ
        self._expr = _expr
    def __str__(self):
        return "(" + self._typ.__str__() + ") (" + self._expr.__str__() + ")"

class cName(cExpression):
    def __init__(self, _name):
        assert(isinstance(_name, str))
        self.name = _name
    def __str__(self):
        return self.name

class AbstractCgenObject(object):
    def _cgen(self):
        raise NotImplementedError
    def __str__(self):
        return self._cgen().__str__()

class cType(AbstractCgenObject):
    def __init__(self, _typ):
        self.typ = _typ
    def _cgen(self):
        return cgen.POD(self.typ, '').inline(True)

cInt = cType("int32")
cUInt = cType("uint32")
cULong = cType("uint64")
cLong = cType("int64")
cShort = cType("int16")
cUShort = cType("uint16")
cUChar = cType("uint8")
cChar = cType("int8")
cFloat = cType("float32")
cDouble = cType("float64")
cVoid = cType("void")

class TypeMap(object):
    _typeMap = { Void: cVoid, ULong:cULong, Long: cLong, UInt:cUInt, Int:cInt, 
                 UShort:cUShort, Short:cShort, UChar:cUChar, Char:cChar, 
                 Float:cFloat, Double:cDouble }
    @classmethod
    def convert(cls, typ):
        assert typ in cls._typeMap
        return cls._typeMap[typ]

class cPointer(AbstractCgenObject):
    def __init__(self, _ctype, _dim):
        assert(isinstance(_ctype, cType))
        self.typ  = _ctype.typ
        self.dim  = _dim
    def _cgen(self):
        decl = cgen.POD(self.typ, '')
        for i in xrange(0, self.dim):
            decl = cgen.Pointer(decl)
        return decl.inline(True)

class cReference(cPointer):
    def __init__(self, _ctype, _dim):
        cPointer.__init__(self, _ctype, _dim)
    def _cgen(self):
        decl = cgen.POD(self.typ, '')
        decl = cgen.Reference(decl)
        for i in xrange(0, self.dim):
            decl = cgen.Pointer(decl)
        return decl.inline(True)

class cVariable(cName):
    def __init__(self, _typ, _name):
        assert(isinstance(_typ, cType) or isinstance(_typ, cPointer))
        self.typ = _typ
        cName.__init__(self, _name)

class cDeclaration(AbstractCgenObject):
    def __init__(self, _ctyp, _cname, _expr=None):
        assert(isinstance(_ctyp, cPointer) or isinstance(_ctyp, cType))
        assert(isinstance(_cname, cName))
        self.cname = _cname
        self.ctyp = _ctyp
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_expr, AbstractExpression) or
               isinstance(_expr, basestring) or
               _expr is None)
        self.expr = _expr
    def _cgen(self):
        if self.expr is not None:
            valDecl = cgen.Value(self.ctyp.__str__(), self.cname.name).inline(True)
            return cgen.Assign(valDecl, self.expr.__str__())
        else:
            return cgen.Value(self.ctyp.__str__(), self.cname.name)

class cStatement(AbstractCgenObject):
    def __init__(self, _expr):
        self.expr = _expr
    def _cgen(self):
        return cgen.Statement(self.expr.__str__())

class cAssign(AbstractCgenObject):
    def __init__(self, _lvalue, _rvalue):
        self.lvalue = _lvalue
        self.rvalue = _rvalue
    def _cgen(self):
        return cgen.Assign(self.lvalue.__str__(), self.rvalue.__str__())

class cBlock(AbstractCgenObject):
    def __init__(self, _parent):
        assert(isinstance(_parent, AbstractCgenObject) or (_parent is None))
        self.parent = _parent
        self.block  = []
        self._isOpen = False
    def __enter__(self):
        self._isOpen = True
        return self
    def __exit__(self, typ, val, tb):
        self._isOpen = False
    def add(self, obj, checkScope = True):
        if (checkScope):
            assert(self._isOpen)
        assert(isinstance(obj, AbstractCgenObject))
        self.block.append(obj)
    def _cgen(self):
        cgenBlock = [ obj._cgen() for obj in self.block ]
        return cgen.Block(cgenBlock)

class cAbstractFCall(cExpression):
    def __init__(self, _args, _name = None, _lib = None):
        for arg in _args:
            arg = Value.numericToValue(arg)
            assert(isinstance(arg, AbstractExpression) or
                   isinstance(arg, basestring))
        self.args = _args
        self.name = _name
        self.lib  = _lib
    def __str__(self):
        arg_str = ", ".join([arg.__str__() for arg in self.args])
        return self.name + "(" + arg_str + ")"

class cFunction(cName):
    def __init__(self, _retTyp, _name, _argDict):
        cName.__init__(self, _name)
        assert(isinstance(_retTyp, cType) or isinstance(_retTyp, cPointer))
        self.retTyp = _retTyp
        self.argDict = _argDict
    def __call__(self, *args):
        call = cAbstractFCall(args, self.name)
        return call

class cSizeof(cUnaryOp):
    def __init__(self, _ctype):
        assert(isinstance(_ctype, cType) or isinstance(_ctype, cPointer))
        self._op    = 'sizeof'
        self._child = _ctype

class cLibraryFunction(cName):
    def __init__(self, _name, _lib):
        self.name = _name
        self.lib  = _lib
    def __call__(self, *args):
        call = cAbstractFCall(args, self.name, self.lib)
        return call

cMalloc = cLibraryFunction('malloc', 'stdlib.h')
cMemAlign = cLibraryFunction('memalign', 'malloc.h')
cMemSet = cLibraryFunction('memset', 'string.h')
cFree = cLibraryFunction('free', 'stdlib.h')
cPrintf = cLibraryFunction('printf', 'stdio.h')

# TODO clean up this function. 
class cFunctionDecl(AbstractCgenObject):
    def __init__(self, _func, _isExternFunc=False, _areParamsVoidPtrs=False):
        assert(isinstance(_func, cFunction))
        self.func = _func
        self.isExternFunc = _isExternFunc
        self.areParamsVoidPtrs = _areParamsVoidPtrs
    def _cgen(self):
        #argDecls = [ cgen.Value(self.func.argDict[arg].__str__(), arg.__str__())\
        #             for arg in self.func.argDict ]
        argDecls = []

        # areParamsVoidPtrs : if the target is to generate shared library
        # implementation using python ctypes

        if not self.areParamsVoidPtrs:
            for arg in self.func.argDict:
                argDecls.append(cgen.Value(self.func.argDict[arg].__str__(), arg.__str__()))
        else:
            for arg in self.func.argDict:
                # print the variable type of input and output arrays (cPointer) as
                # 'void *', so as to handle it using ctypes.c_void_p()
                if isinstance(arg.typ, cPointer):
                    argDecls.append(cgen.Value('void *', arg.__str__()+'_void_arg'))
                else:
                    argDecls.append(cgen.Value(self.func.argDict[arg].__str__(), arg.__str__()))

        typeStr = self.func.retTyp.__str__()
        if self.isExternFunc:
            typeStr = "extern \"C\" " + typeStr
        return cgen.FunctionDeclaration(cgen.Value(typeStr, 
                                        self.func.name), argDecls)

class cFunctionBody(AbstractCgenObject):
    def __init__(self, _fdecl):
        assert(isinstance(_fdecl, cFunctionDecl))
        self.fdecl = _fdecl
        self.func = _fdecl.func
        self.body  = cBlock(self)
    def _cgen(self):
        return cgen.FunctionBody(self.fdecl._cgen(), self.body._cgen())

class cReturn(AbstractCgenObject):
    def __init__(self, _expr):
        self.expr = _expr
    def _cgen(self):
        return cgen.Statement('return ' + self.expr.__str__())

class cContinue(AbstractCgenObject):
    def _cgen(self):
        return cgen.Statement('continue')

class cMacro(cName):
    def __init__(self, _name, _definition, _value):
        self.name = _name
        self.definition  = _definition
        self.value  = _value
    def __call__(self, *args):
        call = cAbstractFCall(args, self.name, None)
        return call

class cMacroDecl(AbstractCgenObject):
    def __init__(self, _macro):
        assert(isinstance(_macro, cMacro))
        self.macro = _macro
    def _cgen(self):
        return cgen.Define(self.macro.definition, self.macro.value)

cMacroMin = cMacro('isl_min', 'isl_min(x,y)', '((x) < (y) ? (x) : (y))')
cMacroMax = cMacro('isl_max', 'isl_max(x,y)', '((x) > (y) ? (x) : (y))')
cMacroFloord = cMacro('isl_floord', 'isl_floord(n,d)', '(((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))')

class cInclude(AbstractCgenObject):
    def __init__(self, _name):
        self.name = _name
    def _cgen(self):
        return cgen.Include(self.name)

class cPragma(AbstractCgenObject):
    def __init__(self, _value):
        self.value = _value
    def _cgen(self):
        return cgen.Pragma(self.value)

class cComment(AbstractCgenObject):
    def __init__(self, _text):
        self.text = _text
    def _cgen(self):
        return cgen.Comment(self.text)

class cDefine(AbstractCgenObject):
    def __init__(self, _symbol, _value):
        self.symbol = _symbol
        self.value = _value
    def _cgen(self):
        return cgen.Define(self.symbol, self.value)

class cFor(AbstractCgenObject):
    def __init__(self, _start, _cond, _update):
        _start = Value.numericToValue(_start)
        _update = Value.numericToValue(_update)
        assert(isinstance(_start, cStatement) or isinstance(_start, cAssign) 
                          or isinstance(_start, AbstractExpression)
                          or isinstance(_start, cDeclaration))
        assert(isinstance(_cond, cCond))
        assert(isinstance(_update, cStatement) or isinstance(_update, cAssign) 
                          or isinstance(_update, AbstractExpression)
                          or isinstance(_start, cDeclaration))
        self.start  = _start
        self.cond   = _cond
        self.update = _update
        self.body   = cBlock(self)
    def _cgen(self):
        start = self.start.__str__().strip(';')
        cond = self.cond.__str__()
        update = self.update.__str__().strip(';')
        return cgen.For(start, cond, update, self.body._cgen())

class cIfThen(AbstractCgenObject):
    def __init__(self, _cond):
        assert(isinstance(_cond, cCond))
        self.cond = _cond
        self.ifBlock   = cBlock(self)
    def _cgen(self):
            return cgen.If(self.cond.__str__(), self.ifBlock._cgen())

class cIfThenElse(AbstractCgenObject):
    def __init__(self, _cond):
        assert(isinstance(_cond, cCond))
        self.cond = _cond
        self.ifBlock   = cBlock(self)
        self.elseBlock = cBlock(self)
    def _cgen(self):
        return cgen.If(self.cond.__str__(), self.ifBlock._cgen(), 
                       self.elseBlock._cgen())

class cModule(AbstractCgenObject):
    def __init__(self, _name):
        self.name     = _name
        self.includes = cBlock(self)
        self.defines  = cBlock(self)
        self.decls    = cBlock(self)
        self.funcs    = cBlock(self)
    def _cgen(self):
        incsgen = [obj._cgen() for obj in self.includes.block]
        defsgen = [obj._cgen() for obj in self.defines.block]
        declsgen = [obj._cgen() for obj in self.decls.block]
        funcsgen = [obj._cgen() for obj in self.funcs.block]
        return cgen.Module(incsgen + defsgen + declsgen + funcsgen)

class cArrayAccess(cExpression):
    def __init__(self, _array, _dims):
        assert len(_array.dims) >= len(_dims)
        self.array = _array
        self.dims =  _dims
    def __str__(self):
        accessStr = ""
        if self.array.layout == 'multidim':
            for dim in self.dims:
                accessStr = accessStr + '[' + dim.__str__() + ']'
        elif self.array.layout == 'contigous':
            expr = None
            for i in xrange(0, len(self.dims)):
                multiplier = None
                for j in xrange(i+1, len(self.array.dims)):
                    if multiplier is None:
                        multiplier = self.array.dims[j]
                    else:
                        multiplier = multiplier * self.array.dims[j]
                product = self.dims[i]        
                if multiplier is not None:
                    product = product * multiplier
                if expr is None:
                    expr = product
                else:
                    expr = expr + product
            accessStr = '[' + expr.__str__() + ']'
        else:
            assert False, self.array.name
        return self.array.name + accessStr

class cArray(cName):
    def __init__(self, _typ, _name, _dims, _layout):
        assert(len(_dims) > 0)
        assert isinstance(_typ, cType)
        assert(_layout == 'contigous' or _layout == 'multidim')
        cName.__init__(self, _name) 
        self.typ = _typ
        self.dims = _dims
        self.layout = None

    def __call__(self, *dims):        
        return cArrayAccess(self, dims)

    def isConstantSize(self):
        for dim in self.dims:
            if not isConstantExpr(dim):
               return False
        return True

    def allocate_contigous(self, block):
        if self.layout == 'contigous':
            # Generate a code block which dynamically allocated memory 
            # for the given array. Single chunk allocation.
            sizeExpr = self.dims[0]
            castType = cPointer(self.typ, 1)
            for i in xrange(1, len(self.dims)):
                sizeExpr = sizeExpr * self.dims[i]
            #expr = cCast(castType, cMemAlign(64, cSizeof(self.typ) * sizeExpr))
            expr = cCast(castType, cMalloc(cSizeof(self.typ) * sizeExpr))
            block.add(cAssign(self.name, expr))
            # Counter which might come in handy
            #var = cVariable(cInt, "_c_" + self.name)
            #varDecl = cDeclaration(cInt, var, 0)
            #block.add(varDecl)

            block.add(cStatement(cMemSet(self.name, 0, cSizeof(self.typ) * sizeExpr)))
        elif self.layout == 'multidim':
            # Generate a code block which dynamically allocated memory 
            # for the given array. Multi chunk allocation.
            l = len(self.dims)
            assert(l > 0)
            castType = cPointer(self.typ, l)
            sizeType = cPointer(self.typ, l-1)
            expr = cCast(castType, cMemAlign(32, cSizeof(sizeType) * self.dims[0]))
            block.add(cAssign(self.name, expr))
            arglist = []
            for i in xrange(1, l):
                # allocate for current level
                var = cVariable(cUInt, cNameGen.getIteratorName())
                arglist.append(var)
                varDecl = cDeclaration(var.typ, var, 0)
                cond = cCond(var, '<', self.dims[i-1])
                inc = cAssign(var, var+1)
                loop = cFor(varDecl, cond, inc)
                block.add(loop, False)
               
                castType = cPointer(self.typ, l-i)
                sizeType = cPointer(self.typ, l-i-1)
                expr = cCast(castType, cMemAlign(32, cSizeof(sizeType) * self.dims[i]))
                loop.body.add(cAssign(self(*arglist), expr), False)
               
                block = loop.body

    def deallocate(self, block):
        assert self.layout == 'contigous'
        block.add(cStatement(cFree(self.name)))

class cArrayDecl(AbstractCgenObject):   
    def __init__(self, _carray):
        self.carray = _carray
        self.carray.layout = 'multidim'
    def _cgen(self):    
        decl = cgen.Value(self.carray.typ.__str__(), self.carray.name)
        for dim in self.carray.dims:
            decl = cgen.ArrayOf(decl, dim)
        return decl

# Probably defunct code beyond this point. Check and move to expression
# simplifier.
def opCountCExpr(cexpr):
    if isinstance(cexpr, AbstractBinaryOpNode):
        leftCount = opCountCExpr(cexpr.left)
        rightCount = opCountCExpr(cexpr.right)
        return leftCount + rightCount + 1
    if isinstance(cexpr, AbstractUnaryOpNode):
        childCount = opCountCExpr(cexpr.child)
        return childCount + 1
    if (isinstance(cexpr, Value) or
        isinstance(cexpr, cVariable)):
        return 0
    if isinstance(cexpr, cArrayAccess):
        # Does not count the ops to compute the array index
        opCount = 0
        for dim in cexpr.dims:
            opCount += opCountCExpr(dim)
        return opCount + 1    
    raise TypeError(type(expr))

def splitCExpr(cexpr):
    # The expressions generated by inlining the functions tend to be
    # large. So splitting the expressions into smaller computations 
    # makes the code more readable and enables the underlying compiler
    # to do better optimizations.

    # Compute the size of the expression
    pass

    




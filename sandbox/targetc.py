from __future__ import absolute_import, division, print_function

#import sys
#sys.path.insert(0, '../../cgen')

import cgen
from cexpr import *

# catch the excpetion from python3 and sort it out ahead of time
try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)
else:
    # 'unicode' exists, must be Python 2
    str = str
    unicode = unicode
    bytes = str
    basestring = basestring

class CNameGen(object):
    _iterator_prefix = "_ci"
    _temporary_prefix = "_ct"

    _iterator_count = 0
    _temp_count = 0

    @classmethod
    def get_iterator_name(cls):
        name = cls._iterator_prefix + str(cls._iterator_count)
        cls._iterator_count+=1
        return name

    @classmethod
    def get_temp_var_name(cls):
        name = cls._temporary_prefix + str(cls._temp_count)
        cls._temp_count+=1
        return name

class CValue(Value):
    def __init__(self, _value, _typ):
        Value.__init__(self, _value, _typ)
    def __str__(self):
        if (self._typ == Float):
            return self._value.__str__()
        return self._value.__str__()

class CExpression(AbstractExpression):
    pass

class CBinaryOp(AbstractBinaryOpNode):
    pass

class CUnaryOp(AbstractUnaryOpNode):
    pass

class CCond(Condition):
    pass

class CCast(CExpression):
    def __init__(self, _typ, _expr):
        assert isinstance(_typ, (CType, CPointer))
        self._typ  = _typ
        self._expr = _expr
    def __str__(self):
        return "(" + self._typ.__str__() + ") (" + self._expr.__str__() + ")"

class CName(CExpression):
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

class CType(AbstractCgenObject):
    def __init__(self, _typ):
        self.typ = _typ
    def _cgen(self):
        return cgen.POD(self.typ, '').inline(True)

c_int = CType("int32")
c_uInt = CType("uint32")
c_uLong = CType("uint64")
c_long = CType("int64")
c_short = CType("int16")
c_uShort = CType("uint16")
c_uChar = CType("uint8")
c_char = CType("int8")
c_float = CType("float32")
c_double = CType("float64")
c_void = CType("void")

class TypeMap(object):
    _type_map = { Void: c_void, ULong:c_uLong, Long: c_long, UInt:c_uInt, Int:c_int,
                 UShort:c_uShort, Short:c_short, UChar:c_uChar, Char:c_char,
                 Float:c_float, Double:c_double }
    @classmethod
    def convert(cls, typ):
        assert typ in cls._type_map
        return cls._type_map[typ]

class CPointer(AbstractCgenObject):
    def __init__(self, _ctype, _dim):
        assert(isinstance(_ctype, CType))
        self.typ  = _ctype.typ
        self.dim  = _dim
    def _cgen(self):
        decl = cgen.POD(self.typ, '')
        for i in range(0, self.dim):
            decl = cgen.Pointer(decl)
        return decl.inline(True)

class CReference(CPointer):
    def __init__(self, _ctype, _dim):
        CPointer.__init__(self, _ctype, _dim)
    def _cgen(self):
        decl = cgen.POD(self.typ, '')
        decl = cgen.Reference(decl)
        for i in range(0, self.dim):
            decl = cgen.Pointer(decl)
        return decl.inline(True)

class CVariable(CName):
    def __init__(self, _typ, _name):
        assert(isinstance(_typ, CType) or isinstance(_typ, CPointer))
        self.typ = _typ
        CName.__init__(self, _name)

class CDeclaration(AbstractCgenObject):
    def __init__(self, _ctyp, _cname, _expr=None):
        assert(isinstance(_ctyp, CPointer) or isinstance(_ctyp, CType))
        assert(isinstance(_cname, CName))
        self.cname = _cname
        self.ctyp = _ctyp
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_expr, AbstractExpression) or
               isinstance(_expr, basestring) or
               _expr is None)
        self.expr = _expr
    def _cgen(self):
        if self.expr is not None:
            val_decl = cgen.Value(self.ctyp.__str__(), self.cname.name).inline(True)
            return cgen.Assign(val_decl, self.expr.__str__())
        else:
            return cgen.Value(self.ctyp.__str__(), self.cname.name)

class CStatement(AbstractCgenObject):
    def __init__(self, _expr):
        self.expr = _expr
    def _cgen(self):
        return cgen.Statement(self.expr.__str__())

class CAssign(AbstractCgenObject):
    def __init__(self, _lvalue, _rvalue):
        self.lvalue = _lvalue
        self.rvalue = _rvalue
    def _cgen(self):
        return cgen.Assign(self.lvalue.__str__(), self.rvalue.__str__())

class CBlock(AbstractCgenObject):
    def __init__(self, _parent):
        assert(isinstance(_parent, AbstractCgenObject) or (_parent is None))
        self.parent = _parent
        self.block  = []
        self._is_open = False
    def __enter__(self):
        self._is_open = True
        return self
    def __exit__(self, typ, val, tb):
        self._is_open = False
    def add(self, obj, check_scope = True):
        if (check_scope):
            assert(self._is_open)
        assert(isinstance(obj, AbstractCgenObject))
        self.block.append(obj)
    def _cgen(self):
        cgen_block = [ obj._cgen() for obj in self.block ]
        return cgen.Block(cgen_block)

class CAbstractFCall(CExpression):
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

class CFunction(CName):
    def __init__(self, _ret_typ, _name, _arg_dict):
        CName.__init__(self, _name)
        assert(isinstance(_ret_typ, CType) or isinstance(_ret_typ, CPointer))
        self.ret_typ = _ret_typ
        self.arg_dict = _arg_dict
    def __call__(self, *args):
        call = CAbstractFCall(args, self.name)
        return call

class CSizeof(CUnaryOp):
    def __init__(self, _ctype):
        assert(isinstance(_ctype, CType) or isinstance(_ctype, CPointer))
        self._op    = 'sizeof'
        self._child = _ctype

class CLibraryFunction(CName):
    def __init__(self, _name, _lib):
        self.name = _name
        self.lib  = _lib
    def __call__(self, *args):
        call = CAbstractFCall(args, self.name, self.lib)
        return call

c_malloc = CLibraryFunction('malloc', 'stdlib.h')
c_memalign = CLibraryFunction('memalign', 'malloc.h')
c_memset = CLibraryFunction('memset', 'string.h')
c_free = CLibraryFunction('free', 'stdlib.h')
c_printf = CLibraryFunction('printf', 'stdio.h')

# TODO clean up this function. 
class CFunctionDecl(AbstractCgenObject):
    def __init__(self, _func, _is_extern_func=False, _are_params_void_ptrs=False):
        assert(isinstance(_func, CFunction))
        self.func = _func
        self.is_extern_func = _is_extern_func
        self.are_params_void_ptrs = _are_params_void_ptrs
    def _cgen(self):
        #arg_decls = [ cgen.Value(self.func.arg_dict[arg].__str__(), arg.__str__())\
        #             for arg in self.func.arg_dict ]
        arg_decls = []

        # are_params_void_ptrs : if the target is to generate shared library
        # implementation using python ctypes

        if not self.are_params_void_ptrs:
            for arg in self.func.arg_dict:
                arg_decls.append(cgen.Value(self.func.arg_dict[arg].__str__(), arg.__str__()))
        else:
            for arg in self.func.arg_dict:
                # print the variable type of input and output arrays (CPointer) as
                # 'void *', so as to handle it using ctypes.c_void_p()
                if isinstance(arg.typ, CPointer):
                    arg_decls.append(cgen.Value('void *', arg.__str__()+'_void_arg'))
                else:
                    arg_decls.append(cgen.Value(self.func.arg_dict[arg].__str__(), arg.__str__()))

        type_str = self.func.ret_typ.__str__()
        if self.is_extern_func:
            type_str = "extern \"C\" " + type_str
        return cgen.FunctionDeclaration(cgen.Value(type_str,
                                        self.func.name), arg_decls)

class CFunctionBody(AbstractCgenObject):
    def __init__(self, _fdecl):
        assert(isinstance(_fdecl, CFunctionDecl))
        self.fdecl = _fdecl
        self.func = _fdecl.func
        self.body  = CBlock(self)
    def _cgen(self):
        return cgen.FunctionBody(self.fdecl._cgen(), self.body._cgen())

class CReturn(AbstractCgenObject):
    def __init__(self, _expr):
        self.expr = _expr
    def _cgen(self):
        return cgen.Statement('return ' + self.expr.__str__())

class CContinue(AbstractCgenObject):
    def _cgen(self):
        return cgen.Statement('continue')

class CMacro(CName):
    def __init__(self, _name, _definition, _value):
        self.name = _name
        self.definition  = _definition
        self.value  = _value
    def __call__(self, *args):
        call = CAbstractFCall(args, self.name, None)
        return call

class CMacroDecl(AbstractCgenObject):
    def __init__(self, _macro):
        assert(isinstance(_macro, CMacro))
        self.macro = _macro
    def _cgen(self):
        return cgen.Define(self.macro.definition, self.macro.value)

c_macro_min = CMacro('isl_min', 'isl_min(x,y)', '((x) < (y) ? (x) : (y))')
c_macro_max = CMacro('isl_max', 'isl_max(x,y)', '((x) > (y) ? (x) : (y))')
c_macro_floord = CMacro('isl_floord', 'isl_floord(n,d)', '(((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))')

class CInclude(AbstractCgenObject):
    def __init__(self, _name):
        self.name = _name
    def _cgen(self):
        return cgen.Include(self.name)

class CPragma(AbstractCgenObject):
    def __init__(self, _value):
        self.value = _value
    def _cgen(self):
        return cgen.Pragma(self.value)

class CComment(AbstractCgenObject):
    def __init__(self, _text):
        self.text = _text
    def _cgen(self):
        return cgen.Comment(self.text)

class CDefine(AbstractCgenObject):
    def __init__(self, _symbol, _value):
        self.symbol = _symbol
        self.value = _value
    def _cgen(self):
        return cgen.Define(self.symbol, self.value)

class CFor(AbstractCgenObject):
    def __init__(self, _start, _cond, _update):
        _start = Value.numericToValue(_start)
        _update = Value.numericToValue(_update)
        assert(isinstance(_start, CStatement) or isinstance(_start, CAssign)
                          or isinstance(_start, AbstractExpression)
                          or isinstance(_start, CDeclaration))
        assert(isinstance(_cond, CCond))
        assert(isinstance(_update, CStatement) or isinstance(_update, CAssign)
                          or isinstance(_update, AbstractExpression)
                          or isinstance(_start, CDeclaration))
        self.start  = _start
        self.cond   = _cond
        self.update = _update
        self.body   = CBlock(self)
    def _cgen(self):
        start = self.start.__str__().strip(';')
        cond = self.cond.__str__()
        update = self.update.__str__().strip(';')
        return cgen.For(start, cond, update, self.body._cgen())

class CIfThen(AbstractCgenObject):
    def __init__(self, _cond):
        assert(isinstance(_cond, CCond))
        self.cond = _cond
        self.if_block   = CBlock(self)
    def _cgen(self):
            return cgen.If(self.cond.__str__(), self.if_block._cgen())

class CIfThenElse(AbstractCgenObject):
    def __init__(self, _cond):
        assert(isinstance(_cond, CCond))
        self.cond = _cond
        self.if_block   = CBlock(self)
        self.else_block = CBlock(self)
    def _cgen(self):
        return cgen.If(self.cond.__str__(), self.if_block._cgen(),
                       self.else_block._cgen())

class CModule(AbstractCgenObject):
    def __init__(self, _name):
        self.name     = _name
        self.includes = CBlock(self)
        self.defines  = CBlock(self)
        self.decls    = CBlock(self)
        self.funcs    = CBlock(self)
    def _cgen(self):
        incsgen = [obj._cgen() for obj in self.includes.block]
        defsgen = [obj._cgen() for obj in self.defines.block]
        declsgen = [obj._cgen() for obj in self.decls.block]
        funcsgen = [obj._cgen() for obj in self.funcs.block]
        return cgen.Module(incsgen + defsgen + declsgen + funcsgen)

class CArrayAccess(CExpression):
    def __init__(self, _array, _dims):
        assert len(_array.dims) >= len(_dims)
        self.array = _array
        self.dims =  _dims
    def __str__(self):
        access_str = ""
        if self.array.layout == 'multidim':
            for dim in self.dims:
                access_str = access_str + '[' + dim.__str__() + ']'
        elif self.array.layout == 'contigous':
            expr = None
            for i in range(0, len(self.dims)):
                multiplier = None
                for j in range(i+1, len(self.array.dims)):
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
            access_str = '[' + expr.__str__() + ']'
        else:
            assert False, self.array.name
        return self.array.name + access_str

class CArray(CName):
    def __init__(self, _typ, _name, _dims, _layout=None):
        assert(len(_dims) > 0)
        assert isinstance(_typ, CType)
        assert(_layout == 'contiguous' or _layout == 'multidim' or \
               _layout == None)
        CName.__init__(self, _name)
        self.typ = _typ
        self.dims = _dims
        self.layout = _layout

    def __call__(self, *dims):        
        return CArrayAccess(self, dims)

    def is_constant_size(self):
        for dim in self.dims:
            if not is_constant_expr(dim):
               return False
        return True

    def allocate_contigous(self, block, pooled):
        if self.layout == 'contigous':
            # Generate a code block which dynamically allocated memory 
            # for the given array. Single chunk allocation.
            size_expr = self.dims[0]
            cast_type = CPointer(self.typ, 1)
            for i in range(1, len(self.dims)):
                size_expr = size_expr * self.dims[i]
            #expr = CCast(cast_type, c_memalign(64, CSizeof(self.typ) * size_expr))
            expr = CCast(cast_type, c_malloc(CSizeof(self.typ) * size_expr))
            block.add(CAssign(self.name, expr))
            # Counter which might come in handy
            #var = CVariable(cInt, "_c_" + self.name)
            #var_decl = CDeclaration(cInt, var, 0)
            #block.add(var_decl)

            #block.add(CStatement(c_memset(self.name, 0, CSizeof(self.typ) * size_expr)))
        elif self.layout == 'multidim':
            # Generate a code block which dynamically allocated memory 
            # for the given array. Multi chunk allocation.
            l = len(self.dims)
            assert(l > 0)
            cast_type = CPointer(self.typ, l)
            size_type = CPointer(self.typ, l-1)
            expr = CCast(cast_type, c_memalign(32, CSizeof(size_type) * self.dims[0]))
            block.add(CAssign(self.name, expr))
            arglist = []
            for i in range(1, l):
                # allocate for current level
                var = CVariable(cUInt, CNameGen.get_iterator_name())
                arglist.append(var)
                var_decl = CDeclaration(var.typ, var, 0)
                cond = CCond(var, '<', self.dims[i-1])
                inc = CAssign(var, var+1)
                loop = CFor(var_decl, cond, inc)
                block.add(loop, False)
               
                cast_type = CPointer(self.typ, l-i)
                size_type = CPointer(self.typ, l-i-1)
                expr = CCast(cast_type, c_memalign(32, CSizeof(size_type) * self.dims[i]))
                loop.body.add(CAssign(self(*arglist), expr), False)
               
                block = loop.body

    def deallocate(self, block, pooled):
        assert self.layout == 'contigous'
        block.add(CStatement(c_free(self.name)))

class CArrayDecl(AbstractCgenObject):
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
def op_count_cExpr(cexpr):
    if isinstance(cexpr, AbstractBinaryOpNode):
        left_count = op_count_cExpr(cexpr.left)
        right_count = op_count_cExpr(cexpr.right)
        return left_count + right_count + 1
    if isinstance(cexpr, AbstractUnaryOpNode):
        child_count = op_count_cExpr(cexpr.child)
        return child_count + 1
    if (isinstance(cexpr, Value) or
        isinstance(cexpr, CVariable)):
        return 0
    if isinstance(cexpr, CArrayAccess):
        # Does not count the ops to compute the array index
        op_count = 0
        for dim in cexpr.dims:
            op_count += opcount_cExpr(dim)
        return op_count + 1
    raise TypeError(type(expr))

def split_cExpr(cexpr):
    # The expressions generated by inlining the functions tend to be
    # large. So splitting the expressions into smaller computations 
    # makes the code more readable and enables the underlying compiler
    # to do better optimizations.

    # Compute the size of the expression
    pass

    




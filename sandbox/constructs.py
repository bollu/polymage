#
# Copyright 2014-2016 Vinay Vasista, Ravi Teja Mullapudi, Uday Bondhugula,
# and others from Multicore Computing Lab, Department of Computer Science
# and Automation, Indian Institute of Science
#

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# constructs.py : Front-end design of PolyMage is defined here.
#

from __future__ import absolute_import, division, print_function

# TODO remove this at some point
from expr_ast import *
from expr_types import *
from expression import *
import logging
import targetc as genc
import math

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s")

"""
class Max(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        if typeCheck:
            _leftTyp = getType(_leftExpr)
            _rightTyp = getType(_rightExpr)
            assert _leftTyp == _rightTyp
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
   
    def getType(self):
        return getType(self._args[0])

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Max(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "isl_max(" + leftStr + ", " + rightStr + ")"  

class Min(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        if typeCheck:
            _leftTyp = getType(_leftExpr)
            _rightTyp = getType(_rightExpr)
            assert _leftTyp == _rightTyp
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
    
    def getType(self):
        return getType(self._args[0])

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Min(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "isl_min(" + leftStr + ", " + rightStr + ")"
"""

class Pow(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr):
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
   
    def getType(self):
        return Double

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Pow(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "pow(" + leftStr + ", " + rightStr + ")"

class Powf(InbuiltFunction):
    def __init__(self, _leftExpr, _rightExpr):
        InbuiltFunction.__init__(self, _leftExpr, _rightExpr)
    
    def getType(self):
        return Float

    def clone(self):
        cloneArgs = [ arg.clone() for arg in self._args ]
        return Powf(cloneArgs[0], cloneArgs[1])

    def __str__(self):
        leftStr = self._args[0].__str__()
        rightStr = self._args[1].__str__()
        return "powf(" + leftStr + ", " + rightStr + ")"

class Log(InbuiltFunction): # Natural Log
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return Double

    def clone(self):
        return Log(self._args[0].clone())

    def __str__(self):
        return "log(" +  self._args[0].__str__() +  ")"

class RandomFloat(InbuiltFunction): # Random Float b/w 0.0f and 1.0f
    def __init__(self):
        InbuiltFunction.__init__(self)

    def getType(self):
        return Float

    def clone(self):
        return RandomFloat()

    def __str__(self):
        return "(static_cast<float> (rand()) / static_cast<float> (RAND_MAX))"

class Exp(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return Double

    def clone(self):
        return Exp(self._args[0].clone())

    def __str__(self):
        return "std::exp(" +  self._args[0].__str__() +  ")"

class Sin(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return Double

    def clone(self):
        return Sin(self._args[0].clone())

    def __str__(self):
        return "std::sin(" +  self._args[0].__str__() +  ")"

class Cos(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Double

    def clone(self):
        return Cos(self._args[0].clone())

    def __str__(self):
        return "std::cos(" +  self._args[0].__str__() +  ")"

class Sqrt(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Double

    def clone(self):
        return Sqrt(self._args[0].clone())

    def __str__(self):
        return "std::sqrt(" +  self._args[0].__str__() +  ")"

class Sqrtf(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)
    
    def getType(self):
        return Float

    def clone(self):
        return Sqrtf(self._args[0].clone())

    def __str__(self):
        return "std::sqrtf(" +  self._args[0].__str__() +  ")"

class Abs(InbuiltFunction):
    def __init__(self, _expr):
        InbuiltFunction.__init__(self, _expr)

    def getType(self):
        return getType(self._args[0])

    def clone(self):
        return Abs(self._args[0].clone())

    def __str__(self):
        return "std::abs(" +  self._args[0].__str__() +  ")"

class Cast(AbstractExpression):
    def __init__(self, _typ, _expr):
        _expr = Value.numericToValue(_expr)
        assert _typ in [Float, Double, 
                        UChar, Char, 
                        UShort, Short, 
                        UInt, Int, 
                        ULong, Long]
        assert(isinstance(_expr, AbstractExpression))
        self._typ  = _typ
        self._expr = _expr

    @property
    def typ(self):
        return self._typ

    @property 
    def expression(self):
        return self._expr

    def collect(self, objType):
        objs = []
        objs += self._expr.collect(objType)
        if (type(self) == objType):
            objs += [self]
        return list(set(objs))

    def clone(self):
        return Cast(self._typ, self._expr.clone())
    
    def replace_refs(self, ref_to_expr_map):
        self._expr = substitute_refs(self._expr, ref_to_expr_map)

    def __str__(self):
        exprStr = self._expr.__str__()
        return "(" + str(genc.TypeMap.convert(self._typ)) + ") " + \
               "(" + exprStr + ")"

    def macro_expand(self):
        self._expr = self._expr.macro_expand()
        return self


class Select(AbstractExpression):
    def __init__(self, _cond, _true_expr, _false_expr, typeCheck = True):
        assert(isinstance(_cond, Condition))
        _true_expr = Value.numericToValue(_true_expr)
        _false_expr = Value.numericToValue(_false_expr)
        assert(isinstance(_true_expr, AbstractExpression))
        assert(isinstance(_false_expr, AbstractExpression))
        if typeCheck:
            trueType = getType(_true_expr)
            falseType = getType(_false_expr)
            assert trueType == falseType, str(_true_expr)+" == "+str(_false_expr)
        self._true_expr = _true_expr
        self._false_expr = _false_expr
        self._cond = _cond

    @property 
    def condition(self):
        return self._cond
    @property
    def true_expression(self):
        return self._true_expr     
    @property
    def false_expression(self):
        return self._false_expr

    def collect(self, objType):
        objs = []
        objs += self._cond.collect(objType)
        objs += self._true_expr.collect(objType)
        objs += self._false_expr.collect(objType)
        if (type(self) == objType):
            objs += [self]
        return list(set(objs))

    def replace_refs(self, ref_to_expr_map):
        self._cond.replace_refs(ref_to_expr_map)
        self._true_expr = substitute_refs(self._true_expr, ref_to_expr_map)
        self._false_expr = substitute_refs(self._false_expr, ref_to_expr_map)
    
    def clone(self):
        return Select(self._cond.clone(),
                      self._true_expr.clone(),
                      self._false_expr.clone())

    def __str__(self):
        condStr = self._cond.__str__()
        trueStr = self._true_expr.__str__()
        falseStr = self._false_expr.__str__()
        return "(" + condStr + "? " + trueStr + ": " + falseStr + ")"

    def macro_expand(self):
        self._true_expr = self._true_expr.macro_expand()
        self._false_expr = self._false_expr.macro_expand()
        return self

class Max(Select):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        Select.__init__(self, Condition(_leftExpr, '>', _rightExpr),
                                        _leftExpr, _rightExpr, typeCheck)
class Min(Select):
    def __init__(self, _leftExpr, _rightExpr, typeCheck = True):
        Select.__init__(self, Condition(_leftExpr, '<', _rightExpr), 
                                        _leftExpr, _rightExpr, typeCheck)

class Variable(AbstractExpression):
    def __init__(self, _typ, _name):
        self._name    = _name
        self._typ     = _typ
    
    @property 
    def name(self):
        return self._name

    @property 
    def typ(self):
        return self._typ

    def clone(self):
        return self

    def __str__(self):
        return self._name.__str__()

    def macro_expand(self):
        return self

class Parameter(Variable):
    def __init__(self, _typ, _name):
        Variable.__init__(self, _typ, _name)
    
class Interval(object):
    def __init__(self, _typ,  _lb, _ub):
        _lb   = Value.numericToValue(_lb)
        _ub   = Value.numericToValue(_ub)
        assert(isinstance(_lb, AbstractExpression))
        assert(isinstance(_ub, AbstractExpression))
        self._lb   = _lb
        self._ub   = _ub
        self._typ  = _typ
 
    @property 
    def lowerBound(self):
        return self._lb
    @property 
    def upperBound(self):
        return self._ub
    @property 
    def typ(self):
        return self._typ

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._lb.collect(objType)
        objs += self._ub.collect(objType)
        return list(set(objs))

    def clone(self):
        return Interval(self._typ, 
                        self._lb.clone(), 
                        self._ub.clone())

    def __str__(self):
        return '(' + self._lb.__str__() + ', ' +\
               self._ub.__str__() + ')'

class Reference(AbstractExpression):
    def __init__(self, _obj, _args):
        _args = [ Value.numericToValue(arg) for arg in _args]
        for arg in _args:
            assert(isinstance(arg, AbstractExpression))
        self._obj  = _obj
        self._args = _args

    @property
    def objectRef(self):
        return self._obj

    def _replace_ref_object(self, cloneObj):
        self._obj = cloneObj

    @property
    def arguments(self):
        return self._args

    def clone(self):
        cloneArgs = [arg.clone() for arg in self._args]
        return Reference(self._obj, cloneArgs)

    def collect(self, objType):
        objs = []
        for arg in self._args:
            objs += arg.collect(objType)
        if (type(self) is objType):
            objs += [self]
        return list(set(objs))

    def __str__(self):
        arg_str = ", ".join([arg.__str__() for arg in self._args])
        return self._obj.name + "(" + arg_str + ")"

    def macro_expand(self):
        expanded_args = []
        for arg in self._args:
            expanded_args.append(arg.macro_expand())

        self.args = expanded_args
        return self

# returns lengths of stuff one level deep in the list
# returns 0 if the object is a value
def get_inner_dimensions(lst):
    if isinstance(lst, list):
        dim = []
        for x in lst:
            if isinstance(x, list):
                dim.append(len(x))
            else:
                dim.append(0)
        return dim
    else:
        return 0


def list_elements_equal(lst):
    return not lst or lst.count(lst[0]) == len(lst)


def is_valid_kernel(kernel, num_dimensions):
    """Checks if the given kernel is a valid stencil by making sure
    that length(vardom) = nesting of kernel
    Parameters
    ----------
    kernel: list
    the kernel corresponding to the Stencil. This is usually a nested list

    num_dimensions: int
    number of dimensions the kernel posesses
    Returns
    -------
    is_valid: Bool
    """

    def is_valid_kernel_reucur(kernel, num_dimensions, closure_data):

        inner_dim = get_inner_dimensions(kernel)
        assert list_elements_equal(inner_dim), ("kernel does not have "
                                                "equal dimensions.\n"
                                                "Erroring dimensions: %s\n"
                                                "Kernel: %s" % (inner_dim,
                                                                kernel))
        if isinstance(kernel, list):
            for subkernel in kernel:
                new_closure_data = {
                    "total_dim": closure_data["total_dim"],
                    "parent": [kernel] + closure_data["parent"]
                }
                assert is_valid_kernel_reucur(subkernel,
                                              num_dimensions - 1,
                                              new_closure_data)
        else:
            if num_dimensions < 0:
                error_str = ("kernel has more dimensions than expected")
            elif num_dimensions > 0:
                error_str = ("Kernel has less dimensions than expected")
            else:
                return True

            parent_chain = " →\n\t".join(map(str, closure_data["parent"]))
            assert num_dimensions == 0, ("%s\n"
                                        "Expected Dimensions: %s\n"
                                        "Incorrect Kernel: %s\n" %
                                        (error_str,
                                         closure_data["total_dim"],
                                         parent_chain))
        return True
    
    
    # workaround for python scope madness - scopes get absolutely
    # _wrecked_ in a recursive inner function, so just pass around
    # a closure like the C people we are
    closure_data = {
        "total_dim": num_dimensions,
        "parent": []
    }
    return is_valid_kernel_reucur(kernel, num_dimensions, closure_data)

def get_valid_kernel_sizes(kernel):
    """
    Provides the sizes along the dimensions of the kernel,
    outermost to innermost for a valid kernel

    Parameters
    ----------
    kernel: nested list
    a valid N-dimensional kernel for computation
    
    Returns
    -------
    sizes: list
    1-D list of sizes, dimensions are ordered outermost to innermost
    """

    def kernel_dim_recur(subkernel):
        if isinstance(subkernel, list):
            if len(subkernel) > 0:
                return [len(subkernel)] + kernel_dim_recur(subkernel[0])
            else:
                return [0]

        else:
            return []
    return kernel_dim_recur(kernel)


def check_type(given_variable, expected_type):
    if not isinstance(given_variable, expected_type):
        raise TypeError("Expected {given_value} to be of type {expected_type}."
            "\nGiven Value: {given_value}"
            "\nGiven Type: {given_type}"
            "\nExpected Type: {expected_type}".format({
                "expected_type": str(type(given_variable)),
                "given_type": str(type(given_variable)),
                "given_value": str(given_variable)
            }))


class Stencil(AbstractExpression):
    def __init__(self, _input_fn, _iteration_vars, _kernel, _origin=None):
        check_type(_input_fn, Function)
        self._input_fn = _input_fn

        for v in _iteration_vars:
            check_type(v, Variable)
        self._iteration_vars = _iteration_vars

        assert is_valid_kernel(_kernel, len(_iteration_vars))
        self._kernel = _kernel

        self._sizes = get_valid_kernel_sizes(self._kernel)

        self._origin = _origin
        if self._origin is None:
            self._origin = list(map(lambda x: (x-1) // 2, self._sizes))

    @property
    def input_func(self):
        return self._input_fn
    @property
    def iter_vars(self):
        return self._iteration_vars
    @property
    def kernel(self):
        return self._kernel
    @property
    def sizes(self):
        return self._sizes
    @property
    def origin(self):
        return self._origin

    def __str__(self):
        return ("Stencil object"
                "\n%s"
                "\n\tinput: %s"
                "\n\titeration vars: %s"
                "\n\tdimension sizes: %s"
                "\n\torigin: %s"
                "\n\tkernel: %s" % (str(Stencil.macro_expand(self)),
                                    self._input_fn,
                                    list(map(str, self.iter_vars)),
                                    self.sizes,
                                    self._origin, self.kernel))

    @staticmethod
    def _build_indexed_kernel_recur(origin_vector, iter_vars, chosen_indeces,
                                    to_choose_sizes, subkernel):
        """
        Builds a list [([variable index], kernel weight] by taking the kernel,
        origin offset, and list of variables as parameters

        Parameters
        ----------
        origin_vector: [Int]
        the relative origin of the kernel with respect to the top left.
        If origin is (0, 0), then the kernel is built up as
        (x, y) to (x + w, y + h), since (0, 0) is taken to be the origin of the
        kernel.
        Usually, the origin is (w/2, h/2, ...)

        iter_vars: [Variable]
        Variables that represnt the iteration axes to
        index the source function
        Usually x, y, z, ...

        chosen_indeces: [Expression]
        pass "frozen" incdeces that have already been chosen. The function
        is now expected to generate all sub-indeces for these chosen
        indeces. One way to look at this is that chosen_indeces represents the
        chosen vector components of the final index.

        to_choose_sizes: [Int]
        Represents the sizes of the indeces that are yet to be chosen. Hence,
        these need to be looped over to pick _every_ index in these indeces.

        subkernel: [Int]^k (k-nested list)
        the remaining sub-space of the kernel that is yet to be chosen. Must
        be indexed from to_choose_sizes.

        Invariants
        ----------
        total kernel dimension: K
        subkernel dimension: K_s

        dim(origin_vector) = K
        len(iter_vars) = K
        len(to_choose_sizes) + len(chosen_indeces) = K

        K_s = len(to_choose_sizes)

        Returns
        -------
        indexed_kernel: [([index_expr: Expression], kernel_weight : Int)]

        Returns a list of tuples
        Each tuple has a list of indexing expressions, used to index the
        Kernel outermost to innerpost, along with the corresponding kernel
        weight at that index.
        """
        chosen = []
        for i in range(to_choose_sizes[0]):
            index_wrt_origin = iter_vars[0] + (i - origin_vector[0])
            if len(to_choose_sizes) == 1:
                chosen.append((chosen_indeces + [index_wrt_origin],
                              subkernel[i]))
            else:
                indexed = \
                    Stencil._build_indexed_kernel_recur(origin_vector[1:],
                                                        iter_vars[1:],
                                                        chosen_indeces +
                                                        [index_wrt_origin],
                                                        to_choose_sizes[1:],
                                                        subkernel[i])
                chosen.extend(indexed)
        return chosen

    @staticmethod
    def _build_indexed_kernel(origin, iter_vars, kernel):
        assert is_valid_kernel(kernel, num_dimensions=len(iter_vars))
        kernel_sizes = get_valid_kernel_sizes(kernel)
        return Stencil._build_indexed_kernel_recur(origin,
                                                   iter_vars,
                                                   [],
                                                   kernel_sizes,
                                                   kernel)

    def macro_expand(self):
        indexed_kernel = self._build_indexed_kernel(self._origin,
                                                    self._iteration_vars,
                                                    self.kernel)
        index_expr = 0
        for (indeces, weight) in indexed_kernel:
            ref = Reference(self._input_fn, indeces)
            index_expr += ref * weight

        # do this to force a pair of brackets around the entire indexing
        # expression
        # TODO: check if this is actually essential
        index_expr = 1 * index_expr

        return index_expr

class TStencil(object):
    def __init__(self, _var_domain, _kernel, _name,
                 _origin=None, _timesteps=1):

        self.name = _name
        self.var_domain = _var_domain
        self.timesteps = int(_timesteps)

        assert is_valid_kernel(_kernel, len(_var_domain))
        self.size = get_valid_kernel_sizes(_kernel)
        self.kernel = _kernel

        if _origin is None:
            self.origin = list(map(lambda x: math.floor(x / 2), self.size))

    def getObjects(self, objType):
        objs = []
        for interval in self.var_domain:
            objs += interval.collect(objType)
        return list(set(objs))

    def __str__(self):
        return ("Stencil object (%s)"
                "\n\tdomain: %s"
                "\n\tdimensions: %s"
                "\n\ttimesteps: %s"
                "\n\torigin: %s"
                "\n\tkernel: %s" % (self.name, list(map(str, self.var_domain)),
                                    self.size, self.timesteps,
                                    self.origin, self.kernel))



class Condition(object):
    def __init__(self, _left, _cond, _right):
        _left  = Value.numericToValue(_left)
        _right = Value.numericToValue(_right)
        assert(_cond in ['<', '<=', '>', '>=', '==', '!=', '&&', '||'])
        if _cond in ['<', '<=', '>', '>=', '==', '!=']:
            assert(isinstance(_left, AbstractExpression))
            assert(isinstance(_right, AbstractExpression))
        if _cond in ['&&', '||']:
            assert(isinstance(_left, Condition))
            assert(isinstance(_right, Condition))
        self._left  = _left
        self._right = _right
        self._cond  = _cond

    @property
    def lhs(self):
        return self._left
    @property
    def rhs(self):
        return self._right
    @property
    def conditional(self):
        return self._cond
   
    def clone(self):
        return Condition(self._left.clone(), self._cond, 
                         self._right.clone())

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._left.collect(objType) + self._right.collect(objType)
        return list(set(objs))

    def replace_refs(self, ref_to_expr_map):
        if(isinstance(self._left, Condition)):
            self._left.replace_refs(ref_to_expr_map)
        else:
            self._left = substitute_refs(self._left, ref_to_expr_map)
        if(isinstance(self._right, Condition)):
            self._right.replace_refs(ref_to_expr_map)
        else:
            self._right = substitute_refs(self._right, ref_to_expr_map)

    def split_to_conjuncts(self):
        conjuncts = []
        if self._cond in ['<', '<=', '>', '>=', '==']:
            conjuncts.append([self])
        elif (self._cond  == '!='):
            less_than = Condition(self._left, '<', self._right)
            conjuncts.append([less_than])
            greater_than = Condition(self._left, '>', self._right)
            conjuncts.append([greater_than])
        elif (self._cond == '||'):
            conjuncts = self._left.split_to_conjuncts() + \
                        self._right.split_to_conjuncts()
        elif (self._cond == '&&'):
            left_conjuncts = self._left.split_to_conjuncts()
            right_conjuncts = self._right.split_to_conjuncts()
            for lconjunct in left_conjuncts:
                for rconjunct in right_conjuncts:
                    conjuncts.append(lconjunct + rconjunct)
        else:
            assert False
        return conjuncts

    def __and__(self, other):
        assert(isinstance(other, Condition))
        return Condition(self, '&&', other)

    def __or__(self, other):
        assert(isinstance(other, Condition))
        return Condition(self, '||', other)
    
    def __str__(self):        
        if (self._cond is None):
            assert self._left is None and self._right is None
            return ""
        left_str = self._left.__str__()
        right_str = self._right.__str__()
        return "(" + left_str + " " + self._cond + " " + right_str + ")"


class Case(object):
    def __init__(self, _cond, _expr):
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_cond, Condition))
        assert(isinstance(_expr, (AbstractExpression, Reduce)))
        self._cond = _cond

        if isinstance(_expr, AbstractExpression):
            self._expr = _expr.macro_expand()
        else:
            self._expr = _expr

    @property
    def condition(self):
        return self._cond
    @property
    def expression(self):
        return self._expr

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        objs = self._cond.collect(objType) + self._expr.collect(objType)
        return list(set(objs))

    def replace_refs(self, ref_to_expr_map):
        self._cond.replace_refs(ref_to_expr_map)
        self._expr = substitute_refs(self._expr, ref_to_expr_map)

    def clone(self):
        return Case(self._cond.clone(), self._expr.clone())

    def __str__(self):
        return 'Case(' + self._cond.__str__() + ')' +\
                '{ ' + self._expr.__str__() + ' }'

class Op(object):
    Sum = 0
    Mul = 1
    Min = 2
    Max = 3

class Reduce(object):
    def __init__(self, _red_ref, _expr, _op_typ):
        assert isinstance(_red_ref, Reference)
        assert isinstance(_red_ref.objectRef, Reduction)
        _expr = Value.numericToValue(_expr)
        assert isinstance(_expr, AbstractExpression)
        assert _op_typ in [Op.Sum, Op.Mul, Op.Min, Op.Max]

        self._red_ref = _red_ref
        self._expr = _expr
        self._op_typ = _op_typ

    @property 
    def accumulate_ref(self):
        return self._red_ref
    @property 
    def expression(self):
        return self._expr
    @property 
    def op_type(self):
        return self._expr

    def replace_refs(self, ref_to_expr_map):
        self._expr = substitute_refs(self._expr, ref_to_expr_map)
        self._red_ref = substitute_refs(self._red_ref, ref_to_expr_map)

    def collect(self, objType):
        if (type(self) is objType):
            return [self]
        
        objs = self._red_ref.collect(objType) + self._expr.collect(objType)
        return list(set(objs))  

    def clone(self):
        return Reduce(self._red_ref.clone(), self._expr.clone(), self._op_typ)

    def __str__(self):
        op_str = None
        op_sep = None
        if (self._op_typ == Op.Sum):
            op_str = ''
            op_sep = ' + '
        elif (self._op_typ == Op.Mul):
            op_str = ''
            op_sep = ' * '
        elif (self._op_typ == Op.Min):
            op_str = 'Min'
            op_sep = ', '
        elif (self._op_typ == Op.Max):
            op_str = 'Max'
            op_sep = ', '
        else:
            assert False

        ret_str = 'Reduce [ ' + self._red_ref.__str__() + ' = ' + \
                  op_str + '(' + self._red_ref.__str__() + op_sep + \
                  self._expr.__str__() + ') ]'
        return ret_str

class Function(object):
    def __init__(self, _varDom, _typ, _name, _const=""):
        self._name      = _name
        # Type of the scalar range of the function
        self._typ       = _typ
        # Variables of the function
        self._variables = None
        # Constant function (standalone)
        if _const == "const":
            self._const = True
        else:
            self._const = False

        # Gives the domain of each variable. Domain of each variable is
        # expected to be over integers. Function evaluation in the
        # lexicographic order of the domain is assumed to be valid.

        assert(len(_varDom[0]) == len(_varDom[1]))
        for i in range(0, len(_varDom[0])):
            assert(isinstance(_varDom[0][i], Variable))
            assert(isinstance(_varDom[1][i], Interval))
            assert(_varDom[0][i].typ ==  _varDom[1][i].typ)
        # add check to ensure that upper bound and lower bound expressions
        # for each variable are only defined in terms of the function and
        # global parameters

        # Should bounds be restricted only to parameters or function variables
        # be allowed? No for now

        # Should the domain be restricted to the positive quadrant?
        # Can this be done automatically
        self._variables = _varDom[0]
        self._varDomain = _varDom[1]

        # dimensionality of the Function
        self._ndims = len(self._variables)

        # * Body of a function is composed of Case and Expression constructs.
        # * The Case constructs are expected to be non-overlapping. Therefore,
        #   value at each point in the function domain is uniquely defined.
        self._body      = []

    @property
    def name(self):
        return self._name
    @property
    def typ(self):
        return self._typ
    @property
    def is_const_func(self):
        return self._const

    @property
    def variableDomain(self):
        return (self._variables, self._varDomain)
           
    @property
    def domain(self):
        return self._varDomain

    @property
    def variables(self):
        return self._variables

    @property
    def ndims(self):
        return self._ndims

    @property
    def defn(self):
        return self._body
    @defn.setter
    def defn(self, _def):
        assert(self._body == [])
        assert(len(_def) > 0), str(_def) + " " + str(self._name)
        case_type = 0
        non_case_type = 0
        for case in _def:
            case = Value.numericToValue(case)
            assert(isinstance(case, (Case, AbstractExpression))),str(case)

            # if the function is defined using Case, all the definition parts
            # in the list '_def' must be of the type Case.
            if isinstance(case, Case):
                case_type += 1
            else:
                non_case_type += 1

            assert(non_case_type <= 1)
            assert(case_type * non_case_type == 0)

            # check if the Case and Expression constructs only use
            # function variables and global parameters

            # MOD -> if _def is not a Case, shouldnt it be disallowed after
            # the first definition?
            self._body.append(case)

    def __call__(self, *args):
        assert(len(args) == len(self._variables))
        for arg in args:
            arg = Value.numericToValue(arg)
            assert(isinstance(arg, AbstractExpression))
        return Reference(self, args)

    def replace_refs(self, ref_to_expr_map):
        num_cases = len(self._body)
        for i in range(0, num_cases):
            if isinstance(self._body[i], Case):
                self._body[i].replace_refs(ref_to_expr_map)
            else:
                self._body[i] = substitute_refs(self._body[i], ref_to_expr_map)

    def getObjects(self, objType):
        objs = []
        for case in self._body:
            objs += case.collect(objType)
        for interval in self._varDomain:
            objs += interval.collect(objType)
        return list(set(objs))

    def hasBoundedIntegerDomain(self):
        boundedIntegerDomain = True
        for varDom in self._varDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
                    break
            else:
                boundedIntegerDomain = False
                break

        return boundedIntegerDomain

    def clone(self):
        newBody = [ c.clone() for c in self._body ]
        varDom = ( [ v.clone() for v in self._variables], 
                   [ d.clone() for d in self._varDomain] )
        _const = ""
        if self.is_const_func:
            _const = "const"
        newFunc = Function(varDom, self._typ, self._name, _const)
        newFunc.defn = newBody
        return newFunc
    
    def __str__(self):
        if (self._body):
            var_str = ", ".join([var.__str__() for var in self._variables])
            dom_str = ', '.join([self._variables[i].__str__() + \
                                 self._varDomain[i].__str__()\
                                   for i in range(len(self._varDomain))])
            case_str = "{ " + "\n ".join([case.__str__() \
                                            for case in self._body]) + " }"
            return "Domain: " + dom_str + '\n' + self._name + \
                   "(" + var_str + ") = " + case_str + '\n'
        else:
            return self._name

class Image(Function):
    def __init__(self, _typ, _name, _dims):
        _dims = [ Value.numericToValue(dim) for dim in _dims ]
        # Have to evaluate if a  stronger constraint 
        # can be imposed. Only AbstractExpression in parameters?
        for dim in _dims:
            assert(isinstance(dim, AbstractExpression))
        self._dims = _dims
        intervals = []
        variables = []
        i = 0
        for dim in self._dims:
            # Just assuming it will not be more that UInt
            intervals.append(Interval(UInt, 0, dim-1))
            variables.append(Variable(UInt, "_" + _name + str(i)))
            i = i + 1
        Function.__init__(self, (variables, intervals), _typ, _name)

    @property
    def dimensions(self):
        return tuple(self._dims)

    def __str__(self):
        dim_str = ", ".join([dim.__str__() for dim in self._dims])
        return self._name.__str__() + "(" + dim_str + ")"

class Reduction(Function):
    def __init__(self, _varDom, _redDom, _typ, _name):
        Function.__init__(self, _varDom, _typ, _name)
        # Gives the domain of the reduction. Reduction domain of each variable
        # is expected to be over integers. Reduction evaluation in the
        # lexicographic order of the domain is assumed to be valid.
        assert(len(_redDom[0]) == len(_redDom[1]))
        for i in range(0, len(_redDom[0])):
            assert(isinstance(_redDom[0][i], Variable))
            assert(isinstance(_redDom[1][i], Interval))
            assert(_redDom[0][i].typ ==  _redDom[1][i].typ)
        # add check to ensure that upper bound and lower bound
        # expressions for each variable are only defined in
        # terms of variables of the function and global parameters

        # Should bounds be restricted only to parameters or function
        # variables be allowed? No for now

        # Should the domain be restricted to the positive quadrant?
        # Can this be done automatically
        self._redVariables = _redDom[0]
        self._redDomain = _redDom[1]

        # Intial value of each accumulator cell. Default is set to zero of the
        # given type
        self._default   = Value(0, _typ)

    @property
    def default(self):
        return self._default
    @default.setter
    def default(self, _expr):
        _expr = Value.numericToValue(_expr)
        assert(isinstance(_expr, AbstractExpression))
        self._default = _expr

    @property
    def reductionDomain(self):
        return self._redDomain

    @property
    def reductionVariables(self):
        return self._redVariables
        
    @property
    def defn(self):
        return self._body

    @defn.setter
    def defn(self, _def):
        assert(self._body == [])
        for case in _def:
            case = Value.numericToValue(case)
            assert(isinstance(case, (Case, Reduce))),str(case)
            # check if the Case and Expression constructs only use
            # function variables and global parameters

            # Which way is better Case inside accumulate or accumulate inside
            # Case

            # MOD -> if _def is not a Case, shouldnt it be disallowed after
            # the first definition?
            self._body.append(case)

    def hasBoundedIntegerDomain(self):
        boundedIntegerDomain = True
        for varDom in self._varDomain:
            if isinstance(varDom, Interval):
                if(not isAffine(varDom.lowerBound) or
                   not isAffine(varDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False

        for redDom in self._redDomain:
            if isinstance(redDom, Interval):
                if(not isAffine(redDom.lowerBound) or
                   not isAffine(redDom.upperBound)):
                    boundedIntegerDomain = False
            else:
                boundedIntegerDomain = False

        return boundedIntegerDomain

    def getObjects(self, objType):
        objs = []
        for case in self._body:
            objs += case.collect(objType)
        for interval in self._varDomain:
            objs += interval.collect(objType)
        for interval in self._redDomain:
            objs += interval.collect(objType)
        objs += self._default.collect(objType)
        return list(set(objs))
   
    def clone(self):
        newBody = [ r.clone() for r in self._body ]
        varDom = ( [ v.clone() for v in self._variables], 
                   [ d.clone() for d in self._varDomain] )
        redDom = ( [ r.clone() for r in self._redVariables],
                   [ d.clone() for d in self._redDomain] )
        newRed = Reduction(varDom, redDom, self._typ, self._name)
        newRed.defn = newBody
        newRed.default = self._default.clone()
        return newRed    

    def __str__(self):
        if (self._body):
            varStr = ", ".join([var.__str__() for var in self._variables])
            domStr = ', '.join([self._variables[i].__str__() + \
                                self._varDomain[i].__str__()\
                                  for i in range(len(self._varDomain))])
            redDomStr = ', '.join([self._redVariables[i].__str__() + \
                                   self._redDomain[i].__str__()\
                                     for i in range(len(self._redDomain))])
            caseStr = "{ " + "\n ".join([case.__str__() \
                                           for case in self._body]) + " }"
            return "Domain: " + domStr + '\n' + \
                   "Reduction Domain: " + redDomStr + '\n' +\
                   self._name + "(" + varStr + ") = " +\
                    caseStr + '\n' + "Default: " + self._default.__str__()
        else:
            return self._name


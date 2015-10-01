from __future__ import absolute_import, division, print_function

import sys
sys.path.insert(0, '../frontend')
sys.path.insert(0, '../codegen')

from Constructs import *
import Poly as opt

# More Python 3 vs 2 mojo
try:
    import queue
except ImportError:
    import Queue as queue

import TargetCx86 as Cgen
from collections import OrderedDict

import pygraphviz as pgv

def islExprToCgen(expr, prologueStmts = None):
    if expr.get_type() == opt.isl._isl.ast_expr_type.op:
        if (expr.get_op_type() == opt.isl._isl.ast_op_type.access or
            expr.get_op_type() == opt.isl._isl.ast_op_type.call or
            expr.get_op_type() == opt.isl._isl.ast_op_type.cond or
            expr.get_op_type() == opt.isl._isl.ast_op_type.member or 
            expr.get_op_type() == opt.isl._isl.ast_op_type.cond or
            expr.get_op_type() == opt.isl._isl.ast_op_type.select):
            assert False, expr.get_op_type()
        if expr.get_op_type() == opt.isl._isl.ast_op_type.min:
            numArgs = expr.get_op_n_arg()
            cmin = Cgen.cMin(islExprToCgen(expr.get_op_arg(0), prologueStmts),
                                  islExprToCgen(expr.get_op_arg(1), prologueStmts))
            for i in xrange(2, numArgs):
                cmin = Cgen.cMin(cmin, islExprToCgen(expr.get_op_arg(i), prologueStmts))
            if prologueStmts is not None:
                varType = Cgen.cgenType.get(Int)
                mincVar = Cgen.cVariable(varType, Cgen.cNameGen.getTempVarName())
                decl = Cgen.cDeclaration(varType, mincVar, cmin)
                prologueStmts.append(decl)
                cmin = mincVar
            return cmin
        if expr.get_op_type() == opt.isl._isl.ast_op_type.max:
            numArgs = expr.get_op_n_arg()
            cmax = Cgen.cMax(islExprToCgen(expr.get_op_arg(0), prologueStmts),
                                  islExprToCgen(expr.get_op_arg(1), prologueStmts))
            for i in xrange(2, numArgs):
                cmax = Cgen.cMax(cmax, islExprToCgen(expr.get_op_arg(i), prologueStmts))
            if prologueStmts is not None:
                varType = Cgen.cgenType.get(Int)
                maxcVar = Cgen.cVariable(varType, Cgen.cNameGen.getTempVarName())
                decl = Cgen.cDeclaration(varType, maxcVar, cmax)
                prologueStmts.append(decl)
                cmax = maxcVar
            return cmax
        if expr.get_op_type() == opt.isl._isl.ast_op_type.fdiv_q:
            assert expr.get_op_n_arg() == 2
            return Cgen.cMacroFloord(islExprToCgen(expr.get_op_arg(0), prologueStmts),
                                     islExprToCgen(expr.get_op_arg(1), prologueStmts))
        if expr.get_op_type() == opt.isl._isl.ast_op_type.add:
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) + \
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)
        if expr.get_op_type() == opt.isl._isl.ast_op_type.mul:
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) * \
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)
        if (expr.get_op_type() == opt.isl._isl.ast_op_type.div or
            expr.get_op_type() == opt.isl._isl.ast_op_type.pdiv_q):
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) / \
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)
        if expr.get_op_type() == opt.isl._isl.ast_op_type.pdiv_r:
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) % \
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)
        if expr.get_op_type() == opt.isl._isl.ast_op_type.sub:
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) - \
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)
        if expr.get_op_type() == opt.isl._isl.ast_op_type.minus:
            assert expr.get_op_n_arg() == 1
            return -islExprToCgen(expr.get_op_arg(0), prologueStmts)
        if (expr.get_op_type() == opt.isl._isl.ast_op_type.names['and'] or
            expr.get_op_type() == opt.isl._isl.ast_op_type.and_then):
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) & \
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)
        if (expr.get_op_type() == opt.isl._isl.ast_op_type.names['or'] or
            expr.get_op_type() == opt.isl._isl.ast_op_type.or_else):
            assert expr.get_op_n_arg() == 2
            return islExprToCgen(expr.get_op_arg(0), prologueStmts) |\
                   islExprToCgen(expr.get_op_arg(1), prologueStmts)

    if expr.get_type() == opt.isl._isl.ast_expr_type.int:
        return expr.get_val().to_python()
    if expr.get_type() == opt.isl._isl.ast_expr_type.id:
        return Cgen.cVariable(Cgen.cInt, expr.get_id().get_name()) 

def islCondToCgen(cond, prologueStmts = None):
    
    compDict = { opt.isl._isl.ast_op_type.eq: '==',
                opt.isl._isl.ast_op_type.ge: '>=',
                opt.isl._isl.ast_op_type.gt: '>',
                opt.isl._isl.ast_op_type.le: '<=',
                opt.isl._isl.ast_op_type.lt: '<' ,
                opt.isl._isl.ast_op_type.names['and']: '&&',
                opt.isl._isl.ast_op_type.and_then: '&&',
                opt.isl._isl.ast_op_type.names['or']: '||',
                opt.isl._isl.ast_op_type.or_else: '||' }

    assert cond.get_op_type() in compDict
    comp = compDict[cond.get_op_type()]
    assert cond.get_op_n_arg() == 2
    if (comp == '&&' or comp == '||'):
        left = islCondToCgen(cond.get_op_arg(0), prologueStmts)
        right = islCondToCgen(cond.get_op_arg(1), prologueStmts)
    else:
        left = islExprToCgen(cond.get_op_arg(0), prologueStmts)
        right = islExprToCgen(cond.get_op_arg(1), prologueStmts)

    return Cgen.cCond(left, comp, right)

def isInnerMostParallel(node):
    # Check if there are no parallel loops within the node
    noInnerParallel = True
    if node.get_type() == opt.isl._isl.ast_node_type.block:
        numNodes = (node.block_get_children().n_ast_node())
        for i in xrange(0, numNodes):
            child =  node.block_get_children().get_ast_node(i)
            noInnerParallel = noInnerParallel and isInnerMostParallel(child)
    else:
        if node.get_type() == opt.isl._isl.ast_node_type.names['for']:
            userNodes = getUserNodesInBody(node.for_get_body())
            var = islExprToCgen(node.for_get_iterator())
            if isSchedDimParallel(userNodes, var.name):
                noInnerParallel = False
        elif node.get_type() == opt.isl._isl.ast_node_type.names['if']:
            noInnerParallel = noInnerParallel and \
                              isInnerMostParallel(node.if_get_then())
            hasElse = node.if_has_else()
            if hasElse:
              noInnerParallel = noInnerParallel and \
                                isInnerMostParallel(node.if_get_else())

        elif node.get_type() == opt.isl._isl.ast_node_type.user:
            pass
        else:
            assert False
    return noInnerParallel            

def getUserNodesInBody(body):
    userNodes = []
    if body.get_type() == opt.isl._isl.ast_node_type.block:
        numNodes = (body.block_get_children().n_ast_node())
        for i in xrange(0, numNodes):
            child =  body.block_get_children().get_ast_node(i)
            userNodes += getUserNodesInBody(child)
    else:
        if body.get_type() == opt.isl._isl.ast_node_type.names['for']:
            userNodes += getUserNodesInBody(body.for_get_body())
        elif body.get_type() == opt.isl._isl.ast_node_type.names['if']:
            userNodes += getUserNodesInBody(body.if_get_then())
            hasElse = body.if_has_else()
            if hasElse:
              userNodes += getUserNodesInBody(body.if_get_else())
        elif body.get_type() == opt.isl._isl.ast_node_type.user:
            userNodes += [body]
        else:
            assert False
    return userNodes

def isSchedDimParallel(userNodes, schedDimName):
    isParallel = True
    for node in userNodes:
        part = node.user_get_expr().get_op_arg(0).get_id().get_user()
        if schedDimName not in part.parallelSchedDims:
            isParallel = False
    return isParallel

def isSchedDimVector(userNodes, schedDimName):
    isVector = True
    for node in userNodes:
        part = node.user_get_expr().get_op_arg(0).get_id().get_user()
        if schedDimName not in part.vectorSchedDim:
            isVector = False
    return isVector

def getArraysForUserNodes(userNodes, cfuncMap):
    arrays = []
    for node in userNodes:
        part = node.user_get_expr().get_op_arg(0).get_id().get_user()
        array, scratch = cfuncMap[part.comp]
        if (array not in arrays) and (True in scratch):
            arrays.append(array)
    return arrays          

def generateCNaiveFromIslAst(node, body, cfuncMap, cparamMap):
    
    if node.get_type() == opt.isl._isl.ast_node_type.block:
        numNodes = (node.block_get_children().n_ast_node())
        for i in xrange(0, numNodes):
            child =  node.block_get_children().get_ast_node(i)
            generateCNaiveFromIslAst(child, body, cfuncMap, cparamMap)
    else:
        if node.get_type() == opt.isl._isl.ast_node_type.names['for']:
            # Convert lb and ub expressions to C expressions
            prologue = []
            cond = islCondToCgen(node.for_get_cond(), prologue)
            var = islExprToCgen(node.for_get_iterator())
            inc = Cgen.cAssign(var, var + islExprToCgen(node.for_get_inc()))
            if prologue is not None:
                for s in prologue:
                    body.add(s)

            prologue = []
            init = islExprToCgen(node.for_get_init(), prologue)
            if prologue is not None:
                for s in prologue:
                    body.add(s)
            varDecl =  Cgen.cDeclaration(var.typ, var, init)
            loop = Cgen.cFor(varDecl, cond, inc)
            # Check if the loop is a parallel or a vector dimension by
            # examining the loop body.

            userNodes = getUserNodesInBody(node.for_get_body())

            dimParallel = isSchedDimParallel(userNodes, var.name)
            dimVector = isSchedDimVector(userNodes, var.name)
            arrays = getArraysForUserNodes(userNodes, cfuncMap)

                                
            if dimParallel:
                ompPragma = Cgen.cPragma("omp parallel for schedule(static)")
                body.add(ompPragma)
            if dimVector:
                vecPragma = Cgen.cPragma("ivdep") 
                body.add(vecPragma)

            body.add(loop)
            # Assuming only one parallel dimension and a whole lot
            # of other things. 
            freeList = []
            if dimParallel:
                with loop.body as lbody:
                    for array in arrays:
                        if array.isConstantSize() and False:
                        #if array.isConstantSize():
                            arrayDecl = Cgen.cArrayDecl(array) 
                            lbody.add(arrayDecl)
                        else:    
                            arrayPtr = Cgen.cPointer(array.typ, 1)
                            arrayDecl = Cgen.cDeclaration(arrayPtr, array)
                            lbody.add(arrayDecl)
                            array.allocate_contigous(lbody)
                            freeList.append(array)
            with loop.body as lbody:
                generateCNaiveFromIslAst(node.for_get_body(), lbody, 
                                         cfuncMap, cparamMap)
                # Deallocate storage
                for array in freeList:
                    array.deallocate(lbody)    

        if node.get_type() == opt.isl._isl.ast_node_type.names['if']:
            ifCond = islCondToCgen(node.if_get_cond())
            hasElse = node.if_has_else()
            if hasElse:
                cifElse = Cgen.cIfThenElse(ifCond)
                with cifElse.ifBlock as ifBlock:
                    generateCNaiveFromIslAst(node.if_get_then(), ifBlock, cfuncMap, cparamMap)
                with cifElse.elseBlock as elseBlock:
                    generateCNaiveFromIslAst(node.if_get_else(), elseBlock, cfuncMap, cparamMap)
                body.add(cifElse)
            else:    
                cif = Cgen.cIfThen(ifCond)
                with cif.ifBlock as ifBlock:
                    generateCNaiveFromIslAst(node.if_get_then(), ifBlock, cfuncMap, cparamMap)
                body.add(cif)

        if node.get_type() == opt.isl._isl.ast_node_type.user:
            # The first argument is the computation object. Retrieving the polyPart.
            polyPart = node.user_get_expr().get_op_arg(0).get_id().get_user()
            if isinstance(polyPart.expr, Accumulate):
                generateCNaiveFromAccumlateNode(node, body, cfuncMap, cparamMap)
            elif isinstance(polyPart.expr, AbstractExpression):
                generateCNaiveFromExpressionNode(node, body, cfuncMap, cparamMap)
            else:
                assert False

def cvariablesFromVariablesAndSched(node, variables, sched):
    cvarMap = {}
    for i in xrange(0, len(variables)):
        varName = variables[i].name
        dim = sched.find_dim_by_name(opt.isl._isl.dim_type.in_, 
                                     varName)
        cvarMap[variables[i]] = \
            islExprToCgen(node.user_get_expr().get_op_arg(dim+1))
    return cvarMap

def generateCNaiveFromAccumlateNode(node, body, cfuncMap, cparamMap):
    polyPart = node.user_get_expr().get_op_arg(0).get_id().get_user()
    domLen = len(polyPart.comp.reductionVariables)
    # XXX this has to be changed similar to Expression Node
    cvarMap = cvariablesFromVariablesAndSched(node, 
                 polyPart.comp.reductionVariables, polyPart.sched)
    expr = generateCExpr(polyPart.expr.expression, cparamMap, cvarMap, 
                         cfuncMap)
    prologue = []
    arrayRef = generateCExpr(polyPart.expr.accumulateRef, cparamMap, 
                             cvarMap, cfuncMap, prologueStmts = prologue)
    assign = Cgen.cAssign(arrayRef, arrayRef + expr)

    if prologue is not None:
        for s in prologue:
            body.add(s)

    if polyPart.pred:
        ccond = generateCCond(polyPart.pred, cparamMap,
                              cvarMap, cfuncMap)
        cif = Cgen.cIfThen(ccond)
        with cif.ifBlock as ifblock:
            ifblock.add(assign)
        body.add(cif)
    else:    
        body.add(assign)

def generateCNaiveFromExpressionNode(node, body, cfuncMap, cparamMap):
    polyPart = node.user_get_expr().get_op_arg(0).get_id().get_user()
    domLen = len(polyPart.comp.variables)
    # Get the mapping to the array
    array, scratch = cfuncMap[polyPart.comp]

    accScratch = [ False for i in xrange(0, domLen) ]
    for i in xrange(0, domLen):
        if i in polyPart.dimTileInfo:
            if (polyPart.dimTileInfo[i][0] != 'none'):
                accScratch[i] = True

    cvarMap = cvariablesFromVariablesAndSched(node, polyPart.comp.variables, 
                                              polyPart.sched)
    arglist = []
    scratchMap = {}
    for i in xrange(0, domLen):
        accExpr = polyPart.comp.variables[i] - \
                  polyPart.comp.domain[i].lowerBound
        if accScratch[i]:           
            varName = polyPart.comp.variables[i].name
            #dim = polyPart.sched.find_dim_by_name(opt.isl._isl.dim_type.in_, 
            #                             '_Acc_' + varName)
            #dimRem = polyPart.sched.find_dim_by_name(opt.isl._isl.dim_type.in_, 
            #                             '_Rem_' + varName)
            mulRem = polyPart.sched.find_dim_by_name(opt.isl._isl.dim_type.in_, 
                                         '_Mul_' + varName)
            #orgVar = Variable(Int, '_Acc_' + varName)
            #remVar = Variable(Int, '_Rem_' + varName)
            #cvarMap[orgVar] = islExprToCgen(node.user_get_expr().get_op_arg(dim+1)) 
            #cvarMap[remVar] = islExprToCgen(node.user_get_expr().get_op_arg(dimRem+1))
            mulVar = Variable(Int, '_Mul_' + varName)
            cvarMap[mulVar] = islExprToCgen(node.user_get_expr().get_op_arg(mulRem+1))
            scratchMap[polyPart.comp.variables[i]] = (mulVar)
            if scratch[i]:
                accExpr = (mulVar)
        arglist.append(generateCExpr(accExpr, cparamMap, cvarMap, 
                                     cfuncMap, scratchMap))
    prologue = []
    expr = generateCExpr(polyPart.expr, cparamMap, cvarMap, cfuncMap, 
                         scratchMap, prologueStmts = prologue)
    assign = Cgen.cAssign(array(*arglist), expr)

    if prologue is not None:
        for s in prologue:
            body.add(s)

    if polyPart.pred:
        ccond = generateCCond(polyPart.pred, cparamMap,
                              cvarMap, cfuncMap, scratchMap)
        cif = Cgen.cIfThen(ccond)
        with cif.ifBlock as ifblock:
            ifblock.add(assign)
        body.add(cif)    
        #var = Cgen.cVariable(Cgen.cInt, "_c_" + polyPart.comp.name)
        #inc = Cgen.cAssign(var, var + 1)
        #body.add(inc)
    else:
        body.add(assign)
        #var = Cgen.cVariable(Cgen.cInt, "_c_" + polyPart.comp.name)
        #inc = Cgen.cAssign(var, var + 1)
        #body.add(inc)

def generateCExpr(expr, cparamMap, cvarMap, cfuncMap,
                  scratchMap = {}, prologueStmts = None):
    if isinstance(expr, AbstractBinaryOpNode):
        left = generateCExpr(expr.left, cparamMap, cvarMap, 
                             cfuncMap, scratchMap, prologueStmts)
        right = generateCExpr(expr.right, cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return AbstractBinaryOpNode(left, right, expr.op)
    if isinstance(expr, AbstractUnaryOpNode):
        child = generateCExpr(expr.child, cparamMap, cvarMap, cfuncMap, 
                              scratchMap, prologueStmts)
        return AbstractUnaryOpNode(child, expr.op)
    if isinstance(expr, Value):
        return Cgen.cValue(expr.value, expr.typ)
    if isinstance(expr, Variable):
        if isinstance(expr, Parameter):
            return cparamMap[expr]
        return cvarMap[expr]
    if isinstance(expr, Reference):
        array, scratch = cfuncMap[expr.objectRef]
        numArgs = len(expr.objectRef.domain)
        shiftedArgs = []
        for i in xrange(numArgs):
            scratchArg = (expr.arguments[i] - 
                                expr.objectRef.domain[i].lowerBound)
            if scratch[i]:
                scratchArg = substituteVars(expr.arguments[i], scratchMap)
            shiftedArgs.append(simplifyExpr(scratchArg))
        args = [ generateCExpr(arg, cparamMap, cvarMap, cfuncMap, 
                               scratchMap, prologueStmts) \
                 for arg in shiftedArgs ]     
        return array(*args)
    if isinstance(expr, Select):
        ccond = generateCCond(expr.condition, cparamMap,
                              cvarMap, cfuncMap, scratchMap, prologueStmts)
        trueExpr = generateCExpr(expr.trueExpression, cparamMap, cvarMap, 
                                 cfuncMap, scratchMap, prologueStmts)
        falseExpr = generateCExpr(expr.falseExpression, cparamMap, cvarMap, 
                                  cfuncMap, scratchMap, prologueStmts)
        if prologueStmts is not None:
            varType = Cgen.cgenType.get(getType(expr.trueExpression))
            truecVar = Cgen.cVariable(varType, Cgen.cNameGen.getTempVarName())
            decl = Cgen.cDeclaration(varType, truecVar, trueExpr)
            prologueStmts.append(decl)
            
            varType = Cgen.cgenType.get(getType(expr.falseExpression))
            falsecVar = Cgen.cVariable(varType, Cgen.cNameGen.getTempVarName())
            decl = Cgen.cDeclaration(varType, falsecVar, falseExpr)
            prologueStmts.append(decl)
            
            varType = Cgen.cgenType.get(getType(expr))
            selcVar = Cgen.cVariable(varType, Cgen.cNameGen.getTempVarName())
            decl = Cgen.cDeclaration(varType, selcVar, 
                                     Cgen.cSelect(ccond, truecVar, falsecVar))
            prologueStmts.append(decl)

            return selcVar

        return Cgen.cSelect(ccond, trueExpr, falseExpr)
    if isinstance(expr, Max):
        cexpr1 = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                               cfuncMap, scratchMap, prologueStmts)
        cexpr2 = generateCExpr(expr.arguments[1], cparamMap, cvarMap, 
                               cfuncMap, scratchMap, prologueStmts)
        return Cgen.cMax(cexpr1, cexpr2)
    if isinstance(expr, Min):
        cexpr1 = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                               cfuncMap, scratchMap, prologueStmts)
        cexpr2 = generateCExpr(expr.arguments[1], cparamMap, cvarMap, 
                               cfuncMap, scratchMap, prologueStmts)
        return Cgen.cMin(cexpr1, cexpr2)
    if isinstance(expr, Pow):
        cexpr1 = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                               cfuncMap, scratchMap, prologueStmts)
        cexpr2 = generateCExpr(expr.arguments[1], cparamMap, cvarMap, 
                               cfuncMap, scratchMap, prologueStmts)
        return Cgen.cPow(cexpr1, cexpr2)
    if isinstance(expr, Powf):
        cexpr1 = generateCExpr(expr.arguments[0], cparamMap, cvarMap, cfuncMap)
        cexpr2 = generateCExpr(expr.arguments[1], cparamMap, cvarMap, cfuncMap)
        return Cgen.cPowf(cexpr1, cexpr2)
    if isinstance(expr, Exp):
        cexpr = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cExp(cexpr)
    if isinstance(expr, Sqrt):
        cexpr = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cSqrt(cexpr)
    if isinstance(expr, Sqrtf):
        cexpr = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cSqrtf(cexpr)
    if isinstance(expr, Sin):
        cexpr = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cSin(cexpr)
    if isinstance(expr, Cos):
        cexpr = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cCos(cexpr)
    if isinstance(expr, Abs):
        cexpr = generateCExpr(expr.arguments[0], cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cAbs(cexpr)
    if isinstance(expr, Cast):
        cexpr = generateCExpr(expr.expression, cparamMap, cvarMap, 
                              cfuncMap, scratchMap, prologueStmts)
        return Cgen.cCast(Cgen.cgenType.get(expr.typ), cexpr)
    raise TypeError(type(expr))

def generateCCond(cond, cparamMap, cvarMap, cfuncMap, 
                  scratchMap = {}, prologueStmts = None):
    if (cond.conditional in ['&&', '||']):
        return Cgen.cCond(generateCCond(cond.lhs, cparamMap, cvarMap, 
                                        cfuncMap, scratchMap), 
                          cond.conditional, 
                          generateCCond(cond.rhs, cparamMap, cvarMap, 
                                        cfuncMap, scratchMap, prologueStmts))
    elif(cond.conditional in ['<', '<=', '>', '>=', '==', '!=']):
        return Cgen.cCond(generateCExpr(cond.lhs, cparamMap, cvarMap, 
                                        cfuncMap, scratchMap, prologueStmts), 
                          cond.conditional, 
                          generateCExpr(cond.rhs, cparamMap, cvarMap, 
                                        cfuncMap, scratchMap, prologueStmts))
    assert False

class Stage:
    """ Stage is a part of the pipeline which realizes a set of computation
        objects. Scheduling and storage allocation for the computation objects
        is done at the level of a stage. A stage also maintains a polyhedral 
        representation of the computation objects when possible.
    """
    def __init__(self, _computeObjs, _ctx, _outputs, _paramConstraints):
        for comp in _computeObjs:
            assert(isinstance(comp, Function))

        self._computeObjs  = _computeObjs
        self._childStages = []
        self._parentStages = []
        self._outputs = _outputs
        refs = []

        for comp in self._computeObjs:
            # Filter out self references
            currRefs = comp.getObjects(Reference)
            currRefs = [ ref for ref in currRefs if not ref.objectRef == comp ]
            refs = refs + currRefs

        self._parentObjs = list(set([ref.objectRef for ref in refs \
                                if not isinstance(ref.objectRef, Image)]))
        self._inputs = list(set([ref.objectRef for ref in refs \
                            if isinstance(ref.objectRef, Image)]))
        # Create a polyhedral representation
        self._polyrep = opt.PolyRep(_ctx, self, _paramConstraints,
                                    _paramEstimates, _tileSizes,
                                    _sizeThreshold, _groupSize, _outputs)

    @property
    def childStages(self):
        return self._childStages
    @property
    def parentStages(self):
        return self._parentStages
    @property
    def parentObjs(self):
        return self._parentObjs
    @property
    def computeObjs(self):
        return self._computeObjs
    @property
    def inputs(self):
        return self._inputs
    @property
    def polyRep(self):
        return self._polyrep

    def getParameters(self):
        params = []
        for comp in self._computeObjs:
            params = params + comp.getObjects(Parameter)
        return list(set(params))

        #return list(set([comp.getObjects(Parameter) for comp in self._computeObjs]))

    def isFused(self):
        return len(self._computeObjs) > 1

    def isPolyhedral(self):
        polyhedral = True
        for comp in self._computeObjs:
            if (not comp.hasBoundedIntegerDomain()):
                polyhedral = False
        return polyhedral

    def orderComputeObjs(self):
        order = {}
        for comp in self.computeObjs:
            order[comp] = 0
        # Doing a topological sort the easy way. There is a topological
        # sort else where in the code hopefully the can be unifed at
        # some point.
        change = True
        while(change):
            change = False
            for comp in self.computeObjs:
                # filter self references
                refs = comp.getObjects(Reference)
                refs = [ref for ref in refs if not ref.objectRef == comp]
                parentObjs = list(set([ref.objectRef for ref in refs\
                                  if not isinstance(ref.objectRef, Image)]))
                for pobj in parentObjs:
                    if (pobj in order  and (order[pobj] >= order[comp])):
                        order[comp] += 1
                        change = True
        return order

    def createLoopVariables(self, variables):
        cvarMap = {}
        # Create iterator variables and bind them to DSL variables
        for i in xrange(0, len(variables)):
            varType = Cgen.cgenType.get(variables[i].typ)
            var = Cgen.cVariable(varType, Cgen.cNameGen.getIteratorName())
            # Binding variables to C iterators
            cvarMap[variables[i]] = var
        return cvarMap

    def createPerfectNestedLoop(self, body, variables, domains,
                                cfuncMap, cparamMap, cvarMap):
        lbody = body
        for i in xrange(0, len(variables)):    
            var = cvarMap[variables[i]]
            # Convert lb and ub expressions to C expressions
            lb  = generateCExpr(domains[i].lowerBound, 
                                cparamMap, cvarMap, cfuncMap)                        
            ub  = generateCExpr(domains[i].upperBound, 
                                cparamMap, cvarMap, cfuncMap)                        
            step  = generateCExpr(domains[i].step, 
                                  cparamMap, cvarMap, cfuncMap)                        
            varDecl =  Cgen.cDeclaration(var.typ, var, lb)
            comp = '<='
            if domains[i].step.value < 0:
                comp = '>='
            cond = Cgen.cCond(var, comp, ub)
            inc = Cgen.cAssign(var, var+step)
            loop = Cgen.cFor(varDecl, cond, inc)
            lbody.add(loop, False)
            lbody = loop.body
        return lbody     

    def generateFunctionScanLoops(self, obj, body, cparamMap, cfuncMap):
        # Compute function points in lexicographic order of domain
        cvarMap = self.createLoopVariables(obj.variables)

        # Generate loops. lbody is the body of the innermost loop.
        lbody = self.createPerfectNestedLoop(body, obj.variables, 
                                             obj.domain, cfuncMap, 
                                             cparamMap, cvarMap)
        arglist = obj.variables
        # Convert function definition into a C expression and add it to loop body
        #hasExpr = False
        for case in obj.definition:
            if(isinstance(case, AbstractExpression)):
                caseexpr = generateCExpr(case, cparamMap, 
                                         cvarMap, cfuncMap)
                arrayRef = generateCExpr(obj(*arglist), cparamMap,
                                         cvarMap, cfuncMap)
                assign = Cgen.cAssign(arrayRef, caseexpr)
                lbody.add(assign, False)
                #hasExpr = True
            elif(isinstance(case, Case)):
                ccond = generateCCond(case.condition, cparamMap,
                                      cvarMap, cfuncMap)
                condexpr = generateCExpr(case.expression, cparamMap,
                                         cvarMap, cfuncMap)
                cif = Cgen.cIfThen(ccond)

                if (isinstance(case.expression, AbstractExpression)):
                    arrayRef = generateCExpr(obj(*arglist), cparamMap,
                                         cvarMap, cfuncMap)
                    assign = Cgen.cAssign(arrayRef, condexpr)
                    with cif.ifBlock as ifblock:
                        ifblock.add(assign)
                        #ifblock.add(Cgen.cContinue())
                else:
                    assert False
                lbody.add(cif, False)
            else:
                assert False
        #if not hasExpr:
        #    caseexpr = generateCExpr(obj.default, cparamMap, 
        #                             cvarMap, cfuncMap)
        #    assign = Cgen.cAssign(array(*arglist), caseexpr)
        #    lbody.add(assign, False)

    def generateAccumulatorScanLoops(self, obj, body, cparamMap, cfuncMap):
        # Compute accumulator points in lexicographic order of reduction domain
        cvarMap = self.createLoopVariables(obj.reductionVariables)

        # Generate loops. lbody is the body of the innermost loop.
        lbody = self.createPerfectNestedLoop(body, obj.reductionVariables, 
                                             obj.reductionDomain, cfuncMap, 
                                             cparamMap, cvarMap)

        # Convert function definition into a C expression and add it to loop body
        #hasExpr = False
        for case in obj.definition:
            if(isinstance(case, Accumulate)):
                caseexpr = generateCExpr(case.expression, cparamMap, 
                                         cvarMap, cfuncMap)
                refArgs = case.accumulateRef.arguments
                accumRef = generateCExpr(obj(*refArgs), cparamMap,
                                             cvarMap,cfuncMap)
                assign = Cgen.cAssign(accumRef, accumRef + caseexpr)                         
                lbody.add(assign, False)
                #hasExpr = True
            elif(isinstance(case, Case)):
                ccond = generateCCond(case.condition, cparamMap,
                                      cvarMap, cfuncMap)
                condexpr = generateCExpr(case.expression, cparamMap,
                                         cvarMap, cfuncMap)
                cif = Cgen.cIfThen(ccond)

                if(isinstance(case.expression, Accumulate)):
                    refArgs = case.accumulateRef.arguments
                    accumRef = generateCExpr(obj(*refArgs), cparamMap,
                                             cvarMap,cfuncMap)
                    assign = Cgen.cAssign(accumRef, accumRef + condexpr)                         
                    with cif.ifBlock as ifblock:
                        ifblock.add(assign)
                        #ifblock.add(Cgen.cContinue())
                else:
                    assert False
                lbody.add(cif, False)
            else:
                assert False
        #if not hasExpr:
        #    caseexpr = generateCExpr(obj.default, cparamMap, 
        #                             cvarMap, cfuncMap)
        #    assign = Cgen.cAssign(array(*arglist), caseexpr)
        #    lbody.add(assign, False)    

    def generateCNaive(self, body, cfuncMap, cparamMap, outputs, isExternAlloc=False):
        self._polyrep.generateCode()
        compObjs = self.orderComputeObjs()
        sortedCompObjs = sorted(compObjs.items(), key=lambda s: (s[1], s[0].name))
        freeList = []

        for obj in [item[0] for item in sortedCompObjs]:
            if(isinstance(obj, Image)):
                continue
            # Allocate storage
            arraydim = len(obj.variables)
            # Conservative storage estimation for each dimension
            # - Check if interval step is increasing or decreasing and pick the 
            #   bound expression accordingly
            notLiveOut = not obj in outputs
            if obj in self.polyRep.polyParts:
                for part in self.polyRep.polyParts[obj]:
                   notLiveOut = notLiveOut and not part.liveout
            else:
                notLiveOut = False

            redDims = [ -1 for i in xrange(0, len(obj.domain))]
            scratch = [ False for i in xrange(0, len(obj.domain))]
            if obj in self.polyRep.polyParts and notLiveOut:
                for part in self.polyRep.polyParts[obj]:
                    for i in xrange(0, len(obj.domain)):
                        if i in part.dimScratchSize:
                            redDims[i] = max(redDims[i], part.dimScratchSize[i])
                            scratch[i] = True
            dims = []
            for i in xrange(0, len(obj.domain)):
                interval = obj.domain[i]
                if redDims[i] == -1:
                    if interval.step.value > 0 :
                        dims.append(simplifyExpr(interval.upperBound - 
                                                 interval.lowerBound + 1))
                    else:
                        dims.append(simplifyExpr(interval.lowerBound - 
                                                interval.upperBound + 1))
                else:
                    dims.append(redDims[i])

            # Simplify the expressions for the array dimensions
            arrayType = Cgen.cgenType.get(obj.typ)
            array = Cgen.cArray(arrayType, obj.name, dims)
            arrayPtr = Cgen.cPointer(arrayType, 1)
            cfuncMap[obj] = (array, scratch)

            if obj not in outputs and not notLiveOut:
                arrayDecl = Cgen.cDeclaration(arrayPtr, array)
                body.add(arrayDecl)

            if not notLiveOut:
                # do not allocate for output arrays if they are
                # already allocated
                if obj in outputs and isExternAlloc:
                    # FIXME: accessing a class member directly,
                    # as a makeshift alternative
                    array.layout = 'contigous'
                    pass
                else:
                    array.allocate_contigous(body)

            if obj not in outputs and not notLiveOut:
                freeList.append(array)

            if obj in self.polyRep.polyParts:
                continue
            if type(obj) == Function:
                self.generateFunctionScanLoops(obj, body, cparamMap, cfuncMap)
            elif type(obj) == Accumulator:
                self.generateAccumulatorScanLoops(obj, body, cparamMap, cfuncMap)
            else:
                assert False

        if self.polyRep.polyast != []:
            for ast in self.polyRep.polyast:
                generateCNaiveFromIslAst(ast, body, cfuncMap, cparamMap)

        # Deallocate storage
        for array in freeList:
            array.deallocate(body)

    def __str__(self):
        comp_str  = "\n\n".join([comp.__str__() for comp in self._computeObjs]) + '\n'
        parent_str = "Parent Objects: " + \
                     ", ".join([ item.name for item in self._parentObjs]) + '\n' 
        return comp_str + '\n' + parent_str + '\n' + self._polyrep.__str__()

class PipeLine:
    def __init__(self, _outputs, _ctx, _paramConstraints, 
                 _paramEstimates, _tileSizes, _sizeThreshold,
                 _groupSize, _inlineDirectives):
        _name = ""
        for out in _outputs:
            _name = _name + out.name
        self._name   = _name
        self._ctx    = _ctx
        self._paramConstraints = _paramConstraints
        self._paramEstimates = _paramEstimates
        self._tileSizes = _tileSizes
        self._sizeThreshold = _sizeThreshold
        self._groupSize = _groupSize
        self._inlineDirectives = _inlineDirectives

        self._outputs = _outputs
        self._stages = self.buildStageGraph()
        self._intialGraph = self.drawPipelineGraph()

        # Bounds Checking
        self.boundsCheckPass()
        # Simple Inlining
        self.inlinePass()
        # Storage Reuse
        # Fusion 
        fusedStage = Stage([f for f in self._stages], _ctx, _paramConstraints,
                            _paramEstimates, _tileSizes, _sizeThreshold, 
                            _groupSize, _outputs)
        self._final = fusedStage
        self._stages = {}
        for out in _outputs:
            self._stages[out] = fusedStage

        inputs = []
        for s in self._stages:
            inputs = inputs + self._stages[s].inputs

        self._inputs = list(set(inputs))

    @property
    def stages(self):
        return self._stages

    @property
    def name(self):
        return self._name

    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs

    @property
    def islContext(self):
        return self._ctx

    @property
    def graph(self):
        return self._intialGraph

    def getParameters(self):
        params=[]
        for stage in self._stages.values():
            params = params + stage.getParameters()
        return list(set(params))

    def drawPipelineGraph(self):
        G = pgv.AGraph(strict=False, directed=True)
        for f in self._stages:
            s = self._stages[f] 
            assert len(s.computeObjs) == 1
            for obj in s.computeObjs:
                for pobj in s.parentObjs:
                    G.add_edge(pobj.name, obj.name)

        G.layout(prog='dot')
        return G

    def drawPipelineGraphWithGroups(self):
        G = pgv.AGraph(strict=False, directed=True)
        for f in self._stages:
            s = self._stages[f]
            for i in xrange(0, len(s.polyRep.groups)):
                subGraphNodes = []
                for p in s.polyRep.groups[i]:
                    if p.comp.name not in subGraphNodes:
                        subGraphNodes.append(p.comp.name)
                G.add_nodes_from(subGraphNodes)
                G.add_subgraph(nbunch = subGraphNodes, 
                               name = "cluster_" + str(i))

        # Temporary hack getting edges from the orginal graph
        # this may not reflect the current pipeline due to 
        # inlining.
        G.add_edges_from(self.graph.edges())
        G.layout(prog='dot')
        return G

    def buildStageGraph(self):
        # Find all the computation objects that are required for the pipeline,
        # create stages for the computation and the depedency graph. This step 
        # assumes that there are no cycles in the pipeline stage graph this has 
        # assumption has to be revisited when dealing with more complex pipelines.
        stages = {}
        q = queue()
        for compObj in self._outputs:
            q.put(compObj)
        while not q.empty():
            obj = q.get()
            if obj not in stages:
                stages[obj] = Stage([obj], self._ctx, self._paramConstraints,
                                    self._paramEstimates, self._tileSizes,
                                    self._sizeThreshold, self._groupSize, 
                                    self._outputs)                
                if len(stages[obj].parentObjs) != 0:
                    for r in stages[obj].parentObjs:
                        q.put(r)
        
        for obj in stages:
            for pobj in stages[obj].parentObjs:
                stages[pobj].childStages.append(stages[obj])
                stages[obj].parentStages.append(stages[pobj])
        
        return stages

    def boundsCheckPass(self):
        """ Bounds check pass analyzes if function values used in the compute 
            objects are within the domain of the functions. Static analysis 
            is only possible when the references to function values are regular
            i.e. they are not data dependent. We restrict ourselves to affine
            references."""
        for stage in self._stages.values():
            for child in stage.childStages:
                opt.checkRefs(child, stage)              
            for inp in stage.inputs:
                # Creating a computation stage for an input which is given
                # is meaningless. Ideally it should be done in a clean way
                # currently abusing stage for construction of a polyhedral
                # representation
                inpStage = Stage([inp], self._ctx, self._paramConstraints,
                                 self._paramEstimates, self._tileSizes,
                                 self._sizeThreshold, self._groupSize,
                                 self._outputs)
                opt.checkRefs(stage, inpStage)

    def inlinePass(self):
        """ Inline pass takes all the inlining decisions and inlines functions 
            at their use points in other stages."""
        # Guidelines for inlining
        # -- Inlining functions which have multiple case constructs can quickly 
        #    turn out to be messy code generation nightmare.
        # -- Splitting an use site can occur when the inlined function is defined 
        #    by different expression in different parts of the function domain. 
        #    Functions which are defined by a single expression over the entire 
        #    domain are good candidates for inlining.
        #    Example :
        #    f(x) = g(x-1) + g(x) + g(x+1) when (1 <= x <= N-1)
        #         = 0                      otherwise
        #    h(x) = f(x-1) + f(x) + f(x+1) when (1 <= x <= N-1)
        #         = 0                      otherwise
        #    Inlining f into h will result in splitting as shown below
        #    h(x) = g(x-2) + g(x-1) + g(x) +      when (2 <= x <= N-2) 
        #           g(x-1) + g(x) + g(x+1) +
        #           g(x)   + g(x+1) + g(x+2)
        #         = 2*g(x-1) + 3*g(x) + 2*g(x+1)  when (x == 1) 
        #           + g(3)      
        #         = g(x-2) + 2*g(x-1) +           when (x == N-1)
        #           3*g(x) + 2*g(x+1)
        #         = 0                             otherwise
        #    For multiple dimensions and complex expression it gets ugly fast
        # -- Inlining without splitting is possible even when the function is 
        #    defined by mutliple expressions. When the function references at
        #    the use site can all be proved to be generated from a single 
        #    expression.
        #    Example :
        #    f(x) = g(x-1) + g(x) + g(x+1) when (1 <= x <= N-1)
        #         = 0                      otherwise
        #    h(x) = f(x-1) + f(x) + f(x+1) when (2 <= x <= N-2)
        #         = 0                      otherwise
        #    Inlining f into h will not result in splitting as shown below
        #    h(x) = g(x-2) + g(x-1) + g(x) +      when (2 <= x <= N-2) 
        #           g(x-1) + g(x) + g(x+1) +
        #           g(x)   + g(x+1) + g(x+2)
        #         = 0                             otherwise
        # -- Users should be enouraged to write functions in a form that makes
        #    inlining easy.
        inlinedCompObjs = []
        for directive in self._inlineDirectives:
            # Only function constructs can be inlined for now
            assert isinstance(directive, Function)
            # Does inling into a fused stage cause problems?
            parentStage = self._stages[directive]
            assert parentStage.computeObjs[0] not in self._outputs
            for child in parentStage.childStages:
                refToInlineExprMap = opt.inline(child, parentStage, noSplit = True)
                child.computeObjs[0].inlineRefs(refToInlineExprMap)
            # Recompute stage graph
            self._stages = self.buildStageGraph()

    def getOrderedStages(self):
        # Topological sorting of stages.
        # Quick and dirty way of doing things have to revisit
        stageOrder = {}
        stageList = [ (len(s.parentStages), s) for s in self._stages.values()]
        level = 0
        while len(stageList) > 0:
            levelAssigned = [ s for s in stageList if s[0] == 0]
            for assigned in levelAssigned:
                stageOrder[assigned[1]] = level
                stageList.remove(assigned)
                def decChild(s):
                    if s[1] in assigned[1].childStages:
                        return (s[0]-1, s[1])
                    else: 
                        return s
                stageList = [ decChild(s) for s in stageList ] 
            level = level + 1
        return sorted(stageOrder.items(), key=lambda s: s[1])
   
    def generateCNaive(self, isExternFunc=False, areParamsVoidPtrs=False, isExternAlloc=False):
        # Order stages
        stageOrder = self.getOrderedStages()
        stageOrder = [ stage[0] for stage in stageOrder ]
        # Generate header and define block
        m = Cgen.cModule('Pipeline')
        with m.includes as incblock:
            incblock.add(Cgen.cInclude('stdio.h'))
            incblock.add(Cgen.cInclude('stdlib.h'))
            incblock.add(Cgen.cInclude('malloc.h'))
            incblock.add(Cgen.cInclude('cmath'))
            incblock.add(Cgen.cInclude('string.h'))
            incblock.add(Cgen.cMacroDecl(Cgen.cMacroMin))
            incblock.add(Cgen.cMacroDecl(Cgen.cMacroMax))
            incblock.add(Cgen.cMacroDecl(Cgen.cMacroFloord))
        with m.funcs as funcblock:
            # Collect all the Inputs and Parameters of the pipeline. 
            # Add them as parameters to the pipeline function.
            # TODO:
            # Also set these data as pipeline properties, since they are
            # required again later during autotuning.
            params = []
            for s in stageOrder:
                for obj in s.computeObjs:
                    params = params + obj.getObjects(Parameter)

            pipelineParams = OrderedDict()
            params = list(set(params))

            cparamMap = {}
            cfuncMap  = {}

            # 1. Collecting scalar parameters
            params.sort(key=lambda x: x.name)
            for param in params:
                cvar = Cgen.cVariable(Cgen.cgenType.get(param.typ), param.name)
                if param.definition is None:
                    pipelineParams[cvar] = cvar.typ
                else:
                    assert False
                # Binding parameters to C variables
                cparamMap[param] = cvar

            # 2. Collecting image parameters
            self.inputs.sort(key=lambda x: x.name)
            for img in self.inputs:
                cpoint = Cgen.cPointer(Cgen.cgenType.get(img.typ), 1)
                cvar = Cgen.cVariable(cpoint, img.name)
                pipelineParams[cvar] = cvar.typ
                # Binding array parameters to C array
                carray = Cgen.cArray(Cgen.cgenType.get(img.typ),
                                     img.name, img.dimensions)
                cfuncMap[img] = (carray,
                                 [ False for i in xrange(0, len(img.dimensions))])
                carray.layout = 'contigous'

            # 3. Collecting output parameters
            self.outputs.sort(key=lambda x: x.name)

            # areParamsVoidPtrs : if the target is to generate shared library
            # implementation using python ctypes

            # isExternAlloc : is the result array allocated outside the
            # generated implementation

            if not isExternAlloc:
                for func in self.outputs:
                    cref = Cgen.cReference(Cgen.cgenType.get(func.typ), 1)
                    cvar = Cgen.cVariable(cref, func.name)
                    pipelineParams[cvar] = cvar.typ
            else:
                for func in self.outputs:
                    cpoint = Cgen.cPointer(Cgen.cgenType.get(func.typ), 1)
                    cvar = Cgen.cVariable(cpoint, func.name)
                    pipelineParams[cvar] = cvar.typ

            pipeline = Cgen.cFunction(Cgen.cVoid, 'pipeline_' + self.name, pipelineParams)
            pipelineDecl = Cgen.cFunctionDecl(pipeline, isExternFunc, areParamsVoidPtrs)
            pipelineBody = Cgen.cFunctionBody(pipelineDecl)

            funcblock.add(pipelineBody)

            with pipelineBody.body as pbody:
                if areParamsVoidPtrs:
                    # Add assignment to inputs
                    for inp in self.inputs:
                        # actual input to be used
                        varType = Cgen.cgenType.get(inp.typ)
                        varPtr = Cgen.cPointer(varType, 1)
                        var = Cgen.cVariable(varType, inp.name)
                        varDecl = Cgen.cDeclaration(varPtr, var)
                        pbody.add(varDecl)

                        # dummy void * input taken as argument
                        dummyType = Cgen.cgenType.get(inp.typ)
                        dummyPtr = Cgen.cPointer(dummyType, 1)
                        dummyCast = Cgen.cCast(dummyPtr, inp.name+'_void_arg')
                        varAssign = Cgen.cAssign(var, dummyCast)
                        pbody.add(varAssign)

                    # Add assignment to output
                    for out in self.outputs:
                        # actual output to be used
                        varType = Cgen.cgenType.get(out.typ)
                        varPtr = Cgen.cPointer(varType, 1)
                        var = Cgen.cVariable(varType, out.name)
                        varDecl = Cgen.cDeclaration(varPtr, var)
                        pbody.add(varDecl)

                        # dummy void * output taken as argument
                        dummyType = Cgen.cgenType.get(out.typ)
                        dummyPtr = Cgen.cPointer(dummyType, 1)
                        dummyCast = Cgen.cCast(dummyPtr, out.name+'_void_arg')
                        varAssign = Cgen.cAssign(var, dummyCast)
                        pbody.add(varAssign)

                for s in stageOrder:
                    s.generateCNaive(pbody, cfuncMap, cparamMap, self.outputs, isExternAlloc)

        return m

    def __str__(self):
        return_str = "Final Stage: " + self._name + "\n"
        for s in self._stages:
            return_str = return_str + s.__str__() + "\n"
        return return_str

def buildPipeLine(outputs, paramConstraints = [], paramEstimates = [], 
                  tileSizes = [64, 64, 64, 64, 64, 64, 64],
                  sizeThreshold = 200 * 200 * 2,
                  groupSize = 100,
                  inlineDirectives = []):
    ctx = opt.isl.Context()
    # Both paramConstraints and inlineDirectives should be sanitized here
    # and proper errors should be reported

    # Check that paramEstimates amd paramConstraints do not conflict
    
    return PipeLine(outputs, ctx, paramConstraints, paramEstimates, 
                    tileSizes, sizeThreshold, groupSize, inlineDirectives)

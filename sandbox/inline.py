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
# inline.py : Function inline pass.
#

from __future__ import absolute_import, division, print_function

from constructs import *
from poly import *
import pipe

def inline_pass(pipeline):
    """
    Inline pass takes all the inlining decisions and inlines functions at
    their use points in other groups.
    """
    # Guidelines for inlining
    # -- Inlining functions which have multiple case constructs can quickly
    #    turn out to be messy code generation nightmare.
    # -- Splitting a use site can occur when the inlined function is defined
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
    # -- Users should be encouraged to write functions in a form that makes
    #    inlining easy.
    inlined_funcs = []
    # TODO: _inline_directives flag
    for directive in pipeline._inline_directives:
        # Only function constructs can be inlined for now
        assert isinstance(directive, Function)

        clone = pipeline._clone_map[directive]
        comp = pipeline._func_map[clone]
        group = comp.group  # the only comp of group

        # One simply does not walk into Inlining
        assert group.comps[0].func not in pipeline.outputs

        drop_inlined = True
        inlined = []
        for child in group.children:
            ref_to_inline_expr_map = \
                piecewise_inline_check(child, group, no_split=True)
            if ref_to_inline_expr_map:
                inline_and_update_graph(pipeline, group, child,
                                        ref_to_inline_expr_map)
            else:
                # for this child inlining did not happen, hence do not drop
                # the compute object and its group
                drop_inlined = False

        if drop_inlined:
            # remove comp_obj
            pipeline.drop_comp(comp)
            pipeline.drop_group(group)
            pipeline._clone_map.pop(directive)

    return

def inline_and_update_graph(pipeline, group, child_group, inline_map):
    """
    calls the actual reference replacements methods at the low level and
    updates the compute graph and group graph accordingly
    """
    comp = group.comps[0]
    func = comp.func
    child_comp = child_group.comps[0]
    child_func = child_comp.func

    # replace the references of the child to inline the comp
    child_func.replace_refs(inline_map)
    pipeline.make_func_independent(comp.func, child_comp.func)

    # create new Group for the child
    new_group = pipe.Group(pipeline._ctx, [child_comp],
                           pipeline._param_constraints)

    pipeline.replace_group(child_group, new_group)

    return

def piecewise_inline_check(child_group, parent_group, no_split = False):
    ref_to_inline_expr_map = {}
    # Inling currently only handles non-fused stages
    assert (not child_group.is_fused() and not parent_group.is_fused())
    # Computation object in the parent group can only be a function
    parent_comp = parent_group.comps[0]
    child_comp = child_group.comps[0]

    # Simple scalar functions which are defined on a non-bounded 
    # integer domain can be inlined.
    # TODO

    # Inlining only if both child and parent stage have a polyhedral
    # representation
    if child_group.polyRep.poly_parts and parent_group.polyRep.poly_parts:
        child_parts = child_group.polyRep.poly_parts
        parent_parts = parent_group.polyRep.poly_parts
        parent_doms = parent_group.polyRep.poly_doms
        for child_part in child_parts[child_comp]:
            # - Collect all the refs and pick out refs to only parent_comp
            child_refs = child_part.refs
            child_refs = [ ref for ref in child_refs \
                                 if ref.objectRef == parent_comp.func ]
            # Compute dependence relations between child and parent
            deps = []
            for ref in child_refs:
                deps += extract_value_dependence(child_part, ref,
                            parent_doms[parent_comp])
            # Check if all the values come from the same parent part
            dep_to_part_map = {}
            for dep in deps:
                # total accessible region by this dependence
                access_region = dep.rel.range().copy().reset_tuple_id()
                # holds the remaining available region after each piece's
                # region is cut off from the access region
                left_over = dep.rel.range().copy().reset_tuple_id()

                # analyze for each parent piece
                for parent_part in parent_parts[parent_comp]:
                    # access region of this part
                    part_region = \
                        parent_part.sched.domain().copy().reset_tuple_id()
                    # portion of the access_region touched by this part's
                    # region
                    partdiff = access_region.subtract(part_region)
                    # chop off the part's region from the left-over region
                    left_over = left_over.subtract(part_region)
                    if(partdiff.is_empty()):
                        dep_to_part_map[dep] = parent_part
                if (not left_over.is_empty()):
                    assert False, "Inlining cannot be done."

            parts = list(set(dep_to_part_map.values()))
            single_part = (len(parts) == 1)

            if(single_part):
                parent_expr = parts[0].expr
                if parts[0].pred:
                    inline = False
                else:
                    inline_deps = []
                    for ref in child_refs:
                        ref_to_inline_expr_map[ref] = parent_expr
            elif(no_split):
                pass
            else:
                pass
    else:
        pass

    return ref_to_inline_expr_map


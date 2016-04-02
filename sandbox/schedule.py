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
# schedule.py : Scheduling the pipeline functions, groups for code generation.
#

from __future__ import absolute_import, division, print_function

import logging
from utils import *
from poly import *

# LOG CONFIG #
schedule_logger = logging.getLogger("schedule.py")
schedule_logger.setLevel(logging.INFO)
LOG = schedule_logger.log

def naive_sched_objs(order):
    # get a reverse map for the object order map
    reverse_map = {}
    max_level = 0
    for obj in order:
        l = order[obj]
        if l not in reverse_map:
            reverse_map[l] = [obj]
        else:
            reverse_map[l] += [obj]
        max_level = max(max_level, l)

    naive_order = {}
    time = 0
    for l in range(0, max_level+1):
        # schedule within the current level
        for obj in reverse_map[l]:
            naive_order[obj] = time
            time += 1
        # sort the next level beforehand
        if l != max_level:
            next_level = []
            for obj in reverse_map[l]:
                obj_kids = [kid for kid in reverse_map[l+1] \
                                   if kid in obj.children]
                next_level += [kid for kid in obj_kids \
                                     if kid not in next_level]
            next_level += [obj for obj in reverse_map[l+1] \
                                  if obj not in next_level]
            reverse_map[l+1] = next_level

    return naive_order

def sort_scheduled_objs(schedule):
    '''
    sorts the objects according to the schedule and returns a list of the
    sorted objects
    '''
    sorted_schedule = sorted(schedule.items(), key=lambda x:x[1])
    sorted_objs = [obj[0] for obj in sorted_schedule]

    return sorted_objs

def naive_sched_groups(pipeline):
    '''
    schedule in level order traversal of the group DAG
    '''
    level_order = pipeline.get_ordered_groups

    return naive_sched_objs(level_order)

def schedule_groups(pipeline):
    # naive scheduling
    grp_schedule = naive_sched_groups(pipeline)

    return grp_schedule


def schedule_parts(group, sorted_comps):
    '''
    Computations which have different scale but map to the same time
    generate a lot of conditionals which can hinder performance. This
    step separates all computations in a time step by adding an additional
    dimension.
    '''
    part_comp_map = group.polyRep.poly_parts
    pi = 0
    for comp in sorted_comps:
        for p in part_comp_map[comp]:
            time_dim = \
              p.sched.find_dim_by_name(isl._isl.dim_type.out, '_t')
            p.sched = \
              p.sched.insert_dims(isl._isl.dim_type.out, time_dim + 1, 1)
            p.sched = p.sched.set_dim_name(isl._isl.dim_type.out,
                                           time_dim + 1, '_o')

            eqs = []
            coeff = {}
            coeff[('constant', 0)] = -pi
            coeff[('out', time_dim + 1)] = 1
            eqs.append(coeff)
            p.sched = add_constraints(p.sched, [], eqs)
            pi += 1

    return

def naive_sched_comps(group):
    '''
    schedule in level order traversal of the group comps DAG
    '''
    level_order = group.get_ordered_comps

    return naive_sched_objs(level_order)

def schedule_within_group(group):
    '''
    schedule the compute objects within each group, possibly for a better
    storage mapping
    '''

    # naive scheduling
    comp_schedule = naive_sched_comps(group)

    # get list of comps sorted according to scheduled order
    sorted_comps = get_sorted_objs(comp_schedule)

    # create a sub-time dimension for poly_parts to introduce an order within
    # the group
    schedule_parts(group, sorted_comps)

    return comp_schedule

def schedule_liveouts(pipeline):
    '''
    Schedule of the group liveouts is the group schedule itself
    '''
    liveouts = pipeline.liveouts
    grp_schedule = pipeline.group_schedule
    comps_schedule = {}
    for comp in liveouts:
        comps_schedule[comp] = grp_schedule[comp.group]

    return comps_schedule

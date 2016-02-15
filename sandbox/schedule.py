from __future__ import absolute_import, division, print_function

import logging
from utils import *
from poly import *

# LOG CONFIG #
schedule_logger = logging.getLogger("schedule.py")
schedule_logger.setLevel(logging.INFO)
LOG = schedule_logger.log

def compute_liveness(children_map, schedule):
    '''
    Given a schedule for a DAG of compute objects, this function computes the
    liveness range of each compute object. The output is a mapping from the
    timestamp in the schedule to a list of compute objects that are live only
    upto the corresponding time.
    '''
    liveness_map = {}
    for comp in schedule:
        last_live = -1
        if comp in children_map:
            for child in children_map[comp]:
                t = schedule[child]
                last_live = max(last_live, t)
            if last_live not in liveness_map:
                liveness_map[last_live] = [comp]
            else:
                liveness_map[last_live] += [comp]

    return liveness_map

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
    return level_order

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
    return level_order

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

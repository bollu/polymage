from __future__ import absolute_import, division, print_function

import logging
from utils import *
from poly import *

# LOG CONFIG #
liveness_logger = logging.getLogger("liveness.py")
liveness_logger.setLevel(logging.INFO)
LOG = liveness_logger.log

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

def liveness_for_group_comps(group, children_map, comps_schedule):
    '''
    Compute liveness for Group Compute Objects
    '''
    liveness_map = compute_liveness(group.children_map, group.comps_schedule)

    group.set_liveness_map(liveness_map)

    return

def liveness_for_pipe_outputs(pipeline):
    '''
    Compute liveness for Group Liveout0s
    '''
    def logs(liveness_map):
        # ***
        log_level = logging.DEBUG-1
        LOG(log_level, "\n_______")
        LOG(log_level, "Liveness map for Liveouts:")
        for time in liveness_map:
            LOG(log_level, str(time)+":"+\
                str([comp.func.name for comp in liveness_map[time]]))
        return

    # compute liveness
    liveness_map = \
        compute_liveness(pipeline.liveouts_children_map,
                         pipeline.liveouts_schedule)

    logs(liveness_map)

    return liveness_map


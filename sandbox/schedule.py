from __future__ import absolute_import, division, print_function

import logging
from utils import *

# LOG CONFIG #
schedule_logger = logging.getLogger("schedule.py")
schedule_logger.setLevel(logging.INFO)
LOG = schedule_logger.log

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

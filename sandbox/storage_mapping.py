from __future__ import absolute_import, division, print_function

import logging
from expression import *
from pipe import *

# LOG CONFIG #
pipe_logger = logging.getLogger("storage_mapping.py")
pipe_logger.setLevel(logging.DEBUG-1)
LOG = pipe_logger.log

def compute_liveness(pipeline, schedule):
    '''
    Given a schedule for a DAG of compute objects, this function computes the
    liveness range of each compute object. The output is a mapping from the
    timestamp in the schedule to a list of compute objects that are live only
    upto the corresponding time.
    '''
    children_map = pipeline.comp_objs_children
    liveness_map = {}
    for comp in schedule:
        children = children_map[comp]
        last_live = 0
        for child in children:
            t = schedule[child]
            last_live = max(last_live, t)
        if last_live not in liveness_map:
            liveness_map[last_live] = [comp]
        else:
            liveness_map[last_live].append([comp])

    liveness_map = sorted(liveness_map.items(), key=lambda x:x[1])

    return liveness_map

class Storage:
    class Dimension:
        def __init__(self, size):
            self._param = size[0]
            self._size_expr = size[1]

            coeff_map = get_affine_var_and_param_coeff(self._size_expr)
            self._coeff = coeff_map[self._param]
            self._const = get_constant_from_expr(self._size_expr)

        @property
        def param(self):
            return self._param
        @property
        def size(self):
            return self._size_expr
        @property
        def coeff(self):
            return self._coeff
        @property
        def const(self):
            return self._const

        def is_constant(self):
            return self._param == 0

    def __init__(self, _dims, _dim_sizes):
        self._dims = _dims
        self._dim_sizes = _dim_sizes
        self._dimension = []
        for dim in range(0, dims):
            self._dimension.append(Dimension(self._dim_sizes[dim]))

    @property
    def dims(self):
        return self._dims
    @property
    def dim_sizes(self):
        return self._dim_sizes

    def get_dim(self, dim):
        assert dim < self._dims
        return self._dimension[dim]

def storage_classification(pipeline):
    '''
    Classifies the compute objects into separate groups based on their storage
    sizes.
    '''
    comps = pipeline._comp_objs
    children = pipeline._comp_objs_children

    def compute_sizes(comps):
        '''
        For each dimension of the compute object, find the interval size and
        compute the total size of the compute object
        '''
        # dict 'comp_size' : {comp -> (interval_sizes, total_size)}
        # list 'interval_sizes' : [ interval_size, ]
        # tuple 'interval_size' : (param, size_expr)

        comp_size = {}
        for comp in comps:
            interval_sizes = []
            total_size = 1
            intervals = comp.domain
            dim = 0
            for interval in intervals:
                params = interval.collect(Parameter)
                assert not len(params) > 1
                if len(params) == 1:
                    param = params[0]
                elif len(params) == 0:
                    param = 0
                size = interval.upperBound - interval.lowerBound
                interval_sizes.append((param, size))
                total_size = total_size * size
                dim += 1
            # ***
            log_level = logging.DEBUG
            LOG(log_level, str(total_size))
            # ***
    
            comp_size[comp] = (interval_sizes, total_size)

        return comp_size

    def create_storage_class(comps, comp_size):
        '''
        Create classes of storage with different sizes. Default storage class
        is 'const' class which includes arrays with no Parameters. Coefficients
        of dimension parameters are used as keys to map 
        '''
        # dict 'storage' : {comp -> Storage}

        storage = {}
        for comp in comps:
            dim_sizes = comp_size[0]
            total_size = comp_size[1]
            dims = len(dim_sizes)
            storage[comp] = Storage(dims, dim_sizes)

    # compute the total size of the compute object using the interval
    # information
    comp_size = compute_sizes(comps)

    return

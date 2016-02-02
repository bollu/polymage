from __future__ import absolute_import, division, print_function

import logging
from expression import *
from pipe import *

# LOG CONFIG #
pipe_logger = logging.getLogger("storage_mapping.py")
pipe_logger.setLevel(logging.INFO)
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

class Dimension:
    def __init__(self, size_map):

        _param = size_map[0]
        if _param == 0:  # constant
            self._param = '0'
        else:  # Parameter
            self._param = _param.name

        self._size_expr = size_map[1]

        coeff_map = get_affine_var_and_param_coeff(self._size_expr)
        self._const = get_constant_from_expr(self._size_expr)
        self._coeff = 1
        if not self.is_constant:
            self._coeff = coeff_map[_param]
        else:
            self._coeff = self._const

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
    @property
    def is_constant(self):
        return self._param == '0'

class Storage:
    def __init__(self, _dims, _dim_sizes):
        self._dims = _dims
        self._dim_sizes = _dim_sizes
        self._dimension = []
        for dim in range(0, self._dims):
            self._dimension.append(Dimension(self._dim_sizes[dim]))
        self._lookup_key = self.generate_key()

    @property
    def dims(self):
        return self._dims
    @property
    def dim_sizes(self):
        return self._dim_sizes
    @property
    def lookup_key(self):
        return self._lookup_key

    def get_dim(self, dim):
        assert dim < self._dims
        return self._dimension[dim]

    def generate_key(self):
        key = [self.dims]

        # get (param, coeff) key from each dim
        param_keys = []
        for dim in range(0, dims):
            storage_dim = storage.get_dim(dim)
            param_keys.append((storage_dim.param, storage_dim.coeff))
        param_keys = sorted(param_keys, key=lambda x:x[0])

        key.extend(param_keys)
        # convert to string because list as a dict key is not allowed
        key = str(key)

        return key

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

    def create_storage_object(comps, comp_size):
        '''
        Create classes of storage with different sizes. Default storage class
        is 'const' class which includes arrays with no Parameters. Coefficients
        of dimension parameters are used as keys to map 
        '''
        # dict 'storage' : {comp -> Storage}

        storage_map = {}
        for comp in comps:
            dim_sizes = comp_size[comp][0]
            total_size = comp_size[comp][1]
            dims = len(dim_sizes)
            storage_map[comp] = Storage(dims, dim_sizes)

        return storage_map

    def find_storage_equivalence(comps, storage_map):
        '''
        Create a mapping to the compute object from it's size properties.
        To create mapping, first we need to generate keys this way-
        - Field 0 : dimensionality 'dim' of the compute object
        - The following 'dim' fields are tuples of Parameter names with their
          respective coefficents. The fields are sorted using the parameter
          names.
        '''
        storage_class_map = {}
        key_map = {}
        for comp in comps:
            storage = storage_map[comp]
            key_map[comp] = storage.lookup_key
            if key not in storage_class_map:
                storage_class_map[key] = [comp]
            else:
                storage_class_map[key].append(comp)

        return key_map, storage_class_map

    # compute the total size of the compute object using the interval
    # information
    comp_size = compute_sizes(comps)

    # create storage classes
    storage_map = create_storage_object(comps, comp_size)

    # find equivalence in size between storage objects and create classes of
    # storage objects
    key_map, storage_class_map = find_storage_equivalence(comps, storage_map)

    # ***
    for key in storage_class_map:
        log_level = logging.DEBUG-1
        comp_names = [comp.name for comp in storage_class_map[key]]
        LOG(log_level, comp_names)
    # ***

    return

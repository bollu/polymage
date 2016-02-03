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
        self._orig_param = _param

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
    def orig_param(self):
        return self._orig_param
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
        self._offsets = self.gen_param_offsets()

    @property
    def dims(self):
        return self._dims
    @property
    def dim_sizes(self):
        return self._dim_sizes
    @property
    def lookup_key(self):
        return self._lookup_key
    @property
    def offsets(self):
        return self._offsets

    def get_dim(self, dim):
        assert dim < self._dims
        return self._dimension[dim]

    def generate_key(self):
        '''
        To create class mapping, we generate keys this way -
        - Field 0 : dimensionality 'dim' of the compute object
        - Following 'dim' fields are tuples of Parameter names with their
          respective coefficents. The fields are sorted using the parameter
          names.
        '''

        key = [self.dims]

        # get (param, coeff) key from each dim
        param_keys = []
        for dim in range(0, self.dims):
            storage_dim = self.get_dim(dim)
            param_keys.append((storage_dim.param, storage_dim.coeff))
        param_keys = sorted(param_keys, key=lambda x:x[0])

        key.extend(param_keys)
        # convert to string because list as a dict key is not allowed
        key = str(key)

        return key

    def gen_param_offsets(self):
        # get (param, const) from each dim
        param_offsets = []
        for dim in range(0, self.dims):
            storage_dim = self.get_dim(dim)
            param_offsets.append((storage_dim.param, storage_dim.const))
        param_offsets = sorted(param_offsets, key=lambda x:x[0])

        return param_offsets

    def compute_total_size(self):
        total_size = 1
        for size in self.dim_sizes:
            total_size *= size

        return total_size

def storage_classification(comps):
    '''
    Classifies the compute objects into separate groups based on their storage
    sizes.
    '''
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
            intervals = comp.domain
            dim = 0
            for interval in intervals:
                params = interval.collect(Parameter)
                assert not len(params) > 1
                if len(params) == 1:
                    param = params[0]
                elif len(params) == 0:  # const
                    param = 0
                size = interval.upperBound - interval.lowerBound
                interval_sizes.append((param, size))
                dim += 1
    
            comp_size[comp] = interval_sizes

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
            dim_sizes = comp_size[comp]
            dims = len(dim_sizes)
            storage_map[comp] = Storage(dims, dim_sizes)

        return storage_map

    def find_storage_equivalence(comps, storage_map):
        '''
        Create a mapping to the compute object from it's size properties.
        The classification can be further improved with the knowledge of param
        constraints or estimates, by not differentiating b/w dimensions of
        equal sizes.
        '''
        storage_class_map = {}
        key_map = {}
        for comp in comps:
            storage = storage_map[comp]
            key = storage.lookup_key
            key_map[comp] = key
            if key not in storage_class_map:
                storage_class_map[key] = [comp]
            else:
                storage_class_map[key].append(comp)

        return key_map, storage_class_map

    def maximal_storage(comps, storage_class_map, storage_map):
        '''
        Compute the maximal storage needed at each dimension individually and
        over approximate the total storage to be the product of maximal storage
        of all dimensions. This can further be improved with the knowledge of
        param constraints or estimates which suggests an exact measure of the
        size of each dimension.
        '''
        new_storage_map = {}
        for key in storage_class_map:
            storage_class = storage_class_map[key]  # a list
            # pick a dummpy comp to get the total number of dimensions and the
            # original parameter associated with each dimesnion
            helper_comp = storage_class[0]
            helper_storage = storage_map[helper_comp]
            dims = helper_storage.dims
            # this list holds the maximal offset value for each dimension
            max_offset = [0 for dim in range(0, dims)]
            for comp in storage_class:
                storage = storage_map[comp]
                offsets = storage.offsets
                for dim in range(0, dims):
                    dim_off = offsets[dim][1]  # its a tuple
                    max_offset[dim] = int(max(max_offset[dim], dim_off))

            # collect the dim storage info and update with the new maximal
            # offset
            dim_sizes = []
            for dim in range(0, dims):
                dim_storage = storage.get_dim(dim)
                param = dim_storage.orig_param
                size = (Fraction(dim_storage.coeff)) * dim_storage.orig_param + \
                       max_offset[dim]
                dim_sizes.append((param, size))

            # final maximal storage for this class
            max_storage = \
                Storage(dims, dim_sizes)

            # all comps of this class now have identical storage
            for comp in storage_class:
                new_storage_map[comp] = max_storage

        # clear the temporary mappings
        storage_class_map.clear()

        return new_storage_map

    # compute the total size of the compute object using the interval
    # information
    comp_size = compute_sizes(comps)

    # create storage classes
    storage_map = create_storage_object(comps, comp_size)

    # find equivalence in size between storage objects and create classes of
    # storage objects
    key_map, storage_class_map = find_storage_equivalence(comps, storage_map)

    # compute the maximal offsets in each dimension of the compute objects,
    # and compute the total_size of the storage for each storage class
    storage_map = maximal_storage(comps, storage_class_map, storage_map)

    return

# TESTME
def classify_scratchpad_storage(comps, size_map):
    pass

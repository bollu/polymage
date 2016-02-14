from __future__ import absolute_import, division, print_function

import logging
from expression import *
from pipe import *

# LOG CONFIG #
storage_logger = logging.getLogger("storage_mapping.py")
storage_logger.setLevel(logging.INFO)
LOG = storage_logger.log

class TypeSizeMap(object):
    _type_size_map = { "void":1,
                       "int8":1, "uint8":1,
                       "int16":2, "uint16":2,
                       "int32":4, "uint32":4,
                       "int64":8, "uint64":8,
                       "float":4, "double":8 }

    @classmethod
    def getsize(cls, typ):
        typ_name = typ.c_type_name()
        assert typ_name in cls._type_size_map
        return cls._type_size_map[typ_name]

def compute_liveness(children_map, schedule):
    '''
    Given a schedule for a DAG of compute objects, this function computes the
    liveness range of each compute object. The output is a mapping from the
    timestamp in the schedule to a list of compute objects that are live only
    upto the corresponding time.
    '''
    liveness_map = {}
    for comp in schedule:
        last_live = 0
        for child in comp.children:
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
    def __init__(self, _typ, _dims, _dim_sizes):
        self._typ = _typ
        self._dims = _dims
        self._dim_sizes = _dim_sizes

        self._dimension = []
        for dim in range(0, self._dims):
            self._dimension.append(Dimension(self._dim_sizes[dim]))

        self._lookup_key = self.generate_key()
        self._offsets = self.gen_param_offsets()

    @property
    def typ(self):
        return self._typ
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
        - Field 0 : size, in bytes, of data type of the compute object
        - Field 1 : dimensionality 'dim' of the compute object
        - Following 'dim' fields are tuples of Parameter names with their
          respective coefficents. The fields are sorted using the parameter
          names.
        '''

        key = [TypeSizeMap.getsize(self.typ), self.dims]

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

    def find_storage_equivalence(comps):
        '''
        Create a mapping to the compute object from it's size properties.
        The classification can be further improved with the knowledge of param
        constraints or estimates, by not differentiating b/w dimensions of
        equal sizes.
        NOTE: This module is unaware of whether the pipeline outputs must be
        excluded from classification with other compute objects.
        '''
        storage_class_map = {}
        for comp in comps:
            storage = comp.orig_storage_class
            key = storage.lookup_key
            if key not in storage_class_map:
                storage_class_map[key] = [comp]
            else:
                storage_class_map[key].append(comp)

        return storage_class_map

    def maximal_storage(comps, storage_class_map):
        '''
        Compute the maximal storage needed at each dimension individually and
        over approximate the total storage to be the product of maximal storage
        of all dimensions. This can further be improved with the knowledge of
        param constraints (or estimates) which suggests an exact (or
        approximate) measure of the size of each dimension.
        '''
        # ***
        log_level = logging.DEBUG
        LOG(log_level, "Storage Classes:")
        # ***
        new_storage_class_map = {}
        for key in storage_class_map:
            class_comps = storage_class_map[key]  # a list

            # ***
            log_level = logging.DEBUG
            LOG(log_level, "_______")
            LOG(log_level, key)
            LOG(log_level, [comp.func.name for comp in class_comps])
            # ***

            # pick a dummpy comp to get the total number of dimensions and the
            # original parameter associated with each dimesnion
            helper_comp = class_comps[0]
            typ = helper_comp.func.typ
            dims = helper_comp.func.ndims
            helper_storage = helper_comp.orig_storage_class

            # this list holds the maximal offset value for each dimension
            max_offset = [0 for dim in range(0, dims)]
            for comp in class_comps:
                storage = comp.orig_storage_class
                offsets = storage.offsets
                for dim in range(0, dims):
                    dim_off = offsets[dim][1]  # its a tuple
                    max_offset[dim] = int(max(max_offset[dim], dim_off))

            # collect the dim storage info and update with the new maximal
            # offset
            dim_sizes = []
            for dim in range(0, dims):
                dim_storage = helper_storage.get_dim(dim)
                param = dim_storage.orig_param
                coeff = dim_storage.coeff
                size = Fraction(coeff) * param + max_offset[dim]
                dim_sizes.append((param, size))

            # final maximal storage for this class
            max_storage = Storage(typ, dims, dim_sizes)

            # all comps of this class now have identical storage
            for comp in class_comps:
                comp.set_storage_class(max_storage)
                key = max_storage.lookup_key
                new_storage_class_map[key] = comp

        # clear the temporary mappings
        storage_class_map.clear()

        return new_storage_class_map

    # find equivalence in size between storage objects and create classes of
    # storage objects
    storage_class_map = find_storage_equivalence(comps)

    # compute the maximal offsets in each dimension of the compute objects,
    # and compute the total_size of the storage for each storage class
    storage_class_map = maximal_storage(comps, storage_class_map)

    return storage_class_map

def allocate_physical_arrays(pipeline):
    '''
    Generate a mapping from logical storage object of the comp (assumed to be
    available at this point), to cgen CArrays. The mapping can be switched
    between naive and optimized (with reuse) versions, given a schedule for
    the comps within its group.
    '''

    return

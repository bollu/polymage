from __future__ import absolute_import, division, print_function

import logging
import targetc as genc
from expression import *
from pipe import *

# LOG CONFIG #
storage_logger = logging.getLogger("storage_mapping.py")
storage_logger.setLevel(logging.DEBUG)
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

        self._id = None

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
    def id_(self):
        return self._id
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

    def generate_id(self):
        self._id = IdGen.get_stg_id()

def classify_storage(comps):
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
            max_storage.generate_id()

            # all comps of this class now have identical storage
            for comp in class_comps:
                comp.set_storage_class(max_storage)
                new_storage_class_map[max_storage] = comp

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


def log_schedule(comps, schedule):
    log_level = logging.DEBUG
    LOG(log_level, "\n_______")
    LOG(log_level, "Schedules :")
    for comp in comps:
        LOG(log_level, comp.func.name+" : "+str(schedule[comp]))
    return

def log_storage_mapping(comps, storage_map):
    log_level = logging.DEBUG
    LOG(log_level, "\n_______")
    LOG(log_level, "Storage Mapping :")
    for comp in comps:
        LOG(log_level, comp.func.name+" : "+str(storage_map[comp]))
    return

def remap_storage_for_comps(storage_class_map, schedule,
                            liveness_map, storage_map):

    # sort comps according to their schedule
    sorted_comps = get_sorted_objs(schedule)

    # initialize a pool of arrays for each storage class
    array_pool = {}
    array_count = 0
    for stg_class in storage_class_map:
        array_pool[stg_class] = []

    for comp in sorted_comps:
        stg_class = comp.storage_class
        # if no array of stg_class is free as of now
        if not array_pool[stg_class]:
            array_count += 1
            storage_map[comp] = array_count
        # there is a free array of stg_class in the pool
        else:
            storage_map[comp] = array_pool[stg_class].pop()

        # return free arrays to pool
        time = schedule[comp]
        # if any comp is not live after this point
        if time in liveness_map:
            free_comps = liveness_map[time]
            for free_comp in free_comps:
                comp_stg_class = free_comp.storage_class
                storage_index = storage_map[free_comp]
                array_pool[comp_stg_class].append(storage_index)

    # ***
    log_schedule(sorted_comps, schedule)
    log_storage_mapping(sorted_comps, storage_map)
    # ***

    return

def remap_storage_for_group(group, storage_class_map, storage_map):

    # compute liveness
    # 1. prepare children map for liveness computation
    children_map = {}
    for comp in group.comps:
        children_map[comp] = \
            [child for child in comp.children \
                     if child.group == group]
    # 2. get schedule for compute objects
    comps_schedule = group.comps_schedule

    liveness_map = compute_liveness(children_map, comps_schedule)

    remap_storage_for_comps(storage_class_map, comps_schedule,
                            liveness_map, storage_map)

    return

def remap_storage_for_liveout(pipeline, storage_class_map, storage_map):
    '''
    Separate out the liveout compute objects in the pipeline, create a
    temporary graph using a children_map, compute liveness_map for this graph,
    and remap storage objects for liveout comps.
    '''
    liveouts = [comp for comp in pipeline.comps \
                       if comp.is_liveout]
    grp_schedule = pipeline.group_schedule
    comps_schedule = {}
    children_map = {}
    for comp in liveouts:
        # schedule of the group liveouts is the group schedule itself
        comps_schedule[comp] = grp_schedule[comp.group]

        # find temporary children map involving only the compute objects which
        # are not scratchpads
        # collect groups where comp is livein
        g_liveouts = []
        if comp.children:
            livein_groups = [child.group for child in comp.children]
            # collect liveouts of these groups
            for g in livein_groups:
                g_liveouts += g.liveouts
        children_map[comp] = g_liveouts

    liveness_map = compute_liveness(children_map, comps_schedule)

    remap_storage_for_comps(storage_class_map, comps_schedule,
                            liveness_map, storage_map)

    return

def remap_storage(pipeline):
    '''
    Map logical storage objects to representative physical arrays
    The mapping can be switched between naive and optimized (with reuse)
    versions, given a schedule for the comps within its group.
    '''
    storage_class_map = pipeline.storage_class_map
    # a mapping from comp -> index of array of comp's storage class
    storage_map = {}

    for group in pipeline.groups:
        remap_storage_for_group(group, storage_class_map, storage_map)

    remap_storage_for_liveout(pipeline, storage_class_map, storage_map)

    return storage_map

def create_physical_arrays(pipeline):
    '''
    Create cgen CArrays for compute objects using the storage mapping from
    logical storage object of the comp (assumed to be available at this point).
    '''
    flat_scratch = 'flatten_scratchpad' in pipeline.options
    for group in pipeline.groups:
        # place where created arrays are recorded
        created = {}
        for comp in group.comps:
            stg_class = comp.storage_class
            created[stg_class] = {}

        # create / map CArray objects to comps
        for comp in group.comps:
            stg_class = comp.storage_class
            array_id = pipeline.storage_map[comp]
            if array_id in created[stg_class]:
                array = created[stg_class][array_id]
            else:
                # array attributes
                array_type = genc.TypeMap.convert(comp.func.typ)
                tag = str(stg_class.id_)
                array_name = genc.CNameGen.get_array_name(tag)
                array_sizes = stg_class.dim_sizes
                # create CArray object
                array = genc.CArray(array_type, array_name, array_sizes)
                if comp.is_liveout:  # full array
                    array.layout = 'contiguous'
                else:  # scratchpad
                    if flat_scratch:  # linearized array
                        array.layout = 'contiguous_static'
                # record the array creation
                created[stg_class][array_id] = array
    return array


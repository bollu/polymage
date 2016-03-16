from __future__ import absolute_import, division, print_function

import pipe
import poly
import logging

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s")

def buildPipeline(outputs,
                  param_estimates = [],
                  param_constraints = [],
                  grouping = [],
                  group_size = None,
                  inline_directives = [],
                  tile_sizes = [],
                  size_threshold = None,
                  pipe_name = None,
                  options = []):

    # Create an isl context that will be used for all polyhedral
    # operations during compilation.
    ctx = poly.isl.Context()

    if group_size == None:
        group_size = 5

    if tile_sizes == []:
        tile_sizes = [16, 16, 16]

    if size_threshold == None:
        size_threshold = 200*200

    #options.append('flatten_scratchpad')
    #options.append('optimize_storage')
    #options.append('early_free')
    #options.append('pool_alloc')

    if 'optimize_storage' in options:
        options.append('flatten_scratchpad')

    options = list(set(options))

    return pipe.Pipeline(_ctx = ctx,
                         _outputs = outputs,
                         _param_estimates = param_estimates,
                         _param_constraints = param_constraints,
                         _grouping = grouping,
                         _group_size = group_size,
                         _inline_directives = inline_directives,
                         _tile_sizes = tile_sizes,
                         _size_threshold = size_threshold,
                         _name = pipe_name,
                         _options = options)

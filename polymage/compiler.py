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
# compiler.py : Entry-level of PolyMage compiler. Output function of the
#               expressed pipeline, program parameter information and
#               compilation flags are passed here.
#

from __future__ import absolute_import, division, print_function

from . import pipe
from . import poly
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

    outputs = list(set(outputs))

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
    #options.append('multipar')

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

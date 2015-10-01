from __future__ import absolute_import, division, print_function

import pipe
import poly
import logging

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s")

def buildPipeline(outputs,
                  param_estimates = [],
                  param_constraints = [],
                  grouping = [],
                  pipe_name = None,
                  options = []):

    # Create an isl context that will be used for all polyhedral
    # operations during compilation.
    ctx = poly.isl.Context()

    return pipe.Pipeline(_ctx = ctx,
                         _outputs = outputs,
                         _param_estimates = param_estimates,
                         _param_constraints = param_constraints,
                         _grouping = grouping,
                         _name = pipe_name,
                         _options = options)

from __future__ import absolute_import, division, print_function

import pipe
import poly
import logging

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s")

def buildPipeline(outputs, paramConstraints = [], grouping = [], pipe_name = None, options = []):
    # Create an isl context that will be used for all polyhedral
    # operations during compilation.
    ctx = poly.isl.Context()

    return pipe.Pipeline(ctx, outputs, paramConstraints, grouping, _name=pipe_name, _options = options)

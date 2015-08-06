from __future__ import absolute_import, division, print_function

import pipe
import poly
import logging

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s")

def buildPipeline(outputs, pipe_name = None, paramConstraints = [], grouping = []):
    # Create an isl context that will be used for all polyhedral
    # operations during compilation.
    ctx = poly.isl.Context()

    return pipe.Pipeline(ctx, outputs, paramConstraints, grouping, pipe_name)

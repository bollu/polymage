# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

import sys
sys.path.insert(0, '../frontend')
sys.path.insert(0, '../codegen')

from IR import *
import Poly as opt

def buildPipeline(outputs, pipeName = None, paramConstraints = [], grouping = []):
    # Create an isl context that will be used for all polyhedral
    # operations during compilation.
    ctx = opt.isl.Context()

    return Pipeline(ctx, outputs, paramConstraints, grouping, pipeName)

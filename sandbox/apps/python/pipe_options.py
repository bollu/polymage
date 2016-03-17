import optparse

parser = optparse.OptionParser()

# Storage Optimization takes care of reusing arrays associated with different
# Functions. Enabling this optimization leads to reuse of arrays after a point
# in the program, where the function is no more used.
help_str = 'True : Optimize for storage, by reusing arrays if possible, \
            False: Create separate array for each Function'
parser.add_option('--optimize_storage',
                  action='store_true',
                  dest='optimize_storage',
                  default=True,
                  help=help_str)

# Flattening of scratchpads to a single dimension
# This is helpful for storage optimization, if enabled, since scratchpads are
# of constant size, and sizes across individual dimensions do not matter once
# the array is linearized.
help_str = 'True : Linearize / flatten const scratchpads to one dimension, \
            False: Create multi-dimensional scratchpads'
parser.add_option('--flatten_scratchpad',
                  action='store_true',
                  dest='flatten_scratchpad',
                  default=True,
                  help=help_str)

# Early freeing of the arrays results in lesser memory footprint, by throwing
# away all unwanted arrays from the memory as soon as possible.
help_str =  'True : Free the arrays not in use as soon as possible, \
             False: Procrastinated freeing of arrays'
parser.add_option('--early_free',
                  action='store_true',
                  dest='early_free',
                  default=True,
                  help=help_str)

# Use a pool of memory, so that allocation from pool can return the pointer to
# an already allocated array, thus resulting in avoiding frequent malloc()
# calls across allocations / iterations of program
parser.add_option('--pool_alloc',
                  action='store_true',
                  dest='pool_alloc',
                  default=True,
                  help='True : Use a pool of memory allocations, \
                        False: generate simple malloc function call')

import optparse
import sys

def parse_args():
    parser = optparse.OptionParser()

    help_str = \
        '"new" : from scratch | "existing" : compile and run |  "ready" : just run'
    parser.add_option('-m', '--mode',
                       type='choice',
                       action='store',
                       dest='mode',
                       choices=['new', 'existing', 'ready', 'tune'],
                       default=['new'],
                       help=help_str)
  
    parser.add_option('-c', '--cycle',
                       action='store',
                       dest='cycle',
                       help='W or V cycle')

    parser.add_option('-a', '--nit',
                       action='store',
                       dest='nit',
                       help='nit')

    parser.add_option('-n', '--n',
                       action='store',
                       dest='n',
                       default=255,
                       help='')

    parser.add_option('-l', '--L',
                       action='store',
                       dest='L',
                       default=2,
                       help='')

    parser.add_option('-e', '--nu1',
                       action='store',
                       dest='nu1',
                       default=10,
                       help='')

    parser.add_option('-f', '--nuc',
                       action='store',
                       dest='nuc',
                       default=0,
                       help='')

    parser.add_option('-g', '--nu2',
                       action='store',
                       dest='nu2',
                       default=0,
                       help='')

    parser.add_option('-p', '--problem',
                       action='store',
                       dest='problem',
                       default=1,
                       help='problem')

    parser.add_option('-r', '--runs',
                       action='store',
                       dest='runs',
                       default=1,
                       help='number of runs',)

    parser.add_option('-t', '--timer',
                       action='store_true',
                       dest='timer',
                       default=False,
                       help='True : report execution time, \
                             False: do not collect timing info',)

    parser.add_option('--pool_alloc',
                       action='store_true',
                       dest='pool_alloc',
                       default=False,
                       help='True : Use a pool of memory allocations, \
                             False: generate simple malloc function call',)

    (options, args) = parser.parse_args()

    return options

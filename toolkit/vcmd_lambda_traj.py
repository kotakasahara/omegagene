#!/usr/bin/python2.6

from optparse import OptionParser
import numpy
import sys
import kkpdb
import kkatomgroup
import kkmmconfig 

def option_parser():
    p = OptionParser()
    
    p.add_option('--cfg', dest='fn_config',
                 type="append",
                 help="Omegagene configuration files")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"

    return opts, args

def _main():
    opts, args = option_parse()
    


if __name__ == "__main__":
    _main()



#!/usr/bin/python2.6

from optparse import OptionParser
import numpy
import sys
import kkpdb
import kkatomgroup
import kkmmconfig 

def option_parser():
    p = OptionParser()
    
    p.add_option('--i-cod', dest='fn_cod',
                 help="Omegagene coordinate file")
    p.add_option('--i-lmb', dest='fn_lmb',
                 help="Omegagene lambda file")
    p.add_option('--i-vs', dest='fn_vs',
                 help="Omegagene vs file")
    p.add_option('--interval-vs', dest='itv_vs',
                 help="VS transition interval")
    p.add_option('--interval-lmb', dest='itv_lmb',
                 help="Lambda output interval")
    p.add_option('--interval-cod', dest='itv_cod',
                 help="Coordinate output interval")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"

    return opts, args

def _main():
    opts, args = option_parse()
    


if __name__ == "__main__":
    _main()



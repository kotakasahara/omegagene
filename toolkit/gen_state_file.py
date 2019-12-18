#!/usr/bin/python2.7

from optparse import OptionParser
import kkpdb 

def get_options():
    p = OptionParser()
    p.add_option('-i', dest='fn_in',
                 help="input file")
    p.add_option('-o', dest='fn_out',
                 help="output file")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args


def read_state(fn):
    f = open(fn)
    state = []
    for line in f:
        state = [int(x) for x in line.strip().split()]
    f.close()
    return state

def write_state(fn, state):
    f = open(fn, "w")
    f.write(str(len(state)) + "\n")
    for s in state:
        f.write(str(s) + "\n")
    f.close()
    return

def _main():
    opts, args = get_options()
    state = read_state(opts.fn_in)
    write_state(opts.fn_out, state)
    return


if __name__ == "__main__":
    _main()


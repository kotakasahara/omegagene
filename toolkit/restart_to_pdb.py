#!/usr/bin/python2.7


import kkpresto_restart as rest
import kkstruct as kst
import kkpdb 
import numpy

from optparse import OptionParser

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_restart',
                 help="file name for restart file")
    p.add_option('--i-pdb', dest='fn_i_pdb',
                 help="file name for pdb file")
    p.add_option('-o', dest='fn_out',
                 help="file name for pdb file")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def _main():
    opts, args = get_options()
    restart = rest.PrestoRestartReader(opts.fn_restart).read_restart()    
    template = kkpdb.PDBReader(opts.fn_i_pdb).read_model()
    
    for i, at in enumerate(template.atoms):
        at.crd = numpy.array(restart.crd[i])
    
    kkpdb.PDBWriter(opts.fn_out).write_model(template)

if __name__ == "__main__":
    _main()

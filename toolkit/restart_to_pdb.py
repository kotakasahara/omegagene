#!/usr/bin/python2.7


import kkpresto_restart as rest
import kkstruct as kst
import kkpdb 
import numpy

from optparse import OptionParser

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_restart',
                 action="append",
                 default=[],
                 help="file name for restart file")
    p.add_option('--i-list', dest='fn_list_restart',
                 help="file name for a list of restart files")
    p.add_option('--i-pdb', dest='fn_i_pdb',
                 help="file name for pdb file")
    #p.add_option('-d', "--del", dest='delete_res',
    #             action="append",
    #             default = [],
    #             help="residue names to be deleted")
    p.add_option('-o', dest='fn_out',
                 help="file name for pdb file")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def read_list(fn):
    flist = []
    f = open(fn)
    for line in f:
        flist.append(line.strip().split()[0])
    return flist

def _main():
    opts, args = get_options()

    template = kkpdb.PDBReader(opts.fn_i_pdb).read_model()

    restart_files = opts.fn_restart
    if opts.fn_list_restart:
        restart_files.extend(read_list(opts.fn_list_restart))
    

    writer = kkpdb.PDBWriter(opts.fn_out)
    writer.open()
    for i_f, fn in enumerate(restart_files):
        restart = rest.PrestoRestartReader(fn).read_restart()        
        for i, at in enumerate(template.atoms):
            at.crd = numpy.array(restart.crd[i])
        template.model_id = i_f+1
        writer.add_model(template)
        #kkpdb.PDBWriter(opts.fn_out).write_model(template)
    writer.close()

if __name__ == "__main__":
    _main()

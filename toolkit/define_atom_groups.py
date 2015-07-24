#!/usr/bin/python2.7

from MDAnalysis import Universe
from MDAnalysis.coordinates.PDB import PrimitivePDBWriter  
from optparse import OptionParser
import sys
import kkatomgroup

def get_options():
    p = OptionParser()
    p.add_option('-p', '--pdb', dest='fn_pdb',
                 help="File name for pdb file")
    p.add_option('-o', '--out', dest='fn_out',
                 help="File name for output, atom group definition file")
    p.add_option('-s', '--select', dest='select',
                 action="append",
                 help="Atom selection strings in MDAnalysis syntex")
    p.add_option('--out-pdb', dest='fn_out_pdb',
                 action="append",
                 help="File name for output pdb files including only selected atoms")
    p.add_option('-n', '--name', dest='select_name',
                 action="append",
                 help="Names for each selection")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def _main():
    opts, args = get_options()
    if len(opts.select) != len(opts.select_name):
        print "The number of selection strings (-s, --select) and the number of names (-n or --name) should be the same."
        sys.exit(1)
    u = None

    try:
        u = Universe(opts.fn_pdb)
    except:
        print "File read error: " + opts.fn_pdb
        sys.exit(1)

    f_o = None    
    try:
        f_o = open(opts.fn_out, "w")
    except:
        print "File open error: " + opts.fn_out
        sys.exit(1)
        
    for i, sel_str in enumerate(opts.select):
        sel_name = opts.select_name
        sel_atoms = u.selectAtoms(sel_str)
        f_o.write(sel_name[i] +  " " + " ".join([str(at.number) for at in sel_atoms]) + "\n")
        if len(opts.fn_out_pdb) > i:
            writer = PrimitivePDBWriter(opts.fn_out_pdb[i])
            writer.write(sel_atoms)
    f_o.close()


if __name__ == "__main__":
    _main()


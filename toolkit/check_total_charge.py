#!/usr/bin/python2.6

from optparse import OptionParser
import numpy
import sys
import kkpdb
import kkpresto

def _main():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_pdb',
                 help="file name for input PDB generated by tplgene")
    p.add_option('--i-tpl', dest='fn_tpl',
                 help="file name for input TPL")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"

    if opts.fn_pdb:
        print "READ PDB: " + opts.fn_pdb
        reader = kkpdb.PDBReader(opts.fn_pdb)
        model = reader.read_model()
        sum_tf = 0.0
        for atom in model.atoms:
            sum_tf += atom.tf
            print str(atom.atom_id) + ": " + str(atom.tf)
        print sum_tf

    if opts.fn_tpl:
        print "READ TPL: " + opts.fn_tpl
        reader = kkpresto.TPLReader(opts.fn_tpl)
        tpl = reader.read_tpl()

        sum_charge = 0.0
        for tplmol in tpl.mols:
            molname = tplmol.mol_name
            sum_charge_mol = 0.0
            for tplatom in tplmol.atoms:
                sum_charge_mol += tplatom.charge
            print "MOL [" + tplmol.mol_name + "] :" + str(sum_charge_mol) + " * " + str(tplmol.mol_num)
            sum_charge += sum_charge_mol * tplmol.mol_num
            print "N_AOTMS: "+ str(len(tplmol.atoms))

        print sum_charge


if __name__ == "__main__":
    _main()

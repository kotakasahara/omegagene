#!/usr/bin/python2.7

import kkpdb 

def get_options():
    p = OptionParser()
    p.add_option('--pdb', dest='fn_pdb',
                 # action="append",
                 help="PDB file")
    p.add_option('--heavy', dest='flg_heavy',
                 action="store_true",
                 help="only heavy atoms")

    p.add_option("--type", dest="group_type",
                 type="choice",
                 choices = ["res"],
                 default="res",
                 help = "Unit of each group")
    
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def _main():
    opts, args = get_options()
    self.structure = kkpdb.PDBReader(self.config.get_val("fn-i-initial-pdb")).read_model()
    
    return


if __name__ == "__main__":
    _main()


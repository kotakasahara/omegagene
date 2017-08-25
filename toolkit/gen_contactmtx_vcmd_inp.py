#!/usr/bin/python2.7

from optparse import OptionParser
import kkatomgroup as kkag
import numpy as np

def get_options():
    p = OptionParser()
    p.add_option('--pref-group', dest='pref_group',
                 action="append",
                 help="prefix for the group names")
    p.add_option('--interval', dest='interval',
                 type="int",
                 help="interval steps for the state transition")
    p.add_option('-o', '--out', dest='fn_out',
                 help="File name for output")
    p.add_option('--atom-groups', dest='fn_atom_groups',
                 help="File name for atom group definition")
    p.add_option('--res-neighbor', dest='res_neighbor',
                 type="int", default=3,
                 help="The maximum threshold of the neighbor residues")
    p.add_option('--default', dest='default_q',
                 type="float",
                 default=1.0,
                 help="default value of the parameter qcano")
    p.add_option('-r','--range', dest='range',
                 action="append",
                 help="The min and max values of lambda in each virtual state, separated with ':'.")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def parse_range(range):
    sep = range.index(":")
    return (int(range[:sep]), int(range[(sep+1):]))

def parse_group(grp_name, pref_group):
    grp_grp = ""
    res_num = -1
    for pref in pref_group:
        idx = -1
        try:
            idx = grp_name.index(pref)
            if idx == 0:
                grp_grp = pref
                res_num = int(grp_name[len(pref):])
                break
        except:
            pass
    return grp_grp, res_num
            


def _main():
    opts, args = get_options()

    grp, grp_names = kkag.AtomGroupsReader(opts.fn_atom_groups).read_groups()
    ranges = [ parse_range(x) for x in opts.range ]

    buf = ""
    n_dim = 0
    for i1, grp_name1 in enumerate(grp_names):
        grp_grp1, res_num1 = parse_group(grp_name1, opts.pref_group)
        if res_num1 < 0: continue
        for i2,grp_name2 in enumerate(grp_names):
            if i1 >= i2: continue
            grp_grp2, res_num2 = parse_group(grp_name2, opts.pref_group)
            if res_num2 < 0: continue
            if abs(res_num1 - res_num2) <= opts.res_neighbor: continue
            
            n_dim += 1
            buf += str(len(ranges))+"\t"+grp_name1+"\t"+grp_name2+"\n"
            for rg in ranges:
                buf += str(rg[0])+"\t"+str(rg[1])+"\n"
    
    print n_dim
    print len(ranges)
    n_states = np.power(len(ranges), n_dim)
    print n_states
    default_q = " ".join([ "0" for x in range(n_dim) ]) + "\t" + str(opts.default_q) + "\n"
    

    fo = open(opts.fn_out, "w")
    fo.write(str(opts.interval)+"\n")
    fo.write(str(n_dim)+"\n")    
    fo.write(buf)
    fo.write(default_q)
    fo.write("END\n")
    fo.close()
    

if __name__ == "__main__":
    _main()

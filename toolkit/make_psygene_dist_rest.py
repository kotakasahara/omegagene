#!/usr/bin/python2.6

from optparse import OptionParser
import re
import kkpdb
import kkpresto
import math

def define_options():
    p = OptionParser()
    
    p.add_option('--i-pdb', dest='fn_pdb',
                 help="pdb file")
    p.add_option('--i-tpl', dest='fn_tpl',
                 help="tpl file")
    p.add_option('-o', dest='fn_out',
                 help="output distance restraint file")
    p.add_option('--chain-nums', dest='molids',
                 action="append",   type = "int",
                 default=[],
                 help="Chain ID to be considered (1, 2, ...) ")
    p.add_option('--res-range', dest='res_range',
                 action="append",   type = "str",
                 default=[],
                 help="range of ResID to be considered, (1-10, ...)")
    p.add_option('--atomnames', dest='atomnames',
                 action="append",
                 default=["CA"],
                 help="atom names")
    p.add_option('--min-sequence-dist', dest='min_seq_dist',
                 type="int",
                 default=4,
                 help="minimum sequence length")
    p.add_option('--max-sequence-dist', dest='max_seq_dist',
                 type="int",
                 default=99999,
                 help="maximum sequence length")
    p.add_option('--max-dist', dest='max_dist',
                 type="float",
                 default=6.5,
                 help="maximum distance")
    p.add_option('--force-coef', dest='force_coef',
                 type="float",
                 default=1.0,
                 help="force coefficient")
    p.add_option('--low-bound', dest='lower_bound', 
                 type="float",
                 default=0.5,
                 help="lower bound for distance")
    p.add_option('--up-bound', dest='upper_bound', 
                 type="float",
                 default=3.0,
                 help="lower bound for distance")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args



def _main():
    opts, args = define_options() 
    print "READ PDB: " + opts.fn_pdb
    model = kkpdb.PDBReader(opts.fn_pdb).read_model()
    
    print "READ TPL: " + opts.fn_tpl
    tpl = kkpresto.TPLReader(opts.fn_tpl).read_tpl()
    tpl.convert_units_to_si()
    mollist = []
    head_atomid = [0]
    prev_head = 0
    for i, mol in enumerate(tpl.mols):
        for j in range(0, mol.mol_num):
            mollist.append(mol)
            head_atomid.append(head_atomid[-1] + len(tpl.mols[i].atoms))
            #if mol.mol_name != "WAT":
                #print "mol " + str(len(mollist)) + " " + str(i) + " " + mol.mol_name
     
    res_range = {}
    for i in range(0, len(opts.molids)):
        if len(opts.res_range) < i+1:
            res_range[opts.molids[i]]=()
        else:
            valid_res = set()
            terms1 = re.compile(":").split(opts.res_range[i])
            for t1 in terms1:
                terms2 = re.compile("-").split(t1)
                for t2 in range(int(terms2[0]), int(terms2[1])+1):
                    valid_res.add(t2)
            print "mol: " + str(i) + " res:" + ",".join([str(x) for x in valid_res])
            res_range[opts.molids[i]] = valid_res
            
    print res_range
    print opts.atomnames

    dist_info = []
    
    for i,molid1 in enumerate(opts.molids):
        mol1 = mollist[molid1]
        for j,molid2 in enumerate(opts.molids):
            mol2 = mollist[molid2]
            if i < j: break
            print (molid1, molid2)
            info = get_dist_info_molmol(
                mol1, mol2, 
                molid1, molid2,
                head_atomid[molid1], head_atomid[molid2],
                model,
                opts.atomnames,
                res_range,
                opts.min_seq_dist,
                opts.max_seq_dist,
                opts.max_dist)
            dist_info.extend(info)
    print "write_dist_info"
    write_dist_info(opts.fn_out, dist_info,
                    opts.lower_bound, opts.upper_bound,
                    opts.force_coef)
    return

def get_dist_info_molmol(mol1, mol2, molid1, molid2,
                         head_atomid1, head_atomid2,
                         model,
                         atomnames, res_range,
                         min_seq_dist, max_seq_dist, 
                         max_dist):
    info = []
    for i, tpl_a1 in enumerate(mol1.atoms):
        atomid1 = head_atomid1 + i
        atom1 = model.atoms[atomid1]
        #print (res_range[molid1][0], res_range[molid1][1], tpl_a1.res_id, tpl_a1.atom_name)

        if not tpl_a1.atom_name in atomnames: continue
        
        if len(res_range[molid1]) > 0 and \
                not tpl_a1.res_id in res_range[molid1]: continue 

        for j, tpl_a2 in enumerate(mol2.atoms):
            atomid2 = head_atomid2 + j
            atom2 = model.atoms[atomid2]
            if not tpl_a2.atom_name in atomnames: continue
            #print (res_range[molid2][0], res_range[molid2][1], tpl_a2.res_id)
            
            if len(res_range[molid2]) > 0 and \
                    not tpl_a2.res_id in res_range[molid2]: continue
            seq_dist = math.fabs(tpl_a1.res_id - tpl_a2.res_id)

            if molid1 == molid2 and (seq_dist < min_seq_dist or seq_dist > max_seq_dist):
                continue 
            
            dist = atom1.distTo(atom2)
            #print (molid1, molid2, tpl_a1.res_id, tpl_a2.res_id, dist)

            if dist >= max_dist: continue
            #print "dd"

            if molid1==molid2 and tpl_a1.res_id >= tpl_a2.res_id: continue
            #print "ee"

            info.append((molid1+1, tpl_a1.res_id, tpl_a1.res_name, tpl_a1.atom_name,
                         molid2+1, tpl_a2.res_id, tpl_a2.res_name, tpl_a2.atom_name,
                         dist))

    return info

def write_dist_info(fn_out, dist_info, low, up, force_coef):
    f_out = open(fn_out,"w")
    f_out.write("RDDSTC> LIST\n")
    for info in dist_info:
        line = "%4d%5d %-4s%4s%6d%5d %-4s%4s"%info[0:8]
        line += "%6.2f%6.2f%9.3f%8.3f%3s"\
            %(force_coef, force_coef, info[8]-low, info[8]+up, "NO" )
        f_out.write(line+"\n")
    f_out.close()
    return

if __name__ == "__main__":
    _main()


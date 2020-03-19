#!/usr/bin/python2.7

import sys
import re
import os
import numpy as np

import kkkit
import kkstruct
import kkpdb
import kkpresto

from optparse import OptionParser

def get_options():
    p = OptionParser()
    
    p.add_option('--pdb', dest='fn_in_pdb',
                 help="filename for the input pdb")
    p.add_option('--param', dest='fn_in_param',
                 help="filename for the input parameter file")
    p.add_option('--tpl', dest='fn_out_tpl',
                 help="filename for the output topology file")
    p.add_option('--molname', dest='mol_name',
                 help="name of the molecule")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

class PotBond(object):
    def __init__(self, res1, atom1, res2, atom2, k, r0):
        self.res1 = res1
        self.atom1 = atom1
        self.res2 = res2
        self.atom2 = atom2
        self.k = k
        self.r0 = r0
        return
        
class PotAngle(object):
    def __init__(self, res1, atom1, res2, atom2, res3, atom3, k, theta0):
        self.res1 = res1
        self.atom1 = atom1
        self.res2 = res2
        self.atom2 = atom2
        self.res3 = res3
        self.atom3 = atom3
        self.k = k
        self.theta0 = theta0
        return
        
class PotTorsion(object):
    def __init__(self, res1, atom1, res2, atom2, res3, atom3,
                 res4, atom4, k, overlap, sym, phase, cal14):
        self.res1 = res1
        self.atom1 = atom1
        self.res2 = res2
        self.atom2 = atom2
        self.res3 = res3
        self.atom3 = atom3
        self.res4 = res4
        self.atom4 = atom4
        self.k = k
        ## k ... energy
        self.overlap = overlap
        ## overlap ... the number of torsions overlapping with this
        self.sym = sym
        ## sym ... symmetry
        self.phase = phase
        ## phase (degree)
        self.cal14 = cal14
        ## flag for calculating 1-4 vdw and coulomb (0 or 1)
        return
        
class PotNonbond(object):
    def __init__(self, res1, atom1, sigma,  epsiron,
                 vdw14, ele14, mass, radius, charge,
                 lmb = 0.0):
        self.res1 = res1
        self.atom1 = atom1
        self.epsiron = epsiron
        self.sigma = sigma
        self.vdw14 = vdw14
        self.ele14 = ele14
        self.mass = mass
        self.radius = radius
        self.charge = charge
        self.lmb = lmb
        return

class ForceField(object):
    def __init__(self):
        self.bond = {}
        ## bond[(res1, atom1, res2, atom2)]
        self.angle = {}
        ## angle[(res1, atom1, res2, atom2, res3, atom3)]
        self.torsion = {}
        ## torsion[(res1, atom1, res2, atom2, res3, atom3, res4, atom4)]
        self.nonbond = {}
        ## nonbond[(res1, atom1)]
        return
    def push_bond(self, res1, atom1, res2, atom2, k, r0):
        self.bond[(res1,atom1,res2,atom2)] = PotBond(res1, atom1, res2, atom2, k, r0)
        return
    def push_angle(self, res1, atom1, res2, atom2, res3, atom3,
                   k, theta0):
        self.angle[(res1,atom1,res2,atom2,res3,atom3)] = PotAngle(res1, atom1, res2, atom2, res3, atom3,
                                                             k, theta0)
        return
    def push_torsion(self, res1, atom1, res2, atom2, res3, atom3,
                     res4, atom4, k, overlap, sym, phase, cal14):
        key = (res1,atom1,res2,atom2,res3,atom3,res4,atom4)
        self.torsion[key] = PotTorsion(res1, atom1, res2, atom2, res3, atom3,
                                res4, atom4, k, overlap, sym, phase, cal14)
        return
    def push_nonbond(self, res1, atom1, sigma, epsiron, vdw14, ele14, mass, radius, charge, lmb = 0.0):
        self.nonbond[(res1,atom1)] = PotNonbond(res1, atom1, 
                                                sigma, epsiron, vdw14, ele14,
                                                mass, radius, charge, lmb)
        return

class ParamReader(kkkit.FileI):
    def __init__(self, fn):
        super(ParamReader, self).__init__(fn)
        return
    def read_param(self):
        ff = ForceField()
        self.open()
        line = self.read_line_com()
        while line:
            terms = line.strip().split()
            if len(terms) < 1:
                line = self.read_line_com()                
                continue
            if terms[0] == "BOND":
                ff.push_bond(terms[1], terms[2], terms[3], terms[4],
                             float(terms[5]), float(terms[6]))
            elif terms[0] == "ANGLE":
                ff.push_angle(terms[1], terms[2], terms[3], terms[4],
                              terms[5], terms[6],
                              float(terms[7]), float(terms[8]))
            elif terms[0] == "TORSION":
                ff.push_torsion(terms[1], terms[2], terms[3], terms[4],
                                terms[5], terms[6], terms[7], terms[8],
                                float(terms[9]), int(terms[10]),
                                int(terms[11]), float(terms[12]), int(terms[13]))
            elif terms[0] == "NONBOND":
                lmb = 0.0
                if len(terms) >= 11:
                    lmb = float(terms[10])
                ff.push_nonbond(terms[1], terms[2], 
                                float(terms[3]), float(terms[4]),
                                float(terms[5]), float(terms[6]),
                                float(terms[7]), float(terms[8]),
                                float(terms[9]), lmb)
            else:
                kkkit.err("Unknown keyword: " + terms[0])
            line = self.read_line_com()
            
        self.close()
        return ff

def set_tpl(model, ff, mol_name):
    tpl = kkpresto.TPL()
    mol_head_atom_id = 0
    
    tplmol = kkpresto.TPLMol(mol_name,
                             1, mol_head_atom_id)
    
    ## bond
    bnd = ff.bond[("*","CA","*","CA")]
    for i in range(len(model.atoms)-1):
        tplmol.add_bond( (i+1, i+2), (bnd.k, bnd.r0) )
        
    ## angle
    ang = ff.angle[("*","CA","*","CA","*","CA")]
    for i in range(len(model.atoms)-2):
        tplmol.add_angle( (i+1, i+2, i+3), (ang.k, ang.theta0*np.pi/180.0) )
    ## torsion
    tor = ff.torsion[("*","CA","*","CA","*","CA","*","CA")]
    for i in range(len(model.atoms)-3):
        
        tplmol.add_torsion( (i+1, i+2, i+3, i+4), (tor.k, tor.overlap,
                                                   tor.sym, tor.phase*np.pi/180.0, tor.cal14) )
    

    atom_type_id_name = {}
    atom_type_name_id = {}
   ## nonbond
    for key, nb  in ff.nonbond.items():
        atom_type = nb.res1+nb.atom1
        if not atom_type in atom_type_name_id:
            atom_type_id = len(atom_type_name_id.keys()) + 1
            atom_type_name_id[atom_type] = atom_type_id
            atom_type_id_name[atom_type_id] = atom_type
        tpl.add_nonbond(atom_type_name_id[atom_type], 
                        [0, 1, nb.sigma, nb.epsiron, nb.vdw14, nb.ele14, nb.lmb, atom_type])


        
    ## atoms
    for i,atom in enumerate(model.atoms):
        atom_type = atom.res_name + atom.atom_name
        atom_type_id = atom_type_name_id[atom_type]
        ffnb = ff.nonbond[(atom.res_name, atom.atom_name)]
        atoms_1_2 = []
        atoms_1_3 = []
        atoms_1_4 = []
        if i < len(model.atoms)-1:
            atoms_1_2.append(1)
        if i < len(model.atoms)-2:
            atoms_1_3.append(2)
        if i < len(model.atoms)-3:
            atoms_1_4.append(3)
        tplmol.add_atom(atom.atom_name,    atom_type,
                        atom_type_id,      atom.res_name,
                        atom.res_id,       ffnb.mass,
                        ffnb.radius,       ffnb.charge,
                        atoms_1_2, atoms_1_3, atoms_1_4,
                        0, 0, 0, 0, 0.0, 0.0, 0.0)

 
    tpl.mols.append(tplmol)
    return tpl

def _main():
    opts, args = get_options()
    model = kkpdb.PDBReader(opts.fn_in_pdb).read_model()
    model.reset_chains_from_atom_info()
    ff = ParamReader(opts.fn_in_param).read_param()

    tpl = set_tpl(model, ff, opts.mol_name)

    f = open(opts.fn_out_tpl,"w")
    f.write(tpl.get_text(header="TPL"))
    f.close()
    
    return
    
if __name__ == "__main__":
    _main()



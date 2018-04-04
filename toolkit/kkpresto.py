#!/usr/bin/python2.7
import re
import kkkit
import numpy as np
import kkdefine as kk
DEBUG = True
##DEBUG = False
J_TO_CAL = 1.0 / 4.1868

class TPLAtom(object):
    def __init__(self, atom_id,
                 atom_name,
                 atom_type,
                 interaction_type,
                 res_name, res_id,
                 mass, radius, charge,
                 atoms_1_2, atoms_1_3, atoms_1_4,
                 z_pos_1_2, z_pos_1_3, z_pos_1_4,
                 z_basis, z_dist, z_ang, z_phase):
        self.n_1_2_atoms = len(atoms_1_2)
        self.n_1_3_atoms = len(atoms_1_3)
        self.n_1_4_atoms = len(atoms_1_4)
        self.atom_id = atom_id
        self.atom_name = atom_name 
        self.atom_type = atom_type
        self.interaction_type = interaction_type
        self.res_name = res_name
        self.res_id = res_id
        self.mass = mass
        self.radius = radius
        self.charge = charge
        self.atoms_1_2 = atoms_1_2
        self.atoms_1_3 = atoms_1_3
        self.atoms_1_4 = atoms_1_4
        self.z_pos_1_2 = z_pos_1_2
        self.z_pos_1_3 = z_pos_1_3
        self.z_pos_1_4 = z_pos_1_4
        self.z_basis = z_basis
        self.z_dist = z_dist
        self.z_ang = z_ang
        self.z_phase = z_phase
        return
    def diff(atom):
        if self.atom_id != atom.atom_id:
            print "atom_id: %5d %5d"%(self.atom_id, atom.atom_id)

class TPLMol():
    def __init__(self, name, num, head_atom_id):
        self.mol_name = name
        self.mol_num = num
        self.head_atom_id = head_atom_id
        self.atoms = []

        ## bonds = ((atom_id1, atom_id2), (energy [kcal/mol], distance[A]))
        self.bonds = []
        ## angles = ((atom_id1, 2, 3), (energy [kcal/mol], angle[degree]))
        self.angles = []
        ## torsions = ((atom_id1, 2, 3, 4),
        ##            (energy [kcal/mol],
        ##             number of overlapping torsions,
        ##             symmetry,
        ##             phase [degree],
        ##             flag for calculating 1-4 interactions)
        self.torsions = []
        self.impropers = []
        self.tor14 = {}
    def diff(mol):
        print "======== MOL ========="
        print "name:    %8s %8s"%(self.mol_name, mol.mol_name)
        print "n_atoms: %8d %8d"%(self.mol_num, mol.mol_num)
        for i, a1 in enumerate(self.atoms):
            a1.diff(mol.atoms[i])
        return 
    def import_from_gromacs_mol(self, gromtopol_mol, gromtopol):
        # set of atoms connecting with the atomid=i via N bonds
        ## cnt2[i] = set(id1, 2, ...)
        cnt2 = {}
        cnt3 = {}
        cnt4 = {}
        cnt4_proper = {}
        for i in range(1,len(gromtopol_mol.atoms)+1):
            cnt2[i] = set()
            cnt3[i] = set()
            cnt4[i] = set()

        for pair, bonds in gromtopol_mol.bonds.items():
            if pair[0] > pair[1]: pair = (pair[1], pair[0])
            cnt2[pair[0]].add(pair[1])
            if len(bonds) != 1:
                sys.stderr.write("Duplicated bond potentials: "+str(pair[0])+" - "+str(pair[1]))
                continue

            #if bond[1][0] == 1:
            #    eps = bond[1][1] * 
            #params = (bond[1][1] * J_CAL, ## kJ/(mol nm^2) => kcal/(mol AA^2)
            #          bond[1][2] * 10.0)  ## nm => A
            #self.add_bond(pair, params)
            bond = bonds[0]
            if bond[0] == 1:    # function type = 1
                print "aaaaa"
                print bond[2]
                params = (bond[2] * J_TO_CAL * 0.01 * 0.5,
                          bond[1] * 10.0)
                self.add_bond(pair, params)
                
        for trio, angles in gromtopol_mol.angles.items():
            if trio[0] > trio[2]: trio = (trio[2], trio[1], trio[0])
            cnt3[trio[0]].add(trio[2])
            if len(angles) != 1:
                sys.stderr.write("Duplicated angle potentials: "+str(pair[0])+" - "+str(pair[1]) + " " + str(pair[2]))
                continue
            for angle in angles:
                if angle[0] == 1:    # function type = 1
                    params = (angle[2] * J_TO_CAL,
                              angle[1] * np.pi / 180.0)
                    self.add_angle(trio, params)

        for quad, dihedrals in gromtopol_mol.dihedrals.items():
            if quad[0] > quad[3]: quad = (quad[3], quad[2], quad[1], quad[0])
            cnt4[quad[0]].add(quad[3])
            ## checking existance of 

            for dihedral in dihedrals:
                if dihedral[0] == 4:    # function type = 4 improper
                    params = (dihedral[2] * J_TO_CAL,
                              1,
                              dihedral[3],
                              dihedral[1] * np.pi / 180.0,
                              0)
                    self.add_improper(quad, params)
                elif dihedral[0] == 9:    # function type = 9 torsion
                    flg_nb14 = 1
                    if quad[0] in cnt4_proper and quad[3] in cnt4_proper[quad[0]]:
                        flg_nb14 = 0
                    params = (dihedral[2] * J_TO_CAL,
                              1,
                              dihedral[3],
                              dihedral[1] * np.pi / 180.0,
                              flg_nb14)
                    cnt4_proper[quad[0]].add(quad[3])
                    self.add_torsion(quad, params)

        new_atoms = {}
        for atom_id, atom in gromtopol_mol.atoms.items():
            resname = atom[2]
            resid = atom[1]
            atom_name = atom[3]
            atom_type = atom[0].upper()
            charge = atom[5]
            mass = atom[6]
            atomtype_id = gromtopol.atomtypes[atom_type]
            if mass == None or mass <= 0.0:
                mass = gromtopol.nonbonds[atomtype_id][2]
            self.add_atom(atom_name, atom_type,
                          atomtype_id,
                          resname, resid,
                          mass,
                          gromtopol.nonbonds[atomtype_id][5] * 10.0 * (2.0**(1.0/6.0)) * 0.5,
                          charge,
                          sorted([x-atom_id for x in cnt2[atom_id]]),
                          sorted([x-atom_id for x in cnt3[atom_id]]),
                          sorted([x-atom_id for x in cnt4[atom_id]]),
                          0,0,0,0,0.0,0.0,0.0)
            
        self.bonds = sorted(self.bonds, key=lambda x:x[0][0])
        self.angles = sorted(self.angles, key=lambda x:x[0][0])
        self.torsions = sorted(self.torsions, key=lambda x:x[0][0])
        self.impropers = sorted(self.impropers, key=lambda x:x[0][0])
        return 
    def import_from_amber(self, ambertopol, amberff, neighbor_backbone=True):
        # ambertopol must be an instance of kkamber_ff.Topology
        for resname, res in ambertopol.res.items():
            self.import_from_amber_res(res, amberff, neighbor_backbone)
        return
    def import_from_amber_res(self, amberres, amberff, neighbor_backbone=True):
        ## for tplgene DB file
        # amberres must be an instance of kkamber_ff.Res
        atomid_bias = 0
        atomid_last= len(amberres.atoms.keys())
        ## neighbor backbone
        if neighbor_backbone:
            atomid_bias = 3
            atomid_last += 4
            new_connectivity = set()
            for cnt in amberres.connectivity:
                new_connectivity.add((cnt[0]+3, cnt[1]+3))
            amberres.connectivity = new_connectivity
            amberres.connectivity.add((1,2))
            amberres.connectivity.add((2,3))
            amberres.connectivity.add((2,4))
            amberres.connectivity.add((atomid_last-2, atomid_last))
            
        cnt2 = {}
        cnt3 = {}
        cnt4 = {}
        for i in range(1, atomid_last+1):
            cnt2[i] = set()
            cnt3[i] = set()
            cnt4[i] = set()
        ## enumerating connection 1-2, 1-3, 1-4
        amberres.enumerate_3_4_connectivity()

        for cnt in amberres.connectivity:
            cnt_pair = (cnt[0], cnt[1])
            if cnt[0] > cnt[1]: cnt_pair = (cnt[1], cnt[0])
            cnt2[cnt_pair[0]].add(cnt_pair[1])

        for cnt in amberres.connectivity3:
            cnt_pair = (cnt[0], cnt[2])
            if cnt[0] > cnt[2]: cnt_pair = (cnt[2], cnt[0])
            cnt3[cnt_pair[0]].add(cnt_pair[1])

        for cnt in amberres.connectivity4:
            cnt_pair = (cnt[0], cnt[3])
            if cnt[0] > cnt[3]: cnt_pair = (cnt[3], cnt[0])
            cnt4[cnt_pair[0]].add(cnt_pair[1])

        self.mol_name = amberres.resname
        self.mol_num = 1
        self.head_atom_id = 1
        
        new_atoms = {}
        if neighbor_backbone:
            ## renumber amberres.atoms:
            for atom_id, info in amberres.atoms.items():
                new_atoms[atom_id+atomid_bias] = info

            new_atoms[1] = ("CA","ct",0,0, 6, 0.0337)
            new_atoms[2] = ("C","c",0,0, 6, 0.5973)
            new_atoms[3] = ("O","o",0,0, 8, -0.5679)
            new_atoms[atomid_last] = ("N","n",0,0, 7, -0.4157)
        amberres.atoms = new_atoms
        for atom_id, atom in amberres.atoms.items():
            resname = amberres.resname
            resid = 1
            if neighbor_backbone and atom_id in [1,2,3]:
                resname = "PRE"
                resid = 0
            elif neighbor_backbone and atom_id == atomid_last:
                resname = "PRE"
                resid = 2

            self.add_atom(atom[0], atom[1].upper(),
                          kk.PrestoDef.Atom_type[atom[1].upper()],
                          resname, resid,
                          kk.KKDEF.ELEM[atom[4]].mass,
                          kk.KKDEF.ELEM[atom[4]].radius,
                          atom[5],
                          sorted([x-atom_id for x in cnt2[atom_id]]),
                          sorted([x-atom_id for x in cnt3[atom_id]]),
                          sorted([x-atom_id for x in cnt4[atom_id]]),
                          0,0,0,0,0.0,0.0,0.0)
        #for aid, amat in amberres.atoms.items():
        #    print str(aid) + " : " + " ".join([str(x) for x in amat])


        # BOND
        #print amberff.bond.keys()
        for cnt in amberres.connectivity:
            atomtype0 = amberres.atoms[cnt[0]][1].lower()
            atomtype1 = amberres.atoms[cnt[1]][1].lower()
            bond = None
            if (atomtype0, atomtype1) in amberff.bond:
                bond = amberff.bond[(atomtype0, atomtype1)]
            elif (atomtype1, atomtype0) in amberff.bond:
                bond = amberff.bond[(atomtype1, atomtype0)]
            if bond:
                self.add_bond(cnt, bond)
                #print bond
                
        # ANGLE
        #print amberff.angle.keys()
        for cnt in amberres.connectivity3:
            atomtype0 = amberres.atoms[cnt[0]][1].lower()
            atomtype1 = amberres.atoms[cnt[1]][1].lower()
            atomtype2 = amberres.atoms[cnt[2]][1].lower()
            angle = None
            #print (atomtype0, atomtype1, atomtype2)
            if (atomtype0, atomtype1, atomtype2) in amberff.angle:
                angle = amberff.angle[(atomtype0, atomtype1, atomtype2)]
            elif (atomtype2, atomtype1, atomtype0) in amberff.angle:
                angle = amberff.angle[(atomtype2, atomtype1, atomtype0)]
            if angle:
                self.add_angle(cnt, angle)

        # TORSION
        #print  amberres.connectivity4
        for cnt in amberres.connectivity4:
            atomtype0 = amberres.atoms[cnt[0]][1].lower()
            atomtype1 = amberres.atoms[cnt[1]][1].lower()
            atomtype2 = amberres.atoms[cnt[2]][1].lower()
            atomtype3 = amberres.atoms[cnt[3]][1].lower()
            atomtypes = (atomtype0, atomtype1, atomtype2, atomtype3)
            params = amberff.get_torsion_params(atomtypes, False)
            for param in params:
                periodicity = param[3]
                interact14 = 1
                if periodicity < 0:
                    periodicity *= -1
                    interact14 = 0
                new_param = (param[1], param[0], periodicity, param[2], interact14)
                self.add_torsion(cnt, new_param)

        for cnt in amberres.connectivity4_impr:
            atomtype0 = amberres.atoms[cnt[0]][1].lower()
            atomtype1 = amberres.atoms[cnt[1]][1].lower()
            atomtype2 = amberres.atoms[cnt[2]][1].lower()
            atomtype3 = amberres.atoms[cnt[3]][1].lower()
            atomtypes = (atomtype0, atomtype1, atomtype2, atomtype3)
            params = amberff.get_torsion_params(atomtypes, True)
            for param in params:
                periodicity = param[3]
                interact14 = 1
                if periodicity < 0:
                    periodicity *= -1
                    interact14 = 0
                new_param = (param[1], param[0], periodicity, param[2], interact14)
                self.add_improper(cnt, new_param)
        self.bonds = sorted(self.bonds, key=lambda x:x[0][0])
        self.angles = sorted(self.angles, key=lambda x:x[0][0])
        self.torsions = sorted(self.torsions, key=lambda x:x[0][0])
        self.impropers = sorted(self.impropers, key=lambda x:x[0][0])
        return 
    def add_atom(self, atom_name, atom_type,
                 interaction_type, res_name, res_id,
                 mass, radius, charge,
                 atoms_1_2, atoms_1_3, atoms_1_4,
                 z_pos_1_2, z_pos_1_3, z_pos_1_4,
                 z_basis, z_dist, z_ang, z_phase):
        self.atoms.append(TPLAtom(len(self.atoms)+1,
                                  atom_name, atom_type,
                                  interaction_type, res_name, res_id,
                                  mass, radius, charge,
                                  atoms_1_2, atoms_1_3, atoms_1_4,
                                  z_pos_1_2, z_pos_1_3, z_pos_1_4,
                                  z_basis, z_dist, z_ang, z_phase))

    def add_bond(self, atoms, params):
        assert(len(atoms) == 2)
        assert(len(params) == 2)
        self.bonds.append((tuple(atoms), params))
    def add_angle(self, atoms, params):
        assert(len(atoms) == 3)
        assert(len(params) == 2)
        self.angles.append((tuple(atoms), params))
        return
    def add_torsion(self, atoms, params):
        assert(len(atoms) == 4)
        assert(len(params) == 5)
        self.torsions.append((tuple(atoms),params))
        return
    def add_improper(self, atoms, params):
        assert(len(atoms) == 4)
        assert(len(params) == 5)
        self.impropers.append((tuple(atoms),params))
        return
    def set_scale_14(self, nonbonds):
        for atom, param in self.torsions:
            if param[4] != 1: continue
            type1 = self.atoms[atom[0]-1].interaction_type
            type4 = self.atoms[atom[3]-1].interaction_type
            ## FY14SE
            tor14param0 = min(nonbonds[type1][4], nonbonds[type4][4] ) * float(param[4])
            ## FY14SV FYTVWS
            tor14param1 = min(nonbonds[type1][5], nonbonds[type4][5] ) * float(param[4])
            self.tor14[atom] = (tor14param0, tor14param1)
        return 
    def get_text(self, header="PRE"):
        text = ""
        text += get_text_atoms(header)
        text += get_text_bonds(header)
        text += get_text_angles(header)
        text += get_text_torsions(header)
        text += get_text_impropers(header)
        return text
    def get_text_atoms_tpl(self, header="PRE"): # for .tpl file
        if len(self.atoms) == 0: return ""
        text = header+"> ATOMS\n"
        text += self.mol_name + "\n"
        for atom in self.atoms:
            text_atom = "%-8s %-5s %2d %-5s%8d %9.4f %9.4f %9.4f %2d %2d %2d -> ;  %d\n"\
                %(atom.atom_name, atom.atom_type, atom.interaction_type,\
                      atom.res_name, atom.res_id, atom.mass, atom.radius,\
                      atom.charge, len(atom.atoms_1_2), len(atom.atoms_1_3),\
                      len(atom.atoms_1_4), atom.atom_id)
            for a12 in atom.atoms_1_2:
                text_atom += "%7d"%a12
            for a13 in atom.atoms_1_3:
                text_atom += "%7d"%a13
            for a14 in atom.atoms_1_4:
                text_atom += "%7d"%a14
            text_atom += "        -> \n"
            text_atom += "%7d%7d%7d%7d %12.4f %12.4f %12.4f\n"\
                %(atom.z_pos_1_2, atom.z_pos_1_3, atom.z_pos_1_4,\
                      atom.z_basis, atom.z_dist, atom.z_ang, atom.z_phase)
            text += text_atom
        return text
    def get_text_atoms(self, header="PRE"): # for DB file
        if len(self.atoms) == 0: return ""
        text = header+"> ATOMS\n"
        text += self.mol_name + "\n"
        for atom in self.atoms:
            text_atom = "%-8s %-5s %2d %-5s%4d %8.3f %9.4f %8.4f %3d %3d %3d;  %d\n"\
                %(atom.atom_name, atom.atom_type, atom.interaction_type,\
                      atom.res_name, atom.res_id, atom.mass, atom.radius,\
                      atom.charge, len(atom.atoms_1_2), len(atom.atoms_1_3),\
                      len(atom.atoms_1_4), atom.atom_id)
            text_atom += "<"
            for a12 in atom.atoms_1_2:
                text_atom += "%4d"%a12
            for a13 in atom.atoms_1_3:
                text_atom += "%4d"%a13
            for a14 in atom.atoms_1_4:
                text_atom += "%4d"%a14
            text_atom += "\n<"
            text_atom += "%4d%4d%4d%4d %11.3f %11.3f %11.3f\n"\
                %(atom.z_pos_1_2, atom.z_pos_1_3, atom.z_pos_1_4,\
                      atom.z_basis, atom.z_dist, atom.z_ang, atom.z_phase)
            text += text_atom
        text += "\n"
        return text
    def get_text_bonds_tpl(self):
        if len(self.bonds) == 0: return ""
        text = "TPL> BONDS\n"
        text += self.mol_name + "\n"
        for i, bnd in enumerate(self.bonds):
            text_bond = "%10d%10d %14.7f %14.7f ;  %d\n"\
                %(bnd[0][0], bnd[0][1], bnd[1][0], bnd[1][1], i+1)
            text += text_bond
        text += "\n"
        return text
    def get_text_bonds(self, header="PRE"):        
        if len(self.bonds) == 0: return ""
        text = header+"> BONDS\n"
        text += self.mol_name + "\n"
        for i, bnd in enumerate(self.bonds):

            text_bond = "%6d%6d            %9.2f%9.4f;  %d\n"\
                %(bnd[0][0], bnd[0][1], bnd[1][0], bnd[1][1], i+1)
            text += text_bond
        text += "\n"
        return text
    def get_text_angles_tpl(self, header="PRE"):
        if len(self.angles) == 0: return ""
        text = header+"> ANGLES\n"
        text += self.mol_name + "\n"
        for i, ang in enumerate(self.angles):
            text_angle = "%10d%10d%10d %14.7f %14.7f ;  %d\n"\
                %(ang[0][0], ang[0][1], ang[0][2], ang[1][0], ang[1][1], i+1)
            text += text_angle
        text += "\n"
        return text
    def get_text_angles(self, header="PRE"):        
        if len(self.angles) == 0: return ""
        text = header + "> ANGLES\n"
        text += self.mol_name + "\n"
        for i, ang in enumerate(self.angles):
            text_angle = "%6d%6d%6d      %9.2f%9.4f;  %d\n"\
                %(ang[0][0], ang[0][1], ang[0][2], ang[1][0], ang[1][1], i+1)
            text += text_angle
        text += "\n"
        return text
    def get_text_torsions_tpl(self, header="TPL"):
        if len(self.torsions) == 0: return ""
        text = header+"> TORSIONS\n"
        text += self.mol_name + "\n"
        for i, trs in enumerate(self.torsions):
            text_trs = "%8d%8d%8d%8d %9.4f %4d %4d %9.4f %4d ; %d\n"\
                %(trs[0][0], trs[0][1], trs[0][2], trs[0][3],\
                      trs[1][0], trs[1][1], trs[1][2],\
                      trs[1][3], trs[1][4], i+1)
            text += text_trs
        text += "\n"
        return text
    def get_text_torsions(self, header="PRE"):        
        if len(self.torsions) == 0: return ""
        text = header + ">TORSIONS\n"
        text += self.mol_name + "\n"
        for i, trs in enumerate(self.torsions):
            text_trs = "%6d%6d%6d%6d%11.4f%5d%5d%9.3f%5d;  %d\n"\
                %(trs[0][0], trs[0][1], trs[0][2], trs[0][3],\
                      trs[1][0], trs[1][1], trs[1][2],\
                      trs[1][3], trs[1][4], i+1)
            text += text_trs
        text += "\n"
        return text
    def get_text_impropers_tpl(self, header="PRE"): 
        if len(self.impropers) == 0: return ""
        text = header + "> IMPROPER-TORSIONS\n"
        text += self.mol_name + "\n"
        for i, trs in enumerate(self.impropers):
            text_trs = "%8d%8d%8d%8d %9.4f %4d %4d %9.4f %4d ; %d\n"\
                %(trs[0][0], trs[0][1], trs[0][2], trs[0][3],\
                      trs[1][0], trs[1][1], trs[1][2],\
                      trs[1][3], trs[1][4], i+1)
            text += text_trs
        return text
    def get_text_impropers(self, header="PRE"): 
        if len(self.impropers) == 0: return ""
        text = header + "> IMPROPER-TORSIONS\n"
        text += self.mol_name + "\n"
        for i, trs in enumerate(self.impropers):
            text_trs = "%6d%6d%6d%6d%11.4f%5d%5d%9.3f%5d;  %d\n"\
                %(trs[0][0], trs[0][1], trs[0][2], trs[0][3],\
                      trs[1][0], trs[1][1], trs[1][2],\
                      trs[1][3], trs[1][4], i+1)
            text += text_trs
        return text
    def shift_atomid(self, offset):
        for atom in self.atoms:
            atom.atom_id += offset
        for bond in self.bonds:
            bond[0] = [ bond[0][i]+offset for i in range(0,2) ]
        for angle in self.angles:
            angle[0] = [ angle[0][i]+offset for i in range(0,3)]
        for torsion in self.torsions:
            torsion[0] = [ torsion[0][i]+offset for i in range(0,4) ]
        for impro in self.impropers:
            impro[0] = [ impro[0][i]+offset for i in range(0,4) ]
        for tor14 in self.to14:
            tor14[0] = [ tor14[0][i]+offset for i in range(0,4) ]
        return

class TPL(object):
    def __init__(self):
        self.mols = []  #mols[] = TplMol
        self.title = ""
        self.nonbonds = {}
        self.functions = {}
        self.nb_pair = {}
        self.nb_pair_hps = {}
        self.atom_id_12 = {}
        self.atom_id_13 = {}
        self.atom_id_14 = {}
        self.atom_id_14_imp = {}
        self.atom_id_14nb = {}
        self.atom_pair_non15 = set()
        self.zm_excess_pairs = set()

        self.enumerated = False
        ## Flag for completion of the enumerate_12_13_14()

    def get_mol_by_name(name):
        for m in self.mols:
            if m.mol_name == name: return m
        return None
    def diff(self, tpl):
        mol1 = set([mol.mol_name for mol in self.mols])
        mol2 = set([mol.mol_name for mol in tpl.mols])
        print "======== MOLECULES ========"
        print "== MATCH"
        print ", ".join(mol1.union(mol2))
        print "== Only in Tpl1 "
        print ", ".join([m for m in mol1 if not m in mol2])
        print "== Only in Tpl2 "
        print ", ".join([m for m in mol2 if not m in mol1])
        print ""
        for m in mol1.union(mol2):
            self.get_mol_by_name(m).diff(tpl.get_mol_by_name(m))
        return 
    def n_atoms(self):
        na = 0
        for mol in self.mols:
            na += len(mol.atoms) * mol.mol_num
        return na
    def add_mol(self, name, num):
        head_atom_id = self.n_atoms()
        self.mols.append(TPLMol(name, num, head_atom_id))
        #print "add_mol " + name + " " + str(head_atom_id) 
        return self.mols
    def convert_units_to_si(self):
        for mol in self.mols:
            for atom in mol.atoms:
                atom.mass *= 1e-3
    def add_function(self, function_id, n_params, name):
        self.functions[function_id] = (n_params, name)
    def add_nonbond(self, atom, params):
        assert(len(params) >= 6)
        self.nonbonds[atom] = tuple(params)
    def set_nonbond_pair(self):
        self.nb_pair = {}
        for i,param_i in self.nonbonds.items():
            for j,param_j in self.nonbonds.items():
                if i<j: break
                pair = (j,i)
                #param_x[2] ... sigma
                #param_x[3] ... epsiron

                p3_ave = np.sqrt((param_i[3]*param_j[3]))
                p2_ij6 = (param_i[2]+param_j[2])**6                
                p2_ij12 = p2_ij6 * p2_ij6
                param1 = p3_ave * p2_ij6 * 2.0
                param2 = p3_ave * p2_ij12
                self.nb_pair[pair] = (param1, param2)
                self.nb_pair[(pair[1], pair[0])] = (param1, param2)

                ## for HPS potential
                cutoff = np.power(2.0, 1.0/6.0) * (param_i[3] + param_j[3]) * 0.5
                lmb = (param_i[6] + param_j[6]) * 0.5
                self.nb_pair_hps[(pair[1], pair[0])] = (cutoff, lmb)
                self.nb_pair_hps[(pair[1], pair[0])] = (cutoff, lmb)

#                print "nb_pair: " + str(pair[0]) + " : " + str(pair[1])
#                print str(param1) + ", " + str(param2)
    def set_scale_14(self):
        for mol in self.mols:
            mol.set_scale_14(self.nonbonds)
    def swap(self,a,b):
        if a > b:
            tmp = a
            a = b
            b = tmp
        return (a,b)
    def enumerate_12_13_14(self):
        """
        The bonded information is recorded in each molecule.
        The instances of each bonded information is enumerated.
        For example, when the molecule A consisting of 10 atoms 
        have 3 instances in the system,
        the instances of bond 1-2 in the molecule A should be
          1-2
          11-12
          21-22
        """
    ## (atom_id1, atom_id2, atom_id1_in_mol, atom_id2_in_mol, mol_id)
        atom_id_12 = []
        atom_id_13 = []
        atom_id_14 = []
        atom_id_14_imp = []
        atom_pair_non15 = set()
        atom_id_14nb = []
        for mol_tpl in self.mols:
            for mol_i in range(0, mol_tpl.mol_num):
                head_atom_id = mol_tpl.head_atom_id + (len(mol_tpl.atoms) * mol_i)
                for pair, params in mol_tpl.bonds:
                    atom1 = head_atom_id-1 + pair[0]
                    atom2 = head_atom_id-1 + pair[1]
                    atom_id_12.append(((atom1,atom2), params))
                    atom_pair_non15.add((atom1,atom2))
                for trio, params in mol_tpl.angles:
                    atom1 = head_atom_id-1 + trio[0]
                    atom2 = head_atom_id-1 + trio[1]
                    atom3 = head_atom_id-1 + trio[2]
                    atom_id_13.append(((atom1,atom2,atom3), params))
                    pair = self.swap(atom1,atom3)
                    atom_pair_non15.add(pair)
                for mol_atom_id, params in mol_tpl.torsions:
                    atom1 = head_atom_id-1 + mol_atom_id[0]
                    atom2 = head_atom_id-1 + mol_atom_id[1]
                    atom3 = head_atom_id-1 + mol_atom_id[2]
                    atom4 = head_atom_id-1 + mol_atom_id[3]
                    atom_id_14.append(((atom1,atom2,atom3,atom4), params))
                    if params[4] != 1: continue
                    pair = self.swap(atom1,atom4)
                    atom_pair_non15.add(pair)
                    #params_nb =self.nb_pair[(mol_tpl.atoms[mol_atom_id[0]].interaction_type,
                    #                          mol_tpl.atoms[mol_atom_id[0]].interaction_type)]
                    params_14nb = (mol_tpl.atoms[mol_atom_id[0]-1].interaction_type,
                                   mol_tpl.atoms[mol_atom_id[3]-1].interaction_type,
                                   mol_tpl.tor14[mol_atom_id][0],
                                   mol_tpl.tor14[mol_atom_id][1])
                    
                    atom_id_14nb.append((pair, params_14nb))
                for quad, params in mol_tpl.impropers:
                    atom1 = head_atom_id-1 + quad[0]
                    atom2 = head_atom_id-1 + quad[1]
                    atom3 = head_atom_id-1 + quad[2]
                    atom4 = head_atom_id-1 + quad[3]
                    atom_id_14_imp.append(((atom1,atom2,atom3,atom4), params))
                    pair = self.swap(atom1,atom4)
                    atom_pair_non15.add(pair)
                #for atom in mol_tpl.atoms:
                #    for dest in atom.atoms_1_2:
                #        atom_id_12.append((head_atom_id-1 + atom.atom_id,
                #                           head_atom_id-1 + atom.atom_id + dest,
                #                           atom.atom_id,
                #                           atom.atom_id + dest,
                #                           mol_i))
                #        atom_pair_non15.add((head_atom_id-1 + atom.atom_id,
                #                             head_atom_id-1 + atom.atom_id + dest))
                #    for dest in atom.atoms_1_3:
                #        atom_id_13.append((head_atom_id-1 + atom.atom_id,
                #                           head_atom_id-1 + atom.atom_id + dest[0],
                #                           head_atom_id-1 + atom.atom_id + dest[1],
                #                           atom.atom_id,
                #                           atom.atom_id + dest[0],
                #                           atom.atom_id + dest[1],
                #                           mol_i))
                #        atom_pair_non15.add((head_atom_id-1 + atom.atom_id,
                #                             head_atom_id-1 + atom.atom_id + dest[1]))
                #    for dest in atom.atoms_1_4:
                #        atom_id_14.append((head_atom_id-1 + atom.atom_id,
                #                           head_atom_id-1 + atom.atom_id + dest[0],
                #                           head_atom_id-1 + atom.atom_id + dest[1],
                #                           head_atom_id-1 + atom.atom_id + dest[2],
                #                           atom.atom_id,
                #                           atom.atom_id + dest[0],
                #                           atom.atom_id + dest[1],
                #                           atom.atom_id + dest[2],
                #                           mol_i))
                #        atom_pair_non15.add((head_atom_id-1 + atom.atom_id,
                #                             head_atom_id-1 + atom.atom_id + dest[2]))
        self.atom_id_12 = atom_id_12
        self.atom_id_13 = atom_id_13
        self.atom_id_14 = atom_id_14
        self.atom_id_14_imp = atom_id_14_imp
        self.atom_pair_non15 = atom_pair_non15
        self.atom_id_14nb = atom_id_14nb
        #print "1-2 atoms:"
        #print self.atom_id_12
        self.enumerated = True
        return self.atom_id_12, self.atom_id_13, self.atom_id_14, self.atom_id_14_imp
    def import_from_gromacs(self, gromtopol):
        for atomtype_id, atom in enumerate(gromtopol.nonbonds):
            params = [0, 1,
                      atom[5] * 10.0 * (2.0**(1.0/6.0)) * 0.5,
                      atom[6] * J_TO_CAL,
                      gromtopol.def_fudge_lj,
                      gromtopol.def_fudge_qq]
            self.add_nonbond(atomtype_id, params)
        for molname, n_mol in gromtopol.n_molecules.items():
            self.add_mol(molname, n_mol)
            self.mols[-1].import_from_gromacs_mol(gromtopol.moltypes[molname], gromtopol)
    def get_text(self, header="TPL"):
        text = header + "> TITLE\n"
        if self.title != "":
            text += self.title + "\n"
        else: text += "MD\n"
        text += "\n"

        text += header + "> MOLECULES\n"
        for mol in self.mols:
            if mol.mol_num == 0: continue
            text += "%-42s%-d\n"%(mol.mol_name, mol.mol_num)
        text += "\n"

        for mol in self.mols:
            if mol.mol_num == 0: continue            
            text += mol.get_text_atoms_tpl(header)
        for mol in self.mols:
            if mol.mol_num == 0: continue            
            text += mol.get_text_bonds_tpl(header)
        for mol in self.mols:
            if mol.mol_num == 0: continue            
            text += mol.get_text_angles_tpl(header)
        for mol in self.mols:
            if mol.mol_num == 0: continue            
            text += mol.get_text_torsions_tpl(header)
        for mol in self.mols:
            if mol.mol_num == 0: continue            
            text += mol.get_text_impropers_tpl(header)
        text += "\n"
        text += header + "> FUNCTIONS\n"    
        text += "         1         4          LENNARD-JONES-AMBER\n"
        text += self.get_text_nonbonds(header)
        return text
    def get_text_nonbonds(self, header="TPL"):
        text = header + "> NONBONDS\n"
        #00000000001111111111222222222233333333334444444444555555555566666666667777777777
        #01234567890123456789012345678901234567890123456789012345678901234567890123456789
        #    1    0   1      1.9080000      0.0860000      0.8333333      0.5000000
        for atom_type, param in self.nonbonds.items():
            text += "%5d%5d%4d%15.7f%15.7f%15.7f%15.7f"%(atom_type, param[0],param[1], param[2], param[3], param[4], param[5])
            if len(param)>=7:
                text += " %8.5f"%(param[6])
                if len(param)>=8:
                    text += " ; " + param[7]
            text += "\n"
        return text
    def reduce_nonbond_parameters(self):
        """
        Remove nonbond atomtypes that did not appear in the mol
        """

        ## detect atom types with same paramters
        dup_types = {}
        for atom_type1, params1 in self.nonbonds.items():
            for atom_type2, params2 in self.nonbonds.items():
                if atom_type1 == atom_type2: break
                if params1 == params2:
                    dup_types[atom_type1] = atom_type2
                    break
            if not atom_type1 in dup_types:
                dup_types[atom_type1] = atom_type1

        used_types = {}
        for mol in self.mols:
            for atom in mol.atoms:
                used_types[dup_types[atom.interaction_type]] = -1
        new = 1
        for old in used_types.keys():
            used_types[old] = new
            new += 1

        for mol in self.mols:
            for atom in mol.atoms:
                atom.interaction_type = used_types[dup_types[atom.interaction_type]]
        new_nonbonds = {}
        for old, new in used_types.items():
            param = self.nonbonds[old]
            new_nonbonds[new] = param
        self.nonbonds = new_nonbonds
        return 0
    def remove_bond_angle_constraints(self, pairs):
        #print "remove_bond_angle_constraints"

        n_bonds_orig = len(self.atom_id_12)
        n_angles_orig = len(self.atom_id_13)

        new_bonds = []
        #new_non15 = set()
        print "n pairs : " + str(len(pairs))
        #print pairs
        for idx, at12 in enumerate(self.atom_id_12):
            #print at12[0]
            if not at12[0] in pairs:
                
                new_bonds.append(at12)
                #new_non15.add(at12[0])
        self.atom_id_12 = new_bonds
                #print "idx " + str(idx)
                #rmv_idx_bond.append(idx)
        new_angles = []
        for idx, at13 in enumerate(self.atom_id_13):
            if not (at13[0][0], at13[0][2]) in pairs:
                new_angles.append(at13)
                #new_non15.add(self.swap(at13[0][0], at13[0][2]))
        self.atom_id_13 = new_angles
        #self.atom_pair_non15 = new_non15
        print "n_bonds: " + str(n_bonds_orig) + " -> " + str(len(self.atom_id_12))
        print "n_angles: " + str(n_angles_orig) + " -> " + str(len(self.atom_id_13))

        return 


class PrestoAsciiReader(kkkit.FileI):
    def __init__(self, fn):
        super(PrestoAsciiReader, self).__init__(fn)
    def readline_comment(self):
        try:
            line = self.f.readline()
            if line=="":  return None
        except: return None
        i = len(line)
        try:
            i = line.index(';')
        except ValueError:
            i = len(line) 
        line = line[0:i] + '\n'
        return line
    def readline(self):
        line = ""
        while line == "":
            line = self.readline_comment()
            if not line: return None
            line = line.strip()

        while len(line)>=2 and line.strip()[-2:] == "->":
            line = line[:-2]
            flg = True
            try:    new_line = self.readline_comment()
            except: return line
            if not new_line: return line
            line += ' ' + new_line.strip()
        return line
        
class TPLReader(PrestoAsciiReader):
    def __init__(self, fn):
        super(TPLReader, self).__init__(fn)
    def read_tpl(self):
        self.open()
        tpl = TPL()        

        reading_area = 0
        ## 0 = default
        ## 1 = ATOM
                     
        reading_mol_id = -1
        reading_mol = ""
        mol_name_id = {}
        line = self.readline()
        while line:
            terms = re.compile("\s+").split(line.strip())
            if terms[0] == "TPL>" or terms[0] == "PRE>":
                reading_area = ""
                if terms[1] == "TITLE":
                    title_line = self.readline()
                    tpl.title = title_line
                    if DEBUG: print tpl.title
                elif terms[1] == "MOLECULES":
                    reading_area = terms[1]
                elif terms[1] == "ATOMS":
                    reading_mol = self.readline().strip()
                    reading_mol_id = mol_name_id[reading_mol] 
                    reading_area = terms[1]
                    tpl.mols[reading_mol_id].head_atom_id = tpl.n_atoms()
                elif terms[1] == "BONDS":
                    reading_mol = self.readline().strip()
                    reading_mol_id = mol_name_id[reading_mol] 
                    reading_area = terms[1]
                elif terms[1] == "ANGLES":
                    reading_mol = self.readline().strip()
                    reading_mol_id = mol_name_id[reading_mol] 
                    reading_area = terms[1]
                elif terms[1] == "TORSIONS":
                    reading_mol = self.readline().strip()
                    reading_mol_id = mol_name_id[reading_mol] 
                    reading_area = terms[1]
                elif terms[1] == "IMPROPER-TORSIONS":
                    reading_mol = self.readline().strip()
                    reading_mol_id = mol_name_id[reading_mol] 
                    reading_area = terms[1]
                elif terms[1] == "FUNCTIONS":
                    reading_area = terms[1]
                elif terms[1] == "NONBONDS":
                    reading_area = terms[1]
            elif reading_area == "MOLECULES":
                mol_terms = re.compile("\s+").split(line.strip())
                mol_name_id[mol_terms[0]] = len(tpl.mols)
                tpl.add_mol(mol_terms[0], int(mol_terms[1]))
                if DEBUG: print "New mol:" + mol_terms[0] + ": " + str(mol_terms[1])
            elif reading_area == "ATOMS":
                #print " ".join(terms)
                n_1_2_atoms = int(terms[8]) #n_1_2_atoms
                n_1_3_atoms = int(terms[9]) #n_1_3_atoms
                n_1_4_atoms = int(terms[10]) #n_1_4_atoms
                n_terms = 11 + n_1_2_atoms + n_1_3_atoms + n_1_4_atoms + 7
                while len(terms) < n_terms:
                    line2 = self.readline().strip()
                    if line2[0] == "<":
                        line2 = line2[1:].strip()
                        print "ERROR: the number of terms in the ATOM field is not enough"
                        return
                    terms2 = re.compile("\s+").split(line2)
                    terms.extend(terms2)

                term_1_2 = 11
                term_1_3 = 11+n_1_2_atoms
                term_1_4 = 11+n_1_2_atoms+n_1_3_atoms
                term_z = term_1_4 + n_1_4_atoms
                atoms_1_2 = [int(x) for x in terms[term_1_2:(term_1_2+n_1_2_atoms)]]
                atoms_1_3 = [int(x) for x in terms[term_1_3:(term_1_3+n_1_3_atoms)]]
                atoms_1_4 = [int(x) for x in terms[term_1_4:(term_1_4+n_1_4_atoms)]]
                tpl.mols[reading_mol_id].add_atom(terms[0], #atom_name
                                                  terms[1], #atom_type
                                                  int(terms[2]), #interaction_type
                                                  terms[3], #res_name
                                                  int(terms[4]), #res_id
                                                  float(terms[5]), #mass
                                                  float(terms[6]), #radius
                                                  float(terms[7]), # charge
                                                  atoms_1_2, atoms_1_3, atoms_1_4,
                                                  int(terms[term_z]), #z_pos_1_2
                                                  int(terms[term_z+1]), #z_pos_1_3
                                                  int(terms[term_z+2]), #z_pos_1_4,
                                                  int(terms[term_z+3]), #z_basis
                                                  float(terms[term_z+4]), #z_dist
                                                  float(terms[term_z+5]), #z_ang
                                                  float(terms[term_z+6])) #z_phase
            elif reading_area == "BONDS":
                atoms = [int(x) for x in terms[0:2]]
                params = [float(x) for x in terms[2:4]]
                ## params[0] ... energy kcal/mol
                ## params[1] ... distance A
                tpl.mols[reading_mol_id].add_bond(atoms, params)
            elif reading_area == "ANGLES":
                atoms = [int(x) for x in terms[0:3]]
                param_e = float(terms[3])
                param_ang = float(terms[4]) * np.pi / 180.0
                params = (param_e, param_ang)
                ## params[0] ... energy
                ## params[1] ... angle
                tpl.mols[reading_mol_id].add_angle(atoms, params)
            elif reading_area == "TORSIONS":
                atoms = [int(x) for x in terms[0:4]]
                params = [float(terms[4]), int(terms[5]), int(terms[6]),
                          float(terms[7]), int(terms[8])]
                ## FYFTOR in Presto: params[0] ... energy
                ## IYTDIV: params[1] ... number of torsions overlapping with this
                ## FYFROT: params[2] ... symmetry
                ## FYFPHS: params[3] ... phase (degree)
                ## IYTNBF: params[4] ... flag for calculating 1-4 vdw and electrostatic
                ##               (0 or 1)

                params[3] *= np.pi / 180.0
                # convert from degree to radian

                tpl.mols[reading_mol_id].add_torsion(atoms, params)
            elif reading_area == "IMPROPER-TORSIONS":
                atoms = [int(x) for x in terms[0:4]]
                params = [float(terms[4]), int(terms[5]), int(terms[6]),
                          float(terms[7]), int(terms[8])]
                params[3] *= np.pi / 180.0

                tpl.mols[reading_mol_id].add_improper(atoms, params)
            elif reading_area == "FUNCTIONS":
                # function_id
                # n_params
                # name
                tpl.add_function(int(terms[0]), int(terms[1]), terms[2])
            elif reading_area == "NONBONDS":
                atom = int(terms[0])
                params = [int(terms[1]), int(terms[2]),
                          float(terms[3]), float(terms[4]),
                          float(terms[5]), float(terms[6])]
                ## params[0] ... always 0
                ## params[1] ... always 1
                ## FYVWRD: params[2] ... distance in angstrome
                ## FYVWME: params[3] ... energy in kcal/mol
                ## FY14SE: params[4] ... energy for 1-4 vdw
                ## FY14SV: params[5] ... energy for 1-4 electrostatic

                ## HPS potential
                if len(terms) >= 7:
                    params.append(float(terms[7]))

                tpl.add_nonbond(atom, params)
            line = self.readline()
        self.close()

        tpl.set_nonbond_pair()
        tpl.set_scale_14()
        return tpl

class InputReader(PrestoAsciiReader):
    def __init__(self, fn):
        self.STRING = 0
        self.INT = 1
        self.FLOAT = 2
        self.KEYWORD = 3

        super(TPLReader, self).__init__(fn)
        self.config_def["INPUT"] = {"RESTI":self.STRING, 
                                    "TPL":self.STRING,
                                    "INITCO":self.STRING,
                                    "LXCELL":self.FLOAT,
                                    "LYCELL":self.FLOAT,
                                    "LZCELL":self.FLOAT,
                                    "DUMMY":self.KEYWORD
                                    }
        self.config = {}
    def read_input(self):
        self.open()
        line = self.readline()
        while line:
            terms = re.compile("\s+").split(line.strip())
            if not terms[0] in self.config_def:
                print("The field is not defined. :" + terms[0])
                sys.exit(0)
            self.configs[terms[0]] = terms[1]
            line = self.readline()
        self.close()
        return self.config

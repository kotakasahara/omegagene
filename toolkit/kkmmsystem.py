import numpy as np
import kkpbc
import kkstruct 
import copy
import kkmmff
import kkpresto
import kkmmeneflo
import kkpresto_shake as shk

class MmSystem(object):
    def __init__(self, model,
                 cell_x, cell_y, cell_z):
        ## parameters 
        self.cutoff = 0.0
        self.cutoff_2 = 0.0
        self.cutoff_i2 = 0.0

        self.model = model
        self.shake = None

        if cell_x and cell_y and cell_z:
            self.pbc = kkpbc.PBC(cell_x, cell_y, cell_z)
        else:
            self.pbc = kkpbc.PBC(model.pbc_box[0],
                                 model.pbc_box[1],
                                 model.pbc_box[2])
        self.n_atoms = len(model.atoms)

        self.ff = kkmmff.FF()
        self.eneflo = kkmmeneflo.Eneflo()
        self.eneflo.set_params(len(model.atoms))
        ## for each atom

        self.charge = np.zeros(self.n_atoms, np.float64)
        self.mass = np.zeros(self.n_atoms, np.float64)
        self.atomtype = np.zeros(self.n_atoms, np.int32)

        self.crd = self.init_3d_array(self.n_atoms)
        self.vel = self.init_3d_array(self.n_atoms)
        self.force = self.init_3d_array(self.n_atoms)
        self.crd_int = self.init_3d_array(self.n_atoms, np.int32)
        #self.pair_flags = self.init_matrix(self.n_atoms, self.n_atoms, np.int32)
        self.set_crd_from_model(self.model)
        self.reset_ene()
    def reset_ene(self):
        self.work_bond = self.init_3d_array(self.n_atoms, np.float64)
        self.work_angle = self.init_3d_array(self.n_atoms)
        self.work_torsion = self.init_3d_array(self.n_atoms)
        self.work_improper = self.init_3d_array(self.n_atoms)
        self.work_14vdw = self.init_3d_array(self.n_atoms)
        self.work_14ele = self.init_3d_array(self.n_atoms)
        self.work_vdw = self.init_3d_array(self.n_atoms)
        self.work_ele = self.init_3d_array(self.n_atoms)

        ## for total in system
        self.temperature = 0.0
        self.kinetic = 0.0
        self.pote_bond = 0.0
        self.pote_angle = 0.0
        self.pote_torsion = 0.0
        self.pote_improper = 0.0
        self.pote_14vdw = 0.0
        self.pote_14ele = 0.0
        self.pote_vdw = 0.0
        self.pote_ele = 0.0
        self.pote_ele_self = 0.0
        self.pote_ele_excess12 = 0.0
        self.pote_ele_excess13 = 0.0
        self.pote_ele_excess14 = 0.0
        self.kinetic = 0.0
        ## set crd from model

        return
    def set_params(self, cutoff):
        self.cutoff = cutoff
        self.cutoff_2 = cutoff**2
        cutoff_i = int(cutoff + 1.5)
        self.cutoff_i2 = cutoff_i**2
        return 
    def init_matrix(self, n, m, elem_type=np.float64):
        arr = []
        for i in range(0, n):
            arr.append(np.zeros(m, elem_type))
        return np.array(arr)
    def init_3d_array(self, n, elem_type=np.float64):
        arr = []
        for i in range(0, n):
            arr.append(np.zeros(3, elem_type))
        return np.array(arr)
    def init_array(self, n):
        arr = []
        for i in range(0, n):
            arr.append(0.0)
        return np.array(arr)

    def set_crd_from_model(self, model):
        for i,atom in enumerate(model.atoms):
            self.crd[i] = copy.deepcopy(atom.crd)
        return
    def set_crd_vel_from_restart(self, restart):
        for i, crd in enumerate(restart.crd):
            self.crd[i] = copy.deepcopy(crd)
            self.vel[i] = copy.deepcopy(restart.vel[i])
        return
    def set_atom_info_from_tpl(self, tpl):
        i = -1
        for mol in tpl.mols:
            for i_mol in range(0, mol.mol_num):
                for atom in mol.atoms:
                    i += 1
                    self.mass[i] = atom.mass
                    self.charge[i] = atom.charge
                    self.atomtype[i] = atom.interaction_type
        return 

    def store_self_energy(self, self_energy):
        self.pote_ele_self = np.sum(self_energy)
        self.pote_ele += self.pote_ele_self
        return
    def store_bond_energy(self, pair, ene, work):
        self.pote_bond += ene
        self.force[pair[0]] += work
        self.force[pair[1]] -= work
        return 
    def store_angle_energy(self, trio, ene, work1, work2):
        self.pote_angle += ene
        self.force[trio[0]] += work1
        self.force[trio[1]] -= work1 + work2
        self.force[trio[2]] += work2
        return 
    def store_torsion_energy(self, quad, ene, work1, work2, work3):
        self.pote_torsion += ene
        self.force[quad[0]] += work1
        self.force[quad[1]] += work2
        self.force[quad[2]] -= work1 + work2 + work3
        self.force[quad[3]] += work3
        return 
    def store_improper_energy(self, quad, ene, work1, work2, work3):
        self.pote_improper += ene
        self.force[quad[0]] += work1
        self.force[quad[1]] += work2
        self.force[quad[2]] -= work1 + work2 + work3
        self.force[quad[3]] += work3
        return 
    def store_14vdw_energy(self, pair, ene, work):
        self.pote_14vdw += ene
        self.force[pair[0]] += work
        self.force[pair[1]] -= work
        return 
    def store_14ele_energy(self, pair, ene, work):
        self.pote_14ele += ene
        self.force[pair[0]] += work
        self.force[pair[1]] -= work
        return 
    def store_pair_vdw_energy(self, pair, ene, work):
        self.pote_vdw += ene
        self.force[pair[0]] += work
        self.force[pair[1]] -= work
        return 
    def store_pair_ele_energy(self, pair, ene, work):
        self.pote_ele += ene
        self.force[pair[0]] += work
        self.force[pair[1]] -= work
        return 
    def print_energies(self):
        line = "RESULT\n"
        line += "======\n"
        line += "\n"
        line += "Potential Energy\n"
        line += "----------------\n"
        line += "Potential  : %4.15e\n"%(self.pote_bond+self.pote_angle+self.pote_torsion+self.pote_improper+self.pote_14vdw+self.pote_14ele+self.pote_vdw+self.pote_ele)
        line += "  Bond     : %4.15e\n"%self.pote_bond
        line += "  Angle    : %4.15e\n"%self.pote_angle
        line += "  Torsion  : %4.15e\n"%self.pote_torsion
        line += "  Improper : %4.15e\n"%self.pote_improper
        line += "  14-VDW   : %4.15e\n"%self.pote_14vdw
        line += "  14-Ele   : %4.15e\n"%self.pote_14ele
        line += "  VDW      : %4.15e\n"%self.pote_vdw
        line += "  Ele      : %4.15e\n"%self.pote_ele
        line += "Kinetic    : %4.15e\n"%self.kinetic
        return line
        

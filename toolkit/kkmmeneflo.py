import numpy as np

class Eneflo(object):
    def __init__(self):
        self.n_atoms = 0
        
        self.pair_terms_atoms = []
        ## [(atom_id1, atom_id2), .... ]
        self.pair_terms_ene = []
        ## [ene, ene, ...]
        self.pair_terms_work = []
        ## [work, work, ...]

        self.pair_ene_i = np.zeros(1)
        self.pair_eneflo_kine = np.zeros(1)
        self.pair_eneflo_pote = np.zeros(1)

        self.trio_ene_i = np.zeros(1)
        self.trio_eneflo_kine = np.zeros(1)
        self.trio_eneflo_pote = np.zeros(1)

        self.quad_ene_i = np.zeros(1)
        self.quad_eneflo_kine = np.zeros(1)
        self.quad_eneflo_pote = np.zeros(1)

        return 
    def set_params(self, n_atoms):
        self.n_atoms = n_atoms
        self.reset_eneflo_values()
    def reset_eneflo_values(self):
        self.pair_terms_atoms = []
        self.pair_terms_ene = []
        self.pair_terms_work = []
        self.trio_terms_atoms = []
        self.trio_terms_ene = []
        self.trio_terms_work1 = []
        self.trio_terms_work2 = []
        self.quad_terms_atoms = []
        self.quad_terms_ene = []
        self.quad_terms_work1 = []
        self.quad_terms_work2 = []
        self.quad_terms_work3 = []
        self.pair_ene_i = np.zeros(self.n_atoms)
        self.pair_eneflo_kine = np.zeros(self.n_atoms)
        self.pair_eneflo_pote = np.zeros(self.n_atoms)

        self.trio_ene_i = np.zeros(self.n_atoms)
        self.trio_eneflo_kine = np.zeros(self.n_atoms)
        self.trio_eneflo_pote = np.zeros(self.n_atoms)
        self.quad_ene_i = np.zeros(self.n_atoms)
        self.quad_eneflo_kine = np.zeros(self.n_atoms)
        self.quad_eneflo_pote = np.zeros(self.n_atoms)

        return
    def append_pair_terms(self, atoms, ene, work):
        self.pair_terms_atoms.append(atoms)
        self.pair_terms_ene.append(ene)
        self.pair_terms_work.append(work)
        return
    def append_trio_terms(self, atoms, ene, work1, work2):
        self.trio_terms_atoms.append(atoms)
        self.trio_terms_ene.append(ene)
        self.trio_terms_work1.append(work1)
        self.trio_terms_work2.append(work2)
        return
    def append_quad_terms(self, atoms, ene, work1, work2, work3):
        self.quad_terms_atoms.append(atoms)
        self.quad_terms_ene.append(ene)
        self.quad_terms_work1.append(work1)
        self.quad_terms_work2.append(work2)
        self.quad_terms_work3.append(work3)
        return

    def calc_eneflo(self, vel):
        self.calc_pair(vel, self.pair_terms_atoms,
                       np.array(self.pair_terms_ene) * 0.5,
                       np.array(self.pair_terms_work) * 0.5)
        self.calc_trio(vel, self.trio_terms_atoms,
                       np.array(self.trio_terms_ene) * 1/3,
                       np.array(self.trio_terms_work1) * 1/3,
                       np.array(self.trio_terms_work2) * 1/3)
        self.calc_quad(vel, self.quad_terms_atoms,
                       np.array(self.quad_terms_ene) * 0.25,
                       np.array(self.quad_terms_work1) * 0.25,
                       np.array(self.quad_terms_work2) * 0.25,
                       np.array(self.quad_terms_work3) * 0.25)
        
        
    def calc_pair(self, vel, atoms,ene, work):
        for i, pair in enumerate(atoms):
            atom_id1, atom_id2 = pair
            self.pair_ene_i[atom_id1] += ene[i]
            self.pair_ene_i[atom_id2] += ene[i]
            self.pair_eneflo_kine[atom_id1] += np.dot(vel[atom_id1], work[i])
            self.pair_eneflo_kine[atom_id2] += np.dot(vel[atom_id2], -1.0*work[i])
            self.pair_eneflo_pote[atom_id1] -= np.dot(vel[atom_id2], -work[i])
            self.pair_eneflo_pote[atom_id2] -= np.dot(vel[atom_id1], 1.0 * work[i])
        return 
    def calc_trio(self, vel, atoms, ene, work1, work2):
        
        for i, trio in enumerate(atoms):
            atom_id1, atom_id2, atom_id3 = trio
            self.trio_ene_i[atom_id1] += ene[i]
            self.trio_ene_i[atom_id2] += ene[i]
            self.trio_ene_i[atom_id3] += ene[i]
            #print "dbg:"
            #print work1
            #print vel[atom_id1]
            #print np.dot(work1,vel[atom_id1])
            self.trio_eneflo_kine[atom_id1] += np.dot(work1[i], vel[atom_id1]) * 2.0
            self.trio_eneflo_kine[atom_id2] += np.dot(-work1[i]-work2[i],vel[atom_id2]) * 2.0
            self.trio_eneflo_kine[atom_id3] += np.dot(work2[i],vel[atom_id3]) * 2.0

            self.trio_eneflo_pote[atom_id1] -= np.dot(-work1[i]-work2[i],vel[atom_id2]) + np.dot(work2[i], vel[atom_id3])
            self.trio_eneflo_pote[atom_id2] -= np.dot(work1[i],vel[atom_id1]) + np.dot(work2[i], vel[atom_id3])
            self.trio_eneflo_pote[atom_id3] -= np.dot(work1[i],vel[atom_id1]) + np.dot(-work1[i]-work2[i], vel[atom_id2])
        return 
        

    def calc_quad(self, vel, atoms, ene, work1, work2, work3):
        work_sum = -(work1++work2+work3)
        for i, quad in enumerate(atoms):
            atom_id1, atom_id2, atom_id3, atom_id4 = quad
            self.quad_ene_i[atom_id1] += ene[i]
            self.quad_ene_i[atom_id2] += ene[i]
            self.quad_ene_i[atom_id3] += ene[i]
            self.quad_ene_i[atom_id4] += ene[i]
            
            fv1 = np.dot(work1[i], vel[atom_id1])
            fv2 = np.dot(work2[i], vel[atom_id2])
            fv3 = np.dot(work_sum[i], vel[atom_id3])
            fv4 = np.dot(work3[i], vel[atom_id4])
            self.quad_eneflo_kine[atom_id1] += fv1 * 3.0
            self.quad_eneflo_kine[atom_id2] += fv2 * 3.0
            self.quad_eneflo_kine[atom_id3] += fv3 * 3.0
            self.quad_eneflo_kine[atom_id4] += fv4 * 3.0

            self.quad_eneflo_pote[atom_id1] -= fv2 + fv3 + fv4
            self.quad_eneflo_pote[atom_id2] -= fv1 + fv3 + fv4
            self.quad_eneflo_pote[atom_id3] -= fv1 + fv2 + fv4
            self.quad_eneflo_pote[atom_id4] -= fv1 + fv3 + fv2
        return

    def print_eneflo(self):
        line = ""
        line += "Eneflo.print_eneflo():\n"
        line += "Potential_i               : " + str(np.sum(self.pair_ene_i) + np.sum(self.trio_ene_i) + np.sum(self.quad_ene_i)) + "\n"
        line += "Energyflow Potential term : " + str(np.sum(self.pair_eneflo_pote) + np.sum(self.trio_eneflo_pote) + np.sum(self.quad_eneflo_pote)) + "\n"
        line += "Energyflow Kinetic term   : " + str(np.sum(self.pair_eneflo_kine) + np.sum(self.trio_eneflo_kine) + np.sum(self.quad_eneflo_kine))
        return line

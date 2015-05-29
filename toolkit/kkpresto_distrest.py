#!/usr/bin/python

import re
import kkpdb
import kkpresto

class PrestoDistRest(object):
    def __init__(self, molid1, molid2, res_id1, res_id2, 
                 res_name1, res_name2, atom_name1, atom_name2,
                 coeff_low,coeff_high, dist_low, dist_high,
                 display = "NO"):
        self.molid = (molid1, molid2)
        self.res_id = (res_id1, res_id2)
        self.res_name = (res_name1, res_name2)
        self.atom_name = (atom_name1, atom_name2)
        self.coeff = (coeff_low, coeff_high)
        self.dist = (dist_low, dist_high)
        self.display = display
        self.atom_id = ([], [])
        return 
    def get_text(self):
        line = "%4d%5d %-4s%4s%6d%5d %-4s%4s"%(self.molid[0],
                                               self.res_id[0],
                                               self.res_name[0],
                                               self.atom_name[0],
                                               self.molid[1],
                                               self.res_id[1],
                                               self.res_name[1],
                                               self.atom_name[1])
        line += "%6.2f%6.2f%9.3f%8.3f%3s"%(self.coeff[0], self.coeff[1],
                                           self.dist[0], self.dist[1], self.display)
        return line

    def set_atom_ids(self, tpl):
        aid = [-1, -1]
        for i in [0,1]:
            aid[i] = self.set_atom_id(tpl, self.molid[i],
                                      self.res_id[i],
                                      self.res_name[i],
                                      self.atom_name[i])
        self.atom_id = tuple(aid)
        return
    def set_atom_id(self, tpl,
                    molid, res_id, res_name,
                    atom_name):
        ## tpl is an instance of kkpresto.TPL
        mol_type_id = 0
        mol_id_in_type = 0
        i = 0
        for mol in tpl.mols:
            for i_mol in range(mol.mol_num):
                if i == molid: break
                mol_id_in_type += 1
                i += 1
            if i == molid: break
            mol_type_id += 1
        ## tpl.mols[mol_type_id]
        ## mol_id_in_type -th instance of the mol

        mol = tpl.mols[mol_type_id]
        for i, atom in enumerate(mol.atoms):
            if atom.res_id == res_id and \
                    atom.res_name == res_name and \
                    atom.atom_name == atom_name:
                atom_id = mol.head_atom_id + mol_id_in_type * len(mol.atoms) + i
                return atom_id
        return -1
    
class PrestoDistRestReader(kkpresto.PrestoAsciiReader):
    def __init__(self, fn):
        super(PrestoDistRestReader, self).__init__(fn)
    def read(self):
        restraints = []
        self.open()
        read_mode = 0
        line = "DUM"
        while 1:
            line = self.readline_comment()
            if not line: break
            #print line
            if line[0:len("RDDSTC> LIST")] == "RDDSTC> LIST":
                read_mode = 1                
                continue
            elif read_mode == 1:
                terms = line.strip().split()
                molid1 = int(terms[0]) - 1 
                res_id1 = int(terms[1])
                res_name1 = terms[2]
                atom_name1 = terms[3]
                molid2 = int(terms[4]) - 1
                res_id2 = int(terms[5])
                res_name2 = terms[6]
                atom_name2 = terms[7]
                force_coef_low = float(terms[8])
                force_coef_high = float(terms[9])
                dist_low = float(terms[10])
                dist_high = float(terms[11])
                display = terms[12]
                #molid1 = int(line[0:4])-1
                #res_id1 = int(line[4:9])
                #res_name1 = line[10:14].strip()
                #atom_name1 = line[14:18].strip()

                #molid2 = int(line[20:24]) -1
                #res_id2 = int(line[24:29])
                #res_name2 = line[30:34].strip()
                #atom_name2 = line[34:38].strip()
            
                #force_coef_low = float(line[38:44])
                #force_coef_high = float(line[44:50])
                #dist_low = float(line[51:59])
                #dist_high = float(line[59:67])
                #display = line[67:70]

                pdr = PrestoDistRest(molid1, molid2, res_id1, res_id2, 
                                     res_name1, res_name2, atom_name1, atom_name2,
                                     force_coef_low, force_coef_high,
                                     dist_low, dist_high,
                                     display)
            
                restraints.append(pdr)

        self.close()
        return restraints

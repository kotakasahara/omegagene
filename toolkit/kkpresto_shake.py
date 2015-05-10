#!/usr/bin/python2.7

import sys
import re
import kkkit
import kkpresto
##DEBUG = True
DEBUG = False

class PrestoShake(object):
    def __init__(self,atom_center, atom_ids, dists):
        self.atom_center = atom_center
        self.atom_ids = atom_ids
        self.dists = dists
        return

class SHKReader(kkpresto.PrestoAsciiReader):
    def __init__(self, fn):
        super(SHKReader,self).__init__(fn)
        ### shake[MOL_ID] = [PrestoShake, ... ]
        self.shake = []
        ### shake_sys[n of elements in a shake group (2 or 3 or 4)]
        ###   = [PrestoShake, ... ]
        self.shake_sys = {}
        self.shake_sys[2] = []
        self.shake_sys[3] = []
        self.shake_sys[4] = []
        self.const_pairs = set()
        
    def read_shk(self):
        self.shake = []

        self.open()
        reading_mol_id = -1
        reading_mol = ""
        line = self.readline()
        while line:
            terms = re.compile("\s+").split(line.strip())
            #print terms
            if line[0:6] == "SHAKE>":
                reading_mol_id += 1
                reading_mol = self.readline().strip()
                print "Mol " + str(reading_mol_id) + ": " + reading_mol
                self.shake.append([])
                #print "Len " + str(len(self.shake))
            elif len(terms) >= 4:
                try:
                    a = int(terms[0])
                    b = int(terms[1])
                except:
                    sys.stderr.write("Error: the first and second terms should be integer values, the number of atoms, and id fo the center atom\n")
                    sys.stderr.write(line)
                    sys.exit(0)
                atom_center = int(terms[1])
                atom_ids = []
                dists = []
                n_atoms = int(terms[0])
                
                try: assert(n_atoms >= 2 and n_atoms<=4)
                except:
                    sys.stderr.write("Error: the number of atoms in a shake unit should be 0, 2, 3, or 4.\n")
                    sys.stderr.write(line)
                    sys.exit(0)

                try:
                    if n_atoms == 2:
                        assert(len(terms) == 4)
                        atom_ids = tuple([int(terms[2])])
                        dists = tuple([float(terms[3])])
                    elif n_atoms == 3:
                        assert(len(terms) == 7)
                        atom_ids = tuple([int(terms[2]),int(terms[3])])
                        dists = tuple([float(terms[4]),float(terms[5]),float(terms[6])])
                    elif n_atoms == 4:
                        assert(len(terms) == 11)
                        atom_ids = tuple([int(terms[2]),int(terms[3]),int(terms[4])])
                        dists = tuple([float(terms[5]),float(terms[6]),float(terms[7]),
                                       float(terms[8]),float(terms[9]),float(terms[10])])
                except:
                    sys.stderr.write("Error: syntax of the shk file\n")
                    sys.stderr.write(line)
                    sys.exit(0)
                self.shake[reading_mol_id].append(PrestoShake(atom_center, atom_ids, dists))
                #print "SHAKE " + str(reading_mol_id) + " " + str(atom_center) + " ... " + ",".join([str(x) for x in atom_ids])
            elif len(terms) == 0: pass
            elif len(terms) == 1 and terms[0] == "0": pass
            else:
                sys.stderr.write("Error: syntax of the shk file\n")
                sys.stderr.write(line)
                sys.exit(0)

            line = self.readline()            
        self.close()
        return self.shake

    def expand_shake_info(self, tpl):
        """
        shake file describes shake unit for each molecule class.
        atom ids are relative ids defined in each molecule.
        the information should be expand to each instance in the system.
        
        arguments:
          tpl: an instance of TPL class in kkpresto.py
        """
        for molid, shakes in enumerate(self.shake):
            #molid = i+1
            print "expand_shake _info " + str(molid)
            tplmol = tpl.mols[molid]
            for shk in shakes:
                for mol_i in range(0, tplmol.mol_num):
                    head_atom_id = tplmol.head_atom_id + len(tplmol.atoms) * mol_i -1
                    
                    new_shk = PrestoShake(shk.atom_center + head_atom_id,
                                          [x + head_atom_id for x in shk.atom_ids],
                                          shk.dists)
                    self.shake_sys[len(shk.atom_ids)+1].append(new_shk)
                    
                    # remove constrained pairs from bonding potentials

                    #pairs = [(new_shk.atom_center, i) for i in new_shk.atom_ids]
                    #pairs.extend([(i, new_shk.atom_center) for i in new_shk.atom_ids])
                    #for i, aid1 in enumerate(new_shk.atom_ids):
                    #    for j in range(i):
                    #        aid2 = new_shk.atom_ids[j]
                    #        pairs.append((aid1, aid2))
                    #        pairs.append((aid2, aid1))
                    
                    if len(new_shk.atom_ids) == 1:
                        self.const_pairs.add((new_shk.atom_center, new_shk.atom_ids[0]))
                        self.const_pairs.add((new_shk.atom_ids[0], new_shk.atom_center))
                    elif len(new_shk.atom_ids) == 2:
                        self.const_pairs.add((new_shk.atom_center, new_shk.atom_ids[0]))
                        self.const_pairs.add((new_shk.atom_ids[0], new_shk.atom_center))
                        self.const_pairs.add((new_shk.atom_center, new_shk.atom_ids[1]))
                        self.const_pairs.add((new_shk.atom_ids[1], new_shk.atom_center))
                        self.const_pairs.add((new_shk.atom_ids[0], new_shk.atom_ids[1]))
                        self.const_pairs.add((new_shk.atom_ids[1], new_shk.atom_ids[0]))
                    elif len(new_shk.atom_ids) == 3:
                        self.const_pairs.add((new_shk.atom_center, new_shk.atom_ids[0]))
                        self.const_pairs.add((new_shk.atom_ids[0], new_shk.atom_center))
                        self.const_pairs.add((new_shk.atom_center, new_shk.atom_ids[1]))
                        self.const_pairs.add((new_shk.atom_ids[1], new_shk.atom_center))
                        self.const_pairs.add((new_shk.atom_center, new_shk.atom_ids[2]))
                        self.const_pairs.add((new_shk.atom_ids[2], new_shk.atom_center))
                        self.const_pairs.add((new_shk.atom_ids[0], new_shk.atom_ids[1]))
                        self.const_pairs.add((new_shk.atom_ids[1], new_shk.atom_ids[0]))
                        self.const_pairs.add((new_shk.atom_ids[0], new_shk.atom_ids[2]))
                        self.const_pairs.add((new_shk.atom_ids[2], new_shk.atom_ids[0]))
                        self.const_pairs.add((new_shk.atom_ids[1], new_shk.atom_ids[2]))
                        self.const_pairs.add((new_shk.atom_ids[2], new_shk.atom_ids[1]))
                    else:
                        pass

        #print "tpl.remove_bond_angle_constraint"
        #tpl.remove_bond_angle_constraints(self.const_pairs)
        return tpl

    def print_shake_info(self):
        print "- For input SHAKE info:"
        for mol_i, shk in enumerate(self.shake):
            print "-- mol: %d: "%mol_i
            print "--- the number of constraint: %d"%len(shk)
        print "- For expanded SHAKE info for the system:"
        for n_const, const in self.shake_sys.items():
            print "-- shake const: %d, num: %d"%(n_const, len(const))
        return

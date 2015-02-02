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
    def read_shk(self):
        ### shake[MOL_ID] = [PrestoShake, ... ]
        shake = []

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
                shake.append({})
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
                shake[reading_mol_id][atom_center] = PrestoShake(atom_center, atom_ids, dists)
                #print "SHAKE " + str(reading_mol_id) + " " + str(atom_center) + " ... " + ",".join([str(x) for x in atom_ids])
            elif len(terms) == 0: pass
            elif len(terms) == 1 and terms[0] == "0": pass
            else:
                sys.stderr.write("Error: syntax of the shk file\n")
                sys.stderr.write(line)
                sys.exit(0)

            line = self.readline()            
        self.close()
        return shake

#!/usr/bin/python

import re
import kkpdb
import kkpresto

class CelestePosRest(object):
    def __init__(self, atomid, crd_x, crd_y, crd_z, dist_margin, coef):
        self.atomid = atomid
        self.crd_x = crd_x
        self.crd_y = crd_y
        self.crd_z = crd_z
        self.dist_margin = dist_margin
        self.coef = coef
        return 
    def get_text(self):
        line = "%8d %10.6f %10.6f %10.6f %6.3f %10.6f"%(self.atomid,
                                                        self.crd_x, self.crd_y, self.crd_z,
                                                        self.dist_margin, self.coef)
        return line

    
class CelestePosRestReader(kkpresto.PrestoAsciiReader):
    def __init__(self, fn):
        super(CelestePosRestReader, self).__init__(fn)
    def read(self):
        restraints = []
        self.open()
        read_mode = 0
        line = "DUM"
        while 1:
            line = self.readline_comment()
            if not line: break
            #print line
            terms = line.strip().split()
            atomid = int(terms[0]) -1
            crd_x = float(terms[1])
            crd_y = float(terms[2])
            crd_z = float(terms[3])
            dist_margin = float(terms[4])
            coef = float(terms[5])
            ppr = CelestePosRest(atomid, crd_x, crd_y, crd_z, dist_margin, coef)
            restraints.append(ppr)

        self.close()
        return restraints

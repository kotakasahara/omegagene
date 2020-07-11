#!/usr/bin/python
import sys
import re
import kkpdb
import kkpresto

POSRESUNIT_NORMAL = 0
POSRESUNIT_Z = 1
POSRESUNIT_MULTIWELL01 = 2

class CelestePosRest(object):
    def __init__(self, atomid, crd_x, crd_y, crd_z, dist_margin, coef,
                 rest_type_txt, params):
        self.atomid = atomid
        self.crd_x = crd_x
        self.crd_y = crd_y
        self.crd_z = crd_z
        self.dist_margin = dist_margin
        self.coef = coef
        self.rest_type = -1
        self.n_params = len(params)
        self.params = params
        if rest_type_txt == "normal":
            self.rest_type = POSRESUNIT_NORMAL
        elif rest_type_txt == "z":
            self.rest_type = POSRESUNIT_Z
        elif rest_type_txt == "multiwell01":
            self.rest_type = POSRESUNIT_MULTIWELL01
        else:
            sys.stderr.write("Error: Unknown position restraint type : " + self.rest_type_txt + "\n")
            sys.exit(1)
        return 
    def get_text(self):
        type_txt = "normal"
        if self.rest_type == POSRESUNIT_Z:
            type_txt = "z"
        if self.rest_type == POSRESUNIT_MULTIWELL01:
            type_txt = "multiwell01"
        line = "%8d %10.6f %10.6f %10.6f %6.3f %10.6f %s"%(self.atomid,
                                                           self.crd_x, self.crd_y, self.crd_z,
                                                           self.dist_margin, self.coef, type_txt)
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
            if len(terms) < 1:
                #sys.stderr.write("skip line : " + line + "\n")
                continue
            atomid = int(terms[0]) -1
            crd_x = float(terms[1])
            crd_y = float(terms[2])
            crd_z = float(terms[3])
            dist_margin = float(terms[4])
            coef = float(terms[5])
            rest_type_txt = "normal"
            if len(terms) >= 7:
                rest_type_txt = terms[6]
            extra_params = []
            if len(terms) >= 8:
                extra_params = [ float(x) for x in terms[7:]]
                print(extra_params)
            ppr = CelestePosRest(atomid, crd_x, crd_y, crd_z, dist_margin, coef,
                                 rest_type_txt,
                                 extra_params)
            restraints.append(ppr)

        self.close()
        return restraints

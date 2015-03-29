#!/usr/bin/python2.7

import kkkit
import numpy
##DEBUG = True
DEBUG = False
import kkstruct
import struct

class PrestoRestart(object):
    def __init__(self):
        self.title = "Untitled"
        self.n_atoms = -1
        self.n_vel = -1
        self.n_loop = 0
        self.time = 0.0
        self.e_total = 0.0
        self.e_kine = 0.0
        self.e_pot = 0.0
        self.crd = numpy.array([])
        self.vel = numpy.array([])
    def set_crd_from_model(self,model):
        self.crd = []
        for atom in model.atoms:
            self.crd.append(atom.crd)
        return numpy.array(crd)

class PrestoRestartReader(kkkit.FileBI):
    def __init__(self,fn):
        super(PrestoRestartReader, self).__init__(fn)
    def open(self):
        return super(PrestoRestartReader, self).open()
    def read_restart(self):
        self.open()
        rest = PrestoRestart()
        buf = struct.unpack("@i",self.f.read(4))[0]
        rest.title = str(struct.unpack("@80s", self.f.read(80)))
        buf, buf = struct.unpack("@ii",self.f.read(8))
        rest.n_atoms, rest.n_vel = \
            struct.unpack("@ii",self.f.read(8))
        buf, buf = struct.unpack("@ii",self.f.read(8))
        rest.n_loop = struct.unpack("@i",self.f.read(4))[0]
        rest.time = struct.unpack("@d",self.f.read(8))[0]
        rest.e_total, rest.e_kine, rest.e_pot = \
            struct.unpack("@ddd",self.f.read(24))
        buf,buf = struct.unpack("@ii",self.f.read(8))
        tmp_crd = []
        for i in range(0,rest.n_atoms):
            crd = struct.unpack("@ddd",self.f.read(24))
            tmp_crd.append(numpy.array(crd))
            ## debug
            print crd
        buf,buf = struct.unpack("@ii",self.f.read(8))
        tmp_vel = []
        for i in range(0,rest.n_atoms):
            vel = struct.unpack("@ddd",self.f.read(24))
            tmp_vel.append(numpy.array(vel))
        rest.crd = numpy.array(tmp_crd)
        rest.vel = numpy.array(tmp_vel)
        self.close()
        return rest

class PrestoRestartWriter(kkkit.FileBO):
    def __init__(self,fn):
        super(PrestoRestartWriter, self).__init__(fn)
    def open(self):
        return super(PrestoRestartWriter, self).open()
    def write_restart(self, rest):
        self.open()
        #l = len(rest.title)
        #for i in range(l, 80):
        #    rest.title += " "
        #print "[" + rest.title + "]"

        self.f.write(struct.pack("@i", 80))
        self.f.write(struct.pack("@80s",rest.title))
        self.f.write(struct.pack("@i", 80))

        self.f.write(struct.pack("@i", 8))
        self.f.write(struct.pack("@ii", rest.n_atoms, rest.n_vel))
        self.f.write(struct.pack("@i", 8))

        self.f.write(struct.pack("@i",36))
        self.f.write(struct.pack("@i",rest.n_loop))
        self.f.write(struct.pack("@d",rest.time))
        self.f.write(struct.pack("@d",rest.e_total))
        self.f.write(struct.pack("@d",rest.e_kine))
        self.f.write(struct.pack("@d",rest.e_pot))
        self.f.write(struct.pack("@i",36))

        self.f.write(struct.pack("@i",rest.n_atoms*3*8))
        for crd in rest.crd:
            self.f.write(struct.pack("@d",crd[0]))
            self.f.write(struct.pack("@d",crd[1]))
            self.f.write(struct.pack("@d",crd[2]))
        self.f.write(struct.pack("@i",rest.n_atoms*3*8))
        #print rest.n_atoms*3

        self.f.write(struct.pack("@i",rest.n_vel*3*8))
        for vel in rest.vel:
            self.f.write(struct.pack("@d",vel[0]))
            self.f.write(struct.pack("@d",vel[1]))
            self.f.write(struct.pack("@d",vel[2]))
        self.f.write(struct.pack("@i",rest.n_vel*3*8))

        self.f.close()



        

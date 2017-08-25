#!/usr/bin/python2.7

import kkkit
import numpy
##DEBUG = True
DEBUG = False
import kkstruct
import struct

class PrestoCrdFrame(object):
    def __init__(self,
                 step, time, cpu_time,
                 e_total, e_kine, temperature, e_pot,
                 rmsf, vdw15, hyd15, rmsd, crds):
        self.step = step
        self.time = time
        self.cpu_time = cpu_time
        self.e_total = e_total
        self.e_kine = e_kine
        self.temperature = temperature
        self.e_pot = e_pot
        self.rmsf = rmsd
        self.vdw15 = vdw15
        self.hyd15 = hyd15
        self.rmsd = rmsd
        self.crds = crds
        return

class PrestoCrdWriter(kkkit.FileBO):
    PREC_SINGLE = 4
    PREC_DOUBLE = 8
    def __init__(self,fn):
        super(PrestoCrdWriter, self).__init__(fn)
    def write_frame(self, in_crd, step, time,
                    bin=False, bin_rev=False,
                    ignore=set(),
                    cpu_time=0.0,
                    total_e=0.0,
                    kinetic_e=0.0,
                    temperature=0.0,
                    potential_e=0.0,
                    rmsf=0.0,
                    vdw=0.0,
                    hyd=0.0,
                    rmsd=0.0):
        self.f.write(struct.pack("@i", 44))
        self.f.write(struct.pack("@i", step)) 
        self.f.write(struct.pack("@f", time)) 
        self.f.write(struct.pack("@f", cpu_time)) #cpu_time
        self.f.write(struct.pack("@f", total_e)) #total E
        self.f.write(struct.pack("@f", kinetic_e)) #kinetic E
        self.f.write(struct.pack("@f", temperature)) #temperature
        self.f.write(struct.pack("@f", potential_e)) #potential E
        self.f.write(struct.pack("@f", rmsf)) #rmsf
        self.f.write(struct.pack("@f", vdw)) #vdw
        self.f.write(struct.pack("@f", hyd)) #hyd
        self.f.write(struct.pack("@f", rmsd)) #rmsd
        self.f.write(struct.pack("@i", 44))
        n_atoms = len(in_crd)
        self.f.write(struct.pack("@i", n_atoms*3*4))
        if not bin:
            for i, x in enumerate(in_crd):
                self.f.write(struct.pack("@fff", x[0], x[1], x[2]))
        else:
            for i, x in enumerate(in_crd):
                self.f.write(x[0]+ x[1]+ x[2])
        self.f.write(struct.pack("@i", n_atoms*3*4))
        return 

class PrestoCrdReader(kkkit.FileBI):
    def __init__(self,fn):
        super(PrestoCrdReader, self).__init__(fn)
        self.n_atoms = 0
        self.n_frames = 0
        self.size_frame = 0
    def read_next_frame(self, bin=False):
        frame = []
        try:
            if bin:
                frame = self.read_frame_binary()
            else:
                frame = self.read_frame()
        except:
            return []
        return frame
    def read_information(self):
        self.open()
        frame_cod = self.read_next_frame()
        self.f.seek(52)
        buf_size = struct.unpack("@i",self.f.read(4))[0] ## n_atoms * 3 * 4
        self.n_atoms = buf_size / 12
        self.size_frame = 52 + buf_size + 8 
        #print self.size_frame
        self.close()
        return
    def count_n_frames(self):
        self.open()
        
        self.n_frames = 0
        flg=True
        while flg:
            flg=False
            try:
                #print self.f.tell()
                self.f.seek(self.size_frame-4, 1)
                buf_size = struct.unpack("@i",self.f.read(4))[0]
                self.n_frames += 1
                flg=True
            except:
                break
        self.close()
        return self.n_frames
    def read_frame(self):
        buf_size = struct.unpack("@i",self.f.read(4))[0] ## 44
        step = struct.unpack("@i",self.f.read(4))[0]
        time = struct.unpack("@f",self.f.read(4))[0]

        cpu_time = struct.unpack("@f",self.f.read(4))[0]
        e_total = struct.unpack("@f",self.f.read(4))[0]
        e_kine = struct.unpack("@f",self.f.read(4))[0]
        temperature = struct.unpack("@f",self.f.read(4))[0]
        e_pot = struct.unpack("@f",self.f.read(4))[0]
        rmsf = struct.unpack("@f",self.f.read(4))[0]
        vdw15 = struct.unpack("@i",self.f.read(4))[0]
        hyd15 = struct.unpack("@i",self.f.read(4))[0]
        rmsd = struct.unpack("@f",self.f.read(4))[0]

        buf_size = struct.unpack("@i",self.f.read(4))[0] ## 44
        
        buf_size = struct.unpack("@i",self.f.read(4))[0] ## n_atoms * 3 * 4
        self.n_atoms = buf_size / 12
        #print n_atoms
        crds = []
        for i in range(0,self.n_atoms):
            #print i*12
            crds.append(numpy.array(struct.unpack("@fff",self.f.read(12))))
        buf_size = struct.unpack("@i",self.f.read(4))[0] ## n_atoms * 3
        frame = PrestoCrdFrame(step, time, cpu_time,
                               e_total, e_kine, temperature, e_pot,
                               rmsf, vdw15, hyd15, rmsd, crds)
        return frame
    def read_frame_binary(self):
        buf_size = struct.unpack("@i",self.f.read(4))[0] ## 44
        step = struct.unpack("@i",self.f.read(4))[0]
        time = struct.unpack("@f",self.f.read(4))[0]

        cpu_time = struct.unpack("@f",self.f.read(4))[0]
        e_total = struct.unpack("@f",self.f.read(4))[0]
        e_kine = struct.unpack("@f",self.f.read(4))[0]
        temperature = struct.unpack("@f",self.f.read(4))[0]
        e_pot = struct.unpack("@f",self.f.read(4))[0]
        rmsf = struct.unpack("@f",self.f.read(4))[0]
        vdw15 = struct.unpack("@i",self.f.read(4))[0]
        hyd15 = struct.unpack("@i",self.f.read(4))[0]
        rmsd = struct.unpack("@f",self.f.read(4))[0]

        buf_size = struct.unpack("@i",self.f.read(4))[0] ## 44
        
        buf_size = struct.unpack("@i",self.f.read(4))[0] ## n_atoms * 3
        n_atoms = buf_size / 12

        crds = []
        for i in range(0,n_atoms):
            crds.append([self.f.read(4), self.f.read(4), self.f.read(4)])
        buf_size = self.f.read(4)
        frame = PrestoCrdFrame(step, time, cpu_time,
                               e_total, e_kine, temperature, e_pot,
                               rmsf, vdw15, hyd15, rmsd, crds)

        return frame

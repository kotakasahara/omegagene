#!/usr/bin/python2.7

import kkkit
import numpy as np
##DEBUG = True
DEBUG = False
import kkstruct
import struct as st
import copy
import sys

class GroTrrFrame(object):
    def __init__(self,
                 step, time, box, crds, vels, frcs):
        self.step = step
        self.time = time
        self.box = box
        self.crds = crds
        self.vels = vels
        self.frcs = frcs
        return

class GroTrrReader(kkkit.FileBI):
    PREC_SINGLE = 4
    PREC_DOUBLE = 8
    def __init__(self, fn):
        super(GroTrrReader, self).__init__(fn)
        self.code_real = "f"
        self.size_real = 4
        self.endian = ">" #big endian

        return
    def open(self):
        super(GroTrrReader, self).open()
        self.set_endian()
        self.read_information()
        return 
    def read_frame(self):
        buf = st.unpack("@i")
        return 
    def set_endian(self):
        #Checking endian
        self.f.seek(0)
        buf = self.f.read(4)
        #print buf
        #print len(buf)
        #print str(st.unpack("<i",buf)) + ":"+str(st.unpack(">i",buf))
        if st.unpack("<i",buf)[0] == 1993:
            self.endian = "<"
            print("Little endian")
        elif st.unpack(">i",buf)[0] != 1993:
            sys.stderr.write("Error: The first 4 byte is not the integer 1993")
        else:
            self.endian = ">"            
            print("Big endian")
        return  self.endian
    def set_real(self, flg_double):
        if flg_double:
            self.code_real = "d"
            self.size_real = 8
        else:
            self.code_real = "f"
            self.size_real = 4
    def read_int(self):
        """
        unpack a value from binary string
        """
        read_string = self.endian + "i"
        i = st.unpack(read_string, self.f.read(4))[0]
        return i
    def read_real(self):
        """
        unpack a value from binary string
        """
        read_string = self.endian + self.code_real
        i = st.unpack(read_string, self.f.read(self.size_real))[0]
        return i
    #def write_int(self, val):
    #    write_string = self.endian + "i"
    #    buf = st.pack(write_string, val)
    #    self.fo.write(buf)
    #    return buf
    def get_bin_real(self, val):
        write_string = self.endian + self.code_real
        bin = st.pack(write_string, val)
        return bin
    def write_real(self, val):
        buf = get_bin_real(val)
        self.fo.write(buf)
        return buf
    def read_information(self):
        """
        read basic information of this trajectory file
        self.box_size
        self.x_size
        self.v_size
        self.f_size
        self.n_atoms
        self.n_frames
        """
        self.f.seek(0)
        self.f.read(32)
        self.size_box = self.read_int() #32-36
        self.f.read(16)                 #36-52
        self.size_x = self.read_int()   #52-56
        self.size_v = self.read_int()   #56-60
        self.size_f = self.read_int()   #60-64
        self.n_atoms = self.read_int()  #64-68
        self.f.read(8) #68-76
        
        if self.size_x == self.n_atoms * 3 * 8:
            self.set_real(True)
            self.size_header = 92
        elif self.size_x == self.n_atoms * 3 * 4:
            self.set_real(False)
            self.size_header = 84
        else:
            print "size_x is not n_atoms*3*4 or *8"
            sys.exit()

        self.read_real() #time
        self.read_real() #lambda

        self.size_frame = self.size_header + self.size_box + self.size_x + self.size_v + self.size_f
        
        ### ?????????????  2013.11.22
        ####if self.code_real=="d": self.size_frame -= 8
        ### ?????????????  2013.11.22
        self.f.seek(0)
        return
    def count_n_frames(self):
        self.f.seek(0)
        self.n_frames = -1
        buf = "dummy"
        while buf != "":
            self.n_frames += 1
            self.f.seek(self.n_frames * self.size_frame)
            buf = self.f.read(4)
        self.f.seek(0)
        #print "n_frames: " + str(self.n_frames)
        #print "n_atoms: " + str(self.n_atoms)
        #print "size_x: " + str(self.size_x)
        #print "size_v: " + str(self.size_v)
        #print "size_f: " + str(self.size_f)
        #print "size_frame: " + str(self.size_frame)
        return
    def read_next_frame(self):
        frame = []
        try:
            frame = self.read_frame()
        except:
            return []
        return frame
    def read_nth_frame(self, frame_n):
        self.f.seek(self.size_frame * frame_n)
        try:
            frame = self.read_frame()
        except:
            return None
        return frame
    def read_frame(self):
        buf   = self.f.read(68)
        step = self.read_int()
        buf   = self.f.read(4)
        time = self.read_real()
        buf   = self.read_real()
        
        #box
        box = []

        if self.size_box != 0:
            for i in range(9):
                box.append(self.read_real() * 10.0)
        crds = []
        if self.size_x != 0:
            for i in range(self.n_atoms):
                crds.append([self.read_real(),
                             self.read_real(),
                             self.read_real()])
        vels = []
        if self.size_v != 0:
            for i in range(self.n_atoms):
                vels.append([self.read_real(),
                             self.read_real(),
                             self.read_real()])
        frcs = []
        if self.size_f != 0:
            for i in range(self.n_atoms):
                frcs.append([self.read_real(),
                             self.read_real(),
                             self.read_real()])
        return GroTrrFrame(step, time,
                           np.array(box),
                           np.array(crds),
                           np.array(vels),
                           np.array(frcs))
    #def write_header(self):
        #0-4: magic_number 1993
        #4-8: 8:double or 4:float
        #8-12: n_atoms
        #12-16: n_frames
    #    self.write_int(1993)
    #    self.write_int(self.size_real)
    #    self.write_int(self.n_atoms)
    #    self.write_int(self.n_frames)
    #    return
    #def write_atom(self,crd_x,crd_y,crd_z):
    #    self.fo.write(crd_x)
    #    self.fo.write(crd_y)
    #    self.fo.write(crd_z)
    #    return 
    
class GroTrrWriter(kkkit.FileBO):
    PREC_SINGLE = 4
    PREC_DOUBLE = 8
    def __init__(self, fn):
        super(GroTrrWriter, self).__init__(fn)
        return
    def write_frame(self, in_box, in_crd, in_vel, in_force, step, time,
                    bin=False, bin_rev=True, prec=PREC_SINGLE, ignore=set(),
                    scale_length = 1.0):
        ## bin = True
        ##   receive values as bit strings and 
        ##   directory write down the strings
        
        ## bin_rev = True
        ##   reversing bit strings
        ##   that means changing the endian
        
        box = None
        crd = None
        vel = None
        force = None
        if not bin and in_box != None:
            box = np.array(in_box) * scale_length
        if not bin and in_crd != None:
            crd = np.array(in_crd) * scale_length
        if not bin and in_vel != None:
            vel = np.array(in_vel) * scale_length

        box_size = 0
        if box != None:
            box_size = 9 * prec # 4byte * 9dim
        x_size = 0
        if crd != None:
            x_size = (len(crd) - len(ignore)) * 3 * prec # 4byte * 3dim

        
        
        v_size = 0
        if vel != None:
            v_size = (len(vel) - len(ignore)) * 3 * prec# 4byte * 3dim
        f_size = 0
        if force != None:
            f_size = (len(force) - len(ignore)) * 3 * prec  # 4byte * 3dim

        prec_esc = "f"
        prec_esc3 = "fff"
        if prec == GroTrrWriter.PREC_DOUBLE:
            prec_esc = "d"
            prec_esc3 = "ddd"
        
        self.f.write(st.pack(">i",1993))
        self.f.write(st.pack(">i",13))
        self.f.write(st.pack(">i",12))
        self.f.write(st.pack(">12s",'GMX_trn_file'))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">i",box_size))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">i",x_size))
        self.f.write(st.pack(">i",v_size))
        self.f.write(st.pack(">i",f_size))
        self.f.write(st.pack(">i",len(crd)-len(ignore)))
        self.f.write(st.pack(">i",step))
        self.f.write(st.pack(">i",0))
        self.f.write(st.pack(">"+prec_esc, time))
        self.f.write(st.pack(">"+prec_esc, 0.0))
        if box_size != 0:
            self.f.write(st.pack(">"+prec_esc3, box[0], 0.0,    0.0))
            self.f.write(st.pack(">"+prec_esc3, 0.0,    box[1], 0.0))
            self.f.write(st.pack(">"+prec_esc3, 0.0,    0.0,    box[2]))
        if x_size != 0:
            if bin:
                for i,x in enumerate(crd):
                    if i in ignore: continue
                    if bin_rev:
                        x[0] = x[0][::-1]
                        x[1] = x[1][::-1]
                        x[2] = x[2][::-1]
                    self.f.write(x[0]+x[1]+x[2])
                    #print x[0] #DEBUG
                    #print st.unpack(">f",x[0]) #DEBUG
            else:
                for i, x in enumerate(crd):
                    if i in ignore: continue
		    self.f.write(st.pack(">"+prec_esc3, x[0],x[1],x[2]))
        if v_size != 0:
            if bin:
                for  i,x in enumerate(vel):
                    if i in ignore: continue
                    if bin_rev:
                        x[0] = x[0][::-1]
                        x[1] = x[1][::-1]
                        x[2] = x[2][::-1]
                    self.f.write(x[0]+x[1]+x[2])
            else:
                for i,x in enumerate(vel):
                    if i in ignore: continue
                    self.f.write(st.pack(">"+prec_esc3, x[0],x[1],x[2]))
        if f_size != 0:
            if bin:
                for i,x in enumerate(force):
                    if i in ignore: continue
                    if bin_rev:
                        x[0] = x[0][::-1]
                        x[1] = x[1][::-1]
                        x[2] = x[2][::-1]
                    self.f.write(x[0]+x[1]+x[2])
            else:
                for i,x in enumerate(force):
                    if i in ignore: continue
                    self.f.write(st.pack(">"+prec_esc3, x[0],x[1],x[2]))
        return 0        
        

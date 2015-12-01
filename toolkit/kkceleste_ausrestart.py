#!/usr/bin/python2.7

import kkkit
import numpy
##DEBUG = True
DEBUG = False
import struct as st

class CelesteAUSRestart(object):

    def __init__(self):
        self. AUS_TYPE = {"dummy":0, "type1":1, "type2":2,
                          "dist-mass-centers":3
                          }
        self.header = ""
        self.aus_type_str = "dummy"
        self.aus_type = 0
        self.n_enhance_groups = 0
        self.enhance_groups = []
        self.n_atoms_in_groups = []
        self.crd_groups = []
        return
    def set_aus_type(self, in_type):
        if not in_type in self.AUS_TYPE:
            print "Error : Invalid aus_type = " + in_type
            sys.exit(1)
        self.aus_type_str = in_type
        self.aus_type = self.AUS_TYPE[self.aus_type_str]
        
    def dump_group_coord(self):
        buf = ""
        buf += st.pack("@i", len(self.header)+1)
        buf += self.header+"\0"
        #print "dbg1130 " + self.header
        buf += st.pack("@i", self.aus_type)
        buf += st.pack("@i", self.n_enhance_groups)
        for grp in self.enhance_groups:
            buf += st.pack("@i", grp)
            #print "dbg1130 " + str(grp)
        for n_atoms in self.n_atoms_in_groups:
            buf += st.pack("@i", n_atoms)
        for grp in self.crd_groups:
            for atom in grp:
                buf += st.pack("@ddd", atom[0], atom[1], atom[2])
        return buf
    def generate_aus_restart(self, restart, atom_groups, atom_group_names,
                             aus_group_names):
        self.header="V-AUS"
        self.enhance_groups = []
        for i_grp, name in enumerate(aus_group_names):
            grp_id = atom_group_names.index(name) 
            self.enhance_groups.append(grp_id)
            self.n_atoms_in_groups.append(len(atom_groups[grp_id]))
            #print "dbg1130 " + name + " " + str(len(atom_groups[name]))
        self.crd_groups = []
        for name in aus_group_names:
            #print "dbg1130 generate_aus_restart " + str(name)
            #print name
            crd_group = []
            grp_id = atom_group_names.index(name) 
            for atom_id in atom_groups[grp_id]:
                atom = restart.crd[atom_id]
                crd_group.append(atom)
            self.crd_groups.append(crd_group)
            print "Group : " + name + " group-id: " + str(grp_id) \
                + " enhance-grp: " + str(len(self.crd_groups)-1) \
                + " n_atoms: " + str(len(crd_group)) + " atoms."
        self.n_enhance_groups = len(self.crd_groups)
        return 

class CelesteAUSRestartReader(kkkit.FileBI):
    def __init__(self, fn):
        super(CelesteAUSRestartReader, self).__init__(fn)
    def read_aus_restart(self, atom_groups, atom_group_names):
        print "Read AUS restart file: "
        self.open()
        self.crd_groups = {}
        rest = CelesteAUSRestart()
        buf = st.unpack("@i",self.f.read(4))[0]
        rest.header = self.f.read(buf)
        rest.aus_type = st.unpack("@i",self.f.read(4))[0]
        rest.n_enhance_groups = st.unpack("@i",self.f.read(4))[0]
        for i in range(rest.n_enhance_groups):
            rest.enhance_groups.append(st.unpack("@i",self.f.read(4))[0])
        for i in range(rest.n_enhance_groups):
            rest.n_atoms_in_groups.append(st.unpack("@i",self.f.read(4))[0])
        for enhance_grp_id, grp_id in enumerate(rest.enhance_groups):
            #print "dbg1130 read_aus_restart " + str(enhance_grp_id) + \
            #    " grp_id: " + str(grp_id) + " atom: " + str(rest.n_atoms_in_groups[enhance_grp_id]) \
            #     + " " + str(len(atom_groups[grp_id]))
            crd_group = []
            for i_atm in range(rest.n_atoms_in_groups[enhance_grp_id]):
                buf_crd = st.unpack("@ddd",self.f.read(8*3))
                crd_group.append(buf_crd)
                #print buf_crd
            name = atom_group_names[grp_id]
            rest.crd_groups.append(crd_group)
            print "Group : " + name + " " + str(len(crd_group)) + " atoms."
        self.close()
        return rest
        
        
        

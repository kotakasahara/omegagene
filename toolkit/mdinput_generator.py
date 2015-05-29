#!/usr/bin/python2.7

################################
##  Celeste input
##    version 
##    15032221
################################
MAGIC=66261
#VERSION = 13111501
#VERSION = 14013101 
#VERSION = 14013204  ## PBC origin, shake
#VERSION = 15020801  ## shake, expand_shake_info()
#VERSION = 15030901  ## vmcmd 
#VERSION = 15032221  ## atom_groups
VERSION = "v.0.34.b" ## version_info

import sys
from optparse import OptionParser
import numpy as np
import struct as st

import kkpresto as prst
import kkpresto_restart as prstrst
import kkpdb as pdb
import kkmmsystem as mmsys
import kkmmff as mmff
import kkmmconfig 
import kkpresto_shake as shk
import kkmm_expand
import define_atom_groups as atgrp
import kkpresto_distrest as disres

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_config',
                 help="file name for config file")
    p.add_option('-o', dest='fn_out',
                 help="file name for config file")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def _main():
    opts, args = get_options()
    mdinputgen = MDInputGen(opts.fn_config, opts.fn_out)
    print "read_files()"
    mdinputgen.read_files()
    print "dump_mdinput()"    
    mdinputgen.dump_mdinput()
    return 

class MDInputGen(object):
    def __init__(self, in_fn_config, in_fn_out):
        self.fn_config = in_fn_config
        self.fn_out = in_fn_out
        self.config = None
        self.system = None
        self.structure = None
        self.restart = None
        self.tpl = None
        self.settle = None
        self.expand = None
        self.atom_groups = {}
        self.dist_rest = []
        return 

    def read_files(self):
        print "read_config"
        self.config = kkmmconfig.ConfigReader(self.fn_config).read_config()
        print "read_tpl"
        print self.config.get_val("fn-i-tpl")
        self.tpl = prst.TPLReader(self.config.get_val("fn-i-tpl")).read_tpl()
        self.tpl.enumerate_12_13_14()
        print "read initial pdb"
        self.structure = pdb.PDBReader(self.config.get_val("fn-i-initial-pdb")).read_model()
        print "prepare system"
        self.system = mmsys.MmSystem(self.structure,
                                     self.config.get_val("cell-x"),
                                     self.config.get_val("cell-y"),
                                     self.config.get_val("cell-z"))
        if self.config.get_val("cell-center-x") and \
                self.config.get_val("cell-center-y") and \
                self.config.get_val("cell-center-z"):
            self.system.pbc.set_center(np.array([self.config.get_val("cell-center-x"),
                                                 self.config.get_val("cell-center-y"),
                                                 self.config.get_val("cell-center-z")]))
        if self.config.get_val("cell-origin-x") and self.config.get_val("cell-origin-y") \
                and self.config.get_val("cell-origin-z"):
            self.system.pbc.origin = np.array([self.config.get_val("cell-origin-x"),
                                          self.config.get_val("cell-origin-y"),
                                          self.config.get_val("cell-origin-z")])

            
        self.system.set_atom_info_from_tpl(self.tpl)
        print "read restart"
        self.restart = prstrst.PrestoRestartReader(self.config.get_val("fn-i-restart")).read_restart()
        print "set_crd_vel_from_restart"
        self.system.set_crd_vel_from_restart(self.restart)

    ## zd self energy
        self.system.ff.set_params(self.config.get_val("cutoff"))
    #system.ff.set_zd_params(config.get_val("ele-ewaldalpha"),
    #                        system.charge,
    #                        tpl.atom_id_12,
    #                        tpl.atom_id_13,
    #                        tpl.atom_id_14nb)
    #system.store_self_energy(system.ff.energy_self)
        self.read_shake()
        self.read_settle()
        self.read_expand()

        if self.config.get_val("fn-i-atom-groups"):
            atom_groups_reader = atgrp.AtomGroupsReader(self.config.get_val("fn-i-atom-groups"))
            print self.config.get_val("fn-i-atom-groups")
            print atom_groups_reader.fn
            self.atom_groups = atom_groups_reader.read_groups()

        if self.config.get_val("fn-i-dist-restraint"):
            dist_rest_reader = disres.PrestoDistRestReader(self.config.get_val("fn-i-dist-restraint"))
            print self.config.get_val("fn-i-dist-restraint")
            self.dist_rest = dist_rest_reader.read()
            for d in self.dist_rest:
                d.set_atom_ids(self.tpl)

        return

    def read_shake(self):
        if self.config.get_val("fn-i-shake"):
            shkreader = shk.SHKReader(self.config.get_val("fn-i-shake"))
            mol_settle = set(self.config.get_val("mol-settle"))
            if not mol_settle:  mol_settle = set()
            shkreader.read_shk(exclude=mol_settle)
            self.tpl = shkreader.expand_shake_info(self.tpl)
            self.system.shake = shkreader.shake_sys
            shkreader.print_shake_info()
        return

    def read_settle(self):
        fn_shk = self.config.get_val("fn-i-shake")
        mol_settle = set(self.config.get_val("mol-settle"))
        if fn_shk and mol_settle:
            shkreader = shk.SHKReader(fn_shk)
            if not mol_settle:  mol_settle = set()
            shkreader.read_shk(readonly=mol_settle)
            self.tpl = shkreader.expand_shake_info(self.tpl)
            self.settle = shkreader.shake_sys
            shkreader.print_shake_info()
        return
        
    def read_expand(self):
        self.expand = None
        if self.config.get_val("fn-i-ttp-v-mcmd-inp"):
            self.expand = kkmm_expand.ExpandConf()
            self.expand.read_mcmdparams(self.config.get_val("fn-i-ttp-v-mcmd-inp"))
            if self.config.get_val("fn-i-ttp-v-mcmd-initial"):
                self.expand.read_init(self.config.get_val("fn-i-ttp-v-mcmd-initial"))
                if self.config.get_val("ttp-v-mcmd-initial-vs"):
                    print "Definition conflicts:"
                    print "The options \"--fn-i-ttp-v-mcmd-initial\" and \"--ttp-v-mcmd-initial-vs\" are mutually exclusive."
                    sys.exit(1)
            elif self.config.get_val("ttp-v-mcmd-initial-vs") and \
                    self.config.get_val("ttp-v-mcmd-seed"):
                self.expand.init_vs = self.config.get_val("ttp-v-mcmd-initial-vs")
                self.expand.seed = self.config.get_val("ttp-v-mcmd-seed")
            else:
                print "For mcmd, --ttp-v-mcmd-initial or the two options --ttp-v-mcmd-initial-vs and --ttp-v-mcmd-seed are required."
                sys.exit(1)
        return
    
    def dump_mdinput(self):
        f = open(self.fn_out, "wb")
        ## Magic number 66261
        f.write(st.pack("@i", MAGIC))
        ## Version
        f.write(st.pack("@i", len(VERSION)))
        f.write(VERSION)
        ## Config
        #buf_config = dump_mmconfig(cfg)
        #f.write(struct.pack("@i", len(buf_config)))
        #f.write(buf_config)
        buf_box = self.dump_box(self.system)
        buf_coordinates = self.dump_crdvel(self.system.crd)
        buf_velocities = self.dump_crdvel(self.system.vel)
        buf_topol = self.dump_topol(self.system, self.tpl)
        buf_shake = ""
        buf_settle = ""
        buf_expand = ""
        if self.system.shake:
            buf_shake = self.dump_shake(self.system.model, self.tpl,
                                        self.system.shake)
        if self.settle:
            buf_settle = self.dump_shake(self.system.model, self.tpl,
                                         self.settle)
        if self.expand:
            buf_expand = self.dump_expand(self.expand)

        buf_atom_groups = self.dump_atom_groups(self.atom_groups)

        buf_dist_rest = self.dump_dist_rest(self.dist_rest)

        #if config.get_val("particle-cluster-shake"):
        f.write(st.pack("@i", len(buf_box)))
        f.write(st.pack("@i", len(buf_coordinates)))
        f.write(st.pack("@i", len(buf_velocities)))
        f.write(st.pack("@i", len(buf_topol)))
        f.write(st.pack("@i", len(buf_shake)))
        f.write(st.pack("@i", len(buf_settle)))
        f.write(st.pack("@i", len(buf_expand)))
        f.write(st.pack("@i", len(buf_atom_groups)))
        f.write(st.pack("@i", len(buf_dist_rest)))
        print "size: buf_box        : " + str(len(buf_box))
        print "size: buf_coordinates: " + str(len(buf_coordinates))
        print "size: buf_velocities : " + str(len(buf_velocities))
        print "size: buf_topol      : " + str(len(buf_topol))
        print "size: buf_shake      : " + str(len(buf_shake))
        print "size: buf_settle     : " + str(len(buf_settle))
        print "size: buf_expand     : " + str(len(buf_expand))
        print "size: buf_atom_groups: " + str(len(buf_atom_groups))
        print "size: buf_dist_rest: " + str(len(buf_dist_rest))

        #f.write(st.pack("@i", len(buf_pcluster)))
        f.write(buf_box)
        f.write(buf_coordinates)
        f.write(buf_velocities)
        f.write(buf_topol)
        f.write(buf_shake)
        f.write(buf_settle)
        f.write(buf_expand)
    #f.write(buf_pcluster)
        f.write(buf_atom_groups)
        f.write(buf_dist_rest)
        f.close()
        return

    def dump_box(self,system):
        buf = ""
        buf += st.pack("@ddd", system.pbc.L[0],0.0,0.0)
        buf += st.pack("@ddd", 0.0, system.pbc.L[1],0.0)
        buf += st.pack("@ddd", 0.0, 0.0, system.pbc.L[2])
        buf += st.pack("@ddd", system.pbc.origin[0],
                       system.pbc.origin[1],
                       system.pbc.origin[2])
        return buf

    def dump_crdvel(self, crd):
        buf = ""
        buf += st.pack("@i", len(crd))
        for atom_crd in crd:
            buf += st.pack("@ddd", atom_crd[0], atom_crd[1], atom_crd[2])
        return buf

    def pack_14(self, params_14):
        buf14 = ""
        buf14 += st.pack("@i", len(params_14))

        for params in params_14:
            atom_id1 = params[0][0]
            atom_id2 = params[0][1]
            atom_id3 = params[0][2]
            atom_id4 = params[0][3]
            energy = params[1][0]
            overlaps = params[1][1]
            symmetry = params[1][2]
            phase = params[1][3]
            flag_14nb = params[1][4]
            buf14 += st.pack("@iiii", atom_id1, atom_id2, atom_id3, atom_id4)
            buf14 += st.pack("@d", energy)
            buf14 += st.pack("@ii", overlaps, symmetry)
            buf14 += st.pack("@d", phase)
            buf14 += st.pack("@i", flag_14nb)
        return buf14

    def dump_topol(self, system, tpl):
        buf = ""
        buf += st.pack("@i", system.n_atoms)
        for val in system.charge:
            buf += st.pack("@d", val)
        for val in system.mass:
            buf += st.pack("@d", val)
        for val in system.atomtype:
            buf += st.pack("@i", val)

        buf_nbpair = ""
        nb_types = set()
        for type_pair in tpl.nb_pair.keys():
            nb_types.add(type_pair[0])
            nb_types.add(type_pair[1])
        buf_nbpair += st.pack("@i", len(nb_types))
        buf_nbpair += st.pack("@i", len(tpl.nb_pair.keys()))
        for type_pair, params in tpl.nb_pair.items():
            buf_nbpair += st.pack("@ii", type_pair[0], type_pair[1]) ## atom_type1, 2
            buf_nbpair += st.pack("@dd", params[0], params[1])
        ## parameters for 6-powered term, 12-powered term

        buf_12 = ""
        buf_12 += st.pack("@i", len(tpl.atom_id_12))
        print "# bonds : "+ str(len(tpl.atom_id_12))
        for params in tpl.atom_id_12:
            atom_id1 = params[0][0]
            atom_id2 = params[0][1]
            epsiron = params[1][0]
            r0 = params[1][1]
            buf_12 += st.pack("@ii", atom_id1, atom_id2)
            buf_12 += st.pack("@dd", epsiron, r0)

        buf_13 = ""
        buf_13 += st.pack("@i", len(tpl.atom_id_13))
        print "# angles : "+ str(len(tpl.atom_id_13))
        for params in tpl.atom_id_13:
            atom_id1 = params[0][0]
            atom_id2 = params[0][1]
            atom_id3 = params[0][2]
            epsiron = params[1][0]
            theta0 = params[1][1]
            buf_13 += st.pack("@iii", atom_id1, atom_id2, atom_id3)
            buf_13 += st.pack("@dd", epsiron, theta0)

        buf_14 = self.pack_14(tpl.atom_id_14)
        print "# torsions : "+ str(len(tpl.atom_id_14))

        buf_14imp = self.pack_14(tpl.atom_id_14_imp)
        print "# improper : "+ str(len(tpl.atom_id_14_imp))

        buf_14nb = ""
        buf_14nb += st.pack("@i", len(tpl.atom_id_14nb))
        for params in tpl.atom_id_14nb:
            atom_id1 = params[0][0]
            atom_id2 = params[0][1]
            atom_type1 = params[1][0]
            atom_type2 = params[1][1]
            coeff_vdw = params[1][2]
            coeff_ele = params[1][3]
            buf_14nb += st.pack("@iiii", atom_id1, atom_id2, atom_type1, atom_type2)
            buf_14nb += st.pack("@dd", coeff_vdw, coeff_ele)

        buf_non15 = ""
        buf_non15 += st.pack("@i", len(tpl.atom_pair_non15)*2)
        for params in tpl.atom_pair_non15:
            atom_id1 = params[0]
            atom_id2 = params[1]
            buf_non15 += st.pack("@ii", atom_id1, atom_id2)
            buf_non15 += st.pack("@ii", atom_id2, atom_id1)

        buf += st.pack("@i", len(buf_nbpair))
        buf += buf_nbpair

        buf += st.pack("@i", len(buf_12))
        buf += buf_12
        buf += st.pack("@i", len(buf_13))
        buf += buf_13
        buf += st.pack("@i", len(buf_14))
        buf += buf_14
        buf += st.pack("@i", len(buf_14imp))
        buf += buf_14imp
        buf += st.pack("@i", len(buf_14nb))
        buf += buf_14nb
        buf += st.pack("@i", len(buf_non15))
        buf += buf_non15
        return buf

    def dump_mmconfig(cfg):
        buf = ""
        keys = ["n-steps",
                "dt",
                "cutoff",
                "electrostatic",
                "ele-cutoff",            
                "ele-ewaldapha", 
                "integrator",
                "print-interval-coord",
                "print-interval-velo",
                "print-interval-log",
                "print-interval-energy",
                "print-interval-energyflow",
                "fn-o-coord", "fn-o-log",
                "fn-o-energy", "fn-o-energyflow",
                "thermostat", "temperature",
                "barostat", "pressure",
                "center-of-motion"]
        buf = ""
        for key in keys:
            if cfg.type_def[key] == kkmmconfig.ConfigReader.STR:
                buf += st.pack("@i", len(cfg.get_val[key]))
                buf += st.pack("@s", cfg.get_val[key])
            elif cfg.type_def[key] == kkmmconfig.ConfigReader.INT:
                buf += st.pack("@i", cfg.get_val[key])
            elif cfg.type_def[key] == kkmmconfig.ConfigReader.FLOAT:
                buf += st.pack("@d", cfg.get_val[key])
        return buf

    def convert_shake_info(self, model, tpl, shake):
        shk_atoms_index = []
        shk_atoms = []
        shk_distances_index = []
        shk_distances = []
        
        buf = ""
        mol_id = 0
        mol_num = 0
        mol_atom_id = 0
        atom_id_head = 0
        tmp_one_atom_units = []
        atom_id = -1
        for atom in model.atoms:
            atom_id += 1
            mol_atom_id += 1
            if mol_atom_id == len(tpl.mols[mol_id].atoms)+1:
                mol_atom_id = 1
                mol_num += 1
                atom_id_head = atom_id
                if tpl.mols[mol_id].mol_num == mol_num:
                    mol_num = 0
                    mol_id += 1
            if shk and mol_id < len(shk) and mol_atom_id in shk[mol_id]:
                shk_atoms_index.append(len(shk_atoms))
                shk_atoms.append(atom_id)
                for c_shk_atom_id in shk[mol_id][mol_atom_id].atom_ids:
                    atom_id_global = atom_id_head + c_shk_atom_id - 1
                    shk_atom.append(atom_id_global)
                shk_distances_index.append(len(shk_distances))
                for c_shk_dist in shk.dists:
                    shk_distances.append(c_shk_dist)
        shk_atoms_index.append(len(shk_atoms))
        shk_distances_index.append(len(shk_distances))
        return shk_atoms_index, shk_atoms, shk_distances_index, shk_distances
    
    def dump_shake(self, model, tpl, shake_sys):
    #shk_atoms_index, shk_atoms, \
    #    shk_distances_index, shk_distances \
    #    = convert_shake_info(model, tpl, shake)
    
    #buf = ""
    #buf += st.pack("@i", len(shk_atoms_index))
    #buf += st.pack("@i", len(shk_atoms))
    #buf += st.pack("@i", len(shk_distance))    
    #for idx in shk_atoms_index:
    #    buf += st.pack("@i", idx)
    #for atom_id in shk_atoms:
    #    buf += st.pack("@i", atom_id)
    #for idx in shk_distances_index:
    #    buf += st.pack("@d", dist)
    #for dist in shk_distances:
    #    buf += st.pack("@d", dist)

        buf = ""
        for shk_size in [2,3,4]:
            buf += st.pack("@i", len(shake_sys[shk_size]))
        for shk_size in [2,3,4]:
            for shk_elem in shake_sys[shk_size]:
                buf += st.pack("@i", shk_elem.atom_center) #
                for dest_id in shk_elem.atom_ids:
                    buf += st.pack("@i", dest_id)
                for dist in shk_elem.dists:
                    buf += st.pack("@d", dist) # (N * N - N)/2 
                
        return buf

    def dump_expand(self, expand):
        buf = ""
        n_vs = len(expand.vmcmd_range.keys())
        buf += st.pack("@i", n_vs)
        #print "DBG2: "+ str(expand.temperature)
        buf += st.pack("@i", expand.interval)
        buf += st.pack("@d", expand.temperature)
        for i in range(1,n_vs+1):
            buf += st.pack("@i", len(expand.vmcmd_params[i])-3)
            #print "Expand: " + str(i)
            #print expand.vmcmd_range[i]
            buf += st.pack("@dddd",
                           expand.vmcmd_range[i][0], # lambda_min
                           expand.vmcmd_range[i][1], # lambda_max
                           expand.vmcmd_range[i][2], # prob
                           expand.vmcmd_range[i][3]) # prob
        #for i in range(1,n_vs+1):

            for prm in expand.vmcmd_params[i]:
                buf += st.pack("@d", float(prm))
        buf += st.pack("@ii", expand.init_vs, expand.seed)
        return buf
    def dump_atom_groups(self, atom_groups):
        buf = ""
        buf += st.pack("@i", len(atom_groups.keys()))
        for name, atoms in atom_groups.items():
            buf += st.pack("@i", len(name))
            buf += name
            buf += st.pack("@i", len(atoms))
        for name, atoms in atom_groups.items():
            for atomid in atoms:
                buf += st.pack("@i", atomid)
        return buf
    def dump_dist_rest(self, dist_rest):
        buf = ""
        buf += st.pack("@i", len(dist_rest))
        for dr in dist_rest:
            buf += st.pack("@i", dr.atom_id[0])
            buf += st.pack("@i", dr.atom_id[1])
            buf += st.pack("@f", dr.coeff[0])
            buf += st.pack("@f", dr.coeff[1])
            buf += st.pack("@f", dr.dist[0])
            buf += st.pack("@f", dr.dist[1])
        return buf
if __name__ == "__main__":
    _main()

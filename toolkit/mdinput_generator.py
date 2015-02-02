#!/usr/bin/python2.7

################################
##  myMD input
##    version 
##    14013204
################################
MAGIC=66261
#VERSION = 13111501
#VERSION = 14013101 
VERSION = 14013204  ## PBC origin, shake

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
    print "read_config"
    config = kkmmconfig.ConfigReader(opts.fn_config).read_config()
    print "read_tpl"
    tpl = prst.TPLReader(config.get_val("fn-i-tpl")).read_tpl()
    tpl.enumerate_12_13_14()
    print "read initial pdb"
    structure = pdb.PDBReader(config.get_val("fn-i-initial-pdb")).read_model()
    print "prepare system"
    system = mmsys.MmSystem(structure,
                            config.get_val("cell-x"),
                            config.get_val("cell-y"),
                            config.get_val("cell-z"))
    if config.get_val("cell-center-x") and config.get_val("cell-center-y") \
            and config.get_val("cell-center-z"):
        system.pbc.set_center(np.array([config.get_val("cell-center-x"),
                                        config.get_val("cell-center-y"),
                                        config.get_val("cell-center-z")]))
    if config.get_val("cell-origin-x") and config.get_val("cell-origin-y") \
            and config.get_val("cell-origin-z"):
        system.pbc.origin = np.array([config.get_val("cell-origin-x"),
                                      config.get_val("cell-origin-y"),
                                      config.get_val("cell-origin-z")])


    system.set_atom_info_from_tpl(tpl)
    print "read restart"
    restart = prstrst.PrestoRestartReader(config.get_val("fn-i-restart")).read_restart()
    print "set_crd_vel_from_restart"
    system.set_crd_vel_from_restart(restart)

    ## zd self energy
    system.ff.set_params(config.get_val("cutoff"))
    #system.ff.set_zd_params(config.get_val("ele-ewaldalpha"),
    #                        system.charge,
    #                        tpl.atom_id_12,
    #                        tpl.atom_id_13,
    #                        tpl.atom_id_14nb)
    #system.store_self_energy(system.ff.energy_self)
    
    if config.get_val("fn-i-shake"):
        system.shake = shk.SHKReader(config.get_val("fn-i-shake")).read_shk()
    dump_mdinput(system, tpl, opts.fn_out, config)
    return

def dump_mdinput(system, tpl, fn_out, config):
    f = open(fn_out, "wb")
    ## Magic number 66261
    f.write(st.pack("@i", MAGIC))
    ## Version
    f.write(st.pack("@i", VERSION))
    ## Config
    #buf_config = dump_mmconfig(cfg)
    #f.write(struct.pack("@i", len(buf_config)))
    #f.write(buf_config)
    buf_box = dump_box(system)
    buf_coordinates = dump_crdvel(system.crd)
    buf_velocities = dump_crdvel(system.vel)
    buf_topol = dump_topol(system, tpl)
    buf_shake = ""
    if system.shake:
        buf_shake = dump_shake(system.model, tpl,
                               system.shake)
    #if config.get_val("particle-cluster-shake"):
        #buf_pcluster = dump_pcluster_from_shake(system.model, tpl, system.shake)

    f.write(st.pack("@i", len(buf_box)))
    f.write(st.pack("@i", len(buf_coordinates)))
    f.write(st.pack("@i", len(buf_velocities)))
    f.write(st.pack("@i", len(buf_topol)))
    f.write(st.pack("@i", len(buf_shake)))
    #f.write(st.pack("@i", len(buf_pcluster)))

    f.write(buf_box)
    f.write(buf_coordinates)
    f.write(buf_velocities)
    f.write(buf_topol)
    f.write(buf_shake)
    #f.write(buf_pcluster)

    f.close()
    return

def dump_box(system):
    buf = ""
    buf += st.pack("@ddd", system.pbc.L[0],0.0,0.0)
    buf += st.pack("@ddd", 0.0, system.pbc.L[1],0.0)
    buf += st.pack("@ddd", 0.0, 0.0, system.pbc.L[2])
    buf += st.pack("@ddd", system.pbc.origin[0],
                   system.pbc.origin[1],system.pbc.origin[2])
    return buf

def dump_crdvel(crd):
    buf = ""
    buf += st.pack("@i", len(crd))
    for atom_crd in crd:
        buf += st.pack("@ddd", atom_crd[0], atom_crd[1], atom_crd[2])
    return buf

def dump_topol(system, tpl):
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
    for params in tpl.atom_id_12:
        atom_id1 = params[0][0]
        atom_id2 = params[0][1]
        epsiron = params[1][0]
        r0 = params[1][1]
        buf_12 += st.pack("@ii", atom_id1, atom_id2)
        buf_12 += st.pack("@dd", epsiron, r0)

    buf_13 = ""
    buf_13 += st.pack("@i", len(tpl.atom_id_13))
    for params in tpl.atom_id_13:
        atom_id1 = params[0][0]
        atom_id2 = params[0][1]
        atom_id3 = params[0][2]
        epsiron = params[1][0]
        theta0 = params[1][1]
        buf_13 += st.pack("@iii", atom_id1, atom_id2, atom_id3)
        buf_13 += st.pack("@dd", epsiron, theta0)

    def pack_14(params_14):
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
    buf_14 = pack_14(tpl.atom_id_14)
    buf_14imp = pack_14(tpl.atom_id_14_imp)

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

def convert_shake_info(model, tpl, shake):
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
        aotm_id += 1
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
    
def dump_shake(model, tpl, shake):
    shk_atoms_index, shk_atoms, \
        shk_distances_index, shk_distances \
        = convert_shake_info(model, tpl, shake)
    
    buf = ""
    buf += st.pack("@i", len(shk_atoms_index))
    buf += st.pack("@i", len(shk_atoms))
    buf += st.pack("@i", len(shk_distance))    
    for idx in shk_atoms_index:
        buf += st.pack("@i", idx)
    for atom_id in shk_atoms:
        buf += st.pack("@i", atom_id)
    for idx in shk_distances_index:
        buf += st.pack("@d", dist)
    for dist in shk_distances:
        buf += st.pack("@d", dist)
    return buf

def dump_pcluster_from_shake(model, tpl, shake):
    shk_atoms_index, shk_atoms, \
        shk_distances_index, shk_distances \
        = convert_shake_info(model, tpl, shake)
    
    buf = ""
    buf += st.pack("@i", len(shk_atoms_index))
    buf += st.pack("@i", len(shk_atoms))
    for idx in shk_atoms_index:
        buf += st.pack("@i", idx)
    for atom_id in shk_atoms:
        buf += st.pack("@i", atom_id)
    return buf

if __name__ == "__main__":
    _main()

#!/usr/bin/python2.6

from optparse import OptionParser
import numpy
import sys
import kkpdb
import kkpresto_crd as kkcrd
import kkpresto_restart as kkrst
import re
import copy
import kkatomgroup

def get_options():
    p = OptionParser()
    p.add_option('--i-cod', dest='fn_i_cod',
                 help="Input trajectory file")
    p.add_option('--o-cod', dest='fn_o_cod',
                 help="Output trajectory file")
    #p.add_option('--i-log', dest='fn_log',
    #             help="Input log file")
    p.add_option('--i-group', dest='fn_group',
                 help="Input group file")
    p.add_option('--group-name', dest='group_name',
                 help="Group name to be output")
    p.add_option('--i-restart', dest='fn_restart',
                 help="Input restart file")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args



def _main():
    opts, args = get_options()
    
    ignore_atoms = set()
    grps = None
    grp_naem = None
    if opts.fn_group:
        agreader = kkatomgroup.AtomGroupsReader(opts.fn_group)
        grps, grp_names = agreader.read_groups()
    grp_id = grp_names.index(opts.group_name)
    print str(grp_id) + " : " + opts.group_name +" " + str(len(grps))
    
    crd_writer = kkcrd.PrestoCrdWriter(opts.fn_o_cod)    
    crd_writer.open()

    crd_reader = kkcrd.PrestoCrdReader(opts.fn_i_cod)    
    crd_reader.open()
    read_frame = 0
    frame_crd = crd_reader.read_next_frame(bin=True)
    while frame_crd:
        crd_writer.write_frame_with_frame(frame_crd, bin=True)
        frame_crd = crd_reader.read_next_frame(bin=True)
        read_frame += 1
        
    crd_reader.close()    

    rst_reader = kkrst.PrestoRestartReader(opts.fn_restart)
    rst = rst_reader.read_restart()

    for i in range(rst.n_atoms):
        if not i in grps[grp_id]:
            ignore_atoms.add(i)

    crd_writer.write_frame(rst.crd,
                           rst.n_loop, rst.time,
                           bin = False, bin_rev=False,
                           ignore=ignore_atoms,
                           total_e = rst.e_total,
                           kinetic_e = rst.e_kine,
                           potential_e = rst.e_pot)

    crd_writer.close()
    return


if __name__ == "__main__":
    _main()

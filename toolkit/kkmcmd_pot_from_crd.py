#!/usr/bin/python2.7
import os
import numpy as np
from optparse import OptionParser
import re
import struct as st
import sys
import kkpresto_crd
import math

def get_options():
    p = OptionParser()
    
    p.add_option('--cod', dest='fn_cod',
                 help="File name of the trajectoory file")
    p.add_option('-o', dest='fn_out',
                 default="prob.out",
                 help="Output file name for the probabilities")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = get_options()    
    n_frames = 0
    print "couting the number of frames:"
    cod_reader = kkpresto_crd.PrestoCrdReader(opts.fn_cod)
    cod_reader.read_information()
    n_frames += cod_reader.count_n_frames()
    print opts.fn_cod + " : " + str(cod_reader.n_frames) + " " + str(n_frames)
    print str(n_frames)

    cod_reader.open()

    frame_cod = cod_reader.read_next_frame()
    fo = open(opts.fn_out,"w")
    while frame_cod:
        pot = frame_cod.e_pot
        fo.write(str(pot)+"\n")
        frame_cod = cod_reader.read_next_frame()        
    cod_reader.close()
    

    fo.close()

    return

if __name__ == "__main__":
    _main()

    

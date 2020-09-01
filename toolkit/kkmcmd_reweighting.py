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
    
    p.add_option('--flg-pot', dest='flg_pot',
                 action="store_true",
                 help="flag for reading text files with potential energies")
    p.add_option('-i', dest='fn_list',
                 help="File name of the trajectory file list")
    p.add_option('--i-cano', dest='fn_cano_dist',
                 help="File name of the canonical distribution")
    p.add_option('-o', dest='fn_prob',
                 default="prob.out",
                 help="Output file name for the probabilities")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args
    

def read_filelist(fn_list):
    files = []
    f = open(fn_list)
    for line in f:
        files.append(line.strip().split()[0])
    f.close()
    return files

def read_pdf(fn_pdf):
    pdf = []
    f = open(fn_pdf)
    for line in f:
        terms = line.strip().split()
        pdf.append((float(terms[0]), float(terms[1])))
    f.close()
    return pdf

def reweight_pot(files, pdf):
    n_frames = 0
    for fn_pot in files:
        f=open(fn_pot)
        for line in f:
            n_frames += 1
        f.close()
    print str(n_frames)

    prob_trj = np.zeros(n_frames, dtype=np.float64)
    
    read_frames = 0
    for fn_pot in files:
        print fn_pot
        f=open(fn_pot)
        for line in f:
            pot = float(line.strip().split()[0])
            prob = 0.0
            for i_bin in range(len(pdf)-1):
                if pot >= pdf[i_bin][0] and pot < pdf[i_bin+1][0]:
                    prob = math.exp(max(pdf[i_bin][1], pdf[i_bin+1][1]))
                    break
            prob_trj[read_frames] = prob
            
            read_frames += 1
        f.close()
        
    prob_trj /= prob_trj.sum()    

    return prob_trj
    
def reweight_codlist(files, pdf):
    n_frames = 0
    print "couting the number of frames:"
    for fn_cod in files:
        cod_reader = kkpresto_crd.PrestoCrdReader(fn_cod)
        cod_reader.read_information()
        n_frames += cod_reader.count_n_frames()
        print fn_cod + " : " + str(cod_reader.n_frames) + " " + str(n_frames)
    print str(n_frames)

    prob_trj = np.zeros(n_frames, dtype=np.float64)

    print "Reweighting"

    read_frame = 0
    for fn_cod in files:
        print fn_cod
        cod_reader = kkpresto_crd.PrestoCrdReader(fn_cod)
        cod_reader.read_information()
        cod_reader.open()
        frame_cod = cod_reader.read_next_frame()
        while frame_cod:
            
            pot = frame_cod.e_pot
            prob = 0.0
            for i_bin in range(len(pdf)-1):
                if pot >= pdf[i_bin][0] and pot < pdf[i_bin+1][0]:
                    prob = math.exp(max(pdf[i_bin][1], pdf[i_bin+1][1]))
                    break
            prob_trj[read_frame] = prob
            
            read_frame += 1
            frame_cod = cod_reader.read_next_frame()
        cod_reader.close()
    
    print "Sum probability : " + str(prob_trj.sum())

    prob_trj /= prob_trj.sum()

    return prob_trj


def _main():
    opts, args = get_options()    
    pdf = read_pdf(opts.fn_cano_dist)
    files = read_filelist(opts.fn_list)

    prob_trj = None

    if opts.flg_pot:
        prob_trj = reweight_pot(files, pdf)
    else:
        prob_trj = reweight_codlist(files, pdf)

    f_out = open(opts.fn_prob, "w")
    for p in prob_trj:
        f_out.write("%16.10e\n"%p)
    f_out.close()

    return

if __name__ == "__main__":
    _main()

    

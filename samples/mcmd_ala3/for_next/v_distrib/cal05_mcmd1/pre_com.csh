#!/bin/csh


set n_stages = $1
set n_series = $2
set margin = 0.0

set traj_weight_coef = 0
set binwidth = 50.0

python ~/local/kktools/mdtools/kkmcmd_weighting_traj.py \
  -i dat_some \
  -o inp_c1_all \
  -c ${traj_weight_coef} -d 10000 \
  -b ${binwidth} -n ${n_series} 
#  --celeste-bin

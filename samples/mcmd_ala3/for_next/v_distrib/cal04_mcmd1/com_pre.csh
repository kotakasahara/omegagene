#!/bin/csh

set n_stages = $1
set n_series = $2
set margin = 0.0

set traj_weight_coef = 0
set binwidth = 50.0

cd dat_some
bash ./link_do.bash ${n_stages} ${n_series}
cd ..

python ~/local/kktools/mdtools/kkmcmd_weighting_traj.py \
  -i dat_some \
  -o inp_c1_all \
  -c ${traj_weight_coef} -d 500000 \
  -b ${binwidth} -n ${n_series} #\
#  --celeste-bin

python ~/local/kktools/mdtools/fornext_vdistrib_mainv.py \
  --ene-margin ${margin} \
  --i-ttpv ttp_v_mcmd.inp \
  > log.txt
#  --celeste-bin \

R --vanilla --slave --args 7 v_distrib.png < r_vdistrib.R

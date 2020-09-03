#!/bin/bash

n_stages=$1
n_series=$2
margin=0.0

traj_weight_coef=0
binwidth=$3

rm -r stg_*

stg=1

## KKTOOLS=/home/usr8/14IAW655/local/kktools

while [ $stg -le ${n_stages} ]
do
  mkdir stg_${stg}
  cp -r dat_some stg_${stg}
  cd stg_${stg}/dat_some
  bash ./link_do.bash ${stg} ${n_series}
  cd ..
  mkdir v_pdf
  python ../kkmcmd_weighting_traj.py \
    -i dat_some \
    -o inp_c1_all \
    -c ${traj_weight_coef} -d 10000 \
      --celeste-bin \
    -b ${binwidth} -n ${n_series} #\


  python ../fornext_vdistrib_mainv.py \
    --ene-margin ${margin} \
    --i-ttpv ../../../../cal04_mcmd1/${stg}/n1/ttp_v_mcmd.inp \
      --celeste-bin \
   > log.txt
#    --i-ttpv ../../../../cal10_mcmd1/ttp_v_mcmd.inp \



  stg=`expr ${stg} + 1`
  cd ..
done

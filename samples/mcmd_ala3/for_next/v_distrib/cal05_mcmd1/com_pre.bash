#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=1::
#$ -l s_core=1

. /etc/profile.d/modules.sh
module load intel cuda openmpi python-extension/2.7 r

#source ~/.zshrc
n_stages=$1
n_series=$2
margin=0.0
KKTOOLS=/gs/hs0/hp170025/local/kktools
OMEGATK=${HOME}/local/og0/toolkit
traj_weight_coef=0
binwidth=5.0

head -1 ../../current_situation > aaa
mdn=`more aaa`
rm aaa
head -2 ../../current_situation > aaa
tail -1 aaa > bbb
nvs=`more bbb`
rm aaa bbb
mdn_str=`printf "%02d" ${mdn}`
cd dat_some
bash ./link_do.bash ${n_stages} ${n_series}
cd ..

python2.7 ${OMEGATK}/kkmcmd_weighting_traj.py \
  -i dat_some \
  -o inp_c1_all \
  -c ${traj_weight_coef} -d 10000 \
  -b ${binwidth} -n ${n_series} \
  --celeste-bin

python2.7 ${OMEGATK}/fornext_vdistrib_mainv.py \
  --ene-margin ${margin} \
  --i-ttpv ../../../cal${mdn_str}_mcmd1/ttp_v_mcmd.inp \
   --celeste-bin \
  > log.txt

R --vanilla --slave --args ${nvs} v_distrib.png < r_vdistrib.R

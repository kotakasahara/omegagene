#!/bin/csh

set n_stages = $1
set n_series = $2
set margin = 0.0

set traj_weight_coef = 0
set binwidth = 50.0

head -1 ../../current_situation > aaa
set mdn = `more aaa`
rm aaa
head -2 ../../current_situation > aaa
tail -1 aaa > bbb
set nvs = `more bbb`
rm aaa bbb
if ( $mdn < 10 ) then
  set mdn_str = "0${mdn}"
else
 set mdn_str = "${mdn}"
endif

cd dat_some
bash ./link_do.bash ${n_stages} ${n_series}
cd ..

python ~/local/kktools/mdtools/kkmcmd_weighting_traj.py \
  -i dat_some \
  -o inp_c1_all \
  -c ${traj_weight_coef} -d 50000 \
  -b ${binwidth} -n ${n_series} \
  --celeste-bin

python ~/local/kktools/mdtools/fornext_vdistrib_mainv.py \
  --ene-margin ${margin} \
  --i-ttpv ../../../cal${mdn_str}_mcmd1/ttp_v_mcmd.inp \
  --celeste-bin \
  > log.txt


# csh ./2_com_all

R --vanilla --slave --args ${nvs} v_distrib.png < r_vdistrib.R

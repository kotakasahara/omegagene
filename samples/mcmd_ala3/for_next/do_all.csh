#!/bin/csh

set mdnum = $1
set vstnum = $2
set nstages = $3
set nseries = $4
set iord_min = 3
set iord_max = 9

if ( $mdnum < 10 ) then
  set mdnum_str = "0${mdnum}"
else
  set mdnum_str = "${mdnum}"
endif

set dircal = "cal${mdnum_str}_mcmd1"
set dirout = "figures/${dircal}"

echo ${mdnum} > current_md_num
echo ${vstnum} > current_vst_num

rm -rf $dirout
mkdir $dirout

cd ../${dircal}
R --vanilla --slave --args ${nstages} ${nseries} ../for_next/${dirout}/energy_timecourse.png < ../for_next/r_ene.R

cd ../for_next/v_distrib

cp -r temp ${dircal}
cd ${dircal}/dat_some
csh ./link_do ${nstages} ${nseries}
cd ../
csh ./2_com_all 1

R --vanilla --slave --args ${vstnum} ../../${dirout}/v_distrib.png < r_vdistrib.R

cd ../../fit_pmc_entire

csh ./123_com ${iord_min} ${iord_max}

##cd ${dircal}
R --vanilla --slave --args ${mdnum_str} ${vstnum} ${iord_min} ${iord_max} ../${dirout}/fitpmc.png < r_fitpmc.R

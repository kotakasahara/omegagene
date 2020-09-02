#!/bin/csh

set iord = $1
set fit_surroundings = $2
set next_pre = $3

  head -1 current_situation > aaa
  set mdnum = `more aaa`
  rm aaa
  head -2 current_situation > aaa
  tail -1 aaa > bbb
  set vsnum = `more bbb`
  rm aaa bbb

if ( $mdnum < 10 ) then
  set mdnum_str = "0${mdnum}"
else
  set mdnum_str = "${mdnum}"
endif


echo "mdnum = $mdnum"
echo "vsnum = $vsnum"
echo "iord = $iord"
echo "fit_surroundings = ${fit_surroundings}"
if ( $fit_surroundings == 0 ) then
  echo "The surrounding virtual states are NOT considered for fitting"
else
  echo "The surrounding virtual states are considered for fitting"
endif

# echo "Ok? [y/others]"
# set ok = $<
# if( $ok != "y" ) then
# exit
# endif

# echo "This is the next run to the pre_run? [y/others]"
# set next_pre = $<

set dircal = "cal${mdnum_str}_mcmd1"

cp -r fit_pmc_entire/${dircal}_${iord}/*  fit_pmc_entire/${dircal}/
echo "cp -r fit_pmc_entire/${dircal}_${iord}/*  fit_pmc_entire/${dircal}/"
## rm -r fit_pmc_entire/${dircal}_*

cd fit_mix
csh ./1_sinple_mix 1 ${next_pre}
if ($fit_surroundings == 1) then
  echo "2_prep_refit"
  csh ./2_prep_refit 1
  csh ./3_refit 1
  csh ./4_re_refit 1
endif

echo "R"
cd ${dircal}
R --vanilla --slave --args ${vsnum}  fitmix.png < ../r_fitmix.R

python2.7 ../gen_ttpinput_tmpl.py -i ../../../${dircal}/ttp_v_mcmd.inp -o ttp_v_mcmd.inp -s 10 -t 300 -n ${vsnum} --range-template

@ mdnum_next = $mdnum + 1
echo ${mdnum_next}

if ( $mdnum_next < 10 ) then
  set mdnum_next_str = "0${mdnum_next}"
else
  set mdnum_next_str = "${mdnum_next}"
endif

set dircal_next = "cal${mdnum_next_str}_mcmd1"
echo ${dircal_next}

mkdir ../../../${dircal_next}
cp ttp_v_mcmd.inp ../../../${dircal_next}
cp ../../../${dircal}/md.inp  ../../../${dircal_next}
cp ../../../${dircal}/md.inp.run  ../../../${dircal_next}

echo "diff"
cd ../../
diff fit_pmc_entire/${dircal}_${iord}/s1.derv.poly  fit_pmc_entire/${dircal}/

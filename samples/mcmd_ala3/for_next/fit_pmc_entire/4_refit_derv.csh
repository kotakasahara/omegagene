#!/bin/sch

################################################# 
set iord = $1
set auto = $2

  head -1 ../current_situation > aaa
  set mdn = `more aaa`
  rm aaa
  head -2 ../current_situation > aaa
  tail -1 aaa > bbb
  set nst = `more bbb`
  rm aaa bbb
if ( $mdn < 10 ) then
  set mdn_str = "0${mdn}"
else
  set mdn_str = "${mdn}"
endif

  echo " "
  echo " ********************************"
  echo " * Deribatives are fit gurther. *"
  echo " ********************************"
  echo " "
  echo " MD done = " $mdn
  echo " N. of v-states = " $nst
  echo " "
  set parms = "7.0 "
  set aaa = ` more bin.size `
  echo "  Parameters: " $parms  $aaa   "       <------------------ OK?"
  echo " "

   set irep = 5
   echo "  repetition times are " $irep
   echo "     If you increase it, change this script. <---- NOTE "
   echo "     Very larfe repetition is wrong. "

  if( $auto == 1 ) then
    set ok = 1
  else
   echo " "
   echo "  Input 1, if OK. "
   echo " "
    set ok = $<
  endif

  if( $ok == 1 ) then
    echo " OK. I proceed. "
    goto p1 
  endif

  echo " I stop here. "
  echo " "
  exit


################################################# 
p1:

echo " ---------- "

set dir = cal${mdn_str}_mcmd1
##set dir = "a"
# Repeating.
  @ mrep = 1
  while ( $mrep <= $irep ) 

    echo "  "
    echo "iteration time = " $mrep

#  Prepare input file for the first iteration.
    if( $mrep == 1 ) then
      @ kk = 1
      cp nul refit_rerv.dat
      while ( $kk <= $nst )
        cat  refit_rerv.dat  ${dir}/s${kk}.derv.dat  > aa
        mv aa refit_rerv.dat
        echo "    Integrating derivatives of v-state: " $kk
        @ kk ++
      end
      mv   refit_rerv.dat   ${dir}/refit_rerv.dat
    endif

#  Prepare input file for further iterations.
    if( $mrep > 1 ) then
      @ kk = 1
      cp nul  refit_rerv.dat 
      while ( $kk <= $nst )
        cat  refit_rerv.dat  ${dir}/s${kk}_d.derv.dat  > aa
        mv aa  refit_rerv.dat 
        echo "    Integrating derivatives of v-state: " $kk
        @ kk ++
      end
      mv refit_rerv.dat  ${dir}/refit_rerv.dat
    endif

    set f1 = ${dir}/refit_rerv.dat 
    echo "   "

# Fitting.
    @ ii = 1
  echo " Order of polynomial ? "
  if( $auto != 1 || $iord < 1 ) then
    set iord = $<
  endif
  echo $iord
    while( $ii <= $nst )

#      @ iord = 9
#      if( $ii == 1 ) @ iord = 9
#      if( $ii == 2 ) @ iord = 9
#      if( $ii == 3 ) @ iord = 9
#      if( $ii == 4 ) @ iord = 9
#      if( $ii == 5 ) @ iord = 9
#      if( $ii == 6 ) @ iord = 9
#      if( $ii == 7 ) @ iord = 9

      echo "    " $ii " order = " $iord

      echo $parms  $aaa  "  999.0  0.0 " $iord " 0.0" > aaa1

      set f2 = ` head -1 ${dir}/range.b${ii} `
      set f3 = ` head -1 ../fit_pmc_entire/${dir}/range.s${ii} `
      echo $f2 $f3  > a1
      echo $f3 > a2

      cat aaa1 a1 a2 $f1 > ${dir}/inp_derv_s${ii}_d.dat
      rm a1 a2

      
      src_pre/aho.exe < ${dir}/inp_derv_s${ii}_d.dat
#  Overwrite.
      mv fort.11 ${dir}/dummy 
      mv fort.12 ${dir}/s${ii}_d.derv.dat   
      mv fort.20 ${dir}/s${ii}_d.derv.poly

      rm fort.13 fort.16

      @ ii ++
    end

    rm aaa1

    @ mrep ++
  end
################################################# 
  rm *.tar
  cp -r ${dir} ${dir}_f
  tar cvf  ${dir}_f.tar  ${dir}_f
  rm -r ${dir}_f

  echo " "
  echo "  end of all "
  echo " "

  exit
#######################################################################

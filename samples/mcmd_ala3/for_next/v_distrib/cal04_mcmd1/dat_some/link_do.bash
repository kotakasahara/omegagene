n_stages=$1
n_series=$2


at=`more ../../../current_md_num`
nvst=`more ../../../current_vst_num`

echo "md No. = " $at
echo "N of v states = " $nvst

if [ $at  -lt 10 ]; then
  at="0${at}"
fi

caldir="cal${at}_mcmd1"
echo $caldir

kk=1
jj=1
filenum=1
rm e*
rm v*


##echo "100.0   1000\n-100000.0\nEND" > tmp.dat
echo " 50.00   611" > tmp.dat
echo "-110100.00" >> tmp.dat
echo "END" >> tmp.dat

totalnum=${n_series}
echo $totalnum
echo $totalnum >> tmp.dat

while [ ${kk} -le ${n_stages} ]
do
  while [ ${jj} -le ${n_series} ]
  do
   ##echo "ln -s ../../../../${caldir}/${kk}/n${jj}/mult.ene e${filenum}
   ln -s ../../../../${caldir}/${kk}/n${jj}/mult.ene e${filenum}
   ln -s ../../../../${caldir}/${kk}/n${jj}/ttp_v_mcmd.out  v${filenum}
   echo "${filenum}   50000  990000000  1 " >> tmp.dat
   jj=`expr $jj + 1`
   filenum=`expr $filenum + 1`
  done
  echo ${kk}_$jj
  jj=1
  kk=`expr $kk + 1`
done

echo "END" >> tmp.dat
mv tmp.dat ../inp_c1_all
cd ..
ls -1 dat_some/e* > fil.all
ls -1 dat_some/v* > filv.all


#!/bin/bash
path=`pwd`
echo $path
STAGE=${path##*/}

SERIES=1
N_SERIES=10

#OMEGABIN=${HOME}/local/omegagene/target/bin/omegagene_wons
#OMEGATK=${HOME}/local/omegagene/toolkit

while [ $SERIES -le $N_SERIES ]
do 
    cd n$SERIES
    python ${OMEGATK}/mdinput_generator.py -i md.inp -o md.inp.cls -v v.0.52 > log_inputgen.txt
    ${OMEGABIN} --cfg md.inp.run --inp md.inp.cls > md.out 
    SERIES=$(( $SERIES + 1 ))
    cd ..
    sleep 0.5s
done

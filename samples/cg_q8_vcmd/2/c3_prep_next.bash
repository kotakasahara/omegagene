#!/bin/bash

path=`pwd`
echo $path
STAGE=${path##*/}

PREV_STAGE=$(( $STAGE -1 ))
echo ${PREV_STAGE}

#OMEGABIN=${HOME}/local/omegagene/target/bin/omegagene_wons
#OMEGATK=${HOME}/local/omegagene/toolkit

ls -1 n*/vcmd_q_raw.txt  > qraw_list.txt
python ${OMEGATK}/merge_vcmd_qraw.py -i ../${PREV_STAGE}/vcmd_next.inp -o vcmd_next.inp  --p-count 0   --i-qraw-list qraw_list.txt --o-qraw vcmd_next_qraw.dat


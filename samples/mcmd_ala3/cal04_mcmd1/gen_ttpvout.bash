#!/bin/bash

stg=0
ser=0

for ser in `seq 1 10 `
do
    for stg in `seq 1 10`
    do
	echo ${stg}/n${ser}/ttp_v_mcmd.out
	echo "0      1" > ${stg}/n${ser}/ttp_v_mcmd.out
    done
done

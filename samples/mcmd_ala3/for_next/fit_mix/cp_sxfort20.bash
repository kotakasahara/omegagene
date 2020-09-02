#!/bin/bash

vs=$1

while [ $vs -ne 0 ]
do
    vsnext=`expr $vs + 1`
    echo "cp s${vs}.fort.20 s${vsnext}.fort.20"
    cp s${vs}.fort.20 s${vsnext}.fort.20
    vs=`expr $vs - 1`
done

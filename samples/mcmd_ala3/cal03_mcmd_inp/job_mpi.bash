#!/bin/bash

cd ${CAL_DIR}
echo "CAL_DIR"
echo ${CAL_DIR}

echo ${JOB_SCRIPT}

echo "R1"
echo ${R1}
echo "R2"
echo ${R2}
echo "R3"
echo ${R3}

if [ -n "${R1}" ]
then
  cd ${R1}
  bash ${JOB_SCRIPT} 0 &
fi

if [ -n "${R2}" ]
then
  cd ${R2}
  bash ${JOB_SCRIPT} 1 &
fi

if [ -n "${R3}" ]
then
  cd ${R3}
  bash ${JOB_SCRIPT} 2 &
fi
wait

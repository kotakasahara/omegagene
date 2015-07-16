#!/bin/bash
#$ -cwd
#$ -l short

~/bin/celeste --inp trpc.cls --cfg md_input.cfg > log_cpu.txt

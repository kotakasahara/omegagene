#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l gpu -l month
#$ -l s_vmem=4G -l mem_req=4G

export LD_LIBRARY_PATH=/usr/local/cuda-6.0/lib64:$LD_LIBRARY_PATH

ulimit -H -v 4194304

~/bin/celeste_gpu_035s6 --inp trpc.cls --cfg md_input.cfg   > log_gpu.txt

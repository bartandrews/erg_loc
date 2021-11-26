#!/bin/bash
#$ -l h_rt=01:00:00
#$ -l h_vmem=32G
#$ -cwd
#$ -o stdout.$JOB_ID
#$ -j y
#$ -pe shared 16
#$ -l highp

python code/W_flow.py -mod heisenberg -L 8 -Nup 4 -pauli 0 -bc o -dis 11000 -W_min 0.5 -W_max 12.5 -W_samp 24

#!/bin/bash
#$ -l h_rt=3:00:00:00
#$ -l h_data=16G
#$ -l h_vmem=512G
#$ -cwd
#$ -o stdout.$JOB_ID
#$ -j y
#$ -pe shared 32
#$ -l highp

source /u/home/b/baandr12/.bash_profile
source /u/home/b/baandr12/.bashrc

python code/L_flow.py -mod heisenberg -thr 32 -L_min 12 -L_max 18 -L_samp 4 -bc o -dis 100 -W 0;
python code/L_flow.py -mod heisenberg -thr 32 -L_min 12 -L_max 18 -L_samp 4 -bc o -dis 100 -W 2.5;
python code/L_flow.py -mod heisenberg -thr 32 -L_min 12 -L_max 18 -L_samp 4 -bc o -dis 100 -W 5

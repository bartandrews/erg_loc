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

python code/T_flow.py -mod ponte2015 -thr 32 -L 8 -h0 2 -T0 7 -dis 1000 -bc o -W 0.5 -T_min 0 -T_max 3 -T_samp 31;
python code/T_flow.py -mod ponte2015 -thr 32 -L 8 -h0 2 -T0 7 -dis 1000 -bc o -W 8 -T_min 0 -T_max 3 -T_samp 31

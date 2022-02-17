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

python code/delta_flow.py -mod spin2021 -L 12 -bc o -T0 1 -T1 1 -W 2 -delta_min 0 -delta_max 1 -delta_samp 11 -dis 100;
python code/delta_flow.py -mod spin2021 -L 14 -bc o -T0 1 -T1 1 -W 2 -delta_min 0.4 -delta_max 0.4 -delta_samp 1 -dis 100

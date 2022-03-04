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

python code/eig_U.py -mod spin2021 -L 12 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 1.0
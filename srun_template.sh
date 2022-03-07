#!/bin/bash
#$ -l h_rt=336:00:00  # wall-clock time limit (max 336h)
#$ -pe shared 32  # jobs per node
#$ -l h_data=16G  # memory per job
#$ -l h_vmem=512G  # total memory
#$ -cwd  # use the current working directory
#$ -o stdout.$JOB_ID  # file name of standard output
#$ -j y  # combine stdout and stderr
#$ -l highp  # for jobs up to 14 days

source /u/home/b/baandr12/.bash_profile
source /u/home/b/baandr12/.bashrc

python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.1;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.2;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.3;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.4;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.5;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.6;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.7;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.8;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 0.9;
python code/eig_U.py -mod spin2021 -L 8 -bc o -W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta 1;

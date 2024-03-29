#!/bin/bash

SLOGS="/u/home/b/baandr12/PycharmProjects/erg_loc/slogs/ent_delta_flow/"
SJOBS_FILES="/u/home/b/baandr12/PycharmProjects/erg_loc/sjobs/ent_delta_flow/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k6,6n)  # sort by L
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	qsub -l h_rt=336:00:00,h_data=64G,h_vmem=512G,highp -pe shared 16 -cwd -o $SLOGS/"$FILE".out -j y -m bea -M bandrews@physics.ucla.edu "$FILE_PATH"
done
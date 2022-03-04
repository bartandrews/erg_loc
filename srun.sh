#!/bin/bash

FILES="sjobs/ent_delta_flow/*"
for f in $FILES
do
    qsub "$f"
done

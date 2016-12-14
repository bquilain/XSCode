#!/bin/sh +x
for i in {14000..15000}; do
    for j in {0..300}; do
        if [ -f DST${i}_${j}.sh ]
            then
            condor_submit condorDST${i}_${j}.sh
        fi
    done
done

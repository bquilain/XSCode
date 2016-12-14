#!/bin/zsh +x
for i in {1..100};do
#    condor_submit condorPMMC_Sand${i}.sh
    condor_submit condorPMMC${i}.sh
done

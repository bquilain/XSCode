#!/bin/zsh +x
#for i in condorAnaBeam*.sh;do
for i in {1..100};do
    condor_submit condorAnaBeam${i}.sh
done
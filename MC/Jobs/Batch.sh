#!/bin/zsh +x
for i in {0..49};do
    condor_submit condorPMMC${i}_Systematics0_0.sh
done

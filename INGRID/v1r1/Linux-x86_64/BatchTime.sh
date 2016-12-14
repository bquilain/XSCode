#!/bin/zsh
#stop at 15791-39 included
export k=0
#sleep 4h
for i in {13062..13094}; do
#    sleep 1h;
    for j in {0..300}; do
#    export j=0
	if [ -f DST${i}_${j}.sh ]
	    then
	    condor_submit condorDST${i}_${j}.sh
	    k=$(($k+1))
	    echo $k
	fi
	if [ $k -eq 20 ]
	    then
	    sleep 4h;
	    k=$((0))
	fi
    done
done

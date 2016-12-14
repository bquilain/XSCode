#!/bin/zsh +x
#for i in /home/bquilain/CC0pi_XS/MC/JobsSystematics/condorPM*;do
#    condor_submit ${i}
#done

export k=1
printf $k a
export b=0$k
echo $b

for i in {0..0};do
#    echo $i
    for j in {0..0};do
#	echo $j
	for k in {1..500};do
	    condor_submit /home/bquilain/CC0pi_XS/MC/JobsSystematics/condorPMMC${k}_Systematics${i}_${j}.sh
	done

#	for k in {14880..15999};do
#	    echo $k
#	    for l in {0..300};do
#		export l1=0$l
#		k1=$(printf "%08d" $k)
#		l1=$(printf "%04d" $l)
#		echo /export/scraid2/data/bquilain/calib_root_file/ingrid_${k1}_${l1}_Calib00.root
#		if [ -f /export/scraid2/data/bquilain/calib_root_file/ingrid_${k1}_${l1}_Calib00.root ]
#		then
#		echo $k1
#		echo $l1
#		    condor_submit /home/bquilain/CC0pi_XS/MC/JobsSystematics/condorPMData_Run${k}_SubRun${l}_Systematics${i}_${j}.sh
#		fi
#	    done
	done
    done
done

#for i in {17..34};do
#for i in {2..7};do
#    echo $i
#    for j in {0..500};do
#    for j in {0..10};do
#	echo $j
#	condor_submit /home/bquilain/CC0pi_XS/MC/JobsSystematics/condorPostrequisite_Systematics${i}_${j}.sh
#    done
#done

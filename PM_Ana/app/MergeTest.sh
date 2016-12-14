#! /bin/bash
isample=${1:-1}
case $isample in
    1000)
	k=0
	l=1000
	;;
    2000)
	k=1000
        l=2000
        ;;
    3000)
	k=2000
        l=3000
        ;;
    4000)
        k=3000
        l=4000
        ;;
    5000)
	k=4000
        l=5000
        ;;
    6000)
        k=5000
        l=6000
        ;;
    7000)
        k=6000
        l=7000
        ;;
    8000)
        k=8000
        l=9000
        ;;
    9000)
        k=9000
        l=10000
        ;;
    *)
	echo "errors"
	exit -1
	;;
esac

for i in `ls /export/scraid2/data/bquilain/LV/ingrid_000*.root `; do
    j=$((j+1))
    echo ${i}
    echo ${j}
#    test "$j" -eq "$k"

    if [ "$j" -le "$l" ]
	then 
	if [ "$j" -gt "$k" ]
	    then echo "we are the best"
	    MYFILES="$MYFILES $i"
	#echo $MYFILES
	fi
    fi
done
echo $MYFILES
hadd TestTom${isample}.root $MYFILES
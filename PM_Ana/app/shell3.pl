#!/usr/bin/perl

for($i=1; $i<200; $i++){
    $file = "/home/bronner/KyotoTracking/IngRec/app/AddedNoise/Final/nd3_${i}.root";
    if(-f $file){
	$out1="/home/ingrid/data/checknoise/nd3_${i}_tcls.root";
	system "./IngTClster -f $file -o $out1";
	$out2="/home/ingrid/data/checknoise/nd3_${i}_brec.root";
	system "./IngBasicRecon -f $out1 -o $out2";
	$out3="/home/ingrid/data/checknoise/nd3_${i}_trk.root";
	system "./IngReconMuBasic -f $out2 -o $out3";
    }
}

#!/usr/bin/perl
$SUM_PATH="/online/daqdata";
$BIN="/home/daq/INGRID_analysis/app";
$match_par="Comment";
$match_comment="test for MPPC gain";


for($i=1863;$i<1887;$i++){
    $file = sprintf("${SUM_PATH}/ingrid_%08d_0000.daq.mid.summary",$i);
    if(!(-f $file)){
	print("not exist ${file}\n");
	exit;
    }
    open SUM_FILE,"< $file";

    @array = <SUM_FILE>;

    foreach $line(@array){
	chomp($line);
	if($line eq "Comment : test for MPPC gain"){
	    print("$i\n");
	    system "$BIN/gain_problem -r 1754 -s ${i}";
	}
    }

}

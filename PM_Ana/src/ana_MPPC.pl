#!/usr/bin/perl
$PATH="/home/daq/data/root_file";
$BIN="/home/daq/INGRID_analysis/app";

for($i=1888;$i<2000;$i++){
    $file = sprintf("${PATH}/ingrid_%08d_0000.daq.mid.new.root.ana_MPPC",$i);
    if(!(-f $file)){
	print("not exist ${file}\n");
	system "${BIN}/convert_adc_to_pe -r $i -t 1";
	system "echo ${BIN}/convert_adc_to_pe -r $i -t 1";
	system "touch $file";
    }

}

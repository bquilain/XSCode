#!/usr/bin/perl
#for($i = 12399; $i<12414; $i++){

for($i = 12422; $i<12423; $i++){
##for($i = 12339; $i<12340; $i++){

    for($j = 12; $j<90; $j++){
    #for($j = 0; $j<1; $j++){
        #$file = sprintf("/home/daq/data/calib_table/ingrid_%08d_%04d_MPPCcalib.txt",$i, $j);

	$file = sprintf("/home/daq/data/dst/ingrid_%08d_%04d_processed.root",$i, $j);
        if(-f $file){
	    printf("there is $file\n");
	    #system "./check_gain_new -r ${i} -s ${j}";
	    system "./IngTClster -r ${i} -s ${j}";
	    
	    system "./IngBasicRecon -r ${i} -s ${j}";
	    $file = sprintf("/home/daq/data/dst/ingrid_%08d_%04d_tclster.root",$i, $j);
	    system "rm $file";
	    system "./IngReconMuBasic -r ${i} -s ${j}";
	    $file = sprintf("/home/daq/data/dst/ingrid_%08d_%04d_recon2.root",$i, $j);
	    system "rm $file";

	}
    }
}

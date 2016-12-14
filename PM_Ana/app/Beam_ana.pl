#!/usr/bin/perl


for($i = 1800;$i<2000;$i++){
    $file = sprintf("/home/daq/data/beam/ingrid_%08d_beam_shot_unixtime.txt",$i);
    
    if (-f $file){
	system"echo '${file}\n'";
	system"../Linux-x86_64/Midas_to_Root.exe -b ${i}";
    }#if -f
}#for i

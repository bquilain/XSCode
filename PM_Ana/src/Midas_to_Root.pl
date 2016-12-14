#!/usr/bin/perl
$PATH="/online/daqdata";
$BIN="/home/daq/offline/work/INGRID/v1r1/Linux-x86_64";
system "source /home/daq/offline/work/Run_At_Start.bashrc";
system "source /home/daq/offline/work/INGRID/v1r1/cmt/setup.sh";
for($i=1888;$i<2000;$i++){
    $file = sprintf("${PATH}/ingrid_%08d_0000.daq.mid.Midas_to_Root",$i);
    if(!(-f $file)){
	print("not exist ${file}\n");
	system "echo ${BIN}/Midas_to_Root.exe -r $i";
	system "${BIN}/Midas_to_Root.exe -r $i";
       	system "touch $file";
    }
}
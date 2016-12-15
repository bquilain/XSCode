#!/usr/bin/perl


system "cp INGRIDEVENTSUMMARY.so        /cern/root/root.v5.20.00/lib/libINGRIDEVENTSUMMARY.so";
system "cp IngridHitSummary.so        /cern/root/root.v5.20.00/lib/libIngridHitSummary.so";
system "cp IngridSimHitSummary.so     /cern/root/root.v5.20.00/lib/libIngridSimHitSummary.so";
system "cp IngridSimVertexSummary.so     /cern/root/root.v5.20.00/lib/libIngridSimVertexSummary.so";
system "cp IngridSimParticleSummary.so     /cern/root/root.v5.20.00/lib/libIngridSimParticleSummary.so";
system "cp BeamInfoSummary.so     /cern/root/root.v5.20.00/lib/libBeamInfoSummary.so";
system "cp IngridBasicReconSummary.so     /cern/root/root.v5.20.00/lib/libIngridBasicReconSummary.so";
system "cp Ingrid1stReducSummary.so     /cern/root/root.v5.20.00/lib/libIngrid1stReducSummary.so";
system "cp NeutInfoSummary.so     /cern/root/root.v5.20.00/lib/libNeutInfoSummary.so";
system "cp IngridTrackSummary.so     /cern/root/root.v5.20.00/lib/libIngridTrackSummary.so";

if(0){
system "cp INGRIDEVENTSUMMARY.so        /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libINGRIDEVENTSUMMARY.so";
system "cp IngridHitSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngridHitSummary.so";
system "cp IngridSimHitSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngridSimHitSummary.so";
system "cp IngridSimVertexSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngridSimVertexSummary.so";
system "cp IngridSimParticleSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngridSimParticleSummary.so";
system "cp BeamInfoSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libBeamInfoSummary.so";
system "cp IngridBasicReconSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngridBasicReconSummary.so";
system "cp Ingrid1stReducSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngrid1stReducSummary.so";
system "cp NeutInfoSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libNeutInfoSummary.so";
system "cp IngridTrackSummary.so          /home/daq/offline/work/ROOT/v5r20p00n02/Linux-x86_64/lib/root/libIngridTrackSummary.so";
}
if(0){# sysbolic link for DSTMaker
system "ln -s INGRIDEVENTSUMMARY.cc        ~/offline/work/INGRID/v1r1/src/INGRIDEVENTSUMMARY.cc";
system "ln -s IngridHitSummary.cc          ~/offline/work/INGRID/v1r1/src/IngridHitSummary.cc";
system "ln -s IngridSimHitSummary.cc       ~/offline/work/INGRID/v1r1/src/IngridSimHitSummary.cc";
system "ln -s IngridSimVertexSummary.cc    ~/offline/work/INGRID/v1r1/src/IngridSimVertexSummary.cc";
system "ln -s IngridSimParticleSummary.cc  ~/offline/work/INGRID/v1r1/src/IngridSimParticleSummary.cc";
system "ln -s BeamInfoSummary.cc           ~/offline/work/INGRID/v1r1/src/BeamInfoSummary.cc";
system "ln -s IngridBasicReconSummary.cc   ~/offline/work/INGRID/v1r1/src/IngridBasicReconSummary.cc";
system "ln -s Ingrid1stReducSummary.cc     ~/offline/work/INGRID/v1r1/src/Ingrid1stReducSummary.cc";
system "ln -s NeutInfoSummary.cc           ~/offline/work/INGRID/v1r1/src/NeutInfoSummary.cc";
system "ln -s IngridTrackSummary.cc        ~/offline/work/INGRID/v1r1/src/IngridTrackSummary.cc";



system "ln -s INGRIDEVENTSUMMARY.h        ~/offline/work/INGRID/v1r1/src/INGRIDEVENTSUMMARY.h";
system "ln -s IngridHitSummary.h          ~/offline/work/INGRID/v1r1/src/IngridHitSummary.h";
system "ln -s IngridSimHitSummary.h       ~/offline/work/INGRID/v1r1/src/IngridSimHitSummary.h";
system "ln -s IngridSimVertexSummary.h    ~/offline/work/INGRID/v1r1/src/IngridSimVertexSummary.h";
system "ln -s IngridSimParticleSummary.h  ~/offline/work/INGRID/v1r1/src/IngridSimParticleSummary.h";
system "ln -s BeamInfoSummary.h           ~/offline/work/INGRID/v1r1/src/BeamInfoSummary.h";
system "ln -s IngridBasicReconSummary.h   ~/offline/work/INGRID/v1r1/src/IngridBasicReconSummary.h";
system "ln -s Ingrid1stReducSummary.h     ~/offline/work/INGRID/v1r1/src/Ingrid1stReducSummary.h";
system "ln -s NeutInfoSummary.h           ~/offline/work/INGRID/v1r1/src/NeutInfoSummary.h";
system "ln -s IngridTrackSummary.h        ~/offline/work/INGRID/v1r1/src/IngridTrackSummary.h";
}

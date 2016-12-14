#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=150;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command11[300];
  int Run=14510;
  for(int i=12;i<=NFiles;i++){
    sprintf(Name1,"AnaBeam%d.sh",i);
    sprintf(Name2,"condorAnaBeam%d.sh",i);
    //sprintf(Command1,"./PMrecon_HitInfo -f ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise.root -o ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root",i,i);
    //sprintf(Command1,"PM_Ana/app/IngAddNoiseMC -f Ingrid_Process/Ing_%d_MIP5.root -o Ingrid_Process/Ing_MC_RealBeam%d_NewGainNewMIP5_wNoise.root",i,i);
    //sprintf(Command1,"./PMreconRev -r 14510 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_calib.root -o Test.root",i,i);
    
    sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSINGRIDana_woXTalk.root",Run,i,Run,i,Run,i,Run,i);
    sprintf(Command11,"/home/bquilain/Ingrid_Process/ingrid_format/app/IngAddBSD -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSINGRIDana_woXTalk.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSINGRIDanabsd_woXTalk.root -v",Run,i,Run,i,Run,i,Run,i);

    //sprintf(Command1,"./PMrecon_HitInfo -r 14515 -s %d -f /export/scraid2/data/bquilain/ingrid_00014515_%04d_calib.root -o /export/scraid2/data/bquilain/ingrid_00014515_%04d_recon.root",i,i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14515 -s %d -f /export/scraid2/data/bquilain/ingrid_00014515_%04d_recon.root -o /export/scraid2/data/bquilain/ingrid_00014515_%04d_ana.root",i,i,i);
    
    
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Sand%d_NewGainNewMIP9_Form1_wNoise_recon.root -o ../../Ing_MC_Sand%d_NewGainNewMIP9_Form1_wNoise_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root -o ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_ana.root",i,i);

    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f ../../PM_MC_RealBeam%d_NewGainNewMIP6.root -o ../../PM_MC_RealBeam%d_NewGainNewMIP6_recon.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../PM_MC_RealBeam%d_NewGainNewMIP6_recon.root -o ../../PM_MC_RealBeam%d_NewGainNewMIP6_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_recon.root -o ../../Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_ana.root",i,i);
    ofstream Script(Name1);
    if(Script)
      {
	Script<<"#!/bin/bash +x"<<endl
	      <<Command1<<endl
	      <<"cd ../../ingrid_format/app"<<endl
	      <<Command11<<endl;
      }


    sprintf(Command2,"Executable = AnaBeam%d.sh",i);
    sprintf(Command3,"Output = condor_AnaBeamlog%d.txt",i);
    sprintf(Command4,"Error = condor_AnaBeamerr%d.txt",i);
    ofstream Condor(Name2);
    if(Condor){
      Condor<<Command2<<endl
	    <<"Universe        = vanilla"<<endl
	    <<"Rank            = kflops"<<endl
	    <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
	    <<"Getenv          = True"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<Command3<<endl
	    <<Command4<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }

  }

}

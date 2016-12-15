#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){
  char Command0[300];
  int NFiles=150;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command11[300], Command12[300];
  int Run=14514;
  for(int i=0;i<=NFiles;i++){
    sprintf(Name1,"AnaBeam%d.sh",i);
    sprintf(Name2,"condorAnaBeam%d.sh",i);
    //sprintf(Command1,"./PMrecon_HitInfo -f ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise.root -o ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root",i,i);
    //sprintf(Command1,"PM_Ana/app/IngAddNoiseMC -f Ingrid_Process/Ing_%d_MIP5.root -o Ingrid_Process/Ing_MC_RealBeam%d_NewGainNewMIP5_wNoise.root",i,i);
    //sprintf(Command12,"./IngMergePM -f ../../ingrid_00014514_00%d_calib.root -a ../../ingrid_00014514_00%d_pmcalib.root -o ../../ingrid_00014514_00%d_merged.root -r 14514 -s %d",i,i,i,i);
    //sprintf(Command1,"./PMrecon_HitInfo -r 14514 -s %d -f ../../ingrid_00014514_00%d_calib.root -o ../../ingrid_00014514_00%d_recon.root",i,i,i);
   

    //sprintf(Command0,"./IngMergePM -r %d -s %d -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmcalib.root -a /export/scraid2/data/bquilain/ingrid_%08d_%04d_calib.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmerged.root",Run,i,Run,i,Run,i,Run,i);
    //sprintf(Command1,"./PMreconRev_woXTalk -r %d -s %d -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmerged.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedrecon.root",Run,i,Run,i,Run,i);
    //sprintf(Command11,"./PMAna_ForXsections -r %d -s %d -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedrecon.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedana.root",Run,i,Run,i,Run,i); 

    /*    sprintf(Command1,"./PMreconRev -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmerged.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedoldrecon_woTalk.root",Run,i,Run,i,Run,i,Run,i);
    sprintf(Command11,"./PMAna_HitInfo_Xsec -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedoldrecon_woTalk.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedoldana_woXTalk.root",Run,i,Run,i,Run,i,Run,i);
    */
    sprintf(Command1,"/export/scraid2/data/bquilain/Programs/INGRID_analysis/app/PMreconRev_KS -r %d -s %d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmerged.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSrecon.root",Run,i,Run,i,Run,i,Run,i);
   // /export/scraid2/data/bquilain/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14514 -s 0 -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSrecon.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSINGRIDana.root
   ///export/scraid2/data/bquilain/Programs/ingrid_format/app/IngAddBSD -r 14514 -s 0 -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSINGRIDana.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSINGRIDanabsd.root -v
    sprintf(Command11,"/export/scraid2/data/bquilain/Programs/INGRID_analysis/app/PMAnaRev_KS -r %d -s %d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSrecon.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPMana.root",Run,i,Run,i,Run,i,Run,i);
    /*
    sprintf(Command1,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmerged.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root",Run,i,Run,i,Run,i,Run,i);
    sprintf(Command11,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSana_woXTalk.root",Run,i,Run,i,Run,i,Run,i);
    sprintf(Command12,"/home/bquilain/Ingrid_Process/ingrid_format/app/IngAddBSD -r 14000 -s 0 -f /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSana_woXTalk.root -o /export/scraid2/data/bquilain/ingrid_%08d_%04d_pmmergedKSanabsd_woXTalk.root -v -p",Run,i,Run,i,Run,i,Run,i);
    */

    //sprintf(Command3,"./PMreconRev_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise_recon.root",i,i);
    //  sprintf(Command5,"./PMAna_ForXsections -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise_ana.root",i,i);

    //  sprintf(Command1,"./PMPEinPMsci -r 14514 -s %d -f /export/scraid2/data/bquilain/ingrid_00014514_%04d_pmcalib.root -o /export/scraid2/data/bquilain/Datafiles/ingrid_00014514_%04d_pmrecon.root",i,i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14514 -s %d -f ../../ingrid_00014514_00%d_recon.root -o ../../ingrid_00014514_00%d_ana.root",i,i,i);
    //sprintf(Command11,"./PMAna_woIngrid -r 14514 -s %d -f /export/scraid1/data/bquilain/Datafiles/ingrid_00014514_%04d_pmrecon.root -o /export/scraid1/data/bquilain/Datafiles/ingrid_00014514_%04d_pmana.root",i,i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Sand%d_NewGainNewMIP9_Form1_wNoise_recon.root -o ../../Ing_MC_Sand%d_NewGainNewMIP9_Form1_wNoise_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root -o ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_ana.root",i,i);

    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f ../../PM_MC_RealBeam%d_NewGainNewMIP6.root -o ../../PM_MC_RealBeam%d_NewGainNewMIP6_recon.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../PM_MC_RealBeam%d_NewGainNewMIP6_recon.root -o ../../PM_MC_RealBeam%d_NewGainNewMIP6_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_recon.root -o ../../Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_ana.root",i,i);
    ofstream Script(Name1);
    if(Script)
      {
	Script<<"#!/bin/bash +x"<<endl
	  //<<Command12<<endl
	  //     <<Command0<<endl
	      <<Command1<<endl
	      <<Command11<<endl;
	  //    <<"cd ../../ingrid_format/app"<<endl
	  //    <<Command12<<endl;
      }

    sprintf(Command2,"Executable = AnaBeam%d.sh",i);
    sprintf(Command3,"Output = condor_AnaBeamlog%d.txt",i);
    sprintf(Command4,"Error = condor_AnaBeamerr%d.txt",i);
    ofstream Condor(Name2);
    if(Condor){
      Condor<<Command2<<endl
	    <<"Universe        = vanilla"<<endl
	    <<"Rank            = kflops"<<endl
	    <<"Getenv          = True"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\")"<<endl

	    <<Command3<<endl
	    <<Command4<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }

  }

}

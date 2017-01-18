#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=1000;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command5[300], Command11[300], Command55[300];
  char Command10[300], Command20[300], Command30[300], Command40[300], Command50[300], Command60[300];
  for(int i=1;i<=NFiles;i++){
    sprintf(Name1,"AnaBeam%d.sh",i);
    //sprintf(Command1,"./PMrecon_HitInfo -f /export/scraid1/data/bquilain/MCfiles/Ing_MC_Beam%d_NewGainNewMIP9_wNoise.root -o /export/scraid1/data/bquilain/MCfiles/Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root",i,i);
    //sprintf(Command1,"../../../PM_Ana/app/IngAddNoiseMC -r 14000 -s 0 -f ../../PM_MC_Beam%d_NewGainNewMIP11Xsec.root -o ../../PM_MC_Beam%d_NewGainNewMIP11Xsec_wNoise.root",i,i);
    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f /export/scraid1/data/bquilain/MCfiles/Ing_MC_Sand%d_NewGainNewMIP9_wNoise.root -o /export/scraid1/data/bquilain/MCfiles/Ing_MC_Sand%d_NewGainNewMIP9_wNoise_recon.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid1/data/bquilain/MCfiles/Ing_MC_Sand%d_NewGainNewMIP9_wNoise_recon.root -o /export/scraid1/data/bquilain/MCfiles/Ing_MC_Sand%d_NewGainNewMIP9_wNoise_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid1/data/bquilain/MCfiles/Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root -o /export/scraid1/data/bquilain/MCfiles/Ing_MC_Beam%d_NewGainNewMIP9_wNoise_ana.root",i,i);

    //sprintf(Command5,"./PMAna_HitInfo_Xsec -r 14000 -s 0 -f /export/scraid1/data/bquilain/MCfiles/PM_MC_Beam%d_NewGainNewMIP11Xsec_recon.root -o /export/scraid1/data/bquilain/MCfiles/PM_MC_Beam%d_NewGainNewMIP11Xsec_ana.root",i,i);
    //    sprintf(Command1,"./PMreconRev -r 14000 -s 0 -f ../../PM_MC_Beam%d_NewGainNewMIP11Xsec.root -o ../../PM_MC_Beam%d_NewGainNewMIP11Xsec_recon.root",i,i);
    //  sprintf(Command5,"./PMAna_woIngrid -r 14000 -s 0 -f ../../PM_MC_Beam%d_NewGainNewMIP11Xsec_recon.root -o ../../PM_MC_Beam%d_NewGainNewMIP11Xsec_ana.root",i,i);


    //sprintf(Command3,"./PMreconRev_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_wNoise_recon.root",i,i);
    //sprintf(Command5,"./PMAna_ForXsections -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_wNoise_ana.root",i,i);



    /*        sprintf(Command3,"./PMreconRev -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_oldrecon_woXTalk.root",i,i);
     sprintf(Command5,"./PMAna_HitInfo_Xsec -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_oldrecon_woXTalk.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP41_oldana_woXTalk.root",i,i);
    */

    //sprintf(Command1,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scbn07/data2/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk_Inefficiency.root -o /export/scbn07/data2/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk_Reduced_Inefficiency_KSrecon_woXTalk.root",i,i);
    //sprintf(Command5,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scbn07/data2/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk_Reduced_Inefficiency_KSrecon_woXTalk.root -o /export/scbn07/data2/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk_Reduced_Inefficiency_KSana_woXTalk.root",i,i);

    //sprintf(Command1,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_wNoise_KSrecon_woXTalk.root",i,i);
	//sprintf(Command5,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_wNoise_KSrecon_woXTalk.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_wNoise_KSana_woXTalk.root",i,i);
    /*    sprintf(Command1,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root",i,i);
      sprintf(Command5,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);
    ///export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_wNoise_ana.root
    sprintf(Command10,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_recon_woXTalk.root",i,i);
    sprintf(Command20,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_recon_woXTalk.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_PMana_woXTalk.root",i,i);
    sprintf(Command30,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o/export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root",i,i);
    sprintf(Command40,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o/export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);
    sprintf(Command50,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o/export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root",i,i);
    sprintf(Command60,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o/export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);
    */



    //sprintf(Command1,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXtalk.root",i,i);
    sprintf(Command5,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk.root",i,i);
    sprintf(Command10,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXtalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk.root",i,i);
      sprintf(Command50,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);
    sprintf(Command20,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_recon_woXTalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_PMana_woXTalk.root",i,i);
    sprintf(Command40,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o/export/scraid3/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);
    sprintf(Command60,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o/export/scraid3/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);



    //sprintf(Command1,"./PMreconRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_wNoise_KSrecon.root",i,i);
    //sprintf(Command5,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_wNoise_KSrecon.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_wNoise_KSana.root",i,i);

    //sprintf(Command1,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXtalk.root",i,i);
    //sprintf(Command5,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXtalk.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk.root",i,i);

    //sprintf(Command11,"./PMreconRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_wNoise_KSrecon.root",i,i);
    //sprintf(Command55,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_wNoise_KSrecon.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_wNoise_KSana.root",i,i);

    //sprintf(Command11,"./PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root",i,i);
    //sprintf(Command55,"./PMAnaRev_KS -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXTalk.root",i,i);
  
    //sprintf(Command3,"./PMreconRev -r 14000 -s 0 -f ../../PM_MC_Beam%d_MIP9_NuE.root -o ../../PM_MC_Beam%d_MIP9_NuE_recon.root",i,i);
    //sprintf(Command5,"./PMAna_woIngrid -r 14000 -s 0 -f ../../PM_MC_Beam%d_MIP9_NuE_recon.root -o ../../PM_MC_Beam%d_MIP9_NuE_ana.root",i,i);
    // sprintf(Command3,"./PMreconRev -r 14000 -s 0 -f ../../PM_MC_Beam%d_Final3_ShapeIngrid_wNoise2.root -o ../../PM_MC_Beam%d_Final3_ShapeIngrid_wNoise2_recon.root",i,i);
    // sprintf(Command5,"./PMAna_woIngrid -r 14000 -s 0 -f ../../PM_MC_Beam%d_Final3_ShapeIngrid_wNoise2_recon.root -o ../../PM_MC_Beam%d_Final3_ShapeIngrid_wNoise2_ana.root",i,i);
    //sprintf(Command3,"./PMreconRev -r 14000 -s 0 -f /home/kikawa/scbn07/data/ingrid_00015835_%04d_cls.root -o ../../Kikawa/ingrid_00015835_%04d_recon.root",i,i);
    //sprintf(Command5,"./PMAna_woIngrid -r 14000 -s 0 -f ../../Kikawa/ingrid_00015835_%04d_recon.root -o ../../Kikawa/ingrid_00015835_%04d_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid1/data/bquilain/MCfiles/PM_MC_RealBeam%d_NewGainNewMIP6_recon.root -o /export/scraid1/data/bquilain/MCfiles/PM_MC_RealBeam%d_NewGainNewMIP6_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid1/data/bquilain/MCfiles/Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_recon.root -o /export/scraid1/data/bquilain/MCfiles/Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_ana.root",i,i);
    ofstream Script(Name1);
    if(Script)
      {
	Script<<"#!/bin/bash +x"<<endl
	  //	      <<Command1<<endl
	  //<<Command5<<endl
	<<Command10<<endl
	  //<<Command20<<endl
	  //	<<Command30<<endl
	  //<<Command40<<endl 
	<<Command50<<endl
	<<Command60<<endl;
	  // <<Command3<<endl
	  //<<Command11<<endl
	//<<Command55<<endl;
      }


    sprintf(Command2,"Executable = AnaBeam%d.sh",i);
    sprintf(Command3,"Output = condor_AnaBeamlog%d.txt",i);
    sprintf(Command4,"Error = condor_AnaBeamerr%d.txt",i);
    sprintf(Name2,"condorAnaBeam%d.sh",i);
    ofstream Condor(Name2);
    cout<<Name2<<endl;
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

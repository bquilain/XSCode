#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=1000;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command10[300], Command5[300];
  char Command20[300], Command30[300], Command40[300], Command100[300], Command50[300];
  for(int i=1;i<=NFiles;i++){
    //sprintf(Name1,"IngRealBeam%d.sh",i);
    //sprintf(Name2,"condorBeam%d.sh",i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../Ing_MC_Beam%d_NewGainNewMIP9_Form4.root -i /export/scraid1/data/kikawa/neutfile/11bfluka_nd3_numu_fe_%d.nt -m 3 -f 1",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../IngCarbon_MC_Beam%d.root -i /export/scbn27/data1/kikawa/neutfile_ch/11bfluka_nd3_numu_ch_%d.nt -m 3 -f 1",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../Ing_MC_Beam%d_MIP9_NuE_2.root -i /export/scbn07/data2/kikawa/neutfile_nue/11bfluka_nd3_nue_fe_%d.nt -m 3 -f 3",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -i /export/scraid2/data/scbn07_backup2/kikawa/neutfile_pm/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 1",i,i);
   
    sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -i /export/scraid2/data/scbn07_backup2/kikawa/neutfile_run2/11bfluka_nd2_numubar_ch_%d.nt -m 2 -f 2",i,i);
    sprintf(Command5,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -i /export/scraid2/data/scbn07_backup2/kikawa/neutfile_run2/11bfluka_nd4_numu_fe_%d.nt -m 2 -f 1",i,i);
    //sprintf(Command10,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -i /export/scraid2/data/bquilain/neutfile/11bfluka_nd3_numu_fe_%d.nt -m 2 -f 1",i,i);
    sprintf(Command20,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -i /export/scraid2/data/scbn07_backup2/kikawa/neutfile_nue/11bfluka_nd2_nue_ch_%d.nt -m 2 -f 3",i,i);
    sprintf(Command10,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDCHBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -i /export/scbn27/data1/kikawa/neutfile_ch/11bfluka_nd3_numu_ch_%d.nt -m 2 -f 1",i,i);

    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected.root -i /export/scbn07/data2/kikawa/neutfile_nd7/11b_nd7_numu_ch_%d_%d.nt -m 2 -f 1",i,k,l);

    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_MIP9_FromIngridV.root -i /export/scbn07/data2/kikawa/neutfile/11bfluka_nd4_numu_fe_%d.nt -m 2 -f 1",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_MIP9_NuE.root -i /export/scbn07/data2/kika//wa/neutfile_nue -m 3 -f 1",i,i); 
    sprintf(Command3,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root",i,i);
    sprintf(Command30,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root",i,i);
    sprintf(Command40,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root",i,i);
    //sprintf(Command50,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root",i,i);
    sprintf(Command50,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDCHBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDCHBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise.root",i,i);
   
 //sprintf(Command50,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_SciBar188.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_SciBar188_wNoise.root",i,i);

    //sprintf(Command5,"/home/bquilain/Ingrid_Process/INGRID_analysis/app/PMreconRev_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise_recon.root",i,i);
    //sprintf(Command10,"/home/bquilain/Ingrid_Process/INGRID_analysis/app/PMAna_ForXsections -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrected_wNoise_ana.root",i,i);

    sprintf(Name1,"IngNuE%d.sh",i);                                                                                                 
    sprintf(Name2,"condorNuE%d.sh",i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../Ing_MC_Sand%d_NewGainNewMIP8.root -i /export/scraid1/data/kikawa/neutfile_nd7/10c_nd7_numu_ch_%d.nt -m 3 -f 1",i,i);
    ofstream Script(Name1);
    if(Script)
    {
      Script<<"#!/bin/bash +x"<<endl
	//<<Command1<<endl
	//<<Command5<<endl
	<<Command10<<endl
	//<<Command20<<endl
	//<<Command3<<endl
	//<<Command30<<endl
	//<<Command40<<endl
	    <<Command50<<endl;

	//   <<Command5<<endl
	//    <<Command10<<endl;
    }


    sprintf(Command2,"Executable = IngNuE%d.sh",i);
    sprintf(Command3,"Output = condor_NuElog%d.txt",i);
    sprintf(Command4,"Error = condor_NuEerr%d.txt",i);
    ofstream Condor(Name2);
    if(Condor){
      Condor<<Command2<<endl
	    <<"Universe        = vanilla"<<endl
	    <<"Rank            = kflops"<<endl
	    <<"Getenv          = True"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<Command3<<endl
	    <<Command4<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }

  }

}

#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=1000;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command10[300], Command5[300], Command6[300], Command7[300];
  char Command20[300], Command30[300], Command40[300], Command100[300], Command50[300];
  for(int i=1;i<=NFiles;i++){
    if(i<=1000) sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d.root -i /home/cvson/scraid2/neutfile5d3d2/wall/11b_nd7_numu_o_%d_0.nt -m 2 -f 1",i,i);
    else sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d.root -i /home/cvson/scraid2/neutfile5d3d2/wallset2/11b_nd7_numu_o_%d_1.nt -m 2 -f 1",i,i-1000);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d.root -i /export/scraid2/data/bquilain/neutfile_pm/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 2",i,i);
    sprintf(Command2,"/home/bquilain/CC0pi_XS/PM_Ana/app/IngAddNoisePMMC_new -f /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise.root",i,i);
    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_recon.root",i,i);
    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_ana.root",i,i);
    sprintf(Command5,"/home/bquilain/CC0pi_XS/Reconstruction/app/SandEventReductionPM -i /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_ana.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root",i,i);
    //sprintf(Command6,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Sand_Run1_%d_Plan.root -f 1 -n 5 -m",i,i);
    sprintf(Command6,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Sand_Run1_%d_Plan.root -f 1 -n 5 -m",i,i);
    sprintf(Command7,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan_BNJDistributions -i /export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Sand_Run1_%d_Plan_BNJ.root -f 1 -n 5 -m",i,i);
    
    sprintf(Name1,"../Jobs/PMMC_Sand%d.sh",i);                                                                                                 
    sprintf(Name2,"../Jobs/condorPMMC_Sand%d.sh",i);

    ofstream Script(Name1);
    if(Script)
    {
      Script<<"#!/bin/bash +x"<<endl
	//	    <<Command1<<endl
	//	    <<Command2<<endl
	//	    <<Command3<<endl
	//    <<Command4<<endl
	//    <<Command5<<endl;
	    <<Command6<<endl
	    <<Command7<<endl;
    }


    sprintf(Command10,"Executable = PMMC_Sand%d.sh",i);
    sprintf(Command20,"Output = condor_PMMC_Sandlog%d.txt",i);
    sprintf(Command30,"Error = condor_PMMC_Sanderr%d.txt",i);
    ofstream Condor(Name2);
    if(Condor){
      Condor<<Command10<<endl
	    <<"Universe        = vanilla"<<endl
	    <<"Rank            = kflops"<<endl
	    <<"Getenv          = True"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
	    <<Command20<<endl
	    <<Command30<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }

  }

}

#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=100;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command5[300], Command11[300], Command55[300];
  char Command10[300], Command20[300], Command30[300], Command40[300], Command50[300], Command60[300];

  for(int i=1;i<=NFiles;i++){
    sprintf(Name1,"AnaBeam%d.sh",i);

    sprintf(Command5,"./PMAna_woIngrid -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSanawoINGRID_woXtalk.root",i,i);
    sprintf(Command10,"./PMAna_woIngrid -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXtalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSanawoINGRID_woXtalk.root",i,i);
      sprintf(Command50,"./PMAna_woIngrid -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_BeamNuMuBar%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSanawoINGRID_woXTalk.root",i,i);
    sprintf(Command20,"./PMAna_woIngrid -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_recon_woXTalk.root -o /export/scraid3/data/bquilain/MCfiles/PM_MC_BeamINGRIDBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_PManawoINGRID_woXTalk.root",i,i);
    sprintf(Command40,"./PMAna_woIngrid -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o/export/scraid3/data/bquilain/MCfiles/PM_MC_BeamINGRIDVerticalBkg%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSanawoINGRID_woXTalk.root",i,i);
    sprintf(Command60,"./PMAna_woIngrid -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSrecon_woXTalk.root -o/export/scraid3/data/bquilain/MCfiles/PM_MC_BeamNuE%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSanawoINGRID_woXTalk.root",i,i);


    ofstream Script(Name1);
    if(Script)
      {
	Script<<"#!/bin/bash +x"<<endl
	  //	      <<Command1<<endl
	      <<Command5<<endl
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

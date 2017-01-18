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
  for(int k=0;k<=2;k++){
    for(double m=70;m<=100;m+=10){
      for(double n=50;n<=100;n+=10){
	for(int i=10;i<=NFiles;i++){
	  sprintf(Name1,"AnaBeam%d_Criteria%d_%2.0f_%2.0f.sh",i,k,m,n);
	  sprintf(Name2,"condorAnaBeam%d_Criteria%d_%2.0f_%2.0f.sh",i,k,m,n);
	  sprintf(Command1,"./VetoAndFV -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_CriteriaVetow%d_%2.0f_%2.0f_wNoise_recon.root -k %d -m %2.0f -n %2.0f",i,i,i,k,m,n,k,m,n);
	  sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_CriteriaVetow%d_%2.0f_%2.0f_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_CriteriaVetow%d_%2.0f_%2.0f_wNoise_ana.root",i,i,k,m,n,i,k,m,n);


	  //sprintf(Command1,"./VetoAndFV -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_KillNeutrons5MeVNewMIP_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_CriteriaVetow%d_%2.0f_%2.0f_wNoise_recon.root -k %d -m %2.0f -n %2.0f",i,i,i,k,m,n,k,m,n);
	  //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_CriteriaVetow%d_%2.0f_%2.0f_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_CriteriaVetow%d_%2.0f_%2.0f_wNoise_ana.root",i,i,k,m,n,i,k,m,n);

	  //sprintf(Command1,"./VetoAndFV -r 14000 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_calib.root -o /export/scraid2/data/bquilain/ingrid_00014510_%04d_CriteriaVetow%d_%2.0f_%2.0f_recon.root -k %d -m %2.0f -n %2.0f",i,i,i,k,m,n,k,m,n);
	  //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_CriteriaVetow%d_%2.0f_%2.0f_recon.root -o /export/scraid2/data/bquilain/ingrid_00014510_%04d_CriteriaVetow%d_%2.0f_%2.0f_ana.root",i,i,k,m,n,i,k,m,n);
	 
	  //sprintf(Command1,"./TrackMatching3D -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_Birks10_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_Criteria%d_%2.0f_%2.0f_wNoise_recon.root -k %d -m %d -n %2.0f",i,i,i,k,m,n,k,m,n);
	  //sprintf(Command1,"./VetoAndFV -r 14000 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_recon.root -o /export/scraid2/data/bquilain/ingrid_00014510_%04d_Criteria%d_%d_%2.0f_ana.root -k %d -m %2.0f -n %2.0f",i,i,i,k,m,n,k,m,n);
	  ofstream Script(Name1);
	  if(Script)
	    {
	      Script<<"#!/bin/bash +x"<<endl;
	      Script<<Command1<<endl;
	      Script<<Command11<<endl;
	    }
	  
	  
	  sprintf(Command2,"Executable = AnaBeam%d_Criteria%d_%2.0f_%2.0f.sh",i,k,m,n);
	  sprintf(Command3,"Output = condor_AnaBeam%d_Criteria%d_%2.0f_%2.0flog.txt",i,k,m,n);
	  sprintf(Command4,"Error = condor_AnaBeam%d_Criteria%d_%2.0f_%2.0ferr.txt",i,k,m,n);
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
    }
  }
}

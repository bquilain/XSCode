#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=100;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command11[300];
  for(int k=3;k<=5;k++){
    for(int m=1;m<=3;m++){
      // if(k!=4 && m>1) continue;
      for(double n=150;n<=200;n+=25){
     
	for(int i=1;i<=100;i++){
	  sprintf(Name1,"AnaBeam%d_Criteria%d_%d_%2.0f.sh",i,k,m,n);
	  sprintf(Name2,"condorAnaBeam%d_Criteria%d_%d_%2.0f.sh",i,k,m,n);
	 
	  sprintf(Command11,"./TrackMatching3D -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Criteria%d_%d_%2.0f_wNoise_ana.root -k %d -m %d -n %2.0f",i,i,i,k,m,n,k,m,n);
	 
	  //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_calib.root -o /export/scraid2/data/bquilain/ingrid_00014510_%04d_recon.root",i,i,i);

	  //sprintf(Command11,"./TrackMatching3D -r 14000 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_recon.root -o /export/scraid2/data/bquilain/ingrid_00014510_%04d_Criteria%d_%d_%2.0f_ana.root -k %d -m %d -n %2.0f",i,i,i,k,m,n,k,m,n);
	  //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/ingrid_00014510_%04d_Criteria%d_%d_%2.0f_recon.root -o /export/scraid2/data/bquilain/ingrid_00014510_%04d_Criteria%d_%d_%2.0f_ana.root",i+10,i+10,k,m,n,i+10,k,m,n);

	  ofstream Script(Name1);
	  if(Script)
	    {
	      Script<<"#!/bin/bash +x"<<endl;
		//	    Script<<Command1<<endl;
	           Script<<Command11<<endl;
	    }
	  
	  
	  sprintf(Command2,"Executable = AnaBeam%d_Criteria%d_%d_%2.0f.sh",i,k,m,n);
	  sprintf(Command3,"Output = condor_AnaBeam_Criteria%d_%d_%2.0flog%d.txt",i,k,m,n);
	  sprintf(Command4,"Error = condor_AnaBeam_Criteria%d_%d_%2.0ferr%d.txt",i,k,m,n);
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

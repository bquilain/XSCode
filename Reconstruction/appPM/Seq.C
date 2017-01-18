#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=115;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command11[300];
  for(int i=109;i<=NFiles;i++){
    sprintf(Name1,"AnaBeam%d.sh",i);
    sprintf(Name2,"condorAnaBeam%d.sh",i);

    sprintf(Command1,"./PMreconRev -r 14510 -s %d -f /export/scraid1/data/bquilain/ingrid_00014510_0%d_pmcalib.root -o test_%d_recon.root",i,i,i);
    sprintf(Command11,"./PMAna_woIngrid -r 14510 -s %d -f test_%d_recon.root -o test_%d_ana.root",i,i,i);

    //   sprintf(Command1,"./PMreconRev -r 14000 -s 0 -f ../../PM_MC_Beam%d_Final3_ShapeIngrid_wNoiseTrue.root -o testMC%d_recon.root",i,i);
    //sprintf(Command11,"./PMAna_woIngrid -r 14000 -s 0 -f testMC%d_recon.root -o testMC%d_ana.root",i,i);

    ofstream Script(Name1);
    if(Script)
      {
	Script<<"#!/bin/bash +x"<<endl
	  //	      <<Command1<<endl
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
	    <<"Getenv          = True"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<Command3<<endl
	    <<Command4<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }

  }

}

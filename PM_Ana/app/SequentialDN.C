#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=1000;
  int NNoise=10;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command5[300], Command6[300], Command10[300];
  for(int i=101;i<=NFiles;i++){
    sprintf(Command1,"./IngNoiseRate -f /export/scraid2/data/bquilain/Datafiles/ingrid_00014510_%04d_calib.root -o TestNoise_14510_%04d.root",i,i);
    sprintf(Command1,"./IngNoiseRate -f /export/scraid2/data/bquilain/Datafiles/ingrid_00014510_%04d_calib.root -o TestNoise_14510_%04d.root",i,i);
    sprintf(Name1,"DNSyst%d_DN.sh",i);
    sprintf(Name2,"condorDNSyst%d_DN.sh",i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../Ing_MC_Sand%d_NewGainNewMIP8.root -i /export/scraid1/data/kikawa/neutfile_nd7/10c_nd7_numu_ch_%d.nt -m 3 -f 1",i,i);
    ofstream Script(Name1);
    if(Script)
    {
      Script<<"#!/bin/bash +x"<<endl
	    <<Command1<<endl;
    }


    sprintf(Command2,"Executable = DNSyst%d_DN.sh",i);
    sprintf(Command3,"Output = condor_DNSystlog%d_DN.txt",i);
    sprintf(Command4,"Error = condor_DNSysterr%d_DN.txt",i);
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


#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=500;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command10[300];
  for(int i=1;i<=NFiles;i++){
    //sprintf(Name1,"IngRealNuE%d.sh",i);
    //sprintf(Name2,"condorNuE%d.sh",i);
    sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0.root -i /export/scraid3/data/bquilain/neutfile/11bfluka_nd3_numu_fe_%d.nt -m 3 -f 1",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%dNuE_BirksCorrectedMIP46.root -i /export/scbn07/data2/kikawa/neutfile_nue/11bfluka_nd3_nue_fe_%d.nt -m 3 -f 3",i,i);

    ///export/scbn07/data2/kikawa/neutfile_nd7
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../IngCarbon_MC_Beam%d.root -i /export/scbn27/data1/kikawa/neutfile_ch/11bfluka_nd3_numu_ch_%d.nt -m 3 -f 1",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%dNuE_KillNeutrons5MevNewMIP.root -i /export/scbn07/data2/kikawa/neutfile_nue/11bfluka_nd3_nue_fe_%d.nt -m 3 -f 3",i,i);


    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../PM_MC_Beam%d_MIP9_FromIngridV.root -i /export/scbn07/data2/kikawa/neutfile/11bfluka_nd4_numu_fe_%d.nt -m 2 -f 1",i,i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../Ing_MC_Beam%d_MIP9_.root -i /export/scbn07/data2/kikawa/neutfile_nue -m 3 -f 1",i,i); 
    //sprintf(Command10,"/home/bquilain/PM_Ana/app/IngAddNoisePMMC_True -f ../PM_MC_Sand%d_NewGainNewMIP11Xsec.root -o ../PM_MC_Sand%d_NewGainNewMIP11Xsec_wNoiseTrue.root",i,i);

        sprintf(Command10,"/home/bquilain/PM_Ana/app/IngAddNoiseMC -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise.root ",i,i);

    sprintf(Name1,"IngNuE%d.sh",i);                                                                                                                                 
    sprintf(Name2,"condorNuE%d.sh",i);
    //sprintf(Command1,"bin/Linux-g++/IngMC -o ../Ing_MC_Sand%d_NewGainNewMIP8.root -i /export/scraid1/data/kikawa/neutfile_nd7/10c_nd7_numu_ch_%d.nt -m 3 -f 1",i,i);
    ofstream Script(Name1);
    if(Script)
    {
      Script<<"#!/bin/bash +x"<<endl
	    <<Command1<<endl;
      	    Script<<Command10<<endl;
    }


    sprintf(Command2,"Executable = IngNuE%d.sh",i);
    sprintf(Command3,"Output = condor_NuElog%d.txt",i);
    sprintf(Command4,"Error = condor_NuEerr%d.txt",i);
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

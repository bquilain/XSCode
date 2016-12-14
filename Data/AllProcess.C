
#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  char Name1[100], Name2[100];
  char Command00[300], Command01[300], Command02[300], Command03[300], Command04[300], Command05[300], Command06[300];
  char Command0[300],Command1[300], Command2[300], Command3[300], Command4[300], Command5[300], Command6[300], Command7[300], Command8[300], Command9[300], Command81[300], Command91[300], Command101[300], Command700[300], Command50[300]; 
  char Command10[300], Command11[300], Command12[300];
  char Command70[300];
  char RunNumber[16];
  char SubRunNumber[16];
  int RunNum;
  int SubRunNum;
  char Name[16];

  ifstream List;
  for(int l=14;l<=14;l++){
    sprintf(Name,"%d000.txt",l);
  List.open(Name);
  char MIDAS[150];
  if(List.is_open()){
    while(!List.eof()){
      List>>MIDAS;
      //cout<<MIDAS<<", 20e="<<MIDAS[50]<<endl;
      sprintf(RunNumber,"%c%c%c%c%c",MIDAS[50],MIDAS[51],MIDAS[52],MIDAS[53],MIDAS[54]);
      sprintf(SubRunNumber,"%c%c%c%c",MIDAS[56],MIDAS[57],MIDAS[58],MIDAS[59]);
      cout<<RunNumber<<", SubRun="<<SubRunNumber<<endl;
      RunNum=atoi(RunNumber);
      SubRunNum=atoi(SubRunNumber);
   
	  sprintf(Name1,"Jobs/DST%d_%d.sh",RunNum,SubRunNum);
	  sprintf(Name2,"Jobs/condorDST%d_%d.sh",RunNum,SubRunNum);

	  sprintf(Command0,"/home/bquilain/T2KSoftware2/INGRID/v1r1/Linux-x86_64/Calc_MPPC_new_ci_2.exe -r %d -s %d -t 1",RunNum,SubRunNum);
	  sprintf(Command1,"if [ -f /export/scraid2/data/bquilain/calib_root_file/ingrid_%08d_%04d_Calib00.root ]",RunNum,SubRunNum);
	  sprintf(Command2,"/home/bquilain/T2KSoftware2/INGRID/v1r1/Linux-x86_64/DSTMaker.exe -r %d -s %d -t 1",RunNum,SubRunNum);
	  sprintf(Command3,"/home/bquilain/T2KSoftware2/INGRID/v1r1/Linux-x86_64/IngCalib_ADCNLCorrected.exe -r %d -s %d -f /export/scraid2/data/bquilain/dst/ingrid_%08d_%04d.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_calib.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
	  sprintf(Command4,"/home/bquilain/T2KSoftware2/INGRID/v1r1/Linux-x86_64/DSTMaker.exe -r %d -s %d -t 1 -p",RunNum,SubRunNum);
	  sprintf(Command50,"/home/bquilain/T2KSoftware2/INGRID/v1r1/Linux-x86_64/IngCalib_ADCNLCorrected.exe -r %d -s %d -f /export/scraid2/data/bquilain/PM/ingrid_%08d_%04d.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmcalib.root -p",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
	  sprintf(Command6,"/home/bquilain/CC0pi_XS/Reconstruction/app/IngMergePM -r %d -s %d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmcalib.root -a /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_calib.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmerged.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
	  sprintf(Command7,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRev_KS_woXTalk -r %d -s %d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmerged.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);

	  sprintf(Command8,"/home/bquilain/CC0pi_XS/Reconstruction/app/IngAnaVertex_HitInfo -r %d -s %d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSINGRIDana_woXTalk.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
	  sprintf(Command9,"~/T2KSoftware2/ingrid_format/app/IngAddBSD -r %d -s %d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSINGRIDana_woXTalk.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSINGRIDanabsd_woXTalk.root -v",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);


	  sprintf(Command81,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRev_KS_woXTalk -r 14000 -s 0 -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
	  sprintf(Command91,"~/T2KSoftware2/ingrid_format/app/IngAddBSD -r %08d -s %04d -f /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk.root -o /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root -v -p",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
	  sprintf(Command101,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi -i /export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_%08d_%04d.root -f 1 -n 5",RunNum,SubRunNum,RunNum,SubRunNum);
	  
	  sprintf(Command700,"QUEUE 1",SubRunNum);
	  ofstream Script(Name1);
	  if(Script)
	    {
	      Script<<"#!/bin/bash +x"<<endl
		//<<Command0<<endl
		//<<Command1<<endl
		//<<"then"<<endl
		//<<Command2<<endl
		//    <<Command3<<endl
		//<<Command4<<endl
		//    <<Command50<<endl
		//<<Command6<<endl
		
		//<<Command7<<endl
		//    <<Command8<<endl
		//    <<Command9<<endl
		//    <<Command81<<endl
		//    <<Command91<<endl
		    <<Command101<<endl;
		//<<"fi"<<endl;
		//		    <<Command01<<endl;
	    }


	  // sprintf(Command4,"Executable = DST%d_%d.sh",RunNum,SubRunNum);

	  //  sprintf(Command5,"Output = condor_DST%d_%dlog.txt",RunNum,SubRunNum);
	  // sprintf(Command6,"Error = condor_DST%d_%derr.txt",RunNum,SubRunNum);
	  sprintf(Command10,"Executable = DST%d_%d.sh",RunNum,SubRunNum);

	  sprintf(Command11,"Output = condor_DST%d_%dlog.txt",RunNum,SubRunNum);
	  sprintf(Command12,"Error = condor_DST%d_%derr.txt",RunNum,SubRunNum);
	  ofstream Condor(Name2);
	  if(Condor){
	    Condor<<Command10<<endl
            <<"Universe        = vanilla"<<endl
            <<"Rank            = kflops"<<endl
            <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
            <<"Universe        = vanilla"<<endl
            <<"Rank            = kflops"<<endl
            <<"Getenv          = True"<<endl
            <<"Arguments      =  $(Process)"<<endl
            <<Command11<<endl
            <<Command12<<endl
            <<"Notification    = Never"<<endl
	      //<<"QUEUE 1"<<endl;
		  <<Command700<<endl;
	      }

    }
  }
  List.close();
  }
}

#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <TMath.h>
using namespace std;
    //0. No Error, nominal case
    //1. TO DO 
    //2. Dark noise, variations 
    //3. Hit efficiency, variations within the difference data and MC
    //4. Light yield, variation of PE with angle btw data and MC.
    //5. Light yield, Birks quenching effect.
    //6: Beam related background (in fact, this mainly evaluate sand muons)
    //7: 2D reconstruction error
    //8: VetoUpstreamCriteria: nominal=0 planes, vary from 0->2 per 1 plane step
    //9: VetoEdgeCriteria: nominal=80cm, vary 70->100 per 10cm steps
    //10: FVCriteria: nominal=100cm, vary 50->100 per 10 cms steps
    //11: Vertexing, plane tolerance: nominal=2, vary 2->4 per 1 plane steps
    //12: Vertexing, transverse tolerance: nominal=15cm, vary 15cm->20cm per 2.5cm steps
    //13: Track matching, plane tolerance: nominal=4, vary 3->5 per 1 plane steps
    //14: INGRID/PM tracks angle matching: nominal=35°, vary 30°->40° per 5° steps
    //15: INGRID/PM tracks transverse position matching: nominal=8.5cm, vary 7.5cm->9.5 per 1cm steps

    //16: Xsec
    //17:Flux
    
 
int main(int argc, char **argv){

  int NFiles=1000;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command5[300], Command6[300];
  char  Command10[300], Command20[300], Command30[300];
  sprintf(Command1,"");
  sprintf(Command2,"");
  sprintf(Command3,"");
  sprintf(Command4,"");
  sprintf(Command5,"");
  sprintf(Command6,"");
  sprintf(Command7,"");
  
  sprintf(Command10,"");
  sprintf(Command20,"");
  sprintf(Command30,"");
  sprintf(Command40,"");

  double ErrorValue;
  int NErrors=16;
  int NE[NErrors];
  for(int n=0;n<NErrors;n++){
    NE[n]=1;
    if(n==2) NE[n]=10;//Number of DN value tested
    else if(n==3) NE[n]=100;//Number of Efficiency number tested
    else if(n==5) NE[n]=2;//Birks constant variations
  }
    
  for(int ErrorType=0;ErrorType<16;ErrorType++){
    for(int n=0;n<NE[0];n++){
	for(int ipln=0;ipln<=2;ipln++){
	  ErrorValue=ipln;
	  if(ErrorType!=8) continue;
	  for(double veto=70;veto<=100;veto+=10){
	    ErrorValue=veto;
	    if(ErrorType!=9) continue;
	    for(double fv=50;fv<=100;fv+=10){
	      ErrorValue=fv;
	      if(ErrorValue!=10) continue;
	      for(int vplan=2;vplan<=4;vplan++){
		ErrorValue=vplan;
		if(ErrorValue!=11) continue;
		for(double vtrans=15.;vtrans<=20.;vtrans+=2.5){
		  ErrorValue=vtrans;
		  if(ErrorValue!=12) continue;
		  for(int trkmplan=3;trkmplan<=5;trkmplan++){
		    ErrorValue=trkmplan;
		    if(ErrorValue!=13) continue;
		    for(double angle=30;angle<=40;angle+=5){
		      ErrorValue=angle;
		      if(ErrorValue!=14) continue;
		      for(double trans=7.5;trans<=8.5;trans+=1.){
			ErrorValue=trans;
			if(ErrorValue!=15) continue;
			for(int i=1;i<=NFiles;i++){
			  if(ErrorType==0){//Case w/o error
			    sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -i /home/cvson/scraid2/neutfile5d3d2/run1/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 1",i,i);
			    //sprintf(Command1,"bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -i /export/scraid2/data/bquilain/neutfile_pm/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 2",i,i);
			    sprintf(Command2,"/home/bquilain/CC0pi_XS/PM_Ana/app/IngAddNoisePMMC_new -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root",i,i);
			    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root",i,i);
			    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root",i,i);
			    sprintf(Command5,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root -f 1 -n 5 -m",i,i);
			    
			    //sprintf(Command6,"source $HOME/T2KSoftware2/T2KReWeight/Run_at_Start.sh");
			    //sprintf(Command7,"source $HOME/T2KSoftware2/T2KReWeight/v1r23/app/genWeightsFromINGRID.exe -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_ReWeight.root",i,i);
			    
			    //sprintf(Command6,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root",i,i);
			    //sprintf(Command5,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d.root -f 1 -n 5 -m",i,i);
			    //sprintf(Command6,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan_BNJDistributions -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_BNJ.root -f 1 -n 5 -m",i,i);
			  }
			  else if(ErrorType==1) continue;
			  else if(ErrorType==2){//Dark noise values.
			    ErrorValue=n;
			    sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%2.1f.root -n %2.0f",i,i,ErrorValue,ErrorValue);	
			    sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%2.1f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%2.1f_recon.root",i,ErrorValue,i,ErrorValue);
			    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%2.1f_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%2.1f_ana.root",i,ErrorValue,i,ErrorValue);
			    sprintf(Command4,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%2.1f_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_DN%2.1f_Plan.root -f 1 -n 5 -m",i,ErrorValue,i,ErrorValue);    
			  }
			  else if(ErrorType==3){//One should give in input the file name of the data/MC difference
			    //The best thing would be just to generate everything w/ an option
			    if(n==0){
			      sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/HitEfficiency -m -i 1 -f 100"); //Generate a file containing MC hit efficiency (XS/files_MCDataComparison/MC_CalibrationPM.root )
			      sprintf(Command2,"/home/bquilain/CC0pi_XS/XS/HitEfficiency -r 14510 -t 14570 -i 0 -f 300"); //Generate a file containing Data hit efficiency (/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/MC_CalibrationPM.root )
			    }
			    ErrorValue=n;
			    sprintf(Command3,"/home/bquilain/CC0pi_XS/XS/GenerateHitEfficiencyMask -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_recon.root",i,i,ErrorValue);
			    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_recon.root",i,ErrorValue,i,ErrorValue);
			    sprintf(Command5,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_ana.root",i,ErrorValue,i,ErrorValue);
			    sprintf(Command6,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_HitEfficiencyRandom%d_Plan.root -f 1 -n 5 -m",i,ErrorValue,i,ErrorValue);    	
			  }
	else if(ErrorType==4){//One should give in input the file name of the data/MC difference
	  //The best thing would be just to generate everything w/ an option
	  if(n==0){
	  sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/CompareCalibrationsPM -m -i 1 -f 100"); //-> Generate a file containing each hit info for MC (XS/files_MCDataComparison/MC_CalibrationPM.root )
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/XS/CompareCalibrationsPM -r 14510 -t 14570 -i 0 -f 300");//-> Generate a file containing each hit info for Data (XS/files_MCDataComparison/Data_CalibrationPM.root )
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/XS/GeneratePEAngleDistributions -o /home/bquilain/CC0pi_XS/XS/files/PEXAngle.root");//Read the data and MC files above and create the dependency of PE with angle.
	  }
	  sprintf(Command4,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root -f 1 -n 5 -m -e %d -v /home/bquilain/CC0pi_XS/XS/files/PEXAngle.root",i,i,ErrorType);    
	}
	else if(ErrorType==5){
	  //pre-requisite: generate MC w/ Birks constant variation from 0.0208+-0.0023 cm/MeV
	  //1. Change birks constant to minus value in src/IngridResponse.cc: 0.0208->0.0185
	  //2. Run this MC
	  ErrorValue=185;
	  if(n==0){
	    sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%3.0f.root -i /home/cvson/scraid2/neutfile5d3d2/run1/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 1",i,i);
	    sprintf(Command2,"/home/bquilain/CC0pi_XS/PM_Ana/app/IngAddNoisePMMC_new -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	    sprintf(Command5,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Systematics%d_%3.0f_Plan.root -f 1 -n 5 -m",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	  }
	  if(n==1){
	    //1. Change birks constant to minus value in src/IngridResponse.cc: 0.0208->0.0231
	    //2. Run this MC
	    ErrorValue=231;
	    sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%3.0f.root -i /home/cvson/scraid2/neutfile5d3d2/run1/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 1",i,i);
	    sprintf(Command2,"/home/bquilain/CC0pi_XS/PM_Ana/app/IngAddNoisePMMC_new -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	    sprintf(Command5,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Systematics%d_%3.0f_Plan.root -f 1 -n 5 -m",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	  }
	}
      else if(ErrorType==6){
	double NominalContamination=0.6;
	double Err=TMath::Sqrt(NominalContamination*NominalContamination+0.2*0.2+0.2*0.2);
	NominalContamination+=1;
	for(double conta=NominalContamination-Err;conta<=NominalContamination+Err;conta+=Err){
	  sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/CC0piSelection -s 1 -w %2.2f",conta);
	}
      }
      else if(ErrorType==7){//Not do anything, since all information is already stored in the trees
	//sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root",i,i);
      }
      else if(ErrorType==8){
	for(int ipln=0;ipln<=2;ipln++){
	  ErrorValue=ipln;
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	  sprintf(Command2,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	}
      }
      else if(ErrorType==9){
	for(double veto=70;veto<=100;veto+=10){
	  ErrorValue=veto;
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	  sprintf(Command2,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	}
      }
      else if(ErrorType==10){
	for(double fv=50;fv<=100;fv+=10){
	  ErrorValue=fv;
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	  sprintf(Command2,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
	}
      }
      else if(ErrorType==11){
	for(int vplan=2;vplan<=4;vplan++){
	ErrorValue=vplan;
	sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
      }
    }
    else if(ErrorType==12){
      for(double vtrans=15.;vtrans<=20.;vtrans+=2.5){
	ErrorValue=vtrans;
	sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
      }
    }
    else if(ErrorType==13){
      for(int trkmplan=3;trkmplan<=5;trkmplan++){
	ErrorValue=trkmplan;
	sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
      }
    }
    else if(ErrorType==14){
      for(double angle=30;angle<=40;angle+=5){
	ErrorValue=angle;
	sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
      }
    }
    else if(ErrorType==15){
      for(double trans=7.5;trans<=8.5;trans+=1.){
	ErrorValue=trans;
	sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root -e %d -v %3.0f",i,i,ErrorType,ErrorValue,ErrorType,ErrorValue);
	sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%3.0f.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%3.0f.root",i,ErrorType,ErrorValue,i,ErrorType,ErrorValue);
      }
    }
    	      

    ofstream Script(Name1);
    if(Script)
    {
      Script<<"#!/bin/bash +x"<<endl
	<<Command1<<endl
	<<Command2<<endl
	<<Command3<<endl
	<<Command4<<endl
        <<Command5<<endl
        <<Command6<<endl;
    }


    sprintf(Command10,"Executable = PMMC%d_Systematics%d_%3.0f.sh",i,ErrorType,ErrorValue);
    sprintf(Command20,"Output = condor_PMMClog%d_Systematics%d_%3.0f.txt",i,ErrorType,ErrorValue);
    sprintf(Command30,"Error = condor_PMMCerr%d_Systematics%d_%3.0f.txt",i,ErrorType,ErrorValue);
    ofstream Condor(Name2);
    if(Condor){
      Condor<<Command10<<endl
	    <<"Universe        = vanilla"<<endl
	    <<"Rank            = kflops"<<endl
	    <<"Getenv          = True"<<endl
	    <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<Command20<<endl
	    <<Command30<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }

  }

}

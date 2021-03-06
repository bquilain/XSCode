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

  int StartError=16;
  int EndError=16;
  int NFiles=1000;
  char Name1[100], Name2[100];
  char Command1[300], Command2[300], Command3[300], Command4[300], Command5[300], Command6[300];
  char  Command10[300], Command20[300], Command30[300];

  float ErrorValue;
  int NE[EndError+1];
  double Step[EndError+1];
  double Start[EndError+1];
  double Nominal,Err;
  //double End[EndError];
  
  for(int n=StartError;n<=EndError;n++){
    NE[n]=1;
    Start[n]=1;
    Step[n]=1;
    
    if(n==2) NE[n]=10;//Number of DN value tested
    else if(n==3) NE[n]=100;//Number of Efficiency number tested
    else if(n==5) NE[n]=2;//Birks constant variations
    else if(n==6){//6: Beam related background (in fact, this mainly evaluate sand muons)
      Nominal=0.6;
      Err=TMath::Sqrt(Nominal*Nominal+0.2*0.2+0.2*0.2);
      Start[n]=Nominal-Err;
      Step[n]=Err;
      NE[n]=3;
    }
    else if(n==7) continue;    //7: 2D reconstruction error

    else if(n==8){    //8: VetoUpstreamCriteria: nominal=0 planes, vary from 0->2 per 1 plane step
      Nominal=0;
      Err=1;
      Start[n]=Nominal;
      Step[n]=Err;
      NE[n]=3;
    }
    else if(n==9){//9: VetoEdgeCriteria: nominal=80cm, vary 70->100 per 10cm steps
      Nominal=80;
      Err=10;
      Start[n]=Nominal-Err;
      Step[n]=Err;
      NE[n]=4;
    }
    else if(n==11){//11: Vertexing, plane tolerance: nominal=2, vary 2->4 per 1 plane steps
      Nominal=2;
      Err=1;
      Start[n]=Nominal;
      Step[n]=Err;
      NE[n]=3;
    }
    else if(n==12){//12: Vertexing, transverse tolerance: nominal=15cm, vary 15cm->20cm per 2.5cm steps
      Nominal=15;
      Err=2.5;
      Start[n]=Nominal;
      Step[n]=Err;
      NE[n]=3;
    }
    else if(n==13){//13: Track matching, plane tolerance: nominal=4, vary 3->5 per 1 plane steps
      Nominal=4;
      Err=1;
      Start[n]=Nominal-Err;
      Step[n]=Err;
      NE[n]=3;
    }
    else if(n==14){//14: INGRID/PM tracks angle matching: nominal=35°, vary 30°->40° per 5° steps
      Nominal=35;
      Err=5;
      Start[n]=Nominal-Err;
      Step[n]=Err;
      NE[n]=3;
    }
    else if(n==15){//15: INGRID/PM tracks transverse position matching: nominal=8.5cm, vary 7.5cm->9.5 per 1cm steps
      Nominal=8.5;
      Err=1.;
      Start[n]=Nominal-Err;
      Step[n]=Err;
      NE[n]=3;
    }
  }
  
  for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      sprintf(Command1,"");
  sprintf(Command2,"");
  sprintf(Command3,"");
  sprintf(Command4,"");
  sprintf(Command5,"");
  sprintf(Command6,"");
  //sprintf(Command7,"");
  
  sprintf(Command10,"");
  sprintf(Command20,"");
  sprintf(Command30,"");
  //sprintf(Command40,"");

    cout<<endl<<"Error Type="<<ErrorType<<endl;
    for(int n=0;n<NE[ErrorType];n++){
      ErrorValue=Start[ErrorType]+n*Step[ErrorType];
      cout<<"Error Value="<<ErrorValue<<", files generated:";
      for(int i=1;i<=NFiles;i++){
	cout<<i<<", ";
	
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
	  sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/IngAddNoisePMMC -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%d.root -n %2.1f",i,i,n,ErrorValue);  
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%d_recon.root",i,n,i,n);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%d_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%d_ana.root",i,n,i,n);
	  sprintf(Command4,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_DN%d_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_DN%d_Plan.root -f 1 -n 5 -m",i,n,i,n);    
	}
	else if(ErrorType==3){//One should give in input the file name of the data/MC difference
	  //The best thing would be just to generate everything w/ an option
	  if(n==0){
	    sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/HitEfficiency -m -i 1 -f 100"); //Generate a file containing MC hit efficiency (XS/files_MCDataComparison/MC_CalibrationPM.root )
	    sprintf(Command2,"/home/bquilain/CC0pi_XS/XS/HitEfficiency -r 14510 -t 14570 -i 0 -f 300"); //Generate a file containing Data hit efficiency (/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/MC_CalibrationPM.root )
	  }
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/XS/GenerateHitEfficiencyMask -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_recon.root",i,i,n);
	  sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_recon.root",i,n,i,n);
	  sprintf(Command5,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_ana.root",i,n,i,n);
	  sprintf(Command6,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_HitEfficiencyRandom%d_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_HitEfficiencyRandom%d_Plan.root -f 1 -n 5 -m",i,n,i,n);    	
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
	  //pre-requisite if you change the MC: do 2 copies and change the Birks constant (cf memo): generate MC w/ Birks constant variation from 0.0208+-0.0023 cm/MeV
	  //ErrorValue=185;
	  	  //1. Change birks constant to minus value in src/IngridResponse.cc: 0.0208->0.0185
	  //2. Run this MC
	    //1. Change birks constant to minus value in src/IngridResponse.cc: 0.0208->0.0231
	    //2. Run this MC

	  if(n==0){
	    sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/Birks_Minus/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%d.root -i /home/cvson/scraid2/neutfile5d3d2/run1/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 1",i,i);
	    sprintf(Command2,"/home/bquilain/CC0pi_XS/PM_Ana/app/IngAddNoisePMMC_new -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	    sprintf(Command5,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Systematics%d_%d_Plan.root -f 1 -n 5 -m",i,ErrorType,n,i,ErrorType,n);
	  }
	  if(n==1){
	    //ErrorValue=231;
	    sprintf(Command1,"/home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/Birks_Plus/bin/Linux-g++/IngMC -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%d.root -i /home/cvson/scraid2/neutfile5d3d2/run1/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 1",i,i);
	    sprintf(Command2,"/home/bquilain/CC0pi_XS/PM_Ana/app/IngAddNoisePMMC_new -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	    sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	    sprintf(Command4,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	    sprintf(Command5,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Systematics%d_%d_Plan.root -f 1 -n 5 -m",i,ErrorType,n,i,ErrorType,n);
	  }
	}
	else if(ErrorType==6){
	  sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/CC0piSelection -s 1 -w %2.2f",ErrorValue);
	}
	else if(ErrorType==7){//Not do anything, since all information is already stored in the trees
	  //sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root",i,i);
	}
	else if(ErrorType==8){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==9){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==10){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMreconRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==11){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==12){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==13){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==14){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==15){
	  sprintf(Command1,"home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v %3.0f",i,i,ErrorType,n,ErrorType,ErrorValue);
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/Reconstruction/app/PMAnaRevOfficial -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",i,ErrorType,n,i,ErrorType,n);
	}
	else if(ErrorType==16){	    
	  sprintf(Command1,"source $HOME/T2KSoftware2/T2KReWeight/Run_at_Start.sh");
	  //sprintf(Command2,"$HOME/T2KSoftware2/T2KReWeight/v1r23/app/genWeightsFromINGRID.exe -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d.root -o /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_ReWeight.root",i,i);
	  sprintf(Command3,"/home/bquilain/CC0pi_XS/XS/XS_CC0pi_Plan -i /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root -o /home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_ReWeight.root -f 1 -n 5 -m -x /export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_ReWeight.root",i,i,i);
	}

	/*
	  sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/CC0piSelection -s 1 -x);

	  sprintf(Command1,"/home/bquilain/CC0pi_XS/XS/CC0piSelection -s 1 -f);
	  for(int iflux=0;iflux<=10;iflux++){
	  sprintf(Command2,"/home/bquilain/CC0pi_XS/XS/UnfoldingOptimisation -d /home/bquilain/CC0pi_XS/XS/files/MCSelected.txt -m /home/bquilain/CC0pi_XS/XS/files/MCSelected_Systematics_Flux%d.txt -n 1",iflux);
	  }
	*/	
	sprintf(Name1,"../JobsSystematics/PMMC%d_Systematics%d_%d.sh",i,ErrorType,n);
	sprintf(Name2,"../JobsSystematics/condorPMMC%d_Systematics%d_%d.sh",i,ErrorType,n);
	cout<<"Name="<<Name1<<endl;
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
	
	
	sprintf(Command10,"Executable = PMMC%d_Systematics%d_%d.sh",i,ErrorType,n);
	sprintf(Command20,"Output = condor_PMMClog%d_Systematics%d_%d.txt",i,ErrorType,n);
	sprintf(Command30,"Error = condor_PMMCerr%d_Systematics%d_%d.txt",i,ErrorType,n);
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
      cout<<endl;
    }
  }
}

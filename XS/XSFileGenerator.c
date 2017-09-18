#include<iomanip>  
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include "setup.h"
//#include "Xsec.cc"
//Xsec * xs = new Xsec();
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
//9: VetoEdgeCriteria: nominal=80mm, vary 60->100 per 20mm steps ; ML 2017/03/14
//10: FVCriteria: nominal=100cm, vary 50->100 per 10 cms steps
//11: Vertexing, plane tolerance: nominal=2, vary 2->4 per 1 plane steps
//12: Vertexing, transverse tolerance: nominal=15cm, vary 15cm->20cm per 2.5cm steps
//13: Track matching, plane tolerance: nominal=4, vary 3->5 per 1 plane steps
//14: INGRID/PM tracks angle matching: nominal=35°, vary 30°->40° per 5° steps
//15: INGRID/PM tracks transverse position matching: nominal=8.5cm, vary 7.5cm->9.5 per 1cm steps

//16: Flux
//17: Xsec
    
 
int main(int argc, char **argv){

  /**************************************Get back environment variables*****************************************************/
  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");
    
  int c=-1;
  bool MC=false;bool Data=false;
  bool SelectionOnly=false;// if true, only the latest part of the analysis (selection & unfolding) are applied
  bool PM=true;
  bool ParticleGun=false;

  bool sandOnly=false;
  int channel=2; // channel for analysis; 0:sand 1:cc0pi 2:cc1pi
  bool retuned=false; // retune the nominal NEUT
   
  while ((c = getopt(argc, argv, "mdspwgSc:t")) != -1) {
    switch(c){
    case 'm': 
      MC=true;
      break;
    case 'd':
      Data=true;
      break;
    case 'w':
      PM=false;//watermodule
      break;
    case 'p':
      PM=true;
      break;
    case 's':
      SelectionOnly=true;
      break;
    case 'S':
      sandOnly=true;
      break;
    case 't':
      retuned=true;
      break;
    case 'g':
      ParticleGun=true;
      break;
    case 'c':
      channel=atoi(optarg);
      break;

    }      
  } 

  char Name1[100], Name2[100];
  char Command001[300], Command002[300], Command003[300], Command004[300], Command005[300], Command006[300];
  char Command01[300], Command02[300], Command03[300], Command04[300], Command05[300], Command06[300];
  char Command1_1[300], Command1_2[300], Command1_3[300], Command1_4[300], Command1_5[300], Command1_6[600];
  char Command1[1000], Command2[300], Command3[300], Command4[300], Command5[300], Command6[300], Command7[300], Command8[300], Command9[300], Command11[300], Command12[300], Command13[300],Command14[300], Command15[300], Command16[300], Command17[300], Command18[300];
  char  Command10[300], Command20[300], Command30[300];
  char Name[16];

  //  xs->Xsec::Initialize();
  InitializeGlobal(PM);
  /////////// VARIABLES TO SWITCH FROM PM TO WM ////////////////////
  char execMC[7];sprintf(execMC,(PM?"IngMC":"IngWMMC"));
  char execRecon[25];sprintf(execRecon,(PM?"PMreconRevOfficial":"Lolirecon"));
  char execAna[25];sprintf(execAna,(PM?"PMAnaRevOfficial":"LoliAna"));
  int MCDetID=(PM?2:5);
  char suffix[3];sprintf(suffix,(PM?"":"_WM"));
  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  char ParticleGenerator[5];sprintf(ParticleGenerator,(ParticleGun?"_PG":""));
  cout<<"Selected detector is "<<DetName<<endl;
  char optionTune[100],fileTune[100]; sprintf(fileTune,"-t ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%s_ReWeight2015.root",suffix,DetName,"%d");
  char Sand[15];sprintf(Sand,(sandOnly?"_Wall":""));
  char InputFile[100];
  cout<<"Selected particle generator is"<<(ParticleGun?" a particle gun":" neutrino event")<<endl;

  
  ///////////FOR THE PREREQUISITE, SO DO NOT CLEAN FOR EACH ERROR/////////////////////
  sprintf(Command01,"");
  sprintf(Command02,"");
  sprintf(Command03,"");
  sprintf(Command04,"");
  sprintf(Command05,"");
  sprintf(Command06,"");

  if(!SelectionOnly){
    for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      //for(int ErrorType=2;ErrorType<=6;ErrorType++){//TEMP

    
      cout<<endl<<"Error Type="<<ErrorType<<endl;

      
      if(ErrorType==3){
	if(PM){
	  sprintf(Command01,"${INSTALLREPOSITORY}/XS/HitEfficiency -m -i 1 -f 1000"); //Generate a file containing MC hit efficiency (XS/files_MCDataComparison/MC_CalibrationPM.root )
	  sprintf(Command02,"${INSTALLREPOSITORY}/XS/HitEfficiency -r 14510 -t 14510 -i 0 -f 300"); //Generate a file containing Data hit efficiency (${INSTALLREPOSITORY}/XS/files_MCDataComparison/MC_CalibrationPM.root )
	}
	else {
	  sprintf(Command01,"${INSTALLREPOSITORY}/XS/HitEfficiency_WM -mw -i 1 -f 1000"); //Generate a file containing MC hit efficiency (XS/files_MCDataComparison/MC_CalibrationPM.root )
	  sprintf(Command02,"${INSTALLREPOSITORY}/XS/HitEfficiency_WM -r 29652 -t 29652 -i 0 -f 300 -w"); //Generate a file containing Data hit efficiency (${INSTALLREPOSITORY}/XS/files_MCDataComparison/MC_CalibrationPM.root )
	}
      }
      else if(ErrorType==4){
	if(PM){	
	  sprintf(Command03,"${INSTALLREPOSITORY}/XS/CompareCalibrationsPM -m -i 1 -f 1000"); //-> Generate a file containing each hit info for MC (XS/files_MCDataComparison/MC_CalibrationPM.root )
	  sprintf(Command04,"${INSTALLREPOSITORY}/XS/CompareCalibrationsPM -r 14510 -t 14570 -i 0 -f 300");//-> Generate a file containing each hit info for Data (XS/files_MCDataComparison/Data_CalibrationPM.root )
	  sprintf(Command05,"${INSTALLREPOSITORY}/XS/GeneratePEAngleDistributions");//Read the data and MC files above and create the dependency of PE with angle.
	}
	else{
	  sprintf(Command03,"${INSTALLREPOSITORY}/XS/CompareCalibrationsPM_WM -mw -i 1 -f 1000");
	  sprintf(Command04,"${INSTALLREPOSITORY}/XS/CompareCalibrationsPM_WM -r 29652 -t 29652 -i 0 -f 300 -w"); //tmp, all data have not been processed yet
	  sprintf(Command05,"${INSTALLREPOSITORY}/XS/GeneratePEAngleDistributions -w");
	}
      }
      if(ErrorType==5){
	// ML 2017/08 not needed anymore
	//	sprintf(Command06,"#Don't forget to copy your MC two times and change the Birks Constant. The path and name of this MC can be changed in XSFileGenerator.c");
      }

 
   
      for(int n=0;n<NE[ErrorType];n++){

	double ErrorValue=Start[ErrorType]+n*Step[ErrorType];
	if(MC){
       
	  for(int i=0;i<NMCfiles;i++){
	  //  for(int i=250;i<650;i++){//TEMP
	    sprintf(Command1,"");
	    sprintf(Command1_1,"");
	    sprintf(Command1_2,"");
	    sprintf(Command1_3,"");
	    sprintf(Command1_4,"");
	    sprintf(Command1_5,"");
	    sprintf(Command1_6,"");
	    sprintf(Command2,"");
	    sprintf(Command3,"");
	    sprintf(Command4,"");
	    sprintf(Command5,"");
	    sprintf(Command6,"");
	    sprintf(Command7,"");
	    sprintf(Command8,"");
	    sprintf(Command9,"");
	    sprintf(Command11,"");
	    sprintf(Command12,"");
	
	    sprintf(Command10,"");
	    sprintf(Command20,"");
	    sprintf(Command30,"");
	    //sprintf(Command40,"");

	    sprintf(optionTune,Form(fileTune,i));

	    if(ErrorType==0){//Case w/o error

	  
	      //################################################Create the total input file####################################################################################


	      //Run 1
	      if(ParticleGun){

		if(i<1000) sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/13a_nd2_numu_ch_%d.nt":"${MCINPUTSTORAGE_WM}/13a_nd2_numu_h2o_%d.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}

		sprintf(Command1,"");
		sprintf(Command1_1,"rm ${MCOUTPUTSTORAGE%s}/%sMC%s_Run1_%d.root \n hadd ${MCOUTPUTSTORAGE%s}/%sMC%s_Run1_%d.root",suffix,DetName,ParticleGenerator,i,suffix,DetName,ParticleGenerator,i);

		for(int ip=0;ip<npdg;ip++){
		  if(i==0) cout<<"pdg values tested="<<pdgValues[ip]<<endl;
		  sprintf(Command1+strlen(Command1),"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC%s%d_Numu_Run1_%d.root -i %s -m %d -f 1 -p %d \n",execMC,suffix,DetName,ParticleGenerator,pdgValues[ip],i,InputFile,MCDetID,pdgValues[ip]);
		  sprintf(Command1_1+strlen(Command1_1)," ${MCOUTPUTSTORAGE%s}/%sMC%s%d_Numu_Run1_%d.root",suffix,DetName,ParticleGenerator,pdgValues[ip],i);
		}
	      }
	      else{
		//NuMu
		// *** warning MC-output's filename changed cause I have only numu ***
		if(i<1000) sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/13a_nd2_numu_ch_%d.nt":"${MCINPUTSTORAGE_WM}/13a_nd2_numu_h2o_%d.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}
		sprintf(Command1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d.root -i %s -m %d -f 1",execMC,suffix,DetName,i,InputFile,MCDetID);
		
		//if(i<=500) sprintf(Command1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numu_Run1_%d.root -i ${MCINPUTSTORAGE}/run1/11bfluka_nd2_numu_ch_%d.nt -m %d -f 1",execMC,i,i,MCDetID);
		//else if(i<=1000) sprintf(Command1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numu_Run1_%d.root -i ${MCINPUTSTORAGE}/run1add1/11bfluka_nd2_numu_ch_%d.nt -m %d -f 1",execMC,i,i-500,MCDetID);
		//else if(i<=1500) sprintf(Command1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numu_Run1_%d.root -i ${MCINPUTSTORAGE}/run1add2/11bfluka_nd2_numu_ch_%d.nt -m %d -f 1",execMC,i,i-1000,MCDetID);
		//else sprintf(Command1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numu_Run1_%d.root -i ${MCINPUTSTORAGE}/run1add3/11bfluka_nd2_numu_ch_%d.nt -m %d -f 1",execMC,i,i-1500,MCDetID);
		
		//NuMuBar
		if(i<1000) sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/ingbg_5.3.6/13a_nd2_numubar_ch_%d.nt":"${MCINPUTSTORAGE_WM}/ingbg_5.3.6/13a_nd2_numubar_h2o_%d.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}
		sprintf(Command1_1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_Numubar_Run1_%d.root -i %s -m %d -f 2",execMC,suffix,DetName,i,InputFile,MCDetID);
		
		// if(i<=500) sprintf(Command1_1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numubar_Run1_%d.root -i ${MCINPUTSTORAGE}/numubar/run1/11bfluka_nd2_numubar_ch_%d.nt -m %d -f 2",execMC,i,i,MCDetID);
		//else if(i<=1000) sprintf(Command1_1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numubar_Run1_%d.root -i ${MCINPUTSTORAGE}/run1add1nubar/11bfluka_nd2_numubar_ch_%d.nt -m %d -f 2",execMC,i,i-500,MCDetID);
		//else if(i<=1500) sprintf(Command1_1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numubar_Run1_%d.root -i ${MCINPUTSTORAGE}/run1add2nubar/11bfluka_nd2_numubar_ch_%d.nt -m %d -f 2",execMC,i,i-1000,MCDetID);
		//else if(i<=2000) sprintf(Command1_1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Numubar_Run1_%d.root -i ${MCINPUTSTORAGE}/run1add3nubar/11bfluka_nd2_numubar_ch_%d.nt -m %d -f 2",execMC,i,i-1500,MCDetID);
		
		//NuE
		if(i<1000) sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/ingbg_5.3.6/13a_nd2_nue_ch_%d.nt":"${MCINPUTSTORAGE_WM}/ingbg_5.3.6/13a_nd2_nue_h2o_%d.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}
		sprintf(Command1_2,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_Nue_Run1_%d.root -i %s -m %d -f 3",execMC,suffix,DetName,i,InputFile,MCDetID);
		
		//if(i<=500) sprintf(Command1_2,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Nue_Run1_%d.root -i ${MCINPUTSTORAGE}/nuerun1/11bfluka_nd2_nue_ch_%d.nt -m %d -f 3",execMC,i,i,MCDetID);
		//else if(i<=1000) sprintf(Command1_2,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Nue_Run1_%d.root -i ${MCINPUTSTORAGE}/nuerun1add1/11bfluka_nd2_nue_ch_%d.nt -m %d -f 3",execMC,i,i-500,MCDetID);
		//else if(i<=1500) sprintf(Command1_2,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Nue_Run1_%d.root -i ${MCINPUTSTORAGE}/nuerun1add2/11bfluka_nd2_nue_ch_%d.nt -m %d -f 3",execMC,i,i-1000,MCDetID);
		//else if(i<=2000) sprintf(Command1_2,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Nue_Run1_%d.root -i ${MCINPUTSTORAGE}/nuerun1add3/11bfluka_nd2_nue_ch_%d.nt -m %d -f 3",execMC,i,i-1500,MCDetID);
		
		//Wall Bkg (mainly Sand Muons)
		if(i<1000)sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/wallbg_5.3.6/13a_nd7_numu_o_%d_1.nt":"${MCINPUTSTORAGE_WM}/wallbg_5.3.6/13a_nd7_numu_o_%d_1.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}	      
		sprintf(Command1_3,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_Wall_Run1_%d.root -i %s -m %d -f 1",execMC,suffix,DetName,i,InputFile,MCDetID);
		
		//if(i<=1000) sprintf(Command1_3,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Wall_Run1_%d.root -i ${MCINPUTSTORAGE}/wall/11b_nd7_numu_o_%d_0.nt -m %d -f 1",execMC,i,i,MCDetID);
		//else if(i<=2000) sprintf(Command1_3,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_Wall_Run1_%d.root -i ${MCINPUTSTORAGE}/wallset2/11b_nd7_numu_o_%d_1.nt -m %d -f 1",execMC,i,i-1000,MCDetID);
		
		//Ingrid Bkg
		//Horizontal modules
		if(i<1000) sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/ingbg_5.3.6/13a_nd3_numu_fe_%d.nt":"${MCINPUTSTORAGE_WM}/ingbg_5.3.6/13a_nd3_numu_fe_%d.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}	      
		sprintf(Command1_4,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_INGRIDH_Run1_%d.root -i %s -m %d -f 1",execMC,suffix,DetName,i,InputFile,MCDetID);
		
		//if(i<=500) sprintf(Command1_4,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDH_Run1_%d.root -i ${MCINPUTSTORAGE}/nd3numurun1/11bfluka_nd3_numu_fe_%d.nt -m %d -f 1",execMC,i,i,MCDetID);
		//else if(i<=1000) sprintf(Command1_4,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDH_Run1_%d.root -i ${MCINPUTSTORAGE}/nd3numurun1add1/11bfluka_nd3_numu_fe_%d.nt -m %d -f 1",execMC,i,i-500,MCDetID);
		//else if(i<=1500) sprintf(Command1_4,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDH_Run1_%d.root -i ${MCINPUTSTORAGE}/nd3numurun1add2/11bfluka_nd3_numu_fe_%d.nt -m %d -f 1",execMC,i,i-1000,MCDetID);
		//else if(i<=2000) sprintf(Command1_4,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDH_Run1_%d.root -i ${MCINPUTSTORAGE}/nd3numurun1add3/11bfluka_nd3_numu_fe_%d.nt -m %d -f 1",execMC,i,i-1500,MCDetID);
		
		//Vertical modules
		if(i<1000) sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/ingbg_5.3.6/13a_nd4_numu_fe_%d.nt":"${MCINPUTSTORAGE_WM}/ingbg_5.3.6/13a_nd4_numu_fe_%d.nt"),i);
		else {cout<<"*** only 1000 NEUT files available ***"<<endl; return 0;}	      
		sprintf(Command1_5,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_INGRIDV_Run1_%d.root -i %s -m %d -f 1",execMC,suffix,DetName,i,InputFile,MCDetID);
		
		//if(i<=500) sprintf(Command1_5,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDV_Run1_%d.root -i ${MCINPUTSTORAGE}/nd4numurun1/11bfluka_nd4_numu_fe_%d.nt -m %d -f 1",execMC,i,i,MCDetID);
		//else if(i<=1000) sprintf(Command1_5,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDV_Run1_%d.root -i ${MCINPUTSTORAGE}/nd4numurun1add1/11bfluka_nd4_numu_fe_%d.nt -m %d -f 1",execMC,i,i-500,MCDetID);
		//else if(i<=1500) sprintf(Command1_5,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDV_Run1_%d.root -i ${MCINPUTSTORAGE}/nd4numurun1add2/11bfluka_nd4_numu_fe_%d.nt -m %d -f 1",execMC,i,i-1000,MCDetID);
		//else if(i<=2000) sprintf(Command1_5,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE}/PMMC_INGRIDV_Run1_%d.root -i ${MCINPUTSTORAGE}/nd4numurun1add3/11bfluka_nd4_numu_fe_%d.nt -m %d -f 1",execMC,i,i-1500,MCDetID);

		//Merge Files Together

		sprintf(Command1_6,"${INSTALLREPOSITORY}/XS/FinalMCOutputMaker -a ${MCOUTPUTSTORAGE%s}/%sMC_Numu_Run1_%d.root -b ${MCOUTPUTSTORAGE%s}/%sMC_Numubar_Run1_%d.root -c ${MCOUTPUTSTORAGE%s}/%sMC_Nue_Run1_%d.root -d ${MCOUTPUTSTORAGE%s}/%sMC_Wall_Run1_%d.root -e ${MCOUTPUTSTORAGE%s}/%sMC_INGRIDH_Run1_%d.root -f ${MCOUTPUTSTORAGE%s}/%sMC_INGRIDV_Run1_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d.root",suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i);
	      }
	      //##########################################################################################################################


		// need for Loli_addcrosstalk for WM here: added.
	      if(!PM) sprintf(Command2,"${INSTALLREPOSITORY}/Reconstruction/appWM/Loli_addcrosstalk_slit -f ${MCOUTPUTSTORAGE_WM}/WMMC%s_Run1_%d.root -o ${MCOUTPUTSTORAGE_WM}/WMMC%s_Run1_%d_wXtalk.root",Sand,i,Sand,i);     	  
		
	      //sprintf(Command1,"bin/Linux-g++/IngMC -o ${MCOUTPUTSTORAGE}/PMMC_Run1_%d.root -i /export/scraid2/data/bquilain/neutfile_pm/11bfluka_nd2_numu_ch_%d.nt -m 2 -f 2",i,i);
	      sprintf(Command3,"${INSTALLREPOSITORY}/Reconstruction/app%s/IngAddNoisePMMC_new -f ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d%s.root -o ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d_wNoise.root %s",DetName,suffix,DetName,ParticleGenerator,Sand,i,(PM?"":"_wXtalk"),suffix,DetName,ParticleGenerator,Sand,i,(PM?"":Form("-n %2.1f -w",DNRateWM)));
	      sprintf(Command4,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d_wNoise.root -o ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d_wNoise_recon.root",DetName,execRecon,suffix,DetName,ParticleGenerator,Sand,i,suffix,DetName,ParticleGenerator,Sand,i);
	      sprintf(Command5,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d_wNoise_recon.root -o ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d_wNoise_ana.root %s",DetName, execAna,suffix,DetName,ParticleGenerator,Sand,i,suffix,DetName,ParticleGenerator,Sand,i,(ParticleGun || sandOnly?"-N":""));
		
	      //ML tmp
	      //sprintf(Command5,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh"); 

	      //NEUT 5.3.2  -> now NEUT 5.3.6 (ML 2017/02/01)

	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC%s%s_Run1_%d_wNoise_ana.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s%s%s_Run1_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,ParticleGenerator,Sand,i,DetName,ParticleGenerator,Sand,i,(PM?"":"w"),optionTune);

	      //NEUT 5.1.4.2
	      //sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan -i ${MCOUTPUTSTORAGE}/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_SciBar188_wNoise_KSana_woXtalk.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_Old_%d_Plan.root -f 1 -m",i,i);
	      
	      //sprintf(Command7,"${T2KREWEIGHTREPOSITORY}/app/genWeightsFromINGRID_2015.exe -i ${MCOUTPUTSTORAGE}/PMMC_Run1_%d.root -o ${MCOUTPUTSTORAGE}/PMMC_Run1_%d_ReWeight2015.root",i,i);

	    }
	    
	    else if(ErrorType==1) continue;
	    else if(ErrorType==2){//Dark noise values.
	      sprintf(Command1,"${INSTALLREPOSITORY}/Reconstruction/app%s/IngAddNoisePMMC_new -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_Systematics%d_%d.root -n %2.1f",DetName,suffix,DetName,i,suffix,DetName,i,ErrorType,n,ErrorValue);  
	      sprintf(Command2,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root",DetName,execRecon,suffix,DetName,i,ErrorType,n,suffix,DetName,i,ErrorType,n);
	      sprintf(Command3,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",DetName,execAna,suffix,DetName,i,ErrorType,n,suffix,DetName,i,ErrorType,n);
	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,i,ErrorType,n,DetName,i,ErrorType,n,(PM?"":"w"),optionTune);
	    }
	    else if(ErrorType==3){//One should give in input the file name of the data/MC difference
	      //The best thing would be just to generate everything w/ an option
	      sprintf(Command1,"${INSTALLREPOSITORY}/XS/GenerateHitEfficiencyMask%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_Systematics%d_%d.root %s",suffix,suffix,DetName,i,suffix,DetName,i,ErrorType,n,(PM?"":"-w"));
	      sprintf(Command2,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root",DetName,execRecon,suffix,DetName,i,ErrorType,n,suffix,DetName,i,ErrorType,n);
	      sprintf(Command3,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",DetName,execAna,suffix,DetName,i,ErrorType,n,suffix,DetName,i,ErrorType,n);
	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,i,ErrorType,n,DetName,i,ErrorType,n,(PM?"":"w"),optionTune);
	    }
	    else if(ErrorType==4){//One should give in input the file name of the data/MC difference
	      // for WM, input in LoliAna -- for PM, input inXS_CC0pi_Plan
	      if(!PM){
		sprintf(Command1,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v ${INSTALLREPOSITORY}/XS/files/PEXAngle_WM.root",DetName,execAna,suffix,DetName,i,suffix,DetName,i,ErrorType,n,ErrorType);
		sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,i,ErrorType,n,DetName,i,ErrorType,n,(PM?"":"w"),optionTune);
	      }
	      else {
		sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan -i ${MCOUTPUTSTORAGE}/PMMC_Run1_%d_wNoise_ana.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_PM_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m -e %d -v ${INSTALLREPOSITORY}/XS/files/PEXAngle_PM.root %s",i,i,ErrorType,n,ErrorType,optionTune);
	      }	
	    }
	    else if(ErrorType==5){
	      // the value of the Birks constant is an option
	      // -B (0,1,or 2) for (-1sigma, nominal, or +1sigma) [0.0185,0.0208,0.0231]
	      int birksindex=2*n;

	      // *** warning MC-output's filename changed cause I have only numu ***
	      //sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/ingbg_5.3.6/13a_nd2_numu_ch_%d.nt":"${MCINPUTSTORAGE_WM}/ingbg_5.3.6/13a_nd2_numu_h2o_%d.nt"),i);// NEUT 5.3.6
	      sprintf(InputFile,(PM?"${MCINPUTSTORAGE}/13a_nd2_numu_ch_%d.nt":"${MCINPUTSTORAGE_WM}/13a_nd2_numu_h2o_%d.nt"),i);
	      sprintf(Command1,"${INSTALLREPOSITORY}/MC/bin/Linux-g++/%s -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_Systematics%d_%d.root -i %s -m %d -f 1 -B %d",execMC,suffix,DetName,i,ErrorType,birksindex,InputFile,MCDetID,birksindex);

	      // re-run FinalOutputMaker with background files (unmodified) 
	      sprintf(Command1_6,"${INSTALLREPOSITORY}/XS/FinalMCOutputMaker -a ${MCOUTPUTSTORAGE%s}/%sMC_Numu_Run1_%d_Systematics%d_%d.root -b ${MCOUTPUTSTORAGE%s}/%sMC_Numubar_Run1_%d.root -c ${MCOUTPUTSTORAGE%s}/%sMC_Nue_Run1_%d.root -d ${MCOUTPUTSTORAGE%s}/%sMC_Wall_Run1_%d.root -e ${MCOUTPUTSTORAGE%s}/%sMC_INGRIDH_Run1_%d.root -f ${MCOUTPUTSTORAGE%s}/%sMC_INGRIDV_Run1_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_Systematics%d_%d.root",suffix,DetName,i,ErrorType,birksindex,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,suffix,DetName,i,ErrorType,birksindex);
	
	      if(!PM) sprintf(Command2,"${INSTALLREPOSITORY}/Reconstruction/appWM/Loli_addcrosstalk_slit -f ${MCOUTPUTSTORAGE_WM}/WMMC%s_Run1_%d_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE_WM}/WMMC%s_Run1_%d_wXtalk_Systematics%d_%d.root",Sand,i,ErrorType,birksindex,Sand,i,ErrorType,birksindex);     	  
	      sprintf(Command3,"${INSTALLREPOSITORY}/Reconstruction/app%s/IngAddNoisePMMC_new -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d%s_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_Systematics%d_%d.root %s",DetName,suffix,DetName,i,(PM?"":"_wXtalk"),ErrorType,birksindex,suffix,DetName,i,ErrorType,birksindex,(PM?"":Form("-n %2.1f -w",DNRateWM)));  
	      sprintf(Command4,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root",DetName,execRecon,suffix,DetName,i,ErrorType,birksindex,suffix,DetName,i,ErrorType,birksindex);
	      sprintf(Command5,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",DetName,execAna,suffix,DetName,i,ErrorType,birksindex,suffix,DetName,i,ErrorType,birksindex);
	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,i,ErrorType,birksindex,DetName,i,ErrorType,birksindex,(PM?"":"w"),optionTune);
	    }
	    else if(ErrorType==6){
	      // done at the level of CC0piSelection.c
	      // sprintf(Command1,"${INSTALLREPOSITORY}/XS/CC0piSelection -s 1 -w %2.2f",ErrorValue);
 	    }
	    else if(ErrorType>=7 && ErrorType<=10){
	      // modif of PMrecon - use the pre-existing files as an input
	      sprintf(Command1,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -e %d -v %3.1f",DetName,execRecon,suffix,DetName,i,suffix,DetName,i,ErrorType,n,ErrorType,ErrorValue);
	      sprintf(Command2,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon_Systematics%d_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root",DetName,execAna,suffix,DetName,i,ErrorType,n,suffix,DetName,i,ErrorType,n);
	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,i,ErrorType,n,DetName,i,ErrorType,n,(PM?"":"w"),optionTune);
	    }
	    else if(ErrorType>=11 && ErrorType<=15){ 
	      // modif of PMAna - use the pre-existing files as an input
	      sprintf(Command1,"${INSTALLREPOSITORY}/Reconstruction/app%s/%s -r 14000 -s 0 -f ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_recon.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -e %d -v %3.1f",DetName,execAna,suffix,DetName,i,suffix,DetName,i,ErrorType,n,ErrorType,ErrorValue);
	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_Systematics%d_%d_Plan.root -f 1 -m%s %s",suffix,suffix,DetName,i,ErrorType,n,DetName,i,ErrorType,n,(PM?"":"w"),optionTune);
	    }
	    else if(ErrorType==16) continue;
	    //sprintf(Command2,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan -i ${MCOUTPUTSTORAGE}/PMMC_Run1_%d_wNoise_ana.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_Run1_%d_Plan_ReWeight2015.root -f 1 -m",i,i);

	    else if(ErrorType>=17){
	      if(n!=0 || ErrorType>17) continue;//this is a something only for XS errors. This is because all the error is contained in one file (the files contains nominal tree + Reweight vector) and so, we do not need to generate different files for different sources or for 7 different variations of it (-3,-2,-1,0,1,2,3)sigmas.
	      // so just run at the first iteration
	      sprintf(Command2,"${T2KREWEIGHT}/app/genWeightsFromINGRID_2015ML.exe -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d.root -o ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_ReWeight2015.root",suffix,DetName,i,suffix,DetName,i); 
	      sprintf(Command6,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan%s -i ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_wNoise_ana.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s%s_Run1_%d_ReWeight2015_Plan.root -f 1 -m -x ${MCOUTPUTSTORAGE%s}/%sMC_Run1_%d_ReWeight2015.root %s",suffix,suffix,DetName,i,DetName,ParticleGenerator,i,suffix,DetName,i,(PM?"":"-w"));

	    }

	    sprintf(Name1,"%s/MC/Jobs/%sMC%d%s_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,i,Sand,ErrorType,n);	
	    sprintf(Name2,"%s/MC/Jobs/condor%sMC%d%s_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,i,Sand,ErrorType,n);
	    if(i==1) cout<<Name1<<" is created, and for other neutfiles also (number of neutfile = NMCFiles: see inc/setup.h"<<endl;
	    ofstream ScriptMC(Name1);
	    if(ScriptMC)
	      {
		if(sandOnly){

		ScriptMC<<"#!/bin/bash +x"<<endl
		  	<<Command1_3<<endl
			<<Command2<<endl
			<<Command3<<endl
			<<Command4<<endl
			<<Command5<<endl; 

		}
		else{

		  ScriptMC<<"#!/bin/bash +x"<<endl
		    // 		      <<Command1<<endl
		    /*  <<Command1_1<<endl
			  <<Command1_2<<endl
			  <<Command1_3<<endl
 			  <<Command1_4<<endl
			  <<Command1_5<<endl
			  <<Command1_6<<endl
			  <<Command2<<endl
		      <<Command3<<endl
		      <<Command4<<endl
		      <<Command5<<endl*/ 
			  <<Command6<<endl;
		  //		  <<Command7<<endl;
		}
	      }
      
	
	    sprintf(Command10,"Executable = %s/MC/Jobs/%sMC%d%s_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,i,Sand,ErrorType,n);
	    sprintf(Command20,"Output = %s/MC/Jobs/log_%sMC%d%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,i,Sand,ErrorType,n);
	    sprintf(Command30,"Error = %s/MC/Jobs/err_%sMC%d%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,i,Sand,ErrorType,n);

	    ofstream CondorMC(Name2);
	    if(CondorMC){
	      CondorMC<<Command10<<endl
		      <<"Universe        = vanilla"<<endl
		      <<"Rank            = kflops"<<endl
		      <<"Getenv          = True"<<endl
		      <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
		      <<"Arguments      =  ${Process)"<<endl
		      <<Command20<<endl
		      <<Command30<<endl
		      <<"Notification    = Never"<<endl
		      <<"QUEUE 1"<<endl;
	    }
	  }
	}
      
	if(Data){
	  ///////////////////////////////////////////DATA//////////////////////////////////////////////////////////////////////////
	  char RunNumber[16];
	  char SubRunNumber[16];
	  int RunNum;
	  int SubRunNum;

	  ifstream List;
	  for(int l=StartRunList;l<=EndRunList;l++){
	    sprintf(Name,"%s/Data/%d000.txt",cINSTALLREPOSITORY,l);
	    List.open(Name);
	    char MIDAS[150];
	    const int wm=(PM?0:8);
	    if(List.is_open()){
	      while(!List.eof()){
		List>>MIDAS;
		//cout<<MIDAS<<", 20e="<<MIDAS[50]<<endl;
		sprintf(RunNumber,"%c%c%c%c%c",MIDAS[50+wm],MIDAS[51+wm],MIDAS[52+wm],MIDAS[53+wm],MIDAS[54+wm]);
		sprintf(SubRunNumber,"%c%c%c%c",MIDAS[56+wm],MIDAS[57+wm],MIDAS[58+wm],MIDAS[59+wm]);
		RunNum=atoi(RunNumber);
		SubRunNum=atoi(SubRunNumber);
	
		if(RunNum<StartRun || RunNum>EndRun) continue;
		if(SubRunNum<StartSubRun || SubRunNum>EndSubRun) continue;
		cout<<RunNumber<<", SubRun="<<SubRunNumber<<endl;
	   
		//sprintf(Name1,"${INSTALLREPOSITORY}/MC/Jobs/PMData_%d_%d.sh",RunNum,SubRunNum);
		//sprintf(Name2,"${INSTALLREPOSITORY}/MC/Jobs/condorPMData%d_%d.sh",RunNum,SubRunNum);

		sprintf(Command1,"");
		sprintf(Command1_1,"");
		sprintf(Command1_2,"");
		sprintf(Command1_3,"");
		sprintf(Command1_4,"");
		sprintf(Command1_5,"");
		sprintf(Command1_6,"");
		sprintf(Command2,"");
		sprintf(Command3,"");
		sprintf(Command4,"");
		sprintf(Command5,"");
		sprintf(Command6,"");
		sprintf(Command7,"");
		sprintf(Command8,"");
		sprintf(Command9,"");
		sprintf(Command11,"");
		sprintf(Command12,"");
		sprintf(Command13,"");
		sprintf(Command14,"");
		sprintf(Command15,"");
		sprintf(Command16,"");
		sprintf(Command17,"");
		sprintf(Command18,"");
	
		sprintf(Command10,"");
		sprintf(Command20,"");
		sprintf(Command30,"");
		//sprintf(Command40,"");

		if(PM){
		  if(ErrorType==0){
		    sprintf(Command1,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/Calc_MPPC_new_ci_2.exe -r %d -s %d -t 1",RunNum,SubRunNum);
		    sprintf(Command2,"if [ -f ${DATAOUTPUTSTORAGE_CALIB}/ingrid_%08d_%04d_Calib00.root ]",RunNum,SubRunNum);
		    sprintf(Command3,"then");
		    sprintf(Command4,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/DSTMaker.exe -r %d -s %d -t 1",RunNum,SubRunNum);
		    sprintf(Command5,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/IngCalib_ADCNLCorrected.exe -r %d -s %d -f ${DATAOUTPUTSTORAGE_DST}/ingrid_%08d_%04d.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_calib.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command6,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/DSTMaker.exe -r %d -s %d -t 1 -p",RunNum,SubRunNum);
		    sprintf(Command7,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/IngCalib_ADCNLCorrected.exe -r %d -s %d -f ${DATAOUTPUTSTORAGE_ROOT}/ingrid_%08d_%04d.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmcalib.root -p",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command8,"${INSTALLREPOSITORY}/Reconstruction/appPM/IngMergePM -r %d -s %d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmcalib.root -a ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_calib.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmerged.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    ////////////////////////////////////////////////////////RECONSTRUCTION////////////////////////////////////////////////////////////
		    sprintf(Command9,"if [ -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmerged.root ]",RunNum,SubRunNum);
		    sprintf(Command11,"then");
		    sprintf(Command12,"${INSTALLREPOSITORY}/Reconstruction/appPM/PMreconRevOfficial -r %d -s %d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmerged.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    //////////////////////////////////////////////////////PM ANALYSIS////////////////////////////////////////////////////////////
		    sprintf(Command13,"${INSTALLREPOSITORY}/Reconstruction/appPM/PMAnaRevOfficial -r %d -s %d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command14,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngAddBSD -r %08d -s %04d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root -v -p",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command15,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh");
		    sprintf(Command16,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan -i ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%08d_%04d.root -f 1 -n 5",RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command17,"fi");
		    sprintf(Command18,"fi");
		  }
		  else if(ErrorType>=7 && ErrorType<=10){
		    sprintf(Command1,"${INSTALLREPOSITORY}/Reconstruction/appPM/PMreconRevOfficial -r %d -s %d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmerged.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk_Systematics%d_%d.root -e %d -v %3.1f",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,ErrorType,ErrorValue);
		    sprintf(Command2,"${INSTALLREPOSITORY}/Reconstruction/appPM/PMAnaRevOfficial -r %d -s %d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk_Systematics%d_%d.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk_Systematics%d_%d.root",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    sprintf(Command3,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngAddBSD -r %08d -s %04d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk_Systematics%d_%d.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk_Systematics%d_%d.root -v -p",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    sprintf(Command4,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh");
		    sprintf(Command5,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan -i ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%08d_%04d_Systematics%d_%d.root -f 1 -n 5",RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		  }
		  else if(ErrorType>=11 && ErrorType<=15){
		    sprintf(Command1,"${INSTALLREPOSITORY}/Reconstruction/appPM/PMAnaRevOfficial -r %d -s %d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSrecon_woXTalk_Systematics%d_%d.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk_Systematics%d_%d.root",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    sprintf(Command2,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngAddBSD -r %08d -s %04d -f ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPMana_woXTalk_Systematics%d_%d.root -o ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk_Systematics%d_%d.root -v -p",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    sprintf(Command3,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh");
		    sprintf(Command4,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan -i ${DATAOUTPUTSTORAGE_RECONSTRUCTED}/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%08d_%04d_Systematics%d_%d.root -f 1 -n 5",RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		  }
		  else continue;
		}
		else { //WM
		  if(ErrorType==0){
		    sprintf(Command1,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/Calc_MPPC_new.exe -r %d -s %d -t 1",RunNum,SubRunNum);
		    sprintf(Command2,"if [ -f ${DATAOUTPUTSTORAGE_WM_CALIB}/ingrid_%08d_%04d_Calib00.root ]",RunNum,SubRunNum);
		    sprintf(Command3,"then");
		    sprintf(Command4,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/DSTMaker_WM.exe -r %d -s %d -t 1 -o ${DATAOUTPUTSTORAGE_WM_DST}/ingrid_%08d_%04d.root",RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command5,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/DSTMaker_WM.exe -r %d -s %d -t 1 -w -o ${DATAOUTPUTSTORAGE_WM_ROOT}/ingrid_%08d_%04d.root",RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command6,"${INGRIDSOFTWARE_INGRID}/v1r1/Linux-x86_64/IngCalib_ADCNLCorrected.exe -r %d -s %d -f ${DATAOUTPUTSTORAGE_WM_DST}/ingrid_%08d_%04d.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_calib.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum); // use Ben's macro
		    sprintf(Command7,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngCalib_new_WM -r %d -s %d -f ${DATAOUTPUTSTORAGE_WM_ROOT}/ingrid_%08d_%04d.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_wmcalib.root -w",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command8,"${INSTALLREPOSITORY}/Reconstruction/appWM/IngMergePM -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_wmcalib.root -a ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_calib.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_wmmerged.root",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    ////////////////////////////////////////////////////////RECONSTRUCTION////////////////////////////////////////////////////////////
		    sprintf(Command9,"if [ -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_wmmerged.root ]",RunNum,SubRunNum);
		    sprintf(Command11,"then");
		    sprintf(Command12,"${INSTALLREPOSITORY}/Reconstruction/appWM/Lolirecon -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_wmmerged.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_recon.root",RunNum,SubRunNum,RunNum,SubRunNum);
		    //////////////////////////////////////////////////////WM ANALYSIS////////////////////////////////////////////////////////////
		    sprintf(Command13,"${INSTALLREPOSITORY}/Reconstruction/appWM/LoliAna -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_recon.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_ana.root",RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command14,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngAddBSD -r %08d -s %04d -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_ana.root -q ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_bsd.root -v -w",RunNum,SubRunNum,RunNum,SubRunNum,RunNum,SubRunNum);
		    //		    sprintf(Command15,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh");
		    sprintf(Command16,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan_WM -i ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_bsd.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%08d_%04d.root -f 1 -w",RunNum,SubRunNum,RunNum,SubRunNum);
		    sprintf(Command17,"fi");
		    sprintf(Command18,"fi");
		  }
		  else if(ErrorType>=7 && ErrorType<=10){
		    sprintf(Command12,"${INSTALLREPOSITORY}/Reconstruction/appWM/Lolirecon -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_wmmerged.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_recon_Systematics%d_%d.root -e %d -v %3.1f",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,ErrorType,ErrorValue);
		    sprintf(Command13,"${INSTALLREPOSITORY}/Reconstruction/appWM/LoliAna -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_recon_Systematics%d_%d.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_ana_Systematics%d_%d.root",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    sprintf(Command14,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngAddBSD -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_ana_Systematics%d_%d.root -q ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_bsd_Systematics%d_%d.root -v -w",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    //sprintf(Command15,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh");
		    sprintf(Command16,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan_WM -i ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_bsd_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%08d_%04d_Systematics%d_%d.root -f 1 -w",RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		  }
		  else if(ErrorType>=11 && ErrorType<=15){
		    sprintf(Command13,"${INSTALLREPOSITORY}/Reconstruction/appWM/LoliAna -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_recon.root -o ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_ana_Systematics%d_%d.root",RunNum,SubRunNum,RunNum,SubRunNum,ErrorType,n);
		    sprintf(Command14,"${INGRIDSOFTWARE_INGRIDFORMAT}/app/IngAddBSD  -f ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_ana_Systematics%d_%d.root -q ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_bsd_Systematics%d_%d.root -v -w",RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		    //sprintf(Command15,"source ${INSTALLREPOSITORY}/source_T2KReweight.sh");
		    sprintf(Command16,"${INSTALLREPOSITORY}/XS/XS_CC0pi_Plan_WM -i ${DATAOUTPUTSTORAGE_WM_RECONSTRUCTED}/ingrid_%08d_%04d_bsd_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%08d_%04d_Systematics%d_%d.root -f 1 -w",RunNum,SubRunNum,ErrorType,n,RunNum,SubRunNum,ErrorType,n);
		  }
		  else continue;
		}


		sprintf(Name1,"%s/MC/Jobs/%sData_Run%d_SubRun%d_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,RunNum,SubRunNum,ErrorType,n);
		sprintf(Name2,"%s/MC/Jobs/condor%sData_Run%d_SubRun%d_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,RunNum,SubRunNum,ErrorType,n);
		cout<<Name1<<endl;
		ofstream ScriptData(Name1);
		if(ScriptData)
		  {
		    ScriptData<<"#!/bin/bash +x"<<endl
			      <<Command1<<endl
			      <<Command2<<endl
			      <<Command3<<endl
			      <<Command4<<endl
			      <<Command5<<endl
			      <<Command6<<endl
			      <<Command7<<endl
			      <<Command8<<endl
			      <<Command9<<endl
			      <<Command11<<endl
			      <<Command12<<endl
			      <<Command13<<endl
			      <<Command14<<endl
			      <<Command15<<endl
			      <<Command16<<endl
			      <<Command17<<endl
			      <<Command18<<endl;
		  }
		
		sprintf(Command10,"Executable = %s/MC/Jobs/PMData_Run%d_SubRun%d_Systematics%d_%d.sh",cINSTALLREPOSITORY,RunNum,SubRunNum,ErrorType,n);
		sprintf(Command20,"Output = %s/MC/Jobs/log_PMData_Run%d_SubRun%d_Systematics%d_%d.txt",cINSTALLREPOSITORY,RunNum,SubRunNum,ErrorType,n);
		sprintf(Command30,"Error = %s/MC/Jobs/err_PMData_Run%d_SubRun%d_Systematics%d_%d.txt",cINSTALLREPOSITORY,RunNum,SubRunNum,ErrorType,n);
		ofstream CondorData(Name2);
		if(CondorData){
		  CondorData<<Command10<<endl
			    <<"Universe        = vanilla"<<endl
			    <<"Rank            = kflops"<<endl
			    <<"Getenv          = True"<<endl
			    <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
			    <<"Arguments      =  ${Process)"<<endl
			    <<Command20<<endl
			    <<Command30<<endl
			    <<"Notification    = Never"<<endl
			    <<"QUEUE 1"<<endl;
		}	
	      }
	    }
	    List.close();
	  }
	  //cout<<endl;
	}
      }
    }
  }
  
  
  ///////////////////////////////////PREREQUISITE//////////////////////////////////////////
  sprintf(Name1,"%s/MC/Jobs/Prerequisite%s.sh",cINSTALLREPOSITORY,DetName);
  sprintf(Name2,"%s/MC/Jobs/condorPrerequisite%s.sh",cINSTALLREPOSITORY,DetName);
  //sprintf(Name1,"$(INSTALLREPOSITORY)/MC/Jobs/Prerequisite.sh");
  //sprintf(Name2,"$(INSTALLREPOSITORY)/MC/Jobs/condorPrerequisite.sh");
  cout<<Name1<<endl;
  ofstream Script(Name1);
  if(Script)
    {
      Script<<"#!/bin/bash +x"<<endl
	    <<Command01<<endl
	    <<Command02<<endl
	    <<Command03<<endl
	    <<Command04<<endl
	    <<Command05<<endl;
    }
  
  
  sprintf(Command10,"Executable = %s/MC/Jobs/Prerequisite%s.sh",cINSTALLREPOSITORY,DetName);
  sprintf(Command20,"Output = %s/MC/Jobs/log_Prerequisite%s.txt",cINSTALLREPOSITORY,DetName);
  sprintf(Command30,"Error = %s/MC/Jobs/err_Prerequisite%s.txt",cINSTALLREPOSITORY,DetName);
  ofstream CondorPre(Name2);
  if(CondorPre){
    CondorPre<<Command10<<endl
	     <<"Universe        = vanilla"<<endl
	     <<"Rank            = kflops"<<endl
	     <<"Getenv          = True"<<endl
	     <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
	     <<"Arguments      =  ${Process)"<<endl
	     <<Command20<<endl
	     <<Command30<<endl
	     <<"Notification    = Never"<<endl
	     <<"QUEUE 1"<<endl;
  }



  
  //////////////////////////////////////POSTREQUISITE////////////////////////////////////////////////. Ne prends pas en compte les Xsec pour l'instant
  char optionWM[2];sprintf(optionWM,(PM?"":"-W"));
  for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
    for(int n=0;n<NE[ErrorType];n++){
      double ErrorValue=Start[ErrorType]+n*Step[ErrorType];
      ///////////FOR THE POSTREQUISITE, SO SHOULD CLEAN FOR EACH ERROR/////////////////////
      sprintf(Command001,"");
      sprintf(Command002,"");
      sprintf(Command003,"");
      sprintf(Command004,"");
      sprintf(Command005,"");
      sprintf(Command006,"");

      sprintf(Name1,"%s/MC/Jobs/Postrequisite%s_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,ErrorType,n);
      sprintf(Name2,"%s/MC/Jobs/condorPostrequisite%s_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,ErrorType,n);
      cout<<Name1<<endl;
      ofstream ScriptPost(Name1);
      if(ErrorType==0){
	sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s%s_Run1_%s_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f %s %s",DetName,ParticleGenerator,"%d",DetName,ParticleGenerator,ErrorType,n,channel,ErrorType,ErrorValue,DataPOT,optionWM,(retuned?"-t 146":""));
	sprintf(Command002,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_%s.root -o ${INSTALLREPOSITORY}/XS/files/DataSelected_Systematics%d_%d -s %d -e %d -v %3.3f %s","%08d","%04d",ErrorType,n,channel,ErrorType,ErrorValue,optionWM);
	sprintf(Command004,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/DataSelected_%s%s_Systematics%d_%d -m ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics%d_%d -o ${INSTALLREPOSITORY}/XS/files/DataUnfolded_%s%s_Systematics%d_%d -n 1",DetName,ParticleGenerator,ErrorType,n,DetName,ParticleGenerator,ErrorType,n,DetName,ParticleGenerator,ErrorType,n);
      }
      else if(ErrorType==6) sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%s_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_%s_Systematics%d_%d -s %d -e %d -v %3.3f -m -w %2.2f -p %2.2f %s %s",DetName,"%d",DetName,ErrorType,n,channel,ErrorType,ErrorValue,ErrorValue,DataPOT,optionWM,(retuned?"-t 146":""));
      else if(ErrorType==5) sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%s_Systematics%d_%d_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_%s_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f %s %s",DetName,"%d",ErrorType,2*n,DetName,ErrorType,n,channel,ErrorType,ErrorValue,DataPOT,optionWM,(retuned?"-t 146":""));
      else sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%s_Systematics%d_%d_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_%s_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f %s %s",DetName,"%d",ErrorType,n,DetName,ErrorType,n,channel,ErrorType,ErrorValue,DataPOT,optionWM,(retuned?"-t 146":""));
      if(ErrorType>=7 && ErrorType<=15){
	sprintf(Command002,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_%s_Systematics%d_%d.root -o ${INSTALLREPOSITORY}/XS/files/DataSelected_Systematics%d_%d -s %d -e %d -v %3.3f %s","%08d","%04d",ErrorType,n,ErrorType,n,channel,ErrorType,ErrorValue,optionWM);
	sprintf(Command004,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/DataSelected_%s%s_Systematics%d_%d -m ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics%d_%d -o ${INSTALLREPOSITORY}/XS/files/DataUnfolded_%s%s_Systematics%d_%d -n 1",DetName,ParticleGenerator,ErrorType,n,DetName,ParticleGenerator,ErrorType,n,DetName,ParticleGenerator,ErrorType,n);
	//sprintf(Command004,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/DataSelected_Systematics%d_%d -m ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics%d_%d -o ${INSTALLREPOSITORY}/XS/files/DataUnfolded_Systematics%d_%d -n 1",ErrorType,n,ErrorType,n,ErrorType,n);
      } 
      else if(ErrorType==16){
	sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s%s_Run1_%s_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f %s %s",DetName,ParticleGenerator,"%d",DetName,ParticleGenerator,ErrorType,n,channel,ErrorType,ErrorValue,DataPOT,optionWM,(retuned?"-t 146":"")); 
	cout <<"Flux" << endl;
	//sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_Run1_%s_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f","%d",ErrorType,n,ErrorType,ErrorValue,DataPOT);	
      }
      else if(ErrorType>=17){
	//double XsecVariation=ErrorValue-CenterXsecVariations*(ErrorType-Systematics_Xsec_Start);//The variation of Xsec parameter, in #sigma. A number between 0 and 175 - the center of the current systematic source (nominal). For example, for Xsec error source #10, it starts from 7*(10-1)=63 and ends at 70. from 63 to 70, it contains the variariation of -3,-2,-1,0,1,2,3 sigma respectively. The center is then located at 66. For the example of a 2 sigma variation, the substraction will be therefore equal to: 68-66=2, which gives the number of sigmas!
	sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s%s_Run1_%s_ReWeight2015_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f %s %s",DetName,ParticleGenerator,"%d",DetName,ParticleGenerator,ErrorType,n,channel,ErrorType,ErrorValue,DataPOT,optionWM,(retuned?"-t 146":""));
	//sprintf(Command001,"${INSTALLREPOSITORY}/XS/CC0piSelection -i ${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%s_ReWeight2015_Plan.root -o ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics%d_%d -s %d -e %d -v %3.3f -m -p %2.2f","%d",DetName,ErrorType,n,ErrorType,ErrorValue,DataPOT);
      }
      
      sprintf(Command003,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics0_0 -m ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics%d_%d -o ${INSTALLREPOSITORY}/XS/files/MCUnfolded_%s%s_Systematics%d_%d -n 1",DetName,ParticleGenerator,DetName,ParticleGenerator,ErrorType,n,DetName,ParticleGenerator,ErrorType,n);
      if(ErrorType==1) sprintf(Command003,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics0_0 -m ${INSTALLREPOSITORY}/XS/files/MCSelected_%s%s_Systematics0_0 -o ${INSTALLREPOSITORY}/XS/files/MCUnfolded_%s%s_Systematics%d_%d -n 1 -s %d",DetName,ParticleGenerator,DetName,ParticleGenerator,DetName,ParticleGenerator,ErrorType,n);
      
      //sprintf(Command003,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics%d_%d -m ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics%d_%d -o ${INSTALLREPOSITORY}/XS/files/MCUnfolded_Systematics%d_%d -n 1",ErrorType,n,ErrorType,n,ErrorType,n);// Not for stat. variations, since unfolded is already done differently
      //if(ErrorType==1) sprintf(Command003,"${INSTALLREPOSITORY}/XS/UnfoldingOptimisation_Dvt -d ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics0_0 -m ${INSTALLREPOSITORY}/XS/files/MCSelected_Systematics0_0 -o ${INSTALLREPOSITORY}/XS/files/MCUnfolded_Systematics%d_%d -n 1 -s 1",ErrorType,n);

      
      if(ScriptPost)
	{
	  ScriptPost<<"#!/bin/bash +x"<<endl
		    <<Command001<<endl;
	  /*<<Command002<<endl
		    <<Command003<<endl
		    <<Command004<<endl;*/
	}
      sprintf(Command10,"Executable = %s/MC/Jobs/Postrequisite%s_Systematics%d_%d.sh",cINSTALLREPOSITORY,DetName,ErrorType,n);
      sprintf(Command20,"Output = %s/MC/Jobs/log_Postrequisite%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
      sprintf(Command30,"Error = %s/MC/Jobs/err_Postrequisite%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);

      ofstream CondorPost(Name2);
      if(CondorPost){
	CondorPost<<Command10<<endl
		  <<"Universe        = vanilla"<<endl
		  <<"Rank            = kflops"<<endl
		  <<"Getenv          = True"<<endl
		  <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
		  <<"Arguments      =  ${Process)"<<endl
		  <<Command20<<endl
		  <<Command30<<endl
		  <<"Notification    = Never"<<endl
		  <<"QUEUE 1"<<endl;
	
      }
    }
  }
  


  
}

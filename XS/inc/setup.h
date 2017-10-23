#ifndef setup_h
#define setup_h
/************CLibs********/
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;

#include <iomanip>
#include <sys/stat.h>
#include <cmath>
/********ROOT Libs*******/
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <THStack.h>
/**********INGRID Libs*************/
#include "TApplication.h"
#include "IngridConstants.h"
#include "INGRID_Dimension.cc"
/*#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "INGRID_Dimension.cc"
/*********Hit Libs****************/
/*
#include "Hit.h"
#include "Reconstruction.h"
*/

double INGRIDSCIBAR=2.0;
Int_t Scyc =  4;
Int_t Ncyc = 12;
int NMod=17;
int NPln=11;
int NPlnPM=18;
int NPlnWM=24;
int NView=2;
int NCh=24;
int NChPM=32;
int NChWM=40;
double IronCarbonRatio=7.87/1.03;
double DNRateWM=2.3;
double MeV2PE=46.0;
double MeV2PEPM=38.6;
double Beam[3]={0,-TMath::Sin(3.8*TMath::Pi()/180.),TMath::Cos(3.8*TMath::Pi()/180.)};
double C[18]={-1.85,0.34,0.59,0.74,0.514,-0.37,1.25,-0.06,0.562,0.82,-0.47,0.6,-0.57,-0.45,0,0,0,0};
#ifndef __INGRID_CONSTANTS__
const float Mmu=105.66/1000.;
const float Mp=938.272/1000.;
const float Mn=939.5659/1000.;
#endif
const float Mpi=139.57/1000.;
float EpiBins[5]={0,0.25,0.50,0.75,1};
float EpiReWeight[5]={0.135,0.4,0.294,1.206,1};
const int StartRun=29596;
const int EndRun=29760; //up to Dec 2016
const int StartSubRun=0;
const int EndSubRun=300;//300
const int StartRunList=29;//14 for PM; Which list will be read (29=>29000.txt). Necessary to use only run processed
const int EndRunList=30;//15 for PM; Which list will be read (29=>29000.txt). Necessary to use only run processed
double NMCfiles=500; // up to now only 1000 are available with NEUT 5.3.6

double DataPOTPM=0.58;//In units of 10^21 POT -- runs 234
//double DataPOTPM=0.76;//In units of 10^21 POT -- runs 234  (+56)
double DataPOTWM=0.72;//In units of 10^21 POT -- run 8
double DataPOT;
const int StartError=17;
const int EndError=37;//15 for det, 16 for det+flux, 36 for det+flux+xs, 37 to add NEUT tunings
const int MaxError=37;
int NFluxFiles;
int StartXsec=0;const int EndXsec=19;int NXsecVariations=7; int NXsecTunings=7;
int CenterXsecVariations=(int) (NXsecVariations-1-((double) (NXsecVariations-1)/2));

int NE[MaxError+1];
double Step[EndError+1];
double Start[EndError+1];
double End[EndError+1];
double Nominal; double Err;
const int Systematics_Detector_Start=2;
const int Systematics_Detector_End=15;
const int Systematics_Flux_Start=16;
const int Systematics_Flux_End=16;
const int Systematics_Xsec_Start=17;
const int Systematics_Xsec_End=Systematics_Xsec_Start+EndXsec;
bool EStatistics=true;//if true, estimate stat. error after unfolding
const int NStatisticalVariations=1000;//number of stat. varied toy experiments to evaluate the stat. error.
    
const int NSamples = 6;//number of track samples, see Reconstruction.cc
const int LimitTracks = 5;
const int LimitHits = 25;
double RangeRelativeDistance = 1.;
const int NDials=147; //NXsecVariations*(EndXsec+1) + NXsecTunings
//For particle gun
const int npdg=4;
int pdgValues[npdg]={13,211,-211,2212};


const int NBinsEnergyFlux=20;
double BinningEnergyFlux[NBinsEnergyFlux+1];
const int NBinsTrueEnergy=6;
const int NBinsRecEnergy=6;
double BinningTrueEnergy[NBinsTrueEnergy+1];
double BinningRecEnergy[NBinsRecEnergy+1];
const int NFSIs=13;//cc0pi+0p,cc0pi+1p,cc0pi+morep,cc1pi,cc1pi0,ccother,nc,+all bkg
const int NIntTypes=6;//ccqe,ccmec,ccres,cccoh,ccres,other
const int NBinsTrueMom=5;
const int NBinsTrueAngle=5;
const int NBinsRecMom=11;// was 17 -- ML 2017/08/03
const int NBinsRecAngle=15;// was 30 -- ML 2017/08/03
double BinningTrueMom[NBinsTrueMom+1];
double BinningTrueAngle[NBinsTrueAngle+1];
double BinningRecMom[NBinsRecMom+1];
double BinningRecAngle[NBinsRecAngle+1];

/////////////////// INITIALIZE ERRORS NOW ////////////////////////    
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
//16:Flux
//17: Xsec

void InitializeGlobal(bool PM=true){

  if(!PM)  Initialize_INGRID_Dimension();
  DataPOT=(PM?DataPOTPM:DataPOTWM);

  for(int i=0;i<=NBinsEnergyFlux;i++){
    if(i<=15) BinningEnergyFlux[i]=i*0.2;//in GeV
    else if(i<=16) BinningEnergyFlux[i]=BinningEnergyFlux[15]+(i-15)*1;
    else if(i<=19) BinningEnergyFlux[i]=BinningEnergyFlux[16]+(i-16)*2; // ML corr 2017/10
    else if(i<=20) BinningEnergyFlux[i]=30.;
    else cout<<"Error in binning the flux in energy. Please look at setup.h"<<endl;
    /*
    if(i==0) BinningEnergyFlux[i]=0;
    else if(i>=1 && i<36) BinningEnergyFlux[i]=0.5+(i-1)*0.1;
    else if(i>=36 && i<42) BinningEnergyFlux[i]=4+(i-36)*1;
    BinningEnergyFlux[42]=10; BinningEnergyFlux[NBinsEnergyFlux]=30;*/
  }


  for(int i=0;i<NBinsTrueMom+1;i++){
    BinningTrueMom[0]=0;
    BinningTrueMom[1]=0.5;
    BinningTrueMom[2]=0.7;
    BinningTrueMom[3]=1.0;
    BinningTrueMom[4]=5.0;
    BinningTrueMom[5]=30.0;
  }

  for(int i=0;i<NBinsRecMom+1;i++){
    /* BinningRecMom[0]=0;
    BinningRecMom[1]=10;    
    if(i>1 && i<NBinsRecMom) BinningRecMom[i]=10+5*(i-1);
    if(i==NBinsRecMom) BinningRecMom[i]=150;*/
    if(i<=6) BinningRecMom[i]=10*i;
    else if(i<=9) BinningRecMom[i]=60+5*(i-6);
  }
  BinningRecMom[10]=85;
  BinningRecMom[11]=150;

  BinningTrueAngle[0]=0;
  BinningTrueAngle[1]=10;
  BinningTrueAngle[2]=20;
  BinningTrueAngle[3]=30;
  BinningTrueAngle[4]=60;
  BinningTrueAngle[5]=180;

  for(int i=0;i<NBinsRecAngle+1;i++){
    //BinningRecAngle[i]=3*i;
    BinningRecAngle[i]=6*i;
  }
    

  for(int n=0;n<=37;n++){
    NE[n]=1;
    Start[n]=1;
    Step[n]=1;

    if(Systematics_Detector_End>=Systematics_Flux_Start || Systematics_Flux_Start>=Systematics_Xsec_Start) cout<<"################################################################################################################### STOP THE PROCESS, THERE IS A PROBLEM IN THE NUMBERING OF YOUR SYST. ERROR SOURCES. CHECK SYSTEMATICS_DETECTOR_START/END...###########################################################################"<<endl;
	 
    if(n>=Systematics_Detector_Start && n<=Systematics_Detector_End){
      if(n==1){//1: Stat
	Start[n]=0;
	NE[n]=1000;
	Step[n]=1;
      }
      else if(n==2){//2: DN
	Start[n]=0;
	NE[n]=10;
	if(PM)Step[n]=1;
	else Step[n]=0.5;
      }
      else if(n==3){//3: Hit efficiency
	NE[n]=1;//Only one probability to mask the channel, given by data - MC. So we try to recover data inefficiency from MC distribution. So only one variation, the one needed to recover data.
	//NE[n]=10;//TEMP
      }
      else if(n==4){//4: absolute difference MC/data light yield with track angle
	NE[n]=1;
      }
      else if(n==5){
	NE[n]=2;//Birks constant variations; 0 is BirksMinus, 1 is BirksPlus
	Start[n]=0;
	Step[n]=1;
      }
      else if(n==6){//6: Beam related background (in fact, this mainly evaluate sand muons)
	Nominal=0.6;
	Err=TMath::Sqrt(Nominal*Nominal+0.2*0.2+0.2*0.2);
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==7){//7: 2D reconstruction error - nb of active planes required
	Nominal=3;
	Err=1;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=5;
      }
      else if(n==8){    //8: VetoUpstreamCriteria: nominal=0 planes, vary from 0->2 per 1 plane step
	// numero of the veto plane in INGRID and PM
	// for the WM, it is defined in Lolirecon.hxx
	Nominal=0;
	Err=1;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==9){//9: VetoEdgeCriteria: nominal=80mm, vary 60->100 per 20mm steps
	//maximal distance of a hit to a track to match them
	//matching with hits in the side veto planes; don't exist for WM
	Nominal=80;
	Err=20;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==10){//10: FVCriteria: nominal=100cm, vary 50->100 per 10 cms steps
	if(PM){
	  Nominal=100;
	  Err=10;
	  Start[n]=Nominal-5*Err;
	  Step[n]=Err;
	  NE[n]=6;
	}
	else { // for WM, nominal = 80, vary 50->80 per 10 cm steps
	  Nominal=80;
	  Err=10;
	  Start[n]=Nominal-3*Err;
	  Step[n]=Err;
	  NE[n]=4;
	}
      }
      else if(n==11){//11: Vertexing, plane tolerance: nominal=2, vary 2->4 per 1 plane steps
	// for WM: nominal is 7, vary 4->10 per 3 planes step (3 is for the pln-grid-grid structure)
	if(PM){
	  Nominal=2;
	  Err=1;
	  Start[n]=Nominal;
	  Step[n]=Err;
	  NE[n]=3;
	}
	else {
	  Nominal=7;
	  Err=3;
	  Start[n]=Nominal-Err;
	  Step[n]=Err;
	  NE[n]=3;
	}
      }
      else if(n==12){//12: Vertexing, transverse tolerance: nominal=15cm, vary 15cm->20cm per 2.5cm steps
	// should be given in mm
	Nominal=150;
	Err=25;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==13){//13: Track matching, plane tolerance: nominal=4, vary 3->5 per 1 plane steps
	// for WM: nominal is 10, vary 4->10 per 3 planes step (3 is for the pln-grid-grid structure)
	if(PM){
	  Nominal=4;
	  Err=1;
	  Start[n]=Nominal-Err;
	  Step[n]=Err;
	  NE[n]=3;
	}
	else {
	  Nominal=10;
	  Err=3;
	  Start[n]=Nominal-2*Err;
	  Step[n]=Err;
	  NE[n]=3;
	}
      }
      else if(n==14){//14: INGRID/PM tracks angle matching: nominal=35°, vary 30°->40° per 5° steps
	Nominal=35;
	Err=5;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==15){//15: INGRID/PM tracks transverse position matching: nominal=8.5cm, vary 7.5cm->9.5 per 1cm steps
	// should be given in mm
	Nominal=85;
	Err=10.;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
    }
    else if(n>=Systematics_Flux_Start && n<=Systematics_Flux_End){//16: flux error
      Nominal=0;
      Err=1;
      Start[n]=Nominal;
      Step[n]=Err;
      NE[n]=500;
      NFluxFiles=NE[n];
    }
    else if(n>=Systematics_Xsec_Start && n<=Systematics_Xsec_End){//17: cross section error
      //one should separate the place in the reweight vector, and the real starts of the Xsection.
      //Real start=0. Place in the Xsection vector is (StartXsec+n)*NXsecVariations.
      Start[n]=(StartXsec+(n-Systematics_Xsec_Start))*NXsecVariations;
      Step[n]=1;
      NE[n]=NXsecVariations;
    }
    else if(n==Systematics_Xsec_End+1){ // tunings
      Start[n]=(StartXsec+(n-Systematics_Xsec_Start))*NXsecVariations;
      Step[n]=1;
      NE[n]=NXsecTunings;
    }
    End[n]=Start[n]+(NE[n]+1)*Step[n];
  }
}

/*
const int NEnergyBins=43;
double EnergyBins[NBinsEnergyFlux+1];
const int NInteractionTypes=11;
const int NEnergyBins=43;
const int NBinsEnu=6;
const int NBinsIron=17;
const int NBinsMom=5;
const int NBinsAngle=4;
const int NBinsAngleRec=4;
double DistIronBin[NBinsIron+1];
double DistMomBin[NBinsMom+1];
double DistAngleBin[NBinsAngle+1];
double DistAngleRecBin[NBinsAngle+1];
const int NBinsTrueAngle=90;
double DistTrueAngleBin[NBinsAngle+1];
const int _NBinsMom=5;//2014/11/15
//const int _NBinsMom=10;
const int _NBinsIron=17;
double _DistIronBin[_NBinsIron+1];
double _DistMomBin[_NBinsMom+1];
const int _NBinsAngle=10;
double _DistAngleBin[_NBinsAngle+1];
const int _NBinsTrueAngle=10;
double _DistTrueAngleBin[_NBinsAngle+1];
*/


// ML for Production 2 -- files NEUT 5.3.3 from Koga-san -- 2017/08/02
const int NBadMCFilesPM=40;
const int NBadMCFilesWM=28;
const int BadMCFilesPM[NBadMCFilesPM]={0,1,15,34,74,141,150,190,189,257,277,278,284,292,318,324,335,386,416,462,594,626,639,751,792,833,851,852,880,886,909,93,174,293,436,518,785,861,903,992};
const int BadMCFilesWM[NBadMCFilesWM]={0,14,19,80,176,189,190,284,415,462,594,655,666,673,743,792,814,988,93,174,293,436,518,544,785,861,903,992};

bool isBadFile(int ifile, bool PM){
  if(PM){
    for(int i=0;i<NBadMCFilesPM;i++)
      if(BadMCFilesPM[i]==ifile) return true;
  }
  else{
    for(int i=0;i<NBadMCFilesWM;i++)
      if(BadMCFilesWM[i]==ifile) return true;
  }
  return false;
}

// number of good files in the range [ifile,ffile-1]
int NGoodFiles(int ifile,int ffile,bool PM){
  int good=ffile-ifile;
  for(int i=ifile;i<ffile;i++)
    if(isBadFile(i,PM)) good--;
  return good;
}


#endif

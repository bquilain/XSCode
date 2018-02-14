//#check if error value for xsec is really -3,-2,-1,0,1,2,3
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
//#include <pair>
#include <sys/stat.h>
#include <cmath>
#include <TMinuit.h>
#include <TFitter.h>
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
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TFrame.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <THStack.h>
#include <TBox.h>
#include "setup.h"
#include "Reconstruction.cc"
#include "Xsec.cc"
//Step 0. Select here if we wish to produce the error table and plots on the reconstructed distributions or on the unfolded distributions (comment #define RECONSTRUCTED)
//#define RECONSTRUCTED
#define TEMP
//#define DEBUG
#define XSTABLE
//3 modes should not be mixed:
//1. Predictions MC: (Varied MC - varied bkg) - (Fixed MC - fixed bkg)
//2. Error to scale to data: (Fixed MC - varied bkg) - (Fixed MC - fixed bkg)
//3. Unfolded error: the way I constructed the code provide us:
//U.(Fixed MC - varied bkg) - U.(Fixed MC - fixed bkg)
// So, if we wish to evaluate the variation of MC prediction: this is the mode 0. For the real systematics, select the mode 1
const int MODE=1;

void Evaluate1DError(TH2D * NominalMC, double **** CovarianceMatrix, TH1D * NominalMC_RecMom, TH1D * NominalMC_RecAngle){

  int NBinsMom = NominalMC->GetNbinsX();
  int NBinsAngle = NominalMC->GetNbinsX();
  
  //First, create the momentum histogram value and error
  for(int e0=0;e0<NBinsMom;e0++){
    double Value=0;
    for(int e1=0;e1<NBinsAngle;e1++) Value+=NominalMC->GetBinContent(e0+1,e1+1);
    double ErrorSquared=0;
    
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	ErrorSquared+= ((double) CovarianceMatrix[e0][e1][e0][f1]);//Case of momentum histogram
      }
    }
    double Error=TMath::Sqrt(ErrorSquared);
    double RelativeError=Error / ( Value == 0 ? 1 : Value );
    NominalMC_RecMom->SetBinContent(e0+1,RelativeError);
    cout<<setprecision(3)<<"Momentum = " << NominalMC_RecMom->GetBinCenter(e0+1) << ", Error="<<100*RelativeError<<"%"<<endl;
    //NominalMC_RecMom->SetBinError(e0+1,Error);
  }

  //Second, create the momentum histogram value and error
  for(int e0=0;e0<NBinsAngle;e0++){
    double Value=0;
    for(int e1=0;e1<NBinsMom;e1++) Value+=NominalMC->GetBinContent(e1+1,e0+1);
    double ErrorSquared=0;
    
    for(int e1=0;e1<NBinsMom;e1++){//loop over effect 1
      for(int f1=0;f1<NBinsMom;f1++){//loop over effect 1
	ErrorSquared+=CovarianceMatrix[e1][e0][f1][e0];//Case of angle histogram
      }
    }
    double Error=TMath::Sqrt(ErrorSquared);
    double RelativeError=Error / ( Value == 0 ? 1 : Value ); 
    NominalMC_RecAngle->SetBinContent(e0+1,RelativeError);
    //NominalMC_RecAngle->SetBinError(e0+1,Error);
  }
  
}

int main(int argc, char ** argv){

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  //TApplication theApp("App",0,0);
  gStyle->SetOptFit(kTRUE);
  char * OutputName = new char[256];
  bool UnfoldedPlots=false;

  bool PM=true;  
  bool MC=false;
  bool Reconstructed=false;
  bool SideBand=false;
  
  int c=-1;
  while ((c = getopt(argc, argv, "o:wpmrb")) != -1) {
    switch(c){
    case 'o':
      OutputName=optarg;
      break;      
    case 'p':
      PM=true;
      break;      
    case 'w':
      PM=false;
      break;      
    case 'm':
      MC=true;
      break;      
    case 'r':
      Reconstructed=true;
      break;
    case 'b':
      SideBand=true;
      cout<<"We will use a side band"<<endl;
      break;
    }
  }

  InitializeGlobal();
#ifdef DEBUG
  cout<<"Initialized"<<endl;
#endif
  //To modify the bining if we have a side band
  ReinitializeUnfoldingBinning(SideBand);
#ifdef DEBUG
  cout<<"ReInitialized"<<endl;
#endif
  Xsec * XS = new Xsec();
  
  int NBinsMom;
  int NBinsAngle;
  double * BinningMom; double * BinningAngle;
  
  if(Reconstructed){
    NBinsMom=NBinsRecMom;
    NBinsAngle=NBinsRecAngle;
    BinningMom = BinningRecMom;
    BinningAngle = BinningRecAngle;
  }
  else{
    NBinsMom=NBinsTrueMom;
    NBinsAngle=NBinsTrueAngle;
    BinningMom = BinningTrueMom;
    BinningAngle = BinningTrueAngle;
  }  


  // Step 1. Declare and initialize the variables and tables
  char suffix[3];sprintf(suffix,(PM?"":"_WM"));
  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  cout<<"Selected detector is "<<DetName<<endl;

  cout<<"CAREFUL: TO UNDERSTAND -> When XS error is varied, XS of 0 is different from the nominal MC that I have!!!!"<<endl; 
  char * txtDataName = new char[512];
  char * txtMCName = new char[512];
  char * txtDataNameSB = new char[512];
  char * txtMCNameSB = new char[512];

    /////////////////////////////////////////////
  double **** MCReconstructedEvents_TrueSignal = new double ***[NBinsTrueMom];
  for(int h=0;h<NBinsTrueMom;h++){
    MCReconstructedEvents_TrueSignal[h] = new double **[NBinsTrueAngle];
    for(int i=0;i<NBinsTrueAngle;i++){
      MCReconstructedEvents_TrueSignal[h][i] = new double *[NBinsRecMom];
      for(int j=0;j<NBinsRecMom;j++){
	MCReconstructedEvents_TrueSignal[h][i][j] = new double [NBinsRecAngle];
      }      
    }
  }
  
  //////////////////////////////////////////////////////////////
  double **** CovarianceFlux = new double ***[NBinsMom];
  double **** CorrelationFlux = new double ***[NBinsMom];
  double **** CovarianceStatistics = new double ***[NBinsMom];
  double **** CorrelationStatistics = new double ***[NBinsMom];
  double **** CovarianceXS = new double ***[NBinsMom];
  double **** CorrelationXS = new double ***[NBinsMom];
 
  for(int h=0;h<NBinsMom;h++){
    CovarianceFlux[h] = new double **[NBinsAngle];
    CorrelationFlux[h] = new double **[NBinsAngle];
    CovarianceXS[h] = new double **[NBinsAngle];
    CorrelationXS[h] = new double **[NBinsAngle];
    CovarianceStatistics[h] = new double **[NBinsAngle];
    CorrelationStatistics[h] = new double **[NBinsAngle];
    for(int i=0;i<NBinsAngle;i++){
      CovarianceFlux[h][i] = new double *[NBinsMom];
      CorrelationFlux[h][i] = new double *[NBinsMom];
      CovarianceXS[h][i] = new double *[NBinsMom];
      CorrelationXS[h][i] = new double *[NBinsMom];
      CovarianceStatistics[h][i] = new double *[NBinsMom];
      CorrelationStatistics[h][i] = new double *[NBinsMom];
      for(int j=0;j<NBinsMom;j++){
	CovarianceFlux[h][i][j] = new double [NBinsAngle];
	CorrelationFlux[h][i][j] = new double [NBinsAngle];
	CovarianceXS[h][i][j] = new double [NBinsAngle];
	CorrelationXS[h][i][j] = new double [NBinsAngle];
	CovarianceStatistics[h][i][j] = new double [NBinsAngle];
	CorrelationStatistics[h][i][j] = new double [NBinsAngle];
	for(int k=0;k<NBinsAngle;k++){
	  CovarianceFlux[h][i][j][k] = 0;
	  CorrelationFlux[h][i][j][k] = 0;
	  CovarianceXS[h][i][j][k] = 0;
	  CorrelationXS[h][i][j][k] = 0;
	  CovarianceStatistics[h][i][j][k] = 0;
	  CorrelationStatistics[h][i][j][k] = 0;
	}
      }      
    }
  }

  double ** Efficiency = new double*[NBinsTrueMom];
  for(int i=0;i<NBinsTrueMom;i++){
    Efficiency[i] = new double [NBinsTrueAngle];
    for(int j=0;j<NBinsTrueAngle;j++){
      Efficiency[i][j] = 0;
    }
  }

  double ** DataReconstructedEvents = new double*[NBinsRecMom];
  double ** MCReconstructedBkgEvents = new double*[NBinsRecMom];
  double ** MCReconstructedEvents = new double*[NBinsRecMom];
  double ** MCTrueEvents = new double*[NBinsRecMom];
  
  for(int i=0;i<NBinsRecMom;i++){
    DataReconstructedEvents[i] = new double [NBinsRecAngle];
    MCReconstructedBkgEvents[i] = new double [NBinsRecAngle];
    MCReconstructedEvents[i] = new double [NBinsRecAngle];
    MCTrueEvents[i] = new double [NBinsRecAngle];
    for(int j=0;j<NBinsRecAngle;j++){
      DataReconstructedEvents[i][j] = 0;
      MCReconstructedBkgEvents[i][j] = 0;
      MCReconstructedEvents[i][j] = 0;
      MCTrueEvents[i][j] = 0;
    }
  }

  
  TH2D * NominalMCBkg = new TH2D("NominalMCBkg","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
  double * NumberOfPOT = new double();


  double ***** CovarianceReduced= new double ****[EndError+1];
  double ***** CorrelationReduced= new double ****[EndError+1];
  for(int g=0;g<EndError+1;g++){
    CovarianceReduced[g] = new double ***[NBinsMom];
    CorrelationReduced[g] = new double ***[NBinsMom];
    for(int h=0;h<NBinsMom;h++){
      CovarianceReduced[g][h] = new double **[NBinsAngle];
      CorrelationReduced[g][h] = new double **[NBinsAngle];
      for(int i=0;i<NBinsAngle;i++){
	CovarianceReduced[g][h][i] = new double *[NBinsMom];
	CorrelationReduced[g][h][i] = new double *[NBinsMom];
	for(int j=0;j<NBinsMom;j++){
	  CovarianceReduced[g][h][i][j] = new double [NBinsAngle];
	  CorrelationReduced[g][h][i][j] = new double [NBinsAngle];
	  for(int k=0;k<NBinsAngle;k++){
	    CovarianceReduced[g][h][i][j][k] = 0;
	    CorrelationReduced[g][h][i][j][k] = 0;
	  }
	}
      }
    }
  }
  
  
  TFile * file = new TFile(OutputName,"recreate");
  TDirectory * NoiseEventFunctions = file->mkdir("NoiseEventFunctions");
  TDirectory * HitEfficiencyFunctions = file->mkdir("HitEfficiencytFunctions");
  TDirectory * DataMCErrorFunctions = file->mkdir("DataMCErrorFunctions");
  TDirectory * FluxEventFunctions = file->mkdir("FluxEventFunctions");
  TDirectory * XSEventFunctions = file->mkdir("XSEventFunctions");
  TDirectory * StatisticsEventFunctions = file->mkdir("StatisticsEventFunctions");

  //TH1D * NEvents[EndError+1][NBinsMom][NBinsAngle];
  
  //Selected number of events
  TH2D * NominalSelectedMC = new TH2D("NominalSelectedMC","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
  TH2D * NominalSelectedData = new TH2D("NominalSelectedData","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
  TH1D * NominalSelectedMC_RecMom = new TH1D("NominalSelectedMC_RecMom","",NBinsMom,BinningMom);
  TH1D * NominalSelectedMC_RecAngle = new TH1D("NominalSelectedMC_RecAngle","",NBinsAngle,BinningAngle);
  TH2D * NominalSelectedMC_XSTemp = new TH2D("NominalSelectedMC_XSTemp","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
  TH2D * NominalMCBkg_XSTemp = new TH2D("NominalMCBkg_XSTemp","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
  TH2D * NominalTrueMC = new TH2D("NominalTrueMC","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);

  TH2D * Error_Minus[EndError+1]; TH2D * Error_Plus[EndError+1];
#ifdef TEMP
  TH2D * Error_Minus2[EndError+1]; TH2D * Error_Plus2[EndError+1];
#endif
  TH1D * Error_RecMom[EndError+1];
  TH1D * Error_RecAngle[EndError+1];

  double ErrorTotalStatistics_Plus[NBinsMom][NBinsAngle];
  double ErrorTotalStatistics_Minus[NBinsMom][NBinsAngle];
  double ErrorTotalStatistics_RecMom_Plus[NBinsMom];
  double ErrorTotalStatistics_RecMom_Minus[NBinsMom];
  double ErrorTotalStatistics_RecAngle_Plus[NBinsAngle];
  double ErrorTotalStatistics_RecAngle_Minus[NBinsAngle];

    //ML 2017/07/18 added
  double ErrorTotalDetector_Plus[NBinsMom][NBinsAngle];
  double ErrorTotalDetector_Minus[NBinsMom][NBinsAngle];
  double ErrorTotalDetector_RecMom_Plus[NBinsMom];
  double ErrorTotalDetector_RecMom_Minus[NBinsMom];
  double ErrorTotalDetector_RecAngle_Plus[NBinsAngle];
  double ErrorTotalDetector_RecAngle_Minus[NBinsAngle];

  double ErrorTotalXS_Plus[NBinsMom][NBinsAngle];
  double ErrorTotalXS_Minus[NBinsMom][NBinsAngle];
  double ErrorTotalXS_RecMom_Plus[NBinsMom];
  double ErrorTotalXS_RecMom_Minus[NBinsMom];
  double ErrorTotalXS_RecAngle_Plus[NBinsAngle];
  double ErrorTotalXS_RecAngle_Minus[NBinsAngle];

  double ErrorTotalFlux_Plus[NBinsMom][NBinsAngle];
  double ErrorTotalFlux_Minus[NBinsMom][NBinsAngle];
  double ErrorTotalFlux_RecMom_Plus[NBinsMom];
  double ErrorTotalFlux_RecMom_Minus[NBinsMom];
  double ErrorTotalFlux_RecAngle_Plus[NBinsAngle];
  double ErrorTotalFlux_RecAngle_Minus[NBinsAngle];

  //For flux error
  TH1D * ErrorStatistics[NBinsMom][NBinsAngle];
  TF1 * fErrorStatistics[NBinsMom][NBinsAngle];

  //For Noise error 
  TH1D * ErrorNoise[NBinsMom][NBinsAngle];
  TH1D * ErrorNoise_Norm[NBinsMom][NBinsAngle];
  TF1 * fErrorNoise[NBinsMom][NBinsAngle];

  //For Hit efficiency error
  TH1D * ErrorHitEfficiency[NBinsMom][NBinsAngle];
  TF1 * fErrorHitEfficiency[NBinsMom][NBinsAngle];

  //For flux error
  TH1D * ErrorFlux[NBinsMom][NBinsAngle];
  TF1 * fErrorFlux[NBinsMom][NBinsAngle];

  //For reconstruction error
  TH1D * MCContent[EndError+1][NBinsMom][NBinsAngle];
  TH1D * DataContent[EndError+1][NBinsMom][NBinsAngle];
 
 //For Xsection errors
  TH1D * ErrorXS[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  TH1D * ErrorXS_Norm[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  TSpline3 * sErrorXS[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  TGraph * gErrorXS[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  double xXS[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];
  double yXS[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];

  TAxis * aX_NominalSelectedMC = (TAxis*) NominalSelectedMC->GetXaxis();
  TAxis * aY_NominalSelectedMC = (TAxis*) NominalSelectedMC->GetYaxis();
  TAxis * aX_NominalSelectedData = (TAxis*) NominalSelectedData->GetXaxis();
  TAxis * aY_NominalSelectedData = (TAxis*) NominalSelectedData->GetYaxis();
  TAxis * aX_NominalSelectedMC_RecMom = (TAxis*) NominalSelectedMC_RecMom->GetXaxis();
  TAxis * aY_NominalSelectedMC_RecMom = (TAxis*) NominalSelectedMC_RecMom->GetYaxis();
  TAxis * aX_NominalSelectedMC_RecAngle = (TAxis*) NominalSelectedMC_RecAngle->GetXaxis();
  TAxis * aY_NominalSelectedMC_RecAngle = (TAxis*) NominalSelectedMC_RecAngle->GetYaxis();


  double TotalValue=0;


























  
  // Step 2. Loop over all the errors to determine the variation of number of events for: each error source & bin
  for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
    //Allow to fill the nominal distribution even if there is no detector systematics (and we wish to start error evaluation at 17, where XS starts)
    if(ErrorType>1 && ErrorType<Systematics_Flux_Start) continue;
    
    //if(!(ErrorType==0 || ErrorType==2 || ErrorType==16 || ErrorType>=17 || (ErrorType>=4 && ErrorType<=7))) continue;//TEMP
    //if(!(ErrorType==0 || (ErrorType==3)  /*|| ErrorType==5*/)) continue;//TEMP
    //if(!(ErrorType==0 || (ErrorType>=17))) continue;//TEMP
    //if(! (ErrorType==0 || ErrorType==1 || ErrorType>=16)) continue;//TEMP
    cout<<"The error currently tested is number "<<ErrorType<<endl;

#ifndef DETECTORSYST
    if(ErrorType>=Systematics_Detector_Start && ErrorType<=Systematics_Detector_End) continue;
#endif

       
    //##############################HISTOGRAM INITIALIZATION#########################################
    Error_Minus[ErrorType] = new TH2D(Form("Error_Minus[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    Error_Plus[ErrorType] = new TH2D(Form("Error_Plus[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
#ifdef TEMP
    Error_Minus2[ErrorType] = new TH2D(Form("Error_Minus2[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    Error_Plus2[ErrorType] = new TH2D(Form("Error_Plus2[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
#endif    
    Error_RecMom[ErrorType] = new TH1D(Form("Error_RecMom[%d]",ErrorType),"",NBinsMom,BinningMom);
    Error_RecAngle[ErrorType] = new TH1D(Form("Error_RecAngle[%d]",ErrorType),"",NBinsAngle,BinningAngle);
    
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(ErrorType==1){
	  ErrorStatistics[e0][e1] = new TH1D(Form("ErrorStatistics[%d][%d]",e0,e1),"",600,-3,3);
	  fErrorStatistics[e0][e1] = new TF1(Form("fErrorStatistics[%d][%d]",e0,e1),"gaus",-3,3);
	}
	if(ErrorType==2){
	  ErrorNoise[e0][e1] = new TH1D(Form("ErrorNoise[%d][%d]",e0,e1),"",NE[2]+1,Start[2],End[2]);
	  ErrorNoise_Norm[e0][e1] = new TH1D(Form("ErrorNoise_Norm[%d][%d]",e0,e1),"",NE[2]+1,Start[2],End[2]);
	  fErrorNoise[e0][e1] = new TF1(Form("fErrorNoise[%d][%d]",e0,e1),"pol1",Start[2],End[2]);
	}
	else if(ErrorType==3){
	  // ErrorHitEfficiency[e0][e1] = new TH1D(Form("ErrorHitEfficiency[%d][%d]",e0,e1),"",NE[2]+1,Start[2],End[2]);
	  //	  fErrorHitEfficiency[e0][e1] = new TF1(Form("fErrorHitEfficiency[%d][%d]",e0,e1),"pol1",Start[2],End[2]);
	  ErrorHitEfficiency[e0][e1] = new TH1D(Form("ErrorHitEfficiency[%d][%d]",e0,e1),"",600,-3,3);
	  fErrorHitEfficiency[e0][e1] = new TF1(Form("fErrorHitEfficiency[%d][%d]",e0,e1),"pol1",-3,3);
	}
	else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
	  MCContent[ErrorType][e0][e1] = new TH1D(Form("MCContent[%d][%d][%d]",ErrorType,e0,e1),"",NE[ErrorType]+1,Start[ErrorType],End[ErrorType]);
	  DataContent[ErrorType][e0][e1] = new TH1D(Form("DataContent[%d][%d][%d]",ErrorType,e0,e1),"",NE[ErrorType]+1,Start[ErrorType],End[ErrorType]);
	}
	else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
	  ErrorFlux[e0][e1] = new TH1D(Form("ErrorFlux[%d][%d]",e0,e1),"",600,-3,3);
	  fErrorFlux[e0][e1] = new TF1(Form("fErrorFlux[%d][%d]",e0,e1),"gaus",-3,3);
	}
	else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	  ErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TH1D(Form("ErrorXS[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),"",620,-3.1,3.1);
	  ErrorXS_Norm[ErrorType-Systematics_Xsec_Start][e0][e1] = new TH1D(Form("ErrorXS_Norm[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),"",620,-3.1,3.1);
	  //sErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = ?????new TF1(Form("fErrorFlux[%d][%d]",e0,e1),"gaus",-3,3);
	}
      }
    }
    //##############################END OF INITIALIZATION#########################################


    /*
    //#################################TEMPORARY, SINCE THE XS VARIATION OF 0 SIGMA DOES NOT CORRESPONDS TO THE NOMINAL MC, WE REDEFINE NOMINAL ONLY FOR XS ERROR AS THE 0 SIGMA VARIATION#################################
    if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
      int n=3;

      if(Reconstructed){
	sprintf(txtMCName,"%s/XS/files/MCSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	if(SideBand){
	  sprintf(txtMCNameSB,"%s/XS/files/MCSideBand_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	  XS->Xsec::LoadInputFilesSB(txtMCName,txtMCName,txtMCNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	}
	else XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	//if(MODE == 1){
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    MCReconstructedEvents[e0][e1] -= MCReconstructedBkgEvents[e0][e1];
	  }
	}
	//}
      
      //XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
      }
      else{
	if(SideBand) sprintf(txtMCName,"%s/XS/files/SB/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents);
	//XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	//XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
      }

      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  NominalSelectedMC_XSTemp->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	  NominalMCBkg_XSTemp->SetBinContent(e0+1,e1+1,MCReconstructedBkgEvents[e0][e1]);
	}
      }
    }
      //#################################################################"
      */
    
    
    


    //##############################LOADING: FILL THE HISTOGRAMS#########################################
    cout<<"Error tested="<<ErrorType<<", number of variations of the error="<<NE[ErrorType]<<endl;

      for(int n=0;n<NE[ErrorType];n++){
      double ErrorValue=Start[ErrorType]+n*Step[ErrorType];
      if(Reconstructed){
	sprintf(txtMCName,"%s/XS/files/MCSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	if(SideBand){
	  sprintf(txtMCNameSB,"%s/XS/files/MCSideBand_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	  XS->Xsec::LoadInputFiles_OnlySelectedDataSB(txtMCName,txtMCNameSB,MCReconstructedEvents);
	}
	else XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
	
	if(SideBand){
	  if(ErrorType>=7 && ErrorType<=Systematics_Detector_End) sprintf(txtDataNameSB,"/home/bquilain/CC0pi_XS/XS/files/DataSideBand_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	  else if(ErrorType==0 && EndError>=7 && StartError<=Systematics_Detector_End) sprintf(txtDataNameSB,"/home/bquilain/CC0pi_XS/XS/files/DataSideBand_%s_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	  else sprintf(txtDataNameSB,"%s/XS/files/DataSideBand_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	}
	else{
	  if(ErrorType>=7 && ErrorType<=Systematics_Detector_End) sprintf(txtDataName,"/home/bquilain/CC0pi_XS/XS/files/DataSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	  else if(ErrorType==0 && EndError>=7 && StartError<=Systematics_Detector_End) sprintf(txtDataName,"/home/bquilain/CC0pi_XS/XS/files/DataSelected_%s_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	  else sprintf(txtDataName,"%s/XS/files/DataSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	}
	//else sprintf(txtDataName,"%s/XS/Selection1000_cutBkg%s.txt",cINSTALLREPOSITORY,cINSTALLREPOSITORY,DetName);
	//cout<<"good"<<endl;
	if(MC){
	  //XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,DataReconstructedEvents);
	  //if(ErrorType==0)
	  if(SideBand) XS->Xsec::LoadInputFilesSB(txtMCName,txtMCName,txtMCNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	  else XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	}
	else{
	  //XS->Xsec::LoadInputFiles_OnlySelectedData(txtDataName,DataReconstructedEvents);
	  //if(ErrorType==0)
	  if(SideBand) XS->Xsec::LoadInputFilesSB(txtDataName,txtMCName,txtDataNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	  else XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	}
      }
      else{
	if(SideBand) sprintf(txtMCName,"%s/XS/files/SB/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents);
	
	if(ErrorType>=7 && ErrorType<=Systematics_Detector_End) sprintf(txtDataName,"/home/bquilain/CC0pi_XS/XS/files/DataUnfolded_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else if(ErrorType==0 && EndError>=7 && StartError<=Systematics_Detector_End) sprintf(txtDataName,"/home/bquilain/CC0pi_XS/XS/files/DataUnfolded_%s_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	
	//else sprintf(txtDataName,"%s/XS/Selection1000_cutBkg%s.txt",cINSTALLREPOSITORY,cINSTALLREPOSITORY,DetName);
	//cout<<"good"<<endl;
	if(MC){
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents);
	  //if(ErrorType==0) XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	}
	else{
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtDataName,DataReconstructedEvents,MCTrueEvents);
	  //if(ErrorType==0) XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	}
      }



      
    
      if(Reconstructed){	  	  
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    MCReconstructedEvents[e0][e1] = MCReconstructedEvents[e0][e1] - MCReconstructedBkgEvents[e0][e1];
	    DataReconstructedEvents[e0][e1] = DataReconstructedEvents[e0][e1] - MCReconstructedBkgEvents[e0][e1];
	    //MCReconstructedEvents[e0][e1] = MCReconstructedEvents[e0][e1] - NominalMCBkg->GetBinContent(e0+1,e1+1);
	    //DataReconstructedEvents[e0][e1] = DataReconstructedEvents[e0][e1] - NominalMCBkg->GetBinContent(e0+1,e1+1);
	  }
	}
      }
      /////////////////////////////////	



      
      
      /////////////////////////////////	
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  //cout<<"Error="<<ErrorType<<", "<<e0<<", "<<e1<<", "<<MCReconstructedEvents[e0][e1]<<endl;
	  //For the reconstructed data, MCReconstructedEvents/DataReconstructedEvents contains the number of selected events -> should substract ther bkg
	  //For the unfolded data, MCReconstructedEvents/DataReconstructedEvents contains the number of unfolded(selected - bkg) events -> should NOT substract the bkg
  
	  if(ErrorType==0){
	    cout<<"Value="<<MCReconstructedEvents[e0][e1]<<endl;
	    NominalSelectedMC->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	    NominalTrueMC->SetBinContent(e0+1,e1+1,MCTrueEvents[e0][e1]);
	    NominalSelectedData->SetBinContent(e0+1,e1+1,DataReconstructedEvents[e0][e1]);
	    NominalMCBkg->SetBinContent(e0+1,e1+1,MCReconstructedBkgEvents[e0][e1]);
	    
	    TotalValue+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  }
	  
	  //if( (MODE==1) && ErrorType>1){
	    //TEMP TO CHANGE TO ERRORTYPE>1
	  if( (MODE==1) && ErrorType>Systematics_Xsec_Start && Reconstructed){
	    //MC = SIGNAL Nominal MC - BKG Varied MC. Note that SIGNAL Nominal MC = SELECTED Nominal MC + BKG Nominal MC.  
	    MCReconstructedEvents[e0][e1] = (NominalSelectedMC->GetBinContent(e0+1,e1+1)+NominalMCBkg->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	    DataReconstructedEvents[e0][e1] = (NominalSelectedData->GetBinContent(e0+1,e1+1)+NominalMCBkg->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	    /* 
	       if(MC && ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	      MCReconstructedEvents[e0][e1] = (NominalSelectedMC_XSTemp->GetBinContent(e0+1,e1+1)+NominalMCBkg_XSTemp->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	      //DataReconstructedEvents[e0][e1] = (NominalSelectedData_XSTemp->GetBinContent(e0+1,e1+1)+NominalMCBkg_XSTemp->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	    }
	    */
	  }
	  
	    
	  double RelativeValue = (MCReconstructedEvents[e0][e1] - NominalSelectedMC->GetBinContent(e0+1,e1+1));
	  if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeValue /= NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  //double RelativeValueData = (DataReconstructedEvents[e0][e1] - NominalSelectedData->GetBinContent(e0+1,e1+1));
	  //if(NominalSelectedData->GetBinContent(e0+1,e1+1)!=0) RelativeValue /= NominalSelectedData->GetBinContent(e0+1,e1+1);
	  //cout<< 
	  
	  if(ErrorType==1){
	    //cout<<"Error2="<<ErrorType<<", "<<e0<<", "<<e1<<", "<<MCReconstructedEvents[e0][e1]<<endl;
	    //if(RelativeValue!=0) cout<<e0<<", "<<e1<<", "<<RelativeValue<<endl;
	    ErrorStatistics[e0][e1]->Fill(RelativeValue);
	    //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;

	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		//CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		//cout<<DataReconstructedEvents[e0][e1]-NominalSelectedData->GetBinContent(e0+1,e1+1)<<", "<<DataReconstructedEvents[f0][f1]-NominalSelectedData->GetBinContent(f0,f1)<<endl
		CovarianceStatistics[e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      }
	    }
	  }
	  if(ErrorType==2){
	    ErrorNoise[e0][e1]->Fill(ErrorValue,MCReconstructedEvents[e0][e1]);
	    ErrorNoise_Norm[e0][e1]->Fill(ErrorValue);
	    //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;
	  }
	  else if(ErrorType==3){
	    Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	    Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));

	    ErrorHitEfficiency[e0][e1]->Fill(RelativeValue);//FINALLY USELESS
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(DataReconstructedEvents[e0][e1]-NominalSelectedData->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-Nominal);
		CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(DataReconstructedEvents[e0][e1]-NominalSelectedData->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-Nominal);
	      }
	    }
	  }
	  else if(ErrorType==4){
	    //if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) cout<<"NominalSelectedMC="<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", Varied="<<MCReconstructedEvents[e0][e1]<<endl;
	    Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	    Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		
	      }
	    }	    
	  //NEvents[NBinsMom][NBinsAngle]->Fill(e0,e1,MCReconstructedEvents[e0][e1]);
	  }
	  else if(ErrorType==5){
	    if(MCReconstructedEvents[e0][e1]>=NominalSelectedMC->GetBinContent(e0+1,e1+1)){
	      if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	      Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	    }
	    else Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);

	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      }
	    }	    
	    //cout<<"Birks attenuation length error is changed"
	    //NEvents[NBinsMom][NBinsAngle]->Fill(e0,e1,MCReconstructedEvents[e0][e1]);
	  }
	  else if(ErrorType==6){
	    if(MCReconstructedEvents[e0][e1]>=NominalSelectedMC->GetBinContent(e0+1,e1+1)) Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	    else Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      }
	    }	    
	  }
	  else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
	    double RelativeMC=MCReconstructedEvents[e0][e1];
	    double RelativeData=DataReconstructedEvents[e0][e1];
	    if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    if(NominalSelectedData->GetBinContent(e0+1,e1+1)!=0) RelativeData/=NominalSelectedData->GetBinContent(e0+1,e1+1);
	    if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) cout<<"n="<<n<<", Relative MC="<<RelativeMC<<", Data="<<RelativeData<<endl;
	    //cout<<"Data="<<DataReconstructedEvents[e0][e1]<<", relative="<<RelativeData<<endl;
	    MCContent[ErrorType][e0][e1]->Fill(ErrorValue,RelativeMC);
	    DataContent[ErrorType][e0][e1]->Fill(ErrorValue,RelativeData);
	    //cout<<"Bin is=("<<e0<<","<<e1<<"), MC="<<MCReconstructedEvents[e0][e1]<<", Nominal MC="<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", Data="<<DataReconstructedEvents[e0][e1]<<", Nominal Data="<<NominalSelectedData->GetBinContent(e0+1,e1+1)<<endl;
	    //cout<<"Data="<<DataReconstructedEvents[e0][e1]<<", relative="<<RelativeData<<endl; 
	  }
	else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
	  //cout<<"Bin "<<e0<<", "<<e1<<", Nominal MC="<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", Modified value="<<MCReconstructedEvents[e0][e1]<<", relative="<<RelativeValue<<endl;
	    ErrorFlux[e0][e1]->Fill(RelativeValue);
	    //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		//CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		//cout<<std::scientific<<setprecision(3)<<"("<<e0<<","<<e1<<")     "<<MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1)<<" vs "<<MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1)<<endl;
		CovarianceFlux[e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));

	      }
	    }
	  }
	  else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	    double XsecVariation=ErrorValue-(ErrorType-Systematics_Xsec_Start)*NXsecVariations-CenterXsecVariations;//The variation of Xsec parameter, in #sigma. A number between 0 and 175 - the center of the current systematic source (nominal). For example, for Xsec error source #10, it starts from 7*(10-1)=63 and ends at 70. from 63 to 70, it contains the variariation of -3,-2,-1,0,1,2,3 sigma respectively. The center is then located at 66. For the example of a 2 sigma variation, the substraction will be therefore equal to: 68-66=2, which gives the number of sigmas!
	    //if(XsecVariation==0) NominalSelectedMC_XSTemp->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	    double RelativeMC=MCReconstructedEvents[e0][e1];
	    //if(NominalSelectedMC_XSTemp->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC_XSTemp->GetBinContent(e0+1,e1+1);
	    if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC->GetBinContent(e0+1,e1+1);

	    ErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]->Fill(XsecVariation,RelativeMC);	    
	    ErrorXS_Norm[ErrorType-Systematics_Xsec_Start][e0][e1]->Fill(XsecVariation);
#ifdef DEBUG
	    cout<<XsecVariation<<", variation="<<RelativeMC<<endl;
#endif
	    xXS[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=XsecVariation;
	    yXS[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=RelativeMC;
	    

	    /*	    if(XsecVariation<-1){
	      if(e0==0 && e1==0) cout<<"Temp, keep only variation +-1 sigma in XS"<<endl;
	      yXS[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=1;//TEMP
	      }*/
	    
	  }
	}	  
      }
    }

    
    //cout<<"HERE="<<ErrorType<<endl;
    file->cd();
    if(ErrorType==1){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  StatisticsEventFunctions->cd();
	  ErrorStatistics[e0][e1]->Fit(Form("fErrorStatistics[%d][%d]",e0,e1),"RQ");
	  
	  ErrorStatistics[e0][e1]->Write();
	  	
	  double Error = fErrorStatistics[e0][e1]->GetParameter(2);
	  if(Error != Error) Error=0;//Case where the bin is empty or fit fails and go to NAN  
	  Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-Error);
	  Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,Error);
	}
      }
    }
    if(ErrorType==2){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  NoiseEventFunctions->cd();
	  ErrorNoise[e0][e1]->Divide(ErrorNoise_Norm[e0][e1]);
	  ErrorNoise[e0][e1]->Fit(Form("fErrorNoise[%d][%d]",e0,e1),"RQ");
	  
	  ErrorNoise[e0][e1]->Write();

	  double RelativeValue=fErrorNoise[e0][e1]->Eval(0)-NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeValue/=NominalSelectedMC->GetBinContent(e0+1,e1+1);

	  Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	  Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));

	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      double RelativeValue2=fErrorNoise[f0][f1]->Eval(0)-NominalSelectedMC->GetBinContent(f0+1,f1+1);
	      if(NominalSelectedMC->GetBinContent(f0+1,f1+1)!=0) RelativeValue2/=NominalSelectedMC->GetBinContent(f0+1,f1+1);
	      //Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]*(1+RelativeValue)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]*(1+RelativeValue2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]*(1+RelativeValue)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]*(1+RelativeValue2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	    }
	  }	    

	  
	  //for(int nt=0;nt<NToys;nt++){
	    
	    //for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      //for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      //}
	    //}	    
	  //}
	}
      }
    }
    else if(ErrorType==3){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  HitEfficiencyFunctions->cd();
	  ErrorHitEfficiency[e0][e1]->Fit(Form("fErrorHitEfficiency[%d][%d]",e0,e1),"RQ");
	  
	  ErrorHitEfficiency[e0][e1]->Write();
	  cout<<"CAREFUL: check function before estimation systematic error"<<endl;
	  
	  //Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(fErrorHitEfficiency[e0][e1]
	}
      }
    }
    else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  
	  double ePlus=0;
	  double eMinus=0;
	  for(int ibin=1;ibin<=MCContent[ErrorType][e0][e1]->GetNbinsX();ibin++){
	    double Error=(DataContent[ErrorType][e0][e1]->GetBinContent(ibin)-MCContent[ErrorType][e0][e1]->GetBinContent(ibin));
	    if( Error > ePlus ) ePlus=Error;
	    else if( (-Error) > eMinus ) eMinus=Error;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		double Error2=(DataContent[ErrorType][f0][f1]->GetBinContent(ibin)-MCContent[ErrorType][f0][f1]->GetBinContent(ibin));
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(NominalSelectedMC->GetBinContent(e0+1,e1+1)*(1+Error)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NominalSelectedMC->GetBinContent(f0+1,f1+1)*(1+Error2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(NominalSelectedMC->GetBinContent(e0+1,e1+1)*(1+Error)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NominalSelectedMC->GetBinContent(f0+1,f1+1)*(1+Error2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      }
	    }
	  }
	Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,ePlus);
	Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,eMinus);
	DataMCErrorFunctions->cd();
	MCContent[ErrorType][e0][e1]->Write();
	DataContent[ErrorType][e0][e1]->Write();
	//Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(fErrorHitEfficiency[e0][e1]

	}
      }
    }
    else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  FluxEventFunctions->cd();
	  ErrorFlux[e0][e1]->Fit(Form("fErrorFlux[%d][%d]",e0,e1),"RQ");
	  ErrorFlux[e0][e1]->Write();

	  double Error = fErrorFlux[e0][e1]->GetParameter(2);
	  if(Error != Error) Error=0;//Case where the bin is empty or fit fails and go to NAN
	  Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-Error);
	  Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,Error);
	}
      }
    }
    else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){      
      //create the spline
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  XSEventFunctions->cd();
	  gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations,xXS[ErrorType-Systematics_Xsec_Start][e0][e1],yXS[ErrorType-Systematics_Xsec_Start][e0][e1]);
	  gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("gErrorXS[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1));
	  sErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TSpline3(Form("sErrorXS[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]);
	  sErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("sErrorXS_%d_%d_%d",ErrorType-Systematics_Xsec_Start,e0,e1)); 
	}
      }
    }
    
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    //cout<<Covariance[ErrorType][ErrorType][f0][f1][f0][f1]<<endl;
	    //double Diagonal=Covariance[ErrorType][ErrorType][e0][e1][e0][e1]*Covariance[ErrorType][ErrorType][f0][f1][f0][f1];
	    //Correlation[ErrorType][ErrorType][e0][e1][f0][f1]=Covariance[ErrorType][ErrorType][e0][e1][f0][f1];
	    //if(Diagonal!=0) Correlation[ErrorType][ErrorType][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	    double Diagonal=CovarianceReduced[ErrorType][e0][e1][e0][e1]*CovarianceReduced[ErrorType][f0][f1][f0][f1];
	    
	    CorrelationReduced[ErrorType][e0][e1][f0][f1]=CovarianceReduced[ErrorType][e0][e1][f0][f1];
	    if(Diagonal!=0) CorrelationReduced[ErrorType][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	  }
	}
      }
    }
    
    if(ErrorType<Systematics_Xsec_Start){
      file->cd();
      Error_Plus[ErrorType]->Write();
      Error_Minus[ErrorType]->Write(); 
    }
  }
  

 
  //#########################SPECIAL TREATMENT FOR XSECTION SOURCE: THE CORRELATION IS ALSO CHECKED BETWEEN THE ERROR SOUCRES (NOT ONLY THE BINS)
  double TotalRelxsError[EndXsec+1];
  
    if(EndError>=Systematics_Xsec_Start){
    //cout<<"hello"<<endl;
  //Starts the toy experiments:

  int NBinsTotal=NBinsMom*NBinsAngle;
  double Error;
  double NEventsXS[NBinsTotal];
  int NToysXsec=500;
  TRandom3 * rxs = new TRandom3();
  
  for(int s1=0;s1<=EndError-Systematics_Xsec_Start;s1++){
  //for(int s1=0;s1<=0;s1++){//TEMP
    
      cout<<"source tested="<<s1+Systematics_Xsec_Start<<endl;

      
      for(int nt=0;nt<NToysXsec;nt++){
	if(nt%100==0) cout<<"toy #"<<nt<<endl;
	Error=10;
	
	while(TMath::Abs(Error)>3) Error=rxs->Gaus(0,1);//My interpolation doesn't go further away than -+3sigma. Therefore, only keep toy experiment inside these values.
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    int bin1=NBinsAngle*e0+e1;
	    NEventsXS[bin1]=NominalSelectedMC->GetBinContent(e0+1,e1+1)*sErrorXS[s1][e0][e1]->Eval(Error);
	    
	    
	    
	    
	    //if(e0==4 && e1==4) cout<<"Error="<<Error<<", "<<NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", value="<<sErrorXS[2][4][4]->Eval(Error)<<endl;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		int bin2=NBinsAngle*f0+f1;
		NEventsXS[bin2]=NominalSelectedMC->GetBinContent(f0+1,f1+1)*sErrorXS[s1][f0][f1]->Eval(Error);
		//CovarianceReduced[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS[bin2]-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		
		CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS[bin2]-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		CovarianceXS[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS[bin2]-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		//CovarianceCurrent[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS[bin2]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		//Covariance[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS[bin2]-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		
		//CovarianceReduced[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[bin1]-sErrorXS[s1][e0][e1]->Eval(0))*(NEventsXS[bin2]-sErrorXS[s1][f0][f1]->Eval(0));
#ifdef DEBUG
		if(CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1] != CovarianceXS[e0][e1][f0][f1]) cout<<"Source="<<s1<<", Duel between CovRed="<<CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]<<" and CovXS="<<CovarianceXS[e0][e1][f0][f1]<<endl;
#endif
	      }
	    }	
	  }
	}
      }
  
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	//double err=TMath::Sqrt(CovarianceReduced[s1][e0][e1][e0][e1]);
	double err=TMath::Sqrt(CovarianceXS[e0][e1][e0][e1]);
	if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) err/=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	Error_Minus[Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,-err);
	Error_Plus[Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,err);

#ifdef TEMP
	double err2=TMath::Sqrt(CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][e0][e1]);
	if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) err2/=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	Error_Minus2[Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,-err2);
	Error_Plus2[Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,err2);
#endif
	
	//cout<<"Error="<<1e6*TMath::Sqrt(CovarianceReduced[e0][e1][e0][e1])<<endl;
	
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      TotalRelxsError[s1]+=CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1];
	      double Diagonal=CovarianceXS[e0][e1][e0][e1]*CovarianceXS[f0][f1][f0][f1];
	      if(Diagonal!=0) CorrelationXS[e0][e1][f0][f1]/=pow(Diagonal,1/2);
	      //double Diagonal=Covariance[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][e0][e1][e0][e1]*Covariance[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][f0][f1][f0][f1];
	      //if(Diagonal!=0) Correlation[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	      //double Diagonal=CovarianceReduced[s1][e0][e1][e0][e1]*CovarianceReduced[s1][f0][f1][f0][f1];
	      //if(Diagonal!=0) CorrelationReduced[s1][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	      
	    }
	  }
      }
    }
    if(TotalValue!=0) TotalRelxsError[s1]=TMath::Sqrt(TotalRelxsError[s1])/TotalValue;
    
    file->cd();
    Error_Minus[Systematics_Xsec_Start+s1]->Write();
    Error_Plus[Systematics_Xsec_Start+s1]->Write();
#ifdef TEMP
    Error_Minus2[Systematics_Xsec_Start+s1]->Write();
    Error_Plus2[Systematics_Xsec_Start+s1]->Write();
#endif
    Evaluate1DError(NominalSelectedMC, CovarianceReduced[Systematics_Xsec_Start+s1], Error_RecMom[Systematics_Xsec_Start+s1], Error_RecAngle[Systematics_Xsec_Start+s1]);
    Error_RecMom[Systematics_Xsec_Start+s1]->Write();
    Error_RecAngle[Systematics_Xsec_Start+s1]->Write();
    
#ifdef XSTABLE
    
    cout<<"Total error so far, when adding="<<Systematics_Xsec_Start+s1<<endl;
    for(int e0=0;e0<NBinsMom;e0++){
      for(int e1=0;e1<NBinsAngle;e1++){
	double ErrorSquared=CovarianceXS[e0][e1][e0][e1];
	double Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	//cout<<"Error XS="<<setprecision(3)<<std::scientific<<ErrorSquared<<endl;
	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;
	cout<<setprecision(0)<<std::fixed<<RelativeError*100<<"%, ";
      }
      cout<<endl;
    }
    
#endif
  }
    }

    //###################################################################################
    //###################################################################################
    //###################################################################################
    //###################################################################################
    //###################################################################################
    // From here, we do not have access to XS errors individually. They are all summed in
    // a XS covariance matrix. This is for RAM reason.
    // It also means that any individual treatment of the systematics should be done
    // before this text!
    //###################################################################################
    //###################################################################################
    //###################################################################################
    //###################################################################################

    /*
#ifdef DEBUG
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    cout << CovarianceFlux[e0][e1][f0][f1] << " ";
#endif
*/











    
    
#ifdef XSTABLE
    /*
    for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      if(ErrorType>=1 && ErrorType<Systematics_Xsec_Start) continue;
      cout<<"Error="<<ErrorType<<endl;
      for(int e0=0;e0<NBinsMom;e0++){
	for(int e1=0;e1<NBinsAngle;e1++){
	  cout<<setprecision(0)<<std::fixed<<"+"<<Error_Plus[ErrorType]->GetBinContent(e0+1,e1+1)*100<<"%/"<<Error_Minus[ErrorType]->GetBinContent(e0+1,e1+1)*100<<"%, ";
	}
	cout<<endl;
      }
      cout<<endl<<endl;
      }*/
#endif
  //#########################################DRAWING PART################################
  TBox * boxsliceErrorStat_RecMom[NBinsMom][NBinsAngle];
  TBox * boxsliceErrorFlux_RecMom[NBinsMom][NBinsAngle];
  TBox * boxsliceErrorXS_RecMom[NBinsMom][NBinsAngle];

  TBox * boxsliceErrorStat_RecAngle[NBinsMom][NBinsAngle];
  TBox * boxsliceErrorFlux_RecAngle[NBinsMom][NBinsAngle];
  TBox * boxsliceErrorXS_RecAngle[NBinsMom][NBinsAngle];

  TBox * boxErrorStat_RecMom[NBinsMom];
  TBox * boxErrorFlux_RecMom[NBinsMom];
  TBox * boxErrorXS_RecMom[NBinsMom];

  TBox * boxErrorStat_RecAngle[NBinsAngle];
  TBox * boxErrorFlux_RecAngle[NBinsAngle];
  TBox * boxErrorXS_RecAngle[NBinsAngle];

  
  //###################################CREATE 1D SLICE DISTRIBUTION##################"
  TH1D * sliceNominalSelectedMC_RecMom[NBinsAngle];
  TH1D * sliceNominalSelectedMC_RecAngle[NBinsMom];
  for(int e1=0;e1<NBinsAngle;e1++) sliceNominalSelectedMC_RecMom[e1] = (TH1D*) NominalSelectedMC->ProjectionX(Form("sliceNominalSelectedMC_RecMom[%d]",e1),e1+1,e1+1);
  for(int e0=0;e0<NBinsMom;e0++) sliceNominalSelectedMC_RecAngle[e0] = (TH1D*) NominalSelectedMC->ProjectionX(Form("sliceNominalSelectedMC_RecAngle[%d]",e0),e0+1,e0+1);
  //###################################


  
  //#########################################STAT ERROR################################
  cout<<"Statistical error:"<<endl;

  double TotalErrorStatSquared=0,TotalRelErrorStat=0.;
  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      for(int f0=0;f0<NBinsMom;f0++){
	for(int f1=0;f1<NBinsAngle;f1++){
	  TotalErrorStatSquared+=CovarianceStatistics[e0][e1][f0][f1];
	}
      }
      
      double Value=0;double Error=0;double ErrorSquared=0;
      Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);

      if(Reconstructed) ErrorSquared=Value;
      else ErrorSquared=CovarianceStatistics[e0][e1][e0][e1];
      Error=TMath::Sqrt(ErrorSquared);

      double RelativeError=TMath::Sqrt(ErrorSquared);
      if(Value!=0) RelativeError/=Value;      
      ErrorTotalStatistics_Plus[e0][e1]=RelativeError;
      ErrorTotalStatistics_Minus[e0][e1]=-RelativeError;
      
      //      NominalSelectedMC->SetBinContent(e0+1,e1+1,Value);
      sliceNominalSelectedMC_RecMom[e1]->SetBinContent(e0+1,Value);
      sliceNominalSelectedMC_RecMom[e1]->SetBinError(e0+1,ErrorSquared);
      
      sliceNominalSelectedMC_RecAngle[e0]->SetBinContent(e1+1,Value);
      sliceNominalSelectedMC_RecAngle[e0]->SetBinError(e1+1,ErrorSquared);
      
      boxsliceErrorStat_RecMom[e0][e1] = new TBox(NominalSelectedMC->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
      boxsliceErrorStat_RecMom[e0][e1]->SetFillColor(kRed);

      boxsliceErrorStat_RecAngle[e0][e1] = new TBox(NominalSelectedMC->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxsliceErrorStat_RecAngle[e0][e1]->SetFillColor(kRed);
    }
  }


  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    
    double Value=0;double Error=0;double ErrorSquared=0;
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      Value+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
    }

    if(Reconstructed) ErrorSquared=Value;
    else{
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	
	for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	  //ErrorSquared+=Covariance[ErrorType][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceStatistics[e0][e1][e0][f1];
	  //cout<<ErrorSquared<<endl;
	  if(CovarianceFlux[e0][e1][e0][f1]!=0) cout<<CovarianceFlux[e0][e1][e0][f1]<<endl;
	}
      }
    }
    
    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalStatistics_RecMom_Plus[e0]=RelativeError;
    ErrorTotalStatistics_RecMom_Minus[e0]=-RelativeError;

    Error=TMath::Sqrt(ErrorSquared);

    cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<Error<<endl;
    NominalSelectedMC_RecMom->SetBinContent(e0+1,Value);
    NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);

    boxErrorStat_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
    boxErrorStat_RecMom[e0]->SetFillColor(kRed);
  }

  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0

    double Value=0;double Error=0;double ErrorSquared=0;
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
      Value+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
    }
    if(Reconstructed) ErrorSquared=Value;
    else{
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	  //ErrorSquared+=Covariance[ErrorType][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceStatistics[e0][e1][f0][e1];
	  //cout<<ErrorSquared<<endl;
	  //if(CovarianceFlux[e0][e1][e0][f1]!=0) cout<<CovarianceFlux[e0][e1][e0][f1]<<endl;
	}
      }
    }

    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalStatistics_RecAngle_Plus[e1]=RelativeError;
    ErrorTotalStatistics_RecAngle_Minus[e1]=-RelativeError;
    
    Error=TMath::Sqrt(ErrorSquared);

    cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
    NominalSelectedMC_RecAngle->SetBinContent(e1+1,Value);
    NominalSelectedMC_RecAngle->SetBinError(e1+1,ErrorSquared);

    boxErrorStat_RecAngle[e1] = new TBox(NominalSelectedMC_RecAngle->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
    boxErrorStat_RecAngle[e1]->SetFillColor(kRed);
    
    if(TotalValue!=0) TotalRelErrorStat=TMath::Sqrt(TotalErrorStatSquared)/TotalValue;
  }
  //#########################################END OF STAT ERROR################################

  //################# ML newly added 2017/07/18  DETECTOR ERROR  #############################
  double TotalErrorDetSquared=0,TotalRelErrorDet=0.;
#ifdef DETECTORSYST
  cout<<"Detector error:"<<endl;
    
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      double Value=NominalSelectedMC->GetBinError(e0+1,e1+1);

      double total_error_plus=0,total_error_minus=0;
      double total_error=0;
      for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	if(ErrorType!=5) continue;
	// I assume uncorrelated error sources 
	total_error_plus+=pow(Error_Plus[ErrorType]->GetBinContent(e0+1,e1+1),2);
	total_error_minus+=pow(Error_Minus[ErrorType]->GetBinContent(e0+1,e1+1),2);
	total_error+=CovarianceReduced[ErrorType][e0][e1][e0][e1];
      }
      ErrorTotalDetector_Plus[e0][e1]=TMath::Sqrt(total_error_plus);// relative error
      ErrorTotalDetector_Minus[e0][e1]=TMath::Sqrt(total_error_minus); // relative error	
      //	cout<<e0<<" "<<e1<<" "<<Value<<" "<<total_error_plus<<" "<<total_error_minus<<" "<<total_error/Value/Value<<endl;

      double ErrorSquared=NominalSelectedMC->GetBinError(e0+1,e1+1);
	
      /*double ErrorPlus=TMath::Sqrt(ErrorSquared+total_error_plus*Value*Value);// absolute error
	double ErrorMinus=TMath::Sqrt(ErrorSquared+total_error_minus*Value*Value);//absolute error
	double Error=(ErrorPlus+ErrorMinus)/2;// I can only store symmetric errors in the TH1
	ErrorSquared=pow(Error,2.);
      */
      ErrorSquared+=total_error;
      double Error=TMath::Sqrt(ErrorSquared);

      sliceNominalSelectedMC_RecMom[e1]->SetBinError(e0+1,ErrorSquared);
      sliceNominalSelectedMC_RecAngle[e0]->SetBinError(e1+1,ErrorSquared);
      
      boxsliceErrorDet_RecMom[e0][e1] = new TBox(NominalSelectedMC->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
      boxsliceErrorDet_RecMom[e0][e1]->SetFillColor(kOrange);

      boxsliceErrorDet_RecAngle[e0][e1] = new TBox(NominalSelectedMC->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxsliceErrorDet_RecAngle[e0][e1]->SetFillColor(kOrange);
    }
  }


  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    
    double Value=0;double total_error_plus=0,total_error_minus=0,total_error=0;
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
      // total_error_plus / _minus are wrong (need the bin-to-bin correlations)
      total_error_plus+=pow(Value*ErrorTotalDetector_Plus[e0][e1],2); // abs squared error
      total_error_minus+=pow(Value*ErrorTotalDetector_Minus[e0][e1],2); // abs squared error
      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	  if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	  if(ErrorType!=5) continue;
	  total_error+=CovarianceReduced[ErrorType][e0][e1][e0][f1];
	}
      }
    }
    double ErrorSquared=NominalSelectedMC_RecMom->GetBinError(e0+1);
    /*      double ErrorPlus=TMath::Sqrt(ErrorSquared+total_error_plus);// absolute error 
	    double ErrorMinus=TMath::Sqrt(ErrorSquared+total_error_minus);//absolute error
	    double Error=(ErrorPlus+ErrorMinus)/2;// I can only store symmetric errors in the TH1
	    ErrorSquared=pow(Error,2.);
    */
    ErrorSquared+=total_error;
    double Error=TMath::Sqrt(ErrorSquared);

    Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);
    //cout<<e0<<" "<<Value<<" "<<total_error_plus/Value/Value<<" "<<total_error_minus/Value/Value<<" "<<total_error/Value/Value<<endl;

    ErrorTotalDetector_RecMom_Plus[e0]=(Value>0?TMath::Sqrt(total_error)/Value:0);// relative error
    ErrorTotalDetector_RecMom_Minus[e0]=(Value>0?TMath::Sqrt(total_error)/Value:0); // relative error	
      
    cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", +-"<<TMath::Sqrt(total_error)<<endl;
    NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);

    boxErrorDet_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
    boxErrorDet_RecMom[e0]->SetFillColor(kOrange);
  }

  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
    
    double Value=0;double total_error_plus=0,total_error_minus=0,total_error=0;
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
      total_error_plus+=pow(Value*ErrorTotalDetector_Plus[e0][e1],2); // abs squared error
      total_error_minus+=pow(Value*ErrorTotalDetector_Minus[e0][e1],2); // abs squared error
      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	  if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	  if(ErrorType!=5) continue;
	  total_error+=CovarianceReduced[ErrorType][e0][e1][f0][e1];
	}
      }
    }
    double ErrorSquared=NominalSelectedMC_RecAngle->GetBinError(e1+1);

    /*      double ErrorPlus=TMath::Sqrt(ErrorSquared+total_error_plus);// absolute error 
	    double ErrorMinus=TMath::Sqrt(ErrorSquared+total_error_minus);//absolute error
	    double Error=(ErrorPlus+ErrorMinus)/2;// I can only store symmetric errors in the TH1
	    ErrorSquared=pow(Error,2.);
    */
    ErrorSquared+=total_error;
    double Error=TMath::Sqrt(ErrorSquared);
      
    Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);
    ErrorTotalDetector_RecAngle_Plus[e1]=(Value>0?TMath::Sqrt(total_error)/Value:0);// relative error
    ErrorTotalDetector_RecAngle_Minus[e1]=(Value>0?TMath::Sqrt(total_error)/Value:0); // relative error	

    cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", +"<<total_error_plus<<", -"<<total_error_minus<<endl;
    NominalSelectedMC_RecAngle->SetBinError(e1+1,ErrorSquared);

    boxErrorDet_RecAngle[e1] = new TBox(NominalSelectedMC_RecAngle->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
    boxErrorDet_RecAngle[e1]->SetFillColor(kOrange);
  }

  for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
    if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
    if(ErrorType!=5) continue;
    for(int e0=0;e0<NBinsMom;e0++){
      for(int e1=0;e1<NBinsAngle;e1++){
	for(int f0=0;f0<NBinsMom;f0++){
	  for(int f1=0;f1<NBinsAngle;f1++){
	    TotalErrorDetSquared+=CovarianceReduced[ErrorType][e0][e1][f0][f1];
	  }
	}
      }
    }
  }
  if(TotalValue!=0) TotalRelErrorDet=TMath::Sqrt(TotalErrorDetSquared)/TotalValue;
  
  //#########################################END OF DETECTOR ERROR################################
#endif

  

  //#########################################XS ERROR################################
  double TotalErrorXSSquared=0,TotalRelErrorXS=0.;
  
  if(EndError>=Systematics_Xsec_Start){
    cout<<"Xsec error:"<<endl;

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){
	  for(int f1=0;f1<NBinsAngle;f1++){
	    TotalErrorXSSquared+=CovarianceXS[e0][e1][f0][f1];
	  }
	}
	
	double Value=0;double Error=0;double ErrorSquared=0;
	Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	ErrorSquared=CovarianceXS[e0][e1][e0][e1];
	/*
	ErrorSquared=0;
	for(int s1=Systematics_Xsec_Start;s1<=EndError;s1++){
#ifdef DEBUG
	  cout<<"source="<<s1<<", err="<<CovarianceReduced[s1][e0][e1][e0][e1]<<endl;
#endif
	  ErrorSquared+=CovarianceReduced[s1][e0][e1][e0][e1];
	  }*/
	double RelativeError=TMath::Sqrt(ErrorSquared);
	//cout<<"Value="<<Value<<", Err="<<
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalXS_Plus[e0][e1]=RelativeError;
	ErrorTotalXS_Minus[e0][e1]=-RelativeError;

	ErrorSquared=ErrorSquared+NominalSelectedMC->GetBinContent(e0+1,e1+1);
	Error=TMath::Sqrt(ErrorSquared);
	
	NominalSelectedMC->SetBinError(e0+1,e1+1,ErrorSquared);
	sliceNominalSelectedMC_RecMom[e1]->SetBinError(e0+1,ErrorSquared);
	sliceNominalSelectedMC_RecAngle[e0]->SetBinError(e1+1,ErrorSquared);
	
	boxsliceErrorXS_RecMom[e0][e1] = new TBox(NominalSelectedMC->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
	boxsliceErrorXS_RecMom[e0][e1]->SetFillColor(kGreen);
	
	boxsliceErrorXS_RecAngle[e0][e1] = new TBox(NominalSelectedMC->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
	boxsliceErrorXS_RecAngle[e0][e1]->SetFillColor(kGreen);
      }
    }




  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);    
    ErrorSquared=0;    
   /*
    //TEMP
    double ErrTemp=0;
    for(int s1=Systematics_Xsec_Start;s1<=EndError;s1++){
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	  //ErrorSquared+=CovarianceXS[e0][e1][e0][f1];      
	  //cout<<"source="<<s1<<", err="<<CovarianceReduced[s1][e0][e1][e0][e1]<<endl;
	  ErrorSquared+=CovarianceReduced[s1][e0][e1][e0][f1];
	  ErrTemp+=CovarianceReduced[s1][e0][e1][e0][f1];
	}
      }
      cout<<"Mom="<<BinningRecMom[e0]<<", s1="<<s1<<", value="<<Value<<", "<<100*TMath::Sqrt(ErrTemp) / ( Value == 0 ? 1 : Value )<<endl;
    }
    */
    //BEFORE
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	ErrorSquared+=CovarianceXS[e0][e1][e0][f1];      
      }
    }
    
    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalXS_RecMom_Plus[e0]=RelativeError;
    ErrorTotalXS_RecMom_Minus[e0]=-RelativeError;
    
    ErrorSquared=ErrorSquared+NominalSelectedMC_RecMom->GetBinError(e0+1);
    Error=TMath::Sqrt(ErrorSquared);
    cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<RelativeError*100<<endl;
    //NominalSelectedMC_RecMom->SetBinContent(e0+1,Value);
    NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);

    boxErrorXS_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),NominalSelectedMC_RecMom->GetBinContent(e0+1)-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),NominalSelectedMC_RecMom->GetBinContent(e0+1)+Error);
    boxErrorXS_RecMom[e0]->SetFillColor(kGreen);
  }






  
  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);    
    
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	//ErrorSquared+=CovarianceXS[e0][e1][f0][e1];      
	for(int s1=Systematics_Xsec_Start;s1<=EndError;s1++){
	  //cout<<"source="<<s1<<", err="<<CovarianceReduced[s1][e0][e1][e0][e1]<<endl;
	  ErrorSquared+=CovarianceReduced[s1][e0][e1][f0][e1];
	}
      }
    }
    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalXS_RecAngle_Plus[e1]=RelativeError;
    ErrorTotalXS_RecAngle_Minus[e1]=-RelativeError;
    ErrorSquared=ErrorSquared+NominalSelectedMC_RecAngle->GetBinError(e1+1);
    Error=TMath::Sqrt(ErrorSquared);
    cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
    //NominalSelectedMC_RecAngle->SetBinContent(e1+1,Value);
    NominalSelectedMC_RecAngle->SetBinError(e1+1,ErrorSquared);

    boxErrorXS_RecAngle[e1] = new TBox(NominalSelectedMC_RecAngle->GetXaxis()->GetBinLowEdge(e1+1),NominalSelectedMC_RecAngle->GetBinContent(e1+1)-Error,NominalSelectedMC_RecAngle->GetXaxis()->GetBinUpEdge(e1+1),NominalSelectedMC_RecAngle->GetBinContent(e1+1)+Error);
    boxErrorXS_RecAngle[e1]->SetFillColor(kGreen);
  }
  if(TotalValue!=0) TotalRelErrorXS=TMath::Sqrt(TotalErrorXSSquared)/TotalValue;
  }
  //#########################################END OF XS ERROR################################

  //#########################################FLUX ERROR################################
  double TotalErrorFluxSquared=0,TotalRelErrorFlux=0.;
  
  if(EndError>=Systematics_Flux_Start){
    cout<<"Flux error:"<<endl;
    //int ErrorType=ErrorType>=Systematics_Flux_Start;

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){
	  for(int f1=0;f1<NBinsAngle;f1++){
	    TotalErrorFluxSquared+=CovarianceFlux[e0][e1][f0][f1];
	  }
	}
	
	double Value=0;double Error=0;double ErrorSquared=0;
	Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	ErrorSquared=CovarianceFlux[e0][e1][e0][e1];
	
	//cout<<setprecision(3)<<std::scientific<<e0<<", "<<e1<<", "<<ErrorSquared<<endl;
	
	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalFlux_Plus[e0][e1]=RelativeError;
	ErrorTotalFlux_Minus[e0][e1]=-RelativeError;

	ErrorSquared+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	Error=TMath::Sqrt(ErrorSquared);
	
	NominalSelectedMC->SetBinError(e0+1,e1+1,ErrorSquared);
	sliceNominalSelectedMC_RecMom[e1]->SetBinError(e0+1,ErrorSquared);
	sliceNominalSelectedMC_RecAngle[e0]->SetBinError(e1+1,ErrorSquared);
	
	boxsliceErrorFlux_RecMom[e0][e1] = new TBox(NominalSelectedMC->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
	boxsliceErrorFlux_RecMom[e0][e1]->SetFillColor(kBlue);
	
	boxsliceErrorFlux_RecAngle[e0][e1] = new TBox(NominalSelectedMC->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
	boxsliceErrorFlux_RecAngle[e0][e1]->SetFillColor(kBlue);
	
      }
    }

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      double Value=0;double ErrorSquared=0;double Error=0;
      Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);

      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1

	for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	  //ErrorSquared+=Covariance[ErrorType][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceFlux[e0][e1][e0][f1];
	  //cout<<ErrorSquared<<endl;
	  //if(CovarianceFlux[e0][e1][e0][f1]!=0) cout<<CovarianceFlux[e0][e1][e0][f1]<<endl;
	}
      }
    double RelativeError=TMath::Sqrt(ErrorSquared);
    double temp=ErrorSquared;
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalFlux_RecMom_Plus[e0]=RelativeError;
    ErrorTotalFlux_RecMom_Minus[e0]=-RelativeError;
    ErrorSquared+=NominalSelectedMC_RecMom->GetBinError(e0+1);

    Error=TMath::Sqrt(ErrorSquared);

    cout<<"Flux error, Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<Error<<", relative="<<RelativeError<<"%"<<", squared="<<temp<<endl;
    boxErrorFlux_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
    boxErrorFlux_RecMom[e0]->SetFillColor(kBlue);
    NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);
    }
  
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
      double Value=0;double ErrorSquared=0;double Error=0;
      Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);
      //ErrorSquared=NominalSelectedMC_RecAngle->GetBinError(e1+1);
      
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	  //ErrorSquared+=Covariance[ErrorType][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceFlux[e0][e1][f0][e1];
	  //cout<<ErrorSquared<<endl;
	  //if(CovarianceFlux[e0][e1][e0][f1]!=0) cout<<CovarianceFlux[e0][e1][e0][f1]<<endl;
	}
      }
      double RelativeError=TMath::Sqrt(ErrorSquared);
      if(Value!=0) RelativeError/=Value;      
      ErrorTotalFlux_RecAngle_Plus[e1]=RelativeError;
      ErrorTotalFlux_RecAngle_Minus[e1]=-RelativeError;
      ErrorSquared+=NominalSelectedMC_RecAngle->GetBinError(e1+1);
      
      Error=TMath::Sqrt(ErrorSquared);

      cout<<"Flux error, Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
      boxErrorFlux_RecAngle[e1] = new TBox(NominalSelectedMC_RecAngle->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxErrorFlux_RecAngle[e1]->SetFillColor(kBlue);
      NominalSelectedMC_RecAngle->SetBinError(e1+1,ErrorSquared);
    }
    
    if(TotalValue!=0) TotalRelErrorFlux=TMath::Sqrt(TotalErrorFluxSquared)/TotalValue;
  }

  //#########################################END FLUX ERROR################################



  
  //#####################################PUT BACK TOTAL ERROR TO SQRT (COVARIANCE)################################
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      double Value=0;double ErrorSquared=0;double Error=0;
      ErrorSquared=NominalSelectedMC->GetBinError(e0+1,e1+1);
      Error=TMath::Sqrt(ErrorSquared);
      //Error=ErrorSquared;
      NominalSelectedMC->SetBinError(e0+1,e1+1,Error);
      sliceNominalSelectedMC_RecMom[e1]->SetBinError(e0+1,Error);
      sliceNominalSelectedMC_RecAngle[e0]->SetBinError(e1+1,Error);
    }
  }

  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    ErrorSquared=NominalSelectedMC_RecMom->GetBinError(e0+1);
    Error=TMath::Sqrt(ErrorSquared);
    //Error=ErrorSquared;
    NominalSelectedMC_RecMom->SetBinError(e0+1,Error);
  }
  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    ErrorSquared=NominalSelectedMC_RecAngle->GetBinError(e1+1);
    Error=TMath::Sqrt(ErrorSquared);
    //Error=ErrorSquared;
    NominalSelectedMC_RecAngle->SetBinError(e1+1,Error);
  }
  //#####################################END OF PUT BACK TOTAL ERROR TO SQRT (COVARIANCE)################################


  
  //#####################################DRAWING OF 1D MOM & ANGLE PLOTS################################
  
  TCanvas * csliceRecMom[NBinsAngle];

  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
    csliceRecMom[e1]= new TCanvas(Form("csliceRecMom[%d]",e1));
    sliceNominalSelectedMC_RecMom[e1]->Draw("E1");
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(EndError>=Systematics_Flux_Start) boxsliceErrorFlux_RecMom[e0][e1]->Draw("same");
      if(EndError>=Systematics_Xsec_Start) boxsliceErrorXS_RecMom[e0][e1]->Draw("same");
      boxsliceErrorStat_RecMom[e0][e1]->Draw("same");
    }
    sliceNominalSelectedMC_RecMom[e1]->Draw("E1same");

    if(Reconstructed) sliceNominalSelectedMC_RecMom[e1]->GetXaxis()->SetTitle("d_{#mu} (cm)");
    else sliceNominalSelectedMC_RecMom[e1]->GetXaxis()->SetTitle("p_{#mu} (GeV)");
   
    sliceNominalSelectedMC_RecMom[e1]->GetYaxis()->SetTitle("Number of events");
    sliceNominalSelectedMC_RecMom[e1]->SetLineColor(1);
    sliceNominalSelectedMC_RecMom[e1]->SetLineWidth(2);
    sliceNominalSelectedMC_RecMom[e1]->GetYaxis()->SetTitleOffset(1.3);
    csliceRecMom[e1]->Write(Form("Canvas_Nominal_RecMom_slice%d",e1));
  }

  TCanvas * cRecMom = new TCanvas("cRecMom");
  NominalSelectedMC_RecMom->Draw("E1");
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    if(EndError>=Systematics_Flux_Start) boxErrorFlux_RecMom[e0]->Draw("same");
    if(EndError>=Systematics_Xsec_Start) boxErrorXS_RecMom[e0]->Draw("same");
    boxErrorStat_RecMom[e0]->Draw("same");
  }
  NominalSelectedMC_RecMom->Draw("E1same");

  if(Reconstructed) NominalSelectedMC_RecMom->GetXaxis()->SetTitle("d_{#mu} (cm)");
  else NominalSelectedMC_RecMom->GetXaxis()->SetTitle("p_{#mu} (GeV)");
  
  NominalSelectedMC_RecMom->GetYaxis()->SetTitle("Number of events");
  NominalSelectedMC_RecMom->SetLineColor(1);
  NominalSelectedMC_RecMom->SetLineWidth(2);
  NominalSelectedMC_RecMom->GetYaxis()->SetTitleOffset(1.3);
  cRecMom->Write("Canvas_Nominal_RecMom");
  
  TCanvas * cRecAngle = new TCanvas("cRecAngle");
  NominalSelectedMC_RecAngle->Draw("E1");
  for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
    if(EndError>=Systematics_Flux_Start) boxErrorFlux_RecAngle[e0]->Draw("same");
      if(EndError>=Systematics_Xsec_Start) boxErrorXS_RecAngle[e0]->Draw("same");
    boxErrorStat_RecAngle[e0]->Draw("same");
  }
  NominalSelectedMC_RecAngle->Draw("E1same");
  NominalSelectedMC_RecAngle->SetLineColor(1);
  NominalSelectedMC_RecAngle->SetLineWidth(2);
  NominalSelectedMC_RecAngle->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  NominalSelectedMC_RecAngle->GetYaxis()->SetTitle("Number of events");
  NominalSelectedMC_RecAngle->GetYaxis()->SetTitleOffset(1.3);
  
  cRecAngle->Write("Canvas_Nominal_RecAngle");
  //#####################################END OF DRAWING################################

  //#########################################ERROR TABLE####################################
  cout<<"2D error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      cout<<setprecision(2)<<std::fixed<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_Plus[e0][e1])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_Plus[e0][e1])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_Plus[e0][e1])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_Plus[e0][e1])<<"(stat),   ";
    }
    cout<<endl;
  }
  cout<<endl<<"1D Momentum error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecMom->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecMom_Plus[e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecMom_Plus[e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_RecMom_Plus[e0])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_RecMom_Plus[e0])<<"(stat)"<<endl;
  }
  
  cout<<endl<<"1D Angle error table"<<endl;  
  for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecAngle->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecAngle_Plus[e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecAngle_Plus[e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_RecAngle_Plus[e0])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_RecAngle_Plus[e0])<<"(stat)"<<endl;
  }

  cout<<endl<<"1bin result"<<endl;
  cout<<setprecision(2)<<std::fixed<<TotalValue<<"+-"<<100*TMath::Abs(TotalRelErrorFlux)<<"%(Flux)+-"<<100*TMath::Abs(TotalRelErrorXS)<<"%(XS)+-"<<100*TMath::Abs(TotalRelErrorDet)<<"%(Det)+-"<<100*TMath::Abs(TotalRelErrorStat)<<"%(stat)"<<endl;

  if(EndError>=Systematics_Xsec_Start){
    cout<<"xscn systematics breakdown:"<<endl; 
    for(int i=StartXsec;i<=EndXsec;i++)
      cout<<setprecision(2)<<"param "<<i<<" +-"<<100*TMath::Abs(TotalRelxsError[i])<<"%"<<endl;
  }
  //#########################################END ERROR TABLE####################################


  //#########################################ERROR TABLE####################################
  cout<<"2D error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      cout<<setprecision(2)<<std::fixed<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_Plus[e0][e1])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_Plus[e0][e1])<<"(XS)+-"<<TMath::Abs(ErrorTotalStatistics_Plus[e0][e1])<<"(stat),   ";
    }
    cout<<endl;
  }
  cout<<endl<<"1D Momentum error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecMom->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecMom_Plus[e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecMom_Plus[e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalStatistics_RecMom_Plus[e0])<<"(stat)"<<endl;
  }
  
  cout<<endl<<"1D Angle error table"<<endl;  
  for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecAngle->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecAngle_Plus[e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecAngle_Plus[e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalStatistics_RecAngle_Plus[e0])<<"(stat)"<<endl;
  }

  //#########################################END ERROR TABLE####################################

    
  NominalSelectedMC_RecMom->Write();
  NominalSelectedMC_RecAngle->Write();
  NominalSelectedData->Write();
  NominalSelectedMC->Write();
  NominalTrueMC->Write();
  NominalMCBkg->Write();
  TH1D * NominalMCBkg_RecMom = (TH1D*) NominalMCBkg->ProjectionX("NominalMCBkg_RecMom",0,NominalMCBkg->GetNbinsX());
  NominalMCBkg_RecMom->Write();
  //NominalMCBkg->Sumw2();
  //NominalSelectedMC->Sumw2();
  NominalMCBkg->Divide(NominalSelectedMC);
  NominalMCBkg->Write("ContaminationMC");
  file->Close();
  return 0;
}

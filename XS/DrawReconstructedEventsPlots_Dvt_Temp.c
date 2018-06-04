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
#define TEMP
//#define DEBUG
#define XSTABLE
#define ANTIBUG
//3 modes should not be mixed:
//1. Predictions MC: (Varied MC - varied bkg) - (Fixed MC - fixed bkg)
//2. Error to scale to data: (Fixed MC - varied bkg) - (Fixed MC - fixed bkg)
//3. Unfolded error: the way I constructed the code provide us:
//U.(Fixed MC - varied bkg) - U.(Fixed MC - fixed bkg)
// So, if we wish to evaluate the variation of MC prediction: this is the mode 0. For the real systematics, select the mode 1
const int MODE=2;

string NameXSParam[EndXsec+1] = {"p$_{F}$ $^{12}$C",
				 "p$_{F}$ $^{16}$O",
				 "E$_{B}$ $^{12}$C",
				 "E$_{B}$ $^{16}$O",
				 "MEC norm. $^{12}$C",
				 "MEC norm. $^{16}$O",
				 "$M_{A}^{QE}$",
				 "$C_{A}$",
				 "$M_{A}^{RES}$",
				 "$I_{bkg}$",
				 "$CC$-other shape",
				 "$CC$-coherent norm.",
				 "$NC$-coherent norm.",
				 "$NC$-other shape",
				 "FSI $\\pi$ absorption",
				 "FSI charge exchange, LE",
				 "FSI charge exchange, HE",
				 "FSI $\\pi$ QE scattering, LE",
				 "FSI $\\pi$ QE scattering, HE",
				 "FSI $\\pi$ production"};


void Evaluate1DError(TH2D * NominalMC, double **** CovarianceMatrix, TH1D * NominalMC_RecMom, TH1D * NominalMC_RecAngle){

  int NBinsMom = NominalMC->GetNbinsX();
  int NBinsAngle = NominalMC->GetNbinsY();

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
  bool BadTuning = true;
  int NDetectors=1;
  
  int c=-1;
  while ((c = getopt(argc, argv, "o:wpmra")) != -1) {
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
    case 'a'://Stands for "all": Compute both detectors.
      NDetectors=3;
      cout<<"We will analyze both PM and WM"<<endl;
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
  //TO DO BNJ
  char suffix[3];sprintf(suffix,(PM?"":"_WM"));
  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  cout<<"Selected detector is "<<DetName<<endl;

  cout<<"CAREFUL: TO UNDERSTAND -> When XS error is varied, XS of 0 is different from the nominal MC that I have!!!!"<<endl; 
  char * txtDataName = new char[512];
  char * txtMCName = new char[512];
  char * txtDataNameSB = new char[512];
  char * txtMCNameSB = new char[512];

    /////////////////////////////////////////////
  //Array for temporary reading the files
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

  
  //////////////////////////////////////////////////////////////
  // Covariance and errors
  double ***** CovarianceFlux = new double ****[2];
  double ***** CorrelationFlux = new double ****[2];
  double ***** CovarianceStatistics = new double ****[2];
  double ***** CorrelationStatistics = new double ****[2];
  double ***** CovarianceXS = new double ****[2];
  double ***** CorrelationXS = new double ****[2];
 

  for(int d=0;d<NDetectors;d++){
    CovarianceFlux[d] = new double ***[NBinsMom];
    CorrelationFlux[d] = new double ***[NBinsMom];
    CovarianceXS[d] = new double ***[NBinsMom];
    CorrelationXS[d] = new double ***[NBinsMom];
    CovarianceStatistics[d] = new double ***[NBinsMom];
    CorrelationStatistics[d] = new double ***[NBinsMom];
    
    for(int h=0;h<NBinsMom;h++){
      CovarianceFlux[d][h] = new double **[NBinsAngle];
      CorrelationFlux[d][h] = new double **[NBinsAngle];
      CovarianceXS[d][h] = new double **[NBinsAngle];
      CorrelationXS[d][h] = new double **[NBinsAngle];
      CovarianceStatistics[d][h] = new double **[NBinsAngle];
      CorrelationStatistics[d][h] = new double **[NBinsAngle];
      for(int i=0;i<NBinsAngle;i++){
	CovarianceFlux[d][h][i] = new double *[NBinsMom];
	CorrelationFlux[d][h][i] = new double *[NBinsMom];
	CovarianceXS[d][h][i] = new double *[NBinsMom];
	CorrelationXS[d][h][i] = new double *[NBinsMom];
	CovarianceStatistics[d][h][i] = new double *[NBinsMom];
	CorrelationStatistics[d][h][i] = new double *[NBinsMom];
	for(int j=0;j<NBinsMom;j++){
	  CovarianceFlux[d][h][i][j] = new double [NBinsAngle];
	  CorrelationFlux[d][h][i][j] = new double [NBinsAngle];
	  CovarianceXS[d][h][i][j] = new double [NBinsAngle];
	  CorrelationXS[d][h][i][j] = new double [NBinsAngle];
	  CovarianceStatistics[d][h][i][j] = new double [NBinsAngle];
	  CorrelationStatistics[d][h][i][j] = new double [NBinsAngle];
	  for(int k=0;k<NBinsAngle;k++){
	    CovarianceFlux[d][h][i][j][k] = 0;
	    CorrelationFlux[d][h][i][j][k] = 0;
	    CovarianceXS[d][h][i][j][k] = 0;
	    CorrelationXS[d][h][i][j][k] = 0;
	    CovarianceStatistics[d][h][i][j][k] = 0;
	    CorrelationStatistics[d][h][i][j][k] = 0;
	  }
	}    
      }
    }
  }

  double ****** CovarianceReduced= new double *****[2];
  double ****** CorrelationReduced= new double *****[2];
  for(int d=0;d<NDetectors;d++){
    CovarianceReduced[d] = new double ****[EndError+1];
    CorrelationReduced[d] = new double ****[EndError+1];
    for(int g=0;g<EndError+1;g++){
      CovarianceReduced[d][g] = new double ***[NBinsMom];
      CorrelationReduced[d][g] = new double ***[NBinsMom];
      for(int h=0;h<NBinsMom;h++){
	CovarianceReduced[d][g][h] = new double **[NBinsAngle];
	CorrelationReduced[d][g][h] = new double **[NBinsAngle];
	for(int i=0;i<NBinsAngle;i++){
	  CovarianceReduced[d][g][h][i] = new double *[NBinsMom];
	  CorrelationReduced[d][g][h][i] = new double *[NBinsMom];
	  for(int j=0;j<NBinsMom;j++){
	    CovarianceReduced[d][g][h][i][j] = new double [NBinsAngle];
	    CorrelationReduced[d][g][h][i][j] = new double [NBinsAngle];
	    for(int k=0;k<NBinsAngle;k++){
	      CovarianceReduced[d][g][h][i][j][k] = 0;
	      CorrelationReduced[d][g][h][i][j][k] = 0;
	    }
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

  //Selected number of events
  TH2D * NominalSelectedMC[NDetectors];
  TH1D * NominalSelectedMC_RecMom[NDetectors];
  TH1D * NominalSelectedMC_RecAngle[NDetectors];
  TH2D * NominalSelectedData[NDetectors];
  TH2D * NominalSelectedMC_XSTemp[NDetectors];
  TH2D * NominalMCBkg_XSTemp[NDetectors];
  TH2D * NominalTrueMC[NDetectors];
  TH2D * NominalMCBkg[NDetectors];
  double ** NumberOfPOT = new double *[NDetectors];

  TH2D * Error_Minus[NDetectors][EndError+1]; TH2D * Error_Plus[NDetectors][EndError+1];
#ifdef TEMP
  TH2D * Error_Minus2[NDetectors][EndError+1]; TH2D * Error_Plus2[NDetectors][EndError+1];
#endif
  TH1D * Error_RecMom[NDetectors][EndError+1];
  TH1D * Error_RecAngle[NDetectors][EndError+1];

  double ErrorTotalStatistics_Plus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalStatistics_Minus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalStatistics_RecMom_Plus[NDetectors][NBinsMom];
  double ErrorTotalStatistics_RecMom_Minus[NDetectors][NBinsMom];
  double ErrorTotalStatistics_RecAngle_Plus[NDetectors][NBinsAngle];
  double ErrorTotalStatistics_RecAngle_Minus[NDetectors][NBinsAngle];

    //ML 2017/07/18 added
  double ErrorTotalDetector_Plus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalDetector_Minus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalDetector_RecMom_Plus[NDetectors][NBinsMom];
  double ErrorTotalDetector_RecMom_Minus[NDetectors][NBinsMom];
  double ErrorTotalDetector_RecAngle_Plus[NDetectors][NBinsAngle];
  double ErrorTotalDetector_RecAngle_Minus[NDetectors][NBinsAngle];

  double ErrorTotalXS_Plus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalXS_Minus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalXS_RecMom_Plus[NDetectors][NBinsMom];
  double ErrorTotalXS_RecMom_Minus[NDetectors][NBinsMom];
  double ErrorTotalXS_RecAngle_Plus[NDetectors][NBinsAngle];
  double ErrorTotalXS_RecAngle_Minus[NDetectors][NBinsAngle];

  double ErrorTotalFlux_Plus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalFlux_Minus[NDetectors][NBinsMom][NBinsAngle];
  double ErrorTotalFlux_RecMom_Plus[NDetectors][NBinsMom];
  double ErrorTotalFlux_RecMom_Minus[NDetectors][NBinsMom];
  double ErrorTotalFlux_RecAngle_Plus[NDetectors][NBinsAngle];
  double ErrorTotalFlux_RecAngle_Minus[NDetectors][NBinsAngle];

  //For flux error
  TH1D * ErrorStatistics[NDetectors][NBinsMom][NBinsAngle];
  TF1 * fErrorStatistics[NDetectors][NBinsMom][NBinsAngle];

  //For Noise error
  TH1D * ErrorNoise[NDetectors][NBinsMom][NBinsAngle];
  TH1D * ErrorNoise_Norm[NDetectors][NBinsMom][NBinsAngle];
  TF1 * fErrorNoise[NDetectors][NBinsMom][NBinsAngle];

  //For Hit efficiency error
  TH1D * ErrorHitEfficiency[NDetectors][NBinsMom][NBinsAngle];
  TF1 * fErrorHitEfficiency[NDetectors][NBinsMom][NBinsAngle];

  //For flux error
  TH1D * ErrorFlux[NBinsMom][NDetectors][NBinsAngle];
  TF1 * fErrorFlux[NBinsMom][NDetectors][NBinsAngle];

  //For reconstruction error
  TH1D * MCContent[EndError+1][NDetectors][NBinsMom][NBinsAngle];
  TH1D * DataContent[EndError+1][NDetectors][NBinsMom][NBinsAngle];
 
 //For Xsection errors
  TH1D * ErrorXS[NDetectors][Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  TSpline3 * sErrorXS[NDetectors][Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  TGraph * gErrorXS[NDetectors][Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
  double xXS[NDetectors][Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];
  double yXS[NDetectors][Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];

  double * TotalValue = new double[NDetectors];

  int NBinsTotal=NBinsMom*NBinsAngle;
  double NEventsXS[NDetectors][NBinsTotal];
  TRandom3 * rxs = new TRandom3();
  double TotalRelxsError[NDetectors][EndXsec+1];
 

  //For DRAWING, initialized in the loop
  TBox * boxsliceErrorStat_RecMom[NDetectors][NBinsMom][NBinsAngle];
  TBox * boxsliceErrorFlux_RecMom[NDetectors][NBinsMom][NBinsAngle];
  TBox * boxsliceErrorXS_RecMom[NDetectors][NBinsMom][NBinsAngle];

  TBox * boxsliceErrorStat_RecAngle[NDetectors][NBinsMom][NBinsAngle];
  TBox * boxsliceErrorFlux_RecAngle[NDetectors][NBinsMom][NBinsAngle];
  TBox * boxsliceErrorXS_RecAngle[NDetectors][NBinsMom][NBinsAngle];

  TBox * boxErrorStat_RecMom[NDetectors][NBinsMom];
  TBox * boxErrorFlux_RecMom[NDetectors][NBinsMom];
  TBox * boxErrorXS_RecMom[NDetectors][NBinsMom];

  TBox * boxErrorStat_RecAngle[NDetectors][NBinsAngle];
  TBox * boxErrorFlux_RecAngle[NDetectors][NBinsAngle];
  TBox * boxErrorXS_RecAngle[NDetectors][NBinsAngle];
  TH1D * sliceNominalSelectedMC_RecMom[NDetectors][NBinsAngle];
  TH1D * sliceNominalSelectedMC_RecAngle[NDetectors][NBinsMom];

  TCanvas * csliceRecMom[NDetectors][NBinsAngle];
  TCanvas * cRecMom[NDetectors];
  TCanvas * cRecAngle[NDetectors];

  for(int d=0;d<NDetectors;d++){

    NominalSelectedMC[d] = new TH2D(Form("NominalSelectedMC%d",d),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    NominalSelectedMC_RecMom[d] = new TH1D(Form("NominalSelectedMC_RecMom%d",d),"",NBinsMom,BinningMom);
    NominalSelectedMC_RecAngle[d] = new TH1D(Form("NominalSelectedMC_RecAngle%d",d),"",NBinsAngle,BinningAngle);
    NominalSelectedData[d] = new TH2D(Form("NominalSelectedData%d",d),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    NominalSelectedMC_XSTemp[d] = new TH2D(Form("NominalSelectedMC_XSTemp%d",d),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    NominalMCBkg_XSTemp[d] = new TH2D(Form("NominalMCBkg_XSTemp%d",d),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    NominalTrueMC[d] = new TH2D(Form("NominalTrueMC%d",d),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    NominalMCBkg[d] = new TH2D(Form("NominalMCBkg%d",d),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    NumberOfPOT[d] = new double();
    TotalValue[d] = 0;
  }



  
  for(int d=0;d<NDetectors;d++){

    if(NDetectors > 1){//Case where we study all detectors.
      if(d==0) PM=true;
      else if(d==1) PM=false;
      sprintf(suffix,(PM?"":"_WM"));
      sprintf(DetName,(PM?"PM":"WM"));
      cout<<"Selected detector is "<<DetName<<endl;
      cout<<"Key for maintaining the correlations between PM and WM: have the same seed for flux and XS parameter variation. For the flux, it is enough to use the same Covariance flux file in CC0piSelection. For XS though, here we should use the same seed"<<endl;
      rxs->SetSeed(1);
    }

      
    // Step 2. Loop over all the errors to determine the variation of number of events for: each error source & bin
    for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      if(ErrorType == Systematics_Flux_Start+1) continue;
      //if(ErrorType>1 && ErrorType<Systematics_Flux_Start ) continue;
      if(ErrorType>1 && ErrorType<Systematics_Xsec_Start) continue;

      cout<<"The error currently tested is number "<<ErrorType<<endl;
      
#ifndef DETECTORSYST
      if(ErrorType>=Systematics_Detector_Start && ErrorType<=Systematics_Detector_End) continue;
#endif

       
      //##############################HISTOGRAM INITIALIZATION#########################################
      Error_Minus[d][ErrorType] = new TH2D(Form("Error_Minus[%d][%d]",d,ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      Error_Plus[d][ErrorType] = new TH2D(Form("Error_Plus[%d][%d]",d,ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
#ifdef TEMP
      Error_Minus2[d][ErrorType] = new TH2D(Form("Error_Minus2[%d][%d]",d,ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      Error_Plus2[d][ErrorType] = new TH2D(Form("Error_Plus2[%d][%d]",d,ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
#endif    
      Error_RecMom[d][ErrorType] = new TH1D(Form("Error_RecMom[%d][%d]",d,ErrorType),"",NBinsMom,BinningMom);
      Error_RecAngle[d][ErrorType] = new TH1D(Form("Error_RecAngle[%d][%d]",d,ErrorType),"",NBinsAngle,BinningAngle);
      
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(ErrorType==1){
	    ErrorStatistics[d][e0][e1] = new TH1D(Form("ErrorStatistics[%d][%d][%d]",d,e0,e1),"",600,-3,3);
	    fErrorStatistics[d][e0][e1] = new TF1(Form("fErrorStatistics[%d][%d][%d]",d,e0,e1),"gaus",-3,3);
	  }
	  if(ErrorType==2){
	    ErrorNoise[d][e0][e1] = new TH1D(Form("ErrorNoise[%d][%d][%d]",d,e0,e1),"",NE[2]+1,Start[2],End[2]);
	    ErrorNoise_Norm[d][e0][e1] = new TH1D(Form("ErrorNoise_Norm[%d][%d][%d]",d,e0,e1),"",NE[2]+1,Start[2],End[2]);
	    fErrorNoise[d][e0][e1] = new TF1(Form("fErrorNoise[%d][%d][%d]",d,e0,e1),"pol1",Start[2],End[2]);
	  }
	  else if(ErrorType==3){
	    // ErrorHitEfficiency[e0][e1] = new TH1D(Form("ErrorHitEfficiency[%d][%d]",e0,e1),"",NE[2]+1,Start[2],End[2]);
	    //	  fErrorHitEfficiency[e0][e1] = new TF1(Form("fErrorHitEfficiency[%d][%d]",e0,e1),"pol1",Start[2],End[2]);
	    ErrorHitEfficiency[d][e0][e1] = new TH1D(Form("ErrorHitEfficiency[%d][%d][%d]",d,e0,e1),"",600,-3,3);
	    fErrorHitEfficiency[d][e0][e1] = new TF1(Form("fErrorHitEfficiency[%d][%d][%d]",d,e0,e1),"pol1",-3,3);
	  }
	  else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
	    MCContent[d][ErrorType][e0][e1] = new TH1D(Form("MCContent[%d][%d][%d][%d]",d,ErrorType,e0,e1),"",NE[ErrorType]+1,Start[ErrorType],End[ErrorType]);
	    DataContent[d][ErrorType][e0][e1] = new TH1D(Form("DataContent[%d][%d][%d][%d]",d,ErrorType,e0,e1),"",NE[ErrorType]+1,Start[ErrorType],End[ErrorType]);
	  }
	  else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
	    ErrorFlux[d][e0][e1] = new TH1D(Form("ErrorFlux[%d][%d][%d]",d,e0,e1),"",600,-3,3);
	    fErrorFlux[d][e0][e1] = new TF1(Form("fErrorFlux[%d][%d][%d]",d,e0,e1),"gaus",-3,3);
	  }
	  else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	    ErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1] = new TH1D(Form("ErrorXS[%d][%d][%d][%d]",d,ErrorType-Systematics_Xsec_Start,e0,e1),"",620,-3.1,3.1);
	    //sErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1] = ?????new TF1(Form("fErrorFlux[%d][%d]",e0,e1),"gaus",-3,3);
	  }
	}
      }
      //##############################END OF INITIALIZATION#########################################
    
        
      if(d != 2){//d==2 <=> ratio: no need to refill
	
    //#################################TEMPORARY, SINCE THE XS VARIATION OF 0 SIGMA DOES NOT CORRESPONDS TO THE NOMINAL MC, WE REDEFINE NOMINAL ONLY FOR XS ERROR AS THE 0 SIGMA VARIATION#################################
      if(BadTuning){
	if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	  int n=3;
	  
	  if(Reconstructed){
	    sprintf(txtMCName,"%s/XS/files/MCSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	    if(SideBand){
	      sprintf(txtMCNameSB,"%s/XS/files/MCSideBand_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	      XS->Xsec::LoadInputFilesSB(txtMCName,txtMCName,txtMCNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	    }
	    else XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
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
	    //XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	    //XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
	  }
	  
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	      NominalSelectedMC_XSTemp[d]->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	      NominalMCBkg_XSTemp[d]->SetBinContent(e0+1,e1+1,MCReconstructedBkgEvents[e0][e1]);
	    }
	  }
	}
      }
      //#################################################################"
            
    
    
    


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
	  if(SideBand) XS->Xsec::LoadInputFilesSB(txtMCName,txtMCName,txtMCNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	  else XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	}
	else{
	  //XS->Xsec::LoadInputFiles_OnlySelectedData(txtDataName,DataReconstructedEvents);
	  //if(ErrorType==0)
	  if(SideBand) XS->Xsec::LoadInputFilesSB(txtDataName,txtMCName,txtDataNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	  else XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	}
      }
      
      else{
	if(SideBand) sprintf(txtMCName,"%s/XS/files/SB/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	//XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents);
	
	if(ErrorType>=7 && ErrorType<=Systematics_Detector_End) sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else if(ErrorType==0 && EndError>=7 && StartError<=Systematics_Detector_End) sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	
	//else sprintf(txtDataName,"%s/XS/Selection1000_cutBkg%s.txt",cINSTALLREPOSITORY,cINSTALLREPOSITORY,DetName);
	//cout<<"good"<<endl;
	if(MC){
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents);
	  //if(ErrorType==0) XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
	}
	else{
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtDataName,DataReconstructedEvents,MCTrueEvents);
	  cout<<"Yo"<<endl;
	  //if(ErrorType==0) XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT[d]);
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
	    NominalSelectedMC[d]->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	    NominalTrueMC[d]->SetBinContent(e0+1,e1+1,MCTrueEvents[e0][e1]);
	    NominalSelectedData[d]->SetBinContent(e0+1,e1+1,DataReconstructedEvents[e0][e1]);
	    NominalMCBkg[d]->SetBinContent(e0+1,e1+1,MCReconstructedBkgEvents[e0][e1]);
	    TotalValue[d]+=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	  }
	  
	  //if( (MODE==1) && ErrorType>1){
	    //TEMP TO CHANGE TO ERRORTYPE>1
	  if( (MODE==1) && ErrorType>Systematics_Xsec_Start && Reconstructed){
	    //MC = SIGNAL Nominal MC - BKG Varied MC. Note that SIGNAL Nominal MC = SELECTED Nominal MC + BKG Nominal MC.  
	    MCReconstructedEvents[e0][e1] = (NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)+NominalMCBkg[d]->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	    DataReconstructedEvents[e0][e1] = (NominalSelectedData[d]->GetBinContent(e0+1,e1+1)+NominalMCBkg[d]->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	     
	    if(BadTuning && MC && ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	      MCReconstructedEvents[e0][e1] = (NominalSelectedMC_XSTemp[d]->GetBinContent(e0+1,e1+1)+NominalMCBkg_XSTemp[d]->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	      //DataReconstructedEvents[e0][e1] = (NominalSelectedData_XSTemp->GetBinContent(e0+1,e1+1)+NominalMCBkg[d]_XSTemp->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	    }
	    
	  }
	  
	    
	  double RelativeValue = (MCReconstructedEvents[e0][e1] - NominalSelectedMC[d]->GetBinContent(e0+1,e1+1));
	  if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) RelativeValue /= NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	  //double RelativeValueData = (DataReconstructedEvents[e0][e1] - NominalSelectedData[d]->GetBinContent(e0+1,e1+1));
	  //if(NominalSelectedData[d]->GetBinContent(e0+1,e1+1)!=0) RelativeValue /= NominalSelectedData[d]->GetBinContent(e0+1,e1+1);
	  //cout<< 
	  
	  if(ErrorType==1){
	    //cout<<"Error2="<<ErrorType<<", "<<e0<<", "<<e1<<", "<<MCReconstructedEvents[e0][e1]<<endl;
	    //if(RelativeValue!=0) cout<<e0<<", "<<e1<<", "<<RelativeValue<<endl;
	    ErrorStatistics[d][e0][e1]->Fill(RelativeValue);
	    //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;

	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		//CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		//cout<<DataReconstructedEvents[e0][e1]-NominalSelectedData[d]->GetBinContent(e0+1,e1+1)<<", "<<DataReconstructedEvents[f0][f1]-NominalSelectedData[d]->GetBinContent(f0,f1)<<endl
		CovarianceStatistics[d][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
	      }
	    }
	  }
	  if(ErrorType==2){
	    ErrorNoise[d][e0][e1]->Fill(ErrorValue,MCReconstructedEvents[e0][e1]);
	    ErrorNoise_Norm[d][e0][e1]->Fill(ErrorValue);
	    //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;
	  }
	  else if(ErrorType==3){
	    Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	    Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));

	    ErrorHitEfficiency[d][e0][e1]->Fill(RelativeValue);//FINALLY USELESS
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(DataReconstructedEvents[e0][e1]-NominalSelectedData[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-Nominal);
		CovarianceReduced[d][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(DataReconstructedEvents[e0][e1]-NominalSelectedData[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-Nominal);
	      }
	    }
	  }
	  else if(ErrorType==4){
	    //if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) cout<<"NominalSelectedMC[d]="<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<", Varied="<<MCReconstructedEvents[e0][e1]<<endl;
	    Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	    Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		CovarianceReduced[d][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		
	      }
	    }	    
	  //NEvents[NBinsMom][NBinsAngle]->Fill(e0,e1,MCReconstructedEvents[e0][e1]);
	  }
	  else if(ErrorType==5){
	    if(MCReconstructedEvents[e0][e1]>=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)){
	      if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	      Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	    }
	    else Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);

	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		CovarianceReduced[d][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
	      }
	    }	    
	    //cout<<"Birks attenuation length error is changed"
	    //NEvents[NBinsMom][NBinsAngle]->Fill(e0,e1,MCReconstructedEvents[e0][e1]);
	  }
	  else if(ErrorType==6){
	    if(MCReconstructedEvents[e0][e1]>=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)) Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	    else Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		CovarianceReduced[d][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
	      }
	    }	    
	  }
	  else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
	    double RelativeMC=MCReconstructedEvents[e0][e1];
	    double RelativeData=DataReconstructedEvents[e0][e1];
	    if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	    if(NominalSelectedData[d]->GetBinContent(e0+1,e1+1)!=0) RelativeData/=NominalSelectedData[d]->GetBinContent(e0+1,e1+1);
	    //if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) cout<<"n="<<n<<", Relative MC="<<RelativeMC<<", Data="<<RelativeData<<endl;
	    //cout<<"Data="<<DataReconstructedEvents[e0][e1]<<", relative="<<RelativeData<<endl;
	    MCContent[d][ErrorType][e0][e1]->Fill(ErrorValue,RelativeMC);
	    DataContent[d][ErrorType][e0][e1]->Fill(ErrorValue,RelativeData);
	    //cout<<"Bin is=("<<e0<<","<<e1<<"), MC="<<MCReconstructedEvents[e0][e1]<<", Nominal MC="<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<", Data="<<DataReconstructedEvents[e0][e1]<<", Nominal Data="<<NominalSelectedData[d]->GetBinContent(e0+1,e1+1)<<endl;
	    //cout<<"Data="<<DataReconstructedEvents[e0][e1]<<", relative="<<RelativeData<<endl; 
	  }
	  else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
	  //cout<<"Bin "<<e0<<", "<<e1<<", Nominal MC="<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<", Modified value="<<MCReconstructedEvents[e0][e1]<<", relative="<<RelativeValue<<endl;
	    ErrorFlux[d][e0][e1]->Fill(RelativeValue);
	    
	    //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		//CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		//cout<<std::scientific<<setprecision(3)<<"("<<e0<<","<<e1<<")     "<<MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<" vs "<<MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)<<endl;
		CovarianceFlux[d][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));

	      }
	    }
	  }
	  else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	    double XsecVariation=ErrorValue-(ErrorType-Systematics_Xsec_Start)*NXsecVariations-CenterXsecVariations;//The variation of Xsec parameter, in #sigma. A number between 0 and 175 - the center of the current systematic source (nominal). For example, for Xsec error source #10, it starts from 7*(10-1)=63 and ends at 70. from 63 to 70, it contains the variariation of -3,-2,-1,0,1,2,3 sigma respectively. The center is then located at 66. For the example of a 2 sigma variation, the substraction will be therefore equal to: 68-66=2, which gives the number of sigmas!
	    //if(XsecVariation==0) NominalSelectedMC[d]_XSTemp->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	    double RelativeMC=MCReconstructedEvents[e0][e1];
	    //if(NominalSelectedMC[d]_XSTemp->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC[d]_XSTemp->GetBinContent(e0+1,e1+1);
	    if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0){
	      if(BadTuning){
		if(NominalSelectedMC_XSTemp[d]->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC_XSTemp[d]->GetBinContent(e0+1,e1+1);
		else RelativeMC = 1;//Only for bad tuning: case where the nominal MC may be not 0 in one bin, but the neut tuned can be 0 => The relative MC = 0 instead of 1(=no variation). TO REMOVE WHEN BADTUNING IS NOT HERE
	      }
	      else if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);

	      
	    }
#ifdef ANTIBUG
	    if(TMath::Abs(RelativeMC)>1e2) RelativeMC=1;
#endif
	    
	    ErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1]->Fill(XsecVariation,RelativeMC);	    
#ifdef DEBUG
	    cout<<XsecVariation<<", variation="<<RelativeMC<<", Entry="<<MCReconstructedEvents[e0][e1]<<", Nominal="<<NominalSelectedMC_XSTemp[d]->GetBinContent(e0+1,e1+1)<<endl;
#endif
	    xXS[d][ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=XsecVariation;
	    yXS[d][ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=RelativeMC;
	    

	    
	  }
	}	  
      }
      }
      }
      else{//The ratio case:
	if(ErrorType==0){
	  NominalSelectedMC[2] = (TH2D*)NominalSelectedMC[1]->Clone(NominalSelectedMC[2]->GetName());
	  NominalSelectedMC[2]->Divide(NominalSelectedMC[0]);
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	      //NominalSelectedMC[2]->SetBinContent(e0+1,e1+1,NominalSelectedMC[0]->GetBinContent(e0+1,e1+1) != 0 ? NominalSelectedMC[1]->GetBinContent(e0+1,e1+1)/NominalSelectedMC[0]->GetBinContent(e0+1,e1+1) : 1);
	      //cout<<"bin2="<<NominalSelectedMC[2]->GetBinContent(e0+1,e1+1)<<endl;
	    }
	  }
	  NominalSelectedMC_RecMom[2] = (TH1D*) NominalSelectedMC_RecMom[1]->Clone(NominalSelectedMC_RecMom[2]->GetName());
	  NominalSelectedMC_RecMom[2]->Divide(NominalSelectedMC_RecMom[0]);
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    //NominalSelectedMC_RecMom[2]->SetBinContent(e0+1,NominalSelectedMC_RecMom[0]->GetBinContent(e0+1) != 0 ? NominalSelectedMC_RecMom[1]->GetBinContent(e0+1)/NominalSelectedMC_RecMom[0]->GetBinContent(e0+1) : 1);
	    //cout<<"bin="<<NominalSelectedMC_RecMom[2]->GetBinContent(e0+1)<<endl;
	  }
	  NominalSelectedMC_RecAngle[2] = (TH1D*) NominalSelectedMC_RecAngle[1]->Clone(NominalSelectedMC_RecAngle[2]->GetName());
	  NominalSelectedMC_RecAngle[2]->Divide(NominalSelectedMC_RecAngle[0]);
	  NominalSelectedData[2] = (TH2D*) NominalSelectedData[1]->Clone(NominalSelectedData[2]->GetName());
	  NominalSelectedData[2]->Divide(NominalSelectedData[0]);
	  
	  NominalTrueMC[2] = (TH2D*) NominalTrueMC[1]->Clone(NominalTrueMC[2]->GetName());
	  NominalTrueMC[2]->Divide(NominalTrueMC[0]);
	  NominalMCBkg[2] = (TH2D*) NominalMCBkg[1]->Clone(NominalMCBkg[2]->GetName());
	  NominalMCBkg[2]->Divide(NominalMCBkg[0]);
	  NominalTrueMC[2] = (TH2D*) NominalTrueMC[1]->Clone(NominalTrueMC[2]->GetName());
	  NominalTrueMC[2]->Divide(NominalTrueMC[0]);
	  
	  TotalValue[2] = TotalValue[1]; if(TotalValue[1] != 0) TotalValue[2] /= TotalValue[0];
	}
	if(BadTuning){
	  if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	    int n=3;
	    NominalSelectedMC_XSTemp[2] = (TH2D*)NominalSelectedMC_XSTemp[1]->Clone(NominalSelectedMC_XSTemp[2]->GetName());
	    NominalSelectedMC_XSTemp[2]->Divide(NominalSelectedMC_XSTemp[0]);
	    NominalMCBkg_XSTemp[2] = (TH2D*)NominalMCBkg_XSTemp[1]->Clone(NominalMCBkg_XSTemp[2]->GetName());
	    NominalMCBkg_XSTemp[2]->Divide(NominalMCBkg_XSTemp[0]);
	  }
	}

	
	//Clone errors
	Error_Minus[2][ErrorType] = (TH2D*) Error_Minus[1][ErrorType]->Clone(Error_Minus[2][ErrorType]->GetName());
	Error_Minus[2][ErrorType]->Divide(Error_Minus[0][ErrorType]);
	Error_Plus[2][ErrorType] = (TH2D*) Error_Plus[1][ErrorType]->Clone(Error_Plus[2][ErrorType]->GetName());
	Error_Plus[2][ErrorType]->Divide(Error_Plus[0][ErrorType]);
#ifdef TEMP
	Error_Minus2[2][ErrorType] = (TH2D*) Error_Minus2[1][ErrorType]->Clone(Error_Minus2[2][ErrorType]->GetName());
	Error_Minus2[2][ErrorType]->Divide(Error_Minus2[0][ErrorType]);
	Error_Plus2[2][ErrorType] = (TH2D*) Error_Plus2[1][ErrorType]->Clone(Error_Plus2[2][ErrorType]->GetName());
	Error_Plus2[2][ErrorType]->Divide(Error_Plus2[0][ErrorType]);
#endif
	Error_RecMom[2][ErrorType] = (TH1D*) Error_RecMom[1][ErrorType]->Clone(Error_RecMom[2][ErrorType]->GetName());
	Error_RecMom[2][ErrorType]->Divide(Error_RecMom[0][ErrorType]);
	Error_RecAngle[2][ErrorType] = (TH1D*) Error_RecAngle[1][ErrorType]->Clone(Error_RecAngle[2][ErrorType]->GetName());
	Error_RecAngle[2][ErrorType]->Divide(Error_RecAngle[0][ErrorType]);
  
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    //Remark: no need to divise by norm since already done at previous step
	    if(ErrorType==1){
	    }
	    else if(ErrorType==2){
	      ErrorNoise[2][e0][e1] = (TH1D*) ErrorNoise[1][e0][e1]->Clone(ErrorNoise[2][e0][e1]->GetName());
	      ErrorNoise[2][e0][e1]->Divide(ErrorNoise[0][e0][e1]);
	    }
	    else if(ErrorType==3){
	      ErrorHitEfficiency[2][e0][e1] = (TH1D*) ErrorHitEfficiency[1][e0][e1]->Clone(ErrorHitEfficiency[2][e0][e1]->GetName());
	      ErrorHitEfficiency[2][e0][e1]->Divide(ErrorHitEfficiency[0][e0][e1]);
	    }
	    else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
	      MCContent[2][ErrorType][e0][e1] = (TH1D*) MCContent[1][ErrorType][e0][e1]->Clone(MCContent[2][ErrorType][e0][e1]->GetName());
	      MCContent[2][ErrorType][e0][e1]->Divide(MCContent[0][ErrorType][e0][e1]);
	      DataContent[2][ErrorType][e0][e1] = (TH1D*) DataContent[1][ErrorType][e0][e1]->Clone(DataContent[2][ErrorType][e0][e1]->GetName());
	      DataContent[2][ErrorType][e0][e1]->Divide(DataContent[0][ErrorType][e0][e1]);
	    }
	    else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
	      ErrorFlux[2][e0][e1] = (TH1D*) ErrorFlux[1][e0][e1]->Clone(ErrorFlux[2][e0][e1]->GetName());
	      ErrorFlux[2][e0][e1]->Divide(ErrorFlux[0][e0][e1]);
	    }
	    else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){	    
	      ErrorXS[2][ErrorType-Systematics_Xsec_Start][e0][e1] = (TH1D*) ErrorXS[1][ErrorType-Systematics_Xsec_Start][e0][e1]->Clone(ErrorXS[2][ErrorType-Systematics_Xsec_Start][e0][e1]->GetName());
	      ErrorXS[2][ErrorType-Systematics_Xsec_Start][e0][e1]->Divide(ErrorXS[0][ErrorType-Systematics_Xsec_Start][e0][e1]);
	      cout<<"Name: "<<ErrorXS[2][ErrorType-Systematics_Xsec_Start][e0][e1]->GetName()<<", CH: "<<ErrorXS[1][ErrorType-Systematics_Xsec_Start][e0][e1]->GetBinContent(2)<<", H20: "<<ErrorXS[1][ErrorType-Systematics_Xsec_Start][e0][e1]->GetBinContent(2)<<", ratio: "<<ErrorXS[2][ErrorType-Systematics_Xsec_Start][e0][e1]->GetBinContent(2)<<endl;
	      for(int i=0;i<NXsecVariations;i++){
		xXS[2][ErrorType-Systematics_Xsec_Start][e0][e1][i] = xXS[0][ErrorType-Systematics_Xsec_Start][e0][e1][i];
		yXS[2][ErrorType-Systematics_Xsec_Start][e0][e1][i] = yXS[1][ErrorType-Systematics_Xsec_Start][e0][e1][i];
		if(yXS[0][ErrorType-Systematics_Xsec_Start][e0][e1][i] != 0) yXS[2][ErrorType-Systematics_Xsec_Start][e0][e1][i] /= yXS[0][ErrorType-Systematics_Xsec_Start][e0][e1][i];
	      }
	    }
	  }
	}
      }
    
      //###############################################################################################################
      cout<<"Well done buddy. You analyze all the files related to this error and stored into 'Error'='Response functions' histograms. Now, the evaluation starts. It can be either just reading the 'Response function' histograms, or drawing randomly the parameter and use the 'Response functions'"<<endl;
      //CHECKED
        
      
    //cout<<"HERE="<<ErrorType<<endl;
    file->cd();
    if(ErrorType==1){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  StatisticsEventFunctions->cd();
	  ErrorStatistics[d][e0][e1]->Fit(Form("fErrorStatistics[%d][%d][%d]",d,e0,e1),"RQ");
	  ErrorStatistics[d][e0][e1]->Write();
	  double Error = fErrorStatistics[d][e0][e1]->GetParameter(2);
	  if(Error != Error) Error=0;//Case where the bin is empty or fit fails and go to NAN  
	  Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-Error);
	  Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,Error);
	}
      }
    }
    if(ErrorType==2){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  NoiseEventFunctions->cd();
	  if(d!=2) ErrorNoise[d][e0][e1]->Divide(ErrorNoise_Norm[d][e0][e1]);
	  ErrorNoise[d][e0][e1]->Fit(Form("fErrorNoise[%d][%d][%d]",d,e0,e1),"RQ");
	  
	  ErrorNoise[d][e0][e1]->Write();

	  double RelativeValue=fErrorNoise[d][e0][e1]->Eval(0)-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	  if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) RelativeValue/=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);

	  Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	  Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));

	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      double RelativeValue2=fErrorNoise[d][f0][f1]->Eval(0)-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1);
	      if(NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)!=0) RelativeValue2/=NominalSelectedMC[d]->GetBinContent(f0+1,f1+1);
	      //Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]*(1+RelativeValue)-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]*(1+RelativeValue2)-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
	      CovarianceReduced[d][ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]*(1+RelativeValue)-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]*(1+RelativeValue2)-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
	    }
	  }	    

	  
	  //for(int nt=0;nt<NToys;nt++){
	    
	    //for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      //for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
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
	  ErrorHitEfficiency[d][e0][e1]->Fit(Form("fErrorHitEfficiency[%d][%d][%d]",d,e0,e1),"RQ");
	  
	  ErrorHitEfficiency[d][e0][e1]->Write();
	  cout<<"CAREFUL: check function before estimation systematic error"<<endl;
	  
	  //Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(fErrorHitEfficiency[e0][e1]
	}
      }
    }
    else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  
	  double ePlus=0;
	  double eMinus=0;
	  for(int ibin=1;ibin<=MCContent[d][ErrorType][e0][e1]->GetNbinsX();ibin++){
	    double Error=(DataContent[d][ErrorType][e0][e1]->GetBinContent(ibin)-MCContent[d][ErrorType][e0][e1]->GetBinContent(ibin));
	    if( Error > ePlus ) ePlus=Error;
	    else if( (-Error) > eMinus ) eMinus=Error;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		double Error2=(DataContent[d][ErrorType][f0][f1]->GetBinContent(ibin)-MCContent[d][ErrorType][f0][f1]->GetBinContent(ibin));
		//Covariance[d][ErrorType][d][ErrorType][e0][e1][f0][f1]+=(1./(NE[d][ErrorType]-1.))*(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)*(1+Error)-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)*(1+Error2)-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		CovarianceReduced[d][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)*(1+Error)-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)*(1+Error2)-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
	      }
	    }
	  }
	Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,ePlus);
	Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,eMinus);
	DataMCErrorFunctions->cd();
	MCContent[d][ErrorType][e0][e1]->Write();
	DataContent[d][ErrorType][e0][e1]->Write();
	//Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(fErrorHitEfficiency[e0][e1]

	}
      }
    }
    else if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  FluxEventFunctions->cd();
	  ErrorFlux[d][e0][e1]->Fit(Form("fErrorFlux[%d][%d][%d]",d,e0,e1),"RQ");
	  ErrorFlux[d][e0][e1]->Write();

	  double Error = fErrorFlux[d][e0][e1]->GetParameter(2);
	  if(Error != Error) Error=0;//Case where the bin is empty or fit fails and go to NAN
	  Error_Minus[d][ErrorType]->SetBinContent(e0+1,e1+1,-Error);
	  Error_Plus[d][ErrorType]->SetBinContent(e0+1,e1+1,Error);
	}
      }
    }
    else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){      
      //create the spline
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  XSEventFunctions->cd();
	  gErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations,xXS[d][ErrorType-Systematics_Xsec_Start][e0][e1],yXS[d][ErrorType-Systematics_Xsec_Start][e0][e1]);
	  gErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("gErrorXS[%d][%d][%d][%d]",d,ErrorType-Systematics_Xsec_Start,e0,e1));
	  sErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1] = new TSpline3(Form("sErrorXS[%d][%d][%d][%d]",d,ErrorType-Systematics_Xsec_Start,e0,e1),gErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1]);
	  if(e0 == 0 && e1 == 0) cout<<"###################################################"<<endl<<ErrorType-Systematics_Xsec_Start<<endl<<endl;
	  sErrorXS[d][ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("sErrorXS_%d_%d_%d_%d",d,ErrorType-Systematics_Xsec_Start,e0,e1)); 
	}
      }
    }
    
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    //cout<<Covariance[d][ErrorType][d][ErrorType][f0][f1][f0][f1]<<endl;
	    //double Diagonal=Covariance[d][ErrorType][d][ErrorType][e0][e1][e0][e1]*Covariance[d][ErrorType][d][ErrorType][f0][f1][f0][f1];
	    //Correlation[d][ErrorType][d][ErrorType][e0][e1][f0][f1]=Covariance[d][ErrorType][d][ErrorType][e0][e1][f0][f1];
	    //if(Diagonal!=0) Correlation[d][ErrorType][d][ErrorType][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	    double Diagonal=CovarianceReduced[d][ErrorType][e0][e1][e0][e1]*CovarianceReduced[d][ErrorType][f0][f1][f0][f1];
	    
	    CorrelationReduced[d][ErrorType][e0][e1][f0][f1]=CovarianceReduced[d][ErrorType][e0][e1][f0][f1];
	    if(Diagonal!=0) CorrelationReduced[d][ErrorType][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	  }
	}
      }
    }
    
    if(ErrorType<Systematics_Xsec_Start){
      file->cd();
      Error_Plus[d][ErrorType]->Write();
      Error_Minus[d][ErrorType]->Write(); 
    }
  }
	//CHECKED
   
  //#########################SPECIAL TREATMENT FOR XSECTION SOURCE: THE CORRELATION IS ALSO CHECKED BETWEEN THE ERROR SOUCRES (NOT ONLY THE BINS)
  
    if(EndError>=Systematics_Xsec_Start){
      //cout<<"hello"<<endl;
      //Starts the toy experiments:
      
      double Error;
      int NToysXsec=500;
  
      for(int s1=0;s1<=EndError-Systematics_Xsec_Start;s1++){
    
      cout<<"source tested="<<s1+Systematics_Xsec_Start<<endl;

      
      for(int nt=0;nt<NToysXsec;nt++){
	if(nt%100==0) cout<<"toy #"<<nt<<endl;
	Error=10;
	
	while(TMath::Abs(Error)>3) Error=rxs->Gaus(0,1);//My interpolation doesn't go further away than -+3sigma. Therefore, only keep toy experiment inside these values.
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    int bin1=NBinsAngle*e0+e1;
	    //cout<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<endl;
	    //cout<<"Bin XS: "<<e0<<","<<e1<<", content: "<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<", err:"<<sErrorXS[d][s1][e0][e1]->Eval(Error)<<endl;
	    NEventsXS[d][bin1]=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)*sErrorXS[d][s1][e0][e1]->Eval(Error);
	    
	    //if(e0==4 && e1==4) cout<<"Error="<<Error<<", "<<NEventsXS[d][bin1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<", value="<<sErrorXS[2][4][4]->Eval(Error)<<endl;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		int bin2=NBinsAngle*f0+f1;
		NEventsXS[d][bin2]=NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)*sErrorXS[d][s1][f0][f1]->Eval(Error);
		//CovarianceReduced[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[d][bin1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NEventsXS[d][bin2]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)); 
		
		CovarianceReduced[d][Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[d][bin1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NEventsXS[d][bin2]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)); 
		CovarianceXS[d][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[d][bin1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NEventsXS[d][bin2]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)); 
		//CovarianceCurrent[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[d][bin1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NEventsXS[d][bin2]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1));
		//Covariance[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[d][bin1]-NominalSelectedMC[d]->GetBinContent(e0+1,e1+1))*(NEventsXS[d][bin2]-NominalSelectedMC[d]->GetBinContent(f0+1,f1+1)); 
		
		//CovarianceReduced[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS[d][bin1]-sErrorXS[s1][e0][e1]->Eval(0))*(NEventsXS[d][bin2]-sErrorXS[s1][f0][f1]->Eval(0));
#ifdef DEBUG
		//		if(CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1] != CovarianceXS[e0][e1][f0][f1]) cout<<"Source="<<s1<<", Duel between CovRed="<<CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]<<" and CovXS="<<CovarianceXS[e0][e1][f0][f1]<<endl;
#endif
	      }
	    }	
	  }
	}
      }
  
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	//double err=TMath::Sqrt(CovarianceReduced[s1][e0][e1][e0][e1]);
	double err=TMath::Sqrt(CovarianceXS[d][e0][e1][e0][e1]);
	if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) err/=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	Error_Minus[d][Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,-err);
	Error_Plus[d][Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,err);

#ifdef TEMP
	double err2=TMath::Sqrt(CovarianceReduced[d][Systematics_Xsec_Start+s1][e0][e1][e0][e1]);
	if(NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)!=0) err2/=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	Error_Minus2[d][Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,-err2);
	Error_Plus2[d][Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,err2);
#endif
	
	//cout<<"Error="<<1e6*TMath::Sqrt(CovarianceReduced[e0][e1][e0][e1])<<endl;
	
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      TotalRelxsError[d][s1]+=CovarianceReduced[d][Systematics_Xsec_Start+s1][e0][e1][f0][f1];
	      double Diagonal=CovarianceXS[d][e0][e1][e0][e1]*CovarianceXS[d][f0][f1][f0][f1];
	      if(Diagonal!=0) CorrelationXS[d][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	      //double Diagonal=Covariance[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][e0][e1][e0][e1]*Covariance[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][f0][f1][f0][f1];
	      //if(Diagonal!=0) Correlation[Systematics_Xsec_Start+s1][Systematics_Xsec_Start+s1][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	      //double Diagonal=CovarianceReduced[s1][e0][e1][e0][e1]*CovarianceReduced[s1][f0][f1][f0][f1];
	      //if(Diagonal!=0) CorrelationReduced[s1][e0][e1][f0][f1]/=pow(Diagonal,1/2);
	      
	    }
	  }
      }
    }
    if(TotalValue[d]!=0) TotalRelxsError[d][s1]=TMath::Sqrt(TotalRelxsError[d][s1])/TotalValue[d];
      
    file->cd();
    Error_Minus[d][Systematics_Xsec_Start+s1]->Write();
    Error_Plus[d][Systematics_Xsec_Start+s1]->Write();
#ifdef TEMP
    Error_Minus2[d][Systematics_Xsec_Start+s1]->Write();
    Error_Plus2[d][Systematics_Xsec_Start+s1]->Write();
#endif
    Evaluate1DError(NominalSelectedMC[d], CovarianceReduced[d][Systematics_Xsec_Start+s1], Error_RecMom[d][Systematics_Xsec_Start+s1], Error_RecAngle[d][Systematics_Xsec_Start+s1]);
    Error_RecMom[d][Systematics_Xsec_Start+s1]->Write();
    Error_RecAngle[d][Systematics_Xsec_Start+s1]->Write();
    
#ifdef XSTABLE
    
    cout<<"Total error so far, when adding="<<Systematics_Xsec_Start+s1<<endl;
    for(int e0=0;e0<NBinsMom;e0++){
      for(int e1=0;e1<NBinsAngle;e1++){
	double ErrorSquared=CovarianceXS[d][e0][e1][e0][e1];
	double Value=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
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
  
    //CHECKED
    

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


     
  //#########################################DRAWING PART################################

  
  //###################################CREATE 1D SLICE DISTRIBUTION##################"
  for(int e1=0;e1<NBinsAngle;e1++) sliceNominalSelectedMC_RecMom[d][e1] = (TH1D*) NominalSelectedMC[d]->ProjectionX(Form("sliceNominalSelectedMC_RecMom[%d][%d]",d,e1),e1+1,e1+1);
  for(int e0=0;e0<NBinsMom;e0++) sliceNominalSelectedMC_RecAngle[d][e0] = (TH1D*) NominalSelectedMC[d]->ProjectionY(Form("sliceNominalSelectedMC_RecAngle[%d][%d]",d,e0),e0+1,e0+1);
  //###################################


  
  //#########################################STAT ERROR################################
  cout<<"Statistical error:"<<endl;

  double TotalErrorStatSquared=0,TotalRelErrorStat=0.;
  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      for(int f0=0;f0<NBinsMom;f0++){
	for(int f1=0;f1<NBinsAngle;f1++){
	  TotalErrorStatSquared+=CovarianceStatistics[d][e0][e1][f0][f1];
	}
      }
      
      double Value=0;double Error=0;double ErrorSquared=0;
      Value=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);

      if(Reconstructed) ErrorSquared=Value;
      else ErrorSquared=CovarianceStatistics[d][e0][e1][e0][e1];
      Error=TMath::Sqrt(ErrorSquared);

      double RelativeError=TMath::Sqrt(ErrorSquared);
      if(Value!=0) RelativeError/=Value;      
      ErrorTotalStatistics_Plus[d][e0][e1]=RelativeError;
      ErrorTotalStatistics_Minus[d][e0][e1]=-RelativeError;
      
      //      NominalSelectedMC[d]->SetBinContent(e0+1,e1+1,Value);
      sliceNominalSelectedMC_RecMom[d][e1]->SetBinContent(e0+1,Value);
      sliceNominalSelectedMC_RecMom[d][e1]->SetBinError(e0+1,ErrorSquared);
      
      sliceNominalSelectedMC_RecAngle[d][e0]->SetBinContent(e1+1,Value);
      sliceNominalSelectedMC_RecAngle[d][e0]->SetBinError(e1+1,ErrorSquared);
      
      boxsliceErrorStat_RecMom[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
      boxsliceErrorStat_RecMom[d][e0][e1]->SetFillColor(kRed);

      boxsliceErrorStat_RecAngle[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC[d]->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxsliceErrorStat_RecAngle[d][e0][e1]->SetFillColor(kRed);
    }
  }


  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    
    double Value=0;double Error=0;double ErrorSquared=0;
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      Value+=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
    }

    if(Reconstructed) ErrorSquared=Value;
    else{
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	
	for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	  //ErrorSquared+=Covariance[d][ErrorType][d][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceStatistics[d][e0][e1][e0][f1];
	  //cout<<ErrorSquared<<endl;
	  if(CovarianceFlux[d][e0][e1][e0][f1]!=0) cout<<CovarianceFlux[d][e0][e1][e0][f1]<<endl;
	}
      }
    }
    
    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalStatistics_RecMom_Plus[d][e0]=RelativeError;
    ErrorTotalStatistics_RecMom_Minus[d][e0]=-RelativeError;

    Error=TMath::Sqrt(ErrorSquared);

    cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<Error<<endl;
    NominalSelectedMC_RecMom[d]->SetBinContent(e0+1,Value);
    NominalSelectedMC_RecMom[d]->SetBinError(e0+1,ErrorSquared);

    boxErrorStat_RecMom[d][e0] = new TBox(NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
    boxErrorStat_RecMom[d][e0]->SetFillColor(kRed);
  }

  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0

    double Value=0;double Error=0;double ErrorSquared=0;
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
      Value+=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
    }
    if(Reconstructed) ErrorSquared=Value;
    else{
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	  //ErrorSquared+=Covariance[d][ErrorType][d][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceStatistics[d][e0][e1][f0][e1];
	  //cout<<ErrorSquared<<endl;
	  //if(CovarianceFlux[d][e0][e1][e0][f1]!=0) cout<<CovarianceFlux[d][e0][e1][e0][f1]<<endl;
	}
      }
    }

    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalStatistics_RecAngle_Plus[d][e1]=RelativeError;
    ErrorTotalStatistics_RecAngle_Minus[d][e1]=-RelativeError;
    
    Error=TMath::Sqrt(ErrorSquared);

    cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
    NominalSelectedMC_RecAngle[d]->SetBinContent(e1+1,Value);
    NominalSelectedMC_RecAngle[d]->SetBinError(e1+1,ErrorSquared);

    boxErrorStat_RecAngle[d][e1] = new TBox(NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
    boxErrorStat_RecAngle[d][e1]->SetFillColor(kRed);
    
    if(TotalValue[d]!=0) TotalRelErrorStat=TMath::Sqrt(TotalErrorStatSquared)/TotalValue[d];
  }
  //#########################################END OF STAT ERROR################################

  //################# ML newly added 2017/07/18  DETECTOR ERROR  #############################
  double TotalErrorDetSquared=0,TotalRelErrorDet=0.;
#ifdef DETECTORSYST
  cout<<"Detector error:"<<endl;
    
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      double Value=NominalSelectedMC[d]->GetBinError(e0+1,e1+1);

      double total_error_plus=0,total_error_minus=0;
      double total_error=0;
      for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	if(ErrorType!=5) continue;
	// I assume uncorrelated error sources 
	total_error_plus+=pow(Error_Plus[d][ErrorType]->GetBinContent(e0+1,e1+1),2);
	total_error_minus+=pow(Error_Minus[d][ErrorType]->GetBinContent(e0+1,e1+1),2);
	total_error+=CovarianceReduced[d][ErrorType][e0][e1][e0][e1];
      }
      ErrorTotalDetector_Plus[d][e0][e1]=TMath::Sqrt(total_error_plus);// relative error
      ErrorTotalDetector_Minus[d][e0][e1]=TMath::Sqrt(total_error_minus); // relative error	
      //	cout<<e0<<" "<<e1<<" "<<Value<<" "<<total_error_plus<<" "<<total_error_minus<<" "<<total_error/Value/Value<<endl;

      double ErrorSquared=NominalSelectedMC[d]->GetBinError(e0+1,e1+1);
	
      ErrorSquared+=total_error;
      double Error=TMath::Sqrt(ErrorSquared);

      sliceNominalSelectedMC_RecMom[d][e1]->SetBinError(e0+1,ErrorSquared);
      sliceNominalSelectedMC_RecAngle[d][e0]->SetBinError(e1+1,ErrorSquared);
      
      boxsliceErrorDet_RecMom[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
      boxsliceErrorDet_RecMom[d][e0][e1]->SetFillColor(kOrange);

      boxsliceErrorDet_RecAngle[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC[d]->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxsliceErrorDet_RecAngle[d][e0][e1]->SetFillColor(kOrange);
    }
  }


  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    
    double Value=0;double total_error_plus=0,total_error_minus=0,total_error=0;
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      Value=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
      // total_error_plus / _minus are wrong (need the bin-to-bin correlations)
      total_error_plus+=pow(Value*ErrorTotalDetector_Plus[d][e0][e1],2); // abs squared error
      total_error_minus+=pow(Value*ErrorTotalDetector_Minus[d][e0][e1],2); // abs squared error
      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	  if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	  if(ErrorType!=5) continue;
	  total_error+=CovarianceReduced[d][ErrorType][e0][e1][e0][f1];
	}
      }
    }
    double ErrorSquared=NominalSelectedMC_RecMom[d]->GetBinError(e0+1);

    ErrorSquared+=total_error;
    double Error=TMath::Sqrt(ErrorSquared);

    Value=NominalSelectedMC_RecMom[d]->GetBinContent(e0+1);
    //cout<<e0<<" "<<Value<<" "<<total_error_plus/Value/Value<<" "<<total_error_minus/Value/Value<<" "<<total_error/Value/Value<<endl;

    ErrorTotalDetector_RecMom_Plus[d][e0]=(Value>0?TMath::Sqrt(total_error)/Value:0);// relative error
    ErrorTotalDetector_RecMom_Minus[d][e0]=(Value>0?TMath::Sqrt(total_error)/Value:0); // relative error	
      
    cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", +-"<<TMath::Sqrt(total_error)<<endl;
    NominalSelectedMC_RecMom[d]->SetBinError(e0+1,ErrorSquared);

    boxErrorDet_RecMom[d][e0] = new TBox(NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
    boxErrorDet_RecMom[d][e0]->SetFillColor(kOrange);
  }

  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
    
    double Value=0;double total_error_plus=0,total_error_minus=0,total_error=0;
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      Value=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
      total_error_plus+=pow(Value*ErrorTotalDetector_Plus[d][e0][e1],2); // abs squared error
      total_error_minus+=pow(Value*ErrorTotalDetector_Minus[d][e0][e1],2); // abs squared error
      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	  if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	  if(ErrorType!=5) continue;
	  total_error+=CovarianceReduced[d][ErrorType][e0][e1][f0][e1];
	}
      }
    }
    double ErrorSquared=NominalSelectedMC_RecAngle[d]->GetBinError(e1+1);

    ErrorSquared+=total_error;
    double Error=TMath::Sqrt(ErrorSquared);
      
    Value=NominalSelectedMC_RecAngle[d]->GetBinContent(e1+1);
    ErrorTotalDetector_RecAngle_Plus[d][e1]=(Value>0?TMath::Sqrt(total_error)/Value:0);// relative error
    ErrorTotalDetector_RecAngle_Minus[d][e1]=(Value>0?TMath::Sqrt(total_error)/Value:0); // relative error	

    cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", +"<<total_error_plus<<", -"<<total_error_minus<<endl;
    NominalSelectedMC_RecAngle[d]->SetBinError(e1+1,ErrorSquared);

    boxErrorDet_RecAngle[d][e1] = new TBox(NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
    boxErrorDet_RecAngle[d][e1]->SetFillColor(kOrange);
  }

  for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
    if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
    if(ErrorType!=5) continue;
    for(int e0=0;e0<NBinsMom;e0++){
      for(int e1=0;e1<NBinsAngle;e1++){
	for(int f0=0;f0<NBinsMom;f0++){
	  for(int f1=0;f1<NBinsAngle;f1++){
	    TotalErrorDetSquared+=CovarianceReduced[d][ErrorType][e0][e1][f0][f1];
	  }
	}
      }
    }
  }
  if(TotalValue[d]!=0) TotalRelErrorDet=TMath::Sqrt(TotalErrorDetSquared)/TotalValue[d];
  
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
	    TotalErrorXSSquared+=CovarianceXS[d][e0][e1][f0][f1];
	  }
	}
	
	double Value=0;double Error=0;double ErrorSquared=0;
	Value=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	ErrorSquared=CovarianceXS[d][e0][e1][e0][e1];

	double RelativeError=TMath::Sqrt(ErrorSquared);
	//cout<<"Value="<<Value<<", Err="<<
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalXS_Plus[d][e0][e1]=RelativeError;
	ErrorTotalXS_Minus[d][e0][e1]=-RelativeError;

	ErrorSquared=ErrorSquared+NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	Error=TMath::Sqrt(ErrorSquared);
	
	NominalSelectedMC[d]->SetBinError(e0+1,e1+1,ErrorSquared);
	sliceNominalSelectedMC_RecMom[d][e1]->SetBinError(e0+1,ErrorSquared);
	sliceNominalSelectedMC_RecAngle[d][e0]->SetBinError(e1+1,ErrorSquared);
	
	boxsliceErrorXS_RecMom[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
	boxsliceErrorXS_RecMom[d][e0][e1]->SetFillColor(kGreen);
	
	boxsliceErrorXS_RecAngle[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC[d]->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
	boxsliceErrorXS_RecAngle[d][e0][e1]->SetFillColor(kGreen);
      }
    }




  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    Value=NominalSelectedMC_RecMom[d]->GetBinContent(e0+1);    
    ErrorSquared=0;    

    //BEFORE
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	ErrorSquared+=CovarianceXS[d][e0][e1][e0][f1];      
      }
    }
    
    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalXS_RecMom_Plus[d][e0]=RelativeError;
    ErrorTotalXS_RecMom_Minus[d][e0]=-RelativeError;
    
    ErrorSquared=ErrorSquared+NominalSelectedMC_RecMom[d]->GetBinError(e0+1);
    Error=TMath::Sqrt(ErrorSquared);
    cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<RelativeError*100<<endl;
    //NominalSelectedMC_RecMom[d]->SetBinContent(e0+1,Value);
    NominalSelectedMC_RecMom[d]->SetBinError(e0+1,ErrorSquared);

    boxErrorXS_RecMom[d][e0] = new TBox(NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinLowEdge(e0+1),NominalSelectedMC_RecMom[d]->GetBinContent(e0+1)-Error,NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinUpEdge(e0+1),NominalSelectedMC_RecMom[d]->GetBinContent(e0+1)+Error);
    boxErrorXS_RecMom[d][e0]->SetFillColor(kGreen);
  }






  
  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    Value=NominalSelectedMC_RecAngle[d]->GetBinContent(e1+1);    
    
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	//ErrorSquared+=CovarianceXS[e0][e1][f0][e1];      
	for(int s1=Systematics_Xsec_Start;s1<=EndError;s1++){
	  //cout<<"source="<<s1<<", err="<<CovarianceReduced[s1][e0][e1][e0][e1]<<endl;
	  ErrorSquared+=CovarianceReduced[d][s1][e0][e1][f0][e1];
	}
      }
    }
    double RelativeError=TMath::Sqrt(ErrorSquared);
    if(Value!=0) RelativeError/=Value;      
    ErrorTotalXS_RecAngle_Plus[d][e1]=RelativeError;
    ErrorTotalXS_RecAngle_Minus[d][e1]=-RelativeError;
    ErrorSquared=ErrorSquared+NominalSelectedMC_RecAngle[d]->GetBinError(e1+1);
    Error=TMath::Sqrt(ErrorSquared);
    cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
    //NominalSelectedMC_RecAngle[d]->SetBinContent(e1+1,Value);
    NominalSelectedMC_RecAngle[d]->SetBinError(e1+1,ErrorSquared);

    boxErrorXS_RecAngle[d][e1] = new TBox(NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinLowEdge(e1+1),NominalSelectedMC_RecAngle[d]->GetBinContent(e1+1)-Error,NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinUpEdge(e1+1),NominalSelectedMC_RecAngle[d]->GetBinContent(e1+1)+Error);
    boxErrorXS_RecAngle[d][e1]->SetFillColor(kGreen);
  }
  if(TotalValue[d]!=0) TotalRelErrorXS=TMath::Sqrt(TotalErrorXSSquared)/TotalValue[d];
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
	    //cout<<CovarianceFlux[d][e0][e1][f0][f1]<<endl;
	    TotalErrorFluxSquared += CovarianceFlux[d][e0][e1][f0][f1];
	  }
	}

	double Value=0;double Error=0;double ErrorSquared=0;
	Value=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	ErrorSquared=CovarianceFlux[d][e0][e1][e0][e1];
	
	//cout<<setprecision(3)<<std::scientific<<e0<<", "<<e1<<", "<<ErrorSquared<<endl;
	
	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalFlux_Plus[d][e0][e1]=RelativeError;
	ErrorTotalFlux_Minus[d][e0][e1]=-RelativeError;

	ErrorSquared+=NominalSelectedMC[d]->GetBinContent(e0+1,e1+1);
	Error=TMath::Sqrt(ErrorSquared);


	NominalSelectedMC[d]->SetBinError(e0+1,e1+1,ErrorSquared);
	sliceNominalSelectedMC_RecMom[d][e1]->SetBinError(e0+1,ErrorSquared);
	sliceNominalSelectedMC_RecAngle[d][e0]->SetBinError(e1+1,ErrorSquared);


	boxsliceErrorFlux_RecMom[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
	boxsliceErrorFlux_RecMom[d][e0][e1]->SetFillColor(kBlue);
	
	boxsliceErrorFlux_RecAngle[d][e0][e1] = new TBox(NominalSelectedMC[d]->GetYaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC[d]->GetYaxis()->GetBinUpEdge(e1+1),Value+Error);
	boxsliceErrorFlux_RecAngle[d][e0][e1]->SetFillColor(kBlue);
	
      }
    }
  
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      double Value=0;double ErrorSquared=0;double Error=0;
      Value=NominalSelectedMC_RecMom[d]->GetBinContent(e0+1);
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	  //ErrorSquared+=Covariance[d][ErrorType][d][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared += CovarianceFlux[d][e0][e1][e0][f1];
	  //cout<<ErrorSquared<<endl;
	  //if(CovarianceFlux[e0][e1][e0][f1]!=0) cout<<CovarianceFlux[e0][e1][e0][f1]<<endl;
	}
      }
      double RelativeError=TMath::Sqrt(ErrorSquared);
      double temp=ErrorSquared;
      if(Value!=0) RelativeError/=Value;      
      ErrorTotalFlux_RecMom_Plus[d][e0]=RelativeError;
      ErrorTotalFlux_RecMom_Minus[d][e0]=-RelativeError;
      ErrorSquared+=NominalSelectedMC_RecMom[d]->GetBinError(e0+1);
      
      Error=TMath::Sqrt(ErrorSquared);
      
      cout<<"Flux error, Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<Error<<", relative="<<RelativeError<<"%"<<", squared="<<temp<<endl;
      boxErrorFlux_RecMom[d][e0] = new TBox(NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom[d]->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
      boxErrorFlux_RecMom[d][e0]->SetFillColor(kBlue);
      NominalSelectedMC_RecMom[d]->SetBinError(e0+1,ErrorSquared);
    }
  
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
      double Value=0;double ErrorSquared=0;double Error=0;
      Value=NominalSelectedMC_RecAngle[d]->GetBinContent(e1+1);
      //ErrorSquared=NominalSelectedMC_RecAngle[d]->GetBinError(e1+1);
      
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	  //ErrorSquared+=Covariance[d][ErrorType][d][ErrorType][e0][e1][e0][f1];      
	  ErrorSquared+=CovarianceFlux[d][e0][e1][f0][e1];
	  //cout<<ErrorSquared<<endl;
	  //if(CovarianceFlux[e0][e1][e0][f1]!=0) cout<<CovarianceFlux[e0][e1][e0][f1]<<endl;
	}
      }
      double RelativeError=TMath::Sqrt(ErrorSquared);
      if(Value!=0) RelativeError/=Value;      
      ErrorTotalFlux_RecAngle_Plus[d][e1]=RelativeError;
      ErrorTotalFlux_RecAngle_Minus[d][e1]=-RelativeError;
      ErrorSquared+=NominalSelectedMC_RecAngle[d]->GetBinError(e1+1);
      
      Error=TMath::Sqrt(ErrorSquared);

      cout<<"Flux error, Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
      boxErrorFlux_RecAngle[d][e1] = new TBox(NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle[d]->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxErrorFlux_RecAngle[d][e1]->SetFillColor(kBlue);
      NominalSelectedMC_RecAngle[d]->SetBinError(e1+1,ErrorSquared);
    }
    
    if(TotalValue[d]!=0) TotalRelErrorFlux=TMath::Sqrt(TotalErrorFluxSquared)/TotalValue[d];
  }

  //#########################################END FLUX ERROR################################



  
  //#####################################PUT BACK TOTAL ERROR TO SQRT (COVARIANCE)################################
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      double Value=0;double ErrorSquared=0;double Error=0;
      ErrorSquared=NominalSelectedMC[d]->GetBinError(e0+1,e1+1);
      Error=TMath::Sqrt(ErrorSquared);
      //Error=ErrorSquared;
      NominalSelectedMC[d]->SetBinError(e0+1,e1+1,Error);
      sliceNominalSelectedMC_RecMom[d][e1]->SetBinError(e0+1,Error);
      sliceNominalSelectedMC_RecAngle[d][e0]->SetBinError(e1+1,Error);
    }
  }

  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    ErrorSquared=NominalSelectedMC_RecMom[d]->GetBinError(e0+1);
    Error=TMath::Sqrt(ErrorSquared);
    //Error=ErrorSquared;
    NominalSelectedMC_RecMom[d]->SetBinError(e0+1,Error);
  }
  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
    double Value=0;double ErrorSquared=0;double Error=0;
    ErrorSquared=NominalSelectedMC_RecAngle[d]->GetBinError(e1+1);
    Error=TMath::Sqrt(ErrorSquared);
    //Error=ErrorSquared;
    NominalSelectedMC_RecAngle[d]->SetBinError(e1+1,Error);
  }
  //#####################################END OF PUT BACK TOTAL ERROR TO SQRT (COVARIANCE)################################


  
  //#####################################DRAWING OF 1D MOM & ANGLE PLOTS################################
  

  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
    csliceRecMom[d][e1]= new TCanvas(Form("csliceRecMom[%d][%d]",d,e1));
    sliceNominalSelectedMC_RecMom[d][e1]->Draw("E1");
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(EndError>=Systematics_Flux_Start) boxsliceErrorFlux_RecMom[d][e0][e1]->Draw("same");
      if(EndError>=Systematics_Xsec_Start) boxsliceErrorXS_RecMom[d][e0][e1]->Draw("same");
      boxsliceErrorStat_RecMom[d][e0][e1]->Draw("same");
    }
    sliceNominalSelectedMC_RecMom[d][e1]->Draw("E1same");

    if(Reconstructed) sliceNominalSelectedMC_RecMom[d][e1]->GetXaxis()->SetTitle("d_{#mu} (cm)");
    else sliceNominalSelectedMC_RecMom[d][e1]->GetXaxis()->SetTitle("p_{#mu} (GeV)");
   
    sliceNominalSelectedMC_RecMom[d][e1]->GetYaxis()->SetTitle("Number of events");
    sliceNominalSelectedMC_RecMom[d][e1]->SetLineColor(1);
    sliceNominalSelectedMC_RecMom[d][e1]->SetLineWidth(2);
    sliceNominalSelectedMC_RecMom[d][e1]->GetYaxis()->SetTitleOffset(1.3);
    csliceRecMom[d][e1]->Write(Form("Canvas_Nominal_RecMom%d_slice%d",d,e1));
  }

  cRecMom[d] = new TCanvas(Form("cRecMom[%d]",d));
  NominalSelectedMC_RecMom[d]->Draw("E1");
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    if(EndError>=Systematics_Flux_Start) boxErrorFlux_RecMom[d][e0]->Draw("same");
    if(EndError>=Systematics_Xsec_Start) boxErrorXS_RecMom[d][e0]->Draw("same");
    boxErrorStat_RecMom[d][e0]->Draw("same");
  }
  NominalSelectedMC_RecMom[d]->Draw("E1same");

  if(Reconstructed) NominalSelectedMC_RecMom[d]->GetXaxis()->SetTitle("d_{#mu} (cm)");
  else NominalSelectedMC_RecMom[d]->GetXaxis()->SetTitle("p_{#mu} (GeV)");
  
  NominalSelectedMC_RecMom[d]->GetYaxis()->SetTitle("Number of events");
  NominalSelectedMC_RecMom[d]->SetLineColor(1);
  NominalSelectedMC_RecMom[d]->SetLineWidth(2);
  NominalSelectedMC_RecMom[d]->GetYaxis()->SetTitleOffset(1.3);
  cRecMom[d]->Write(Form("Canvas_Nominal_RecMom_%d",d));
  
  cRecAngle[d] = new TCanvas(Form("cRecAngle[%d]",d));
  NominalSelectedMC_RecAngle[d]->Draw("E1");
  for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
    if(EndError>=Systematics_Flux_Start) boxErrorFlux_RecAngle[d][e0]->Draw("same");
      if(EndError>=Systematics_Xsec_Start) boxErrorXS_RecAngle[d][e0]->Draw("same");
    boxErrorStat_RecAngle[d][e0]->Draw("same");
  }
  NominalSelectedMC_RecAngle[d]->Draw("E1same");
  NominalSelectedMC_RecAngle[d]->SetLineColor(1);
  NominalSelectedMC_RecAngle[d]->SetLineWidth(2);
  NominalSelectedMC_RecAngle[d]->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  NominalSelectedMC_RecAngle[d]->GetYaxis()->SetTitle("Number of events");
  NominalSelectedMC_RecAngle[d]->GetYaxis()->SetTitleOffset(1.3);
  
  cRecAngle[d]->Write(Form("Canvas_Nominal_RecAngle_%d",d));
  //#####################################END OF DRAWING################################

  //#########################################ERROR TABLE####################################
  cout<<"2D error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      cout<<setprecision(2)<<std::fixed<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_Plus[d][e0][e1])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_Plus[d][e0][e1])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_Plus[d][e0][e1])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_Plus[d][e0][e1])<<"(stat),   ";
    }
    cout<<endl;
  }
  cout<<endl<<"1D Momentum error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecMom[d]->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecMom_Plus[d][e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecMom_Plus[d][e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_RecMom_Plus[d][e0])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_RecMom_Plus[d][e0])<<"(stat)"<<endl;
  }
  
  cout<<endl<<"1D Angle error table"<<endl;  
  for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecAngle[d]->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecAngle_Plus[d][e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecAngle_Plus[d][e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_RecAngle_Plus[d][e0])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_RecAngle_Plus[d][e0])<<"(stat)"<<endl;
  }

  cout<<endl<<"1bin result"<<endl;
  cout<<setprecision(2)<<std::fixed<<TotalValue[d]<<"+-"<<100*TMath::Abs(TotalRelErrorFlux)<<"%(Flux)+-"<<100*TMath::Abs(TotalRelErrorXS)<<"%(XS)+-"<<100*TMath::Abs(TotalRelErrorDet)<<"%(Det)+-"<<100*TMath::Abs(TotalRelErrorStat)<<"%(stat)"<<endl;

  if(EndError>=Systematics_Xsec_Start){
    cout<<"xscn systematics breakdown:"<<endl; 
    for(int i=StartXsec;i<=EndXsec;i++)
      cout<<setprecision(2)<<NameXSParam[i]<<" & "<<100*TMath::Abs(TotalRelxsError[d][i])<<"\\% \\\\"<<endl;
  }
  //#########################################END ERROR TABLE####################################


  //#########################################ERROR TABLE####################################
  cout<<"2D error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      cout<<setprecision(2)<<std::fixed<<NominalSelectedMC[d]->GetBinContent(e0+1,e1+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_Plus[d][e0][e1])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_Plus[d][e0][e1])<<"(XS)+-"<<TMath::Abs(ErrorTotalStatistics_Plus[d][e0][e1])<<"(stat),   ";
    }
    cout<<endl;
  }
  cout<<endl<<"1D Momentum error table"<<endl;  
  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecMom[d]->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecMom_Plus[d][e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecMom_Plus[d][e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalStatistics_RecMom_Plus[d][e0])<<"(stat)"<<endl;
  }
  
  cout<<endl<<"1D Angle error table"<<endl;  
  for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
    cout<<setprecision(2)<<std::fixed<<NominalSelectedMC_RecAngle[d]->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecAngle_Plus[d][e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecAngle_Plus[d][e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalStatistics_RecAngle_Plus[d][e0])<<"(stat)"<<endl;
  }
    
  //#########################################END ERROR TABLE####################################

    
  NominalSelectedMC_RecMom[d]->Write();
  NominalSelectedMC_RecAngle[d]->Write();
  NominalSelectedData[d]->Write();
  NominalSelectedMC[d]->Write();
  NominalTrueMC[d]->Write();
  NominalMCBkg[d]->Write();
  NominalMCBkg[d]->Divide(NominalSelectedMC[d]);
  NominalMCBkg[d]->Write("ContaminationMC");
#ifdef PLOTS
  //TCanvas * 
#endif
    
  }   

  if(NDetectors > 1){
    cout<<"xscn systematics breakdown:"<<endl; 
    for(int i=StartXsec;i<=EndXsec;i++)
      cout<<setprecision(2)<<NameXSParam[i]<<" & $"<<100*TMath::Abs(TotalRelxsError[0][i])<<" \\%$ & $"<<100*TMath::Abs(TotalRelxsError[1][i])<<"\\%$ \\\\"<<endl;
  }

     
  file->Close();
  return 0;
}

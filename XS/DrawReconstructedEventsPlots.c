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
#include <TMatrixDEigen.h>
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

//Step 0. comment #define DETECTORSYST if you wish to study straight XS and Flux systematics
#define DETECTORSYST 

//#define XSTABLE
//3 modes should not be mixed:
//1. Predictions MC: (Varied MC - varied bkg) - (Fixed MC - fixed bkg)
//2. Error to scale to data: (Fixed MC - varied bkg) - (Fixed MC - fixed bkg)
//3. Unfolded error: the way I constructed the code provide us:
//U.(Fixed MC - varied bkg) - U.(Fixed MC - fixed bkg)
// So, if we wish to evaluate the variation of MC prediction: this is the mode 0. For the real systematics, select the mode 1
const int MODE=1;

inline bool IsReducedXSSyst(int s){ // for model
  return (s==4 || s==5 || s==14 || s==15 || s==17 || s==19);
}

inline bool IsReadyDetSyst(int s){
  return (s==3 || s==4 || s==6 || s==5);
}

inline bool IsReadySyst(int s){
  return (IsReadyDetSyst(s) || s>Systematics_Detector_End || s<Systematics_Detector_Start);
}


void Evaluate1DError(TH2D * NominalMC, double **** CovarianceMatrix, TH1D * NominalMC_RecMom, TH1D * NominalMC_RecAngle){

  // activitates the use of trash bins in the true phase space
  bool TrueBinning=true;
#ifdef RECONSTRUCTED
  TrueBinning=false;
#endif

  int NBinsMom = NominalMC->GetNbinsX();
  int NBinsAngle = NominalMC->GetNbinsX();
  
  //First, create the momentum histogram value and error
  for(int e0=0;e0<NBinsMom;e0++){
    if(TrueBinning && IsTrueMomTrashBin[e0]) continue;
    double Value=0;
    for(int e1=0;e1<NBinsAngle;e1++) {
      if(TrueBinning && IsTrueAngleTrashBin[e1]) continue;
      Value+=NominalMC->GetBinContent(e0+1,e1+1);
    }
    double ErrorSquared=0;
    
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      if(TrueBinning && IsTrueAngleTrashBin[e1]) continue;
      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	if(TrueBinning && IsTrueAngleTrashBin[f1]) continue;
	ErrorSquared+= ((double) CovarianceMatrix[e0][e1][e0][f1]);//Case of momentum histogram
	//	cout<<CovarianceMatrix[e0][e1][e0][f1]<<endl;
      }
    }
    double Error=TMath::Sqrt(ErrorSquared);
    double RelativeError=Error / ( Value == 0 ? 1 : Value );
    NominalMC_RecMom->SetBinContent(e0+1,RelativeError);
    cout<<setprecision(3)<<"Momentum = " << NominalMC_RecMom->GetBinCenter(e0+1) << ", Error="<<RelativeError<<endl;
    //NominalMC_RecMom->SetBinError(e0+1,Error);
  }
 
  //Second, create the angle histogram value and error
  for(int e0=0;e0<NBinsAngle;e0++){
    double Value=0;
    if(TrueBinning && IsTrueAngleTrashBin[e0]) continue;
    for(int e1=0;e1<NBinsMom;e1++) {
      if(TrueBinning && IsTrueMomTrashBin[e1]) continue; 
      Value+=NominalMC->GetBinContent(e1+1,e0+1);
    }
    double ErrorSquared=0;
    
    for(int e1=0;e1<NBinsMom;e1++){//loop over effect 1
      if(TrueBinning && IsTrueMomTrashBin[e1]) continue; 
      for(int f1=0;f1<NBinsMom;f1++){//loop over effect 1
	if(TrueBinning && IsTrueMomTrashBin[f1]) continue; 
	ErrorSquared+=(double)CovarianceMatrix[e1][e0][f1][e0];//Case of angle histogram
      }
    }
    double Error=TMath::Sqrt(ErrorSquared);
    double RelativeError=Error / ( Value == 0 ? 1 : Value ); 
    NominalMC_RecAngle->SetBinContent(e0+1,RelativeError);
    cout<<setprecision(3)<<"Angle = " << NominalMC_RecAngle->GetBinCenter(e0+1) << ", Error="<<RelativeError<<endl;
    //NominalMC_RecAngle->SetBinError(e0+1,Error);
  }

  //Finally cout the total 1-bin error
  double Value=0, ErrorSquared=0;
  for(int e0=0;e0<NBinsMom;e0++){
    if(TrueBinning && IsTrueMomTrashBin[e0]) continue; 
    for(int e1=0;e1<NBinsAngle;e1++){
      if(TrueBinning && IsTrueAngleTrashBin[e1]) continue; 
      Value+=NominalMC->GetBinContent(e0+1,e1+1);
      for(int f0=0;f0<NBinsMom;f0++){
	if(TrueBinning && IsTrueMomTrashBin[f0]) continue; 
	for(int f1=0;f1<NBinsAngle;f1++){
	  if(TrueBinning && IsTrueAngleTrashBin[f1]) continue; 
	  ErrorSquared+=(double)CovarianceMatrix[e0][e1][f0][f1];  
	  //  cout<<CovarianceMatrix[e0][e1][f0][f1]<<endl;
	}
      }
    }
  }
  double Error=TMath::Sqrt(ErrorSquared);
  double RelativeError=Error / (Value==0 ? 1 : Value);
  cout<<"Total 1-bin relative error = "<<RelativeError<<" "<<ErrorSquared<<endl;
}



int main(int argc, char ** argv){

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  //TApplication theApp("App",0,0);
  gStyle->SetOptFit(kTRUE);
  char * OutputName = new char[256];
  bool UnfoldedPlots=false;

  bool PM=true;  
  bool ratio=false;
  bool water=false;
  
  bool MC=false;
  bool useKogaDetError=true;

  bool UseTrashBins=false;
  bool normalized=false;

  int FakeDataSet=0;
  bool FakeData=false;

  int nIter=1; 
  bool allIters=false;
      
  // 4 possible modes:
  //   -p -> PM
  //   -w -> WM
  //   -W -> water
  //   -r -> ratio

  int c=-1;
  while ((c = getopt(argc, argv, "o:wWrpmtnF:i:I")) != -1) {
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
    case 't':
      UseTrashBins=true;
      break;
    case 'n':
      normalized=true;
      break;
    case 'r':
      ratio=true;
      break;
    case 'W':
      water=true;
      break;
    case 'F':
      FakeData=true;
      FakeDataSet=atoi(optarg);
      break;
    case 'i':
      nIter=atoi(optarg);
      break;
    case 'I':
      allIters=true;
      break;
    }
  }

  if(ratio || water) {
    normalized=true;
    PM=false; //WM is used as first file
  }
  if(!FakeData) allIters=false;
  else UseTrashBins=true;

  // I only use the Load methods so I can initialize with any detector
  Xsec * XS = new Xsec(PM);
  XS->Xsec::Initialize();

  // ################### DECIDE WHETHER TO USE TRASH BINS ################
  //boolean to activate the use  of trash bin (in true phase space only)
#ifdef RECONSTRUCTED
  UseTrashBins=false;
  normalized=false;
  ratio=false;
  water=false;
  FakeData=false;
#endif


#ifdef RECONSTRUCTED
  const int NBinsMom=NBinsRecMom;
  const int NBinsAngle=NBinsRecAngle;
  double BinningMom[NBinsMom+1];
  double BinningAngle[NBinsAngle+1];
  for(int i=0;i<NBinsMom+1;i++) BinningMom[i]=BinningRecMom[i];
  for(int i=0;i<NBinsAngle+1;i++) BinningAngle[i]=BinningRecAngle[i];
  
  //double BinningMom[NBinsMom+1]=BinningRecMom;
  //double BinningAngle[NBinsAngle+1]=BinningRecAngle;
  
#else
   const int NBinsMom=NBinsTrueMom;
   const int NBinsAngle=NBinsTrueAngle;
  double BinningMom[NBinsMom+1];
  double BinningAngle[NBinsAngle+1];
  for(int i=0;i<NBinsMom+1;i++) BinningMom[i]=BinningTrueMom[i];
  for(int i=0;i<NBinsAngle+1;i++) BinningAngle[i]=BinningTrueAngle[i];
  //double BinningMom[NBinsMom+1]=BinningTrueMom;
  //double BinningAngle[NBinsAngle+1]=BinningTrueAngle;
#endif

  

 
  // Step 1. Declare and initialize the variables and tables
  bool merged=true;
  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  char merge[12]; sprintf(merge,(!PM && merged?"_merged":""));
  cout<<"Selected detector is "<<DetName<<endl;

  cout<<"CAREFUL: TO UNDERSTAND -> When XS error is varied, XS of 0 is different from the nominal MC that I have!!!!"<<endl; 
  char * txtDataName = new char[512];
  char * txtMCName = new char[512];
  double DataReconstructedEvents[NBinsMom][NBinsAngle];
  double MCReconstructedEvents[NBinsMom][NBinsAngle];//either rec or unfolded
  double MCTrueEvents[NBinsMom][NBinsAngle];//REALLY TRUE MC INFORMATION

  // for combined analysis (water or ratio)
  double MCReconstructedEvents_h2o[NBinsMom][NBinsAngle];
  double MCTrueEvents_h2o[NBinsMom][NBinsAngle];
  double MCReconstructedEvents_ch[NBinsMom][NBinsAngle];
  double MCTrueEvents_ch[NBinsMom][NBinsAngle];

  //For check of MC purity etc...
  double Efficiency[NBinsTrueMom][NBinsTrueAngle];
  double NominalEfficiency[NBinsTrueMom][NBinsTrueAngle];
  double NumberOfPOT;
  double MCReconstructedBkgEvents[NBinsMom][NBinsAngle];//whether rec, whether unfolded
  double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsMom][NBinsAngle];
  TH2D * NominalMCBkg = new TH2D("NominalMCBkg","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);

  TH1D * Chi2Stat=new TH1D("Chi2Stat","#chi^{2} stat",nIter+1,0,nIter+1);
  TH1D * Chi2Bias=new TH1D("Chi2Bias","#chi^{2} bias",nIter+1,0,nIter+1);
  TH1D * Chi2BiasCorr=new TH1D("Chi2BiasCorr","#chi^{2} bias+corr",nIter+1,0,nIter+1);
  TH1D * Chi2Rel=new TH1D("Chi2Relative","#chi^{2} relative",nIter+1,0,nIter+1);
  TH1D * Chi2RelCorr=new TH1D("Chi2RelativeCorr","#chi^{2} rel.+corr",nIter+1,0,nIter+1);
  TH1D * Chi2Data=new TH1D("Chi2Data","",nIter+1,0,nIter+1);
  TH1D * Chi2DataCorrr=new TH1D("Chi2DataCorr","",nIter+1,0,nIter+1);

  double PreviousIterUnfolded[NBinsTrueMom][NBinsTrueAngle];
  double InfiniteIterUnfolded[NBinsTrueMom][NBinsTrueAngle];

  int iIter=nIter,fIter=nIter+1; 
  if(allIters) iIter=0;

  for(int iter=iIter;iter<fIter;iter++){

    // all Covariance matrices named ...Targets are global Cov. matrices with h2o and ch mixed.
    double ***** CovarianceReduced;
    CovarianceReduced = new double ****[EndError+1];

    for(int h=0;h<EndError+1;h++) {
      CovarianceReduced[h] = new double ***[NBinsMom];
        }

    for(int h=0;h<EndError+1;h++){ 
      for(int i=0;i<NBinsMom;i++){
	CovarianceReduced[h][i] = new double **[NBinsAngle];
      }
    }
    for(int h=0;h<EndError+1;h++){ 
      for(int i=0;i<NBinsMom;i++){
	for(int j=0;j<NBinsAngle;j++){
	  CovarianceReduced[h][i][j] = new double *[NBinsMom];
	}
      }
    }
    for(int h=0;h<EndError+1;h++){ 
      for(int i=0;i<NBinsMom;i++){
	for(int j=0;j<NBinsAngle;j++){
	  for(int k=0;k<NBinsMom;k++){
	    CovarianceReduced[h][i][j][k] = new double [NBinsAngle];
	  }
	}
      }
    }
  
    double CovarianceReducedTargets[EndError+1][NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{{0}}}}}}};

    double CovarianceStatistics[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};
    double CovarianceFlux[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};
    double CovarianceBirksPlus[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};
    double CovarianceBirksMinus[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};
    double CovarianceTotal[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};

    double CovarianceStatisticsTargets[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{0}}}}}};
    double CovarianceFluxTargets[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{0}}}}}};
    double CovarianceBirksPlusTargets[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{0}}}}}};
    double CovarianceBirksMinusTargets[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{0}}}}}};
    double CovarianceTotalTargets[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{0}}}}}};

    double CorrelationReduced[EndError+1][NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{{0}}}}};
    double CorrelationStatistics[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};
    double CorrelationFlux[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle]={{{{0}}}};

    double **** CovarianceXS;
    double CovarianceXSTargets[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle][2][2]={{{{{{0}}}}}};;
    CovarianceXS = new double ***[NBinsMom];
    for(int i=0;i<NBinsMom;i++) {
      CovarianceXS[i] = new double **[NBinsAngle];
    }
    for(int i=0;i<NBinsMom;i++){
      for(int j=0;j<NBinsAngle;j++){
	CovarianceXS[i][j] = new double *[NBinsMom];
      }
    }
    for(int i=0;i<NBinsMom;i++){
      for(int j=0;j<NBinsAngle;j++){
	for(int k=0;k<NBinsMom;k++){
	  CovarianceXS[i][j][k] = new double [NBinsAngle];
	}
      }
    }
    //double CovarianceXS[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle];
    double CorrelationXS[NBinsMom][NBinsAngle][NBinsMom][NBinsAngle];
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	int bin1=NBinsAngle*e0+e1;
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    int bin2=NBinsAngle*f0+f1;
	    CovarianceXS[e0][e1][f0][f1]=0;
	  }
	}
      }
    }

  

    for(int s0=0;s0<=EndError;s0++){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  int bin1=NBinsAngle*e0+e1;
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      int bin2=NBinsAngle*f0+f1;
	      CovarianceReduced[s0][e0][e1][f0][f1]=0;
	      CorrelationReduced[s0][e0][e1][f0][f1]=0;
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
    TH2D * NominalSelectedMC_ch = new TH2D("NominalSelectedMC_ch","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    TH2D * NominalSelectedMC_h2o = new TH2D("NominalSelectedMC_h2o","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    TH2D * NominalSelectedData = new TH2D("NominalSelectedData","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    TH1D * NominalSelectedMC_RecMom = new TH1D("NominalSelectedMC_RecMom","d_{#mu} distribution of selected events",NBinsMom,BinningMom);
#ifndef RECONSTRUCTED
    NominalSelectedMC_RecMom->SetTitle("p_{#mu} distribution of selected events");
#endif
    TH1D * NominalSelectedMC_RecAngle = new TH1D("NominalSelectedMC_RecAngle","#theta_{#mu} distribution of selected events",NBinsAngle,BinningAngle);
    TH2D * NominalSelectedMC_XSTemp = new TH2D("NominalSelectedMC_XSTemp","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);

    TH2D * NominalTrueMC = new TH2D("NominalTrueMC","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    TH2D * NominalTrueMC_ch = new TH2D("NominalTrueMC_ch","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
    TH2D * NominalTrueMC_h2o = new TH2D("NominalTrueMC_h2o","",NBinsMom,BinningMom,NBinsAngle,BinningAngle);

    TH2D * Error_Minus[EndError+1]; TH2D * Error_Plus[EndError+1];
    TH2D * Error_Minus_ch[EndError+1]; TH2D * Error_Plus_ch[EndError+1];
    TH2D * Error_Minus_h2o[EndError+1]; TH2D * Error_Plus_h2o[EndError+1];

    TH1D * Error_RecMom[EndError+1];
    TH1D * Error_RecAngle[EndError+1];

    double ErrorTotalStatistics_Plus[NBinsMom][NBinsAngle];
    double ErrorTotalStatistics_Minus[NBinsMom][NBinsAngle];
    double ErrorTotalStatistics_RecMom_Plus[NBinsMom];
    double ErrorTotalStatistics_RecMom_Minus[NBinsMom];
    double ErrorTotalStatistics_RecAngle_Plus[NBinsAngle];
    double ErrorTotalStatistics_RecAngle_Minus[NBinsAngle];

    //ML 2017/07/18 added
    double ErrorTotalDetector_Plus[NBinsMom][NBinsAngle]={{0.}};
    double ErrorTotalDetector_Minus[NBinsMom][NBinsAngle]={{0.}};
    double ErrorTotalDetector_RecMom_Plus[NBinsMom]={0};
    double ErrorTotalDetector_RecMom_Minus[NBinsMom]={0};
    double ErrorTotalDetector_RecAngle_Plus[NBinsAngle]={0};
    double ErrorTotalDetector_RecAngle_Minus[NBinsAngle]={0};

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
    TH1D * ErrorEffFlux[NBinsTrueMom][NBinsTrueAngle];
    TF1 * fErrorEffFlux[NBinsTrueMom][NBinsTrueAngle];

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
    // ML (tmp?) for errors with some invalid points 
    double xXS_reduced[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations-2];
    double yXS_reduced[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations-2];

    TSpline3 * sErrorXS_ch[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
    TGraph * gErrorXS_ch[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
    TSpline3 * sErrorXS_h2o[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
    TGraph * gErrorXS_h2o[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle];
    double xXS_ch[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];
    double yXS_ch[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];
    double xXS_h2o[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];
    double yXS_h2o[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations];
    double xXS_reduced_ch[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations-2];
    double yXS_reduced_ch[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations-2];
    double xXS_reduced_h2o[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations-2];
    double yXS_reduced_h2o[Systematics_Xsec_End-Systematics_Xsec_Start+1][NBinsMom][NBinsAngle][NXsecVariations-2];


    TAxis * aX_NominalSelectedMC = (TAxis*) NominalSelectedMC->GetXaxis();
    TAxis * aY_NominalSelectedMC = (TAxis*) NominalSelectedMC->GetYaxis();
    TAxis * aX_NominalSelectedData = (TAxis*) NominalSelectedData->GetXaxis();
    TAxis * aY_NominalSelectedData = (TAxis*) NominalSelectedData->GetYaxis();
    TAxis * aX_NominalSelectedMC_RecMom = (TAxis*) NominalSelectedMC_RecMom->GetXaxis();
    TAxis * aY_NominalSelectedMC_RecMom = (TAxis*) NominalSelectedMC_RecMom->GetYaxis();
    TAxis * aX_NominalSelectedMC_RecAngle = (TAxis*) NominalSelectedMC_RecAngle->GetXaxis();
    TAxis * aY_NominalSelectedMC_RecAngle = (TAxis*) NominalSelectedMC_RecAngle->GetYaxis();

    double TotalValue=0., TotalValue_ch=0., TotalValue_h2o=0.;

    double RelativeErrorWM[Systematics_Detector_End+1]={0,0,0.10,0,0,0,0,0.6,0.82,0.41,0.31,0.13,0,0.35,.60,0};
    double RelativeErrorPM[Systematics_Detector_End+1]={0,0,0.12,0,0,0,0,0.58,0.58,0,0.18,0.12,0,0.97,.27,0};

    TH1D* h_test=new TH1D("gaus_test","gaus_test",100,-5,5);
    TH1D* h_test_red=new TH1D("gaus_test_red","gaus_test",100,-5,5);

    // Step 2. Loop over all the errors to determine the variation of number of events for: each error source & bin
    for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      //Allow to fill the nominal distribution even if there is no detector systematics (and we wish to start error evaluation at 17, where XS starts)
      cout<<"The error currently tested is number "<<ErrorType<<endl;

#ifndef DETECTORSYST
      if(ErrorType>=Systematics_Detector_Start && ErrorType<=Systematics_Detector_End) continue;
      //    if(ErrorType>=Systematics_Detector_Start && ErrorType<=16) continue;
#endif

        
      //##############################HISTOGRAM INITIALIZATION#########################################
      //Error_Minus and _Plus are used for simple cov matrix; filled with RELATIVE ERRORS
      Error_Minus[ErrorType] = new TH2D(Form("Error_Minus[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      Error_Plus[ErrorType] = new TH2D(Form("Error_Plus[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      //The following are used to build the 2targets cov matrix; filled with ABSOLUTE ERRORS
      Error_Minus_ch[ErrorType] = new TH2D(Form("Error_Minus_ch[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      Error_Plus_ch[ErrorType] = new TH2D(Form("Error_Plus_ch[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      Error_Minus_h2o[ErrorType] = new TH2D(Form("Error_Minus_h2o[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);
      Error_Plus_h2o[ErrorType] = new TH2D(Form("Error_Plus_h2o[%d]",ErrorType),"",NBinsMom,BinningMom,NBinsAngle,BinningAngle);

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
	  else if(ErrorType==Systematics_Flux_Start/* && ErrorType<=Systematics_Flux_End*/){
	    ErrorFlux[e0][e1] = new TH1D(Form("ErrorFlux%d_%d",e0,e1),"",600,-3,3);
	    fErrorFlux[e0][e1] = new TF1(Form("fErrorFlux[%d][%d]",e0,e1),"gaus",-3,3);
	  }
	  else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	    ErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TH1D(Form("ErrorXS[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),"",620,-3.1,3.1);
	    ErrorXS_Norm[ErrorType-Systematics_Xsec_Start][e0][e1] = new TH1D(Form("ErrorXS_Norm[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),"",620,-3.1,3.1);
	    //sErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = ?????new TF1(Form("fErrorFlux[%d][%d]",e0,e1),"gaus",-3,3);
	  }
	}
      }
#ifdef RECONSTRUCTED    
      if(ErrorType==0){
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over effect 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over effect 1
	    ErrorEffFlux[c0][c1] = new TH1D(Form("ErrorEffFlux%d_%d",c0,c1),"",600,-3,3);
	    fErrorEffFlux[c0][c1] = new TF1(Form("fErrorEffFlux[%d][%d]",c0,c1),"gaus",-3,3);
	  }
	} 
      }
#endif
      //##############################END OF INITIALIZATION#########################################
   


      //#################################TEMPORARY, SINCE THE XS VARIATION OF 0 SIGMA DOES NOT CORRESPONDS TO THE NOMINAL MC, WE REDEFINE NOMINAL ONLY FOR XS ERROR AS THE 0 SIGMA VARIATION#################################
      // ----------------------------------- ML 2017/08/31 for me it corresponds ------------------------
      /*
	if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	int n=3;

	#ifdef RECONSTRUCTED
	sprintf(txtMCName,"%s/XS/files/MCSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	//if(MODE == 1){
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	MCReconstructedEvents[e0][e1] -= MCReconstructedBkgEvents[e0][e1];
	}
	}
	//}
      
	//XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
	#else
	sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s_Systematics%d_%d.root",cINSTALLREPOSITORY,DetName,ErrorType,n);
	XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents);
	//XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	//XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
	#endif

	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	NominalSelectedMC_XSTemp->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	}
	}
	}
      */
      //#################################################################"


      //##############################SELECTED ERRORS#########################################
      if (FakeData && ErrorType>1) continue;
      //      if (ErrorType>=18 && ErrorType<=22) continue;
#ifdef RECONSTRUCTED
      if (ErrorType==1) continue;
 
#ifdef DETECTORSYST
      if(useKogaDetError && !IsReadySyst(ErrorType)){
	double RelativeBkgError=(PM?RelativeErrorPM[ErrorType]:RelativeErrorWM[ErrorType])/100.;
	for(int e0=0;e0<NBinsMom;e0++){
	  for(int e1=0;e1<NBinsAngle;e1++){
	    //MC = SIGNAL Nominal MC - BKG Varied MC. Note that SIGNAL Nominal MC = SELECTED Nominal MC + BKG Nominal MC.  
	    double Value1 = (NominalSelectedMC->GetBinContent(e0+1,e1+1)+NominalMCBkg->GetBinContent(e0+1,e1+1)) - NominalMCBkg->GetBinContent(e0+1,e1+1)*(1+RelativeBkgError);
	    double Error1=Value1-NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    //double Value1=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-fabs(Error1));
	    Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,+fabs(Error1));
	    for(int f0=0;f0<NBinsMom;f0++){
	      for(int f1=0;f1<NBinsAngle;f1++){
		double Value2 = (NominalSelectedMC->GetBinContent(f0+1,f1+1)+NominalMCBkg->GetBinContent(f0+1,f1+1)) - NominalMCBkg->GetBinContent(f0+1,f1+1)*(1+RelativeBkgError);
		double Error2=Value2-NominalSelectedMC->GetBinContent(f0+1,f1+1);
		CovarianceReduced[ErrorType][e0][e1][f0][f1]=Error1*Error2;
	      } 
	    }
	  }
	}
	continue;
      }
#endif

#else
#ifdef DETECTORSYST
      if(ErrorType==2 && !FakeData) NE[ErrorType]=2; // Scaled errors put in this slot
      else if(!IsReadySyst(ErrorType)) continue;
#else
      if(!IsReadySyst(ErrorType)) continue;
#endif
#endif

      //##############################LOADING: FILL THE HISTOGRAMS#########################################
      for(int n=0;n<NE[ErrorType];n++){
	double ErrorValue=Start[ErrorType]+n*Step[ErrorType];

	//if(ErrorType==5 && n==0) continue; // to select the highest 1-bin error
	// this selection should be done on a bin-by-bin basis though

	//The variation of Xsec parameter, in #sigma:
	double XsecVariation=ErrorValue-(ErrorType-Systematics_Xsec_Start)*NXsecVariations-CenterXsecVariations; //. A number between 0 and 175 - the center of the current systematic source (nominal). For example, for Xsec error source #10, it starts from 7*(10-1)=63 and ends at 70. from 63 to 70, it contains the variariation of -3,-2,-1,0,1,2,3 sigma respectively. The center is then located at 66. For the example of a 2 sigma variation, the substraction will be therefore equal to: 68-66=2, which gives the number of sigmas!      
	if(ErrorType>=Systematics_Xsec_Start){
	  // some troublesome Xsec values where only +-1 sigma is meaningfull
	  if(IsReducedXSSyst(ErrorType-Systematics_Xsec_Start) && abs(XsecVariation)>1) continue;
	}
      
#ifdef RECONSTRUCTED
	if(ErrorType==5) sprintf(txtMCName,"%s/XS/files/MCSelected_%s%s_Systematics%d_%d_rescaled.txt",cINSTALLREPOSITORY,DetName,merge,ErrorType,n);
	else if(ErrorType==4) sprintf(txtMCName,"%s/XS/files/MCSelected_%s%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,merge,ErrorType,n);
	else sprintf(txtMCName,"%s/XS/files/MCSelected_%s%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,merge,ErrorType,n);

	//      sprintf(txtMCName,"%s/XS/files/MCSelected_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,ErrorType,n,DetName);
	//sprintf(txtMCName,"%s/XS/Selection1000_cutBkg%s.txt",cINSTALLREPOSITORY,cINSTALLREPOSITORY,DetName);
	XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,MCReconstructedEvents);
      
	if(ErrorType>=7 && ErrorType<=Systematics_Detector_End) sprintf(txtDataName,"/home/bquilain/CC0pi_XS/XS/files/DataSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else if(ErrorType==0 && EndError>=7 && StartError<=Systematics_Detector_End) sprintf(txtDataName,"/home/bquilain/CC0pi_XS/XS/files/DataSelected_%s_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	else sprintf(txtDataName,"%s/XS/files/DataSelected_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	//else sprintf(txtDataName,"%s/XS/Selection1000_cutBkg%s.txt",cINSTALLREPOSITORY,cINSTALLREPOSITORY,DetName);
	//cout<<"good"<<endl;
	if(MC){
	  //XS->Xsec::LoadInputFiles_OnlySelectedData(txtMCName,DataReconstructedEvents);
	  //if(ErrorType==0)
	  XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,&NumberOfPOT);
	}
	else{
	  //XS->Xsec::LoadInputFiles_OnlySelectedData(txtDataName,DataReconstructedEvents);
	  //if(ErrorType==0)
	  XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,&NumberOfPOT);
	}
      
	if(ErrorType==0) {
	  for(int c0=0;c0<NBinsTrueMom;c0++){
	    for(int c1=0;c1<NBinsTrueAngle;c1++){
	      NominalEfficiency[c0][c1] = NumberOfPOT*Efficiency[c0][c1];
	    }
	  }
	}

#else 
	char FD[5];sprintf(FD,FakeData?Form("FD%d_",FakeDataSet):"");
	if(water || ratio){
	  // the first file loaded is WM
	  if(ErrorType==5) sprintf(txtMCName,"%s/XS/files/MCUnfolded_%sWM_merged_Systematics%d_%d_rescaled.root",cINSTALLREPOSITORY,FD,ErrorType,n);
	  else if(ErrorType==2 && useKogaDetError) sprintf(txtMCName,"%s/XS/files/MCUnfolded_%sWM_merged_SystematicsScaled1.34_%d.root",cINSTALLREPOSITORY,FD,n);
	  else sprintf(txtMCName,"%s/XS/files/MCUnfolded_%sWM_merged_Systematics%d_%d.root",cINSTALLREPOSITORY,FD,ErrorType,n);
	  
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents_h2o,MCTrueEvents_h2o,normalized,iter);

	  // the second file loaded is PM
	  if(ErrorType==5) sprintf(txtMCName,"%s/XS/files/MCUnfolded_PM_Systematics%d_%d_rescaled.root",cINSTALLREPOSITORY,ErrorType,n);
	  else if(ErrorType==2 && useKogaDetError) sprintf(txtMCName,"%s/XS/files/MCUnfolded_PM_SystematicsScaled%1.2f_%d.root",cINSTALLREPOSITORY,1.32,n);
	  else sprintf(txtMCName,"%s/XS/files/MCUnfolded_PM_Systematics%d_%d.root",cINSTALLREPOSITORY,ErrorType,n);

	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents_ch,MCTrueEvents_ch,normalized,iter);
	}
	else{
	  if(ErrorType==5) sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s%s%s_Systematics%d_%d_rescaled.root",cINSTALLREPOSITORY,FD,DetName,merge,ErrorType,n);
	  else if(ErrorType==2 && useKogaDetError) sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s%s%s_SystematicsScaled%1.2f_%d.root",cINSTALLREPOSITORY,FD,DetName,merge,(PM?1.32:1.34),n);
	  else sprintf(txtMCName,"%s/XS/files/MCUnfolded_%s%s%s_Systematics%d_%d.root",cINSTALLREPOSITORY,FD,DetName,merge,ErrorType,n);
	  
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents,normalized,iter);
	}

	

	/* ML not used yet    
	   if(ErrorType>=7 && ErrorType<=Systematics_Detector_End) sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	   else if(ErrorType==0 && EndError>=7 && StartError<=Systematics_Detector_End) sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d%s.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	   else sprintf(txtDataName,"%s/XS/files/DataUnfolded_%s_Systematics%d_%d.txt",cINSTALLREPOSITORY,DetName,ErrorType,n);
	*/
	//else sprintf(txtDataName,"%s/XS/Selection1000_cutBkg%s.txt",cINSTALLREPOSITORY,cINSTALLREPOSITORY,DetName);
	//cout<<"good"<<endl;
	/*if(MC){
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtMCName,MCReconstructedEvents,MCTrueEvents,normalized);
	  //if(ErrorType==0) XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	  }
	  else{
	  XS->Xsec::LoadInputFiles_OnlyUnfoldedData(txtDataName,DataReconstructedEvents,MCTrueEvents,normalized);
	  //if(ErrorType==0) XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);
	  }*/
#endif
 


      
    
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
#ifdef RECONSTRUCTED
	    //cout<<"mom="<<e0<<" ang="<<e1<<" data="<<MCReconstructedEvents[e0][e1]<<" bkg="<<MCReconstructedBkgEvents[e0][e1]<<" signal="<<MCReconstructedEvents[e0][e1] - MCReconstructedBkgEvents[e0][e1]<<endl;
	    MCReconstructedEvents[e0][e1] = MCReconstructedEvents[e0][e1] - MCReconstructedBkgEvents[e0][e1];
	    DataReconstructedEvents[e0][e1] = DataReconstructedEvents[e0][e1] - MCReconstructedBkgEvents[e0][e1];
	    //MCReconstructedEvents[e0][e1] = MCReconstructedEvents[e0][e1] - NominalMCBkg->GetBinContent(e0+1,e1+1);
	    //DataReconstructedEvents[e0][e1] = DataReconstructedEvents[e0][e1] - NominalMCBkg->GetBinContent(e0+1,e1+1);
	  
	    if( (MODE==1) && ErrorType>1 ){
	      //MC = SIGNAL Nominal MC - BKG Varied MC. Note that SIGNAL Nominal MC = SELECTED Nominal MC + BKG Nominal MC.  
	      MCReconstructedEvents[e0][e1] = (NominalSelectedMC->GetBinContent(e0+1,e1+1)+NominalMCBkg->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	      DataReconstructedEvents[e0][e1] = (NominalSelectedData->GetBinContent(e0+1,e1+1)+NominalMCBkg->GetBinContent(e0+1,e1+1)) - MCReconstructedBkgEvents[e0][e1];
	    }
#else
	    if(water || ratio){
	      // WM --> water
	      MCReconstructedEvents_h2o[e0][e1] = (NTargetWM*MCReconstructedEvents_h2o[e0][e1]-NTargetWM_sci*MCReconstructedEvents_ch[e0][e1])/NTargetWM_h2o;
	      MCTrueEvents_h2o[e0][e1] = (NTargetWM*MCTrueEvents_h2o[e0][e1]-NTargetWM_sci*MCTrueEvents_ch[e0][e1])/NTargetWM_h2o;

	      if(ratio){
		MCReconstructedEvents[e0][e1]=MCReconstructedEvents_h2o[e0][e1]/MCReconstructedEvents_ch[e0][e1];
		MCTrueEvents[e0][e1]=MCTrueEvents_h2o[e0][e1]/MCTrueEvents_ch[e0][e1];
	      }
	      if(water){
		MCReconstructedEvents[e0][e1]=MCReconstructedEvents_h2o[e0][e1];
		MCTrueEvents[e0][e1]=MCTrueEvents_h2o[e0][e1];
	      }
	    }
#endif
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
	      NominalSelectedMC->SetBinContent(e0+1,e1+1,MCReconstructedEvents[e0][e1]);
	      //ML for consistency, the error should be initialized to ErrorSquared=Value; (reco case)
	      NominalSelectedMC->SetBinError(e0+1,e1+1,MCReconstructedEvents[e0][e1]);

	      NominalSelectedMC_ch->SetBinContent(e0+1,e1+1,MCReconstructedEvents_ch[e0][e1]);
	      NominalSelectedMC_h2o->SetBinContent(e0+1,e1+1,MCReconstructedEvents_h2o[e0][e1]);

	      NominalTrueMC->SetBinContent(e0+1,e1+1,MCTrueEvents[e0][e1]);//not filled ifdef RECONSTRUCTED

	      NominalTrueMC_ch->SetBinContent(e0+1,e1+1,MCTrueEvents_ch[e0][e1]);
	      NominalTrueMC_h2o->SetBinContent(e0+1,e1+1,MCTrueEvents_h2o[e0][e1]);

	      NominalSelectedData->SetBinContent(e0+1,e1+1,DataReconstructedEvents[e0][e1]);
	      NominalMCBkg->SetBinContent(e0+1,e1+1,MCReconstructedBkgEvents[e0][e1]);

	      if(!UseTrashBins ||( !IsTrueMomTrashBin[e0] && !IsTrueAngleTrashBin[e1])){
		if(ratio){
		  TotalValue_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
		  TotalValue_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
		}
		else TotalValue+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	      }
	    }
	  
	    double RelativeValue = (MCReconstructedEvents[e0][e1] - NominalSelectedMC->GetBinContent(e0+1,e1+1));
	    if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeValue /= NominalSelectedMC->GetBinContent(e0+1,e1+1);

	    double AbsoluteValue_ch = (MCReconstructedEvents_ch[e0][e1] - NominalSelectedMC_ch->GetBinContent(e0+1,e1+1));
	    double AbsoluteValue_h2o = (MCReconstructedEvents_h2o[e0][e1] - NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1));
	    double RelativeValue_ch = AbsoluteValue_ch/(NominalSelectedMC_ch->GetBinContent(e0+1,e1+1)!=0?NominalSelectedMC_ch->GetBinContent(e0+1,e1+1):1);
	    double RelativeValue_h2o = AbsoluteValue_h2o/(NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1)!=0?NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1):1);
	    //double RelativeValueData = (DataReconstructedEvents[e0][e1] - NominalSelectedData->GetBinContent(e0+1,e1+1));
	    //if(NominalSelectedData->GetBinContent(e0+1,e1+1)!=0) RelativeValue /= NominalSelecteOAdData->GetBinContent(e0+1,e1+1);
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
		    if(ratio){
		      CovarianceStatisticsTargets[e0][e1][f0][f1][0][0]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceStatisticsTargets[e0][e1][f0][f1][0][1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		      CovarianceStatisticsTargets[e0][e1][f0][f1][1][0]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceStatisticsTargets[e0][e1][f0][f1][1][1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		    }
		    CovarianceStatistics[e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		}
	      }
	    }
	    if(ErrorType==2){
#ifdef RECONSTRUCTED
	      ErrorNoise[e0][e1]->Fill(ErrorValue,MCReconstructedEvents[e0][e1]);
	      ErrorNoise_Norm[e0][e1]->Fill(ErrorValue);
#else
	      if(MCReconstructedEvents[e0][e1]>=NominalSelectedMC->GetBinContent(e0+1,e1+1)) Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	      else Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-RelativeValue);

	      if(ratio){
		if(AbsoluteValue_ch>0) Error_Plus_ch[ErrorType]->SetBinContent(e0+1,e1+1,AbsoluteValue_ch);
		else Error_Minus_ch[ErrorType]->SetBinContent(e0+1,e1+1,-AbsoluteValue_ch);

		if(AbsoluteValue_h2o>0) Error_Plus_h2o[ErrorType]->SetBinContent(e0+1,e1+1,AbsoluteValue_h2o);
		else Error_Minus_h2o[ErrorType]->SetBinContent(e0+1,e1+1,-AbsoluteValue_h2o);
	      }

#endif
	      //cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;
	    }
	    else if(ErrorType==3){
	      Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	      Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));

	      ErrorHitEfficiency[e0][e1]->Fill(RelativeValue);//FINALLY USELESS
	      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
		for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		  //Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(DataReconstructedEvents[e0][e1]-NominalSelectedData->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-Nominal);
		  //CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(DataReconstructedEvents[e0][e1]-NominalSelectedData->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-Nominal);cout<<Nominal<<"f"<<endl;
		  CovarianceReduced[ErrorType][e0][e1][f0][f1]=(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1)); // ML new 2017/07/20 -- I put here an expression similar to 4
		  if(ratio){
		    CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][0]=(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		    CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][1]=(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		  }
		}
	      }
	    }
	    else if(ErrorType==4){
	      //if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) cout<<"NominalSelectedMC="<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", Varied="<<MCReconstructedEvents[e0][e1]<<endl;
	      Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	      Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));
	      //  cout<<e0<<" "<<e1<<" "<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<" "<<MCReconstructedEvents[e0][e1]<<" "<<RelativeValue<<endl;
	      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
		for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		  //Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		  // WARNING *** ml 2017/07/20 NE[ErrorType] is 1 here !!!
		  CovarianceReduced[ErrorType][e0][e1][f0][f1]=(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		  if(ratio){
		    CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][0]=(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		    CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][1]=(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		  }
		
		}
	      }	    
	      //NEvents[NBinsMom][NBinsAngle]->Fill(e0,e1,MCReconstructedEvents[e0][e1]);
	    }
	    else if(ErrorType==5){
	      if(MCReconstructedEvents[e0][e1]>=NominalSelectedMC->GetBinContent(e0+1,e1+1)){
		//if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
		Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	      }
	      else Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-RelativeValue);

	      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
		for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		  //Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));

		  //CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		  
		  if(ratio){
		    if(n==1){
		      CovarianceBirksPlusTargets[e0][e1][f0][f1][0][0]=(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceBirksPlusTargets[e0][e1][f0][f1][0][1]=(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		      CovarianceBirksPlusTargets[e0][e1][f0][f1][1][0]=(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceBirksPlusTargets[e0][e1][f0][f1][1][1]=(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		    }
		    else {
		      CovarianceBirksMinusTargets[e0][e1][f0][f1][0][0]=(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceBirksMinusTargets[e0][e1][f0][f1][0][1]=(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		      CovarianceBirksMinusTargets[e0][e1][f0][f1][1][0]=(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceBirksMinusTargets[e0][e1][f0][f1][1][1]=(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		    }
		  }
		  else {
		    if(n==1) CovarianceBirksPlus[e0][e1][f0][f1]=(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		    else CovarianceBirksMinus[e0][e1][f0][f1]=(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		  }

		}
	      }	    
	      //cout<<"Birks attenuation length error is changed"
	      //NEvents[NBinsMom][NBinsAngle]->Fill(e0,e1,MCReconstructedEvents[e0][e1]);
	    }
	    else if(ErrorType==6){
	      if(MCReconstructedEvents[e0][e1]>NominalSelectedMC->GetBinContent(e0+1,e1+1)) Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
	      else if(MCReconstructedEvents[e0][e1]<NominalSelectedMC->GetBinContent(e0+1,e1+1))Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-RelativeValue);

	      if(ratio){
		if(AbsoluteValue_ch>0) Error_Plus_ch[ErrorType]->SetBinContent(e0+1,e1+1,AbsoluteValue_ch);
		else Error_Minus_ch[ErrorType]->SetBinContent(e0+1,e1+1,-AbsoluteValue_ch);

		if(AbsoluteValue_h2o>0) Error_Plus_h2o[ErrorType]->SetBinContent(e0+1,e1+1,AbsoluteValue_h2o);
		else Error_Minus_h2o[ErrorType]->SetBinContent(e0+1,e1+1,-AbsoluteValue_h2o);
	      }
	      /*for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
		for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		}
		}*/	    
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
	      if(ErrorType==16){
		//cout<<"Bin "<<e0<<", "<<e1<<", Nominal MC="<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", Modified value="<<MCReconstructedEvents[e0][e1]<<", relative="<<RelativeValue<<endl;
		ErrorFlux[e0][e1]->Fill(RelativeValue);
		//cout<<e0<<","<<e1<<"="<<MCReconstructedEvents[e0][e1]<<endl;
		for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
		  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		    //Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		    //CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		    //cout<<std::scientific<<setprecision(3)<<"("<<e0<<","<<e1<<")     "<<MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1)<<" vs "<<MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0,f1)<<endl;
		    if(ratio){
		      CovarianceFluxTargets[e0][e1][f0][f1][0][0]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceFluxTargets[e0][e1][f0][f1][0][1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_h2o[e0][e1]-NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		      CovarianceFluxTargets[e0][e1][f0][f1][1][0]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_h2o[f0][f1]-NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1));
		      CovarianceFluxTargets[e0][e1][f0][f1][1][1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents_ch[e0][e1]-NominalSelectedMC_ch->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents_ch[f0][f1]-NominalSelectedMC_ch->GetBinContent(f0+1,f1+1));
		    }
		    CovarianceFlux[e0][e1][f0][f1]+=(1./(NE[ErrorType]-1.))*(MCReconstructedEvents[e0][e1]-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		  }
		}
	      }
	      else if(ErrorType==17){ // Antinu, nue bkg
		if(MCReconstructedEvents[e0][e1]>NominalSelectedMC->GetBinContent(e0+1,e1+1)) Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,RelativeValue);
		else if (MCReconstructedEvents[e0][e1]<NominalSelectedMC->GetBinContent(e0+1,e1+1)) Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-RelativeValue);

		if(ratio){
		  if(AbsoluteValue_ch>0) Error_Plus_ch[ErrorType]->SetBinContent(e0+1,e1+1,AbsoluteValue_ch);
		  else Error_Minus_ch[ErrorType]->SetBinContent(e0+1,e1+1,-AbsoluteValue_ch);

		  if(AbsoluteValue_h2o>0) Error_Plus_h2o[ErrorType]->SetBinContent(e0+1,e1+1,AbsoluteValue_h2o);
		  else Error_Minus_h2o[ErrorType]->SetBinContent(e0+1,e1+1,-AbsoluteValue_h2o);
		}		
	      }
	    }

	    else if(ErrorType>=Systematics_Xsec_Start && ErrorType<=Systematics_Xsec_End){
	      double RelativeMC=MCReconstructedEvents[e0][e1];

	      //if(NominalSelectedMC_XSTemp->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC_XSTemp->GetBinContent(e0+1,e1+1);
	      if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeMC/=NominalSelectedMC->GetBinContent(e0+1,e1+1);
 
	      ErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]->Fill(XsecVariation,RelativeMC);	    
	      ErrorXS_Norm[ErrorType-Systematics_Xsec_Start][e0][e1]->Fill(XsecVariation);
	      //cout<<XsecVariation<<endl;
	      xXS[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=XsecVariation;
	      yXS[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=RelativeMC;
	      // cout<<ErrorType<<" "<<e0<<" "<<e1<<" "<<XsecVariation<<" "<<RelativeMC<<endl;
	      xXS_reduced[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations-2))]=XsecVariation;
	      yXS_reduced[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations-2))]=RelativeMC;
	      //   if(abs(XsecVariation)<=1)   cout<<XsecVariation<<" "<<e0<<" "<<e1<<" "<<MCReconstructedEvents[e0][e1]<<" "<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<" "<<RelativeMC<<endl;

	      if(ratio){
		double RelativeMC_ch=MCReconstructedEvents_ch[e0][e1];
		if(NominalSelectedMC_ch->GetBinContent(e0+1,e1+1)!=0) RelativeMC_ch/=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
		xXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=XsecVariation;
		yXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=RelativeMC_ch;
		xXS_reduced_ch[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations-2))]=XsecVariation;
		yXS_reduced_ch[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations-2))]=RelativeMC_ch;
		
		double RelativeMC_h2o=MCReconstructedEvents_h2o[e0][e1];
		if(NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1)!=0) RelativeMC_h2o/=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
		xXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=XsecVariation;
		yXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations))]=RelativeMC_h2o;
		xXS_reduced_h2o[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations-2))]=XsecVariation;
		yXS_reduced_h2o[ErrorType-Systematics_Xsec_Start][e0][e1][((int) (XsecVariation+CenterXsecVariations-2))]=RelativeMC_h2o;
	      }
	    }
	  }	  
	}
	
	if(ratio) TotalValue=TotalValue_h2o/TotalValue_ch;

#ifdef RECONSTRUCTED
	if(ErrorType==Systematics_Flux_Start){
	  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over effect 0
	    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over effect 1
	      ErrorEffFlux[c0][c1]->Fill(Efficiency[c0][c1]*NumberOfPOT/NominalEfficiency[c0][c1]-1.);
	    }
	  }
	}
#endif
      }
    
      cout<<"HERE="<<ErrorType<<endl;
    
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
      else if(ErrorType==2){
#ifdef RECONSTRUCTED
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    NoiseEventFunctions->cd();
	    // cout<<e0<<" "<<e1<<" "<<NominalSelectedMC->GetBinContent(e0+1,e1+1)<<endl;
	    ErrorNoise[e0][e1]->Divide(ErrorNoise_Norm[e0][e1]);
	    ErrorNoise[e0][e1]->Fit(Form("fErrorNoise[%d][%d]",e0,e1),"R");
	  
	    ErrorNoise[e0][e1]->Write();

	    double RelativeValue=fErrorNoise[e0][e1]->Eval(0)-NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) RelativeValue/=NominalSelectedMC->GetBinContent(e0+1,e1+1);

	    Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(RelativeValue));
	    Error_Plus[ErrorType]->SetBinContent(e0+1,e1+1,TMath::Abs(RelativeValue));


	    // ML WARNING: at this point MCReconstructedEvents contains the data of Systematics2_9.txt ...
	    // 2017/07/20 changed to NominalSelectedMC (cf treatment of Error>=7)
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		double RelativeValue2=fErrorNoise[f0][f1]->Eval(0)-NominalSelectedMC->GetBinContent(f0+1,f1+1);
		if(NominalSelectedMC->GetBinContent(f0+1,f1+1)!=0) RelativeValue2/=NominalSelectedMC->GetBinContent(f0+1,f1+1);
		//Covariance[ErrorType][ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]*(1+RelativeValue)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]*(1+RelativeValue2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		//CovarianceReduced[ErrorType][e0][e1][f0][f1]+=(MCReconstructedEvents[e0][e1]*(1+RelativeValue)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(MCReconstructedEvents[f0][f1]*(1+RelativeValue2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
		CovarianceReduced[ErrorType][e0][e1][f0][f1]=(NominalSelectedMC->GetBinContent(e0+1,e1+1)*(1+RelativeValue)-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NominalSelectedMC->GetBinContent(f0+1,f1+1)*(1+RelativeValue2)-NominalSelectedMC->GetBinContent(f0+1,f1+1));
	      }
	    }	    
	  }
	}
#else
	// case of unfolded events : used as free spot for Koga's det errors
	for(int e0=0;e0<NBinsMom;e0++){
	  for(int e1=0;e1<NBinsAngle;e1++){
	    double Value1=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    double Error1=fmax(Error_Plus[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus[ErrorType]->GetBinContent(e0+1,e1+1));
	    Error1*=Value1;

	    double Error1_ch=fmax(Error_Plus_ch[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus_ch[ErrorType]->GetBinContent(e0+1,e1+1));	
	    double Error1_h2o=fmax(Error_Plus_h2o[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus_h2o[ErrorType]->GetBinContent(e0+1,e1+1));

	    for(int f0=0;f0<NBinsMom;f0++){
	      for(int f1=0;f1<NBinsAngle;f1++){
		double Value2=NominalSelectedMC->GetBinContent(f0+1,f1+1);
		double Error2=fmax(Error_Plus[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus[ErrorType]->GetBinContent(f0+1,f1+1));
		Error2*=Value2;

		CovarianceReduced[ErrorType][e0][e1][f0][f1]=Error1*Error2;

		if(ratio){
		  double Error2_ch=fmax(Error_Plus_ch[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus_ch[ErrorType]->GetBinContent(f0+1,f1+1));	
		  double Error2_h2o=fmax(Error_Plus_h2o[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus_h2o[ErrorType]->GetBinContent(f0+1,f1+1));
		  
		  CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][0]=Error1_h2o*Error2_h2o;
		  CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][1]=Error1_ch*Error2_ch;
		  // they are uncorrelated
		}
	      }
	    }
	  }
	}
#endif
      }
      /*else if(ErrorType==3){
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	HitEfficiencyFunctions->cd();
	cout<<e0<<" "<<e1<<endl;
	ErrorHitEfficiency[e0][e1]->Fit(Form("fErrorHitEfficiency[%d][%d]",e0,e1),"RQ");
	ErrorHitEfficiency[e0][e1]->Write();
	cout<<"CAREFUL: check function before estimation systematic error"<<endl;
	  
	//Error_Minus[ErrorType]->SetBinContent(e0+1,e1+1,-TMath::Abs(fErrorHitEfficiency[e0][e1]
	}
	}
	}*/
      else if(ErrorType==5){
	for(int e0=0;e0<NBinsMom;e0++){
	  for(int e1=0;e1<NBinsAngle;e1++){
	    double Value1=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    double Error1=fmax(Error_Plus[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus[ErrorType]->GetBinContent(e0+1,e1+1));
	    Error1*=Value1;
	    for(int f0=0;f0<NBinsMom;f0++){
	      for(int f1=0;f1<NBinsAngle;f1++){
		double Value2=NominalSelectedMC->GetBinContent(f0+1,f1+1);
		double Error2=fmax(Error_Plus[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus[ErrorType]->GetBinContent(f0+1,f1+1));
		Error2*=Value2;
		//   CovarianceReduced[ErrorType][e0][e1][f0][f1]=Error1*Error2;
	      } 
	    }
	  }
	}
      }
      else if(ErrorType==6){// ML added 2017/07/21
	for(int e0=0;e0<NBinsMom;e0++){
	  for(int e1=0;e1<NBinsAngle;e1++){
	    double Value1=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    double Error1=fmax(Error_Plus[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus[ErrorType]->GetBinContent(e0+1,e1+1));
	    Error1*=Value1;
	    double Error1_ch, Error1_h2o;
	    if(ratio){
	      Error1_ch=fmax(Error_Plus_ch[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus_ch[ErrorType]->GetBinContent(e0+1,e1+1));
	      Error1_h2o=fmax(Error_Plus_h2o[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus_h2o[ErrorType]->GetBinContent(e0+1,e1+1));
	    }
	    for(int f0=0;f0<NBinsMom;f0++){
	      for(int f1=0;f1<NBinsAngle;f1++){
		double Value2=NominalSelectedMC->GetBinContent(f0+1,f1+1);
		double Error2=fmax(Error_Plus[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus[ErrorType]->GetBinContent(f0+1,f1+1));
		Error2*=Value2;
		CovarianceReduced[ErrorType][e0][e1][f0][f1]=Error1*Error2;
		
		if(ratio){
		  double Error2_ch=fmax(Error_Plus_ch[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus_ch[ErrorType]->GetBinContent(f0+1,f1+1));
		  double Error2_h2o=fmax(Error_Plus_h2o[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus_h2o[ErrorType]->GetBinContent(f0+1,f1+1));
		  CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][0]+=Error1_h2o*Error2_h2o;
		  CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][1]+=Error1_h2o*Error2_ch;
		  CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][0]+=Error1_ch*Error2_h2o;
		  CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][1]+=Error1_ch*Error2_ch;
		}	      

	      }
	    }
	  }
	}
      }
      else if(ErrorType>=7 && ErrorType<=Systematics_Detector_End){
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  
	    double ePlus=0;
	    double eMinus=0;
	    for(int ibin=1;ibin<=MCContent[ErrorType][e0][e1]->GetNbinsX();ibin++){ // the x axis is ErrorValue (e.g. nb of veto planes...)
	      double Error=(DataContent[ErrorType][e0][e1]->GetBinContent(ibin)-MCContent[ErrorType][e0][e1]->GetBinContent(ibin));
	      if( Error > ePlus ) ePlus=Error;
	      else if( (-Error) > (-eMinus) ) eMinus=Error; //ML 2017/07/20 add the - sign
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
      if(ErrorType>=Systematics_Flux_Start && ErrorType<=Systematics_Flux_End){
	cout<<"create the spline... "<<ErrorType<<endl;
	if(ErrorType==16){
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
#ifdef RECONSTRUCTED
	  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over effect 0
	    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over effect 1
	      FluxEventFunctions->cd();
	      ErrorEffFlux[c0][c1]->Fit(Form("fErrorEffFlux[%d][%d]",c0,c1),"RQ");
	      ErrorEffFlux[c0][c1]->Write();
	    }
	  }
#endif
	}
	else if(ErrorType==17){
	  for(int e0=0;e0<NBinsMom;e0++){
	    for(int e1=0;e1<NBinsAngle;e1++){
	      double Value1=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	      double Error1=fmax(Error_Plus[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus[ErrorType]->GetBinContent(e0+1,e1+1));
	      Error1*=Value1;
	      double Error1_ch, Error1_h2o;
	      if(ratio){
		Error1_ch=fmax(Error_Plus_ch[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus_ch[ErrorType]->GetBinContent(e0+1,e1+1));
		Error1_h2o=fmax(Error_Plus_h2o[ErrorType]->GetBinContent(e0+1,e1+1),Error_Minus_h2o[ErrorType]->GetBinContent(e0+1,e1+1));
	      }
	      for(int f0=0;f0<NBinsMom;f0++){
		for(int f1=0;f1<NBinsAngle;f1++){
		  double Value2=NominalSelectedMC->GetBinContent(f0+1,f1+1);
		  double Error2=fmax(Error_Plus[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus[ErrorType]->GetBinContent(f0+1,f1+1));
		  Error2*=Value2;
		  CovarianceFlux[e0][e1][f0][f1]+=Error1*Error2;

		  if(ratio){
		    double Error2_ch=fmax(Error_Plus_ch[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus_ch[ErrorType]->GetBinContent(f0+1,f1+1));
		    double Error2_h2o=fmax(Error_Plus_h2o[ErrorType]->GetBinContent(f0+1,f1+1),Error_Minus_h2o[ErrorType]->GetBinContent(f0+1,f1+1));
		    CovarianceFluxTargets[e0][e1][f0][f1][0][0]+=Error1_h2o*Error2_h2o;
		    CovarianceFluxTargets[e0][e1][f0][f1][0][1]+=Error1_h2o*Error2_ch;
		    CovarianceFluxTargets[e0][e1][f0][f1][1][0]+=Error1_ch*Error2_h2o;
		    CovarianceFluxTargets[e0][e1][f0][f1][1][1]+=Error1_ch*Error2_ch;
		  }	      


		}
	      }
	    }
	  }
	}
      }
      else if(ErrorType>=Systematics_Xsec_Start/* && ErrorType<=Systematics_Xsec_End*/){      
	/*    int a=18;
	      if(a>=18 && a<=37){      */
	cout<<"create the spline... "<<ErrorType<<endl;
	//create the spline
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    XSEventFunctions->cd();
	    if(IsReducedXSSyst(ErrorType-Systematics_Xsec_Start)){
	      gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations-4,xXS_reduced[ErrorType-Systematics_Xsec_Start][e0][e1],yXS_reduced[ErrorType-Systematics_Xsec_Start][e0][e1]);
	      if(ratio){
		gErrorXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations-4,xXS_reduced_ch[ErrorType-Systematics_Xsec_Start][e0][e1],yXS_reduced_ch[ErrorType-Systematics_Xsec_Start][e0][e1]);
		gErrorXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations-4,xXS_reduced_h2o[ErrorType-Systematics_Xsec_Start][e0][e1],yXS_reduced_h2o[ErrorType-Systematics_Xsec_Start][e0][e1]);
	      }
	    }
	    else {
	      gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations,xXS[ErrorType-Systematics_Xsec_Start][e0][e1],yXS[ErrorType-Systematics_Xsec_Start][e0][e1]);
	      if(ratio){
		gErrorXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations,xXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1],yXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1]);
		gErrorXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1] = new TGraph(NXsecVariations,xXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1],yXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1]);
	      }
	      
	    }
	    gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("gErrorXS_%d_%d_%d",ErrorType-Systematics_Xsec_Start,e0,e1));
	    sErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1] = new TSpline3(Form("sErrorXS[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),gErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]);
	    sErrorXS[ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("sErrorXS_%d_%d_%d",ErrorType-Systematics_Xsec_Start,e0,e1)); 

	    if(ratio){
	      gErrorXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("gErrorXS_ch_%d_%d_%d",ErrorType-Systematics_Xsec_Start,e0,e1));
	      gErrorXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1]->Write(Form("gErrorXS_h2o_%d_%d_%d",ErrorType-Systematics_Xsec_Start,e0,e1));
	      sErrorXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1] = new TSpline3(Form("sErrorXS_ch[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),gErrorXS_ch[ErrorType-Systematics_Xsec_Start][e0][e1]);
	      sErrorXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1] = new TSpline3(Form("sErrorXS_h2o[%d][%d][%d]",ErrorType-Systematics_Xsec_Start,e0,e1),gErrorXS_h2o[ErrorType-Systematics_Xsec_Start][e0][e1]);
	    }
	  }
	}
      }
      ErrorType=ErrorType;
      //    cout<<ErrorType<<" xsec?"<<(ErrorType>17)<<" "<<(ErrorType<=37)<<" "<<ErrorType<<" "<<((ErrorType>17) && (ErrorType<=37))<<endl;

      // covariance reduced not filled here for Flux or XS errors    
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
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
    //--- although here no correlations is assumed in the computation of the total CovarianceXS ... -----
    double TotalRelxsError[EndXsec+1]={0};
    bool CorrelateMEC=true;
    double MEC_C_O_corr=1;
    if(EndError>=Systematics_Xsec_Start && !FakeData){
      //cout<<"hello"<<endl;
      //Starts the toy experiments:
 
      double Error;
      double NEventsXS1,NEventsXS2;
      double NEventsXS1b,NEventsXS2b;
      double NEventsXS1_ch,NEventsXS2_ch,NEventsXS1_h2o,NEventsXS2_h2o;
      double NEventsXS1b_ch,NEventsXS2b_ch,NEventsXS1b_h2o,NEventsXS2b_h2o;

      int NToysXsec=5000;
      TRandom3 * rxs = new TRandom3();
      rxs->SetSeed(0);

      for(int s1=0;s1<=EndError-Systematics_Xsec_Start;s1++){
	//for(int s1=0;s1<=0;s1++){//TEMP

	if(CorrelateMEC && s1==5) continue;
	cout<<"source tested="<<s1<<endl;
      
	for(int nt=0;nt<NToysXsec;nt++){
	  if(nt%1000==0) cout<<"toy #"<<nt<<endl;
	  Error=10;
	
	  //My interpolation doesn't go further away than -+3sigma. Therefore, only keep toy experiment inside these values.
	  //while((IsReducedXSSyst(s1)&&TMath::Abs(Error)>1) || (!IsReducedXSSyst(s1)&&TMath::Abs(Error)>3)) Error=rxs->Gaus(0,1);
	  while(TMath::Abs(Error)>3) Error=rxs->Gaus(0,1);
	  if(IsReducedXSSyst(s1)){
	    if(Error<-1) Error=-1.;
	    if(Error>1) Error=1.;
	  }
	  if(s1==4) h_test_red->Fill(Error);
	  if(s1==6) h_test->Fill(Error);

	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	      //	    cout<<e0<<" "<<e1<<" "<<bin1<<" "<<Error<<endl;
	      NEventsXS1=NominalSelectedMC->GetBinContent(e0+1,e1+1)*sErrorXS[s1][e0][e1]->Eval(Error);
	      if(CorrelateMEC && s1==4) NEventsXS1b=NominalSelectedMC->GetBinContent(e0+1,e1+1)*sErrorXS[s1+1][e0][e1]->Eval(Error);
	      if(ratio){
 		NEventsXS1_ch=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1)*(sErrorXS_ch[s1][e0][e1]->Eval(Error)-1);
		NEventsXS1_h2o=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1)*(sErrorXS_h2o[s1][e0][e1]->Eval(Error)-1);

		if(CorrelateMEC && s1==4){
		  NEventsXS1b_ch=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1)*(sErrorXS_ch[s1+1][e0][e1]->Eval(Error)-1);
		  NEventsXS1b_h2o=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1)*(sErrorXS_h2o[s1+1][e0][e1]->Eval(Error)-1);
		}
	      }
	      
	      //if(e0==4 && e1==4) cout<<"Error="<<Error<<", "<<NEventsXS[bin1]-NominalSelectedMC->GetBinContent(e0+1,e1+1)<<", value="<<sErrorXS[2][4][4]->Eval(Error)<<endl;
	      for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
		for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		  NEventsXS2=NominalSelectedMC->GetBinContent(f0+1,f1+1)*sErrorXS[s1][f0][f1]->Eval(Error);
		  if(CorrelateMEC && s1==4) NEventsXS2b=NominalSelectedMC->GetBinContent(f0+1,f1+1)*sErrorXS[s1+1][f0][f1]->Eval(Error);
		  //cout<<bin1<<" "<<NEventsXS[bin1]<<" "<<sErrorXS[s1][e0][e1]->Eval(Error)<<" "<<bin2<<" "<<NEventsXS[bin2]<<" "<<sErrorXS[s1][f0][f1]->Eval(Error)<<endl;
		  CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS1-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		  if(CorrelateMEC && s1==4){
		    CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS1b-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2b-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		    CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=MEC_C_O_corr*(1./(NToysXsec-1.))*(NEventsXS1-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2b-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		    CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1]+=MEC_C_O_corr*(1./(NToysXsec-1.))*(NEventsXS1b-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		  }
		  // assuming no correlations between error sources:
		  CovarianceXS[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS1-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		  if(CorrelateMEC && s1==4){
		    CovarianceXS[e0][e1][f0][f1]+=(1./(NToysXsec-1.))*(NEventsXS1b-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2b-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		    CovarianceXS[e0][e1][f0][f1]+=MEC_C_O_corr*(1./(NToysXsec-1.))*(NEventsXS1-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2b-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		    CovarianceXS[e0][e1][f0][f1]+=MEC_C_O_corr*(1./(NToysXsec-1.))*(NEventsXS1b-NominalSelectedMC->GetBinContent(e0+1,e1+1))*(NEventsXS2-NominalSelectedMC->GetBinContent(f0+1,f1+1)); 
		  }
		  
		  if(ratio){
		    NEventsXS2_ch=NominalSelectedMC_ch->GetBinContent(f0+1,f1+1)*(sErrorXS_ch[s1][f0][f1]->Eval(Error)-1);
		    NEventsXS2_h2o=NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1)*(sErrorXS_h2o[s1][f0][f1]->Eval(Error)-1);

		    CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][0]+=NEventsXS1_h2o*NEventsXS2_h2o/((float)NToysXsec-1.);
		    CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][1]+=NEventsXS1_h2o*NEventsXS2_ch/((float)NToysXsec-1.);
		    CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][0]+=NEventsXS1_ch*NEventsXS2_h2o/((float)NToysXsec-1.);
		    CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][1]+=NEventsXS1_ch*NEventsXS2_ch/((float)NToysXsec-1.);

		    CovarianceXSTargets[e0][e1][f0][f1][0][0]+=NEventsXS1_h2o*NEventsXS2_h2o/((float)NToysXsec-1.);
		    CovarianceXSTargets[e0][e1][f0][f1][0][1]+=NEventsXS1_h2o*NEventsXS2_ch/((float)NToysXsec-1.);
		    CovarianceXSTargets[e0][e1][f0][f1][1][0]+=NEventsXS1_ch*NEventsXS2_h2o/((float)NToysXsec-1.);
		    CovarianceXSTargets[e0][e1][f0][f1][1][1]+=NEventsXS1_ch*NEventsXS2_ch/((float)NToysXsec-1.);

		    if(CorrelateMEC && s1==4){
		      NEventsXS2b_ch=NominalSelectedMC_ch->GetBinContent(f0+1,f1+1)*(sErrorXS_ch[s1+1][f0][f1]->Eval(Error)-1);
		      NEventsXS2b_h2o=NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1)*(sErrorXS_h2o[s1+1][f0][f1]->Eval(Error)-1);

		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][0]+=NEventsXS1b_h2o*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][1]+=NEventsXS1b_h2o*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][0]+=NEventsXS1b_ch*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][1]+=NEventsXS1b_ch*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][0]+=MEC_C_O_corr*NEventsXS1_h2o*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][1]+=MEC_C_O_corr*NEventsXS1_h2o*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][0]+=MEC_C_O_corr*NEventsXS1_ch*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][1]+=MEC_C_O_corr*NEventsXS1_ch*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][0]+=MEC_C_O_corr*NEventsXS1b_h2o*NEventsXS2_h2o/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][1]+=MEC_C_O_corr*NEventsXS1b_h2o*NEventsXS2_ch/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][0]+=MEC_C_O_corr*NEventsXS1b_ch*NEventsXS2_h2o/((float)NToysXsec-1.);
		      CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][1]+=MEC_C_O_corr*NEventsXS1b_ch*NEventsXS2_ch/((float)NToysXsec-1.);

		      CovarianceXSTargets[e0][e1][f0][f1][0][0]+=NEventsXS1b_h2o*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][0][1]+=NEventsXS1b_h2o*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][1][0]+=NEventsXS1b_ch*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][1][1]+=NEventsXS1b_ch*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][0][0]+=MEC_C_O_corr*NEventsXS1_h2o*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][0][1]+=MEC_C_O_corr*NEventsXS1_h2o*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][1][0]+=MEC_C_O_corr*NEventsXS1_ch*NEventsXS2b_h2o/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][1][1]+=MEC_C_O_corr*NEventsXS1_ch*NEventsXS2b_ch/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][0][0]+=MEC_C_O_corr*NEventsXS1b_h2o*NEventsXS2_h2o/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][0][1]+=MEC_C_O_corr*NEventsXS1b_h2o*NEventsXS2_ch/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][1][0]+=MEC_C_O_corr*NEventsXS1b_ch*NEventsXS2_h2o/((float)NToysXsec-1.);
		      CovarianceXSTargets[e0][e1][f0][f1][1][1]+=MEC_C_O_corr*NEventsXS1b_ch*NEventsXS2_ch/((float)NToysXsec-1.);
		    }	
		  }		  
		}
	      }	
	    }
	  }
	}

	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    double err=TMath::Sqrt(CovarianceReduced[s1][e0][e1][e0][e1]);
	    if(NominalSelectedMC->GetBinContent(e0+1,e1+1)!=0) err/=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    Error_Minus[Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,-err);
	    Error_Plus[Systematics_Xsec_Start+s1]->SetBinContent(e0+1,e1+1,err);
	  
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	      if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	      for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
		if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
		if(ratio){
		  TotalRelxsError[s1] +=CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
		  TotalRelxsError[s1] +=CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
		  TotalRelxsError[s1] += -2*CovarianceReducedTargets[Systematics_Xsec_Start+s1][e0][e1][f0][f1][1][0]/(TotalValue_ch*TotalValue_h2o);
		}
		else TotalRelxsError[s1]+=CovarianceReduced[Systematics_Xsec_Start+s1][e0][e1][f0][f1];

		/*CorrelationXS[e0][e1][f0][f1]=CovarianceXS[e0][e1][e0][e1];		
		double Diagonal=CovarianceXS[e0][e1][e0][e1]*CovarianceXS[f0][f1][f0][f1];
		if(Diagonal!=0) CorrelationXS[e0][e1][f0][f1]/=pow(Diagonal,1/2);*/
	      }
	    }
	  }
	}
	if(ratio) TotalRelxsError[s1]*=pow(TotalValue,2.);
	file->cd();
	Error_Minus[Systematics_Xsec_Start+s1]->Write();
	Error_Plus[Systematics_Xsec_Start+s1]->Write();
	
	Evaluate1DError(NominalSelectedMC, CovarianceReduced[Systematics_Xsec_Start+s1], Error_RecMom[Systematics_Xsec_Start+s1], Error_RecAngle[Systematics_Xsec_Start+s1]);
	Error_RecMom[Systematics_Xsec_Start+s1]->Write();
	Error_RecAngle[Systematics_Xsec_Start+s1]->Write();
	
	if(TotalValue!=0) TotalRelxsError[s1]=TMath::Sqrt(TotalRelxsError[s1])/TotalValue;
    
#ifdef XSTABLE
    
	cout<<"Total error so far, when adding="<<Systematics_Xsec_Start+s1<<endl;
	for(int e0=0;e0<NBinsMom;e0++){
	  for(int e1=0;e1<NBinsAngle;e1++){
	    double ErrorSquared=CovarianceXS[e0][e1][e0][e1];
	    double Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	    cout<<"Error XS="<<setprecision(3)<<std::scientific<<ErrorSquared<<endl;
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
    TBox * boxsliceErrorDet_RecMom[NBinsMom][NBinsAngle];
    TBox * boxsliceErrorFlux_RecMom[NBinsMom][NBinsAngle];
    TBox * boxsliceErrorXS_RecMom[NBinsMom][NBinsAngle];

    TBox * boxsliceErrorStat_RecAngle[NBinsMom][NBinsAngle];
    TBox * boxsliceErrorDet_RecAngle[NBinsMom][NBinsAngle];
    TBox * boxsliceErrorFlux_RecAngle[NBinsMom][NBinsAngle];
    TBox * boxsliceErrorXS_RecAngle[NBinsMom][NBinsAngle];

    TBox * boxErrorStat_RecMom[NBinsMom];
    TBox * boxErrorDet_RecMom[NBinsMom];
    TBox * boxErrorFlux_RecMom[NBinsMom];
    TBox * boxErrorXS_RecMom[NBinsMom];

    TBox * boxErrorStat_RecAngle[NBinsAngle];
    TBox * boxErrorDet_RecAngle[NBinsAngle];
    TBox * boxErrorFlux_RecAngle[NBinsAngle];
    TBox * boxErrorXS_RecAngle[NBinsAngle];

  
    //###################################CREATE 1D SLICE DISTRIBUTION##################"
    TH1D * sliceNominalSelectedMC_RecMom[NBinsAngle];
    TH1D * sliceNominalSelectedMC_RecAngle[NBinsMom];
    // --------  ML 2017/12/14 - because of trash bins the projection is not that simple; done just after
    //  for(int e1=0;e1<NBinsAngle;e1++) sliceNominalSelectedMC_RecMom[e1] = (TH1D*) NominalSelectedMC->ProjectionX(Form("sliceNominalSelectedMC_RecMom[%d]",e1),e1+1,e1+1);
    //  for(int e0=0;e0<NBinsMom;e0++) sliceNominalSelectedMC_RecAngle[e0] = (TH1D*) NominalSelectedMC->ProjectionX(Form("sliceNominalSelectedMC_RecAngle[%d]",e0),e0+1,e0+1);
    for(int e1=0;e1<NBinsAngle;e1++) sliceNominalSelectedMC_RecMom[e1] = new TH1D(Form("sliceNominalSelectedMC_RecMom[%d]",e1),Form("%d #circ < #theta_{#mu} < %d #circ",(int)BinningAngle[e1],(int)BinningAngle[e1+1]),NBinsMom,BinningMom);
    for(int e0=0;e0<NBinsMom;e0++) sliceNominalSelectedMC_RecAngle[e0] = new TH1D(Form("sliceNominalSelectedMC_RecAngle[%d]",e0),Form("%f GeV < p_{#mu} < %f GeV",BinningMom[e0],BinningMom[e0+1]),NBinsAngle,BinningAngle);
    //###################################


    //############################ PRODUCE HISTOS FOR TRUE MC ######################
    TH1D * sliceNominalTrueMC_RecMom[NBinsAngle];
    TH1D * sliceNominalTrueMC_RecAngle[NBinsMom];
    for(int e1=0;e1<NBinsAngle;e1++) sliceNominalTrueMC_RecMom[e1] = new TH1D(Form("sliceNominalTrueMC_RecMom[%d]",e1),Form("%d #circ < #theta_{#mu} < %d #circ",(int)BinningAngle[e1],(int)BinningAngle[e1+1]),NBinsMom,BinningMom);
    for(int e0=0;e0<NBinsMom;e0++) sliceNominalTrueMC_RecAngle[e0] = new TH1D(Form("sliceNominalTrueMC_RecAngle[%d]",e0),Form("%f GeV < p_{#mu} < %f GeV",BinningMom[e0],BinningMom[e0+1]),NBinsAngle,BinningAngle);
    TH1D * NominalTrueMC_RecMom = new TH1D("NominalTrueMC_RecMom","p_{#mu} distribution of selected events",NBinsMom,BinningMom);
    TH1D * NominalTrueMC_RecAngle = new TH1D("NominalTrueMC_RecAngle","#theta_{#mu} distribution of selected events",NBinsAngle,BinningAngle);
  
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	sliceNominalTrueMC_RecMom[e1]->SetBinContent(e0+1,NominalTrueMC->GetBinContent(e0+1,e1+1));
	sliceNominalTrueMC_RecAngle[e0]->SetBinContent(e1+1,NominalTrueMC->GetBinContent(e0+1,e1+1));
      }
    }

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;

      double TrueValue=0, TrueValue_ch=0., TrueValue_h2o=0.;
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	TrueValue+=NominalTrueMC->GetBinContent(e0+1,e1+1);
	if(ratio){
	  TrueValue_ch+=NominalTrueMC_ch->GetBinContent(e0+1,e1+1);
	  TrueValue_h2o+=NominalTrueMC_h2o->GetBinContent(e0+1,e1+1);
	}
      }
      if(ratio) NominalTrueMC_RecMom->SetBinContent(e0+1,TrueValue);
      else NominalTrueMC_RecMom->SetBinContent(e0+1,TrueValue_h2o/TrueValue_ch);
    }

    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
      if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      double TrueValue=0., TrueValue_ch=0., TrueValue_h2o=0.;
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	TrueValue+=NominalTrueMC->GetBinContent(e0+1,e1+1);
	if(ratio){
	  TrueValue_ch+=NominalTrueMC_ch->GetBinContent(e0+1,e1+1);
	  TrueValue_h2o+=NominalTrueMC_h2o->GetBinContent(e0+1,e1+1);
	}
      }
      if(ratio) NominalTrueMC_RecAngle->SetBinContent(e1+1,TrueValue);
      else NominalTrueMC_RecAngle->SetBinContent(e1+1,TrueValue_h2o/TrueValue_ch);
    }
    //###########################################################################################    


    //####################################### TOTAL COVARIANCE ################################
    cout<<"\nTotal error from total covariance matrix:"<<endl;

    int NBins=NBinsAngle*NBinsMom;
    if(UseTrashBins) NBins=(NBinsAngle-1)*(NBinsMom-2);

    double data[NBins*NBins],datacorr[NBins*NBins];

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    /*	    CovarianceTotal[e0][e1][f0][f1]=CovarianceFlux[e0][e1][f0][f1]+CovarianceXS[e0][e1][f0][f1];
#ifdef DETECTORSYST
	    for(int ErrorType=Systematics_Detector_Start;ErrorType<=Systematics_Detector_End;ErrorType++){
	      if(!useKogaDetError && !IsReadyDetSyst(ErrorType)) continue;
	      if(ErrorType==5){//tmp I use only Birks +1sigma
		CovarianceTotal[e0][e1][f0][f1]+=CovarianceBirksPlus[e0][e1][f0][f1];
	      }
	      else CovarianceTotal[e0][e1][f0][f1]+=CovarianceReduced[ErrorType][e0][e1][f0][f1];
	    }
#endif*/
#ifndef RECONSTRUCTED
	    CovarianceTotal[e0][e1][f0][f1]+=CovarianceStatistics[e0][e1][f0][f1];
#endif

	    int BinX=(e0-1)*(NBinsAngle-1)+e1;
	    int BinY=(f0-1)*(NBinsAngle-1)+f1;
	    if(!UseTrashBins) {BinX=e0*NBinsAngle+e1; BinY=f0*NBinsAngle+f1;}
	    else if(IsTrueMomTrashBin[e0] || IsTrueAngleTrashBin[e1] || IsTrueMomTrashBin[f0] || IsTrueAngleTrashBin[f1]) continue;

	    data[BinX+BinY*NBins]=CovarianceTotal[e0][e1][f0][f1];


	  }
	}
	cout<<"("<<e0<<","<<e1<<")=\t"<<TMath::Sqrt(CovarianceTotal[e0][e1][e0][e1])/NominalSelectedMC->GetBinContent(e0+1,e1+1)<<endl;
      }
    }

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    int BinX=(e0-1)*(NBinsAngle-1)+e1;
	    int BinY=(f0-1)*(NBinsAngle-1)+f1;
	    if(!UseTrashBins) {BinX=e0*NBinsAngle+e1; BinY=f0*NBinsAngle+f1;}
	    else if(IsTrueMomTrashBin[e0] || IsTrueAngleTrashBin[e1] || IsTrueMomTrashBin[f0] || IsTrueAngleTrashBin[f1]) continue;
	    datacorr[BinX+BinY*NBins]=CovarianceTotal[e0][e1][f0][f1]/TMath::Sqrt(CovarianceTotal[e0][e1][e0][e1]*CovarianceTotal[f0][f1][f0][f1]);
	  }
	}
      }
    }
	    
    TMatrixTSym<double> CovMatrix(NBins,data);
    TMatrixTSym<double> InvCovMatrix(CovMatrix); 
    InvCovMatrix.Invert(); 
    TMatrixTSym<double> CorrMatrix(NBins,datacorr);
    TMatrixTSym<double> InvCorrMatrix(CorrMatrix); 
    InvCorrMatrix.Invert(); 
    CovMatrix.Print();
    InvCovMatrix.Print();
    /*TMatrixDEigen CorrEigen(CorrMatrix);
    CorrEigen.GetEigenValues().Print();*/
    //###########################################################################################    


    //######################### FRACTIONNAL COVARIANCE MATRICES (2 TARGETS) ##############################
    if(ratio){
      TH2D* FracCovFlux=new TH2D("FracCovFlux","Fractionnal covariance matrix",NBins*2,0,NBins*2,NBins*2,0,NBins*2);
      TH2D* FracCovXS=new TH2D("FracCovXS","Fractionnal covariance matrix",NBins*2,0,NBins*2,NBins*2,0,NBins*2);
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      int BinX=(e0-1)*(NBinsAngle-1)+e1+1;
	      int BinY=(f0-1)*(NBinsAngle-1)+f1+1;
	      if(!UseTrashBins) {BinX=e0*NBinsAngle+e1+1; BinY=f0*NBinsAngle+f1+1;}
	      else if(IsTrueMomTrashBin[e0] || IsTrueAngleTrashBin[e1] || IsTrueMomTrashBin[f0] || IsTrueAngleTrashBin[f1]) continue;
	    
	      double valueX_ch=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	      double valueX_h2o=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	      double valueY_ch=NominalSelectedMC_ch->GetBinContent(f0+1,f1+1);
	      double valueY_h2o=NominalSelectedMC_h2o->GetBinContent(f0+1,f1+1);
	      /*
	      if(CovarianceXSTargets[e0][e1][f0][f1][0][0]<0) cout<<"**Negative cov value ww**"<<endl;
	      if(CovarianceXSTargets[e0][e1][f0][f1][0][1]<0) cout<<"**Negative cov value wc**"<<endl;
	      if(CovarianceXSTargets[e0][e1][f0][f1][1][0]<0) cout<<"**Negative cov value cw**"<<endl;
	      if(CovarianceXSTargets[e0][e1][f0][f1][1][1]<0) cout<<"**Negative cov value cc**"<<endl;*/

	      FracCovFlux->SetBinContent(BinX,BinY,CovarianceFluxTargets[e0][e1][f0][f1][0][0]/sqrt(fabs(CovarianceFluxTargets[e0][e1][f0][f1][0][0]))/sqrt(valueX_h2o*valueY_h2o));
	      FracCovFlux->SetBinContent(BinX,BinY+NBins,CovarianceFluxTargets[e0][e1][f0][f1][0][1]/sqrt(fabs(CovarianceFluxTargets[e0][e1][f0][f1][0][1]))/sqrt(valueX_h2o*valueY_ch));
	      FracCovFlux->SetBinContent(BinX+NBins,BinY,CovarianceFluxTargets[e0][e1][f0][f1][1][0]/sqrt(fabs(CovarianceFluxTargets[e0][e1][f0][f1][1][0]))/sqrt(valueX_ch*valueY_h2o));
	      FracCovFlux->SetBinContent(BinX+NBins,BinY+NBins,CovarianceFluxTargets[e0][e1][f0][f1][1][1]/sqrt(fabs(CovarianceFluxTargets[e0][e1][f0][f1][1][1]))/sqrt(valueX_ch*valueY_ch));

	      FracCovXS->SetBinContent(BinX,BinY,CovarianceXSTargets[e0][e1][f0][f1][0][0]/sqrt(fabs(CovarianceXSTargets[e0][e1][f0][f1][0][0]))/sqrt(valueX_h2o*valueY_h2o));
	      FracCovXS->SetBinContent(BinX,BinY+NBins,CovarianceXSTargets[e0][e1][f0][f1][0][1]/sqrt(fabs(CovarianceXSTargets[e0][e1][f0][f1][0][1]))/sqrt(valueX_h2o*valueY_ch));
	      FracCovXS->SetBinContent(BinX+NBins,BinY,CovarianceXSTargets[e0][e1][f0][f1][1][0]/sqrt(fabs(CovarianceXSTargets[e0][e1][f0][f1][1][0]))/sqrt(valueX_ch*valueY_h2o));
	      FracCovXS->SetBinContent(BinX+NBins,BinY+NBins,CovarianceXSTargets[e0][e1][f0][f1][1][1]/sqrt(fabs(CovarianceXSTargets[e0][e1][f0][f1][1][1]))/sqrt(valueX_ch*valueY_ch));
	    }
	  }
	}
      }
      file->cd();
      FracCovFlux->Write();
      FracCovXS->Write();
    }
    //###########################################################################################    


    //####################################### STAT ERROR ################################
    cout<<"\nStatistical error:"<<endl;

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	double Value=0;double Error=0;double ErrorSquared=0;
	Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
#ifdef RECONSTRUCTED    
	ErrorSquared=Value;
#else
	ErrorSquared=CovarianceStatistics[e0][e1][e0][e1];
#endif
	Error=TMath::Sqrt(ErrorSquared);

	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalStatistics_Plus[e0][e1]=RelativeError;
	ErrorTotalStatistics_Minus[e0][e1]=-RelativeError;

	cout<<e0<<" "<<e1<<" "<<Value<<" "<<Error<<" "<<ErrorSquared<<endl;      
	NominalSelectedMC->SetBinError(e0+1,e1+1,ErrorSquared);

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
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;

      double Value=0;double Error=0;double ErrorSquared=0;
      double Value_ch=0,Value_h2o=0;
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	Value+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	if(ratio){
	  Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	  Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	}
      }
      if(ratio) Value=Value_h2o/Value_ch;

#ifndef RECONSTRUCTED    
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	  if(ratio){
	    ErrorSquared+=CovarianceStatisticsTargets[e0][e1][e0][f1][0][0]/pow(Value_h2o,2.);
	    ErrorSquared+=CovarianceStatisticsTargets[e0][e1][e0][f1][1][1]/pow(Value_ch,2.);
	    ErrorSquared+= -CovarianceStatisticsTargets[e0][e1][e0][f1][1][0]/(Value_ch*Value_h2o);
	    ErrorSquared+= -CovarianceStatisticsTargets[e0][e1][e0][f1][0][1]/(Value_ch*Value_h2o);
	  }
	  else  ErrorSquared+=CovarianceStatistics[e0][e1][e0][f1];
	}
      }
      if(ratio) ErrorSquared=ErrorSquared*pow(Value,2.);
#else      
      ErrorSquared=Value;
#endif
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
      if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      double Value=0;double Error=0;double ErrorSquared=0;
      double Value_ch=0,Value_h2o=0;
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	Value+=NominalSelectedMC->GetBinContent(e0+1,e1+1)/**(BinningTrueMom[e0+1]-BinningTrueMom[e0])*/;
	if(ratio){
	  Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	  Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	}
      }
      if(ratio) Value=Value_h2o/Value_ch;

#ifndef RECONSTRUCTED
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	  if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	  if(ratio){
	    ErrorSquared+=CovarianceStatisticsTargets[e0][e1][f0][e1][0][0]/pow(Value_h2o,2.);
	    ErrorSquared+=CovarianceStatisticsTargets[e0][e1][f0][e1][1][1]/pow(Value_ch,2.);
	    ErrorSquared+= -CovarianceStatisticsTargets[e0][e1][f0][e1][1][0]/(Value_ch*Value_h2o);
	    ErrorSquared+= -CovarianceStatisticsTargets[e0][e1][f0][e1][0][1]/(Value_ch*Value_h2o);
	  }
	  else  ErrorSquared+=CovarianceStatistics[e0][e1][f0][e1];
	}
      }
      if(ratio) ErrorSquared=ErrorSquared*pow(Value,2.);
#else
      ErrorSquared=Value;
#endif
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
    }

    double TotalErrorStatSquared=0,TotalRelErrorStat=0.;
#ifndef RECONSTRUCTED
    for(int e0=0;e0<NBinsMom;e0++){
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      for(int e1=0;e1<NBinsAngle;e1++){
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	for(int f0=0;f0<NBinsMom;f0++){
	  if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	  for(int f1=0;f1<NBinsAngle;f1++){
	    if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	    if(ratio){
	      TotalErrorStatSquared+=CovarianceStatisticsTargets[e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
	      TotalErrorStatSquared+=CovarianceStatisticsTargets[e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
	      TotalErrorStatSquared+= -CovarianceStatisticsTargets[e0][e1][f0][f1][1][0]/(TotalValue_ch*TotalValue_h2o);
	      TotalErrorStatSquared+= -CovarianceStatisticsTargets[e0][e1][f0][f1][0][1]/(TotalValue_ch*TotalValue_h2o);
	    }
	    else TotalErrorStatSquared+=CovarianceStatistics[e0][e1][f0][f1];
	  }
	}
      }
    }
    if(ratio) TotalErrorStatSquared=TotalErrorStatSquared*pow(TotalValue,2.);
#else
    TotalErrorStatSquared=TotalValue;
#endif
    if(TotalValue!=0) TotalRelErrorStat=TMath::Sqrt(TotalErrorStatSquared)/TotalValue;
    
 
    //#########################################END OF STAT ERROR################################

    //################# ML newly added 2017/07/18  DETECTOR ERROR  #############################
    double TotalErrorDetSquared=0.,TotalRelErrorDet=0.;
    double TotaldetErrorSquared[Systematics_Detector_End+1]={0};
    double TotalReldetError[Systematics_Detector_End+1]={0};
#ifdef DETECTORSYST
    cout<<"Detector error:"<<endl;
    
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	double Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	double Value_ch=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	double Value_h2o=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);

	double total_error_plus=0,total_error_minus=0;
	double total_error=0;
	for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	  if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	  if(!useKogaDetError && !IsReadyDetSyst(ErrorType)) continue;
	  if(!ratio) {
	    if(ErrorType==5){
	      total_error+=fmax(CovarianceBirksPlus[e0][e1][e0][e1],CovarianceBirksMinus[e0][e1][e0][e1]);
	    }
	    else total_error+=CovarianceReduced[ErrorType][e0][e1][e0][e1];
	  }
	  else {
	    if(ErrorType==5){
	      double totalBirksPlus=CovarianceBirksPlusTargets[e0][e1][e0][e1][0][0]/pow(Value_h2o,2.)
		+CovarianceBirksPlusTargets[e0][e1][e0][e1][1][1]/pow(Value_ch,2.)
		-2*CovarianceBirksPlusTargets[e0][e1][e0][e1][0][1]/(Value_h2o*Value_ch);

	      double totalBirksMinus=CovarianceBirksMinusTargets[e0][e1][e0][e1][0][0]/pow(Value_h2o,2.)
		+CovarianceBirksMinusTargets[e0][e1][e0][e1][1][1]/pow(Value_ch,2.)
		-2*CovarianceBirksMinusTargets[e0][e1][e0][e1][0][1]/(Value_h2o*Value_ch);

	      total_error+=fmax(totalBirksPlus,totalBirksMinus);
	    }
	    else {
	      total_error+=CovarianceReducedTargets[ErrorType][e0][e1][e0][e1][0][0]/pow(Value_h2o,2.);
	      total_error+=CovarianceReducedTargets[ErrorType][e0][e1][e0][e1][1][1]/pow(Value_ch,2.);
	      total_error+=-2*CovarianceReducedTargets[ErrorType][e0][e1][e0][e1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
	    }
	  }
	}
	if(ratio) total_error*=pow(Value,2.);

	ErrorTotalDetector_Plus[e0][e1]=TMath::Sqrt(total_error);// relative error
	ErrorTotalDetector_Minus[e0][e1]=TMath::Sqrt(total_error); // relative error	
	//	cout<<e0<<" "<<e1<<" "<<Value<<" "<<total_error_plus<<" "<<total_error_minus<<" "<<total_error/Value/Value<<endl;

	double ErrorSquared=NominalSelectedMC->GetBinError(e0+1,e1+1);
	
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
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      double Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);

      double Value_ch=0,Value_h2o=0;
      if(ratio){
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	  Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	}
      }

      double total_error=0;
      for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	if(!useKogaDetError && !IsReadyDetSyst(ErrorType)) continue;
	if(ErrorType==5){
	  double errorPlus=0,errorMinus=0;
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	      if(ratio){
		errorPlus+=CovarianceBirksPlusTargets[e0][e1][e0][f1][0][0]/pow(Value_h2o,2.);
		errorPlus+=CovarianceBirksPlusTargets[e0][e1][e0][f1][1][1]/pow(Value_ch,2.);
		errorPlus+=-CovarianceBirksPlusTargets[e0][e1][e0][f1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
		errorPlus+=-CovarianceBirksPlusTargets[e0][e1][e0][f1][1][0]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors

		errorMinus+=CovarianceBirksMinusTargets[e0][e1][e0][f1][0][0]/pow(Value_h2o,2.);
		errorMinus+=CovarianceBirksMinusTargets[e0][e1][e0][f1][1][1]/pow(Value_ch,2.);
		errorMinus+=-CovarianceBirksMinusTargets[e0][e1][e0][f1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
		errorMinus+=-CovarianceBirksMinusTargets[e0][e1][e0][f1][1][0]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
	      }
	      else{
		errorPlus+=CovarianceBirksPlus[e0][e1][e0][f1];
		errorMinus+=CovarianceBirksMinus[e0][e1][e0][f1];
	      }
	    }
	  }
	  total_error+=fmax(errorPlus,errorMinus);
	}
	else {
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	      if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	      
	      if(ratio){
		total_error+=CovarianceReducedTargets[ErrorType][e0][e1][e0][f1][0][0]/pow(Value_h2o,2.);
		total_error+=CovarianceReducedTargets[ErrorType][e0][e1][e0][f1][1][1]/pow(Value_ch,2.);
		total_error+=-CovarianceReducedTargets[ErrorType][e0][e1][e0][f1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
		total_error+=-CovarianceReducedTargets[ErrorType][e0][e1][e0][f1][1][0]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
	      }
	      else total_error+=CovarianceReduced[ErrorType][e0][e1][e0][f1];
	    }
	  }
	}
      }
      if(ratio) total_error*=pow(Value,2.);

      double ErrorSquared=NominalSelectedMC_RecMom->GetBinError(e0+1);

      ErrorSquared+=total_error;
      double Error=TMath::Sqrt(ErrorSquared);

      //cout<<e0<<" "<<Value<<" "<<total_error_plus/Value/Value<<" "<<total_error_minus/Value/Value<<" "<<total_error/Value/Value<<endl;

      ErrorTotalDetector_RecMom_Plus[e0]=(Value>0?TMath::Sqrt(total_error)/Value:0);// relative error
      ErrorTotalDetector_RecMom_Minus[e0]=(Value>0?TMath::Sqrt(total_error)/Value:0); // relative error	
      
      cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", +-"<<TMath::Sqrt(total_error)<<endl;
      NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);

      boxErrorDet_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
      boxErrorDet_RecMom[e0]->SetFillColor(kOrange);
    }

    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
      if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      double Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);

      double Value_ch=0,Value_h2o=0; 
      if(ratio){
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	  if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	  Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	  Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	}
      }
	
      double total_error=0;
      for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
	if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
	if(!useKogaDetError && !IsReadyDetSyst(ErrorType)) continue;
	if(ErrorType==5){
	  double errorMinus=0,errorPlus=0;
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	      if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	      if(ratio){
		errorPlus+=CovarianceBirksPlusTargets[e0][e1][f0][e1][0][0]/pow(Value_h2o,2.);
		errorPlus+=CovarianceBirksPlusTargets[e0][e1][f0][e1][1][1]/pow(Value_ch,2.);
		errorPlus+=-CovarianceBirksPlusTargets[e0][e1][f0][e1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
		errorPlus+=-CovarianceBirksPlusTargets[e0][e1][f0][e1][1][0]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors

		errorMinus+=CovarianceBirksMinusTargets[e0][e1][f0][e1][0][0]/pow(Value_h2o,2.);
		errorMinus+=CovarianceBirksMinusTargets[e0][e1][f0][e1][1][1]/pow(Value_ch,2.);
		errorMinus+=-CovarianceBirksMinusTargets[e0][e1][f0][e1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
		errorMinus+=-CovarianceBirksMinusTargets[e0][e1][f0][e1][1][0]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
	      }
	      else{
		errorPlus+=CovarianceBirksPlus[e0][e1][f0][e1];
		errorMinus+=CovarianceBirksMinus[e0][e1][f0][e1];
	      }
	    }
	  }
	  total_error+=fmax(errorPlus,errorMinus);
	}
	else {
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	    for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	      if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	      if(ratio){
		total_error+=CovarianceReducedTargets[ErrorType][e0][e1][f0][e1][0][0]/pow(Value_h2o,2.);
		total_error+=CovarianceReducedTargets[ErrorType][e0][e1][f0][e1][1][1]/pow(Value_ch,2.);
		total_error+=-CovarianceReducedTargets[ErrorType][e0][e1][f0][e1][0][1]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
		total_error+=-CovarianceReducedTargets[ErrorType][e0][e1][f0][e1][1][0]/(Value_h2o*Value_ch);// it's 0 for non-correlated errors
	      }
	      else total_error+=CovarianceReduced[ErrorType][e0][e1][f0][e1];
	    }
	  }
	}
      }
      if(ratio) total_error*=pow(Value,2.); 

      double ErrorSquared=NominalSelectedMC_RecAngle->GetBinError(e1+1);
      ErrorSquared+=total_error;
      double Error=TMath::Sqrt(ErrorSquared);      

      ErrorTotalDetector_RecAngle_Plus[e1]=(Value>0?TMath::Sqrt(total_error)/Value:0);// relative error
      ErrorTotalDetector_RecAngle_Minus[e1]=(Value>0?TMath::Sqrt(total_error)/Value:0); // relative error	

      NominalSelectedMC_RecAngle->SetBinError(e1+1,ErrorSquared);

      boxErrorDet_RecAngle[e1] = new TBox(NominalSelectedMC_RecAngle->GetXaxis()->GetBinLowEdge(e1+1),Value-Error,NominalSelectedMC_RecAngle->GetXaxis()->GetBinUpEdge(e1+1),Value+Error);
      boxErrorDet_RecAngle[e1]->SetFillColor(kOrange);
    }

    for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      if(ErrorType>Systematics_Detector_End || ErrorType<Systematics_Detector_Start) continue;
      if(!useKogaDetError && !IsReadyDetSyst(ErrorType)) continue;
      if(ErrorType==5){
	double errorPlus=0,errorMinus=0;
	for(int e0=0;e0<NBinsMom;e0++){
	  if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	  for(int e1=0;e1<NBinsAngle;e1++){
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    for(int f0=0;f0<NBinsMom;f0++){
	      if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	      for(int f1=0;f1<NBinsAngle;f1++){
		if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
		if(ratio){
		  errorPlus+=CovarianceBirksPlusTargets[e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
		  errorPlus+=CovarianceBirksPlusTargets[e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
		  errorPlus+=-CovarianceBirksPlusTargets[e0][e1][f0][f1][0][1]/(TotalValue_h2o*TotalValue_ch);// it's 0 for non-correlated errors
		  errorPlus+=-CovarianceBirksPlusTargets[e0][e1][f0][f1][1][0]/(TotalValue_h2o*TotalValue_ch);// it's 0 for non-correlated errors

		  errorMinus+=CovarianceBirksMinusTargets[e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
		  errorMinus+=CovarianceBirksMinusTargets[e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
		  errorMinus+=-CovarianceBirksMinusTargets[e0][e1][f0][f1][0][1]/(TotalValue_h2o*TotalValue_ch);// it's 0 for non-correlated errors
		  errorMinus+=-CovarianceBirksMinusTargets[e0][e1][f0][f1][1][0]/(TotalValue_h2o*TotalValue_ch);// it's 0 for non-correlated errors
		}
		else{
		  errorPlus+=CovarianceBirksPlus[e0][e1][f0][f1];
		  errorMinus+=CovarianceBirksMinus[e0][e1][f0][f1];
		}
	      }
	    }
	  }
	}

	cout<<ErrorType<<" "<<TotalErrorDetSquared<<" "<<fmax(errorPlus,errorMinus)<<" "<<sqrt(fmax(errorPlus,errorMinus))<<" "<<errorPlus<<" "<<errorMinus<<endl; 

	TotalErrorDetSquared+=fmax(errorPlus,errorMinus);
	TotaldetErrorSquared[5]=fmax(errorPlus,errorMinus);
      }
      else {
	double error=0.;
	for(int e0=0;e0<NBinsMom;e0++){
	  if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	  for(int e1=0;e1<NBinsAngle;e1++){
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    for(int f0=0;f0<NBinsMom;f0++){
	      if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	      for(int f1=0;f1<NBinsAngle;f1++){
		if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
		if(ratio){
		  error+=CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
		  error+=CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
		  error+=-CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][0][1]/(TotalValue_h2o*TotalValue_ch);// it's 0 for non-correlated errors
		  error+=-CovarianceReducedTargets[ErrorType][e0][e1][f0][f1][1][0]/(TotalValue_h2o*TotalValue_ch);// it's 0 for non-correlated errors
		}
		else error+=CovarianceReduced[ErrorType][e0][e1][f0][f1];
	      }
	    }
	  }
	}

	cout<<ErrorType<<" "<<TotalErrorDetSquared<<" "<<error<<" "<<sqrt(error)<<endl; 
	
	TotalErrorDetSquared+=error;
	TotaldetErrorSquared[ErrorType]=error;		
      }

      if(ratio)	TotaldetErrorSquared[ErrorType]*=pow(TotalValue,2.); 
      if(TotalValue!=0) TotalReldetError[ErrorType]=TMath::Sqrt(TotaldetErrorSquared[ErrorType])/TotalValue;
    }
    if(ratio) 	TotalErrorDetSquared*=pow(TotalValue,2.);    
    if(TotalValue!=0) TotalRelErrorDet=TMath::Sqrt(TotalErrorDetSquared)/TotalValue;
    
  
    //#########################################END OF DETECTOR ERROR################################
#endif
 

    //#########################################XS ERROR################################
    double TotalErrorXSSquared=0,TotalRelErrorXS=0.;

    if(EndError>=Systematics_Xsec_Start){
      cout<<"\nXsec error:"<<endl;

      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  double Value=0;double Error=0;double ErrorSquared=0;
	  Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  ErrorSquared=CovarianceXS[e0][e1][e0][e1];

	  double RelativeError=TMath::Sqrt(ErrorSquared);
	  if(Value!=0) RelativeError/=Value;      
	  ErrorTotalXS_Plus[e0][e1]=RelativeError;
	  ErrorTotalXS_Minus[e0][e1]=-RelativeError;

	  //	ErrorSquared+=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  ErrorSquared+=NominalSelectedMC->GetBinError(e0+1,e1+1); // ML modified 2017/07/18
	  Error=TMath::Sqrt(ErrorSquared);
	  cout<<e0<<" "<<e1<<" "<<Value<<" "<<Error<<" "<<ErrorSquared<<endl;      
	
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
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;

	double Value=0;double ErrorSquared=0;double Error=0;
	Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);    
	ErrorSquared=NominalSelectedMC_RecMom->GetBinError(e0+1);    

	double Value_ch=0,Value_h2o=0;
	if(ratio){
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	    Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	  }
	}

    	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	    if(ratio){
 	      ErrorSquared+=CovarianceXSTargets[e0][e1][e0][f1][0][0]/pow(Value_h2o,2.);
	      ErrorSquared+=CovarianceXSTargets[e0][e1][e0][f1][1][1]/pow(Value_ch,2.);
	      ErrorSquared+= -CovarianceXSTargets[e0][e1][e0][f1][1][0]/(Value_ch*Value_h2o);
	      ErrorSquared+= -CovarianceXSTargets[e0][e1][e0][f1][0][1]/(Value_ch*Value_h2o);
	    }
	    else  ErrorSquared+=CovarianceXS[e0][e1][e0][f1];
	  }
	}
	if(ratio) ErrorSquared=ErrorSquared*pow(Value,2.);
	
	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalXS_RecMom_Plus[e0]=RelativeError;
	ErrorTotalXS_RecMom_Minus[e0]=-RelativeError;

	Error=TMath::Sqrt(ErrorSquared);
	cout<<"Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<Error<<endl;
	//NominalSelectedMC_RecMom->SetBinContent(e0+1,Value);
	NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);

	boxErrorXS_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),NominalSelectedMC_RecMom->GetBinContent(e0+1)-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),NominalSelectedMC_RecMom->GetBinContent(e0+1)+Error);
	boxErrorXS_RecMom[e0]->SetFillColor(kGreen);
      }


      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      
	double Value=0;double ErrorSquared=0;double Error=0;
	Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);    

	double Value_ch=0,Value_h2o=0;
	if(ratio){
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	    Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	    Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	  }
	}
    
	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	  if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 0
	    if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	    if(ratio){
	      ErrorSquared+=CovarianceXSTargets[e0][e1][f0][e1][0][0]/pow(Value_h2o,2.);
	      ErrorSquared+=CovarianceXSTargets[e0][e1][f0][e1][1][1]/pow(Value_ch,2.);
	      ErrorSquared+= -CovarianceXSTargets[e0][e1][f0][e1][1][0]/(Value_ch*Value_h2o);
	      ErrorSquared+= -CovarianceXSTargets[e0][e1][f0][e1][0][1]/(Value_ch*Value_h2o);
	    }
	    else  ErrorSquared+=CovarianceXS[e0][e1][f0][e1];
	  }
	}
	if(ratio) ErrorSquared=ErrorSquared*pow(Value,2.);

	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalXS_RecAngle_Plus[e1]=RelativeError;
	ErrorTotalXS_RecAngle_Minus[e1]=-RelativeError;
	ErrorSquared+=NominalSelectedMC_RecAngle->GetBinError(e1+1);
 
	Error=TMath::Sqrt(ErrorSquared);
	cout<<"Angle="<<BinningRecAngle[e1]<<", value="<<Value<<", "<<Error<<endl;
	//NominalSelectedMC_RecAngle->SetBinContent(e1+1,Value);
	NominalSelectedMC_RecAngle->SetBinError(e1+1,ErrorSquared);

	boxErrorXS_RecAngle[e1] = new TBox(NominalSelectedMC_RecAngle->GetXaxis()->GetBinLowEdge(e1+1),NominalSelectedMC_RecAngle->GetBinContent(e1+1)-Error,NominalSelectedMC_RecAngle->GetXaxis()->GetBinUpEdge(e1+1),NominalSelectedMC_RecAngle->GetBinContent(e1+1)+Error);
	boxErrorXS_RecAngle[e1]->SetFillColor(kGreen);
      }
  

      for(int e0=0;e0<NBinsMom;e0++){
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	for(int e1=0;e1<NBinsAngle;e1++){
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  for(int f0=0;f0<NBinsMom;f0++){
	    if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	    for(int f1=0;f1<NBinsAngle;f1++){
	      if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	      if(ratio){
		TotalErrorXSSquared+=CovarianceXSTargets[e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
		TotalErrorXSSquared+=CovarianceXSTargets[e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
		TotalErrorXSSquared+= -CovarianceXSTargets[e0][e1][f0][f1][1][0]/(TotalValue_ch*TotalValue_h2o);
		TotalErrorXSSquared+= -CovarianceXSTargets[e0][e1][f0][f1][0][1]/(TotalValue_ch*TotalValue_h2o);
	      }
	      else TotalErrorXSSquared+=CovarianceXS[e0][e1][f0][f1];
	    }
	  }
	}
      }
      if(ratio) TotalErrorXSSquared=TotalErrorXSSquared*pow(TotalValue,2.);
      if(TotalValue!=0) TotalRelErrorXS=TMath::Sqrt(TotalErrorXSSquared)/TotalValue;
    }
    //#########################################END OF XS ERROR################################

    //#########################################FLUX ERROR################################
    double TotalErrorFluxSquared=0,TotalRelErrorFlux=0.;

    if(EndError>=Systematics_Flux_Start){
      cout<<"\nFlux error:"<<endl;

      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  double Value=0;double Error=0;double ErrorSquared=0;
	  Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  ErrorSquared=CovarianceFlux[e0][e1][e0][e1];
	  // antinu, nue component is already included in the Covariance Matrix
	
	  //cout<<setprecision(3)<<std::scientific<<e0<<", "<<e1<<", "<<ErrorSquared<<endl;
	
	  double RelativeError=TMath::Sqrt(ErrorSquared);
	  if(Value!=0) RelativeError/=Value;      
	  ErrorTotalFlux_Plus[e0][e1]=RelativeError;
	  ErrorTotalFlux_Minus[e0][e1]=-RelativeError;


	  // ErrorSquared+=NominalSelectedMC->GetBinContent(e0+1,e1+1); 
	  ErrorSquared+=NominalSelectedMC->GetBinError(e0+1,e1+1); // ML modified 2017/07/18
	  Error=TMath::Sqrt(ErrorSquared);
	  cout<<e0<<" "<<e1<<" "<<Value<<" +/-"<<RelativeError*Value<<" ("<<RelativeError*100<<"%) and from fit "<<Error_Plus[Systematics_Flux_Start]->GetBinContent(e0+1,e1+1)*100<<"%"<<endl;      
	
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
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	double Value=0;double ErrorSquared=0;double Error=0;
	Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);

	double Value_ch=0,Value_h2o=0;
	if(ratio){
	  for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	    Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	    Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	  }
	}
	
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  for(int f1=0;f1<NBinsAngle;f1++){//loop over effect 1
	    if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	    if(ratio){
	      ErrorSquared+=CovarianceFluxTargets[e0][e1][e0][f1][0][0]/pow(Value_h2o,2.);
	      ErrorSquared+=CovarianceFluxTargets[e0][e1][e0][f1][1][1]/pow(Value_ch,2.);
	      ErrorSquared+= -CovarianceFluxTargets[e0][e1][e0][f1][1][0]/(Value_ch*Value_h2o);
	      ErrorSquared+= -CovarianceFluxTargets[e0][e1][e0][f1][0][1]/(Value_ch*Value_h2o);
	    }
	    else ErrorSquared+=CovarianceFlux[e0][e1][e0][f1];
	  }
	}
	if(ratio) ErrorSquared=ErrorSquared*pow(Value,2.);

	double RelativeError=TMath::Sqrt(ErrorSquared);
	if(Value!=0) RelativeError/=Value;      
	ErrorTotalFlux_RecMom_Plus[e0]=RelativeError;
	ErrorTotalFlux_RecMom_Minus[e0]=-RelativeError;
	ErrorSquared+=NominalSelectedMC_RecMom->GetBinError(e0+1);

	Error=TMath::Sqrt(ErrorSquared);

	cout<<"Flux error, Mom="<<BinningRecMom[e0]<<", value="<<Value<<", "<<Error<<endl;
	boxErrorFlux_RecMom[e0] = new TBox(NominalSelectedMC_RecMom->GetXaxis()->GetBinLowEdge(e0+1),Value-Error,NominalSelectedMC_RecMom->GetXaxis()->GetBinUpEdge(e0+1),Value+Error);
	boxErrorFlux_RecMom[e0]->SetFillColor(kBlue);
	NominalSelectedMC_RecMom->SetBinError(e0+1,ErrorSquared);
      }
  
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	double Value=0;double ErrorSquared=0;double Error=0;
	Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);
      
	double Value_ch=0,Value_h2o=0;
	if(ratio){
	  for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	    if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	    Value_ch+=NominalSelectedMC_ch->GetBinContent(e0+1,e1+1);
	    Value_h2o+=NominalSelectedMC_h2o->GetBinContent(e0+1,e1+1);
	  }
	}

	for(int e0=0;e0<NBinsMom;e0++){//loop over effect 1
	  if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	  for(int f0=0;f0<NBinsMom;f0++){//loop over effect 1
	    if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	    if(ratio){
	      ErrorSquared+=CovarianceFluxTargets[e0][e1][f0][e1][0][0]/pow(Value_h2o,2.);
	      ErrorSquared+=CovarianceFluxTargets[e0][e1][f0][e1][1][1]/pow(Value_ch,2.);
	      ErrorSquared+= -CovarianceFluxTargets[e0][e1][f0][e1][1][0]/(Value_ch*Value_h2o);
	      ErrorSquared+= -CovarianceFluxTargets[e0][e1][f0][e1][0][1]/(Value_ch*Value_h2o);
	    }
	    else ErrorSquared+=CovarianceFlux[e0][e1][f0][e1];
	  }
	}
	if(ratio) ErrorSquared=ErrorSquared*pow(Value,2.);

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
    
  
      for(int e0=0;e0<NBinsMom;e0++){
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	for(int e1=0;e1<NBinsAngle;e1++){
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  for(int f0=0;f0<NBinsMom;f0++){
	    if(UseTrashBins && IsTrueMomTrashBin[f0]) continue;
	    for(int f1=0;f1<NBinsAngle;f1++){
	      if(UseTrashBins && IsTrueAngleTrashBin[f1]) continue;
	      if(ratio){
		TotalErrorFluxSquared+=CovarianceFluxTargets[e0][e1][f0][f1][0][0]/pow(TotalValue_h2o,2.);
		TotalErrorFluxSquared+=CovarianceFluxTargets[e0][e1][f0][f1][1][1]/pow(TotalValue_ch,2.);
		TotalErrorFluxSquared+= -CovarianceFluxTargets[e0][e1][f0][f1][1][0]/(TotalValue_ch*TotalValue_h2o);
		TotalErrorFluxSquared+= -CovarianceFluxTargets[e0][e1][f0][f1][0][1]/(TotalValue_ch*TotalValue_h2o);
	      }
	      else TotalErrorFluxSquared+=CovarianceFlux[e0][e1][f0][f1];
	    }
	  }
	}
      }
      if(ratio) TotalErrorFluxSquared=TotalErrorFluxSquared*pow(TotalValue,2.);
      if(TotalValue!=0) TotalRelErrorFlux=TMath::Sqrt(TotalErrorFluxSquared)/TotalValue;
    }
    //#########################################END FLUX ERROR################################



  
    //#####################################PUT BACK TOTAL ERROR TO SQRT (COVARIANCE)################################
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
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
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      double Value=0;double ErrorSquared=0;double Error=0;
      ErrorSquared=NominalSelectedMC_RecMom->GetBinError(e0+1);
      Error=TMath::Sqrt(ErrorSquared);
      //Error=ErrorSquared;
      NominalSelectedMC_RecMom->SetBinError(e0+1,Error);
    }
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
      if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      double Value=0;double ErrorSquared=0;double Error=0;
      ErrorSquared=NominalSelectedMC_RecAngle->GetBinError(e1+1);
      Error=TMath::Sqrt(ErrorSquared);
      //Error=ErrorSquared;
      NominalSelectedMC_RecAngle->SetBinError(e1+1,Error);
    }
    //#####################################END OF PUT BACK TOTAL ERROR TO SQRT (COVARIANCE)################################


    //#################################### DIVIDE BIN CONTENT AND ERROR BY THE BINWIDTH #######################
    // ML added 2017/12/13 -- only at true level
#ifndef RECONSTRUCTED
    if(normalized && !ratio){
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	  if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	  double BinWidth=(BinningTrueMom[e0+1]-BinningTrueMom[e0])*(BinningTrueAngle[e1+1]-BinningTrueAngle[e1]); 
	  double Value=NominalSelectedMC->GetBinContent(e0+1,e1+1);
	  double Error=NominalSelectedMC->GetBinError(e0+1,e1+1);
	  NominalSelectedMC->SetBinContent(e0+1,e1+1,Value/BinWidth);
	  NominalSelectedMC->SetBinError(e0+1,e1+1,Error/BinWidth);

	  sliceNominalSelectedMC_RecMom[e1]->SetBinContent(e0+1,Value/BinWidth);
	  sliceNominalSelectedMC_RecAngle[e0]->SetBinContent(e1+1,Value/BinWidth);
	  sliceNominalSelectedMC_RecMom[e1]->SetBinError(e0+1,Error/BinWidth);
	  sliceNominalSelectedMC_RecAngle[e0]->SetBinError(e1+1,Error/BinWidth);

	  double TrueValue=NominalTrueMC->GetBinContent(e0+1,e1+1);
	  sliceNominalTrueMC_RecMom[e1]->SetBinContent(e0+1,TrueValue/BinWidth);
	  sliceNominalTrueMC_RecAngle[e0]->SetBinContent(e1+1,TrueValue/BinWidth);

	  // now resize the error boxes
	  boxsliceErrorStat_RecMom[e0][e1]->SetY1(boxsliceErrorStat_RecMom[e0][e1]->GetY1()/BinWidth);
	  boxsliceErrorStat_RecMom[e0][e1]->SetY2(boxsliceErrorStat_RecMom[e0][e1]->GetY2()/BinWidth);
	  boxsliceErrorStat_RecAngle[e0][e1]->SetY1(boxsliceErrorStat_RecAngle[e0][e1]->GetY1()/BinWidth);
	  boxsliceErrorStat_RecAngle[e0][e1]->SetY2(boxsliceErrorStat_RecAngle[e0][e1]->GetY2()/BinWidth);
#ifdef DETECTORSYST
	  boxsliceErrorDet_RecMom[e0][e1]->SetY1(boxsliceErrorDet_RecMom[e0][e1]->GetY1()/BinWidth);
	  boxsliceErrorDet_RecMom[e0][e1]->SetY2(boxsliceErrorDet_RecMom[e0][e1]->GetY2()/BinWidth);
	  boxsliceErrorDet_RecAngle[e0][e1]->SetY1(boxsliceErrorDet_RecAngle[e0][e1]->GetY1()/BinWidth);
	  boxsliceErrorDet_RecAngle[e0][e1]->SetY2(boxsliceErrorDet_RecAngle[e0][e1]->GetY2()/BinWidth);
#endif
	  if(EndError>=Systematics_Xsec_Start){
	    boxsliceErrorXS_RecMom[e0][e1]->SetY1(boxsliceErrorXS_RecMom[e0][e1]->GetY1()/BinWidth);
	    boxsliceErrorXS_RecMom[e0][e1]->SetY2(boxsliceErrorXS_RecMom[e0][e1]->GetY2()/BinWidth);
	    boxsliceErrorXS_RecAngle[e0][e1]->SetY1(boxsliceErrorXS_RecAngle[e0][e1]->GetY1()/BinWidth);
	    boxsliceErrorXS_RecAngle[e0][e1]->SetY2(boxsliceErrorXS_RecAngle[e0][e1]->GetY2()/BinWidth);
	  }
	  if(EndError>=Systematics_Flux_Start){
	    boxsliceErrorFlux_RecMom[e0][e1]->SetY1(boxsliceErrorFlux_RecMom[e0][e1]->GetY1()/BinWidth);
	    boxsliceErrorFlux_RecMom[e0][e1]->SetY2(boxsliceErrorFlux_RecMom[e0][e1]->GetY2()/BinWidth);
	    boxsliceErrorFlux_RecAngle[e0][e1]->SetY1(boxsliceErrorFlux_RecAngle[e0][e1]->GetY1()/BinWidth);
	    boxsliceErrorFlux_RecAngle[e0][e1]->SetY2(boxsliceErrorFlux_RecAngle[e0][e1]->GetY2()/BinWidth);
	  }
	}
      }
  

      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	double BinWidth=(BinningTrueMom[e0+1]-BinningTrueMom[e0]);
	double Value=NominalSelectedMC_RecMom->GetBinContent(e0+1);
	double Error=NominalSelectedMC_RecMom->GetBinError(e0+1);

	NominalSelectedMC_RecMom->SetBinContent(e0+1,Value/BinWidth);
	NominalSelectedMC_RecMom->SetBinError(e0+1,Error/BinWidth);

	double TrueValue=NominalTrueMC_RecMom->GetBinContent(e0+1);
	NominalTrueMC_RecMom->SetBinContent(e0+1,TrueValue/BinWidth);

	// now resize the error boxes
	boxErrorStat_RecMom[e0]->SetY1(boxErrorStat_RecMom[e0]->GetY1()/BinWidth);
	boxErrorStat_RecMom[e0]->SetY2(boxErrorStat_RecMom[e0]->GetY2()/BinWidth);
#ifdef DETECTORSYST
	boxErrorDet_RecMom[e0]->SetY1(boxErrorDet_RecMom[e0]->GetY1()/BinWidth);
	boxErrorDet_RecMom[e0]->SetY2(boxErrorDet_RecMom[e0]->GetY2()/BinWidth);
#endif
	if(EndError>=Systematics_Xsec_Start){
	  boxErrorXS_RecMom[e0]->SetY1(boxErrorXS_RecMom[e0]->GetY1()/BinWidth);
	  boxErrorXS_RecMom[e0]->SetY2(boxErrorXS_RecMom[e0]->GetY2()/BinWidth);
	}
	if(EndError>=Systematics_Flux_Start){
	  boxErrorFlux_RecMom[e0]->SetY1(boxErrorFlux_RecMom[e0]->GetY1()/BinWidth);
	  boxErrorFlux_RecMom[e0]->SetY2(boxErrorFlux_RecMom[e0]->GetY2()/BinWidth);
	}
      }


      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	double BinWidth=(BinningTrueAngle[e1+1]-BinningTrueAngle[e1]);
	double Value=NominalSelectedMC_RecAngle->GetBinContent(e1+1);
	double Error=NominalSelectedMC_RecAngle->GetBinError(e1+1);

	NominalSelectedMC_RecAngle->SetBinContent(e1+1,Value/BinWidth);
	NominalSelectedMC_RecAngle->SetBinError(e1+1,Error/BinWidth);

	double TrueValue=NominalTrueMC_RecAngle->GetBinContent(e1+1);
	NominalTrueMC_RecAngle->SetBinContent(e1+1,TrueValue/BinWidth);

	// now resize the error boxes
	boxErrorStat_RecAngle[e1]->SetY1(boxErrorStat_RecAngle[e1]->GetY1()/BinWidth);
	boxErrorStat_RecAngle[e1]->SetY2(boxErrorStat_RecAngle[e1]->GetY2()/BinWidth);
#ifdef DETECTORSYST
	boxErrorDet_RecAngle[e1]->SetY1(boxErrorDet_RecAngle[e1]->GetY1()/BinWidth);
	boxErrorDet_RecAngle[e1]->SetY2(boxErrorDet_RecAngle[e1]->GetY2()/BinWidth);
#endif
	if(EndError>=Systematics_Xsec_Start){
	  boxErrorXS_RecAngle[e1]->SetY1(boxErrorXS_RecAngle[e1]->GetY1()/BinWidth);
	  boxErrorXS_RecAngle[e1]->SetY2(boxErrorXS_RecAngle[e1]->GetY2()/BinWidth);
	}
	if(EndError>=Systematics_Flux_Start){
	  boxErrorFlux_RecAngle[e1]->SetY1(boxErrorFlux_RecAngle[e1]->GetY1()/BinWidth);
	  boxErrorFlux_RecAngle[e1]->SetY2(boxErrorFlux_RecAngle[e1]->GetY2()/BinWidth);
	}
      }
    }
#endif
    //############################### END OF DIVIDE BIN CONTENT AND ERROR BY THE BINWIDTH #######################


  
    //#####################################DRAWING OF 1D MOM & ANGLE PLOTS################################
  
    TCanvas * csliceRecMom[NBinsAngle];
    TCanvas * csliceRecAngle[NBinsMom];

    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
      if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      csliceRecMom[e1]= new TCanvas(Form("csliceRecMom[%d]",e1));
      sliceNominalSelectedMC_RecMom[e1]->Draw("E1");
      for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
	if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
	if(EndError>=Systematics_Flux_Start) boxsliceErrorFlux_RecMom[e0][e1]->Draw("same");
	if(EndError>=Systematics_Xsec_Start) boxsliceErrorXS_RecMom[e0][e1]->Draw("same");
#ifdef DETECTORSYST
	if(EndError>=Systematics_Detector_Start) boxsliceErrorDet_RecMom[e0][e1]->Draw("same");
#endif
	boxsliceErrorStat_RecMom[e0][e1]->Draw("same");
      }
      sliceNominalTrueMC_RecMom[e1]->SetLineStyle(2);
      sliceNominalTrueMC_RecMom[e1]->SetLineWidth(2);
      sliceNominalTrueMC_RecMom[e1]->Draw("same");
      sliceNominalSelectedMC_RecMom[e1]->Draw("E1same");
#ifdef RECONSTRUCTED
      sliceNominalSelectedMC_RecMom[e1]->GetXaxis()->SetTitle("d_{#mu} (cm)");
      sliceNominalSelectedMC_RecMom[e1]->GetYaxis()->SetTitle("Number of events");
#else
      sliceNominalSelectedMC_RecMom[e1]->GetXaxis()->SetTitle("p_{#mu} (GeV)");
      sliceNominalSelectedMC_RecMom[e1]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{#mu} d#theta_{#mu}} (cm^{2}/nucleon/GeV/#circ)");
#endif    
      sliceNominalSelectedMC_RecMom[e1]->SetLineColor(1);
      sliceNominalSelectedMC_RecMom[e1]->SetLineWidth(2);
      sliceNominalSelectedMC_RecMom[e1]->GetYaxis()->SetTitleOffset(1.3);
      csliceRecMom[e1]->Write(Form("Canvas_Nominal_RecMom_slice%d",e1));
    }

    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      csliceRecAngle[e0]= new TCanvas(Form("csliceRecAngle[%d]",e0));
      sliceNominalSelectedMC_RecAngle[e0]->Draw("E1");
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
	if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
	if(EndError>=Systematics_Flux_Start) boxsliceErrorFlux_RecAngle[e0][e1]->Draw("same");
	if(EndError>=Systematics_Xsec_Start) boxsliceErrorXS_RecAngle[e0][e1]->Draw("same");
#ifdef DETECTORSYST
	if(EndError>=Systematics_Detector_Start) boxsliceErrorDet_RecAngle[e0][e1]->Draw("same");
#endif
	boxsliceErrorStat_RecAngle[e0][e1]->Draw("same");
      }
      sliceNominalTrueMC_RecAngle[e0]->SetLineStyle(2);
      sliceNominalTrueMC_RecAngle[e0]->SetLineWidth(2);
      sliceNominalTrueMC_RecAngle[e0]->Draw("same");
      sliceNominalSelectedMC_RecAngle[e0]->Draw("E1same");
#ifdef RECONSTRUCTED
      sliceNominalSelectedMC_RecAngle[e0]->GetYaxis()->SetTitle("Number of events");
#else
      sliceNominalSelectedMC_RecAngle[e0]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{#mu} d#theta_{#mu}} (cm^{2}/nucleon/GeV/#circ)");
#endif    
      sliceNominalSelectedMC_RecAngle[e0]->GetXaxis()->SetTitle("#theta_{#mu} (GeV)");
      sliceNominalSelectedMC_RecAngle[e0]->SetLineColor(1);
      sliceNominalSelectedMC_RecAngle[e0]->SetLineWidth(2);
      sliceNominalSelectedMC_RecAngle[e0]->GetYaxis()->SetTitleOffset(1.3);
      csliceRecAngle[e0]->Write(Form("Canvas_Nominal_RecAngle_slice%d",e0));
    }

    TCanvas * cRecMom = new TCanvas("cRecMom");
    NominalSelectedMC_RecMom->Draw("E1");
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      if(UseTrashBins && IsTrueMomTrashBin[e0]) continue;
      if(EndError>=Systematics_Flux_Start) boxErrorFlux_RecMom[e0]->Draw("same");
      if(EndError>=Systematics_Xsec_Start) boxErrorXS_RecMom[e0]->Draw("same");
#ifdef DETECTORSYST
      if(EndError>=Systematics_Detector_Start) boxErrorDet_RecMom[e0]->Draw("same");
#endif
      boxErrorStat_RecMom[e0]->Draw("same");
    }
    NominalTrueMC_RecMom->SetLineStyle(2);
    NominalTrueMC_RecMom->SetLineWidth(2);
    NominalTrueMC_RecMom->Draw("same");
    NominalSelectedMC_RecMom->Draw("E1same");
#ifdef RECONSTRUCTED
    NominalSelectedMC_RecMom->GetXaxis()->SetTitle("d_{#mu} (cm)");
    NominalSelectedMC_RecMom->GetYaxis()->SetTitle("Number of events");
#else
    NominalSelectedMC_RecMom->GetXaxis()->SetTitle("p_{#mu} (GeV)");
    NominalSelectedMC_RecMom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} (cm^{2}/nucleon/GeV)");
#endif    

    NominalSelectedMC_RecMom->SetLineColor(1);
    NominalSelectedMC_RecMom->SetLineWidth(2);
    NominalSelectedMC_RecMom->GetYaxis()->SetTitleOffset(1.3);
    cRecMom->Write("Canvas_Nominal_RecMom");
  
    TCanvas * cRecAngle = new TCanvas("cRecAngle");
    NominalSelectedMC_RecAngle->Draw("E1");
    for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 0
      if(UseTrashBins && IsTrueAngleTrashBin[e1]) continue;
      if(EndError>=Systematics_Flux_Start) boxErrorFlux_RecAngle[e1]->Draw("same");
      if(EndError>=Systematics_Xsec_Start) boxErrorXS_RecAngle[e1]->Draw("same");
#ifdef DETECTORSYST
      if(EndError>=Systematics_Detector_Start) boxErrorDet_RecAngle[e1]->Draw("same");
#endif
      boxErrorStat_RecAngle[e1]->Draw("same");
    }
    NominalTrueMC_RecAngle->SetLineStyle(2);
    NominalTrueMC_RecAngle->SetLineWidth(2);
    NominalTrueMC_RecAngle->Draw("same");
    NominalSelectedMC_RecAngle->Draw("E1same");
    NominalSelectedMC_RecAngle->SetLineColor(1);
    NominalSelectedMC_RecAngle->SetLineWidth(2);
    NominalSelectedMC_RecAngle->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
#ifdef RECONSTRUCTED
    NominalSelectedMC_RecAngle->GetYaxis()->SetTitle("Number of events");
#else
    NominalSelectedMC_RecAngle->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta_{#mu}} (cm^{2}/nucleon/#circ)");
#endif
    NominalSelectedMC_RecAngle->GetYaxis()->SetTitleOffset(1.3);
  
    cRecAngle->Write("Canvas_Nominal_RecAngle");
    //#####################################END OF DRAWING################################

    //#########################################END OF DRAWING PART################################


    //################################ CHI2 COMPUTATIONS #################################
    double chi2stat=0., chi2bias=0., chi2biascorr=0., chi2rel=0., chi2relcorr=0.;
    if(FakeData){
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  if(UseTrashBins && (IsTrueMomTrashBin[c0] || IsTrueAngleTrashBin[c1])) continue;
	  double TrueC=NominalTrueMC->GetBinContent(c0+1,c1+1);
	  double BiasC=NominalSelectedMC->GetBinContent(c0+1,c1+1)-TrueC;
	  double IterBiasC=NominalSelectedMC->GetBinContent(c0+1,c1+1)-PreviousIterUnfolded[c0][c1];
	  int BinX=(c0-1)*(NBinsAngle-1)+c1;
	  
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over dause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      if(UseTrashBins && (IsTrueMomTrashBin[d0] || IsTrueAngleTrashBin[d1])) continue;
	      double TrueD=NominalTrueMC->GetBinContent(d0+1,d1+1);
	      double BiasD=NominalSelectedMC->GetBinContent(d0+1,d1+1)-TrueD;
	      double IterBiasD=NominalSelectedMC->GetBinContent(d0+1,d1+1)-PreviousIterUnfolded[d0][d1];
	      int BinY=(d0-1)*(NBinsAngle-1)+d1;
	      // if(BinX==0)  cout<<"BinX=0: "<<d0<<" "<<d1<<" "<<BinY<<" "<<BiasD<<" "<<TrueD<<" "<<chi2stat<<endl;

	      chi2stat+=(BiasC)*InvCovMatrix(BinX,BinY)*(BiasD);
	      chi2biascorr+=BiasC/sqrt(TrueC)*InvCorrMatrix(BinX,BinY)*BiasD/sqrt(TrueD);
	      if(iter>0) chi2relcorr+=IterBiasC/sqrt(PreviousIterUnfolded[c0][c1])*InvCorrMatrix(BinX,BinY)*IterBiasD/sqrt(PreviousIterUnfolded[d0][d1]);
	    }
	  }
	  // cout<<c0<<" "<<c1<<" "<<BinX<<" "<<BiasC<<" "<<TrueC<<" "<<chi2stat<<endl;
 
	  chi2bias+=BiasC*BiasC/TrueC;
	  if(iter>0) chi2rel+=IterBiasC*IterBiasC/PreviousIterUnfolded[c0][c1];
	}
      }
      chi2stat/=NBins;  //cout<<chi2stat<<" for "<<NBins<<" bins"<<endl;
      chi2bias/=NBins;  
      chi2biascorr/=NBins;  
      chi2rel/=NBins;  
      chi2relcorr/=NBins;  

      Chi2Stat->SetBinContent(iter+1,chi2stat);
      Chi2Bias->SetBinContent(iter+1,chi2bias);
      Chi2BiasCorr->SetBinContent(iter+1,chi2biascorr);
      Chi2Rel->SetBinContent(iter+1,chi2rel);
      Chi2RelCorr->SetBinContent(iter+1,chi2relcorr);
    }

    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	PreviousIterUnfolded[c0][c1]=NominalSelectedMC->GetBinContent(c0+1,c1+1);
      }
    }

    //#########################################ERROR TABLE####################################
    double power=(normalized && !ratio?1e42:1);
#ifdef RECONSTRUCTED
    power=1.;
#endif 

    cout<<"2D error table"<<endl;  
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsAngle;e1++){//loop over effect 1
	cout<<setprecision(2)<<std::fixed<<power*NominalSelectedMC->GetBinContent(e0+1,e1+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_Plus[e0][e1])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_Plus[e0][e1])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_Plus[e0][e1])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_Plus[e0][e1])<<"(stat),   ";
      }
      cout<<endl;
    }
    cout<<endl<<"1D Momentum error table"<<endl;  
    for(int e0=0;e0<NBinsMom;e0++){//loop over effect 0
      cout<<setprecision(2)<<std::fixed<<power*NominalSelectedMC_RecMom->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecMom_Plus[e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecMom_Plus[e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_RecMom_Plus[e0])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_RecMom_Plus[e0])<<"(stat)"<<endl;
    }
  
    cout<<endl<<"1D Angle error table"<<endl;  
    for(int e0=0;e0<NBinsAngle;e0++){//loop over effect 0
      cout<<setprecision(2)<<std::fixed<<power*NominalSelectedMC_RecAngle->GetBinContent(e0+1)<<"+-"<<TMath::Abs(ErrorTotalFlux_RecAngle_Plus[e0])<<"(Flux)+-"<<TMath::Abs(ErrorTotalXS_RecAngle_Plus[e0])<<"(XS)+-"<<TMath::Abs(ErrorTotalDetector_RecAngle_Plus[e0])<<"(Det)+-"<<TMath::Abs(ErrorTotalStatistics_RecAngle_Plus[e0])<<"(stat)"<<endl;
    }

    cout<<endl<<"1bin result"<<endl;
    cout<<setprecision(2)<<std::fixed<<power*TotalValue<<"+-"<<100*TMath::Abs(TotalRelErrorFlux)<<"%(Flux)+-"<<100*TMath::Abs(TotalRelErrorXS)<<"%(XS)+-"<<100*TMath::Abs(TotalRelErrorDet)<<"%(Det)+-"<<100*TMath::Abs(TotalRelErrorStat)<<"%(stat)"<<endl;
    cout<<setprecision(4)<<"chi2/ndf = "<<chi2<<" ("<<NBins<<" bins, "<<iter<<" iterations)"<<endl;

#ifdef DETECTORSYST
    cout<<"\ndetector systematics breakdown:"<<endl; 
    for(int ierr=StartError;ierr<=EndError;ierr++){
      if(ierr>=Systematics_Detector_Start && ierr<=Systematics_Detector_End)
	cout<<setprecision(2)<<"error "<<ierr<<" +-"<<100*TMath::Abs(TotalReldetError[ierr])<<"%"<<endl;
    }
#endif

    if(EndError>=Systematics_Xsec_Start){
      cout<<"\nxscn systematics breakdown:"<<endl; 
      for(int i=StartXsec;i<=EndXsec;i++)
	cout<<setprecision(2)<<"param "<<i<<" +-"<<100*TMath::Abs(TotalRelxsError[i])<<"%"<<endl;
      if(CorrelateMEC) cout<<"MEC-C and -O are treated as "<<(int)100*MEC_C_O_corr<<"% correlated in error #4"<<endl;
    }
    //#########################################END ERROR TABLE####################################

    
    NominalSelectedMC_RecMom->Write();
    NominalSelectedMC_RecAngle->Write();
    NominalSelectedData->Write();
    NominalSelectedMC->Write();
    NominalTrueMC->Write();
    NominalTrueMC_RecMom->Write();
    NominalTrueMC_RecAngle->Write();
    NominalMCBkg->Write();
    TH1D * NominalMCBkg_RecMom = (TH1D*) NominalMCBkg->ProjectionX("NominalMCBkg_RecMom",0,NominalMCBkg->GetNbinsX());
    NominalMCBkg_RecMom->Write();
    //NominalMCBkg->Sumw2();
    //NominalSelectedMC->Sumw2();
    NominalMCBkg->Divide(NominalSelectedMC);
    NominalMCBkg->Write("ContaminationMC");
    h_test->Write();
    h_test_red->Write();
    file->Close();

  }

  TFile * globalout=new TFile(Form("chi2_FD%d_%s.root",FakeDataSet,DetName),"recreate");
  globalout->cd();
  Chi2Stat->GetYaxis()->SetRangeUser(0,10);
  Chi2Stat->GetYaxis()->SetTitle("#chi^{2}/NBins");
  Chi2Stat->GetXaxis()->SetTitle("Number of iterations");
  Chi2Stat->SetMarkerColor(kBlue);
  Chi2Stat->SetMarkerStyle(20);
  Chi2Stat->Write();

  Chi2Bias->GetYaxis()->SetRangeUser(0,10);
  Chi2Bias->GetYaxis()->SetTitle("#chi^{2}/NBins");
  Chi2Bias->GetXaxis()->SetTitle("Number of iterations");
  Chi2Bias->SetMarkerColor(kRed);
  Chi2Bias->SetMarkerStyle(2);
  Chi2Bias->Write();

  Chi2BiasCorr->GetYaxis()->SetRangeUser(0,10);
  Chi2BiasCorr->GetYaxis()->SetTitle("#chi^{2}/NBins");
  Chi2BiasCorr->GetXaxis()->SetTitle("Number of iterations");
  Chi2BiasCorr->SetMarkerColor(kRed);
  Chi2BiasCorr->SetMarkerStyle(20);
  Chi2BiasCorr->Write();

  Chi2Rel->GetYaxis()->SetRangeUser(0,10);
  Chi2Rel->GetYaxis()->SetTitle("#chi^{2}/NBins");
  Chi2Rel->GetXaxis()->SetTitle("Number of iterations");
  Chi2Rel->SetMarkerColor(kBlack);
  Chi2Rel->SetMarkerStyle(2);
  Chi2Rel->Write();

  Chi2RelCorr->GetYaxis()->SetRangeUser(0,10);
  Chi2RelCorr->GetYaxis()->SetTitle("#chi^{2}/NBins");
  Chi2RelCorr->GetXaxis()->SetTitle("Number of iterations");
  Chi2RelCorr->SetMarkerColor(kBlack);
  Chi2RelCorr->SetMarkerStyle(20);
  Chi2RelCorr->Write();

  globalout->Close();

  return 0;
}


//To correct: covariance should be 1/(n-1) & add a plot for each bin for each different source, not only for the total
//How to set the number of files of flux? In a card, that would be ideal! Call this NFluxFiles, add also NXsecSources. For now, the latter should be <25.
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
#include <TF2.h>
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
//#include "Reconstruction.cc"
#include "Xsec.cc"
//#define DEBUG
#define PLOTS
bool SideBand = false;
const int NBinsIteration=1001;
const int NBinsSystToys=1;
bool UseXS=false;//true:study XS, false: study number of events
bool Fit=false;
//Keep in mind that value is divided by sqrt, so to obtain  the absolute number of event, multiply by sart(nevt)

//////////////////////////////////////////////

int main(int argc, char ** argv){

  //TApplication theApp("App",0,0);
  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat(kFALSE);
  //  Xsec * XS = new Xsec();
  //XS->Xsec::Initialize();
  InitializeGlobal();
#ifdef DEBUG
  cout<<"Initialized"<<endl;
#endif
  //To modify the bining if we have a side band
  ReinitializeUnfoldingBinning(SideBand);
#ifdef DEBUG
  cout<<"ReInitialized"<<endl;
#endif
  int NBinsTrueMomPlots=0;int NBinsTrueAnglePlots=0;
  int FirstBinMomPlots=1; int FirstBinAnglePlots=0;  
//int StartBinTrueMomPlots, StartBinTrueAnglePlots;
  if(SideBand){//We do not care to watch both trash bin of signal, the CC1pi and the Other interaction result.
    //StartBinTrueMomPlots = 1;
    //StartBinTrueAnglePlots = 1;
    NBinsTrueMomPlots = NBinsTrueMom-2;
    NBinsTrueAnglePlots = NBinsTrueAngle-2;
  }
  else{
    //StartBinTrueMomPlots = 1;
    //StartBinTrueAnglePlots = 1;
    NBinsTrueMomPlots = NBinsTrueMomSignal-1;
    NBinsTrueAnglePlots = NBinsTrueAngleSignal-1;
  }

  int ListBinPlots[NBinsTrueMomPlots][NBinsTrueAnglePlots][2];
  bool UsedBinPlots[NBinsTrueMom][NBinsTrueAngle];
  bool UsedBinPlotsMom[NBinsTrueMom];
  bool UsedBinPlotsAngle[NBinsTrueAngle];

  int CurrentBinMom=0;
  
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    bool MomUsed=false;
    
    int CurrentBinAngle=0;
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      bool AngleUsed=false;
      if(c0 < FirstBinMomPlots || c0 >= (NBinsTrueMomSignal -1) || c1 < FirstBinAnglePlots || c1 >= (NBinsTrueAngleSignal-1)){// first mom bin is throw, last signal mom bin is throw, SB is throw, last angle bin is throw
	UsedBinPlots[c0][c1] = false;
	continue;
      }
      else{
	UsedBinPlots[c0][c1] = true;
	ListBinPlots[CurrentBinMom][CurrentBinAngle][0] = c0;
	ListBinPlots[CurrentBinMom][CurrentBinAngle][1] = c1;
	MomUsed=true;
	AngleUsed=true;
      }
      if(AngleUsed) CurrentBinAngle++;
    }
    
    if(MomUsed) CurrentBinMom++;
  }

  for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
    UsedBinPlotsMom[c0] = false;
    for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
      if(UsedBinPlots[c0][c1]){
	UsedBinPlotsMom[c0] = true;
	cout<<"Mom bin #"<<c0<<" is used"<<endl;
	break;
      }
    }
  }
  for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
    UsedBinPlotsAngle[c1] = false;
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      if(UsedBinPlots[c0][c1]){
	UsedBinPlotsAngle[c1] = true;
	cout<<"Ang bin #"<<c1<<" is used"<<endl;
	break;
      }
    }
  }

  static const int ccolor[] = {632, 800, 400, 416+2, 416, 416-8, 600, 616+2, 432, 1, 920, 616};
  std::vector <int> color (ccolor, ccolor + sizeof(ccolor) / sizeof(ccolor[0]));
  int ic=0;
  while(color.size() <NBinsTrueMom*NBinsTrueAngle){
    color.push_back(921+ic);
    ic++;
  }

  Xsec * XS = new Xsec();
  gStyle->SetPaintTextFormat("2.2f");
  
  double vLikelihood[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double vUnfolding[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double vInitialPriorMC[NBinsTrueMom][NBinsTrueAngle];
  double vInitialPrior[NBinsTrueMom][NBinsTrueAngle];
  double vPrior[NBinsTrueMom][NBinsTrueAngle];
  double vPriorNormalised[NBinsTrueMom][NBinsTrueAngle];
  double vPosterior[NBinsTrueMom][NBinsTrueAngle];

  int c=-1;
  bool Systematics_Flux=false;
  bool Systematics_Xsec=false;
  bool Systematics_Detector=false;
  bool Reference=false;
 char * InputNameMC=new char[256];char * InputNameData=new char[256];
 char * OutputNameMC=new char[256];
 
  while ((c = getopt(argc, argv, "fxdm:r:o:")) != -1) {
    switch(c){
    case 'm':
      InputNameMC=optarg;
      break;
    case 'r':
      Reference=true;//For comparing the data with another data set: relative chi2 between data sets instead of with the true value of a model.
      break;
    case 'o':
      OutputNameMC=optarg;
      break;
    case 'f':
      Systematics_Flux=true;
      break; 
    case 'x':
      Systematics_Xsec=true;
      break; 
    case 'd':
      Systematics_Detector=true;
      break; 
    }
  }


  //###############################################################
  //###############################################################
  //################# 0. OUTPUT FILE PREPARATION ##################
  //###############################################################
  //###############################################################
  TFile * ofile = new TFile(OutputNameMC,"recreate");
  TDirectory * General = ofile->mkdir("General");//Contains general information about convergence
  TDirectory * NoStatOrSystVariation = ofile->mkdir("NoStatOrSystVariation");
  TDirectory * StatAndSystVariations = ofile->mkdir("StatAndSystVariations");//Contains 
  TDirectory * StatisticalVariations = ofile->mkdir("StatisticalVariations");
  TDirectory * SystematicVariations = ofile->mkdir("SystematicVariations");
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
  











  //###############################################################
  //###############################################################
  //################## 1. HISTOGRAM PREPARATION ###################
  //###############################################################
  //###############################################################

  //################## 0. BASIC HISTOGRAMS ###################
  TCanvas * canEvents[NBinsTrueMom][NBinsTrueAngle];
  TH2D * NumberOfEvents;
  TH2D * NumberOfEvents_TrueRec;
  
  //############# 1. NO STAT. OR SYST. VARIATION #############
  TH2D * UnfoldedDistribution = new TH2D("UnfoldedDistribution","Events distribution in true binning after unfolding applied",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * UnfoldedFinalIterationDistribution = new TH2D("UnfoldedFinalIterationDistribution",Form("Events distribution in true binning after unfolding applied %d times",NBinsIteration),NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * TrueDistribution = new TH2D("TrueDistribution","Events distribution in true binning for the true MC information",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * EventReconstructed = new TH2D("EventReconstructed","Events distribution in true binning for the true MC information",NBinsRecMom,BinningRecMom,NBinsRecAngle,BinningRecAngle);
  
  //The sigma here is the expected error from the true bin
  TH1D * EventXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2MCStatErrorXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * EventTotalXIterations;
  TH1D * Chi2MCStatErrorTotalXIterations;
 
  
  //############## 2. STAT. + SYST. VARIATIONS ###############
  //Simply the 1D pull distribution in all bins
  TH1D * pullEvents[NBinsTrueMom][NBinsTrueAngle];
  //Events pull distributions with the number of iterations.
  TH2D * pullEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullEventsmeanXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullEventssigmaXIterations[NBinsTrueMom][NBinsTrueAngle];
  //Events pull distributions with the systematic varied toy experiments.
  TH2D * pullstatEventsXsystEvents[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventsmeanXsystEvents[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventssigmaXsystEvents[NBinsTrueMom][NBinsTrueAngle];
  TH3D * pullstatEventsXsystEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  //
  TH2D * BiasDistribution = new TH2D("BiasDistribution","Bias=Mean of the gaussian we used to fit the pull distribution, given in true variables",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * ErrorDistribution = new TH2D("ErrorDistribution","Error=width of the gaussian we used to fit the pull distribution, given in true variables",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  //Events pull distributions and chi2, only for stat variations!!!!.
  TH1D * BiasXIterations[NBinsSystToys][NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2XIterations[NBinsSystToys][NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2TotalXIterations[NBinsSystToys];

 //To create some chi2 compared with Data Final
  TH1D * Chi2MCStatErrorFinalIteration;//Look the difference only for data after the maximal number of iterations
  TH3D * pullstatEventsXsystEventsXIterationsFinalIteration[NBinsTrueMom][NBinsTrueAngle];
  TH1D * BiasXIterationsFinalIteration[NBinsSystToys][NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2XIterationsFinalIteration[NBinsSystToys][NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2TotalXIterationsFinalIteration[NBinsSystToys];
  TH1D * Chi2DataStatErrorXIterationsFinalIteration[NBinsSystToys][NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2DataStatErrorTotalXIterationsFinalIteration[NBinsSystToys];
  TH1D * CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[NBinsSystToys];

  //To create a chi2 as the standard one, but just for the nominal experiement
  TH1D * Chi2DataStatErrorXIterations[NBinsSystToys][NBinsTrueMom][NBinsTrueAngle];
  TH1D * Chi2DataStatErrorTotalXIterations[NBinsSystToys];
  TH1D * CorrelatedChi2DataStatErrorTotalXIterations[NBinsSystToys];

  
  //############### 3. STAT. VARIATIONS ONLY #################
  TCanvas * canstatEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  //Events pull distributions only on stat variations with the number of iterations.
  TH2D * pullstatEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventsmeanXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventssigmaXIterations[NBinsTrueMom][NBinsTrueAngle];
  //0. Preparation for covariance estimation
  TH2D * CovarianceStat_Iteration1;
  TH2D * CorrelationStat_Iteration1;
  TH2D * CovarianceStat_Iteration3;
  TH2D * CorrelationStat_Iteration3;
  TH2D * CovarianceStat;
  TH2D * CorrelationStat;
  TH2D * NVariationsStat;
 
  TMatrixDSym CovarianceXIterations[NBinsIteration];
  TMatrixDSym CovarianceInvertXIterations[NBinsIteration];
  TMatrixDSym CovarianceXIterationsFinalIteration[NBinsIteration];
  TMatrixDSym CovarianceInvertXIterationsFinalIteration[NBinsIteration];
  
  for(int it=0;it<NBinsIteration;it++){
    CovarianceXIterations[it].ResizeTo((NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
    CovarianceInvertXIterations[it].ResizeTo((NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
    CovarianceXIterationsFinalIteration[it].ResizeTo((NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
    CovarianceInvertXIterationsFinalIteration[it].ResizeTo((NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
  }
  //= new TH3D("CovarianceXIterations","Covariance natrix wrt the number of iterations",NBinsIteration,0,NBinsIteration,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots)); 
  //TH3D * CovarianceXIterationsFinalIterations = new TH3D("CovarianceXIterationsFinalIteration","Covariance natrix wrt the number of iterations",NBinsIteration,0,NBinsIteration,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
  
  //TH3D * CovarianceXIterations = new TH3D("CovarianceXIterations","Covariance natrix wrt the number of iterations",NBinsIteration,0,NBinsIteration,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots)); 
  //TH3D * CovarianceXIterationsFinalIterations = new TH3D("CovarianceXIterationsFinalIteration","Covariance natrix wrt the number of iterations",NBinsIteration,0,NBinsIteration,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots)); 

  
  //############### 4. SYST. VARIATIONS ONLY #################
  TCanvas * cansystEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  //Events pull distributions only on syst variations with the number of iterations.
  TH2D * pullsystEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullsystEventsmeanXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullsystEventssigmaXIterations[NBinsTrueMom][NBinsTrueAngle];


  
  
  
  //################## 4. INITIALIZATION ####################
  int NBins=(NBinsTrueMomPlots)*(NBinsTrueAnglePlots);

  NumberOfEvents = new TH2D("NumberOfEvents","Number of events (true MC information) for each true bin",NBinsTrueMomPlots,0,NBinsTrueMomPlots,NBinsTrueAngle,0,NBinsTrueAngle);
  NumberOfEvents_TrueRec = new TH2D("NumberOfEvents_TrueRec","Number of events (true MC information) for each true & reconstructed bin",NBinsRecMom*NBinsRecAngle,0,NBinsRecMom*NBinsRecAngle,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));


  CovarianceStat = new TH2D("CovarianceStat","CovarianceStat matrix",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
  CorrelationStat = new TH2D("CorrelationStat","CorrelationStat matrix",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
  NVariationsStat = new TH2D("NVariationsStat","",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));

  CovarianceStat_Iteration1 = new TH2D("CovarianceStat_Iteration1","CovarianceStat_Iteration1 matrix",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
  CorrelationStat_Iteration1 = new TH2D("CorrelationStat_Iteration1","CorrelationStat_Iteration1 matrix",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));

  CovarianceStat_Iteration3 = new TH2D("CovarianceStat_Iteration3","CovarianceStat_Iteration3 matrix",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
  CorrelationStat_Iteration3 = new TH2D("CorrelationStat_Iteration3","CorrelationStat_Iteration3 matrix",(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),0,(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));

  
  for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
      //int imom = ListBinPlots[c0][c1][0];
      //int iang = ListBinPlots[c0][c1][1];
      if(!UsedBinPlots[c0][c1]) continue;
      
      NumberOfEvents->SetBinContent(c0+1,c1+1,0);
      for(int d0=0;d0<NBinsTrueMomPlots;d0++){//loop over cause 0
	for(int d1=0;d1<NBinsTrueAnglePlots;d1++){//loop over cause 1
	  //int imom2 = ListBinPlots[c0][c1][0];
	  //int iang2 = ListBinPlots[c0][c1][1];
	  if(!UsedBinPlots[d0][d1]) continue;
	  int BinX=c0*(NBinsTrueAnglePlots)+c1+1;
	  int BinY=d0*(NBinsTrueAnglePlots)+d1+1;
	  CovarianceStat->SetBinContent(BinX,BinY,0);
	  CorrelationStat->SetBinContent(BinX,BinY,0);
	  NVariationsStat->SetBinContent(BinX,BinY,0);

	  CovarianceStat_Iteration1->SetBinContent(BinX,BinY,0);
	  CorrelationStat_Iteration1->SetBinContent(BinX,BinY,0);

	  CovarianceStat_Iteration3->SetBinContent(BinX,BinY,0);
	  CorrelationStat_Iteration3->SetBinContent(BinX,BinY,0);
	}
      }
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  int BinX=e0*(NBinsRecAngle-1)+e1+1;
	  int BinY=c0*(NBinsTrueAngle-1)+c1+1;
	  NumberOfEvents_TrueRec->SetBinContent(BinX,BinY,0);
	}
      }
    }
  }
  //
  

  
  Chi2MCStatErrorTotalXIterations = new TH1D("Chi2MCStatErrorTotalXIterations","Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
  EventTotalXIterations = new TH1D("EventTotalXIterations","Total number of event distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);


  //To create some chi2 compared with Data Final
  Chi2MCStatErrorFinalIteration = new TH1D("Chi2MCStatErrorFinalIteration","Total chi2/NDF distribution as a function of the unfolding iterations",3*NBinsTrueMom*NBinsTrueAngle,0,3*NBinsTrueMom*NBinsTrueAngle);


  //  for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0 
  //for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
  //  if(!UsedBinPlots[c0][c1]) continue;
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0 
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
  

      EventXIterations[c0][c1] = new TH1D(Form("EventXIterations[%d][%d]",c0,c1),"Chi2 distribution as a function of the unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration);
      Chi2MCStatErrorXIterations[c0][c1] = new TH1D(Form("Chi2MCStatErrorXIterations[%d][%d]",c0,c1),"Chi2 distribution as a function of the unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration);

      pullEvents[c0][c1] = new TH1D(Form("pullEvents[%d][%d]",c0,c1),"",100,-10,10);

      pullEventsXIterations[c0][c1] = new TH2D(Form("pullEventsXIterations[%d][%d]",c0,c1),"Pull distribution X unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration,100,-10,10);
      pullEventsmeanXIterations[c0][c1] = new TH1D(Form("pullEventsmeanXIterations[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration);
      pullEventssigmaXIterations[c0][c1] = new TH1D(Form("pullEventssigmaXIterations[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration);
 
      
      pullstatEventsXIterations[c0][c1] = new TH2D(Form("pullstatEventsXIterations[%d][%d]",c0,c1),"Pull distribution X unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration,100,-10,10);
      pullstatEventsmeanXIterations[c0][c1] = new TH1D(Form("pullstatEventsmeanXIterations[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);
      pullstatEventssigmaXIterations[c0][c1] = new TH1D(Form("pullstatEventssigmaXIterations[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);

      pullsystEventsXIterations[c0][c1] = new TH2D(Form("pullsystEventsXIterations[%d][%d]",c0,c1),"Pull distribution X unfolding iterations, for true binning (Syst. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration,100,-10,10);
      pullsystEventsmeanXIterations[c0][c1] = new TH1D(Form("pullsystEventsmeanXIterations[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X unfolding iterations, for true binning (Syst. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);
      pullsystEventssigmaXIterations[c0][c1] = new TH1D(Form("pullsystEventssigmaXIterations[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X unfolding iterations, for true binning (Syst. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);      


      pullstatEventsXsystEvents[c0][c1] = new TH2D(Form("pullstatEventsXsystEvents[%d][%d]",c0,c1),"Pull distribution X systematically varied toy experiment, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsSystToys,0,NBinsSystToys,100,-10,10);
      pullstatEventsmeanXsystEvents[c0][c1] = new TH1D(Form("pullstatEventsmeanXsystEvents[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X systematically varied toy experiment, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsSystToys,0,NBinsSystToys);
      pullstatEventssigmaXsystEvents[c0][c1] = new TH1D(Form("pullstatEventssigmaXsystEvents[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X systematically varied toy experiment, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsSystToys,0,NBinsSystToys);
      pullstatEventsXsystEventsXIterations[c0][c1] = new TH3D(Form("pullstatEventsXsystEventsXIterations[%d][%d]",c0,c1),"Pull distribution X systematically varied toy experiment X Number of unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration,NBinsSystToys,0,NBinsSystToys,100,-10,10);

      pullstatEventsXsystEventsXIterationsFinalIteration[c0][c1] = new TH3D(Form("pullstatEventsXsystEventsXIterationsFinalIteration[%d][%d]",c0,c1),"Pull distribution X systematically varied toy experiment X Number of unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration,NBinsSystToys,0,NBinsSystToys,100,-10,10);
    }
  }

  for(int isys=0;isys<NBinsSystToys;isys++){
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	//if(!UsedBinPlots[c0][c1]) continue;
  	BiasXIterations[isys][c0][c1] = new TH1D(Form("BiasXIterations[%d][%d][%d]",isys,c0,c1),Form("For the toy XP %d: Bias distribution as a function of the unfolding iterations, for true binning",isys),NBinsIteration,0,NBinsIteration);
	Chi2XIterations[isys][c0][c1] = new TH1D(Form("Chi2XIterations[%d][%d][%d]",isys,c0,c1),Form("For the toy XP %d: Chi2 distribution as a function of the unfolding iterations, for true binning",isys),NBinsIteration,0,NBinsIteration);
	Chi2DataStatErrorXIterations[isys][c0][c1] = new TH1D(Form("Chi2DataStatErrorXIterations[%d][%d][%d]",isys,c0,c1),Form("For the toy XP %d: Chi2DataStatError distribution as a function of the unfolding iterations, for true binning",isys),NBinsIteration,0,NBinsIteration);

	BiasXIterationsFinalIteration[isys][c0][c1] = new TH1D(Form("BiasXIterationsFinalIteration[%d][%d][%d]",isys,c0,c1),Form("For the toy XP %d: Bias distribution as a function of the unfolding iterations, for true binning",isys),NBinsIteration,0,NBinsIteration);
	Chi2XIterationsFinalIteration[isys][c0][c1] = new TH1D(Form("Chi2XIterationsFinalIteration[%d][%d][%d]",isys,c0,c1),Form("For the toy XP %d: Chi2 distribution as a function of the unfolding iterations, for true binning",isys),NBinsIteration,0,NBinsIteration);	
	Chi2DataStatErrorXIterationsFinalIteration[isys][c0][c1] = new TH1D(Form("Chi2DataStatErrorXIterationsFinalIteration[%d][%d][%d]",isys,c0,c1),Form("For the toy XP %d: Chi2DataStatError distribution as a function of the unfolding iterations, for true binning",isys),NBinsIteration,0,NBinsIteration);
      }
    }
    Chi2DataStatErrorTotalXIterations[isys]= new TH1D(Form("Chi2DataStatErrorTotalXIterations[%d]",isys),"For the toy XP %d: Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
    CorrelatedChi2DataStatErrorTotalXIterations[isys]= new TH1D(Form("CorrelatedChi2DataStatErrorTotalXIterations[%d]",isys),"For the toy XP %d: Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
    Chi2TotalXIterations[isys]= new TH1D(Form("Chi2TotalXIterations[%d]",isys),"For the toy XP %d: Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
    
    Chi2TotalXIterationsFinalIteration[isys]= new TH1D(Form("Chi2TotalXIterationsFinalIteration[%d]",isys),"For the toy XP %d: Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
    Chi2DataStatErrorTotalXIterationsFinalIteration[isys]= new TH1D(Form("Chi2DataStatErrorTotalXIterationsFinalIteration[%d]",isys),"For the toy XP %d: Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[isys]= new TH1D(Form("CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[%d]",isys),"For the toy XP %d: Total chi2/NDF distribution as a function of the unfolding iterations",NBinsIteration,0,NBinsIteration);
  }


  TAxis * aX_UnfoldedDistribution = (TAxis*) UnfoldedDistribution->GetXaxis();
  aX_UnfoldedDistribution->LabelsOption("d");
  //aY_UnfoldedDistribution->LabelsOption("uh");
  TAxis * aX_TrueDistribution = (TAxis*) TrueDistribution->GetXaxis();
  TAxis * aX_BiasDistribution = (TAxis*) BiasDistribution->GetXaxis();
  TAxis * aX_ErrorDistribution = (TAxis*) ErrorDistribution->GetXaxis();
  
  TAxis * aY_UnfoldedDistribution = (TAxis*) UnfoldedDistribution->GetYaxis();
  TAxis * aY_TrueDistribution = (TAxis*) TrueDistribution->GetYaxis();
  TAxis * aY_BiasDistribution = (TAxis*) BiasDistribution->GetYaxis();
  TAxis * aY_ErrorDistribution = (TAxis*) ErrorDistribution->GetYaxis();
  
  for(int ibin=1;ibin<=NBinsTrueMom;ibin++){
    aX_UnfoldedDistribution->SetBinLabel(ibin,Form("%4.1f",BinningTrueMom[ibin-1]));
    aX_TrueDistribution->SetBinLabel(ibin,Form("%4.1fGeV-%4.1fGeV",BinningTrueMom[ibin-1],BinningTrueMom[ibin]));
    aX_BiasDistribution->SetBinLabel(ibin,Form("%4.1fGeV-%4.1fGeV",BinningTrueMom[ibin-1],BinningTrueMom[ibin]));
    aX_ErrorDistribution->SetBinLabel(ibin,Form("%4.1fGeV-%4.1fGeV",BinningTrueMom[ibin-1],BinningTrueMom[ibin]));
  }

  for(int ibin=1;ibin<=NBinsTrueAngle;ibin++){
    aY_TrueDistribution->SetBinLabel(ibin,Form("%2.0f#circ-%2.0f#circ",BinningTrueAngle[ibin-1],BinningTrueAngle[ibin]));
    aY_BiasDistribution->SetBinLabel(ibin,Form("%2.0f#circ-%2.0f#circ",BinningTrueAngle[ibin-1],BinningTrueAngle[ibin]));
    aY_ErrorDistribution->SetBinLabel(ibin,Form("%2.0f#circ-%2.0f#circ",BinningTrueAngle[ibin-1],BinningTrueAngle[ibin]));
  }



  TH1D * hTemp; TF1 * gausTemp = new TF1("gausTemp","gaus",-3,3);
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
  









  
  //###############################################################
  //###############################################################
  //###################### 2. LOAD INPUT FILE #####################
  //###############################################################
  //###############################################################
  TFile * file = new TFile(InputNameMC,"read");
  //TFile * file = new TFile(".root","recreate");
  file->cd();
  TChain * wtree = new TChain("wtree");
  wtree->Add(InputNameMC);
  cout<<"opening file "<<InputNameMC<<" for data input"<<endl;

  int Iterations, Systematics, Statistics;
  double Events[NBinsTrueMom][NBinsTrueAngle];
  double EventsRec[NBinsRecMom][NBinsRecAngle];
  double EventsAll[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double TrueEvents[NBinsTrueMom][NBinsTrueAngle];
  double XSection[NBinsTrueMom][NBinsTrueAngle];
  double TrueXSection[NBinsTrueMom][NBinsTrueAngle];

  TBranch* Br_Iterations = wtree->GetBranch("Iterations");
  Br_Iterations->SetAddress(&Iterations);
  wtree->SetBranchAddress("Iterations",&Iterations);

  TBranch* Br_Systematics = wtree->GetBranch("Systematics");
  Br_Systematics->SetAddress(&Systematics);
  wtree->SetBranchAddress("Systematics",&Systematics);

  TBranch* Br_Statistics = wtree->GetBranch("Statistics");
  Br_Statistics->SetAddress(&Statistics);
  wtree->SetBranchAddress("Statistics",&Statistics);

  TBranch* Br_Events = wtree->GetBranch("Events");
  Br_Events->SetAddress(Events);
  wtree->SetBranchAddress("Events",Events);

  TBranch* Br_TrueEvents = wtree->GetBranch("TrueEvents");
  Br_TrueEvents->SetAddress(TrueEvents);
  wtree->SetBranchAddress("TrueEvents",TrueEvents);

  TBranch* Br_EventsAll = wtree->GetBranch("EventsAll");
  Br_EventsAll->SetAddress(EventsAll);
  wtree->SetBranchAddress("EventsAll",EventsAll);

  TBranch* Br_EventsRec = wtree->GetBranch("EventsRec");
  Br_EventsRec->SetAddress(EventsRec);
  wtree->SetBranchAddress("EventsRec",EventsRec);

  TBranch* Br_XSection = wtree->GetBranch("XSection");
  Br_XSection->SetAddress(XSection);
  wtree->SetBranchAddress("XSection",XSection);

  TBranch* Br_TrueXSection = wtree->GetBranch("TrueXSection");
  Br_TrueXSection->SetAddress(TrueXSection);
  wtree->SetBranchAddress("TrueXSection",TrueXSection);
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################




  

  




  //###############################################################
  //###############################################################
  //################### 3. LOOP OVER INPUT FILE ###################
  //###############################################################
  //###############################################################
  int nevt=wtree->GetEntries();
  int NMaxIterations=0;
  int NStatVariations=0; 

    int ToyID[2] = {-1,-1};//To identify if the loop concerns the same toy (with just one more iteration) or if we are chaning of toys.
    double EventsFinalIteration[NBinsTrueMom][NBinsTrueAngle];
    double XSectionFinalIteration[NBinsTrueMom][NBinsTrueAngle];
     
    for(int ievt=0;ievt<nevt;ievt++){
      wtree->GetEntry(ievt);
      if(ievt%100 == 0) cout<<ievt<<"/"<<nevt<<endl;
#ifdef DEBUG
    cout<<"Number of iteration="<<Iterations<<", number of mom bins="<<NBinsTrueMom<<", angle bins="<<NBinsTrueAngle<<endl;
#endif
    //if(Statistics != 0) continue;

    
      //cout<<"Syst:"<<Systematics<<", Stat:"<<Statistics<<", Current Syst:"<<ToyID[0]<<", Current Stat:"<<ToyID[1]<<endl;
      //TEST IF WE ARE STUDYING THE SAME TOY XP THAN PREVIOUS ITERATION.
    if(ToyID[0] != Systematics /*|| ToyID[1] != Statistics*/){//i.e. NEW TOY
	//We just changed of events, so fill the data final of previous final
	double chi2FinalIteration=0;double chi2Bin = 0; 
	for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	    if(!UsedBinPlots[c0][c1]) continue;
	    if(UseXS){
	      chi2Bin = XSectionFinalIteration[c0][c1]-TrueXSection[c0][c1];
	      if(TrueXSection[c0][c1] != 0) chi2Bin /= TMath::Sqrt(TrueXSection[c0][c1]);
	    }
	    else{
	      //cout<<EventsFinalIteration[c0][c1]<<", "<<TrueEvents[c0][c1]<<endl;
	      chi2Bin = EventsFinalIteration[c0][c1]-TrueEvents[c0][c1];
	      if(TrueEvents[c0][c1] != 0) chi2Bin /= TMath::Sqrt(TrueEvents[c0][c1]);
	    }	      
	    chi2FinalIteration += pow(chi2Bin,2.);
	  }
	}
	//cout<<"Chi2 = "<<chi2FinalIteration<<endl;
	Chi2MCStatErrorFinalIteration->Fill(chi2FinalIteration);
	  
	//Then update the event ID
	ToyID[0] = Systematics;
	ToyID[1] = Statistics;

	//Search for the final iteration step of this toy
	for(int ievt2=ievt+1;ievt2<nevt;ievt2++){
      
	  wtree->GetEntry(ievt2);
	  if(ievt2%100 == 0) cout<<"first event loop:"<<ievt2<<"/"<<nevt<<", stat:"<<Statistics<<endl;

	  if(ToyID[0] == Systematics && Statistics == 0 /*&& ToyID[1] == Statistics*/){

	    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
		if(!UsedBinPlots[c0][c1]) continue;
		EventsFinalIteration[c0][c1] = Events[c0][c1];
		XSectionFinalIteration[c0][c1] = XSection[c0][c1];
	      }
	    }
	  }
	  else break;
	}
	//Set the ievt back to its initial value
	wtree->GetEntry(ievt);
      }
      
      //Then we wish to calculate 2 chi2. The first one is the fake data vs final fake data with the unfolding stat error. The second is the final fake data with the true vs error from true number of events. (chi2 modified)




    if(Statistics > NStatVariations) NStatVariations = Statistics;
    if(Iterations > NMaxIterations) NMaxIterations=Iterations;
    if(Iterations>=NBinsIteration){cout<<"Reset manually the number of iterations in this code to match the input file"<<endl; break;}
    if(Systematics>=NBinsSystToys){cout<<"Reset manually the number of syst toy XP in this code to match the input file"<<endl;break;}

    if(Statistics==0 && Systematics ==0){
      for(int e0=0;e0<NBinsRecMomSignal;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngleSignal;e1++){//loop over effect 1
	  EventReconstructed->SetBinContent(e0,e1,EventsRec[e0][e1]);
	}
      }
    }
    
      for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	  //if(!UsedBinPlots[c0][c1]) continue;
  
#ifdef DEBUG
	  if(Iterations == 0 && Statistics == 0){
	    cout<<setprecision(4)<<"Systematic toy="<<Systematics<<", Mom bin="<<c0<<", ang bin="<<c1<<", Number of Events="<<Events[c0][c1]<<", True events="<<TrueEvents[c0][c1]<<", XS="<<XSection[c0][c1]<<", True XS="<<TrueXSection[c0][c1]<<endl;
	  }
#endif

	  //######### 0. FIRST EVENT: NO STAT OR ERROR VARIATION ##########
	  if(ievt==0){//Only for the non-varied distribution. We do not want to pile-up the information of all the varied distributions!
	    double nEvts=NumberOfEvents->GetBinContent(c0+1,c1+1);
	    for(int e0=0;e0<NBinsRecMomSignal;e0++){//loop over effect 0
	      for(int e1=0;e1<NBinsRecAngleSignal;e1++){//loop over effect 1
		int BinX=e0*(NBinsRecAngleSignal-1)+e1+1;
		int BinY=c0*(NBinsTrueAnglePlots)+c1+1;
		double nEvts_TrueRec=NumberOfEvents_TrueRec->GetBinContent(BinX,BinY);
		nEvts_TrueRec+=EventsAll[c0][c1][e0][e1];
		nEvts+=EventsAll[c0][c1][e0][e1];
		NumberOfEvents_TrueRec->SetBinContent(BinX,BinY,nEvts_TrueRec);
	      }
	    }
	    //cout<<"nevts="<<nEvts<<endl;
	    NumberOfEvents->SetBinContent(c0+1,c1+1,nEvts);
	  }
	  //###############################################################
	  

	  //############### 2. NO SYST OR STAT. VARIATION #################
	  //if(Statistics==0 && Systematics ==0) EventXIterations[c0][c1]->Fill(Iterations,XSection[c0][c1]);
	  if(Statistics==0 && Systematics ==0){
	    EventXIterations[c0][c1]->Fill(Iterations,Events[c0][c1]);
	    //cout<<"("<<c0<<","<<c1<<"), Iterations="<<Iterations<<", Events="<<Events[c0][c1]<<endl;
	  }
	  //###############################################################



	    
	  //######### 2. GENERAL VARIATIONS, WITH SYST AND STAT. ##########
	  double Value=0.; double ValueFinalIteration=0.;
	  if(UseXS){
	    Value=XSection[c0][c1]-TrueXSection[c0][c1];
	    if(TrueXSection[c0][c1]!=0) Value/=TMath::Sqrt(TrueXSection[c0][c1]);
	    ValueFinalIteration=XSection[c0][c1]-XSectionFinalIteration[c0][c1];
	    if(XSectionFinalIteration[c0][c1]!=0) Value/=TMath::Sqrt(XSectionFinalIteration[c0][c1]);
	  }
	  else{
	    Value=Events[c0][c1]-TrueEvents[c0][c1];
	    if(TrueEvents[c0][c1] != 0) Value /= TMath::Sqrt(TrueEvents[c0][c1]);
	    ValueFinalIteration=Events[c0][c1]-EventsFinalIteration[c0][c1];
	    if(EventsFinalIteration[c0][c1]!=0) ValueFinalIteration /= TMath::Sqrt(EventsFinalIteration[c0][c1]);
#ifdef DEBUG
	    //if(Iterations == 0)
	    cout<<"Events="<<Events[c0][c1]<<", true="<<TrueEvents[c0][c1]<<endl;
	    cout<<"Val="<<Value*100<<"%, Value ref="<<ValueFinalIteration*100<<"%"<<endl;
#endif
	  }
	  UnfoldedDistribution->SetBinContent(c0+1,c1+1,Events[c0][c1]);
	  UnfoldedFinalIterationDistribution->SetBinContent(c0+1,c1+1,EventsFinalIteration[c0][c1]);
	  TrueDistribution->SetBinContent(c0+1,c1+1,TrueEvents[c0][c1]);
	  pullEventsXIterations[c0][c1]->Fill(Iterations,Value);
	  pullstatEventsXsystEvents[c0][c1]->Fill(Systematics,Value);
	  pullEvents[c0][c1]->Fill(Value);
	  pullstatEventsXsystEventsXIterations[c0][c1]->Fill(Iterations,Systematics,Value);
	  pullstatEventsXsystEventsXIterationsFinalIteration[c0][c1]->Fill(Iterations,Systematics,ValueFinalIteration);
	  //###############################################################
	  


	  
	  //######## 2. SYSTEMATIC ONLY VARIATIONS, NO STAT VAR  ##########
	  if(Statistics==0){
	    pullsystEventsXIterations[c0][c1]->Fill(Iterations,Value);
	  }
	  //###############################################################
	



	  //########### 3. STAT ONLY VARIATIONS, NO SYST VAR  #############
	  if(Systematics==0){//no syst variation
	    pullstatEventsXIterations[c0][c1]->Fill(Iterations,Value);
	    for(int d0=0;d0<NBinsTrueMomPlots;d0++){//loop over dause 0
	      for(int d1=0;d1<NBinsTrueAnglePlots;d1++){//loop over cause 1
		//if(!UsedBinPlots[d0][d1]) continue;
		int BinX=c0*(NBinsTrueAnglePlots)+c1+1;
		int BinY=d0*(NBinsTrueAnglePlots)+d1+1;

		double Cova,CovaFinalIteration;
		if(UseXS){
		  Cova=(XSection[c0][c1]-TrueXSection[c0][c1])*(XSection[d0][d1]-TrueXSection[d0][d1]);
		  CovaFinalIteration=(XSection[c0][c1]-XSectionFinalIteration[c0][c1])*(XSection[d0][d1]-XSectionFinalIteration[d0][d1]);
		}
		else{
		  Cova=(Events[c0][c1]-TrueEvents[c0][c1])*(Events[d0][d1]-TrueEvents[d0][d1]);
		  CovaFinalIteration=(Events[c0][c1]-EventsFinalIteration[c0][c1])*(Events[d0][d1]-EventsFinalIteration[d0][d1]);
		} 
		if(Iterations == 0) CovarianceStat_Iteration1->SetBinContent(BinX,BinY,CovarianceStat_Iteration1->GetBinContent(BinX,BinY)+Cova); 
		else if(Iterations == 2) CovarianceStat_Iteration3->SetBinContent(BinX,BinY,CovarianceStat_Iteration3->GetBinContent(BinX,BinY)+Cova); 
		CovarianceStat->SetBinContent(BinX,BinY,CovarianceStat->GetBinContent(BinX,BinY)+Cova);
		
		CovarianceXIterations[Iterations](BinX-1,BinY-1) = CovarianceXIterations[Iterations](BinX-1,BinY-1)+Cova;
		CovarianceInvertXIterations[Iterations](BinX-1,BinY-1) = CovarianceInvertXIterations[Iterations](BinX-1,BinY-1)+Cova;
		CovarianceXIterationsFinalIteration[Iterations](BinX-1,BinY-1) = CovarianceXIterationsFinalIteration[Iterations](BinX-1,BinY-1)+CovaFinalIteration;
		CovarianceInvertXIterationsFinalIteration[Iterations](BinX-1,BinY-1) = CovarianceInvertXIterationsFinalIteration[Iterations](BinX-1,BinY-1)+CovaFinalIteration;

		int nvar=NVariationsStat->GetBinContent(BinX,BinY);
		nvar++;
		NVariationsStat->SetBinContent(BinX,BinY,nvar);
		
	      }
	    }
	  }
	  //###############################################################

	}
      }
    }
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
  //###############################################################
    //General->cd();
    cout<<"number of stat variations="<<NStatVariations<<endl;
    for(int it=0;it<NBinsIteration;it++){
      for(int a=0;a<CovarianceXIterations[it].GetNrows();a++){
	for(int b=0;b<CovarianceXIterations[it].GetNcols();b++){
	  CovarianceInvertXIterations[it](a,b) = CovarianceInvertXIterations[it](a,b) / NStatVariations;
	  CovarianceXIterations[it](a,b) = CovarianceXIterations[it](a,b) / NStatVariations;
	  CovarianceInvertXIterationsFinalIteration[it](a,b) = CovarianceInvertXIterationsFinalIteration[it](a,b) / NStatVariations;
	  CovarianceXIterationsFinalIteration[it](a,b) = CovarianceXIterationsFinalIteration[it](a,b) / NStatVariations;
	}
      }
      CovarianceInvertXIterations[it].Invert();
      CovarianceInvertXIterationsFinalIteration[it].Invert();
      //CovarianceXIterations[it].Write(Form("CovarianceXIterations_%d",it));
      //CovarianceInvertXIterations[it].Write(Form("CovarianceInvertXIterations_%d",it));  
    } 





    

  //###############################################################
  //###############################################################
  //###################### 4. FILL HISTOGRAMS #####################
  //###############################################################
  //###############################################################
#ifdef PLOTS	      
    cout<<"Write plots"<<endl;

    //################## 0. BASIC HISTOGRAMS ###################
    General->cd();
    TCanvas * canTrue = new TCanvas("cNumberOfEventsTrue","");
    TrueDistribution->GetXaxis()->SetTitle("Momentum bin");
    TrueDistribution->GetYaxis()->SetTitle("Momentum angle");
    TrueDistribution->Draw("colztext");
    canTrue->Write();
    
    EventReconstructed->GetXaxis()->SetTitle("d_{#mu} (cm)");
    EventReconstructed->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
    EventReconstructed->Write();
    TH1D * EventReconstructedMom = (TH1D*) EventReconstructed->ProjectionX("EventReconstructedMom",0,EventReconstructed->GetNbinsY());EventReconstructedMom->Write();
    TH1D * EventReconstructedAngle = (TH1D*) EventReconstructed->ProjectionY("EventReconstructedAngle",0,EventReconstructed->GetNbinsX());EventReconstructedAngle->Write();
    
    //############# 1. NO STAT. OR SYST. VARIATION #############
    NoStatOrSystVariation->cd();
    TCanvas * canTrueXUnfolded = new TCanvas("cNumberOfEvents_UnfoldedAndTrue","");
    UnfoldedDistribution->GetXaxis()->SetTitle("Momentum bin");
    UnfoldedDistribution->GetYaxis()->SetTitle("Momentum angle");
    TrueDistribution->GetXaxis()->SetTitle("Momentum bin");
    TrueDistribution->GetYaxis()->SetTitle("Momentum angle");
  
    canTrueXUnfolded->Divide(1,2);
    canTrueXUnfolded->cd(1);
    UnfoldedDistribution->Draw("colztext");
    canTrueXUnfolded->cd(2);
    TrueDistribution->Draw("colztext");
    canTrueXUnfolded->Write();

    TLine * lBins[NBinsTrueMomPlots][2];
    TLatex * tBins[NBinsTrueMomPlots][2];
    cout.precision(2);
    cout.setf(ios::fixed);
    
    for(int imom=0;imom<NBinsTrueMomPlots;imom++){
      cout<<"imom="<<imom<<endl;
      lBins[imom][0] = new TLine(imom*(NBinsTrueAnglePlots),0,imom*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots));
      lBins[imom][1] = new TLine(0,imom*(NBinsTrueAnglePlots),(NBinsTrueMomPlots)*(NBinsTrueAnglePlots),imom*(NBinsTrueAnglePlots));
      
      //tBins[imom][0] = new TLatex(imom*(NBinsTrueAnglePlots)+NBinsTrueAngle/4,-1,Form("Bin%d p_{#mu}",imom+1));
      //tBins[imom][1] = new TLatex(-1,imom*(NBinsTrueAnglePlots)+NBinsTrueAngle/4,Form("Bin%d p_{#mu}",imom+1));
      //tBins[imom][0]->SetTextSize(0.04);
      //tBins[imom][1]->SetTextSize(0.04);

      tBins[imom][0] = new TLatex(imom*(NBinsTrueAnglePlots)+NBinsTrueAngle/7,-1,Form("%4.1f-%4.1fGeV",BinningTrueMom[imom],BinningTrueMom[imom+1]));      
      tBins[imom][1] = new TLatex(-1,imom*(NBinsTrueAnglePlots)+NBinsTrueAngle/15,Form("%4.1f-%4.1fGeV",BinningTrueMom[imom],BinningTrueMom[imom+1]));
      tBins[imom][0]->SetTextSize(0.03);
      tBins[imom][1]->SetTextSize(0.03);

      tBins[imom][1]->SetTextAngle(90);
      
      lBins[imom][0]->SetLineColor(1);
      lBins[imom][0]->SetLineWidth(4);
      lBins[imom][0]->SetLineStyle(2);
      
      lBins[imom][1]->SetLineColor(1);
      lBins[imom][1]->SetLineWidth(4);
      lBins[imom][1]->SetLineStyle(2);
      
      tBins[imom][0]->Draw("same");
      tBins[imom][1]->Draw("same");
      lBins[imom][0]->Draw("same");
      lBins[imom][1]->Draw("same");
    }


    
    for(int iit=0;iit<EventTotalXIterations->GetNbinsX();iit++){
      double evtTotal = 0;
      double chi2Total = 0;
      int NBinsTotal = 0;
      for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	  if(!UsedBinPlots[c0][c1]) continue;
	  double evtUnfolded = EventXIterations[c0][c1]->GetBinContent(iit+1);
	  double evtTrue = TrueDistribution->GetBinContent(c0+1,c1+1);
	  double error = TMath::Sqrt(evtTrue);
	  double chi2 = pow(evtUnfolded-evtTrue,2.);
	  if( error != 0) chi2 /= pow(error,2.);
    //TEMP: TO FIX WITH MODULABLE NUMBER OF BINS: C0 should start at NStartBin and change c0 -1 -> c0 below
	  if(c0>0){
	    Chi2MCStatErrorXIterations[c0][c1]->SetBinContent(iit+1,chi2);
	    chi2Total += chi2;
	    NBinsTotal ++;
	  }	 
	  evtTotal += evtUnfolded;
	}
      }
      EventTotalXIterations->SetBinContent(iit+1,evtTotal);
      chi2Total /= NBinsTotal;
      Chi2MCStatErrorTotalXIterations->SetBinContent(iit+1,chi2Total);
    }
    EventTotalXIterations->Write();
    Chi2MCStatErrorTotalXIterations->Write();


    TCanvas * cChi2MCStatError = new TCanvas("cChi2MCStatError","");
    TLegend * lModified = new TLegend(0.75,0.2,0.99,0.8);
    lModified->SetFillColor(0.);
    lModified->SetLineColor(0.);    
    Chi2MCStatErrorTotalXIterations->SetLineWidth(3);
    Chi2MCStatErrorTotalXIterations->SetLineStyle(2);
    Chi2MCStatErrorTotalXIterations->GetXaxis()->SetTitle("Number of iterations");
    Chi2MCStatErrorTotalXIterations->GetYaxis()->SetTitle("#chi^{2} / NDF");
    Chi2MCStatErrorTotalXIterations->GetYaxis()->SetRangeUser(1e-10,1e4);
    Chi2MCStatErrorTotalXIterations->Draw();
    lModified->AddEntry(Chi2MCStatErrorTotalXIterations,"Total");
    TLine * line1 = new TLine(0.,1.,Chi2MCStatErrorTotalXIterations->GetXaxis()->GetXmax(),1.);
    line1->SetLineColor(1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(10);
    line1->Draw("same");
    
    //TEMP: TO FIX WITH MODULABLE NUMBER OF BINS: C0 should start at NStartBin and change c0 -1 -> c0 below
    for(int c0=1;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	
	if(!UsedBinPlots[c0][c1]) continue;
  	Chi2MCStatErrorXIterations[c0][c1]->SetLineWidth(2);
	Chi2MCStatErrorXIterations[c0][c1]->SetLineColor(color[(c0-FirstBinMomPlots)*NBinsTrueAnglePlots+(c1-FirstBinAnglePlots)]);
	lModified->AddEntry(Chi2MCStatErrorXIterations[c0][c1] ,Form("p_{#mu} %d, #theta_{#mu} %d",c0,c1));
	Chi2MCStatErrorXIterations[c0][c1]->Draw("same");
	EventXIterations[c0][c1]->Write();
	Chi2MCStatErrorXIterations[c0][c1]->Write();
      }
    }
    lModified->Draw("same");
    cChi2MCStatError->SetGridx();
    cChi2MCStatError->SetGridy();
    cChi2MCStatError->SetLogy();
    cChi2MCStatError->Write();

  
    //############## 2. STAT. + SYST. VARIATIONS ###############
    StatAndSystVariations->cd();
    
    TF1 * fPoisson[NBinsTrueMom][NBinsTrueAngle];
    //General variations: both systematic and statictal uncertainties are varied
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	if(!UsedBinPlots[c0][c1]) continue;
	fPoisson[c0][c1] = new TF1(Form("fPoisson[%d][%d]",c0,c1),"[0]*(TMath::Poisson((x*[1]+[1]),[1]))",-10,10);
	fPoisson[c0][c1]->SetParameter(1,TrueDistribution->GetBinContent(c0+1,c1+1));
	fPoisson[c0][c1]->SetParameter(0,pullEvents[c0][c1]->Integral());
	//Draw the pull distribution for all the variations in the input file
	canEvents[c0][c1] = new TCanvas(Form("cPullDistribution_%d_%d",c0,c1),Form("True distribution of events after unfolding in the bin %d_%d",c0,c1));
	pullEvents[c0][c1]->GetXaxis()->SetTitle("Pull");
	pullEvents[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	pullEvents[c0][c1]->Draw();
	pullEvents[c0][c1]->Fit(Form("fPoisson[%d][%d]",c0,c1),"RQ");
	pullEvents[c0][c1]->Fit("gaus","Q");
	fPoisson[c0][c1]->SetLineColor(1);
	fPoisson[c0][c1]->Draw("same");
	canEvents[c0][c1]->Write();
	//Draw the pull distribution with iterations & behaviour of the mean/error of the pull, for all the variations in the input file
	
	pullEventsXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding iterations");
	pullEventsXIterations[c0][c1]->GetYaxis()->SetTitle("Pull");
	
	
	for(int ibinx=1;ibinx<=pullEventsXIterations[c0][c1]->GetNbinsX();ibinx++){
	  gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	  if(Fit){
	    pullEventsXIterations[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");     
	    pullEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(1));
	    pullEventsmeanXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(1));
	    pullEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(2));
	    pullEventssigmaXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(2));
	    
	    if(ibinx==pullEventsXIterations[c0][c1]->GetNbinsX()){//only for the last bin of the number of iteration. In the case that there is no variations on the number of iterations, this is equal to the number of iteration set.
	      BiasDistribution->SetBinContent(c0+1,c1+1,gausTemp->GetParameter(1));
	      ErrorDistribution->SetBinContent(c0+1,c1+1,gausTemp->GetParameter(2));
	    }
	  }
	  else{
	    TH1D * Slice = (TH1D*) pullEventsXIterations[c0][c1]->ProjectionY("Slice",ibinx,ibinx);
	    pullEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,Slice->GetMean());
	    pullEventsmeanXIterations[c0][c1]->SetBinError(ibinx,Slice->GetMeanError());
	    pullEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,Slice->GetRMS());
	    pullEventssigmaXIterations[c0][c1]->SetBinError(ibinx,Slice->GetRMSError());
	    
	    if(ibinx==pullEventsXIterations[c0][c1]->GetNbinsX()){//only for the last bin of the number of iteration. In the case that there is no variations on the number of iterations, this is equal to the number of iteration set.
	      BiasDistribution->SetBinContent(c0+1,c1+1,Slice->GetMean());
	      ErrorDistribution->SetBinContent(c0+1,c1+1,Slice->GetRMS());
	    }
	  } 
	  
	}
	pullEventsmeanXIterations[c0][c1]->SetLineColor(1);
	pullEventssigmaXIterations[c0][c1]->SetLineColor(kGray);
	pullEventsmeanXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding Iterations");
	//pullEventsmeanXIterations[c0][c1]->Draw("E1");
	//pullEventssigmaXIterations[c0][c1]->Draw("E1same");
	pullEventsXIterations[c0][c1]->Write();
	
      }
    }
    cout<<"Start to watch toy XP / toy XP"<<endl;
    //Want to build a plot chi2 statistics vs iterations for each syst toy XP. First, obtain the Iteration vs stat pull for each syst toy XP.
    TH2D * TempstatPullXIterations;
    double Chi2Total[NBinsSystToys][NBinsIteration]={{0.}};
    double Chi2NominalTotal[NBinsSystToys][NBinsIteration]={{0.}};
    double CorrelatedChi2NominalTotal[NBinsSystToys][NBinsIteration]={{0.}};
    int NDFTotal[NBinsSystToys][NBinsIteration]={{0.}};

    //ProjectAndEvaluateChi2(pullstatEventsXsystEventsXIterations,)
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	if(!UsedBinPlots[c0][c1]) continue;
	cout<<"Bin "<<c0<<", "<<c1<<endl;
	for(int isys=1;isys<=pullstatEventsXsystEventsXIterations[c0][c1]->GetNbinsY();isys++){
	  cout<<"Toy XP #"<<isys<<endl;
	  pullstatEventsXsystEventsXIterations[c0][c1]->GetYaxis()->SetRange(isys,isys);
	  TempstatPullXIterations = (TH2D*) pullstatEventsXsystEventsXIterations[c0][c1]->Project3D("zx");
	  TempstatPullXIterations->Write(Form("StatPullXIterationsForToy[%d][%d][%d]",isys-1,c0,c1));
	  //then, for each iteration, wants to fit the pull:

	  for(int ibinx=1;ibinx<=TempstatPullXIterations->GetNbinsX();ibinx++){
	    double Mean;
	    double RMS;
	    TH1D * Slice;

	    if(Fit){
	      gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	      TempstatPullXIterations->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	      Mean = gausTemp->GetParameter(1);
	      RMS = gausTemp->GetParameter(2);
	    }
	    else{
	      Slice = (TH1D*) TempstatPullXIterations->ProjectionY("Slice",ibinx,ibinx);
	      Mean = Slice->GetMean();
	      RMS = Slice->GetRMS();
	    }
	    //Then, put this into a chi2
	    double Chi2=pow(Mean,2.);
	    if(RMS !=0 ) Chi2 /= pow(RMS,2.);
	    BiasXIterations[isys-1][c0][c1]->SetBinContent(ibinx,Mean);
	    Chi2XIterations[isys-1][c0][c1]->SetBinContent(ibinx,Chi2);
	    Chi2Total[isys-1][ibinx-1] += Chi2;
	    NDFTotal[isys-1][ibinx-1] ++;
	    double Chi2Nominal=pow(EventXIterations[c0][c1]->GetBinContent(ibinx) - TrueDistribution->GetBinContent(c0+1,c1+1),2.);
	    double ErrorNominal=TMath::Sqrt(TrueDistribution->GetBinContent(c0+1,c1+1))*RMS;//Multiply the spread of the pull by the number of event to have the spread in number of events
	    if(ErrorNominal) Chi2Nominal /= pow(ErrorNominal,2.);
	    //cout<<EventXIterations[c0][c1]->GetBinContent(ibinx)<<", "<<TrueDistribution->GetBinContent(c0+1,c1+1)<<", chi2="<<Chi2Nominal<<endl;
	    Chi2DataStatErrorXIterations[isys-1][c0][c1]->SetBinContent(ibinx,Chi2Nominal);
	    Chi2NominalTotal[isys-1][ibinx-1] += Chi2Nominal;

	    for(int d0=0;d0<NBinsTrueMomPlots;d0++){//loop over cause 0
	      for(int d1=0;d1<NBinsTrueAnglePlots;d1++){//loop over cause 1
		if(!UsedBinPlots[d0][d1]) continue;
		double Chi2Correlated=(EventXIterations[c0][c1]->GetBinContent(ibinx) - TrueDistribution->GetBinContent(c0+1,c1+1))*(EventXIterations[d0][d1]->GetBinContent(ibinx) - TrueDistribution->GetBinContent(d0+1,d1+1));
		int BinX=c0*(NBinsTrueAnglePlots)+c1;
		int BinY=d0*(NBinsTrueAnglePlots)+d1;
		double ErrorInvertSquaredCorrelated=CovarianceInvertXIterations[ibinx-1](BinX,BinY);
		CorrelatedChi2NominalTotal[isys-1][ibinx-1] += Chi2Correlated*ErrorInvertSquaredCorrelated;
		//if(c0 == d0 && c1 == d1) cout<<setprecision(6)<<TMath::Sqrt(CovarianceXIterations[ibinx-1](BinX,BinY))<<", "<<Chi2Correlated<<", "<<1./TMath::Sqrt(ErrorInvertSquaredCorrelated)<<endl;
	      }
	    }
	    //cout<<"Iteration #"<<ibinx<<", chi2="<<Chi2NominalTotal[isys-1][ibinx-1]<<", correlated chi2="<<CorrelatedChi2NominalTotal[isys-1][ibinx-1]<<endl;
	  }
	  BiasXIterations[isys-1][c0][c1]->Write();
	  Chi2XIterations[isys-1][c0][c1]->Write();
	}	
      }	
    }

    TCanvas * cChi2FittedDataStatError = new TCanvas("cChi2FittedDataStatError","");
    for(int isys=1;isys<=pullstatEventsXsystEventsXIterations[0][0]->GetNbinsY();isys++){
      for(int ibinx=1;ibinx<=pullstatEventsXsystEventsXIterations[0][0]->GetNbinsX();ibinx++){
	cout<<ibinx<<endl;
	//if(!UsedBinPlotsMom[ibinx-1]) continue;
	cout<<ibinx<<" vs "<<NBinsTrueMomPlots<<", chi2="<<Chi2Total[isys-1][ibinx-1]<<", NDF="<<NDFTotal[isys-1][ibinx-1]<<endl;
	Chi2Total[isys-1][ibinx-1] /= NDFTotal[isys-1][ibinx-1];
	Chi2TotalXIterations[isys-1]->SetBinContent(ibinx,Chi2Total[isys-1][ibinx-1]);
      }
      Chi2TotalXIterations[isys-1]->Write();
  
      //Chi2TotalXIterations->Draw();
      //if(isys==1) l->AddEntry(Chi2TotalXIterations[isys-1],"Total");
    
      Chi2TotalXIterations[isys-1]->Draw();
      Chi2TotalXIterations[isys-1]->GetXaxis()->SetTitle("Number of iterations");
      Chi2TotalXIterations[isys-1]->GetYaxis()->SetTitle("#chi^{2} / NDF");
      Chi2TotalXIterations[isys-1]->SetLineWidth(3);
      Chi2TotalXIterations[isys-1]->SetLineStyle(2);
      Chi2TotalXIterations[isys-1]->GetYaxis()->SetRangeUser(1e-5,100.);
      for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	
	for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	  if(!UsedBinPlots[c0][c1]) continue;
	  Chi2XIterations[isys-1][c0][c1]->SetLineWidth(2);
	  Chi2XIterations[isys-1][c0][c1]->SetLineColor(color[(c0-FirstBinMomPlots)*NBinsTrueAnglePlots+(c1-FirstBinAnglePlots)]);
	  //if(isys==1 && c0==0) l->AddEntry(Chi2XIterations[isys-1][c0][c1] , Form("#theta_{#mu} %d",c1));
	  Chi2XIterations[isys-1][c0][c1]->Draw("same");
	  //int BinY=c0*(NBinsTrueAnglePlots)+c1+1;
	  //Chi2XIterations[c0][c1]->SetLineColor(BinY);
	  //Chi2XIterations[c0][c1]->Draw("same");
	  //l->AddEntry( Chi2XIterations[c0][c1] , Form("Mom %d, Angle %d",c0,c1));
	}
	//l->Draw("same");
	//cChi2[isys-1][c0]->SetGridx();
	//cChi2[isys-1][c0]->SetGridy();
	//cChi2[isys-1][c0]->SetLogy();
       //cChi2[isys-1][c0]->Write();
      }
    }
    line1->Draw("same");
    lModified->Draw("same");
    cChi2FittedDataStatError->SetGridx();
    cChi2FittedDataStatError->SetGridy();
    cChi2FittedDataStatError->SetLogy();
    cChi2FittedDataStatError->Write();

  TCanvas * canBias = new TCanvas("cPull_Mean");
  gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(kFALSE);
  BiasDistribution->Draw("colztext");
  canBias->Write();
  TCanvas * canError = new TCanvas("cPull_Error");
  gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(kFALSE);
  ErrorDistribution->Draw("colztext");
  canError->Write();

  
  TCanvas * cChi2DataStatError = new TCanvas("cChi2DataStatError");
  //TCanvas * cChi2DataStatError[NBinsSystToys][NBinsTrueMomPlots];
  for(int isys=1;isys<=pullstatEventsXsystEventsXIterations[0][0]->GetNbinsY();isys++){
    for(int ibinx=1;ibinx<=pullstatEventsXsystEventsXIterations[0][0]->GetNbinsX();ibinx++){
      Chi2NominalTotal[isys-1][ibinx-1] /= NDFTotal[isys-1][ibinx-1];
      Chi2DataStatErrorTotalXIterations[isys-1]->SetBinContent(ibinx,Chi2NominalTotal[isys-1][ibinx-1]);
      CorrelatedChi2NominalTotal[isys-1][ibinx-1] /= NDFTotal[isys-1][ibinx-1];
      CorrelatedChi2DataStatErrorTotalXIterations[isys-1]->SetBinContent(ibinx,CorrelatedChi2NominalTotal[isys-1][ibinx-1]);
      //cout<<"Iteration #"<<ibinx<<", Chi2="<<Chi2NominalTotal[isys-1][ibinx-1]<<endl;
    }
    Chi2DataStatErrorTotalXIterations[isys-1]->GetXaxis()->SetTitle("Number of iterations");
    Chi2DataStatErrorTotalXIterations[isys-1]->GetYaxis()->SetTitle("#chi^{2} / NDF");
    Chi2DataStatErrorTotalXIterations[isys-1]->Write();
    Chi2DataStatErrorTotalXIterations[isys-1]->Draw();
    Chi2DataStatErrorTotalXIterations[isys-1]->SetLineWidth(3);
    Chi2DataStatErrorTotalXIterations[isys-1]->SetLineStyle(2);
    Chi2DataStatErrorTotalXIterations[isys-1]->GetYaxis()->SetRangeUser(1e-5,100.);

    CorrelatedChi2DataStatErrorTotalXIterations[isys-1]->GetXaxis()->SetTitle("Number of iterations");
    CorrelatedChi2DataStatErrorTotalXIterations[isys-1]->GetYaxis()->SetTitle("#chi^{2} / NDF");
    CorrelatedChi2DataStatErrorTotalXIterations[isys-1]->Write();
    
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      //cChi2DataStatError[isys-1][c0] = new TCanvas(Form("cChi2DataStatError for toy XP %d and mom %d",isys-1,c0),"");
      
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	if(!UsedBinPlots[c0][c1]) continue;
	Chi2DataStatErrorXIterations[isys-1][c0][c1]->SetLineWidth(2);
	Chi2DataStatErrorXIterations[isys-1][c0][c1]->SetLineColor(color[(c0-FirstBinMomPlots)*NBinsTrueAnglePlots+(c1-FirstBinAnglePlots)]);
	Chi2DataStatErrorXIterations[isys-1][c0][c1]->Draw("same");
	//int BinY=c0*(NBinsTrueAnglePlots)+c1+1;
	//Chi2DataStatErrorXIterations[c0][c1]->SetLineColor(BinY);
	//Chi2DataStatErrorXIterations[c0][c1]->Draw("same");
	//l->AddEntry( Chi2DataStatErrorXIterations[c0][c1] , Form("Mom %d, Angle %d",c0,c1));
      }
    }
  }
  line1->Draw("same");
  lModified->Draw("same");
  cChi2DataStatError->SetGridx();
  cChi2DataStatError->SetGridy();
  cChi2DataStatError->SetLogy();
  cChi2DataStatError->Write();


  TH2D * TempstatPullXIterationsFinalIteration;
  double Chi2TotalFinalIteration[NBinsSystToys][NBinsIteration]={{0.}};
  double Chi2NominalTotalFinalIteration[NBinsSystToys][NBinsIteration]={{0.}};
  double CorrelatedChi2NominalTotalFinalIteration[NBinsSystToys][NBinsIteration]={{0.}};
  int NDFTotalFinalIteration[NBinsSystToys][NBinsIteration]={{0.}};
 
    //ProjectAndEvaluateChi2(pullstatEventsXsystEventsXIterationsFinalIteration,)
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	if(!UsedBinPlots[c0][c1]) continue;
	cout<<"Bin "<<c0<<", "<<c1<<endl;
	
	for(int isys=1;isys<=pullstatEventsXsystEventsXIterationsFinalIteration[c0][c1]->GetNbinsY();isys++){
	  cout<<"Toy XP #"<<isys<<endl;
	  
	  pullstatEventsXsystEventsXIterationsFinalIteration[c0][c1]->GetYaxis()->SetRange(isys,isys);
	  TempstatPullXIterationsFinalIteration = (TH2D*) pullstatEventsXsystEventsXIterationsFinalIteration[c0][c1]->Project3D("zx");
	  TempstatPullXIterationsFinalIteration->Write(Form("StatPullXIterationsFinalIterationForToy[%d][%d][%d]",isys-1,c0,c1));
	  //then, for each iteration, wants to fit the pull:
	  for(int ibinx=1;ibinx<=TempstatPullXIterationsFinalIteration->GetNbinsX();ibinx++){
	    double Mean;
	    double RMS;
	    TH1D * Slice;

	    if(Fit){
	      gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	      TempstatPullXIterationsFinalIteration->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	      Mean = gausTemp->GetParameter(1);
	      RMS = gausTemp->GetParameter(2);
	    }
	    else{
	      Slice = (TH1D*) TempstatPullXIterationsFinalIteration->ProjectionY("Slice",ibinx,ibinx);
	      Mean = Slice->GetMean();
	      RMS = Slice->GetRMS();
	    }

	    //Then, put this into a chi2
	    double Chi2=pow(Mean,2.);
	    if(RMS !=0 ) Chi2 /= pow(RMS,2.);
	    //BiasXIterationsFinalIteration[isys-1][c0][c1]->SetBinContent(ibinx,Mean);
	    Chi2XIterationsFinalIteration[isys-1][c0][c1]->SetBinContent(ibinx,Chi2);
	    Chi2TotalFinalIteration[isys-1][ibinx-1] += Chi2;
	    NDFTotalFinalIteration[isys-1][ibinx-1] ++;

	    double Chi2Nominal=pow(EventXIterations[c0][c1]->GetBinContent(ibinx) - UnfoldedFinalIterationDistribution->GetBinContent(c0+1,c1+1),2.);
	    double ErrorNominal=TMath::Sqrt(UnfoldedFinalIterationDistribution->GetBinContent(c0+1,c1+1))*RMS;//Multiply the spread of the pull by the number of event to have the spread in number of events
	    if(ErrorNominal) Chi2Nominal /= pow(ErrorNominal,2.);
	    //cout<<EventXIterations[c0][c1]->GetBinContent(ibinx)<<", "<<TrueDistribution->GetBinContent(c0+1,c1+1)<<", chi2="<<Chi2Nominal<<endl;
	    Chi2DataStatErrorXIterationsFinalIteration[isys-1][c0][c1]->SetBinContent(ibinx,Chi2Nominal);
	    Chi2NominalTotalFinalIteration[isys-1][ibinx-1] += Chi2Nominal;

	    for(int d0=0;d0<NBinsTrueMomPlots;d0++){//loop over cause 0
	      for(int d1=0;d1<NBinsTrueAnglePlots;d1++){//loop over cause 1
		if(!UsedBinPlots[d0][d1]) continue;
		double Chi2Correlated=(EventXIterations[c0][c1]->GetBinContent(ibinx) - UnfoldedFinalIterationDistribution->GetBinContent(c0+1,c1+1))*(EventXIterations[d0][d1]->GetBinContent(ibinx) - UnfoldedFinalIterationDistribution->GetBinContent(d0+1,d1+1));
		int BinX=c0*(NBinsTrueAnglePlots)+c1;
		int BinY=d0*(NBinsTrueAnglePlots)+d1;
		double ErrorInvertSquaredCorrelated=CovarianceInvertXIterationsFinalIteration[ibinx-1](BinX,BinY);
		CorrelatedChi2NominalTotalFinalIteration[isys-1][ibinx-1] += Chi2Correlated*ErrorInvertSquaredCorrelated;
		//if(c0 == d0 && c1 == d1) cout<<setprecision(6)<<TMath::Sqrt(CovarianceXIterations[ibinx-1](BinX,BinY))<<", "<<Chi2Correlated<<", "<<1./TMath::Sqrt(ErrorInvertSquaredCorrelated)<<endl;
	      }
	    }

	  }
	  //BiasXIterationsFinalIteration[isys-1][c0][c1]->Write();
	  Chi2XIterationsFinalIteration[isys-1][c0][c1]->Write();
	}	
      }	
    }
   
    TCanvas * cChi2FinalIteration = new TCanvas("cChi2FinalIteration");
    //TCanvas * cChi2FinalIteration[NBinsSystToys][NBinsTrueMomPlots];
    for(int isys=1;isys<=pullstatEventsXsystEventsXIterationsFinalIteration[0][0]->GetNbinsY();isys++){
      for(int ibinx=1;ibinx<=pullstatEventsXsystEventsXIterationsFinalIteration[0][0]->GetNbinsX();ibinx++){
	Chi2TotalFinalIteration[isys-1][ibinx-1] /= NDFTotalFinalIteration[isys-1][ibinx-1];
	Chi2TotalXIterationsFinalIteration[isys-1]->SetBinContent(ibinx,Chi2TotalFinalIteration[isys-1][ibinx-1]);
      }
      //Chi2TotalXIterationsFinalIteration[isys-1]->Write();
      Chi2TotalXIterationsFinalIteration[isys-1]->GetXaxis()->SetTitle("Number of iterations");
      Chi2TotalXIterationsFinalIteration[isys-1]->GetYaxis()->SetTitle("#chi^{2} / NDF");
      Chi2TotalXIterationsFinalIteration[isys-1]->Draw();
      Chi2TotalXIterationsFinalIteration[isys-1]->SetLineWidth(3);
      Chi2TotalXIterationsFinalIteration[isys-1]->SetLineStyle(2);
      Chi2TotalXIterationsFinalIteration[isys-1]->GetYaxis()->SetRangeUser(1e-5,100.);
      
      for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	//cChi2FinalIteration[isys-1][c0] = new TCanvas(Form("cChi2FinalIteration for toy XP %d and mom %d",isys-1,c0),"");

	for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	  if(!UsedBinPlots[c0][c1]) continue;
	  Chi2XIterationsFinalIteration[isys-1][c0][c1]->SetLineWidth(2);
	  Chi2XIterationsFinalIteration[isys-1][c0][c1]->SetLineColor(color[(c0-FirstBinMomPlots)*NBinsTrueAnglePlots+(c1-FirstBinAnglePlots)]);
	  Chi2XIterationsFinalIteration[isys-1][c0][c1]->Draw("same");
	  //int BinY=c0*(NBinsTrueAnglePlots)+c1+1;->SetLineColor((c0-FirstBinMomPlots)*NBinsTrueAnglePlots+(c1-FirstBinAnglePlots
	  //Chi2XIterationsFinalIteration[c0][c1]->SetLineColor(BinY);
	  //Chi2XIterationsFinalIteration[c0][c1]->Draw("same");
	  //l->AddEntry( Chi2XIterationsFinalIteration[c0][c1] , Form("Mom %d, Angle %d",c0,c1));
	}
      }
    }
    lModified->Draw("same");
    line1->Draw("same");
    cChi2FinalIteration->SetGridx();
    cChi2FinalIteration->SetGridy();
    cChi2FinalIteration->SetLogy();
    cChi2FinalIteration->Write();
    Chi2MCStatErrorFinalIteration->Write();

    TCanvas * cChi2DataStatErrorFinalIteration = new TCanvas("cChi2DataStatErrorFinalIteration");
    //TCanvas * cChi2DataStatErrorFinalIteration[NBinsSystToys][NBinsTrueMomPlots];
    for(int isys=1;isys<=pullstatEventsXsystEventsXIterationsFinalIteration[0][0]->GetNbinsY();isys++){
      for(int ibinx=1;ibinx<=pullstatEventsXsystEventsXIterationsFinalIteration[0][0]->GetNbinsX();ibinx++){
	Chi2NominalTotalFinalIteration[isys-1][ibinx-1] /= NDFTotalFinalIteration[isys-1][ibinx-1];
	Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->SetBinContent(ibinx,Chi2NominalTotalFinalIteration[isys-1][ibinx-1]);
      CorrelatedChi2NominalTotalFinalIteration[isys-1][ibinx-1] /= NDFTotal[isys-1][ibinx-1];
      CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->SetBinContent(ibinx,CorrelatedChi2NominalTotalFinalIteration[isys-1][ibinx-1]);
      //cout<<"Iteration #"<<ibinx<<", Chi2="<<Chi2NominalTotal[isys-1][ibinx-1]<<endl;
    }



      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->GetXaxis()->SetTitle("Number of iterations");
      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->GetYaxis()->SetTitle("#chi^{2} / NDF");
      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->Write();
      
      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->Draw();
      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->SetLineWidth(3);
      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->SetLineStyle(2);
      Chi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->GetYaxis()->SetRangeUser(1e-5,100.);

      CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->GetXaxis()->SetTitle("Number of iterations");
    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->GetYaxis()->SetTitle("#chi^{2} / NDF");
    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[isys-1]->Write();


      for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
	//cChi2DataStatErrorFinalIteration[isys-1][c0] = new TCanvas(Form("cChi2DataStatErrorFinalIteration for toy XP %d and mom %d",isys-1,c0),"");
	
	for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	  if(!UsedBinPlots[c0][c1]) continue;
	  Chi2DataStatErrorXIterationsFinalIteration[isys-1][c0][c1]->SetLineWidth(2);
	  Chi2DataStatErrorXIterationsFinalIteration[isys-1][c0][c1]->SetLineColor(color[(c0-FirstBinMomPlots)*NBinsTrueAnglePlots+(c1-FirstBinAnglePlots)]);
	  Chi2DataStatErrorXIterationsFinalIteration[isys-1][c0][c1]->Draw("same");
	  //int BinY=c0*(NBinsTrueAnglePlots)+c1+1;
	  //Chi2DataStatErrorXIterationsFinalIteration[c0][c1]->SetLineColor(BinY);
	  //Chi2DataStatErrorXIterationsFinalIteration[c0][c1]->Draw("same");
	  //l->AddEntry( Chi2DataStatErrorXIterationsFinalIteration[c0][c1] , Form("Mom %d, Angle %d",c0,c1));
	}
      }
    }
    lModified->Draw("same");
    line1->Draw("same");
    cChi2DataStatErrorFinalIteration->SetGridx();
    cChi2DataStatErrorFinalIteration->SetGridy();
    cChi2DataStatErrorFinalIteration->SetLogy();
    cChi2DataStatErrorFinalIteration->Write();

    TCanvas * cChi2DataStatErrorComparison = new TCanvas("cChi2DataStatErrorComparison");
    Chi2DataStatErrorTotalXIterations[0]->SetLineWidth(2);    
    Chi2DataStatErrorTotalXIterations[0]->SetLineStyle(1);    
    Chi2DataStatErrorTotalXIterations[0]->SetLineColor(kBlue+2);    
    Chi2DataStatErrorTotalXIterations[0]->Draw();
    Chi2DataStatErrorTotalXIterations[0]->GetXaxis()->SetRangeUser(0,50);
    TLine * line2 = new TLine(0.,1.,50,1.);
    line2->SetLineColor(1);
    line2->SetLineWidth(2);
    line2->SetLineStyle(10);
    line2->Draw("same");

    CorrelatedChi2DataStatErrorTotalXIterations[0]->SetLineWidth(2);    
    CorrelatedChi2DataStatErrorTotalXIterations[0]->SetLineStyle(1);    
    CorrelatedChi2DataStatErrorTotalXIterations[0]->SetLineColor(kGreen+2);    
    CorrelatedChi2DataStatErrorTotalXIterations[0]->Draw("same");
    
    Chi2DataStatErrorTotalXIterationsFinalIteration[0]->SetLineWidth(2);    
    Chi2DataStatErrorTotalXIterationsFinalIteration[0]->SetLineStyle(1);    
    Chi2DataStatErrorTotalXIterationsFinalIteration[0]->SetLineColor(kRed);    
    Chi2DataStatErrorTotalXIterationsFinalIteration[0]->Draw("same");    

    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[0]->SetLineWidth(2);    
    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[0]->SetLineStyle(1);    
    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[0]->SetLineColor(kMagenta);    
    CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[0]->Draw("same");

    TLegend * lComparison = new TLegend(0.75,0.2,0.89,0.89);
    lComparison->SetFillColor(0.);
    lComparison->SetLineColor(0.);    
    lComparison->AddEntry(Chi2DataStatErrorTotalXIterations[0],"#chi^{2}_{NC} / NDF");
    lComparison->AddEntry(CorrelatedChi2DataStatErrorTotalXIterations[0],"#chi^{2} / NDF");
    lComparison->AddEntry(Chi2DataStatErrorTotalXIterationsFinalIteration[0],"#chi^{2}_{NC,DD} / NDF");
    lComparison->AddEntry(CorrelatedChi2DataStatErrorTotalXIterationsFinalIteration[0],"#chi^{2}_{DD} / NDF");
    lComparison->Draw("same");
    
    cChi2DataStatErrorComparison->SetGridx();
    cChi2DataStatErrorComparison->SetGridy();
    cChi2DataStatErrorComparison->SetLogy();
    cChi2DataStatErrorComparison->Write();
  ////////////////////STAT VS SYST///////////////////////////////
  
    //Draw the pull distribution with iterations & behaviour of the mean/error of the pull, BUT ONLY FOR STATISTICAL VARIATIONS in the input file	
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	if(!UsedBinPlots[c0][c1]) continue;
	canstatEventsXIterations[c0][c1] = new TCanvas(Form("cstatPullXIterations%d_%d",c0,c1),Form("Pull maximum and 1#sigma error of events after unfolding in the bin %d_%d",c0,c1));
	pullstatEventsXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding iterations");
	pullstatEventsXIterations[c0][c1]->GetYaxis()->SetTitle("Pullstat");
	for(int ibinx=1;ibinx<=pullstatEventsXIterations[c0][c1]->GetNbinsX();ibinx++){
	  double Mean;
	  double MeanError;
	  double RMS;
	  double RMSError;
	  TH1D * Slice;
	  
	  if(Fit){
	    gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	    pullstatEventsXIterations[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	    Mean = gausTemp->GetParameter(1);
	    RMS = gausTemp->GetParameter(2);
	    MeanError = gausTemp->GetParError(1);
	    RMSError = gausTemp->GetParError(2);
	  }
	  else{
	    Slice = (TH1D*) pullstatEventsXIterations[c0][c1]->ProjectionY("Slice",ibinx,ibinx);
	    Mean = Slice->GetMean();
	    RMS = Slice->GetRMS();
	    MeanError = Slice->GetMeanError();
	    RMSError = Slice->GetRMSError();
	  }
	  pullstatEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,Mean);
	  pullstatEventsmeanXIterations[c0][c1]->SetBinError(ibinx,MeanError);
	  pullstatEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,RMS);
	  pullstatEventssigmaXIterations[c0][c1]->SetBinError(ibinx,RMSError);
	  //cout<<"center="<<Mean<<", error"<<RMS<<endl;	  
	}
	pullstatEventsmeanXIterations[c0][c1]->SetLineColor(kRed);
	pullstatEventssigmaXIterations[c0][c1]->SetLineColor(kOrange);
	pullstatEventsmeanXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding Iterations");
	pullstatEventsmeanXIterations[c0][c1]->Draw("E1");
	pullstatEventssigmaXIterations[c0][c1]->Draw("E1same");
	pullstatEventsmeanXIterations[c0][c1]->GetYaxis()->SetRangeUser(-0.2,0.2);
	canstatEventsXIterations[c0][c1]->Write();

	pullstatEventsXsystEvents[c0][c1]->Write();
	pullstatEventsmeanXsystEvents[c0][c1]->Write();	
	pullstatEventssigmaXsystEvents[c0][c1]->Write();	
      }
    }
    //StatisticalVariations->Close();

    
    //############### 3. STAT. VARIATIONS ONLY #################
    //4. For the total error, build the correlation matrix from the covariant one
    StatisticalVariations->cd();
    CovarianceStat->Divide(NVariationsStat);
    CovarianceStat_Iteration1->Divide(NVariationsStat);
    CovarianceStat_Iteration3->Divide(NVariationsStat);

    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	//if(!UsedBinPlots[c0][c1]) continue;
	for(int d0=0;d0<NBinsTrueMomPlots;d0++){//loop over cause 0
	  for(int d1=0;d1<NBinsTrueAnglePlots;d1++){//loop over cause 1
	    //if(!UsedBinPlots[d0][d1]) continue;
	    int BinX=c0*(NBinsTrueAnglePlots)+c1+1;
	    int BinY=d0*(NBinsTrueAnglePlots)+d1+1;

	    double Corr=CovarianceStat->GetBinContent(BinX,BinY);
	    double NormCorr=TMath::Sqrt(CovarianceStat->GetBinContent(BinX,BinX)*CovarianceStat->GetBinContent(BinY,BinY));
	    if(NormCorr) Corr/=NormCorr;
	    CorrelationStat->SetBinContent(BinX,BinY,Corr);

	    Corr=CovarianceStat_Iteration1->GetBinContent(BinX,BinY);
	    NormCorr=TMath::Sqrt(CovarianceStat_Iteration1->GetBinContent(BinX,BinX)*CovarianceStat_Iteration1->GetBinContent(BinY,BinY));
	    if(NormCorr) Corr/=NormCorr;
	    CorrelationStat_Iteration1->SetBinContent(BinX,BinY,Corr);

	    Corr=CovarianceStat_Iteration3->GetBinContent(BinX,BinY);
	    NormCorr=TMath::Sqrt(CovarianceStat_Iteration3->GetBinContent(BinX,BinX)*CovarianceStat_Iteration3->GetBinContent(BinY,BinY));
	    if(NormCorr) Corr/=NormCorr;
	    CorrelationStat_Iteration3->SetBinContent(BinX,BinY,Corr);	    
	  }
	}
      }
    }
    //
    TCanvas * canCorrelationStat = new TCanvas("cCorrelationStat");
    gStyle->SetPaintTextFormat("2.2f");
    gStyle->SetOptStat(kFALSE);
    CorrelationStat->GetXaxis()->SetLabelSize(0.);
    CorrelationStat->GetYaxis()->SetLabelSize(0.);
    CorrelationStat->Draw("colztext");
    for(int imom=0;imom<NBinsTrueMomPlots;imom++){
      tBins[imom][0]->Draw("same");
      tBins[imom][1]->Draw("same");
      lBins[imom][0]->Draw("same");
      lBins[imom][1]->Draw("same");
    }
    canCorrelationStat->Write();
   
    TCanvas * canCorrelationStat_Iteration1 = new TCanvas("cCorrelationStat_Iteration1");
    gStyle->SetPaintTextFormat("2.2f");
    gStyle->SetOptStat(kFALSE);
    CorrelationStat_Iteration1->GetXaxis()->SetLabelSize(0.);
    CorrelationStat_Iteration1->GetYaxis()->SetLabelSize(0.);
    CorrelationStat_Iteration1->Draw("colztext");
    for(int imom=0;imom<NBinsTrueMomPlots;imom++){
      tBins[imom][0]->Draw("same");
      tBins[imom][1]->Draw("same");
      lBins[imom][0]->Draw("same");
      lBins[imom][1]->Draw("same");
    }
    canCorrelationStat_Iteration1->Write();

    TCanvas * canCorrelationStat_Iteration3 = new TCanvas("cCorrelationStat_Iteration3");
    gStyle->SetPaintTextFormat("2.2f");
    gStyle->SetOptStat(kFALSE);
    CorrelationStat_Iteration3->GetXaxis()->SetLabelSize(0.);
    CorrelationStat_Iteration3->GetYaxis()->SetLabelSize(0.);
    CorrelationStat_Iteration3->Draw("colztext");
    for(int imom=0;imom<NBinsTrueMomPlots;imom++){
      tBins[imom][0]->Draw("same");
      tBins[imom][1]->Draw("same");
      lBins[imom][0]->Draw("same");
      lBins[imom][1]->Draw("same");
    }
    canCorrelationStat_Iteration3->Write();

    //////////////////////////////////////////////////////////////////////////////



    

    //############### 4. SYST. VARIATIONS ONLY #################
    SystematicVariations->cd();
    for(int c0=0;c0<NBinsTrueMomPlots;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAnglePlots;c1++){//loop over cause 1
	if(!UsedBinPlots[c0][c1]) continue;
	//Systematic variation
	//Draw the pull distribution with iterations & behaviour of the mean/error of the pull, BUT ONLY FOR SYSTEMATICS VARIATIONS in the input file	
	cansystEventsXIterations[c0][c1] = new TCanvas(Form("csystPullXIterations%d_%d",c0,c1),Form("Pull maximum and 1#sigma error of events after unfolding in the bin %d_%d",c0,c1));
	pullsystEventsXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding iterations");
	pullsystEventsXIterations[c0][c1]->GetYaxis()->SetTitle("Pullsyst");
	for(int ibinx=1;ibinx<=pullsystEventsXIterations[c0][c1]->GetNbinsX();ibinx++){

	  double Mean;
	  double MeanError;
	  double RMS;
	  double RMSError;
	  TH1D * Slice;
	 
 	  if(Fit){
	    gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	    pullsystEventsXIterations[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	    Mean = gausTemp->GetParameter(1);
	    RMS = gausTemp->GetParameter(2);
	    MeanError = gausTemp->GetParError(1);
	    RMSError = gausTemp->GetParError(2);
	  }
	  else{
	    Slice = (TH1D*) pullsystEventsXIterations[c0][c1]->ProjectionY("Slice",ibinx,ibinx);
	    Mean = Slice->GetMean();
	    RMS = Slice->GetRMS();
	    MeanError = Slice->GetMeanError();
	    RMSError = Slice->GetRMSError();
	  }

	  pullsystEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,Mean);
	  pullsystEventsmeanXIterations[c0][c1]->SetBinError(ibinx,MeanError);
	  pullsystEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,RMS);
	  pullsystEventssigmaXIterations[c0][c1]->SetBinError(ibinx,RMSError);
	  
	}
	pullsystEventsmeanXIterations[c0][c1]->SetLineColor(kRed);
	pullsystEventssigmaXIterations[c0][c1]->SetLineColor(kOrange);
	pullsystEventsmeanXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding Iterations");
	pullsystEventsmeanXIterations[c0][c1]->Write();
	pullsystEventsmeanXIterations[c0][c1]->Draw("E1same");
	pullsystEventssigmaXIterations[c0][c1]->Draw("E1same");
	pullsystEventsmeanXIterations[c0][c1]->GetYaxis()->SetRangeUser(-0.05,0.08);
	cansystEventsXIterations[c0][c1]->Write();
	
	
	//Draw the pull distribution with systematic bias & behaviour of the mean/error of the pull. PUll is only done FOR STATISTICAL VARIATIONS in the input file, but x axis is the systematic variation!!	
	pullstatEventsXsystEvents[c0][c1]->GetXaxis()->SetTitle("Systematically varied exp.");
	pullstatEventsXsystEvents[c0][c1]->GetYaxis()->SetTitle("Pullstat");
	for(int ibinx=1;ibinx<=pullstatEventsXsystEvents[c0][c1]->GetNbinsX();ibinx++){
	  double Mean;
	  double MeanError;
	  double RMS;
	  double RMSError;
	  TH1D * Slice;
	  if(Fit){
	    gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	    pullstatEventsXsystEvents[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	    Mean = gausTemp->GetParameter(1);
	    RMS = gausTemp->GetParameter(2);
	    MeanError = gausTemp->GetParError(1);
	    RMSError = gausTemp->GetParError(2);
	  }
	  else{
	    Slice = (TH1D*) pullstatEventsXsystEvents[c0][c1]->ProjectionY("Slice",ibinx,ibinx);
	    Mean = Slice->GetMean();
	    RMS = Slice->GetRMS();
	    MeanError = Slice->GetMeanError();
	    RMSError = Slice->GetRMSError();
	  }
	  pullstatEventsmeanXsystEvents[c0][c1]->SetBinContent(ibinx,Mean);
	  pullstatEventsmeanXsystEvents[c0][c1]->SetBinError(ibinx,MeanError);
	  pullstatEventssigmaXsystEvents[c0][c1]->SetBinContent(ibinx,RMS);
	  pullstatEventssigmaXsystEvents[c0][c1]->SetBinError(ibinx,RMSError);
	  //cout<<"center="<<Mean<<", error"<<RMS<<endl;
	}
	//leg->Draw("lsame");
      }
    }
    //SystematicVariations->Close();

    General->cd();
    TCanvas * canNumberOfEvents = new TCanvas("cNumberOfEvents_True");
    gStyle->SetPaintTextFormat("2.2f");
    gStyle->SetOptStat(kFALSE);
    NumberOfEvents->Draw("colztext");
    canNumberOfEvents->Write();

    TVectorD * VReconstructedEvents = (TVectorD*) file->Get("VReconstructedEvents");
    cout<<"#########################"<<endl<<endl;
    cout<<NBinsTrueMomSignal<<endl;


    TH2D * Prior = (TH2D*) file->Get(Form("Prior%d",0));
    TH1D * PriorMom = (TH1D*) Prior->ProjectionX("PriorMom",0,Prior->GetNbinsY());
    TH1D * PriorAngle = (TH1D*) Prior->ProjectionY("PriorAngle",0,Prior->GetNbinsX());
    PriorMom->Scale(1./PriorMom->Integral());
    PriorAngle->Scale(1./PriorAngle->Integral());
    Prior->Write(); PriorMom->Write(); PriorAngle->Write();

 

    TH2D * NEvents[NMaxIterations+1];
    TH1D * NEventsMom[NMaxIterations+1];
    TH1D * NEventsAngle[NMaxIterations+1];
    for(int it=0;it<=NMaxIterations;it++){ 
      NEvents[it] = new TH2D(Form("NEvents_Iteration%d",it),"",NBinsTrueMomPlots,0,NBinsTrueMomPlots,NBinsTrueAnglePlots,0,NBinsTrueAnglePlots);
      for(int imom=0;imom<NBinsTrueMomPlots;imom++){
	for(int iang=0;iang<NBinsTrueAnglePlots;iang++){
	  if(!UsedBinPlots[imom][iang]) continue;
  	  NEvents[it]->SetBinContent(imom+1,iang+1,EventXIterations[imom][iang]->GetBinContent(it+1));
	}
      }
      
      NEventsMom[it] = (TH1D*) NEvents[it]->ProjectionX(Form("NEventsMom_Iteration%d",it),0,NEvents[it]->GetNbinsY());
      NEventsAngle[it] = (TH1D*) NEvents[it]->ProjectionY(Form("NEventsAngle_Iteration%d",it),0,NEvents[it]->GetNbinsX());
    }
    
    TLegend * legNEvents = new TLegend(0.7,0.7,1.,1.);
    legNEvents->SetFillColor(0.);
    legNEvents->SetLineColor(0.);    

    /*    
    TH2D * Prior[NMaxIterations+1];
    TH1D * PriorMom[NMaxIterations+1];
    TH1D * PriorAngle[NMaxIterations+1];

    for(int it=0;it<=NMaxIterations;it++){ 
      Prior[it] = (TH2D*) file->Get(Form("Prior%d",it));
      NEvents[it] = new TH2D(Form("NEvents_Iteration%d",it),"",NBinsTrueMomPlots,0,NBinsTrueMomPlots,NBinsTrueAnglePlots,0,NBinsTrueAnglePlots);
      for(int imom=0;imom<NBinsTrueMomPlots;imom++){
	for(int iang=0;iang<NBinsTrueAnglePlots;iang++){
	  NEvents[it]->SetBinContent(imom+1,iang+1,EventXIterations[imom][iang]->GetBinContent(it+1));
	}
      }
      //Prior[it]->Write();
      //
      PriorMom[it] = new TH1D(Form("PriorMom_Iteration%d",it),"",NBinsTrueMomPlots,0,NBinsTrueMomPlots);
      NEventsMom[it] = new TH1D(Form("NEventsMom_Iteration%d",it),"",NBinsTrueMomPlots,0,NBinsTrueMomPlots);
      for(int imom=0;imom<NBinsTrueMomPlots;imom++){
	double prior1D=0;
	double nevents1D=0;
	for(int iang=0;iang<NBinsTrueAnglePlots;iang++){
	  prior1D += Prior[it]->GetBinContent(imom+1,iang+1);
	  nevents1D += NEvents[it]->GetBinContent(imom+1,iang+1);
	}
	PriorMom[it]->SetBinContent(imom+1,prior1D);
	NEventsMom[it]->SetBinContent(imom+1,nevents1D);
      }
      //PriorMom[it]->Write();
      
      PriorAngle[it] = new TH1D(Form("PriorAngle_Iteration%d",it),"",NBinsTrueAnglePlots,0,NBinsTrueAnglePlots);
      NEventsAngle[it] = new TH1D(Form("NEventsAngle_Iteration%d",it),"",NBinsTrueAnglePlots,0,NBinsTrueAnglePlots);
      for(int iang=0;iang<NBinsTrueAnglePlots;iang++){
	double prior1D=0; 
	double nevents1D=0; 
	for(int imom=0;imom<NBinsTrueMomPlots;imom++){
	  prior1D += Prior[it]->GetBinContent(imom+1,iang+1);
	  nevents1D += NEvents[it]->GetBinContent(imom+1,iang+1);
	}
	PriorAngle[it]->SetBinContent(iang+1,prior1D);
	NEventsAngle[it]->SetBinContent(iang+1,nevents1D);
      }
      //PriorAngle[it]->Write();
    }
    
    Prior[0]->Write();
    NEvents[0]->Write();
    TLegend * legPrior = new TLegend(0.7,0.7,1.,1.);
    legPrior->SetFillColor(0.);
    legPrior->SetLineColor(0.);    

    TCanvas * cPriorMom = new TCanvas("cPriorMom");
    TH1D * TruePriorMomFull = (TH1D*) TrueDistribution->ProjectionX("TruePriorMomFull",0,TrueDistribution->GetNbinsY());
    TH1D * TruePriorMom = new TH1D("TruePriorMom","",NBinsTrueMomPlots,0,NBinsTrueMomPlots);
    for(int imom=0;imom<NBinsTrueMomPlots;imom++) TruePriorMom->SetBinContent(imom+1,TruePriorMomFull->GetBinContent(imom+1));

    TruePriorMom->Scale(1./TruePriorMom->Integral());
    TruePriorMom->SetLineColor(kBlue);
    TruePriorMom->SetLineStyle(2);
    TruePriorMom->SetLineWidth(3);
    TruePriorMom->Draw();
    legPrior->AddEntry(TruePriorMom,"MC true");
    //TruePriorMom->SetYaxis()->SetRangeUser(0,;
    for(int it=0;it<=NMaxIterations;it++){
      
      PriorMom[it]->Scale(1./PriorMom[it]->Integral());
      
      if(it == 0 || it == 1 || it == 2 || it == 5 || it == 10 || it == 100 || it == 1000){
	if(it == 0) PriorMom[it]->SetLineColor(1);
	else if(it == 1) PriorMom[it]->SetLineColor(kRed);
	else if(it == 2) PriorMom[it]->SetLineColor(kGreen+2);
	else if(it == 5) PriorMom[it]->SetLineColor(kBlue);
	else if(it == 10) PriorMom[it]->SetLineColor(kMagenta);
	else if(it == 100) PriorMom[it]->SetLineColor(kCyan);
	else if(it == 1000) PriorMom[it]->SetLineColor(kGray);
	PriorMom[it]->SetLineWidth(2);
	PriorMom[it]->Draw("same");
	legPrior->AddEntry(PriorMom[it],Form("%d iterations",it));
      }
    }
    legPrior->Draw("same");
    TAxis * aX_PriorMom = (TAxis*) TruePriorMom->GetXaxis();
    for(int ibin=1;ibin<=NBinsTrueMomPlots;ibin++) aX_PriorMom->SetBinLabel(ibin,Form(Form("%4.1fGeV-%4.1fGeV",BinningTrueMom[ibin-1],BinningTrueMom[ibin])));
    cPriorMom->Write();
    
    TCanvas * cPriorAngle = new TCanvas("cPriorAngle");
    TH1D * TruePriorAngleFull = (TH1D*) TrueDistribution->ProjectionY("TruePriorAngleFull",0,TrueDistribution->GetNbinsX());
    TH1D * TruePriorAngle = new TH1D("TruePriorAngle","",NBinsTrueAnglePlots,0,NBinsTrueAnglePlots);
    for(int imom=0;imom<NBinsTrueAnglePlots;imom++) TruePriorAngle->SetBinContent(imom+1,TruePriorAngleFull->GetBinContent(imom+1));

    TruePriorAngle->Scale(1./TruePriorAngle->Integral());
    TruePriorAngle->SetLineColor(kBlue);
    TruePriorAngle->SetLineWidth(3);
    TruePriorAngle->SetLineStyle(2);
    TruePriorAngle->Draw();
    
    for(int it=0;it<=NMaxIterations;it++){
      
      PriorAngle[it]->Scale(1./PriorAngle[it]->Integral());
      
      if(it == 0 || it == 1 || it == 2 || it == 5 || it == 10 || it == 100 || it == 1000){
	if(it == 0) PriorAngle[it]->SetLineColor(1);
	else if(it == 1) PriorAngle[it]->SetLineColor(kRed);
	else if(it == 2) PriorAngle[it]->SetLineColor(kGreen+2);
	else if(it == 5) PriorAngle[it]->SetLineColor(kBlue);
	else if(it == 10) PriorAngle[it]->SetLineColor(kMagenta);
	else if(it == 100) PriorAngle[it]->SetLineColor(kCyan);
	else if(it == 1000) PriorAngle[it]->SetLineColor(kGray);
	PriorAngle[it]->SetLineWidth(2);
	PriorAngle[it]->Draw("same");
      }
    }
    TAxis * aX_PriorAngle = (TAxis*) TruePriorAngle->GetXaxis();
    for(int ibin=1;ibin<=NBinsTrueAnglePlots;ibin++) aX_PriorAngle->SetBinLabel(ibin,Form(Form("%2.0f#circ-%2.0f#circ",BinningTrueAngle[ibin-1],BinningTrueAngle[ibin])));
    legPrior->Draw("same");
    cPriorAngle->Write();
    */
    
    //TruePriorMom->SetYaxis()->SetRangeUser(0,;

    ////////////////////////////////////////////////////
    TCanvas * cNEventsMom = new TCanvas("cNEventsMom");
    cNEventsMom->Divide(1,2);
    cNEventsMom->cd(1);
    TH1D * TrueNEventsMomFull = (TH1D*) TrueDistribution->ProjectionX("TrueNEventsMomFull",0,TrueDistribution->GetNbinsY());
    TH1D * TrueNEventsMom = new TH1D("TrueNEventsMom","",NBinsTrueMomPlots,0,NBinsTrueMomPlots);
    for(int imom=0;imom<NBinsTrueMomPlots;imom++){
      if(!UsedBinPlots[imom][0]) continue;
      TrueNEventsMom->SetBinContent(imom+1,TrueNEventsMomFull->GetBinContent(imom+1));
    }
    TrueNEventsMom->SetLineColor(kBlue);
    TrueNEventsMom->SetLineStyle(2);
    TrueNEventsMom->SetLineWidth(3);
    TrueNEventsMom->Draw();
    legNEvents->AddEntry(TrueNEventsMom,"MC true");
    PriorMom->Scale(TrueNEventsMom->Integral());
    PriorMom->SetLineColor(kBlack);
    PriorMom->SetLineWidth(2);
    PriorMom->Draw("same");
    legNEvents->AddEntry(PriorMom,"0 iterations");
    //TrueNEventsMom->SetYaxis()->SetRangeUser(0,;
    for(int it=0;it<=NMaxIterations;it++){  
      if(it == 0 || it == 1 || it == 4 || it == 9 || it == 99 || it == 999){
	if(it == 0) NEventsMom[it]->SetLineColor(kRed);
	else if(it == 1) NEventsMom[it]->SetLineColor(kBlue);
	else if(it == 4) NEventsMom[it]->SetLineColor(kGreen+2);
	else if(it == 9) NEventsMom[it]->SetLineColor(kMagenta);
	else if(it == 99) NEventsMom[it]->SetLineColor(kCyan);
	else if(it == 999) NEventsMom[it]->SetLineColor(kGray);
	NEventsMom[it]->SetLineWidth(2);
	NEventsMom[it]->Draw("same");
	legNEvents->AddEntry(NEventsMom[it],Form("%d iterations",it+1));
      }
    }
    legNEvents->Draw("same");
    TAxis * aX_NEventsMom = (TAxis*) TrueNEventsMom->GetXaxis();
    for(int ibin=1;ibin<=NBinsTrueMomPlots;ibin++){
      if(!UsedBinPlots[ibin][0]) continue;
      aX_NEventsMom->SetBinLabel(ibin,Form(Form("%4.1fGeV-%4.1fGeV",BinningTrueMom[ibin-1],BinningTrueMom[ibin])));
    }
    aX_NEventsMom->SetLabelSize(0.05);
    aX_NEventsMom->SetTitleSize(0.05);
    TrueNEventsMom->GetYaxis()->SetLabelSize(0.05);
    TrueNEventsMom->GetYaxis()->SetTitleSize(0.05);
      
    cNEventsMom->cd(2);
    TH1D * TrueNEventsMomRelative = (TH1D*) TrueNEventsMom->Clone("TrueNEventsMomRelative");
    TrueNEventsMomRelative->Divide(TrueNEventsMom);
    TrueNEventsMomRelative->GetYaxis()->SetTitle("#frac{Unfolded}{True}");
    TrueNEventsMomRelative->GetYaxis()->CenterTitle();
    TrueNEventsMomRelative->Draw();
    TH1D * NEventsMomRelative[NMaxIterations+1];
    
    for(int it=0;it<=NMaxIterations;it++){
      NEventsMomRelative[it] = (TH1D*) NEventsMom[it]->Clone(Form("NEventsMomRelative[%d]",it));
      NEventsMomRelative[it]->Divide(TrueNEventsMom);
      if(it == 0 || it == 1 || it == 4 || it == 9 || it == 99 || it == 999) NEventsMomRelative[it]->Draw("same");
	//legNEvents->AddEntry(NEventsMom[it],Form("%d iterations",it));
    }
    cNEventsMom->Write();


    TCanvas * cNEventsAngle = new TCanvas("cNEventsAngle");
    cNEventsAngle->Divide(1,2);
    cNEventsAngle->cd(1);
    TH1D * TrueNEventsAngleFull = (TH1D*) TrueDistribution->ProjectionY("TrueNEventsAngleFull",0,TrueDistribution->GetNbinsX());
    TH1D * TrueNEventsAngle = new TH1D("TrueNEventsAngle","",NBinsTrueAnglePlots,0,NBinsTrueAnglePlots);
    for(int iang=0;iang<NBinsTrueAnglePlots;iang++){
      if(!UsedBinPlotsAngle[iang]) continue;
      TrueNEventsAngle->SetBinContent(iang+1,TrueNEventsAngleFull->GetBinContent(iang+1));
    }
    TrueNEventsAngle->SetLineColor(kBlue);
    TrueNEventsAngle->SetLineStyle(2);
    TrueNEventsAngle->SetLineWidth(3);
    TrueNEventsAngle->Draw();
    PriorAngle->Scale(TrueNEventsAngle->Integral());
    PriorAngle->SetLineColor(kBlack);
    PriorAngle->SetLineWidth(2);
    PriorAngle->Draw("same");
    //TrueNEventsAngle->SetYaxis()->SetRangeUser(0,;
    for(int it=0;it<=NMaxIterations;it++){  
      if(it == 0 || it == 1 || it == 4 || it == 9 || it == 99 || it == 999){
	if(it == 0) NEventsAngle[it]->SetLineColor(kRed);
	else if(it == 1) NEventsAngle[it]->SetLineColor(kBlue);
	else if(it == 4) NEventsAngle[it]->SetLineColor(kGreen+2);
	else if(it == 9) NEventsAngle[it]->SetLineColor(kMagenta);
	else if(it == 99) NEventsAngle[it]->SetLineColor(kCyan);
	else if(it == 999) NEventsAngle[it]->SetLineColor(kGray);
	NEventsAngle[it]->SetLineWidth(2);
	NEventsAngle[it]->Draw("same");
      }
    }
    legNEvents->Draw("same");
    TAxis * aX_NEventsAngle = (TAxis*) TrueNEventsAngle->GetXaxis();
    for(int ibin=1;ibin<=NBinsTrueAnglePlots;ibin++){
      if(!UsedBinPlotsAngle[ibin]) continue;
      aX_NEventsAngle->SetBinLabel(ibin,Form(Form("%4.1fGeV-%4.1fGeV",BinningTrueAngle[ibin-1],BinningTrueAngle[ibin])));
    }
    aX_NEventsAngle->SetLabelSize(0.05);
    aX_NEventsAngle->SetTitleSize(0.05);
    TrueNEventsAngle->GetYaxis()->SetLabelSize(0.05);
    TrueNEventsAngle->GetYaxis()->SetTitleSize(0.05);
      
    cNEventsAngle->cd(2);
    TH1D * TrueNEventsAngleRelative = (TH1D*) TrueNEventsAngle->Clone("TrueNEventsAngleRelative");
    TrueNEventsAngleRelative->Divide(TrueNEventsAngle);
    TrueNEventsAngleRelative->GetYaxis()->SetTitle("#frac{Unfolded}{True}");
    TrueNEventsAngleRelative->GetYaxis()->CenterTitle();
    TrueNEventsAngleRelative->Draw();
    TH1D * NEventsAngleRelative[NMaxIterations+1];
    
    for(int it=0;it<=NMaxIterations;it++){
      NEventsAngleRelative[it] = (TH1D*) NEventsAngle[it]->Clone(Form("NEventsAngleRelative[%d]",it));
      NEventsAngleRelative[it]->Divide(TrueNEventsAngle);
      if(it == 0 || it == 1 || it == 4 || it == 9 || it == 99 || it == 999) NEventsAngleRelative[it]->Draw("same");
	//legNEvents->AddEntry(NEventsAngle[it],Form("%d iterations",it));
    }
    cNEventsAngle->Write();




    
    TMatrixD * MLikelihood = (TMatrixD*) file->Get("MLikelihood");

    TCanvas * cLikelihood = new TCanvas("cLikelihood");
    MLikelihood->Draw("colztext");
    cLikelihood->Write();
    TH2D * LikelihoodMom = new TH2D("LikelihoodMom","",NBinsRecMom,BinningRecMom,NBinsTrueMom,BinningTrueMom);//It is equal to sum(L2D)overRecAndTrueangles / Sum(prior)_overTrueAngle
    LikelihoodMom->GetXaxis()->SetTitle("d_{#mu} (cm)");
    LikelihoodMom->GetYaxis()->SetTitle("p_{#mu} (GeV/c)");
    for(int rmom=0;rmom<NBinsRecMom;rmom++){
      for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	//First, estimate the numerator of LikelihoodMom
	double L1D=0;
	for(int rang=0;rang<NBinsRecAngle;rang++){
	  for(int tang=0;tang<NBinsTrueAngle;tang++){
	    int bintrue=tmom*NBinsTrueAngle+tang;
	    int binrec=rmom*NBinsRecAngle+rang;
	    double L2D=(*MLikelihood)(bintrue,binrec);
	    L1D += L2D*Prior->GetBinContent(tmom+1,tang+1);
	  }
	}
	//Second, the denominator
	double prior1D=0;
	for(int tang=0;tang<NBinsTrueAngle;tang++){
	  double prior2D = Prior->GetBinContent(tmom+1,tang+1);
	  prior1D += prior2D;
	}
	//Divide
	if(prior1D !=0 ) L1D /= prior1D;
	LikelihoodMom->SetBinContent(rmom+1,tmom+1,L1D);
      }
    }

    for(int rmom=0;rmom<NBinsRecMom;rmom++){
      double LikelihoodNorm = 0;
      for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	LikelihoodNorm += LikelihoodMom->GetBinContent(rmom+1,tmom+1);
      }
   
      for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	if(LikelihoodNorm != 0) LikelihoodMom->SetBinContent(rmom+1,tmom+1,LikelihoodMom->GetBinContent(rmom+1,tmom+1)/LikelihoodNorm);
      }
    }
    LikelihoodMom->Write();

    TH2D * LikelihoodAngle = new TH2D("LikelihoodAngle","",NBinsRecAngle,BinningRecAngle,NBinsTrueAngle,BinningTrueAngle);//It is equal to sum(L2D)overRecAndTrueangles / Sum(prior)_overTrueMom
    LikelihoodAngle->GetXaxis()->SetTitle("rec #theta_{#mu} (#circ)");
    LikelihoodAngle->GetYaxis()->SetTitle("true #theta_{#mu} (#circ)");
    for(int rang=0;rang<NBinsRecAngle;rang++){
      for(int tang=0;tang<NBinsTrueAngle;tang++){
	//First, estimate the numerator of LikelihoodAngle
	double L1D=0;
	for(int rmom=0;rmom<NBinsRecMom;rmom++){
	  for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	    int bintrue=tmom*NBinsTrueAngle+tang;
	    int binrec=rmom*NBinsRecAngle+rang;
	    double L2D=(*MLikelihood)(bintrue,binrec);
	    L1D += L2D*Prior->GetBinContent(tmom+1,tang+1);
	  }
	}
	//Second, the denominator
	double prior1D=0;
	for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	  double prior2D = Prior->GetBinContent(tmom+1,tang+1);
	  prior1D += prior2D;
	}
	//Divide
	if(prior1D !=0 ) L1D /= prior1D;
	LikelihoodAngle->SetBinContent(rang+1,tang+1,L1D);
      }
    }
    LikelihoodAngle->Write();
      
    TMatrixD * MUnfolding = (TMatrixD*) file->Get("MUnfolding");
    TCanvas * cUnfolding = new TCanvas("cUnfolding");
    MUnfolding->Draw("colztext");
    cUnfolding->Write();

    TH2D * UnfoldingMom = new TH2D("UnfoldingMom","",NBinsRecMom,BinningRecMom,NBinsTrueMom,BinningTrueMom);//It is equal to sum(L2D)overRecAndTrueangles / Sum(prior)_overTrueAngle
    UnfoldingMom->GetXaxis()->SetTitle("d_{#mu} (cm)");
    UnfoldingMom->GetYaxis()->SetTitle("p_{#mu} (GeV/c)");
    for(int rmom=0;rmom<NBinsRecMom;rmom++){
      for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	//First, estimate the numerator of UnfoldingMom
	double L1D=0;
	for(int rang=0;rang<NBinsRecAngle;rang++){
	  for(int tang=0;tang<NBinsTrueAngle;tang++){
	    int bintrue=tmom*NBinsTrueAngle+tang;
	    int binrec=rmom*NBinsRecAngle+rang;
	    double L2D=(*MUnfolding)(bintrue,binrec);
	    L1D += L2D*Prior->GetBinContent(tmom+1,tang+1);
	  }
	}
	//Second, the denominator
	double prior1D=0;
	for(int tang=0;tang<NBinsTrueAngle;tang++){
	  double prior2D = Prior->GetBinContent(tmom+1,tang+1);
	  prior1D += prior2D;
	}
	//Divide
	if(prior1D !=0 ) L1D /= prior1D;
	UnfoldingMom->SetBinContent(rmom+1,tmom+1,L1D);
      }
    }

    for(int rmom=0;rmom<NBinsRecMom;rmom++){
      double UnfoldingNorm = 0;
      for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	UnfoldingNorm += UnfoldingMom->GetBinContent(rmom+1,tmom+1);
      }
   
      for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	if(UnfoldingNorm != 0) UnfoldingMom->SetBinContent(rmom+1,tmom+1,UnfoldingMom->GetBinContent(rmom+1,tmom+1)/UnfoldingNorm);
      }
    }
    UnfoldingMom->Write();

    TH2D * UnfoldingAngle = new TH2D("UnfoldingAngle","",NBinsRecAngle,BinningRecAngle,NBinsTrueAngle,BinningTrueAngle);//It is equal to sum(L2D)overRecAndTrueangles / Sum(prior)_overTrueMom
    UnfoldingAngle->GetXaxis()->SetTitle("rec #theta_{#mu} (#circ)");
    UnfoldingAngle->GetYaxis()->SetTitle("true #theta_{#mu} (#circ)");
    for(int rang=0;rang<NBinsRecAngle;rang++){
      for(int tang=0;tang<NBinsTrueAngle;tang++){
	//First, estimate the numerator of UnfoldingAngle
	double L1D=0;
	for(int rmom=0;rmom<NBinsRecMom;rmom++){
	  for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	    int bintrue=tmom*NBinsTrueAngle+tang;
	    int binrec=rmom*NBinsRecAngle+rang;
	    double L2D=(*MUnfolding)(bintrue,binrec);
	    L1D += L2D*Prior->GetBinContent(tmom+1,tang+1);
	  }
	}
	//Second, the denominator
	double prior1D=0;
	for(int tmom=0;tmom<NBinsTrueMom;tmom++){
	  double prior2D = Prior->GetBinContent(tmom+1,tang+1);
	  prior1D += prior2D;
	}
	//Divide
	if(prior1D !=0 ) L1D /= prior1D;
	UnfoldingAngle->SetBinContent(rang+1,tang+1,L1D);
      }
    }
    UnfoldingAngle->Write();

    
    

    //////////////////2D case///////////////////////////////////-> Not used (cf Memo from the 2015/12/09 
    
    TH2D * ReconstructedLikelihood[NBinsTrueMom][NBinsTrueAngle];
    TF2 * Gaussian2D[NBinsTrueMom][NBinsTrueAngle];
    TH1D * ReconstructedLikelihood_RecMom[NBinsTrueMom][NBinsTrueAngle];
    TH1D * ReconstructedLikelihood_RecAngle[NBinsTrueMom][NBinsTrueAngle];
    TF1 * Gaussian1D_RecMom[NBinsTrueMom][NBinsTrueAngle];
    TF1 * Gaussian1D_RecAngle[NBinsTrueMom][NBinsTrueAngle];
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	ReconstructedLikelihood[c0][c1]= new TH2D(Form("ReconstructedLikelihood[%d][%d]",c0,c1),"",NBinsRecMom,BinningRecMom,NBinsRecAngle,BinningRecAngle);

	
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    int bintrue=c0*NBinsTrueAngle+c1;
	    int binrec=e0*NBinsRecAngle+e1;
	    //cout<<bintrue<<", "<<binrec<<endl;
	    //cout<<(*MLikelihood)(bintrue,binrec)<<endl;
	    //cout<<(*MLikelihood)[bintrue][binrec]<<endl;
	    double correlation=(*MLikelihood)(bintrue,binrec);
	    ReconstructedLikelihood[c0][c1]->SetBinContent(e0+1,e1+1,correlation);
	  }
	}
	double minrecmom=BinningRecMom[0];
	double maxrecmom=BinningRecMom[NBinsRecMom];
	double minrecangle=BinningRecAngle[0];
	double maxrecangle=BinningRecAngle[NBinsRecAngle];
	
	Gaussian2D[c0][c1] = new TF2(Form("Gaussian2D[%d][%d]",c0,c1),"[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",minrecmom,maxrecmom,minrecangle,maxrecangle);
	Gaussian2D[c0][c1]->SetParameters(1,10,10,10,10);
	ReconstructedLikelihood[c0][c1]->Fit(Gaussian2D[c0][c1],"RQ0");
	ReconstructedLikelihood[c0][c1]->Write(Form("ReconstructedLikelihood[%d][%d]",c0,c1));


	
	ReconstructedLikelihood_RecMom[c0][c1] = (TH1D*) ReconstructedLikelihood[c0][c1]->ProjectionX(Form("ReconstructedLikelihood_RecMom[%d][%d]",c0,c1),1,NBinsRecAngle);
	ReconstructedLikelihood_RecAngle[c0][c1] = (TH1D*) ReconstructedLikelihood[c0][c1]->ProjectionY(Form("ReconstructedLikelihood_RecAngle[%d][%d]",c0,c1),1,NBinsRecMom);
	/*	Gaussian1D_RecMom[c0][c1] = new TF1(Form("Gaussian1D_RecMom[%d][%d]",c0,c1),"[0]*TMath::Gaus(x,[1],[2])",minrecmom,maxrecmom);
	Gaussian1D_RecMom[c0][c1]->SetParameters(ReconstructedLikelihood_RecMom[c0][c1]->Integral(),ReconstructedLikelihood_RecMom[c0][c1]->GetMaximum(),1);
	Gaussian1D_RecAngle[c0][c1] = new TF1(Form("Gaussian1D_RecAngle[%d][%d]",c0,c1),"[0]*TMath::Gaus(x,[1],[2])",minrecangle,maxrecangle);
	Gaussian1D_RecAngle[c0][c1]->SetParameters(ReconstructedLikelihood_RecAngle[c0][c1]->Integral(),ReconstructedLikelihood_RecAngle[c0][c1]->GetMaximum(),1);
	*/
	Gaussian1D_RecMom[c0][c1] = new TF1(Form("Gaussian1D_RecMom[%d][%d]",c0,c1),"gaus",minrecmom,maxrecmom);
	Gaussian1D_RecAngle[c0][c1] = new TF1(Form("Gaussian1D_RecAngle[%d][%d]",c0,c1),"gaus",minrecangle,maxrecangle);

	ReconstructedLikelihood_RecMom[c0][c1]->Fit(Gaussian1D_RecMom[c0][c1],"RQ0");
	ReconstructedLikelihood_RecMom[c0][c1]->Write(Form("ReconstructedLikelihood_RecMom[%d][%d]",c0,c1));
	ReconstructedLikelihood_RecAngle[c0][c1]->Fit(Gaussian1D_RecAngle[c0][c1],"RQ0");
	ReconstructedLikelihood_RecAngle[c0][c1]->Write(Form("ReconstructedLikelihood_RecAngle[%d][%d]",c0,c1));

	//cout<<setprecision(1)<<"True=("<<BinningTrueMom[c0]<<"-"<<BinningTrueMom[c0+1]<<"GeV, "<<BinningTrueAngle[c1]<<"-"<<BinningTrueAngle[c1+1]<<"°), Fitted=("<<Gaussian1D_RecMom[c0][c1]->GetParameter(1)-Gaussian1D_RecMom[c0][c1]->GetParameter(2)<<"-"<<Gaussian1D_RecMom[c0][c1]->GetParameter(1)+Gaussian1D_RecMom[c0][c1]->GetParameter(2)<<"GeV, "<<Gaussian1D_RecAngle[c0][c1]->GetParameter(1)-Gaussian1D_RecAngle[c0][c1]->GetParameter(2)<<"-"<<Gaussian1D_RecAngle[c0][c1]->GetParameter(1)+Gaussian1D_RecAngle[c0][c1]->GetParameter(2)<<"°)"<<endl;
	cout<<setprecision(1)<<"True=("<<BinningTrueMom[c0]<<"-"<<BinningTrueMom[c0+1]<<"GeV, "<<BinningTrueAngle[c1]<<"-"<<BinningTrueAngle[c1+1]<<"°), Fitted=("<<Gaussian1D_RecMom[c0][c1]->GetParameter(1)<<"cm, "<<Gaussian1D_RecAngle[c0][c1]->GetParameter(1)<<"°)"<<", check chi2/ndf="<<(Gaussian1D_RecMom[c0][c1]->GetChisquare()/Gaussian1D_RecMom[c0][c1]->GetNDF())<<(Gaussian1D_RecMom[c0][c1]->GetChisquare()/Gaussian1D_RecMom[c0][c1]->GetNDF())<<endl;
      }
    }
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	double Diagonal=0;double NonDiagonal=0;
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    int bintrue=c0*NBinsTrueAngle+c1;
	    int binrec=e0*NBinsRecAngle+e1;
	    //cout<<bintrue<<", "<<binrec<<endl;
	    double unfolding=(*MUnfolding)(bintrue,binrec);
	    double centralrecmom=BinningRecMom[e0]+(BinningRecMom[e0+1]-BinningRecMom[e0])/2;//Central bin value in dmu
	    double centralrecangle=BinningRecAngle[e1]+(BinningRecAngle[e1+1]-BinningRecAngle[e1])/2;//Central bin value in thetamu
	    
 	    //if( TMath::Abs(centralrecmom-Gaussian1D_RecMom[c0][c1]->GetParameter(1))<Gaussian1D_RecMom[c0][c1]->GetParameter(2) && TMath::Abs(centralrecangle-Gaussian1D_RecAngle[c0][c1]->GetParameter(1))<Gaussian1D_RecAngle[c0][c1]->GetParameter(2)) Diagonal+=unfolding*(*VReconstructedEvents)[binrec];
	    /*
	    //We want: the central bin value to be within max+-1sigma
	    if(TMath::Abs(centralrecmom-Gaussian1D_RecMom[c0][c1]->GetParameter(1))<Gaussian1D_RecMom[c0][c1]->GetParameter(2) && TMath::Abs(centralrecangle-Gaussian1D_RecAngle[c0][c1]->GetParameter(1))<Gaussian1D_RecAngle[c0][c1]->GetParameter(2)){
	      Diagonal+=unfolding*(*VReconstructedEvents)[binrec]; //the condition is just: if the central bin value to be within max+-1sigma , the bin is diagonal
	    }
	    else NonDiagonal+=unfolding*(*VReconstructedEvents)[binrec];
	  }
	}
	*/
	//different condition here than above: if any part of the bin is within mean gaus+-sigma, then this bin is diagonal.
	    if( ((centralrecmom<Gaussian1D_RecMom[c0][c1]->GetParameter(1) && TMath::Abs(Gaussian1D_RecMom[c0][c1]->GetParameter(1)-BinningRecMom[e0+1])<Gaussian1D_RecMom[c0][c1]->GetParameter(2)) || (centralrecmom>Gaussian1D_RecMom[c0][c1]->GetParameter(1) && TMath::Abs(Gaussian1D_RecMom[c0][c1]->GetParameter(1)-BinningRecMom[e0])<Gaussian1D_RecMom[c0][c1]->GetParameter(2))) && ((centralrecangle<Gaussian1D_RecAngle[c0][c1]->GetParameter(1) && TMath::Abs(Gaussian1D_RecAngle[c0][c1]->GetParameter(1)-BinningRecAngle[e1+1])<Gaussian1D_RecAngle[c0][c1]->GetParameter(2)) || (centralrecangle>Gaussian1D_RecAngle[c0][c1]->GetParameter(1) && TMath::Abs(Gaussian1D_RecAngle[c0][c1]->GetParameter(1)-BinningRecAngle[e1])<Gaussian1D_RecAngle[c0][c1]->GetParameter(2)))) Diagonal+=unfolding*(*VReconstructedEvents)[binrec]; //the condition is just: if any part of the bin is within mean gaus+-sigma, then this bin is diagonal.
	    else NonDiagonal+=unfolding*(*VReconstructedEvents)[binrec];
	  }
	  }
	
	cout<<"Bin ("<<c0<<","<<c1<<"), ratio="<<NonDiagonal/(Diagonal)<<", Diagonal="<<Diagonal<<", Non diagonal="<<NonDiagonal<<endl; 
      }
    }


    ///////////////////////////////////////////////////////////////////////////////////////


    
    //canCorrelationStat->Write();
    NumberOfEvents->Write();
    NumberOfEvents_TrueRec->Write();
    BiasDistribution->Write();
    ErrorDistribution->Write();
    cUnfolding->Write();
   //ofile->Write();
    ofile->Close();
#endif    
    //theApp.Run();
  
  return 0;
}

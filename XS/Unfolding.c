
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
#define DEBUG


int main(int argc, char ** argv){

  TApplication theApp("App",0,0);

  double vLikelihood[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double vUnfolding[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double vInitialPriorMC[NBinsTrueMom][NBinsTrueAngle];
  double vInitialPrior[NBinsTrueMom][NBinsTrueAngle];
  double vPrior[NBinsTrueMom][NBinsTrueAngle];
  double vPriorNormalised[NBinsTrueMom][NBinsTrueAngle];
  double vPosterior[NBinsTrueMom][NBinsTrueAngle];


  bool PriorMC=true;
  int NIterations=2;

  int c=-1;
  int NPriors=0;
  char * fDataName = new char[256];
  char * fMCName = new char[256];
  sprintf(fDataName,"DistributionsData.root");
  sprintf(fMCName,"DistributionsMC.root");

  
  while ((c = getopt(argc, argv, "d:m:")) != -1) {
    switch(c){
    case 'd':
      fDataName=optarg;
      break; 
   case 'm':
      fMCName=optarg;
      break;
    }
  }
  
  cout<<"looking for files named: "<<fDataName<<" and "<<fMCName<<endl;
  Xsec * XS = new Xsec();

  double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double DataReconstructedEvents[NBinsRecMom][NBinsRecAngle], MCReconstructedBkgEvents[NBinsRecMom][NBinsRecAngle], MCReconstructedEvents[NBinsRecMom][NBinsRecAngle];
  double MCEfficiency[NBinsTrueMom][NBinsTrueAngle];
  double TrueEventsDistribution[NBinsTrueMom][NBinsTrueAngle];

  //0. Load the input distributions
  XS->Xsec::LoadInputFiles(fDataName,fMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents, MCEfficiency);

  //1. Build the likelihood matrix & intial prior from the MC
  XS->Xsec::BuildLikelihood(vLikelihood, vInitialPriorMC, TrueEventsDistribution, MCReconstructedEvents_TrueSignal);

  for(int it=0;it<NIterations;it++){

    //2. Build the prior used in the unfolding
    XS->Xsec::SetPrior(vPriorNormalised,vPrior,vInitialPriorMC,vInitialPrior,vPosterior,PriorMC,it);

    //3. Build the unfolding matrix
    XS->Xsec::BuildUnfolding(vUnfolding,vLikelihood,vPriorNormalised);

    //4. Apply the unfolding on data
    XS->Xsec::ApplyUnfolding(vPosterior, vUnfolding, DataReconstructedEvents, MCReconstructedBkgEvents);

  }

  //5. Correct the effects of detector efficiency, flux and bin width to obtain the cross section

  theApp.Run();

  return 0;
}

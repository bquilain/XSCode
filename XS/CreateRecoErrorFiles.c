// ML 2017/11/27
//   This macro is used to mimic a +1sigma and -1 sigma global variation of the number of events (normalization error)
//   Input: the nominal (Syst0_0) output of CC0piSelection
//   Output: a .txt file with scaled number of events, efficiency etc --> to feed the Unfolding macro


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
#include "Xsec.cc"

int main(int argc, char ** argv){

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  bool PM=true;  
  bool MC=true;
  char * OutputName = new char[256];

  double RelativeError=0.;
  bool UsersChoice=false;

  int c=-1;
  while ((c = getopt(argc, argv, "wpv:")) != -1) {
    switch(c){
    case 'p':
      PM=true;
      break;      
    case 'w':
      PM=false;
      break;      
    case 'v':
      RelativeError=atof(optarg)/100.;//enter a value in %
      UsersChoice=true;
    }
  }
  

  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  cout<<"Selected detector is "<<DetName<<endl;
  Xsec * XS = new Xsec(PM);
  XS->Xsec::Initialize();

  if(!UsersChoice) RelativeError=(PM? 1.32:1.34)*0.01; // ML update 2017/11/20

  char * txtMCName = new char[512];
  sprintf(txtMCName, "%s/XS/files/MCSelected_%s%s_Systematics0_0.txt",cINSTALLREPOSITORY,DetName,(!PM?"_merged":""));

  char * txtOutName = new char[512];
  ofstream fEvent;
  
  // input is in reconstructed variables
  const int NBinsMom=NBinsRecMom;
  const int NBinsAngle=NBinsRecAngle;
  double BinningMom[NBinsMom+1];
  double BinningAngle[NBinsAngle+1];
  for(int i=0;i<NBinsMom+1;i++) BinningMom[i]=BinningRecMom[i];
  for(int i=0;i<NBinsAngle+1;i++) BinningAngle[i]=BinningRecAngle[i];
  //

  double MCReconstructedEvents[NBinsMom][NBinsAngle];//whether rec, whether unfolded
  double MCReconstructedBkgEvents[NBinsMom][NBinsAngle];//whether rec, whether unfolded
  double Efficiency[NBinsTrueMom][NBinsTrueAngle];
  double * NumberOfPOT = new double();
  double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsMom][NBinsAngle];

  XS->Xsec::LoadInputFiles(txtMCName,txtMCName,MCReconstructedEvents_TrueSignal,MCReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,Efficiency,NumberOfPOT);

  double MCReconstructedEvents_out[NBinsMom][NBinsAngle]={{0.}};//whether rec, whether unfolded
  double MCReconstructedBkgEvents_out[NBinsMom][NBinsAngle]={{0.}};//whether rec, whether unfolded
  double Efficiency_out[NBinsTrueMom][NBinsTrueAngle]={{0.}};
  double MCReconstructedEvents_TrueSignal_out[NBinsTrueMom][NBinsTrueAngle][NBinsMom][NBinsAngle]={{{{0.}}}};

  double sigma=0.;

  for(int var=0;var<2;var++){
    if(var==0) sigma=-RelativeError;
    else sigma=+RelativeError;

    // Scale the data
    cout<<"Now scaling all events by "<<sigma*100.<<" %"<<endl;
    for(int e0=0;e0<NBinsMom;e0++){
      for(int e1=0;e1<NBinsAngle;e1++){
	if(MCReconstructedEvents[e0][e1]>0) MCReconstructedEvents_out[e0][e1]=MCReconstructedEvents[e0][e1]*(1+sigma);
	if(MCReconstructedBkgEvents[e0][e1]>0) MCReconstructedBkgEvents_out[e0][e1]=MCReconstructedBkgEvents[e0][e1]*(1+sigma);
	for(int c0=0;c0<NBinsTrueMom;c0++){
	  for(int c1=0;c1<NBinsTrueAngle;c1++){
	    if(MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]>0) MCReconstructedEvents_TrueSignal_out[c0][c1][e0][e1]=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]*(1+sigma);
	  }
	}
      }
    }
    for(int c0=0;c0<NBinsTrueMom;c0++){
      for(int c1=0;c1<NBinsTrueAngle;c1++){
	if(Efficiency[c0][c1]>0) Efficiency_out[c0][c1]=Efficiency[c0][c1]*(1+sigma);
      }
    }


    // Write the output txt file
    sprintf(txtOutName, "%s/XS/files/MCSelected_%s%s_SystematicsScaled%1.2f_%d.txt",cINSTALLREPOSITORY,DetName,(!PM?"_merged":""),100.*RelativeError,var);
    cout<<"opening file "<<txtOutName<<" to write output"<<endl;
    fEvent.open(txtOutName,ios::out);
    
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	fEvent<<MCReconstructedEvents_out[e0][e1]<<" ";	
      }
    }
    fEvent<<666<<" ";
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    fEvent<<MCReconstructedEvents_TrueSignal_out[c0][c1][e0][e1]<<" ";
	  }
	}
      }
    }
    fEvent<<666<<" ";
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	fEvent<<MCReconstructedBkgEvents_out[e0][e1]<<" ";
      }
    }
    fEvent<<666<<" ";
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	fEvent<<Efficiency_out[c0][c1]<<" ";
      }
    }

    fEvent<<666;
    fEvent<<" "<<(*NumberOfPOT)<<" ";
    fEvent<<666;
    
    fEvent.close();
    cout<<"file "<<txtOutName<<" is closed"<<endl;
    
  }

}

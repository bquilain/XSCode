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
#include "Reconstruction.cc"
#include "Xsec.cc"
#define DEBUG

int main(int argc, char ** argv){

  //TApplication theApp("App",0,0);
  gStyle->SetOptFit(kTRUE);
  Xsec * XS = new Xsec();
  XS->Xsec::Initialize();
  
  double vLikelihood[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double vUnfolding[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double vInitialPriorMC[NBinsTrueMom][NBinsTrueAngle];
  double vInitialPrior[NBinsTrueMom][NBinsTrueAngle];
  double vPrior[NBinsTrueMom][NBinsTrueAngle];
  double vPriorNormalised[NBinsTrueMom][NBinsTrueAngle];
  double vPosterior[NBinsTrueMom][NBinsTrueAngle];

  int c=-1;
  bool Systematics=false;
 char * InputNameMC=new char[256];char * InputNameData=new char[256];
 
  while ((c = getopt(argc, argv, "fxdm:")) != -1) {
    switch(c){
    case 'm':
      InputNameMC=optarg;
      break;
    case 's':
      Systematics=true;
      break; 
    }
  }

  //00.Output file preparation
  TFile * ofile = new TFile("plots/UnfoldingPlots.root","recreate");
  
  //0. Preparation for covariance estimation
  int NBins=(NBinsTrueMom-1)*(NBinsTrueAngle-1);
  TH2D * Covariance = new TH2D("Covariance","Covariance matrix",(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1),(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1));
  TH2D * Correlation = new TH2D("Correlation","Correlation matrix",(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1),(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1));
  TH2D * NVariations = new TH2D("NVariations","",(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1),(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1));

  TH2D * NumberOfEvents = new TH2D("NumberOfEvents","Number of events (true MC information) for each true bin",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * NumberOfEvents_TrueRec = new TH2D("NumberOfEvents_TrueRec","Number of events (true MC information) for each true & reconstructed bin",NBinsRecMom*NBinsRecAngle,0,NBinsRecMom*NBinsRecAngle,(NBinsTrueMom-1)*(NBinsTrueAngle-1),0,(NBinsTrueMom-1)*(NBinsTrueAngle-1));

  for(int c0=0;c0<NBinsTrueMom-1;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle-1;c1++){//loop over cause 1
      NumberOfEvents->SetBinContent(c0+1,c1+1,0);
      for(int d0=0;d0<NBinsTrueMom-1;d0++){//loop over cause 0
	for(int d1=0;d1<NBinsTrueAngle-1;d1++){//loop over cause 1
	  int BinX=c0*(NBinsTrueAngle-1)+c1+1;
	  int BinY=d0*(NBinsTrueAngle-1)+d1+1;
	  Covariance->SetBinContent(BinX,BinY,0);
	  Correlation->SetBinContent(BinX,BinY,0);
	  NVariations->SetBinContent(BinX,BinY,0);
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
  
  
  ifstream fEvent_Default;
  ifstream fEvent_Varied;
  double Value_Default[NBinsTrueMom][NBinsTrueAngle];
  double Value_Default_Data[NBinsTrueMom][NBinsTrueAngle];
  double Sigma[NBinsTrueMom][NBinsTrueAngle];
  double Value_Varied[NBinsTrueMom][NBinsTrueAngle];
  double Value_Varied_Data[NBinsTrueMom][NBinsTrueAngle];
  double MaximalDifference[NBinsTrueMom][NBinsTrueAngle];
  int check=0;

  //1. Load the default Unfolding & intialize bin histograms
  TChain * wtree = new TChain("wtree");
  wtree->Add(InputNameMC);
  for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
    for(int n=0;n<NE[ErrorType];n++){
      double ErrorValue=Start[ErrorType]+n*Step[ErrorType];
      wtree->Add("");
      //plots/PriorMC_Sample5.root
    }
  }
    
  int Iterations, Systematics, Statistics;
  double Events[NBinsTrueMom][NBinsTrueAngle];
  double EventsAll[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double TrueEvents[NBinsTrueMom][NBinsTrueAngle];

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

  TBranch* Br_EventsAll = wtree->GetBranch("EventsAll");
  Br_EventsAll->SetAddress(EventsAll);
  wtree->SetBranchAddress("EventsAll",EventsAll);  

  TBranch* Br_TrueEvents = wtree->GetBranch("TrueEvents");
  Br_TrueEvents->SetAddress(TrueEvents);
  wtree->SetBranchAddress("TrueEvents",TrueEvents);  
////////////////////////////////////////////////////////////////////////////////
  
  int nevt=wtree->GetEntries();

  TCanvas * canEvents[NBinsTrueMom][NBinsTrueAngle];
  TCanvas * canEventsXIterations[NBinsTrueMom][NBinsTrueAngle];

  //2D true variable distributions
  TH2D * UnfoldedDistribution = new TH2D("UnfoldedDistribution","Events distribution in true binning after unfolding applied",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * TrueDistribution = new TH2D("TrueDistribution","Events distribution in true binning for the true MC information",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * BiasDistribution = new TH2D("BiasDistribution","Bias=Mean of the gaussian we used to fit the pull distribution, given in true variables",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * ErrorDistribution = new TH2D("ErrorDistribution","Error=width of the gaussian we used to fit the pull distribution, given in true variables",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
  
  //Events pull distributions, only for stat variations!!!!.
  TH1D * pullEvents[NBinsTrueMom][NBinsTrueAngle];

  //Events pull distributions with the number of iterations.
  TH2D * pullEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullEventsmeanXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullEventssigmaXIterations[NBinsTrueMom][NBinsTrueAngle];

  //Events pull distributions only on stat variations with the number of iterations.
  TH2D * pullstatEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventsmeanXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventssigmaXIterations[NBinsTrueMom][NBinsTrueAngle];

  //Events pull distributions only on syst variations with the number of iterations.
  TH2D * pullsystEventsXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullsystEventsmeanXIterations[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullsystEventssigmaXIterations[NBinsTrueMom][NBinsTrueAngle];

  //Events pull distributions with the systematic varied toy experiments.
  TH2D * pullstatEventsXsystEvents[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventsmeanXsystEvents[NBinsTrueMom][NBinsTrueAngle];
  TH1D * pullstatEventssigmaXsystEvents[NBinsTrueMom][NBinsTrueAngle];

  
  TH1D * hTemp; TF1 * gausTemp = new TF1("gausTemp","gaus",-3,3);
  int NBinsIteration=2;
  
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      pullEvents[c0][c1] = new TH1D(Form("pullEvents[%d][%d]",c0,c1),"",100,-2,2);

      pullEventsXIterations[c0][c1] = new TH2D(Form("pullEventsXIterations[%d][%d]",c0,c1),"Pull distribution X unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration,100,-2,2);
      pullEventsmeanXIterations[c0][c1] = new TH1D(Form("pullEventsmeanXIterations[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration);
      pullEventssigmaXIterations[c0][c1] = new TH1D(Form("pullEventssigmaXIterations[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X unfolding iterations, for true binning",NBinsIteration,0,NBinsIteration);

      pullstatEventsXIterations[c0][c1] = new TH2D(Form("pullstatEventsXIterations[%d][%d]",c0,c1),"Pull distribution X unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration,100,-2,2);
      pullstatEventsmeanXIterations[c0][c1] = new TH1D(Form("pullstatEventsmeanXIterations[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);
      pullstatEventssigmaXIterations[c0][c1] = new TH1D(Form("pullstatEventssigmaXIterations[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X unfolding iterations, for true binning (Stat. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);

      pullsystEventsXIterations[c0][c1] = new TH2D(Form("pullsystEventsXIterations[%d][%d]",c0,c1),"Pull distribution X unfolding iterations, for true binning (Syst. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration,100,-2,2);
      pullsystEventsmeanXIterations[c0][c1] = new TH1D(Form("pullsystEventsmeanXIterations[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X unfolding iterations, for true binning (Syst. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);
      pullsystEventssigmaXIterations[c0][c1] = new TH1D(Form("pullsystEventssigmaXIterations[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X unfolding iterations, for true binning (Syst. error variations only between toy experiment to build pull)",NBinsIteration,0,NBinsIteration);      


      pullstatEventsXsystEvents[c0][c1] = new TH2D(Form("pullstatEventsXsystEvents[%d][%d]",c0,c1),"Pull distribution X systematically varied toy experiment, for true binning (Stat. error variations only between toy experiment to build pull)",100,0,100,100,-2,2);
      pullstatEventsmeanXsystEvents[c0][c1] = new TH1D(Form("pullstatEventsmeanXsystEvents[%d][%d]",c0,c1),"Mean (gaussian fit) pull distribution X systematically varied toy experiment, for true binning (Stat. error variations only between toy experiment to build pull)",100,0,100);
      pullstatEventssigmaXsystEvents[c0][c1] = new TH1D(Form("pullstatEventssigmaXsystEvents[%d][%d]",c0,c1),"Error (gaussian fit) pull distribution X systematically varied toy experiment, for true binning (Stat. error variations only between toy experiment to build pull)",100,0,100);
																								   
    }
  }
  
  for(int ievt=0;ievt<nevt;ievt++){
    wtree->GetEntry(ievt);
    //cout<<"Statistics="<<Statistics<<endl;
    for(int c0=0;c0<NBinsTrueMom-1;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle-1;c1++){//loop over cause 1
	UnfoldedDistribution->SetBinContent(c0+1,c1+1,Events[c0][c1]);
	TrueDistribution->SetBinContent(c0+1,c1+1,TrueEvents[c0][c1]);
	double Value=Events[c0][c1]-TrueEvents[c0][c1];
	if(TrueEvents[c0][c1]!=0) Value/=TrueEvents[c0][c1];
	pullEventsXIterations[c0][c1]->Fill(Iterations,Value);
	pullstatEventsXsystEvents[c0][c1]->Fill(Systematics,Value);
	
	if(Statistics==0) pullsystEventsXIterations[c0][c1]->Fill(Iterations,Value);
	if(Systematics==0){
	  pullstatEventsXIterations[c0][c1]->Fill(Iterations,Value);
	  pullEvents[c0][c1]->Fill(Value);
	  for(int d0=0;d0<NBinsTrueMom-1;d0++){//loop over dause 0
	    for(int d1=0;d1<NBinsTrueAngle-1;d1++){//loop over cause 1
	      int BinX=c0*(NBinsTrueAngle-1)+c1+1;
	      int BinY=d0*(NBinsTrueAngle-1)+d1+1;
	      double Cova=Covariance->GetBinContent(BinX,BinY);
	      Cova+=(Events[c0][c1]-TrueEvents[c0][c1])*(Events[d0][d1]-TrueEvents[d0][d1]);
	      Covariance->SetBinContent(BinX,BinY,Cova);

	      int nvar=NVariations->GetBinContent(BinX,BinY);
	      nvar++;
	      NVariations->SetBinContent(BinX,BinY,1);;
	    }
	  }
	  if(ievt==0){//Only for the non-varied distribution. We do not want to pile-up the information of all the varied distributions!
	    double nEvts=NumberOfEvents->GetBinContent(c0+1,c1+1);
	    for(int e0=0;e0<NBinsRecMom-1;e0++){//loop over effect 0
	      for(int e1=0;e1<NBinsRecAngle-1;e1++){//loop over effect 1
		int BinX=e0*(NBinsRecAngle-1)+e1+1;
		int BinY=c0*(NBinsTrueAngle-1)+c1+1;
		double nEvts_TrueRec=NumberOfEvents_TrueRec->GetBinContent(BinX,BinY);
		nEvts_TrueRec+=EventsAll[c0][c1][e0][e1];
		nEvts+=EventsAll[c0][c1][e0][e1];
		NumberOfEvents_TrueRec->SetBinContent(BinX,BinY,nEvts_TrueRec);
	      }
	    }
	    cout<<"nevts="<<nEvts<<endl;
	    NumberOfEvents->SetBinContent(c0+1,c1+1,nEvts);
	  }
	}
      }
    }
  }

  //4. For the total error, build the correlation matrix from the covariant one
  Covariance->Divide(NVariations);
  for(int c0=0;c0<NBinsTrueMom-1;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle-1;c1++){//loop over cause 1
      for(int d0=0;d0<NBinsTrueMom-1;d0++){//loop over cause 0
	for(int d1=0;d1<NBinsTrueAngle-1;d1++){//loop over cause 1
	  int BinX=c0*(NBinsTrueAngle-1)+c1+1;
	  int BinY=d0*(NBinsTrueAngle-1)+d1+1;
	  double Corr=Covariance->GetBinContent(BinX,BinY);
	  double NormCorr=TMath::Sqrt(Covariance->GetBinContent(BinX,BinX)*Covariance->GetBinContent(BinY,BinY));
	  if(NormCorr) Corr/=NormCorr;
	  Correlation->SetBinContent(BinX,BinY,Corr);
	}
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  
	      
  ofile->cd();
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

  TCanvas * canCorrelations = new TCanvas("cCorrelations");
  gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(kFALSE);
  Correlation->Draw("colztext");
  
  TLine * lBins[NBinsTrueMom-1][2];
  TLatex * tBins[NBinsTrueMom-1][2];
  cout.precision(2);
  cout.setf(ios::fixed);
  
  for(int imom=0;imom<NBinsTrueMom-1;imom++){
    cout<<"imom="<<imom<<endl;
    lBins[imom][0] = new TLine(imom*(NBinsTrueAngle-1),0,imom*(NBinsTrueAngle-1),(NBinsTrueMom-1)*(NBinsTrueAngle-1));
    lBins[imom][1] = new TLine(0,imom*(NBinsTrueAngle-1),(NBinsTrueMom-1)*(NBinsTrueAngle-1),imom*(NBinsTrueAngle-1));

    tBins[imom][0] = new TLatex(imom*(NBinsTrueAngle-1)+NBinsTrueAngle/4,-1,Form("Bin%d p_{#mu}",imom+1));
    tBins[imom][1] = new TLatex(-1,imom*(NBinsTrueAngle-1)+NBinsTrueAngle/4,Form("Bin%d p_{#mu}",imom+1));
    tBins[imom][1]->SetTextAngle(90);
    tBins[imom][0]->SetTextSize(0.04);
    tBins[imom][1]->SetTextSize(0.04);

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



  
    for(int c0=0;c0<NBinsTrueMom-1;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle-1;c1++){//loop over cause 1
	//Draw the pull distribution for all the variations in the input file
	canEvents[c0][c1] = new TCanvas(Form("cPullDistribution_%d_%d",c0,c1),Form("True distribution of events after unfolding in the bin %d_%d",c0,c1));
	pullEvents[c0][c1]->GetXaxis()->SetTitle("Pull");
	pullEvents[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	pullEvents[c0][c1]->Draw();
	pullEvents[c0][c1]->Fit("gaus","RQ");
	canEvents[c0][c1]->Write();

	//Draw the pull distribution with iterations & behaviour of the mean/error of the pull, for all the variations in the input file
       
	canEventsXIterations[c0][c1] = new TCanvas(Form("cPullXIterations%d_%d",c0,c1),Form("Pull maximum and 1#sigma error of events after unfolding in the bin %d_%d",c0,c1));
	pullEventsXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding iterations");
	pullEventsXIterations[c0][c1]->GetYaxis()->SetTitle("Pull");
	for(int ibinx=1;ibinx<=pullEventsXIterations[c0][c1]->GetNbinsX();ibinx++){
	  gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	  pullEventsXIterations[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	  pullEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(1));
	  pullEventsmeanXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(1));
	  pullEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(2));
	  pullEventssigmaXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(2));
	}
	pullEventsmeanXIterations[c0][c1]->SetLineColor(1);
	pullEventssigmaXIterations[c0][c1]->SetLineColor(kGray);
	pullEventsmeanXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding Iterations");
	//pullEventsmeanXIterations[c0][c1]->Draw("E1");
	//pullEventssigmaXIterations[c0][c1]->Draw("E1same");
	
	//Draw the pull distribution with iterations & behaviour of the mean/error of the pull, BUT ONLY FOR STATISTICAL VARIATIONS in the input file	
	pullstatEventsXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding iterations");
	pullstatEventsXIterations[c0][c1]->GetYaxis()->SetTitle("Pullstat");
	for(int ibinx=1;ibinx<=pullstatEventsXIterations[c0][c1]->GetNbinsX();ibinx++){
	  gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	  pullstatEventsXIterations[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	  pullstatEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(1));
	  pullstatEventsmeanXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(1));
	  pullstatEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(2));
	  pullstatEventssigmaXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(2));
	  //cout<<"center="<<gausTemp->GetParameter(1)<<", error"<<gausTemp->GetParameter(2)<<endl;

	  if(ibinx==pullsystEventsXIterations[c0][c1]->GetNbinsX()){//only for the last bin of the number of iteration. In the case that there is no variations on the number of iterations, this is equal to the number of iteration set.
	    BiasDistribution->SetBinContent(c0+1,c1+1,gausTemp->GetParameter(1));
	    ErrorDistribution->SetBinContent(c0+1,c1+1,gausTemp->GetParameter(2));
	  }
	}
	pullstatEventsmeanXIterations[c0][c1]->SetLineColor(kRed);
	pullstatEventssigmaXIterations[c0][c1]->SetLineColor(kOrange);
	pullstatEventsmeanXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding Iterations");
	pullstatEventsmeanXIterations[c0][c1]->Draw("E1");
	pullstatEventssigmaXIterations[c0][c1]->Draw("E1same");
	
	//Draw the pull distribution with iterations & behaviour of the mean/error of the pull, BUT ONLY FOR SYSTEMATICS VARIATIONS in the input file	
	pullsystEventsXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding iterations");
	pullsystEventsXIterations[c0][c1]->GetYaxis()->SetTitle("Pullsyst");
	for(int ibinx=1;ibinx<=pullsystEventsXIterations[c0][c1]->GetNbinsX();ibinx++){
	  gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	  pullsystEventsXIterations[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	  pullsystEventsmeanXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(1));
	  pullsystEventsmeanXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(1));
	  pullsystEventssigmaXIterations[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(2));
	  pullsystEventssigmaXIterations[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(2));

	}
	pullsystEventsmeanXIterations[c0][c1]->SetLineColor(kBlue);
	pullsystEventssigmaXIterations[c0][c1]->SetLineColor(kCyan);
	pullsystEventsmeanXIterations[c0][c1]->GetXaxis()->SetTitle("Unfolding Iterations");
	//pullsystEventsmeanXIterations[c0][c1]->Draw("E1same");
	//pullsystEventssigmaXIterations[c0][c1]->Draw("E1same");


	  //Draw the pull distribution with systematic bias & behaviour of the mean/error of the pull. PUll is only done FOR STATISTICAL VARIATIONS in the input file, but x axis is the systematic variation!!	
	pullstatEventsXsystEvents[c0][c1]->GetXaxis()->SetTitle("Systematically varied exp.");
	pullstatEventsXsystEvents[c0][c1]->GetYaxis()->SetTitle("Pullstat");
	for(int ibinx=1;ibinx<=pullstatEventsXsystEvents[c0][c1]->GetNbinsX();ibinx++){
	  gausTemp->SetParameter(0,1.);gausTemp->SetParameter(0,0.);gausTemp->SetParameter(0,1.);
	  pullstatEventsXsystEvents[c0][c1]->FitSlicesY(gausTemp,ibinx,ibinx,0,"RQ0");
	  pullstatEventsmeanXsystEvents[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(1));
	  pullstatEventsmeanXsystEvents[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(1));
	  pullstatEventssigmaXsystEvents[c0][c1]->SetBinContent(ibinx,gausTemp->GetParameter(2));
	  pullstatEventssigmaXsystEvents[c0][c1]->SetBinError(ibinx,gausTemp->GetParError(2));
	  //cout<<"center="<<gausTemp->GetParameter(1)<<", error"<<gausTemp->GetParameter(2)<<endl;
	}
	//leg->Draw("lsame");
	canEventsXIterations[c0][c1]->Write();
	pullEventsXIterations[c0][c1]->Write();
      }
    }
    for(int c0=0;c0<NBinsTrueMom-1;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle-1;c1++){//loop over cause 1  
	pullstatEventsXsystEvents[c0][c1]->Write();
	pullstatEventsmeanXsystEvents[c0][c1]->Write();	
	pullstatEventssigmaXsystEvents[c0][c1]->Write();	
      }
    }
    
    TCanvas * canNumberOfEvents = new TCanvas("cNumberOfEvents_True");
    gStyle->SetPaintTextFormat("2.2f");
    gStyle->SetOptStat(kFALSE);
    NumberOfEvents->Draw("colztext");
    canNumberOfEvents->Write();
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

    

    /////////////////Here we treat the optimisation of the true binning////////////////////
    TMatrixD * MUnfolding = (TMatrixD*) file->Get("MUnfolding");
    TMatrixD * MLikelihood = (TMatrixD*) file->Get("MLikelihood");
    TVectorD * VReconstructedEvents = (TVectorD*) file->Get("VReconstructedEvents");
    /*
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    LikelihoodMatrix->SetBinContent();
	    
	    */
    
    TCanvas * cLikelihood = new TCanvas("cLikelihood");
    MLikelihood->Draw("colztext");
    cLikelihood->Write();
    
    TCanvas * cUnfolding = new TCanvas("cUnfolding");
    MUnfolding->Draw("colztext");
    cUnfolding->Write();
    
    
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
    /*
    TH2D * ReconstructedLikelihood[NBinsTrueMom][NBinsTrueAngle];
    TF2 * Gaussian2D[NBinsTrueMom][NBinsTrueAngle];
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	ReconstructedLikelihood[c0][c1]= new TH2D(Form("ReconstructedLikelihood[%d][%d]",c0,c1),"",NBinsRecMom,BinningRecMom,NBinsRecAngle,BinningRecAngle);
	double minrecmom=BinningRecMom[0];
	double maxrecmom=BinningRecMom[NBinsRecMom];
	double minrecangle=BinningRecAngle[0];
	double maxrecangle=BinningRecAngle[NBinsRecAngle];
	
	Gaussian2D[c0][c1] = new TF2(Form("Gaussian2D[%d][%d]",c0,c1),"[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",minrecmom,maxrecmom,minrecangle,maxrecangle);
	Gaussian2D[c0][c1]->SetParameters(1,10,10,10,10);
	
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
	ReconstructedLikelihood[c0][c1]->Fit(Gaussian2D[c0][c1],"R0");
	ReconstructedLikelihood[c0][c1]->Write(Form("ReconstructedLikelihood[%d][%d]",c0,c1));
      }
    }
   
    TH1D * ReconstructedLikelihood_RecMom[NBinsTrueMom][NBinsTrueAngle];
    TH1D * ReconstructedLikelihood_RecAngle[NBinsTrueMom][NBinsTrueAngle];
    TF1 * Gaussian1D_RecMom[NBinsTrueMom][NBinsTrueAngle];
    TF1 * Gaussian1D_RecAngle[NBinsTrueMom][NBinsTrueAngle];
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	ReconstructedLikelihood_RecMom[c0][c1]= new TH1D(Form("ReconstructedLikelihood_RecMom[%d][%d]",c0,c1),"",NBinsRecMom,BinningRecMom);
	ReconstructedLikelihood_RecAngle[c0][c1]= new TH1D(Form("ReconstructedLikelihood_RecAngle[%d][%d]",c0,c1),"",NBinsRecAngle,BinningRecAngle);
	
	double minrecmom=BinningRecMom[0];
	double maxrecmom=BinningRecMom[NBinsRecMom];
	double minrecangle=BinningRecAngle[0];
	double maxrecangle=BinningRecAngle[NBinsRecAngle];
	
	Gaussian1D_RecMom[c0][c1] = new TF1(Form("Gaussian1D_RecMom[%d][%d]",c0,c1),"[0]*TMath::Gaus(x,[1],[2])",minrecmom,maxrecmom);
	Gaussian1D_RecMom[c0][c1]->SetParameters(1,10,10);
	Gaussian1D_RecAngle[c0][c1] = new TF1(Form("Gaussian1D_RecAngle[%d][%d]",c0,c1),"[0]*TMath::Gaus(x,[1],[2])",minrecangle,maxrecangle);
	Gaussian1D_RecAngle[c0][c1]->SetParameters(1,10,10);
	
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    int bintrue=c0*NBinsTrueAngle+c1;
	    int binrec=e0*NBinsRecAngle+e1;
	    //cout<<bintrue<<", "<<binrec<<endl;
	    //cout<<(*MLikelihood)(bintrue,binrec)<<endl;
	    //cout<<(*MLikelihood)[bintrue][binrec]<<endl;
	    double correlation=(*MLikelihood)(bintrue,binrec);
	    //cout<<setprecision(5)<<"Bin Rec Mom="<<e0+1<<", Rec Angle="<<e1+1<<", Correlation="<<correlation<<endl;
	    ReconstructedLikelihood_RecMom[c0][c1]->SetBinContent(e0+1,correlation);
	    ReconstructedLikelihood_RecAngle[c0][c1]->SetBinContent(e1+1,correlation);
	  }
	}
	ReconstructedLikelihood_RecMom[c0][c1]->Fit(Gaussian1D_RecMom[c0][c1],"R0");
	ReconstructedLikelihood_RecMom[c0][c1]->Write(Form("ReconstructedLikelihood_RecMom[%d][%d]",c0,c1));
	ReconstructedLikelihood_RecAngle[c0][c1]->Fit(Gaussian1D_RecAngle[c0][c1],"R0");
	ReconstructedLikelihood_RecAngle[c0][c1]->Write(Form("ReconstructedLikelihood_RecAngle[%d][%d]",c0,c1));
      }
    }
*/

    ///////////////////////////////////////////////////////////////////////////////////////


    
    canCorrelations->Write();
    NumberOfEvents->Write();
    NumberOfEvents_TrueRec->Write();
    BiasDistribution->Write();
    ErrorDistribution->Write();
    cUnfolding->Write();
   //ofile->Write();
    ofile->Close();
    
    //theApp.Run();
  
  return 0;
}

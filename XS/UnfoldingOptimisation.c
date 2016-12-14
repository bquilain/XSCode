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



  bool PriorMC=true;
  int c=-1;
  int NPriors=0;
  char * fDataName = new char[256];
  char * fMCName = new char[256];
  sprintf(fDataName,"DistributionsData.root");
  sprintf(fMCName,"DistributionsMC.root");
  char * OutputName = new char[256];
  int NIterations=1;
  int NToys=1;
  bool StatisticalFluctuations=false;
  bool SystematicFluctuations=false;
  bool NIterationFluctuations=false;
  int NTarget=8.299e28;
  
  while ((c = getopt(argc, argv, "d:m:n:t:syio:a:")) != -1) {
    switch(c){
    case 'd':
      fDataName=optarg;
      break; 
   case 'm':
      fMCName=optarg;
      break;
   case 'o':
      OutputName=optarg;
      break;
    case 'n':
      NIterations=atoi(optarg);
      break;
    case 't':
      NToys=atoi(optarg);
      break;
    case 's':
      StatisticalFluctuations=true;
      break;
    case 'y':
      SystematicFluctuations=true;
      break;
    case 'i':
      NIterationFluctuations=true;
      break;
    case 'a':
      NTarget=atoi(optarg);
      break;
    }
  }

  //////////////////////////////////PREPARE FOR READING OF TXT AND ROOT FILES/////////////////////////////////////
  char * rootDataName = new char[256];
  char * rootMCName = new char[256];
  sprintf(rootMCName,"%s.root",fMCName);
  sprintf(rootDataName,"%s.root",fDataName);
  
  char * txtDataName = new char[256];
  char * txtMCName = new char[256]; 
  sprintf(txtMCName,"%s.txt",fMCName);
  sprintf(txtDataName,"%s.txt",fDataName);

  char * txtOutputName = new char[256];
  char * rootOutputName = new char[256];
  sprintf(txtOutputName,"%s.txt",OutputName);
  sprintf(rootOutputName,"%s.root",OutputName);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"looking for files named: "<<txtDataName<<" and "<<txtMCName<<endl;

  double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double DataReconstructedEvents[NBinsRecMom][NBinsRecAngle], MCReconstructedBkgEvents[NBinsRecMom][NBinsRecAngle], MCReconstructedEvents[NBinsRecMom][NBinsRecAngle];
  double MCEfficiency[NBinsTrueMom][NBinsTrueAngle];
  double TrueEventsDistribution[NBinsTrueMom][NBinsTrueAngle];
  double * NumberOfNeutrino=new double();
  
  double UnfoldedData[NToys][NBinsTrueMom][NBinsTrueAngle];
  double table_TrueEventsDistribution[NToys][NBinsTrueMom][NBinsTrueAngle];

  for(int nt=0;nt<NToys;nt++){//loop on toy experiments
    
    //0. Load the input distributions
    XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,MCEfficiency,NumberOfNeutrino);
    
#ifdef DEBUG
    cout<<"DEBUG, LOAD FILES/////////////////////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"DEBUGGING, check number of neutrino chosen (read in the -d file):"<<*NumberOfNeutrino<<endl;
        
    cout<<"DEBUGGING, check data file events:"<<endl;
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	cout<<DataReconstructedEvents[e0][e1]<<" ";
      }
      cout<<endl;
    }

    cout<<endl<<"DEBUGGING, check MC file events:"<<endl;
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	cout<<MCReconstructedEvents[e0][e1]<<" ";
      }
      cout<<endl;
    }
    cout<<endl<<"DEBUGGING, check MC background events:"<<endl;
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	cout<<MCReconstructedBkgEvents[e0][e1]<<" ";
      }
      cout<<endl;
    }
    
    cout<<endl<<"DEBUGGING, True distribution of events (pmu, thetamu, ironmu,recthetamu):"<<endl;
   for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    cout<<MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]<<" ";
	  }
	}
      }
    }
    cout<<endl;
    
    double TrueRecBins[NBinsRecMom][NBinsRecAngle];
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    TrueRecBins[e0][e1]=0;
	    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
		TrueRecBins[e0][e1]+=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	    }
	}
      }
    }

    cout<<endl<<"test data - background in rec bins vs sum of true distribution (in case of MC w/o error, should be the same"<<endl;
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	cout<<DataReconstructedEvents[e0][e1]-MCReconstructedBkgEvents[e0][e1]<<"vs"<<TrueRecBins[e0][e1]<<" ";
      }
      cout<<endl;
    }
    
#endif

    if(StatisticalFluctuations){
      //Toy0. Statistically vary the DataReconstructedEvents distribution.
      XS->Xsec::GenerateStatisticalFluctuations(DataReconstructedEvents);
    }

    if(SystematicFluctuations){
      //Toy1. Fluctuation of the MCpart: prior, bkg, unfolding matrix (through likelihood) -> Looks systematics detector variation
      //Step 1: the user should specify variation he wants:
      double RelativeSigma[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	      RelativeSigma[c0][c1][e0][e1]=0.1;
	    }
	  }
	}
      }
    //Step 2: deduce the variation
      XS->Xsec::GenerateMCFluctuations(MCReconstructedEvents_TrueSignal,RelativeSigma);
    }

    if(NIterationFluctuations){
      //Toy 3: vary the number of iterations.
      NIterations=NToys;
    }

    //1. Build the likelihood matrix & intial prior from the MC
    XS->Xsec::BuildLikelihood(vLikelihood, vInitialPriorMC, TrueEventsDistribution, MCReconstructedEvents_TrueSignal);
    
#ifdef DEBUG    
    cout<<endl<<"Likelihood matrix (pmu, thetamu, ironmu,recthetamu):"<<endl;
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    cout<<vLikelihood[c0][c1][e0][e1]<<" ";
	  }
	}
      }
    }
    cout<<endl;
    
    cout<<endl<<"Initial prior (pmu, thetamu):"<<endl;
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    cout<<vInitialPriorMC[c0][c1]<<" ";
      }
    }
    cout<<endl;
    
#endif

    for(int it=0;it<NIterations;it++){      
      //2. Build the prior used in the unfolding
      XS->Xsec::SetPrior(vPriorNormalised,vPrior,vInitialPriorMC,vInitialPrior,vPosterior,PriorMC,it);

      
      //3. Build the unfolding matrix
      XS->Xsec::BuildUnfolding(vUnfolding,vLikelihood,vPriorNormalised);
      

      //4. Apply the unfolding on data
      XS->Xsec::ApplyUnfolding(vPosterior, vUnfolding, DataReconstructedEvents, MCReconstructedBkgEvents);

#ifdef DEBUG
      cout<<"Iteration Step #"<<it+1<<" ////////////////////////////////////////////////////////////////////////"<<endl;

      cout<<endl<<"Prior chosen (pmu, thetamu):"<<endl;
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  cout<<vPriorNormalised[c0][c1]<<" ";
	}
      }
      cout<<endl;

      cout<<endl<<"Unfolding matrix (pmu, thetamu, ironmu,recthetamu):"<<endl;
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    cout<<vUnfolding[c0][c1][e0][e1]<<" ";
	  }
	}
      }
    }
    cout<<endl;
    
    cout<<endl<<"Posterior (pmu, thetamu):"<<endl;
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	cout<<vPosterior[c0][c1]<<" ";
      }
    }
    cout<<endl;
    
#endif

    }

    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	UnfoldedData[nt][c0][c1]=vPosterior[c0][c1];
	table_TrueEventsDistribution[nt][c0][c1]=TrueEventsDistribution[c0][c1];
      }
    }
  }

#ifdef DEBUG
  TMatrixD LMom(NBinsTrueMom,NBinsRecMom);
  TMatrixD UMom(NBinsTrueMom,NBinsRecMom);
  TH1D * PriorMom = new TH1D("PriorMom","Prior of the 1D momentum",NBinsTrueMom,BinningTrueMom);
  for(int c0=0;c0<NBinsTrueMom;c0++){
    for(int e0=0;e0<NBinsRecMom;e0++){
      LMom[c0][e0]=0;
      UMom[c0][e0]=0;
    }
  }
  for(int c0=0;c0<NBinsTrueMom;c0++){
    for(int c1=0;c1<NBinsTrueAngle;c1++){
      PriorMom->Fill(PriorMom->GetBinCenter(c0+1),vPriorNormalised[c0][c1]);
      for(int e0=0;e0<NBinsRecMom;e0++){
	for(int e1=0;e1<NBinsRecAngle;e1++){
	  LMom[c0][e0]+=vLikelihood[c0][c1][e0][e1];
	  UMom[c0][e0]+=vUnfolding[c0][c0][e0][e1];
	}
      }
    }
  }

  TCanvas * cLMom = new TCanvas("cLMom","Likelihood matrix for 1D momentumXdistance");
  LMom.Draw("colz");
  TCanvas * cUMom = new TCanvas("cUMom","Unfolding matrix for 1D momentumXdistance");
  UMom.Draw("colz");
  TCanvas * cPriorMom = new TCanvas("cPriorMom","Prior for 1D momentumXdistance");
  PriorMom->Draw();
  
#endif
  
  TCanvas * canEvents[NBinsTrueMom][NBinsTrueAngle];
  TH1D * Events[NBinsTrueMom][NBinsTrueAngle];
  TCanvas * canxsection[NBinsTrueMom][NBinsTrueAngle];
  TH1D * hxsection[NBinsTrueMom][NBinsTrueAngle];
  char Name[256];char Title[256];
  double xsection[NBinsTrueMom][NBinsTrueAngle];

  TFile * ROOT_file = new TFile(rootMCName);
  TH2D * hTrueDistribution = (TH2D*) ROOT_file->Get("MCTrueEvents");
  TH2D * htruedistribution = (TH2D*) hTrueDistribution->Clone("htruedistribution");
  double truedistribution[NBinsTrueMom][NBinsTrueAngle];

  TFile * file = new TFile(rootOutputName,"recreate");
  file->cd();

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
#ifdef DEBUG
	  cout<<"Bin pmu="<<c0<<", Angle="<<c1<<endl;
#endif
	sprintf(Name,"Eventsbin%d_%d",c0,c1);
	sprintf(Title,"True distribution after unfolding in the bin %d_%d",c0,c1);
	Events[c0][c1] = new TH1D(Name,Title,40,-0.1,0.1);

	sprintf(Name,"xsectionbin%d_%d",c0,c1);
	sprintf(Title,"Double differential xsection after unfolding in the bin %d_%d",c0,c1);
	hxsection[c0][c1] = new TH1D(Name,Title,40,-0.1,0.1);

	
	for(int nt=0;nt<NToys;nt++){
	  double Value=UnfoldedData[nt][c0][c1]-table_TrueEventsDistribution[0][c0][c1];
	  if(table_TrueEventsDistribution[0][c0][c1]!=0) Value/=table_TrueEventsDistribution[0][c0][c1];
	  Events[c0][c1]->Fill(Value);

	  double BinningWidth=(BinningTrueMom[c0+1]-BinningTrueMom[c0])*(BinningTrueAngle[c1+1]-BinningTrueAngle[c1]);
	  double Flux=*NumberOfNeutrino;
	  double Efficiency=MCEfficiency[c0][c1];
	  double Correction=Flux*BinningWidth*Efficiency*NTarget;
	  //if(Correction!=0) Value/=Correction;
	  xsection[c0][c1]=UnfoldedData[nt][c0][c1];
	  truedistribution[c0][c1]=hTrueDistribution->GetBinContent(c0+1,c1+1);
								   
	  hxsection[c0][c1]->Fill(xsection[c0][c1]);
	  htruedistribution->SetBinContent(c0+1,c1+1,truedistribution[c0][c1]);
	  
#ifdef DEBUG
	  cout<<"Number of events of the Toy="<<UnfoldedData[nt][c0][c1]<<", true="<<table_TrueEventsDistribution[0][c0][c1]<<", value="<<Value<<endl;
#endif

	}
	  
	sprintf(Name,"canEvents%d_%d",c0,c1);
	sprintf(Title,"True distribution of events after unfolding in the bin %d_%d",c0,c1);
	canEvents[c0][c1] = new TCanvas(Name,Title);
	Events[c0][c1]->GetXaxis()->SetTitle("Pull");
	Events[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	Events[c0][c1]->Draw();
	Events[c0][c1]->Fit("gaus","R");

	sprintf(Name,"canxsection%d_%d",c0,c1);
	sprintf(Title,"Double differential xsection after unfolding in the bin %d_%d",c0,c1);
	canxsection[c0][c1] = new TCanvas(Name,Title);
	hxsection[c0][c1]->GetXaxis()->SetTitle("Pull");
	hxsection[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	hxsection[c0][c1]->Draw();
	hxsection[c0][c1]->Fit("gaus","R");
	
	canEvents[c0][c1]->Write();
	canxsection[c0][c1]->Write();
      }
  }
  if(NToys==1){
    ofstream fEvent;
    fEvent.open(txtOutputName,ios::out);
   
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	fEvent<<xsection[c0][c1]<<" ";
      }
    }	  
    fEvent<<666;
    fEvent.close();
    
     htruedistribution->Write();
  }
  file->Close();
  //theApp.Run();

  return 0;
}

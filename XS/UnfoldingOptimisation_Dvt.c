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
//#define DEBUG
//#define DEBUG2
//#define FLUXCONTROL

int main(int argc, char ** argv){

  //TApplication theApp("App",0,0);
  gStyle->SetOptFit(kTRUE);

  bool PriorMC=true;
  int c=-1;
  int NPriors=0;
  //Signal
  char * fDataName = new char[256];
  char * fMCName = new char[256];
  sprintf(fDataName,"DistributionsData.root");
  sprintf(fMCName,"DistributionsMC.root");
  //Will be only used if side band
  char * fDataNameSB = new char[256];
  char * fMCNameSB = new char[256];
  sprintf(fDataNameSB,"DistributionsDataSB.root");
  sprintf(fMCNameSB,"DistributionsMCSB.root");
  //
  char * OutputName = new char[256];
  const int NVariations=3;
  int NToys[NVariations]={1,1,1};//Number of toy experiement varying respectively number of iterations (0), number of systematics variation (1) and number of stat. variations (2). 
  bool StatisticalFluctuations=false;
  bool SystematicFluctuations=false;
  bool NIterationFluctuations=false;
  //int NTarget=8.299e28;
  double NTarget = NTargetPM;
  bool SideBand=false; bool SideBandData=false; bool SideBandMC=false;
  bool FakeData=false;
  bool BkgSub=false;
  bool StatisticalUncertaintyEstimation = false;
  
  while ((c = getopt(argc, argv, "d:m:n:a:s:y:io:b:c:pfBSw")) != -1) {
    switch(c){
    case 'd':
      fDataName=optarg;
      break; 
   case 'm':
      fMCName=optarg;
      break;
    case 'n':
      NToys[0]=atoi(optarg);
      break;
   case 'o':
      OutputName=optarg;
      break;
    case 's':
      StatisticalFluctuations=true;
      NToys[2]=atoi(optarg);
      break;
    case 'y':
      SystematicFluctuations=true;
      NToys[1]=atoi(optarg);
      break;
    case 'i':
      NIterationFluctuations=true;
      //NToys[0]=atoi(optarg);
      break;
    case 'a':
      NTarget=atof(optarg);
      break;
    case 'b':
      SideBandData=true;
      fDataNameSB=optarg;
      break;
    case 'c':
      SideBandMC=true;
      fMCNameSB=optarg;
      break;
    case 'p':
      PriorMC=false;
      cout<<"For now, a flat prior is chosen if the PriorMC is not chosen as in this case!"<<endl;
      break;
    case 'f':
      FakeData=true;
      cout<<"You are using fake data set, not real data"<<endl;
      break;
    case 'B':
      BkgSub=true;
      cout<<"You are using a background substraction method"<<endl;
      break;
    case 'S':
      StatisticalUncertaintyEstimation=true;
      cout<<"You are really estimating the stat error."<<endl;
      break;
    case 'w':
      NTarget = NTargetWM;
      break;
    }
  }
  if(SideBandData && SideBandMC){
    SideBand=true;
    cout<<"We will use a side band"<<endl;
  }

  //InitializeGlobal();
  InitializeGlobal(true,1);
#ifdef DEBUG
  cout<<"Initialized"<<endl;
  cout<<"Nbins true mom="<<NBinsTrueMom<<", true angle="<<NBinsTrueAngle<<", rec mom="<<NBinsRecMom<<", rec angle="<<NBinsRecAngle<<endl;
  cout<<"Signal: Nbins true mom="<<NBinsTrueMomSignal<<", true angle="<<NBinsTrueAngleSignal<<", rec mom="<<NBinsRecMomSignal<<", rec angle="<<NBinsRecAngleSignal<<endl;
  cout<<"SB: Nbins true mom="<<NBinsTrueMomSB<<", true angle="<<NBinsTrueAngleSB<<", rec mom="<<NBinsRecMomSB<<", rec angle="<<NBinsRecAngleSB<<endl;
#endif
  Xsec * XS = new Xsec();
  //To modify the bining if we have a side band
  ReinitializeUnfoldingBinning(SideBand);
#ifdef DEBUG
  cout<<"ReInitialized"<<endl;
  cout<<"Nbins true mom="<<NBinsTrueMom<<", true angle="<<NBinsTrueAngle<<", rec mom="<<NBinsRecMom<<", rec angle="<<NBinsRecAngle<<endl;
  cout<<"Signal: Nbins true mom="<<NBinsTrueMomSignal<<", true angle="<<NBinsTrueAngleSignal<<", rec mom="<<NBinsRecMomSignal<<", rec angle="<<NBinsRecAngleSignal<<endl;
  cout<<"SB: Nbins true mom="<<NBinsTrueMomSB<<", true angle="<<NBinsTrueAngleSB<<", rec mom="<<NBinsRecMomSB<<", rec angle="<<NBinsRecAngleSB<<endl;
#endif
  
  /////////////////////////////////////////////
  double **** vLikelihood = new double ***[NBinsTrueMom];
  double **** vUnfolding = new double ***[NBinsTrueMom];
  double **** MCReconstructedEvents_TrueSignal = new double ***[NBinsTrueMom];
  double **** MCReconstructedEvents_TrueSignal_Default = new double ***[NBinsTrueMom];
  double **** DataReconstructedEvents_TrueSignal = new double ***[NBinsTrueMom];
  double **** DataReconstructedEvents_TrueSignal_Default = new double ***[NBinsTrueMom];
  double **** RelativeSigma = new double ***[NBinsTrueMom];

  for(int h=0;h<NBinsTrueMom;h++){
    vLikelihood[h] = new double **[NBinsTrueAngle];
    vUnfolding[h] = new double **[NBinsTrueAngle];
    MCReconstructedEvents_TrueSignal[h] = new double **[NBinsTrueAngle];
    MCReconstructedEvents_TrueSignal_Default[h] = new double **[NBinsTrueAngle];
    DataReconstructedEvents_TrueSignal[h] = new double **[NBinsTrueAngle];
    DataReconstructedEvents_TrueSignal_Default[h] = new double **[NBinsTrueAngle];
    RelativeSigma[h] = new double **[NBinsTrueAngle];
    for(int i=0;i<NBinsTrueAngle;i++){
      vLikelihood[h][i] = new double *[NBinsRecMom];
      vUnfolding[h][i] = new double *[NBinsRecMom];
      MCReconstructedEvents_TrueSignal[h][i] = new double *[NBinsRecMom];
      MCReconstructedEvents_TrueSignal_Default[h][i] = new double *[NBinsRecMom];
      DataReconstructedEvents_TrueSignal[h][i] = new double *[NBinsRecMom];
      DataReconstructedEvents_TrueSignal_Default[h][i] = new double *[NBinsRecMom];
      RelativeSigma[h][i] = new double *[NBinsRecMom];
      for(int j=0;j<NBinsRecMom;j++){
	vLikelihood[h][i][j] = new double [NBinsRecAngle];
	vUnfolding[h][i][j] = new double [NBinsRecAngle];
	MCReconstructedEvents_TrueSignal[h][i][j] = new double [NBinsRecAngle];
	MCReconstructedEvents_TrueSignal_Default[h][i][j] = new double [NBinsRecAngle];
	DataReconstructedEvents_TrueSignal[h][i][j] = new double [NBinsRecAngle];
	DataReconstructedEvents_TrueSignal_Default[h][i][j] = new double [NBinsRecAngle];
	RelativeSigma[h][i][j] = new double [NBinsRecAngle];
      }      
    }
  }

  double ** vInitialPriorMC = new double*[NBinsTrueMom];
  double ** vInitialPrior = new double*[NBinsTrueMom];
  double ** vPrior = new double*[NBinsTrueMom];
  double ** vPriorNormalised = new double*[NBinsTrueMom];
  double ** vPosterior = new double*[NBinsTrueMom];
  double ** vPosterior_SignalOnly = new double*[NBinsTrueMom];
  double ** DataEfficiency = new double*[NBinsTrueMom];
  double ** MCEfficiency_Default = new double*[NBinsTrueMom];
  double ** MCEfficiency = new double*[NBinsTrueMom];
  double ** UnfoldedData_SignalOnly = new double*[NBinsTrueMom];
  double ** TrueEventsDistribution_SignalOnly = new double*[NBinsTrueMom];
  double ** table_TrueEventsDistribution_SignalOnly = new double*[NBinsTrueMom];
  double ** XSection_SignalOnly = new double*[NBinsTrueMom];
  double ** TrueXSection_SignalOnly = new double*[NBinsTrueMom];//not used
  for(int i=0;i<NBinsTrueMom;i++){
    vInitialPriorMC[i] = new double[NBinsTrueAngle];
    vInitialPrior[i] = new double[NBinsTrueAngle];
    vPrior[i] = new double[NBinsTrueAngle];
    vPriorNormalised[i] = new double[NBinsTrueAngle];
    vPosterior[i] = new double[NBinsTrueAngle];
    vPosterior_SignalOnly[i] = new double[NBinsTrueAngle];
    DataEfficiency[i] = new double [NBinsTrueAngle];
    MCEfficiency_Default[i] = new double [NBinsTrueAngle];
    MCEfficiency[i] = new double [NBinsTrueAngle];
    UnfoldedData_SignalOnly[i] = new double [NBinsTrueAngle];
    TrueEventsDistribution_SignalOnly[i] = new double [NBinsTrueAngle];
    table_TrueEventsDistribution_SignalOnly[i] = new double [NBinsTrueAngle];
    XSection_SignalOnly[i] = new double [NBinsTrueAngle];
    TrueXSection_SignalOnly[i] = new double [NBinsTrueAngle];//not used
  }

  double ** DataReconstructedEvents = new double*[NBinsRecMom];
  double ** DataReconstructedEvents_Default = new double*[NBinsRecMom];
  double ** MCReconstructedBkgEvents = new double*[NBinsRecMom];
  double ** MCReconstructedEvents = new double*[NBinsRecMom];
  double ** TrueRecBins = new double*[NBinsRecMom];
  for(int i=0;i<NBinsRecMom;i++){
    DataReconstructedEvents[i] = new double [NBinsRecAngle];
    DataReconstructedEvents_Default[i] = new double [NBinsRecAngle];
    MCReconstructedBkgEvents[i] = new double [NBinsRecAngle];
    MCReconstructedEvents[i] = new double [NBinsRecAngle];
    TrueRecBins[i] = new double [NBinsRecAngle];
  }
  double * NumberOfNeutrino = new double;

  TMatrixD * MUnfolding = new TMatrixD(NBinsTrueMom*NBinsTrueAngle,NBinsRecMom*NBinsRecAngle);
  TMatrixD * MLikelihood = new TMatrixD(NBinsTrueMom*NBinsTrueAngle,NBinsRecMom*NBinsRecAngle);
  TH2D * hPrior[NToys[0]];
  for(int it=0;it<NToys[0];it++){
    hPrior[it]= new TH2D(Form("hPrior%d",it),"",NBinsTrueMom,BinningTrueMom,NBinsTrueAngle,BinningTrueAngle);
  }
  TVectorD * VReconstructedEvents = new TVectorD(NBinsRecMom*NBinsRecAngle);

  //cout<<NToys[0]<<", "<<NToys[1]<<", "<<NToys[2]<<endl;
  //////////////////////////////////PREPARE FOR READING OF TXT AND ROOT FILES/////////////////////////////////////
  char * rootDataName = new char[256];
  char * rootMCName = new char[256];
  sprintf(rootMCName,"%s.root",fMCName);
  sprintf(rootDataName,"%s.root",fDataName);
  
  char * txtDataName = new char[256];
  char * txtMCName = new char[256]; 
  sprintf(txtMCName,"%s.txt",fMCName);
  sprintf(txtDataName,"%s.txt",fDataName);

  cout<<"looking for files named: "<<txtDataName<<" and "<<txtMCName<<endl;

  //Side band. Used only if a side band option is activated
  char * rootDataNameSB = new char[256];
  char * rootMCNameSB = new char[256];
  sprintf(rootMCNameSB,"%s.root",fMCNameSB);
  sprintf(rootDataNameSB,"%s.root",fDataNameSB);
  
  char * txtDataNameSB = new char[256];
  char * txtMCNameSB = new char[256]; 
  sprintf(txtMCNameSB,"%s.txt",fMCNameSB);
  sprintf(txtDataNameSB,"%s.txt",fDataNameSB);
  
  cout<<"Side band: looking for files named: "<<txtDataNameSB<<" and "<<txtMCNameSB<<endl;
  
  ///////////////////////////////////OUTPUT TREE/////////////
  char * txtOutputName = new char[256];
  char * rootOutputName = new char[256];
  sprintf(txtOutputName,"%s.txt",OutputName);
  sprintf(rootOutputName,"%s.root",OutputName);
  int nt0;int nt1;int nt2;
  TFile * file = new TFile(rootOutputName,"recreate");
  file->cd();

  /*
  double UnfoldedData_SignalOnly_Fill[NBinsTrueMomMax][NBinsTrueAngleMax];
  double XSection_SignalOnly_Fill[NBinsTrueMomMax][NBinsTrueAngleMax];
  double table_TrueEventsDistribution_SignalOnly_Fill[NBinsTrueMomMax][NBinsTrueAngleMax];
  double MCReconstructedEvents_TrueSignal_Fill[NBinsTrueMomMax][NBinsTrueAngleMax][NBinsRecMomMax][NBinsRecAngleMax];
  */
  double UnfoldedData_SignalOnly_Fill[NBinsTrueMom][NBinsTrueAngle];
  double XSection_SignalOnly_Fill[NBinsTrueMom][NBinsTrueAngle];
  double TrueXSection_SignalOnly_Fill[NBinsTrueMom][NBinsTrueAngle];
  double table_TrueEventsDistribution_SignalOnly_Fill[NBinsTrueMom][NBinsTrueAngle];
  double DataReconstructedEvents_TrueSignal_Fill[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double DataReconstructedEvents_Fill[NBinsRecMom][NBinsRecAngle];
  
  TTree*              wtree    = new TTree("wtree","wtree");
  wtree->SetDirectory(file);
  wtree              -> Branch   ("Iterations",&nt0,"Iterations/I");
  wtree              -> Branch   ("Systematics",&nt1,"Systematics/I");
  wtree              -> Branch   ("Statistics",&nt2,"Statistics/I");
  
  /*  
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 1
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      wtree              -> Branch   (Form("Events_%d_%d",c0,c1),&(UnfoldedData_SignalOnly[c0][c1]),Form("Events_%d_%d/D",c0,c1));
      wtree              -> Branch   (Form("TrueEvents_%d_%d",c0,c1),&(table_TrueEventsDistribution_SignalOnly[c0][c1]),Form("TrueEvents_%d_%d/D",c0,c1));
      wtree              -> Branch   (Form("XSection_SignalOnly_%d_%d",c0,c1),&(XSection_SignalOnly[c0][c1]),Form("XSection_SignalOnly_%d_%d/D",c0,c1));
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over cause 1
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over cause 1
	  wtree              -> Branch   (Form("EventsAll_%d_%d_%d_%d",c0,c1,e0,e1),&(MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]),Form("EventsAll_%d_%d_%d_%d/D",c0,c1,e0,e1));
	}
      }
    }
  }
  */  
  //wtree              -> Branch   ("Events",&UnfoldedData_SignalOnly);
  //wtree              -> Branch   ("Events",&UnfoldedData_SignalOnly,"Events/D");
  wtree              -> Branch   ("Events",UnfoldedData_SignalOnly_Fill,Form("Events[%d][%d]/D",NBinsTrueMom,NBinsTrueAngle));
  wtree              -> Branch   ("EventsAll",DataReconstructedEvents_TrueSignal_Fill,Form("EventsAll[%d][%d][%d][%d]/D",NBinsTrueMom,NBinsTrueAngle,NBinsRecMom,NBinsRecAngle));
  wtree              -> Branch   ("TrueEvents",table_TrueEventsDistribution_SignalOnly_Fill,Form("TrueEvents[%d][%d]/D",NBinsTrueMom,NBinsTrueAngle));
  wtree              -> Branch   ("XSection",XSection_SignalOnly_Fill,Form("XSection[%d][%d]/D",NBinsTrueMom,NBinsTrueAngle));
  wtree              -> Branch   ("TrueXSection",TrueXSection_SignalOnly_Fill,Form("TrueXSection[%d][%d]/D",NBinsTrueMom,NBinsTrueAngle));
  wtree              -> Branch   ("EventsRec",DataReconstructedEvents_Fill,Form("EventsRec[%d][%d]/D",NBinsRecMom,NBinsRecAngle));

#ifdef FLUXCONTROL
  TFile * fFlux = (TFile*) new TFile("FluxTest.root","update");
  //TH1D * hFlux = new TH1D("hFlux","",1e4,1e11,1e15);
  TH1D * hFlux = (TH1D*) fFlux->Get("hFlux");
#endif

  //00. IF FLAT PRIOR
  double Norm = 0;
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      vInitialPrior[c0][c1]=1;
      Norm += vInitialPrior[c0][c1];
      ///(NBinsTrueMom*NBinsTrueAngle);
    }
  }
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      vInitialPrior[c0][c1] /= Norm;
    }
  }

  //0. Load the input distributions
  if(SideBand) XS->Xsec::LoadInputFilesSB(txtDataName,txtMCName,txtDataNameSB,txtMCNameSB,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,MCEfficiency,NumberOfNeutrino,DataReconstructedEvents_TrueSignal,DataEfficiency,FakeData);
  else XS->Xsec::LoadInputFiles(txtDataName,txtMCName,MCReconstructedEvents_TrueSignal,DataReconstructedEvents,MCReconstructedEvents,MCReconstructedBkgEvents,MCEfficiency,NumberOfNeutrino,DataReconstructedEvents_TrueSignal,DataEfficiency,FakeData);

#ifdef DEBUG
    cout<<"DEBUG, LOAD FILES/////////////////////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"DEBUGGING, check number of neutrino chosen (read in the -d file):"<<*NumberOfNeutrino<<endl;
#endif
#ifdef DEBUG2
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
      cout<<"Mom bin:"<<c0<<endl;
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	cout<<"Angle bin:"<<c1<<endl;
	double NeV=0;
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    NeV += DataReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	    //cout<<MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]<<" ";
	  }
	}
	cout<<NeV<<endl;
      }
    }
    cout<<endl;
    
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	TrueRecBins[e0][e1]=0;
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    TrueRecBins[e0][e1]+=DataReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	  }
	}
      }
    }
    
    cout<<endl<<"test data - background in rec bins vs sum of true distribution (in case of MC w/o error, should be the same"<<endl;
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	//cout<<DataReconstructedEvents[e0][e1]-MCReconstructedBkgEvents[e0][e1]<<"vs"<<TrueRecBins[e0][e1]<<" ";
	cout<<DataReconstructedEvents[e0][e1]<<"vs"<<TrueRecBins[e0][e1]<<" ";
      }
      cout<<endl;
    }
    
#endif
    for(nt1=0;nt1<NToys[1];nt1++){//loop on MC variations (syst)
      cout<<nt1<<endl;

#ifdef DEBUG
      cout<<"######################################################################"<<endl;
      cout<<"#################### Systematic toy #"<<nt1<<" ###############################"<<endl;
      cout<<"######################################################################"<<endl;
#endif
      
      if(nt1==0){//for each change of MC variations, we reset the Default to the original one.
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    DataReconstructedEvents_Default[e0][e1]=DataReconstructedEvents[e0][e1];
	    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
		MCReconstructedEvents_TrueSignal_Default[c0][c1][e0][e1]=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
		DataReconstructedEvents_TrueSignal_Default[c0][c1][e0][e1]=DataReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	      }
	    }
	  }
	}
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    MCEfficiency_Default[c0][c1] = MCEfficiency[c0][c1];
	  }
	}
      }
      else{//after each syst variation, we reset the values to default. Then, we apply variations.
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
		MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]=MCReconstructedEvents_TrueSignal_Default[c0][c1][e0][e1];
		/*if(!FakeData)*/ DataReconstructedEvents_TrueSignal[c0][c1][e0][e1]=DataReconstructedEvents_TrueSignal_Default[c0][c1][e0][e1];
	      }
	    }
	  }
	}
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    MCEfficiency[c0][c1] = MCEfficiency_Default[c0][c1];
	  }
	}	  

	//Toy1. Fluctuation of the MCpart: prior, bkg, unfolding matrix (through likelihood) -> Looks systematics detector variation
	//Step 1: the user should specify variation he wants:
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
		RelativeSigma[c0][c1][e0][e1]=.2;
	      }
	    }
	  }
	}
      }


	
      for(nt2=0;nt2<NToys[2];nt2++){//loop on data variations (stat)
	if(nt2%10 == 0) cout<<"Number of statistical toys = "<<nt2<<endl;
#ifdef DEBUG
      cout<<"#################### Statistics toy #"<<nt2<<" ###############################"<<endl;
#endif
      
	  
      //after each stat variation, we reset the values to default. Then, we apply variations.
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  /*if(!FakeData)*/ DataReconstructedEvents[e0][e1]=DataReconstructedEvents_Default[e0][e1];
	  //MCReconstructedEvents[e0][e1]=MCReconstructedEvents_Default[e0][e1];
	}
      }
      
	if(SystematicFluctuations && nt2==0){
	  //Step 2: deduce the variation
	  /*if(FakeData){
	    XS->Xsec::GenerateDataFluctuations(DataReconstructedEvents,DataReconstructedEvents_TrueSignal,RelativeSigma);
	    cout<<"variation applied"<<endl;
	    }
	    else*/
	  //XS->Xsec::GenerateDetectorMCFluctuations(MCReconstructedEvents_TrueSignal,RelativeSigma);
	  XS->Xsec::GenerateXSModelMCFluctuations(MCEfficiency,MCReconstructedEvents_TrueSignal,RelativeSigma);
	}
	if(StatisticalFluctuations && nt2!=0){
	  //Toy0. Statistically vary the DataReconstructedEvents distribution.
	  XS->Xsec::GenerateStatisticalFluctuations(DataReconstructedEvents);
	}
	
	//1. Build the likelihood matrix & intial prior from the MC
	//TO UNCOMMENT INSTEAD OF BELOW
	XS->Xsec::BuildLikelihood(vLikelihood, vInitialPriorMC, MCReconstructedEvents_TrueSignal);
	//TEMP
	//XS->Xsec::BuildLikelihood(vLikelihood, vInitialPriorMC, DataReconstructedEvents_TrueSignal);
	//Provide the true signal distribution projected only on true bins. The important point is that we only sum over the reconstructed bins of the SIGNAL REGION
	XS->Xsec::ProjectOnTruePhaseSpace_OnlySignal(TrueEventsDistribution_SignalOnly, DataReconstructedEvents_TrueSignal);

	//TEMP
	/*	
	double TotalMCContent=0;
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    double TrueMCContent=0;
	    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1 
		TrueMCContent += DataReconstructedEvents_TrueSignal[c0][c1][e0][e1];
		TotalMCContent += DataReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	      }
	    }
	    vInitialPriorMC[c0][c1] = TrueMCContent;
	  }
	}
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    if(TotalMCContent != 0) vInitialPriorMC[c0][c1] /= TotalMCContent;
	  }
	  }*/
#ifdef DEBUG2
	
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
	if(PriorMC) cout<<vInitialPriorMC[c0][c1]<<" ";
	else cout<<vInitialPrior[c0][c1]<<" ";
      }
    }
    cout<<endl;
    
#endif

    for(nt0=0;nt0<NToys[0];nt0++){
#ifdef DEBUG
      if(nt0%10 == 0) cout<<"Number of toys = "<<nt0<<endl;
#endif
      
      //2. Build the prior used in the unfolding
      XS->Xsec::SetPrior(vPriorNormalised,vPrior,vInitialPriorMC,vInitialPrior,vPosterior,PriorMC,nt0);
      if(nt1==0 && nt2==0){
	for(int ibinx=0;ibinx<NBinsTrueMom;ibinx++){
	  for(int ibiny=0;ibiny<NBinsTrueAngle;ibiny++){
	    hPrior[nt0]->SetBinContent(ibinx+1,ibiny+1,vPriorNormalised[ibinx][ibiny]);
	  }
	}
      }
      
#ifdef DEBUG
      if(nt0%10 == 0) cout<<"Prior is set"<<endl;
#endif      
      //3. Build the unfolding matrix
      XS->Xsec::BuildUnfolding(vUnfolding,vLikelihood,vPriorNormalised);
      
#ifdef DEBUG
      if(nt0%10 == 0) cout<<"Unfolding matrix is built"<<endl;
#endif
      //4. Apply the unfolding on data
      if(BkgSub) XS->Xsec::ApplyUnfoldingBkgSubstraction(vPosterior, vUnfolding, DataReconstructedEvents, MCReconstructedBkgEvents);
      else XS->Xsec::ApplyUnfolding(vPosterior, vPosterior_SignalOnly, vUnfolding, DataReconstructedEvents);
    
#ifdef DEBUG2
      cout<<"Iteration Step #"<<nt0+1<<" ////////////////////////////////////////////////////////////////////////"<<endl;

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
    
#ifdef FLUXCONTROL
    hFlux->Fill(*NumberOfNeutrino);
#endif
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	//Data unfolded that we wish to give to use are not the whole unfolding result. Only the sum over reconstructed bins in the signal reegion. In the example of CC0pi cross section, we do not want the true CC0pi in the reconstructed CC1pi sample to be added to our signal bin result.
	UnfoldedData_SignalOnly[c0][c1]=vPosterior_SignalOnly[c0][c1];
	
	table_TrueEventsDistribution_SignalOnly[c0][c1]=TrueEventsDistribution_SignalOnly[c0][c1];
	//trick for side band, because lower limit of the first bin pf side band is always assumed to be 0.
	double BinLowerLimitMom = BinningTrueMom[c0];
	double BinLowerLimitAngle = BinningTrueAngle[c1];	
	//if(c0==NBinsTrueMomSignal) BinLowerLimitMom = 0;
	//if(c1==NBinsTrueAngleSignal) BinLowerLimitAngle = 0;
	//
	double BinningWidth=(BinningTrueMom[c0+1]-BinLowerLimitMom)*(BinningTrueAngle[c1+1]-BinLowerLimitAngle);
      
	double Flux=*NumberOfNeutrino;
	double Efficiency=MCEfficiency[c0][c1];
	double Correction=Flux*BinningWidth*Efficiency*NTarget;
#ifdef DEBUG
	cout<<c0<<", "<<c1<<", correction = "<<Correction<<", Bin width = "<<BinningWidth<<", flux = "<<Flux<<", NTarget = "<<NTarget<<", Efficiency = "<<Efficiency<<endl;
#endif
	XSection_SignalOnly[c0][c1]=UnfoldedData_SignalOnly[c0][c1]/ ( Correction == 0 ? 1 : Correction );
	if(FakeData){//In this case, the true XS is given by correcting with the true data efficiency
	  double dataEfficiency=DataEfficiency[c0][c1];
	  Correction=Flux*BinningWidth*dataEfficiency*NTarget;
	  TrueXSection_SignalOnly[c0][c1]=TrueEventsDistribution_SignalOnly[c0][c1]/ ( Correction == 0 ? 1 : Correction );
	}
	else{
	  TrueXSection_SignalOnly[c0][c1]=TrueEventsDistribution_SignalOnly[c0][c1]/ ( Correction == 0 ? 1 : Correction );
	}
	//cout<<"bin ("<<c0<<","<<c1<<"), benji says efficiency="<<Efficiency<<", correction="<<Correction<<", Nev="<<XSection_SignalOnly[c0][c1]<<", XS="<<XSection_SignalOnly[c0][c1]<<endl;
	//cout<<"("<<c0<<","<<c1<<")     Correction="<<Correction<<", XS="<<XSection_SignalOnly[c0][c1]<<endl;
	//XSection_SignalOnly[c0][c1]=UnfoldedData_SignalOnly[c0][c1]/Correction;
      }
    }

    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    (*MUnfolding)(c0*NBinsTrueAngle+c1,e0*NBinsRecAngle+e1)=vUnfolding[c0][c1][e0][e1];
	    (*MLikelihood)(c0*NBinsTrueAngle+c1,e0*NBinsRecAngle+e1)=vLikelihood[c0][c1][e0][e1];
	    (*VReconstructedEvents)(e0*NBinsRecAngle+e1)=DataReconstructedEvents[e0][e1];
	    //cout<<(*VReconstructedEvents)[e0*NBinsRecAngle+e1]<<endl;	    
	  }
	}
      }
    }
    
    //cout<<endl<<endl;
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	XSection_SignalOnly_Fill[c0][c1]=XSection_SignalOnly[c0][c1];
	TrueXSection_SignalOnly_Fill[c0][c1]=TrueXSection_SignalOnly[c0][c1];
	UnfoldedData_SignalOnly_Fill[c0][c1]=UnfoldedData_SignalOnly[c0][c1];
	table_TrueEventsDistribution_SignalOnly_Fill[c0][c1]=table_TrueEventsDistribution_SignalOnly[c0][c1];
	//cout<<setprecision(3)<<table_TrueEventsDistribution_SignalOnly[c0][c1]<<"X"<<UnfoldedData_SignalOnly[c0][c1]<<" ";
	 
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    DataReconstructedEvents_TrueSignal_Fill[c0][c1][e0][e1]=DataReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	  }
	}
      }
      //cout<<endl;
    }
    
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	DataReconstructedEvents_Fill[e0][e1] = DataReconstructedEvents[e0][e1];
      }
    }
    
    if(NIterationFluctuations) wtree->Fill();//In the case of check of the unfolding behaviour w/ number of iterations, the number of events is saved for each iterations
    }
    if(!(NIterationFluctuations) && !StatisticalUncertaintyEstimation) wtree->Fill();//In the normal case, the number of events is saved only at the end (after all iterations has been applied)
#ifdef DEBUG
    cout<<endl<<"Output: Number of events in each bin:"<<endl;
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	double Ratio = UnfoldedData_SignalOnly[c0][c1]-TrueEventsDistribution_SignalOnly[c0][c1];
	Ratio = Ratio/ ( TrueEventsDistribution_SignalOnly[c0][c1] == 0 ? 1 : TrueEventsDistribution_SignalOnly[c0][c1] );
	double RatioXS = XSection_SignalOnly[c0][c1]-TrueXSection_SignalOnly[c0][c1];
	RatioXS = RatioXS/ ( TrueXSection_SignalOnly[c0][c1] == 0 ? 1 : TrueXSection_SignalOnly[c0][c1] );
	
	cout<<"Events:"<<UnfoldedData_SignalOnly[c0][c1]<<", true="<<TrueEventsDistribution_SignalOnly[c0][c1]<<", post/true="<<Ratio*100<<"%"<<endl;
	cout<<"Efficiency: data="<<MCEfficiency[c0][c1]<<", mc="<<DataEfficiency[c0][c1]<<endl;
	cout<<"XS:"<<XSection_SignalOnly[c0][c1]<<", true="<<TrueXSection_SignalOnly[c0][c1]<<", post/true="<<RatioXS*100<<"%"<<endl;
      }
      cout<<endl;
    }
#endif
      }//end of loop on toys 1: statistical variations
    }//end of loop on toys 2: systematic variations
    if(!(NIterationFluctuations) && StatisticalUncertaintyEstimation) wtree->Fill();

#ifdef DEBUG
    /*
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
    */
#endif
  /*
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
	  double Value=UnfoldedData_SignalOnly[nt][c0][c1]-table_TrueEventsDistribution_SignalOnly[0][c0][c1];
	  if(table_TrueEventsDistribution_SignalOnly[0][c0][c1]!=0) Value/=table_TrueEventsDistribution_SignalOnly[0][c0][c1];
	  Events[c0][c1]->Fill(Value);

	  double BinningWidth=(BinningTrueMom[c0+1]-BinningTrueMom[c0])*(BinningTrueAngle[c1+1]-BinningTrueAngle[c1]);
	  double Flux=*NumberOfNeutrino;
	  double Efficiency=MCEfficiency[c0][c1];
	  double Correction=Flux*BinningWidth*Efficiency*NTarget;
	  //if(Correction!=0) Value/=Correction;
	  xsection[c0][c1]=UnfoldedData_SignalOnly[nt][c0][c1];
	  truedistribution[c0][c1]=hTrueDistribution->GetBinContent(c0+1,c1+1);
								   
	  hxsection[c0][c1]->Fill(xsection[c0][c1]);
	  htruedistribution->SetBinContent(c0+1,c1+1,truedistribution[c0][c1]);
	  
#ifdef DEBUG
	  cout<<"Number of events of the Toy="<<UnfoldedData_SignalOnly[nt][c0][c1]<<", true="<<table_TrueEventsDistribution_SignalOnly[0][c0][c1]<<", value="<<Value<<endl;
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
  //theApp.Run();
*/
  MUnfolding->Write("MUnfolding");
  MLikelihood->Write("MLikelihood");
  VReconstructedEvents->Write("VReconstructedEvents");
  for(int it=0;it<NToys[0];it++) hPrior[it]->Write(Form("Prior%d",it));
  
  wtree->Write();
  file->Write();
  file->Close();
#ifdef FLUXCONTROL
  fFlux->cd();
  hFlux->Write();
  fFlux->Close();
#endif
  //delete
  for(int h=0;h<NBinsTrueMom;h++){
    for(int i=0;i<NBinsTrueAngle;i++){
      for(int j=0;j<NBinsRecMom;j++){
	delete vLikelihood[h][i][j];
	delete vUnfolding[h][i][j];
	delete MCReconstructedEvents_TrueSignal[h][i][j];
	delete MCReconstructedEvents_TrueSignal_Default[h][i][j];
	delete DataReconstructedEvents_TrueSignal[h][i][j];
	delete DataReconstructedEvents_TrueSignal_Default[h][i][j];
	delete RelativeSigma[h][i][j];
      }      
      delete vLikelihood[h][i];
      delete vUnfolding[h][i];
      delete MCReconstructedEvents_TrueSignal[h][i];
      delete MCReconstructedEvents_TrueSignal_Default[h][i];
      delete DataReconstructedEvents_TrueSignal[h][i];
      delete DataReconstructedEvents_TrueSignal_Default[h][i];
      delete RelativeSigma[h][i];
    }
    delete vLikelihood[h];
    delete vUnfolding[h];
    delete MCReconstructedEvents_TrueSignal[h];
    delete MCReconstructedEvents_TrueSignal_Default[h];
    delete DataReconstructedEvents_TrueSignal[h];
    delete DataReconstructedEvents_TrueSignal_Default[h];
    delete RelativeSigma[h];
  }
  delete vLikelihood;
  delete vUnfolding;
  delete MCReconstructedEvents_TrueSignal;
  delete MCReconstructedEvents_TrueSignal_Default;
  delete DataReconstructedEvents_TrueSignal;
  delete DataReconstructedEvents_TrueSignal_Default;
  delete RelativeSigma;





  for(int i=0;i<NBinsTrueMom;i++){
    delete vInitialPriorMC[i];
    delete vInitialPrior[i];
    delete vPrior[i];
    delete vPriorNormalised[i];
    delete vPosterior[i];
    delete vPosterior_SignalOnly[i];
    delete DataEfficiency[i];
    delete MCEfficiency[i];
    delete TrueEventsDistribution_SignalOnly[i];
    delete UnfoldedData_SignalOnly[i];
    delete table_TrueEventsDistribution_SignalOnly[i];
    delete XSection_SignalOnly[i];
    delete TrueXSection_SignalOnly[i];
  }  
  delete vInitialPriorMC;
  delete vInitialPrior;
  delete vPrior;
  delete vPriorNormalised;
  delete vPosterior;
  delete vPosterior_SignalOnly;
  delete DataEfficiency;
  delete MCEfficiency;
  delete TrueEventsDistribution_SignalOnly;
  delete UnfoldedData_SignalOnly;
  delete table_TrueEventsDistribution_SignalOnly;
  delete XSection_SignalOnly;
  delete TrueXSection_SignalOnly;

 
  for(int i=0;i<NBinsRecMom;i++){
    delete DataReconstructedEvents[i];
    delete DataReconstructedEvents_Default[i];
    delete MCReconstructedBkgEvents[i];
    delete MCReconstructedEvents[i];
    delete TrueRecBins[i];
  }
  delete DataReconstructedEvents;
  delete DataReconstructedEvents_Default;
  delete MCReconstructedBkgEvents;
  delete MCReconstructedEvents;
  delete TrueRecBins;

  delete NumberOfNeutrino;

  return 0;
}

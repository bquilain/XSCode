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
  bool Systematics_Flux=false;
  bool Systematics_Xsec=false;
  bool Systematics_Detector=false;
 char InputNameMC[256];char InputNameData[256];char Name[256];char Title[256];
    
  while ((c = getopt(argc, argv, "fxd")) != -1) {
    switch(c){
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


  //0. Preparation for covariance estimation
  int NBins=NBinsTrueMom*NBinsTrueAngle;
  TH2D * Covariance = new TH2D("Covariance","",NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle,NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle);
  TH2D * NVariations = new TH2D("NVariations","",NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle,NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle);
  TH2D * CoRMS = new TH2D("CoRMS","",NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle,NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle);
  TH2D * Correlation = new TH2D("Correlation","",NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle,NBinsTrueMom*NBinsTrueAngle,0,NBinsTrueMom*NBinsTrueAngle);
  
  ifstream fEvent_Default;
  ifstream fEvent_Varied;
  double Value_Default[NBinsTrueMom][NBinsTrueAngle];
  double Value_Default_Data[NBinsTrueMom][NBinsTrueAngle];
  double Sigma[NBinsTrueMom][NBinsTrueAngle];
  double Value_Varied[NBinsTrueMom][NBinsTrueAngle];
  double Value_Varied_Data[NBinsTrueMom][NBinsTrueAngle];
  double MaximalDifference[NBinsTrueMom][NBinsTrueAngle];
  int check=0;

  TCanvas * canEvents[NBinsTrueMom][NBinsTrueAngle];
  TH1D * Events_MC[NBinsTrueMom][NBinsTrueAngle];
  TH1D * Events_Data[NBinsTrueMom][NBinsTrueAngle];
  TCanvas * canxsection[NBinsTrueMom][NBinsTrueAngle];
  TH1D * xsection_Data[NBinsTrueMom][NBinsTrueAngle];
  TH1D * xsection_MC[NBinsTrueMom][NBinsTrueAngle];
  

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	  int BinX=c0*NBinsTrueAngle+c1+1;
	  int BinY=d0*NBinsTrueAngle+d1+1;	
	  Covariance->SetBinContent(BinX,BinY,0);
	}
      }
    }
  }


  //1. Load the default Unfolding & intialize bin histograms
  sprintf(InputNameMC,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Systematics_Xsec_Error0_0Sigma.txt");
  //sprintf(InputNameMC,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Systematics_Flux0.txt");
  fEvent_Default.open(InputNameMC,ios::in);
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      fEvent_Default>>Value_Default[c0][c1];

      sprintf(Name,"Events_MC_bin%d_%d",c0,c1);
      sprintf(Title,"True distribution after unfolding in the bin %d_%d",c0,c1);
      Events_MC[c0][c1] = new TH1D(Name,Title,300,-3,3);
      
      sprintf(Name,"xsection_MC_bin%d_%d",c0,c1);
      sprintf(Title,"Double differential xsection after unfolding in the bin %d_%d",c0,c1);
      xsection_MC[c0][c1] = new TH1D(Name,Title,300,-3,3);
      
    }
  }
  fEvent_Default>>check;
  if(check!=666) cout<<"problem in reading data file"<<endl;
  fEvent_Default.close();

  sprintf(InputNameData,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Data.txt");
  fEvent_Default.open(InputNameData,ios::in);
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      fEvent_Default>>Value_Default_Data[c0][c1];

      sprintf(Name,"Events_Data_bin%d_%d",c0,c1);
      sprintf(Title,"True distribution after unfolding in the bin %d_%d",c0,c1);
      Events_Data[c0][c1] = new TH1D(Name,Title,300,-3,3);
      
      sprintf(Name,"xsection_Data_bin%d_%d",c0,c1);
      sprintf(Title,"Double differential xsection after unfolding in the bin %d_%d",c0,c1);
      xsection_Data[c0][c1] = new TH1D(Name,Title,300,-3,3);
      
    }
  }
  fEvent_Default>>check;
  if(check!=666) cout<<"problem in reading data file"<<endl;
  fEvent_Default.close();
  
  ////////////////////////////////////////////////////////////////////////////////
  
  
  
  if(Systematics_Flux){

    //2. For each different MC, check the variation => pull. Only one error source here:
    for(int n=0;n<1;n++){
      
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      int BinX=c0*NBinsTrueAngle+c1+1;
	      int BinY=d0*NBinsTrueAngle+d1+1;	
	      NVariations->SetBinContent(BinX,BinY,0);
	    }
	  }
	}
      }
      
      
      for(int m=0;m<NFluxFiles;m++){//number of files for this error source
      sprintf(InputNameMC,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Systematics_Flux%d.txt",m);
      fEvent_Varied.open(InputNameMC,ios::in);
	
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  fEvent_Varied>>Value_Varied[c0][c1];
	  //cout<<Value_Varied[c0][c1]<<endl;
	  Sigma[c0][c1]=(Value_Varied[c0][c1]-Value_Default[c0][c1]);
	  if(Value_Default[c0][c1]!=0) Sigma[c0][c1]/=Value_Default[c0][c1];
	  cout<<"Sigma="<<Sigma[c0][c1]<<endl;
	  xsection_MC[c0][c1]->Fill(Sigma[c0][c1]);
	  //xsection_MC[c0][c1]->Fill(Value);
	}
      }
      fEvent_Varied>>check;
      if(check!=666) cout<<"problem in reading data file"<<endl;
      fEvent_Varied.close();

      //3. For each error source, build the covariance matrix. Here, just add in quadrature
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      int BinX=c0*NBinsTrueAngle+c1+1;
	      int BinY=d0*NBinsTrueAngle+d1+1;	
	      double Cova=Covariance->GetBinContent(BinX,BinY);
	      Cova+=(Value_Varied[c0][c1]-Value_Default[c0][c1])*(Value_Varied[d0][d1]-Value_Default[d0][d1]);
	      Covariance->SetBinContent(BinX,BinY,Cova);

	      int nvar=NVariations->GetBinContent(BinX,BinY);
	      nvar++;
	      NVariations->SetBinContent(BinX,BinY,1);;
	    }
	  }
	}
      }
      }
      Covariance->Divide(NVariations);
    }
      
    //4. For the total error, build the correlation matrix from the covariant one
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	  for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	    int BinX=c0*NBinsTrueAngle+c1+1;
	    int BinY=d0*NBinsTrueAngle+d1+1;	
	    double Corr=Covariance->GetBinContent(BinX,BinY);
	    double NormCorr=TMath::Sqrt(Covariance->GetBinContent(BinX,BinX)*Covariance->GetBinContent(BinY,BinY));
	    if(NormCorr) Corr/=NormCorr;
	    Correlation->SetBinContent(BinX,BinY,Corr);
	    //CoRMS->SetBinContent(BinX,BinY,CoRMS);
	  }
	}
      }
    }

    TCanvas * cCova = new TCanvas();
    Covariance->Draw("colz");
    TCanvas * cCorr = new TCanvas();
    Correlation->Draw("colz");
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	
	sprintf(Name,"canEvents%d_%d",c0,c1);
	sprintf(Title,"True distribution of events after unfolding in the bin %d_%d",c0,c1);
	canEvents[c0][c1] = new TCanvas(Name,Title);
	Events_MC[c0][c1]->GetXaxis()->SetTitle("Pull");
	Events_MC[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	Events_MC[c0][c1]->Draw();
	Events_MC[c0][c1]->Fit("gaus","R");
	
	sprintf(Name,"canxsection%d_%d",c0,c1);
	sprintf(Title,"Double differential xsection after unfolding in the bin %d_%d",c0,c1);
	canxsection[c0][c1] = new TCanvas(Name,Title);
	xsection_MC[c0][c1]->GetXaxis()->SetTitle("Pull");
	xsection_MC[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	xsection_MC[c0][c1]->Draw();
	xsection_MC[c0][c1]->Fit("gaus","R");
      }
    }
  }
  else if(Systematics_Xsec){//asumming no correlation between different systs.
    
    //2. For each different MC, check the variation => pull
    for(int n=StartXsec;n<=EndXsec;n++){
      cout<<endl<<endl<<"Xsection tweak value="<<n<<endl;
      sprintf(InputNameMC,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Systematics_Xsec_Error%d_1Sigma.txt",n);
      fEvent_Varied.open(InputNameMC,ios::in);
	
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  fEvent_Varied>>Value_Varied[c0][c1];
	  Sigma[c0][c1]=(Value_Varied[c0][c1]-Value_Default[c0][c1]);
	  if(Value_Default[c0][c1]!=0) Sigma[c0][c1]/=Value_Default[c0][c1];
	  cout<<"Sigma="<<Sigma[c0][c1]*100<<"  ";
	  //xsection_MC[c0][c1]->Fill(Value);
	}
	cout<<endl;
      }
      fEvent_Varied>>check;
      if(check!=666) cout<<"problem in reading data file"<<endl;
      fEvent_Varied.close();

      //3. For each error source, build the covariance matrix. For now, we assume that error sources are uncorrelated between each other
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      int BinX=c0*NBinsTrueAngle+c1+1;
	      int BinY=d0*NBinsTrueAngle+d1+1;	
	      double Cova=Covariance->GetBinContent(BinX,BinY);
	      Cova+=(Value_Varied[c0][c1]-Value_Default[c0][c1])*(Value_Varied[d0][d1]-Value_Default[d0][d1]);
	      Covariance->SetBinContent(BinX,BinY,Cova);
	    }
	  }
	}
      }
    }

    //4. For the total error, build the correlation matrix from the covariant one
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	  for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	    int BinX=c0*NBinsTrueAngle+c1+1;
	    int BinY=d0*NBinsTrueAngle+d1+1;	
	    double Corr=Covariance->GetBinContent(BinX,BinY);
	    double NormCorr=TMath::Sqrt(Covariance->GetBinContent(BinX,BinX)*Covariance->GetBinContent(BinY,BinY));
	    if(NormCorr) Corr/=NormCorr;
	    Correlation->SetBinContent(BinX,BinY,Corr);
	    //CoRMS->SetBinContent(BinX,BinY,CoRMS);
	  }
	}
      }
    }

    TCanvas * cCova = new TCanvas();
    Covariance->Draw("colz");
    TCanvas * cCorr = new TCanvas();
    Correlation->Draw("colz");
    
  }
  else if(Systematics_Detector){
    //2. For each different MC, check the variation => pull. Only one error source here:
    for(int ErrorType=StartError;ErrorType<=EndError;ErrorType++){
      cout<<endl<<"Error Type="<<ErrorType<<endl;	
      
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      int BinX=c0*NBinsTrueAngle+c1+1;
	      int BinY=d0*NBinsTrueAngle+d1+1;	
	      NVariations->SetBinContent(BinX,BinY,0);
	    }
	  }
	}
      }
      
    
    //Assume non correlated error sources!
    //Don't forget to take the highest difference as an error
    //In the case of errors >6, (varied -nominal) becomes? That's quite the same in fact, except that we are looking no more MC_varied -MC_nominal but: (Data-MC)_varied-(Data-MC)_nominal. So one should only run also data and take the data -MC as difference. Not so simple because 100% is MC and data nominal... Which means that one should do: MC_varied/MC_nominal-Data_varied/Data_Nominal => It gives the relative error (not the absolute as for the other ones). So should be scaled (w/ MC, data events?)? To do!
    //Do things for all errors separately!
    //For errors >=8, 
    //1. For each error and its variation, read the MC and data files
    //2. For each error, among the iteration, take the largest variation. Input it in Covariance.
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	MaximalDifference[c0][c1]=0;
      }
    }

    if(ErrorType==3){//Efficiency error: fill Covariant matrix for each variation within the error
      for(int m=0;m<NE[ErrorType];m++){
	double ErrorValue=Start[ErrorType]+m*Step[ErrorType];
	cout<<"Error Value="<<ErrorValue<<", files generated:";
	
	sprintf(InputNameMC,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Systematics%d_%d.txt",ErrorType,m);
	fEvent_Varied.open(InputNameMC,ios::in);
	
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    fEvent_Varied>>Value_Varied[c0][c1];
	    //cout<<Value_Varied[c0][c1]<<endl;
	    Sigma[c0][c1]=(Value_Varied[c0][c1]-Value_Default[c0][c1]);
	    if(Value_Default[c0][c1]!=0) Sigma[c0][c1]/=Value_Default[c0][c1];
	    cout<<"Sigma="<<Sigma[c0][c1]<<endl;
	    xsection_MC[c0][c1]->Fill(Sigma[c0][c1]);
	    //xsection_MC[c0][c1]->Fill(Value);
	  }
	}
      fEvent_Varied>>check;
      if(check!=666) cout<<"problem in reading data file"<<endl;
      fEvent_Varied.close();

      //3. For each error source, build the covariance matrix. Here, just add in quadrature
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      int BinX=c0*NBinsTrueAngle+c1+1;
	      int BinY=d0*NBinsTrueAngle+d1+1;	
	      double Cova=Covariance->GetBinContent(BinX,BinY);
	      Cova+=(Value_Varied[c0][c1]-Value_Default[c0][c1])*(Value_Varied[d0][d1]-Value_Default[d0][d1]);
	      Covariance->SetBinContent(BinX,BinY,Cova);

	      int nvar=NVariations->GetBinContent(BinX,BinY);
	      nvar++;
	      NVariations->SetBinContent(BinX,BinY,1);;
	    }
	  }
	}
      }
      }
      Covariance->Divide(NVariations);
    }
    else{//Case where we take the maximal variations as the systematics. For error>6, we compare MC vs Data. Else, we compare MC with MC.
      for(int m=0;m<NE[ErrorType];m++){
	double ErrorValue=Start[ErrorType]+m*Step[ErrorType];
	cout<<"Error Value="<<ErrorValue<<", files generated:";
	////////////////////////READ MC////////////////////////////////
	sprintf(InputNameMC,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Systematics%d_%d.txt",ErrorType,m);
	fEvent_Varied.open(InputNameMC,ios::in);
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	    fEvent_Varied>>Value_Varied[c0][c1];
	  }
	}
	fEvent_Varied>>check;
	if(check!=666) cout<<"problem in reading mc file"<<endl;
	fEvent_Varied.close();
	
	if(ErrorType>6){ 
	  //////////////////////////READ DATA////////////////////////////////////
	  sprintf(InputNameData,"/home/bquilain/CC0pi_XS/XS/files/Unfolded_Data_Systematics%d_%d.txt",ErrorType,m);
	  fEvent_Varied.open(InputNameData,ios::in);
	  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	      fEvent_Varied>>Value_Varied_Data[c0][c1];
	    }
	  }
	  fEvent_Varied>>check;
	  if(check!=666) cout<<"problem in reading data file"<<endl;
	  fEvent_Varied.close();
	}
	////////////////////////////CHECK IF THIS IS THE LARGEST DIFFERENCE//////////////////////////////////////
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1

	    if(ErrorType>6){
	      double MCRatio=Value_Varied[c0][c1];
	      if(Value_Default[c0][c1]!=0) MCRatio/=Value_Default[c0][c1];
	      double DataRatio=Value_Varied_Data[c0][c1];
	      if(Value_Default_Data[c0][c1]!=0) DataRatio/=Value_Default_Data[c0][c1];
	      double max=TMath::Abs(DataRatio-MCRatio)*Value_Default[c0][c1];//Because we want not relative but absolute number to create the covariant matrix
	      if(max>MaximalDifference[c0][c1]) MaximalDifference[c0][c1]=max;
	    }
	    else{
	      double max=TMath::Abs(Value_Varied[c0][c1]-Value_Default[c0][c1]);
	      if(max>MaximalDifference[c0][c1]) MaximalDifference[c0][c1]=max;    
	    }
	    
	  }
	}
      }//End of variation of the error
      
      //3. For each error source, build the covariance matrix. Here, just add in quadrature
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	    for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	      int BinX=c0*NBinsTrueAngle+c1+1;
	      int BinY=d0*NBinsTrueAngle+d1+1;	
	      double Cova=Covariance->GetBinContent(BinX,BinY);
	      Cova+=(MaximalDifference[c0][c1]*Value_Default[c0][c1])*(Value_Varied[d0][d1]*Value_Default[d0][d1]);
	      Covariance->SetBinContent(BinX,BinY,Cova);
	      
	      int nvar=NVariations->GetBinContent(BinX,BinY);
	      nvar++;
	      NVariations->SetBinContent(BinX,BinY,1);;
	    }
	  }
	}
      }
    Covariance->Divide(NVariations);
    }
    }
  
    //4. For the total error, build the correlation matrix from the covariant one
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int d0=0;d0<NBinsTrueMom;d0++){//loop over cause 0
	  for(int d1=0;d1<NBinsTrueAngle;d1++){//loop over cause 1
	    int BinX=c0*NBinsTrueAngle+c1+1;
	    int BinY=d0*NBinsTrueAngle+d1+1;	
	    double Corr=Covariance->GetBinContent(BinX,BinY);
	    double NormCorr=TMath::Sqrt(Covariance->GetBinContent(BinX,BinX)*Covariance->GetBinContent(BinY,BinY));
	    if(NormCorr) Corr/=NormCorr;
	    Correlation->SetBinContent(BinX,BinY,Corr);
	    //CoRMS->SetBinContent(BinX,BinY,CoRMS);
	  }
	}
      }
    }

    TCanvas * cCova = new TCanvas();
    Covariance->Draw("colz");
    TCanvas * cCorr = new TCanvas();
    Correlation->Draw("colz");
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	
	sprintf(Name,"canEvents%d_%d",c0,c1);
	sprintf(Title,"True distribution of events after unfolding in the bin %d_%d",c0,c1);
	canEvents[c0][c1] = new TCanvas(Name,Title);
	Events_MC[c0][c1]->GetXaxis()->SetTitle("Pull");
	Events_MC[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	Events_MC[c0][c1]->Draw();
	Events_MC[c0][c1]->Fit("gaus","R");
	
	sprintf(Name,"canxsection%d_%d",c0,c1);
	sprintf(Title,"Double differential xsection after unfolding in the bin %d_%d",c0,c1);
	canxsection[c0][c1] = new TCanvas(Name,Title);
	xsection_MC[c0][c1]->GetXaxis()->SetTitle("Pull");
	xsection_MC[c0][c1]->GetYaxis()->SetTitle("Number of toy experiments");
	xsection_MC[c0][c1]->Draw();
	xsection_MC[c0][c1]->Fit("gaus","R");
      }
    }
  }

    
  theApp.Run();
  
  return 0;
}

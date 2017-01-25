//CAREFUL: ONLY STOPPED SAMPLE RIGHT NOW (SAMPLE==3)
//For no, we have a bug correction in this code, that set CL=1e-30 if this is 0. But try to avoid this if possible -> Check the CUmulative distributions
//For now, isreconstructed is not systematically applied
//For now, I removed the track witdth and matching parameters!

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
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
//#include <TApplication.h>
#include <TBox.h>
#include <TVectorD.h>
#include "setup.h"
//#include "Xsec.cc"
//#include "Reconstruction.cc"
//#define DEBUG

//Xsec * XS = new Xsec();
char Type[32];
char Name[256];
char Name0[256];
double ScalingMC;
double DataEquivalent=1;
double SandReweight=0.6;

void ProduceStack(TH1D * h[NFSIs], THStack * hStack){
  
  for(int fsi=0;fsi<NFSIs;fsi++){
    
    h[fsi]->GetYaxis()->SetTitleOffset(1.3);
    h[fsi]->Scale((5.86e20/1e21)*100/1000);
    
    if(fsi<3){//CC0Pi+/-/0
      if(fsi==0) h[0]->SetFillColor(kRed);
      else if(fsi==1) h[1]->SetFillColor(kOrange+10);
      else h[2]->SetFillColor(kOrange+9);
    }
    else{
      h[3]->SetFillColor(kBlue+2);//CC1Pi+/-
      h[4]->SetFillColor(kAzure+7);//CCNPi+/-
      h[5]->SetFillColor(kAzure+10);//CCpi0
      h[6]->SetFillColor(kGray);//NC
      if(fsi>6 && fsi<11){
	if(fsi==7) h[7]->SetFillColor(kYellow);//Sand
	if(fsi==8) h[8]->SetFillColor(kYellow+4);//Horizontal||Vertical Module Bkg
	if(fsi==9) h[9]->SetFillColor(kMagenta);//AntiNuMu
	if(fsi==10) h[10]->SetFillColor(kGreen);//NuE
      }
    }
    h[fsi]->SetLineColor(1);
    h[fsi]->SetLineWidth(2);
    hStack->Add(h[fsi]);
  }

}

void InitialiseTable(double DataSelected[NBinsRecMom][NBinsRecAngle],double MCSelected[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle],double BkgSelected[NBinsRecMom][NBinsRecAngle],double Efficiency[NBinsTrueMom][NBinsTrueAngle],double TotalCC0piEvent[NBinsTrueMom][NBinsTrueAngle]){
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      Efficiency[c0][c1]=0;
      TotalCC0piEvent[c0][c1]=0;
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  MCSelected[c0][c1][e0][e1]=0;
	  BkgSelected[e0][e1]=0;
	  DataSelected[e0][e1]=0;
	}
      }
    }
  }
}

void CC0piDistributions(TChain * wtree,bool IsData,int Selection,bool Plots,bool Systematics_Flux,int File_Number,TVectorD FluxVector,bool Systematics_Xsec,int dial,bool Systematics_Detector,int ErrorType,int ErrorIteration){
  cout<<"hello"<<endl;
  int nevt=(int) wtree->GetEntries();
  cout<<"number of events="<<nevt<<endl;
  int err;
  int sig;
  TH1D * NeutrinoFlux;
  TRandom3 * rand = new TRandom3();
  ///////////////////////////////LOAD THE NEUTRINO FLUX
  TFile * f = new TFile("files/nd34_tuned_11bv3.1_250ka.root");
  NeutrinoFlux = (TH1D*) f->Get("ing3_tune_numu");
  NeutrinoFlux->SetDirectory(0);//Detach the histogram from the TFile (one can close the file)
  f->Close();

  
  ///////////////////////////////NOW, TAKE CARE OF THE SELECTION
  int FSIInt;
  int Num_Int;
  int nTracks[10];
  double weight;
  bool IsFV;
  bool IsSand;
  bool IsAnti;
  bool IsBkgH;
  bool IsBkgV;
  bool IsNuE;
  bool IsDetected[10];
  bool SelectionFV[10];bool SelectionOV[10];
  double Enu;
  int nIngBasRec;
  double Errorweight;
  double TrueMomentumMuon, TrueAngleMuon;
  double POT;
  double TrackAngle[10][20];
  int TypeOfTrack[10][20];
  double CLMuon[10][20];// = new vector<double> [10][20];
  double CLMuon_Plan[10][20];// = new vector<double> [10][20];
  double CLMuon_Likelihood[10][20];// = new vector<double> [10][20];
  double ProportionHighPE[10][20];
  double MeanHighPE[10][20];
  double HighestPE[10][20];
  double Momentum[10][20];
  double IronDistance[10][20];
  double PlasticDistance[10][20];
  int Sample[10][20];
  bool IsReconstructed[10][20];
  double ReWeight[175];
  double CriteriaAngleX[10][20];
  double CriteriaAngleY[10][20];
  double CriteriaHalfWayX[10][20];
  double CriteriaHalfWayY[10][20];
  double TrackWidth[10][20];
  int Spill;
  int GoodSpill;


  /////////////////////////////////PREPARE FOR READING/////////////////////////////////////////
  TBranch* Br_Spill = wtree->GetBranch("Spill");
  Br_Spill->SetAddress(&Spill);
  wtree->SetBranchAddress("Spill",&Spill);

  TBranch* Br_GoodSpill = wtree->GetBranch("GoodSpill");
  Br_GoodSpill->SetAddress(&GoodSpill);
  wtree->SetBranchAddress("GoodSpill",&GoodSpill);

  TBranch* Br_FSIInt = wtree->GetBranch("FSIInt");
  Br_FSIInt->SetAddress(&FSIInt);
  wtree->SetBranchAddress("FSIInt",&FSIInt);

  TBranch* Br_nIngBasRec = wtree->GetBranch("nIngBasRec");
  Br_nIngBasRec->SetAddress(&nIngBasRec);
  wtree->SetBranchAddress("nIngBasRec",&nIngBasRec);

  TBranch * Br_Num_Int = wtree->GetBranch("InteractionType");
  Br_Num_Int->SetAddress(&Num_Int);
  wtree->SetBranchAddress("InteractionType",&Num_Int);
  
  TBranch * Br_nTracks = wtree->GetBranch("nTracks[10]");
  Br_nTracks->SetAddress(&nTracks);
  wtree->SetBranchAddress("nTracks[10]",&nTracks);
 
  TBranch * Br_weight = wtree->GetBranch("weight");
  Br_weight->SetAddress(&weight);
  wtree->SetBranchAddress("weight",&weight);

  TBranch * Br_Enu = wtree->GetBranch("Enu");
  Br_Enu->SetAddress(&Enu);
  wtree->SetBranchAddress("Enu",&Enu);

  TBranch * Br_TrueMomentumMuon = wtree->GetBranch("TrueMomentumMuon");
  Br_TrueMomentumMuon->SetAddress(&TrueMomentumMuon);
  wtree->SetBranchAddress("TrueMomentumMuon",&TrueMomentumMuon);

  TBranch * Br_TrueAngleMuon = wtree->GetBranch("TrueAngleMuon");
  Br_TrueAngleMuon->SetAddress(&TrueAngleMuon);
  wtree->SetBranchAddress("TrueAngleMuon",&TrueAngleMuon);

  TBranch * Br_IsFV = wtree->GetBranch("IsFV");
  Br_IsFV->SetAddress(&IsFV);
  wtree->SetBranchAddress("IsFV",&IsFV);

  TBranch * Br_IsSand = wtree->GetBranch("IsSand");
  Br_IsSand->SetAddress(&IsSand);
  wtree->SetBranchAddress("IsSand",&IsSand);
  
  TBranch * Br_IsAnti = wtree->GetBranch("IsAnti");
  Br_IsAnti->SetAddress(&IsAnti);
  wtree->SetBranchAddress("IsAnti",&IsAnti);
  
  TBranch * Br_IsBkgH = wtree->GetBranch("IsBkgH");
  Br_IsBkgH->SetAddress(&IsBkgH);
  wtree->SetBranchAddress("IsBkgH",&IsBkgH);

  TBranch * Br_IsBkgV = wtree->GetBranch("IsBkgV");
  Br_IsBkgV->SetAddress(&IsBkgV);
  wtree->SetBranchAddress("IsBkgV",&IsBkgV);

  TBranch * Br_IsNuE = wtree->GetBranch("IsNuE");
  Br_IsNuE->SetAddress(&IsNuE);
  wtree->SetBranchAddress("IsNuE",&IsNuE);

  TBranch * Br_SelectionFV = wtree->GetBranch("SelectionFV[10]");
  Br_SelectionFV->SetAddress(&SelectionFV);
  wtree->SetBranchAddress("SelectionFV[10]",&SelectionFV);

  TBranch * Br_SelectionOV = wtree->GetBranch("SelectionOV[10]");
  Br_SelectionOV->SetAddress(&SelectionOV);
  wtree->SetBranchAddress("SelectionOV[10]",&SelectionOV);

  TBranch * Br_IsDetected = wtree->GetBranch("IsDetected[10]");
  Br_IsDetected->SetAddress(&IsDetected);
  wtree->SetBranchAddress("IsDetected[10]",&IsDetected);
  
  TBranch * Br_TrackAngle = wtree->GetBranch("TrackAngle[10][20]");
  Br_TrackAngle->SetAddress(TrackAngle);
  wtree->SetBranchAddress("TrackAngle[10][20]",TrackAngle);
 
 TBranch * Br_TypeOfTrack = wtree->GetBranch("TypeOfTrack[10][20]");
  Br_TypeOfTrack->SetAddress(&TypeOfTrack);
  wtree->SetBranchAddress("TypeOfTrack[10][20]",&TypeOfTrack);

  TBranch * Br_CLMuon = wtree->GetBranch("CLMuon[10][20]");
  Br_CLMuon->SetAddress(&CLMuon);
  wtree->SetBranchAddress("CLMuon[10][20]",&CLMuon);

  TBranch * Br_CLMuon_Likelihood = wtree->GetBranch("CLMuon_Likelihood[10][20]");
  Br_CLMuon_Likelihood->SetAddress(&CLMuon_Likelihood);
  wtree->SetBranchAddress("CLMuon_Likelihood[10][20]",&CLMuon_Likelihood);

  TBranch * Br_CLMuon_Plan = wtree->GetBranch("CLMuon_Plan[10][20]");
  Br_CLMuon_Plan->SetAddress(&CLMuon_Plan);
  wtree->SetBranchAddress("CLMuon_Plan[10][20]",&CLMuon_Plan);

  TBranch * Br_ProportionHighPE = wtree->GetBranch("ProportionHighPE[10][20]");
  Br_ProportionHighPE->SetAddress(&ProportionHighPE);
  wtree->SetBranchAddress("ProportionHighPE[10][20]",&ProportionHighPE);

    TBranch * Br_MeanHighPE = wtree->GetBranch("MeanHighPE[10][20]");
  Br_MeanHighPE->SetAddress(&MeanHighPE);
  wtree->SetBranchAddress("MeanHighPE[10][20]",&MeanHighPE);

    TBranch * Br_HighestPE = wtree->GetBranch("HighestPE[10][20]");
  Br_HighestPE->SetAddress(&HighestPE);
  wtree->SetBranchAddress("HighestPE[10][20]",&HighestPE); 
    
  TBranch * Br_Momentum = wtree->GetBranch("Momentum[10][20]");
  Br_Momentum->SetAddress(Momentum);
  wtree->SetBranchAddress("Momentum[10][20]",Momentum);

  TBranch * Br_IronDistance = wtree->GetBranch("IronDistance[10][20]");
  Br_IronDistance->SetAddress(IronDistance);
  wtree->SetBranchAddress("IronDistance[10][20]",IronDistance);

  TBranch * Br_PlasticDistance = wtree->GetBranch("PlasticDistance[10][20]");
  Br_PlasticDistance->SetAddress(PlasticDistance);
  wtree->SetBranchAddress("PlasticDistance[10][20]",PlasticDistance);

  TBranch * Br_Sample = wtree->GetBranch("Sample[10][20]");
  Br_Sample->SetAddress(Sample);
  wtree->SetBranchAddress("Sample[10][20]",Sample);

  TBranch * Br_IsReconstructed = wtree->GetBranch("IsReconstructed[10][20]");
  Br_IsReconstructed->SetAddress(IsReconstructed);
  wtree->SetBranchAddress("IsReconstructed[10][20]",IsReconstructed);

  TBranch * Br_CriteriaAngleX = wtree->GetBranch("CriteriaAngleX[10][20]");
  Br_CriteriaAngleX->SetAddress(CriteriaAngleX);
  wtree->SetBranchAddress("CriteriaAngleX[10][20]",CriteriaAngleX);

  TBranch * Br_CriteriaAngleY = wtree->GetBranch("CriteriaAngleY[10][20]");
  Br_CriteriaAngleY->SetAddress(CriteriaAngleY);
  wtree->SetBranchAddress("CriteriaAngleY[10][20]",CriteriaAngleY);

  TBranch * Br_CriteriaHalfWayX = wtree->GetBranch("CriteriaHalfWayX[10][20]");
  Br_CriteriaHalfWayX->SetAddress(CriteriaHalfWayX);
  wtree->SetBranchAddress("CriteriaHalfWayX[10][20]",CriteriaHalfWayX);

  TBranch * Br_CriteriaHalfWayY = wtree->GetBranch("CriteriaHalfWayY[10][20]");
  Br_CriteriaHalfWayY->SetAddress(CriteriaHalfWayY);
  wtree->SetBranchAddress("CriteriaHalfWayY[10][20]",CriteriaHalfWayY);

  TBranch * Br_ReWeight = wtree->GetBranch("ReWeight[175]");
  Br_ReWeight->SetAddress(ReWeight);
  wtree->SetBranchAddress("ReWeight[175]",ReWeight);

  TBranch * Br_POT = wtree->GetBranch("POT");
  Br_POT->SetAddress(&POT);
  wtree->SetBranchAddress("POT",&POT);

  TBranch * Br_TrackWidth = wtree->GetBranch("TrackWidth[10][20]");
  Br_TrackWidth->SetAddress(TrackWidth);
  wtree->SetBranchAddress("TrackWidth[10][20]",TrackWidth);
  /*
  TH1D * FluxDistribution = new TH1D("FluxDistribution","",NBinsEnergyFlux,BinningEnergyFlux);
  for(int i=0;i<NBinsEnergyFlux;i++){
    if(Systematics_Flux) FluxDistribution->SetBinContent(i+1,hNormNuMu_Binned->GetBinContent(i+1)+hNormNuMu_Binned->GetBinContent(i+1)*FluxVector[i]);
    else FluxDistribution->SetBinContent(i+1,0.);
  }
*/
  
  int Counter=0;
  int Counter2=0;
  int Counter3=0;
  double POTCount=0;

  double DataSelected[NBinsRecMom][NBinsRecAngle];
  double MCSelected[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double BkgSelected[NBinsRecMom][NBinsRecAngle];
  double Efficiency[NBinsTrueMom][NBinsTrueAngle];
  double TotalCC0piEvent[NBinsTrueMom][NBinsTrueAngle];
  
  InitialiseTable(DataSelected,MCSelected,BkgSelected,Efficiency,TotalCC0piEvent);
  
  TH2D * MCTrueEvents = new TH2D("MCTrueEvents","",NBinsTrueMom,BinningTrueMom,NBinsTrueAngle,BinningTrueAngle);
  
  TH1D * hMuCL[NFSIs];
  TH1D * hMuCL_Lowest[NFSIs];
  TH2D * hMuCL_2tracks[NFSIs];
  TH1D * hNTracks[NFSIs];
  TH1D * hRecMom[NFSIs];
  TH1D * hRecAngle[NFSIs];

  TH1D * hRecMom_CC0pi[NFSIs];
  TH1D * hRecAngle_CC0pi[NFSIs];

  TH1D * hRecMom_CC1pi[NFSIs];
  TH1D * hRecAngle_CC1pi[NFSIs];

  TH1D * hRecMom_CCNpi[NFSIs];
  TH1D * hRecAngle_CCNpi[NFSIs];

  TH1D * hMuCL_TrueMuon = new TH1D("hMuCL_TrueMuon","Muon confidence level for the sand muons",100,0,1);
  TH1D * hMuCL_TruePion = new TH1D("hMuCL_TruePion","Muon confidence level for the true muons",100,0,1);
  TH1D * hMuCL_TrueProton = new TH1D("hMuCL_TrueProton","Muon confidence level for the true protons",100,0,1);

  TH2D * PE_Lowest_CC0pi = new TH2D("PE_Lowest_CC0pi","",50,0,1,200,0,1000);
  TH2D * PE_Lowest_Other = new TH2D("PE_Lowest_Other","",50,0,1,200,0,1000);
				    
  if(Plots){
      char Name[256];char Title[256];

      for(int i=0;i<NFSIs;i++){
	sprintf(Name,"hMuCL%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi",i);
	hMuCL[i] = new TH1D(Name,Title,100,0,1);
	//hMuCL->Sumw2();
	hMuCL[i]->GetXaxis()->SetTitle("#mu_{CL}");    hMuCL[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hMuCL_Lowest%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi, only for the lowest MuCL tracks if ntracks>1",i);
	hMuCL_Lowest[i] = new TH1D(Name,Title,100,0,1);
	//hMuCL_Lowest->Sumw2();
	hMuCL_Lowest[i]->GetXaxis()->SetTitle("#mu_{CL_Lowest}");    hMuCL_Lowest[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hMuCL_2tracks%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi",i);
	hMuCL_2tracks[i] = new TH2D(Name,Title,100,0,1,100,0,1);
	//hMuCL_2tracks->Sumw2();
	hMuCL_2tracks[i]->GetXaxis()->SetTitle("higher #mu_{CL}");    hMuCL_2tracks[i]->GetYaxis()->SetTitle("lower #mu_{CL}");
	
	sprintf(Name,"hNTracks%d",i);
	sprintf(Title,"Number of tracks for the %dth-fsi",i);
	hNTracks[i] = new TH1D(Name,Title,5,0,5);
	//hNTracks->Sumw2();
	hNTracks[i]->GetXaxis()->SetTitle("Number of reconstructed tracks");
	hNTracks[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecMom%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom[i] = new TH1D(Name,Title,NBinsRecMom,BinningRecMom);
	//hRecMom->Sumw2();
	hRecMom[i]->GetXaxis()->SetTitle("Penetration in iron (cm)");
	hRecMom[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle[i] = new TH1D(Name,Title,NBinsRecAngle,BinningRecAngle);
	//hRecAngle->Sumw2();
	hRecAngle[i]->GetXaxis()->SetTitle("Angle (cm)");
	hRecAngle[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecMom_CC0pi%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CC0pi[i] = new TH1D(Name,Title,NBinsRecMom,BinningRecMom);
	//hRecMom_CC0pi->Sumw2();
	hRecMom_CC0pi[i]->GetXaxis()->SetTitle("Penetration in iron (cm)");
	hRecMom_CC0pi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CC0pi%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CC0pi[i] = new TH1D(Name,Title,NBinsRecAngle,BinningRecAngle);
	//hRecAngle_CC0pi->Sumw2();
	hRecAngle_CC0pi[i]->GetXaxis()->SetTitle("Angle_CC0pi (cm)");
	hRecAngle_CC0pi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecMom_CC1pi%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CC1pi[i] = new TH1D(Name,Title,NBinsRecMom,BinningRecMom);
	//hRecMom_CC1pi->Sumw2();
	hRecMom_CC1pi[i]->GetXaxis()->SetTitle("Penetration in iron (cm)");
	hRecMom_CC1pi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CC1pi%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CC1pi[i] = new TH1D(Name,Title,NBinsRecAngle,BinningRecAngle);
	//hRecAngle_CC1pi->Sumw2();
	hRecAngle_CC1pi[i]->GetXaxis()->SetTitle("Angle_CC1pi (cm)");
	hRecAngle_CC1pi[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hRecMom_CCNpi%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CCNpi[i] = new TH1D(Name,Title,NBinsRecMom,BinningRecMom);
	//hRecMom_CCNpi->Sumw2();
	hRecMom_CCNpi[i]->GetXaxis()->SetTitle("Penetration in iron (cm)");
	hRecMom_CCNpi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CCNpi%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CCNpi[i] = new TH1D(Name,Title,NBinsRecAngle,BinningRecAngle);
	//hRecAngle_CCNpi->Sumw2();
	hRecAngle_CCNpi[i]->GetXaxis()->SetTitle("Angle_CCNpi (cm)");
	hRecAngle_CCNpi[i]->GetYaxis()->SetTitle("Number of events");

      }

  }    
  
  //cout<<"going to read events"<<endl;
  for(int ievt=0;ievt<nevt;ievt++){
    POTCount+=POT;
    //if(ievt%10000==0){cout<<ievt<<", "<<"POT processed="<<POTCount<<", xsec="<<IsXsec<<endl;
    if(ievt==(nevt-1)) cout<<"Final POT="<<POTCount<<endl;
    wtree->GetEntry(ievt);
    if(!IsData) weight*=ScalingMC;//Adjust the MC distribution to the amount of data we process
    if(IsSand) weight=(1+SandReweight);
    ////////////////DETERMINE THE REWEIGHTING OF THE EVENT IF SYSTEMATICS ERROR FLUX
    if(Systematics_Flux){
      for(int i=0;i<NBinsEnergyFlux;i++){
	if(Enu<BinningEnergyFlux[i+1]) {
	  Errorweight=FluxVector[i];
	  break;
	}
      }
      weight=weight+weight*Errorweight;
    }
    else if(Systematics_Xsec) weight=weight*ReWeight[dial];
    //cout<<ReWeight[5]<<endl;
    
    bool IsNuMu=((!IsSand) && (!IsAnti) && (!IsNuE) && (!IsBkgH) && (!IsBkgV));
    if(IsData && (!GoodSpill || !Spill)) continue;
    if((IsFV && !IsData) || IsData){
#ifdef DEBUG
      cout<<"FV"<<", nbasrec="<<nIngBasRec<<endl;
#endif
      int BinTrueMom=0;
      int BinTrueAngle=0;

      for(int i=0;i<=NBinsTrueMom;i++){
	if(TrueMomentumMuon<BinningTrueMom[i+1]){BinTrueMom=i;break;}
      }
      for(int i=0;i<=NBinsTrueAngle;i++){
	if(TrueAngleMuon<BinningTrueAngle[i+1]){BinTrueAngle=i;break;}
      }
      //cout<<FSIInt<<", "<<Num_Int<<endl;
      if(FSIInt<3 && IsNuMu){
	TotalCC0piEvent[BinTrueMom][BinTrueAngle]+=weight;
	MCTrueEvents->Fill(BinTrueMom+1,BinTrueAngle+1,weight);
      }
      for(int irec=0;irec<nIngBasRec;irec++){

	bool MuonFound=false;
	int MuonTrue;int MuonRec=0;
	int LowestMuCL=0;
	bool Trash=false;
	if(IsData) FSIInt=0;
	if(Num_Int==1) FSIInt=1;
	else if(Num_Int==2) FSIInt=2;
	else if(Num_Int>2 && Num_Int<20) FSIInt=3;
	else if(Num_Int<30) FSIInt=4;
	else if(Num_Int>=30) FSIInt=5;
   
	if(Selection==1 && !SelectionFV[irec]) continue;
	if(Selection==0 && !SelectionOV[irec]) continue;
	if(!IsDetected[irec]) continue;

	  ////////////////////DETERMINE MUON TRACK/////////////////////////////////
	  for(int itrk=0;itrk<nTracks[irec];itrk++){
	    //if(CLMuon_Likelihood[irec][itrk]==0) CLMuon_Likelihood[irec][itrk]=1e-30;
	    if(IsReconstructed[irec][itrk]){
	      //CLMuon_Plan[irec][itrk]=TMath::Log(CLMuon_Plan[irec][itrk]);
	      //if(CLMuon_Plan[irec][itrk]<-50) CLMuon_Plan[irec][itrk]=-50;
	     
	      //cout<<TypeOfTrack[irec][itrk]<<endl;
	      if(CLMuon_Likelihood[irec][itrk]!=CLMuon_Likelihood[irec][itrk] || CLMuon_Likelihood[irec][itrk]==-1){ Trash=true; continue;}
	      if(TypeOfTrack[irec][itrk]==13) MuonTrue=itrk;
	      if(CLMuon_Likelihood[irec][itrk]>=CLMuon_Likelihood[irec][MuonRec]) MuonRec=itrk;
	      MuonFound=true;
	      if(CLMuon_Likelihood[irec][itrk]<CLMuon_Likelihood[irec][LowestMuCL]) LowestMuCL=itrk;
	      
	      if(Plots){
		//if(SelectionOV[irec]){//MC case: muon at 99%. If data, one should only take long tracks since gamma contamination
		//cout<<CLMuon_Likelihood[irec][itrk]<<", "<<TMath::Abs(TypeOfTrack[irec][itrk])<<endl;
		if(SelectionFV[irec] && !IsData){
		  if(TMath::Abs(TypeOfTrack[irec][itrk])==13) hMuCL_TrueMuon->Fill(CLMuon_Likelihood[irec][itrk],weight);
		  else if(TMath::Abs(TypeOfTrack[irec][itrk])==211) hMuCL_TruePion->Fill(CLMuon_Likelihood[irec][itrk],weight);
		  else if(TMath::Abs(TypeOfTrack[irec][itrk])==2212) hMuCL_TrueProton->Fill(CLMuon_Likelihood[irec][itrk],weight);
		  //cout<<"hello, particle is="<<TMath::Abs(TypeOfTrack[irec][itrk])<<endl;
		}
	      }
	    }
	  }


	  
	  if(!MuonFound){cout<<"whether NC, whether bug to solve"<<endl;continue;}
	    double MuonLike=0;double ProtonLike=0;
	    for(int itrk=0;itrk<nTracks[irec];itrk++){
	      if(CLMuon_Likelihood[irec][itrk]>0.7) MuonLike++;
	      else if(CLMuon_Likelihood[irec][itrk]<0.3) ProtonLike++;
	    }
	    double rnd=rand->Uniform(0,1);


	    
	    if(Plots){
	      if(Sample[irec][MuonRec]==4 /*|| (Sample[irec][MuonRec]>=4 && rnd<=0.1)*/){
		//cout<<CLMuon_Likelihood[irec][MuonRec]<<endl;
		if(nTracks[irec]==1) hMuCL[FSIInt]->Fill(CLMuon_Likelihood[irec][MuonRec],weight);
		if(nTracks[irec]==2) hMuCL_2tracks[FSIInt]->Fill(CLMuon_Likelihood[irec][MuonRec],CLMuon_Likelihood[irec][LowestMuCL],weight);
		if(nTracks[irec]==2 && CLMuon_Likelihood[irec][MuonRec]>0.7){
		  hMuCL_Lowest[FSIInt]->Fill(CLMuon_Likelihood[irec][LowestMuCL],weight);
		  if(CLMuon_Likelihood[irec][LowestMuCL]<0.1){
		    if(FSIInt<3) PE_Lowest_CC0pi->Fill(ProportionHighPE[irec][LowestMuCL],MeanHighPE[irec][LowestMuCL],weight);
		    else PE_Lowest_Other->Fill(ProportionHighPE[irec][LowestMuCL],MeanHighPE[irec][LowestMuCL],weight);		    
		  }
		}
		if(MuonLike==1){
		  
		}
		hNTracks[FSIInt]->Fill(nTracks[irec],weight);
		hRecMom[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
		hRecAngle[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	      }
	    }
	    //if((nTracks[irec]==1 && CLMuon_Likelihood[irec][0]>-3) || ((nTracks[irec]==2 && TMath::Max(CLMuon_Likelihood[irec][0],CLMuon_Likelihood[irec][1])>-3 && TMath::Min(CLMuon_Likelihood[irec][0],CLMuon_Likelihood[irec][1])<-3) && (TMath::Min(CLMuon_Likelihood[irec][0],CLMuon_Likelihood[irec][1])>=0)) || ((nTracks[irec]==3 && CLMuon_Likelihood[irec][MuonRec]>-3 && CLMuon_Likelihood[irec][(MuonRec+1)%3]<-3 && CLMuon_Likelihood[irec][(MuonRec+2)%3]<-3)) || ((nTracks[irec]==4 && CLMuon_Likelihood[irec][MuonRec]>-3 && CLMuon_Likelihood[irec][(MuonRec+1)%4]<-3 && CLMuon_Likelihood[irec][(MuonRec+2)%4]<-3 && CLMuon_Likelihood[irec][(MuonRec+3)%4]<-3))){
	      
	    //if(TrackWidth[irec][MuonRec]>1.1 && IronDistance[irec][MuonRec]<15) continue;
	    //else if(TrackWidth[irec][MuonRec]>1.2 && IronDistance[irec][MuonRec]<20) continue;
	    //else if(TrackWidth[irec][MuonRec]>1.3 && IronDistance[irec][MuonRec]<30) continue;
	    //else if(TrackWidth[irec][MuonRec]>1.5) continue;
	    //cout<<weight<<endl;
	    int BinRecMom=0;
	    int BinRecAngle=0;
	    bool old=true;
	    
	    for(int i=0;i<=NBinsRecMom;i++){
	      if((IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio))<BinningRecMom[i+1]){BinRecMom=i;break;}
	    }
	    for(int i=0;i<=NBinsRecAngle;i++){
	      if(TrackAngle[irec][MuonRec]<BinningRecAngle[i+1]){BinRecAngle=i;break;}
	    }

	  ///////////////////THE ACTUAL SELECTION///////////////////////////////////////
	    //Note: should I add "reconstructed" to the selection??	    	
	    if(Sample[irec][MuonRec]!=4 /*|| (Sample[irec][MuonRec]>=4 && rnd>0.1)*/) continue;
    
      	    if(MuonLike==1/*&& nTracks[irec]<=2*/){//CC0pi

	      DataSelected[BinRecMom][BinRecAngle]+=weight;
	      if(!IsData){
		if(FSIInt<3 && IsNuMu){
		  MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
		  Efficiency[BinTrueMom][BinTrueAngle]+=weight;
		}
		else{
		  BkgSelected[BinRecMom][BinRecAngle]+=weight;
		}
	      }
	      hRecMom_CC0pi[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	      hRecAngle_CC0pi[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);

	    }
	    else if(MuonLike==2){//Side band CC1pi
	      hRecMom_CC1pi[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	      hRecAngle_CC1pi[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	    }
	    else if(MuonLike>=3){//Side band CCNpi
	      hRecMom_CCNpi[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	      hRecAngle_CCNpi[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	    }	     
	    }
      }
    }
  
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      if(!IsData && TotalCC0piEvent[c0][c1]!=0) Efficiency[c0][c1]/=TotalCC0piEvent[c0][c1];
    }
  }
    
  char OutNameEvent[256];
  char OutNameEvent_root[256];
  ofstream fEvent;
  if(Systematics_Flux){
    sprintf(OutNameEvent,"files/MCSelected_Systematics_Flux%d.txt",File_Number);
    sprintf(OutNameEvent_root,"files/MCSelected_Systematics_Flux%d.root",File_Number);
    //sprintf(OutNameEvent,"files/MCSelected_Systematics%d_%d.txt",ErrorType,ErrorIteration);
    //sprintf(OutNameEvent_root,"files/MCSelected_Systematics%d_%d.root",ErrorType,ErrorIteration);
    for(int i=1;i<=NBinsEnergyFlux;i++){
      double Nneut=NeutrinoFlux->GetBinContent(i);
      NeutrinoFlux->SetBinContent(i,Nneut+Nneut*FluxVector[i-1]);
    }
  }
  
  else if(Systematics_Xsec){
    int Var=(int) (dial/NXsecVariations);
    int Sig=dial%NXsecVariations-CenterXsecVariations;
    //sprintf(OutNameEvent,"files/MCSelected_Systematics%d_%d.txt",ErrorType,ErrorIteration);
    //sprintf(OutNameEvent_root,"files/MCSelected_Systematics%d_%d.root",ErrorType,ErrorIteration);
    sprintf(OutNameEvent,"files/MCSelected_Systematics_Xsec_Error%d_%dSigma.txt",Var,Sig);
    sprintf(OutNameEvent_root,"files/MCSelected_Systematics_Xsec_Error%d_%dSigma.root",Var,Sig);
  }
   else if(Systematics_Detector){
     sprintf(OutNameEvent,"files/MCSelected_Systematics%d_%d.txt",ErrorType,ErrorIteration);
     sprintf(OutNameEvent_root,"files/MCSelected_Systematics%d_%d.root",ErrorType,ErrorIteration);
  }
  else{
    if(IsData) sprintf(OutNameEvent,"files/DataSelected.txt");
    else{
      sprintf(OutNameEvent,"files/MCSelected.txt");
      sprintf(OutNameEvent_root,"files/MCSelected.root");
    }
  }
  
  cout<<"opening file "<<OutNameEvent<<" for writing output"<<endl;
  fEvent.open(OutNameEvent,ios::out);
  if(IsData){
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	fEvent<<DataSelected[e0][e1]<<" ";	
      }
    }
  }
  else{
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	fEvent<<DataSelected[e0][e1]<<" ";	
      }
    }
    fEvent<<666<<" ";
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	    fEvent<<MCSelected[c0][c1][e0][e1]<<" ";
	  }
	}
      }
    }
    fEvent<<666<<" ";
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	fEvent<<BkgSelected[e0][e1]<<" ";
      }
    }
    fEvent<<666<<" ";
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	fEvent<<Efficiency[c0][c1]<<" ";
      }
    }
  }

  fEvent<<666;
  if(IsData) NeutrinoFlux->Scale(POTCount/1e21);
  else NeutrinoFlux->Scale((DataEquivalent*1e21)/1e21);

  fEvent<<" "<<NeutrinoFlux->Integral()<<" ";
  fEvent<<666;
  
  fEvent.close();
  cout<<"file "<<OutNameEvent<<" is closed"<<endl;

  if(!IsData){
    TFile * wfile = new TFile(OutNameEvent_root,"recreate");
    MCTrueEvents->Write();
    wfile->Close();
  }
  
  delete MCTrueEvents;
  
  if(Plots){
    if(!IsData){
      THStack * Stack_MuCL = new THStack("Stack_MuCL","");
      THStack * Stack_MuCL_Lowest = new THStack("Stack_MuCL_Lowest","");
      THStack * Stack_NTracks = new THStack("Stack_NTracks","");
      THStack * Stack_RecMom = new THStack("Stack_RecMom","");
      THStack * Stack_RecAngle = new THStack("Stack_RecAngle","");

      THStack * Stack_RecMom_CC0pi = new THStack("Stack_RecMom_CC0pi","");
      THStack * Stack_RecAngle_CC0pi = new THStack("Stack_RecAngle_CC0pi","");
      THStack * Stack_RecMom_CC1pi = new THStack("Stack_RecMom_CC1pi","");
      THStack * Stack_RecAngle_CC1pi = new THStack("Stack_RecAngle_CC1pi","");
      THStack * Stack_RecMom_CCNpi = new THStack("Stack_RecMom_CCNpi","");
      THStack * Stack_RecAngle_CCNpi = new THStack("Stack_RecAngle_CCNpi","");

      ProduceStack(hMuCL,Stack_MuCL);
      ProduceStack(hMuCL_Lowest,Stack_MuCL_Lowest);
      //No stack produced for 2tracks 2D plot, but should still normalise as done in the produce stack!
      for(int fsi=0;fsi<NFSIs;fsi++) hMuCL_2tracks[fsi]->Scale((5.86e20/1e21)*100/1000);
      
      ProduceStack(hNTracks,Stack_NTracks);
      ProduceStack(hRecMom,Stack_RecMom);
      ProduceStack(hRecAngle,Stack_RecAngle);

      ProduceStack(hRecMom_CC0pi,Stack_RecMom_CC0pi);
      ProduceStack(hRecAngle_CC0pi,Stack_RecAngle_CC0pi);
      ProduceStack(hRecMom_CC1pi,Stack_RecMom_CC1pi);
      ProduceStack(hRecAngle_CC1pi,Stack_RecAngle_CC1pi);
      ProduceStack(hRecMom_CCNpi,Stack_RecMom_CCNpi);
      ProduceStack(hRecAngle_CCNpi,Stack_RecAngle_CCNpi);
      /*
	TLegend * leg = new TLegend(0.7,0.8,1,1);
	if(hRecMom[0]) leg->Add(hRecMom[0],"CC0#pi");
	if(hRecMom[1]) leg->Add(hRecMom[1],"CC0#pi");
      
      if(fsi<3){//CC0Pi+/-/0
	if(fsi==0) h[0]->SetFillColor(kRed);
	else if(fsi==1) h[1]->SetFillColor(kOrange+10);
	else h[2]->SetFillColor(kOrange+9);
      }
      else{
	h[3]->SetFillColor(kBlue+2);//CC1Pi+/-
	h[4]->SetFillColor(kAzure+7);//CCNPi+/-
	h[5]->SetFillColor(kAzure+10);//CCpi0
	h[6]->SetFillColor(kGray);//NC
	if(fsi>6 && fsi<11){
	  if(fsi==7) h[7]->SetFillColor(kYellow);//Sand
	  if(fsi==8) h[8]->SetFillColor(kYellow+4);//Horizontal||Vertical Module Bkg
	  if(fsi==9) h[9]->SetFillColor(kMagenta);//AntiNuMu
	  if(fsi==10) h[10]->SetFillColor(kGreen);//NuE
	}
      }
      */      
      
      TFile * fMC = new TFile("plots/MCPlots_Plan.root","RECREATE");
      Stack_MuCL->Write();
      Stack_MuCL_Lowest->Write();
      Stack_NTracks->Write();
      Stack_RecMom->Write();
      Stack_RecAngle->Write();
      Stack_RecMom_CC0pi->Write();
      Stack_RecAngle_CC0pi->Write();
      Stack_RecMom_CC1pi->Write();
      Stack_RecAngle_CC1pi->Write();
      Stack_RecMom_CCNpi->Write();
      Stack_RecAngle_CCNpi->Write();
      for(int fsi=0;fsi<NFSIs;fsi++){
	hMuCL[fsi]->Write();
	hMuCL_Lowest[fsi]->Write();
	hMuCL_2tracks[fsi]->Write();
	hNTracks[fsi]->Write();
      }
      hMuCL_TrueMuon->Write();
      hMuCL_TruePion->Write();
      hMuCL_TrueProton->Write();
      PE_Lowest_CC0pi->Write();
      PE_Lowest_Other->Write();
      
      fMC->Close();
    }
    else{
      TFile * fData = new TFile("plots/DataPlots.root","RECREATE");
      hMuCL[0]->Write("Data_MuCL");
      hMuCL_Lowest[0]->Write("Data_MuCL_Lowest");
      hNTracks[0]->Write("Data_NTracks");
      hRecMom[0]->Write("Data_RecMom");
      hRecAngle[0]->Write("Data_RecAngle");
      hRecMom_CC0pi[0]->Write("Data_RecMom");
      hRecAngle_CC0pi[0]->Write("Data_RecAngle");
      hRecMom_CC1pi[0]->Write("Data_RecMom");
      hRecAngle_CC1pi[0]->Write("Data_RecAngle");
      hRecMom_CCNpi[0]->Write("Data_RecMom");
      hRecAngle_CCNpi[0]->Write("Data_RecAngle");

      hMuCL_TrueMuon->Write();
      hMuCL_TruePion->Write();
      hMuCL_TrueProton->Write();

      fData->Close();
    }
    //ProduceStack(hNTracks,Stack_NTracks);
    //ProduceStack(hRecMom,Stack_RecMom);
    //ProduceStack(hRecAngle,Stack_RecAngle);
    for(int fsi=0;fsi<NFSIs;fsi++){
      delete hMuCL[fsi];
      delete hMuCL_2tracks[fsi];
      delete hMuCL_Lowest[fsi];
      delete hNTracks[fsi];
      delete hRecMom[fsi];
      delete hRecAngle[fsi];
      delete hRecMom_CC0pi[fsi];
      delete hRecAngle_CC0pi[fsi];
      delete hRecMom_CC1pi[fsi];
      delete hRecAngle_CC1pi[fsi];
      delete hRecMom_CCNpi[fsi];
      delete hRecAngle_CCNpi[fsi];
    
    }
    delete hMuCL_TrueMuon;
    delete hMuCL_TruePion;
    delete hMuCL_TrueProton;
    delete PE_Lowest_CC0pi;
    delete PE_Lowest_Other;  
    }
  delete rand;
}

int main(int argc, char ** argv){

  cout<<"welcome, here is the CC0pi selection"<<endl;
  char fName[256];
  bool MC=false;
  bool Data=true;
  int Sample=0;
  bool IsPlots=true;
  bool Systematics_Flux=false;
  bool Systematics_Xsec=false;
  bool Systematics_Detector=false;
  int ErrorType=0;int n=0;
  int Xsec_dial=0;
  int File_Number=0;
  TVectorD * FluxVector=new TVectorD();
  
  int c=-1;
  while ((c = getopt(argc, argv, "w:ms:fxd")) != -1) {
    switch(c){
    case 'w':
      SandReweight=atof(optarg);
      break;
    case 'm':
      MC=true;
      Data=false;
      break;
    case 's':
      Sample=atoi(optarg);//0 for sand muon, 1 for CC0pi
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
  
  cout<<"welcome"<<endl;
  //  XS->Xsec::Initialize();
  InitializeGlobal();
  
  TChain * chain = new TChain("wtree");
   if(Data){
     for(int i=StartRun;i<=EndRun;i++){
       //if(i>=14337 && i<=14429) continue;
       sprintf(fName,"root_input/CC0piTree_%d.root",i);
       cout<<fName<<endl;
       chain->Add(fName);
     }
   }
   else{
     cout<<"Please enter the amount of data you'd like to mimic (I will adjust MC stat.) in unit of 1e21 POT:"<<endl;
     cin>>DataEquivalent;
     ScalingMC=DataEquivalent/NMCfiles;//one MC file is equivalent to 1e21 POT
       
          for(int i=1;i<=NMCfiles;i++){
       if((i==165 || i==166 || i==167 || i==168) || (i==665 || i==666 || i==667 || i==668) || i==555) continue;
       //sprintf(fName,"root_input/CC0piTree%d.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Likelihood_ReWeight.root",i);
       sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Old_%d_Likelihood.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Likelihood.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_RandomPE.root",i);
       cout<<fName<<endl;
       chain->Add(fName);
       }
     
   }

   
   if(Systematics_Flux && !Data){
     TFile * FluxError = new TFile("files/ErrorVarFinal.root");
     FluxVector = new TVectorD();
     for(int iflux=0;iflux<=NFluxFiles;iflux++){
       cout<<"variation in flux #"<<iflux<<endl;
       sprintf(Name,"Var[%d]",iflux);
       FluxVector = (TVectorD*) FluxError->Get(Name);
       File_Number=iflux;
       CC0piDistributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,n);
     }
     FluxError->Close();
   }
   if(Systematics_Xsec && !Data){
     cout<<"Center="<<CenterXsecVariations<<endl;
     for(int ixsec=StartXsec;ixsec<=EndXsec;ixsec++){
       for(int sig=CenterXsecVariations-1;sig<=CenterXsecVariations+1;sig++){//-1,0,1 sigma
	 Xsec_dial=ixsec*NXsecVariations+sig;
	 File_Number=Xsec_dial;
 	 CC0piDistributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,n);	 
       }
     }
   }
   if(Systematics_Detector){
     for(ErrorType=StartError;ErrorType<=EndError;ErrorType++){
       SandReweight=Start[ErrorType]+Step[ErrorType];//Set SandReweight to the nominal value
       cout<<endl<<"Error Type="<<ErrorType<<endl;

       for(n=0;n<NE[ErrorType];n++){
	 double ErrorValue=Start[ErrorType]+n*Step[ErrorType];
	 cout<<"Error Value="<<ErrorValue<<", files generated:";
	 
	 if(Data){
	   chain->Delete();
	   for(int i=StartRun;i<=EndRun;i++){  
	     for(int j=StartSubRun;j<=EndSubRun;j++){  
	       sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_%08d_%04d_Systematics%d_%d.root",i,j,ErrorType,n);
	     }
	   }
	 }
	 else{
	   chain->Delete();
	   if(ErrorType==6){//Change sand weight
	     SandReweight=ErrorValue;
	     for(int i=1;i<=NMCfiles;i++){
	       if((i==165 || i==166 || i==167 || i==168) || (i==665 || i==666 || i==667 || i==668) || i==555) continue;
	       sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_ReWeight.root",i);
	       cout<<fName<<endl;
	       chain->Add(fName);
	     }
	   }
	   else{
	     for(int i=1;i<=NMCfiles;i++){	
	       if((i==165 || i==166 || i==167 || i==168) || (i==665 || i==666 || i==667 || i==668) || i==555) continue;
	       sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Systematics%d_%d_Plan.root",i,ErrorType,n);
	       cout<<fName<<endl;
	       chain->Add(fName);
	     }
	   }
	 }
	 CC0piDistributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,n);	 
       }
     }

     
   }
   else{
     cout<<"Simple selection on the MC or data"<<endl;
     CC0piDistributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,n);
   }
   
   return 0;
}

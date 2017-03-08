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
#include <TPie.h>
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
#define DEBUG

//Xsec * XS = new Xsec();
char Type[32];
char Name[256];
char Name0[256];
double ScalingMC;
double DataEquivalent=1;
double SandReweight=0.6;
double MuonSample1=3;
double MuonSample2=3;
double MuonCut=-2.5;
double ProtonCut=-3.5;
double MuonCut_WM=0.3;// to be tuned
double ProtonCut_WM=0.3;// to be tuned


void ProduceStack(TH1D * h[NFSIs], THStack * hStack){

  for(int fsi=0;fsi<NFSIs;fsi++){
    
    h[fsi]->GetYaxis()->SetTitleOffset(1.3);
    //h[fsi]->Scale((5.86e20/1e21)*100/1000);
    
    if(fsi<3){//CC0Pi+/-/0
      if(fsi==0) h[0]->SetFillColor(kOrange+10);
      else if(fsi==1) h[1]->SetFillColor(kRed);
      else h[2]->SetFillColor(kOrange+8);
      h[2]->SetLineWidth(3);
    }
    else{
      h[3]->SetFillColor(kAzure+10);//CC1Pi+/-
      h[4]->SetFillColor(kAzure+7);//CCNPi+/-
      h[5]->SetFillColor(kBlue+2);//CCpi0
      h[6]->SetFillColor(kGray);//NC
      if(fsi>6 && fsi<=11){
	if(fsi==7) h[7]->SetFillColor(kYellow);//Sand
	if(fsi==8) h[8]->SetFillColor(kYellow+4);//Horizontal Module Bkg
	if(fsi==9) h[9]->SetFillColor(kYellow+7);//Vertical Module Bkg
	if(fsi==10) h[10]->SetFillColor(kMagenta);//AntiNuMu
	if(fsi==11) h[11]->SetFillColor(kGreen);//NuE
	}
    }
    h[fsi]->SetLineColor(1);
    h[fsi]->SetLineWidth(2);
    hStack->Add(h[fsi]);
  }
  h[0]->SetTitle("CC0#pi (0p)");
  h[1]->SetTitle("CC0#pi (1p)");
  h[2]->SetTitle("MEC");
  h[3]->SetTitle("CC1#pi^{#pm}");
  h[4]->SetTitle("CC1#pi^{0}");
  h[5]->SetTitle("CCother");
  h[6]->SetTitle("NC");
  h[7]->SetTitle("Sand");
  h[8]->SetTitle("INGRID Bkg H");
  h[9]->SetTitle("INGRID Bkg V");
  h[10]->SetTitle("anti-#nu_{#mu}");
  h[11]->SetTitle("#nu_{e}");
}

  
void InitialiseTable(double DataSelected[NBinsRecMom][NBinsRecAngle],double MCSelected[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle],double BkgSelected[NBinsRecMom][NBinsRecAngle],double Efficiency[NBinsTrueMom][NBinsTrueAngle],double TotalCC0piEvent[NBinsTrueMom][NBinsTrueAngle],double TotalCC1piEvent[NBinsTrueMom][NBinsTrueAngle]){
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      Efficiency[c0][c1]=0;
      TotalCC0piEvent[c0][c1]=0;
      TotalCC1piEvent[c0][c1]=0;
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

//void CC0piDistributions(TChain * wtree,bool IsData,int Sample,bool IsPlots,int SelectedError_Source,double SelectedError_Variation,char * OutNameEvent);

void Distributions(TChain * wtree,bool IsData,int Selection,bool Plots,bool Systematics_Flux,int File_Number,TVectorD FluxVector,bool Systematics_Xsec,int dial,bool Systematics_Detector,int ErrorType,char * outnameevent, bool _isPM){
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
  int FSIInt;//Final state true information
  int Num_Int;//Interaction number
  int nTracks[10];//number of tracks reconstructed per vertex
  double weight;//weight of each event
  bool IsFV;//True vertex is in FV
  bool IsSand;
  bool IsAnti;
  bool IsBkgH;
  bool IsBkgV;
  bool IsNuE;
  bool IsDetected[10];//Event is reconstructed (means >=1 vertex reconstructed)
  bool SelectionFV[10];bool SelectionOV[10];//Vertex is reconstructed within/out of FV
  double Enu;//True neutrino energy
  int nIngBasRec;//Number of vertices reconstructed per event
  double Errorweight;
  double TrueMomentumMuon, TrueAngleMuon;//True muon kinematic variables
  double TrueMomentumPion, TrueAnglePion;//True pion kinematic variables
  double POT;//Number of POT
  double TrackAngle[10][20];//Reconstructed angle of the reconstructed track
  double TrackThetaX[10][20];
  double TrackThetaY[10][20];
  double CoplanarityAngle[10];
  double OpeningAngle[10];
  int TypeOfTrack[10][20];//pdg of the reconstructed track. Careful that it is based on an algorithm defined in Reconstruction.cc
  double CLMuon_KS[10][20];// = new vector<double> [10][20];
  double CLMuon_Likelihood[10][20];
  double CLMuon_Plan[10][20];
  double ProportionHighPE[10][20];
  double MeanHighPE[10][20];
  double HighestPE[10][20];
  double Momentum[10][20];
  double IronDistance[10][20];
  double PlasticDistance[10][20];
  int Sample[10][20];//Geometric properties of the track, defined in Reconstruction::SelectTrackSample
  bool IsReconstructed[10][20];
  double ReWeight[175];
  double CriteriaAngleX[10][20];
  double CriteriaAngleY[10][20];
  double CriteriaHalfWayX[10][20];
  double CriteriaHalfWayY[10][20];
  double TrackWidth[10][20];
  double TotalCharge[10][20];
  int Spill;
  int GoodSpill;
  int LastChannelIX[10][20],LastChannelIY[10][20];

  /////////////////////////////////PREPARE FOR READING/////////////////////////////////////////
  wtree->SetBranchAddress("Spill",&Spill);
  wtree->SetBranchAddress("GoodSpill",&GoodSpill);
  wtree->SetBranchAddress("FSIInt",&FSIInt);
  wtree->SetBranchAddress("nIngBasRec",&nIngBasRec);
  wtree->SetBranchAddress("InteractionType",&Num_Int);
  wtree->SetBranchAddress("nTracks[10]",nTracks);
  wtree->SetBranchAddress("weight",&weight);
  wtree->SetBranchAddress("Enu",&Enu);
  wtree->SetBranchAddress("TrueMomentumMuon",&TrueMomentumMuon);
  wtree->SetBranchAddress("TrueAngleMuon",&TrueAngleMuon);
  wtree->SetBranchAddress("TrueMomentumPion",&TrueMomentumPion);
  wtree->SetBranchAddress("TrueAnglePion",&TrueAnglePion);
  wtree->SetBranchAddress("IsFV",&IsFV);
  wtree->SetBranchAddress("IsSand",&IsSand);
  wtree->SetBranchAddress("IsAnti",&IsAnti);
  wtree->SetBranchAddress("IsBkgH",&IsBkgH);
  wtree->SetBranchAddress("IsBkgV",&IsBkgV);
  wtree->SetBranchAddress("IsNuE",&IsNuE);
  wtree->SetBranchAddress("SelectionFV[10]",SelectionFV);
  wtree->SetBranchAddress("SelectionOV[10]",SelectionOV);
  wtree->SetBranchAddress("IsDetected[10]",IsDetected);
  wtree->SetBranchAddress("TrackAngle[10][20]",TrackAngle);
  wtree->SetBranchAddress("OpeningAngle[10]",OpeningAngle);
  wtree->SetBranchAddress("CoplanarityAngle[10]",CoplanarityAngle);
  wtree->SetBranchAddress("TrackThetaX[10][20]",TrackThetaX);
  wtree->SetBranchAddress("TrackThetaY[10][20]",TrackThetaY);
  wtree->SetBranchAddress("TypeOfTrack[10][20]",TypeOfTrack);
  wtree->SetBranchAddress("CLMuon_KS[10][20]",CLMuon_KS);
  wtree->SetBranchAddress("CLMuon_Plan[10][20]",CLMuon_Plan);
  wtree->SetBranchAddress("CLMuon_Likelihood[10][20]",CLMuon_Likelihood);
  wtree->SetBranchAddress("ProportionHighPE[10][20]",ProportionHighPE);
  wtree->SetBranchAddress("MeanHighPE[10][20]",MeanHighPE);
  wtree->SetBranchAddress("HighestPE[10][20]",HighestPE); 
  wtree->SetBranchAddress("Momentum[10][20]",Momentum);
  wtree->SetBranchAddress("IronDistance[10][20]",IronDistance);
  wtree->SetBranchAddress("PlasticDistance[10][20]",PlasticDistance);
  wtree->SetBranchAddress("Sample[10][20]",Sample);
  wtree->SetBranchAddress("IsReconstructed[10][20]",IsReconstructed);
  wtree->SetBranchAddress("CriteriaAngleX[10][20]",CriteriaAngleX);
  wtree->SetBranchAddress("CriteriaAngleY[10][20]",CriteriaAngleY);
  wtree->SetBranchAddress("CriteriaHalfWayX[10][20]",CriteriaHalfWayX);
  wtree->SetBranchAddress("CriteriaHalfWayY[10][20]",CriteriaHalfWayY);
  wtree->SetBranchAddress("ReWeight[175]",ReWeight);
  wtree->SetBranchAddress("POT",&POT);
  wtree->SetBranchAddress("TrackWidth[10][20]",TrackWidth);
  wtree->SetBranchAddress("TotalCharge[10][20]",TotalCharge);
  wtree->SetBranchAddress("LastChannelINGRIDX[10][20]",LastChannelIX);
  wtree->SetBranchAddress("LastChannelINGRIDY[10][20]",LastChannelIY);


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
  double DataSelected_full[NBinsRecMom][NBinsRecAngle];
  double MCSelected_full[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double BkgSelected_full[NBinsRecMom][NBinsRecAngle];
  double Efficiency_full[NBinsTrueMom][NBinsTrueAngle];
  double TotalCC0piEvent[NBinsTrueMom][NBinsTrueAngle]={{0}};
  double TotalCC1piEvent[NBinsTrueMom][NBinsTrueAngle]={{0}};
  double TotalCC1pi=0;

  const int nCuts=7;
  const string myCuts[nCuts]={"Generated","Reconstruted + one INGRID track","FCFV","CC1pi Selection","2 or 3 tracks","MuTrk is INGRID stop/through","MuTrk is INGRID stop"};
  double nEventsInter[nCuts][NFSIs]={{0}};
  double nEvents[nCuts]={0};
  TH1D* CutEfficiency=new TH1D("Efficiency","Efficiency",nCuts,0,nCuts);
  TH1D* CutPurity=new TH1D("Purity","Purity",nCuts,0,nCuts);
  for(int i=0;i<nCuts;i++){
    CutEfficiency->GetXaxis()->SetBinLabel(i+1,myCuts[i].c_str());
    CutPurity->GetXaxis()->SetBinLabel(i+1,myCuts[i].c_str());
  }
  CutEfficiency->SetMarkerColor(kBlue);
  CutPurity->SetMarkerColor(kRed);
  CutEfficiency->SetLineColor(kBlue);
  CutPurity->SetLineColor(kRed);
  CutEfficiency->SetMarkerStyle(2);
  CutPurity->SetMarkerStyle(2);

  InitialiseTable(DataSelected,MCSelected,BkgSelected,Efficiency,TotalCC0piEvent,TotalCC0piEvent);
  InitialiseTable(DataSelected_full,MCSelected_full,BkgSelected_full,Efficiency_full,TotalCC0piEvent,TotalCC0piEvent);

  TH2D * MCTrueEvents = new TH2D("MCTrueEvents","",NBinsTrueMom,BinningTrueMom,NBinsTrueAngle,BinningTrueAngle);
  
  TH1D * hMuCL[NFSIs];
  TH1D * hMuCL_Lowest[NFSIs];
  TH2D * hMuCL_2tracks[NFSIs];
  TH1D * hNTracks[NFSIs];
  TH1D * hRecMom[NFSIs];
  TH1D * hRecAngle[NFSIs];
  TH1D * hFCFVTrueEvents[NFSIs];

  TH1D * hRecMom_CC0pi[NFSIs];
  TH1D * hRecAngle_CC0pi[NFSIs];

  TH1D * hRecMom_CC1pi[NFSIs];
  TH1D * hRecAngle_CC1pi[NFSIs];
  TH1D * hRecMom_CC1pi_restr[NFSIs];
  TH1D * hRecAngle_CC1pi_restr[NFSIs];
  TH1D * hRecMom_CC1pi_full[NFSIs];
  TH1D * hRecAngle_CC1pi_full[NFSIs];

  TH1D * hRecMom_CCNpi[NFSIs];
  TH1D * hRecAngle_CCNpi[NFSIs];

  TH1D * hSampleSecondTrack[NFSIs];

  TH1D * hMuCL_TrueMuon;
  TH1D * hMuCL_TruePion;
  TH1D * hMuCL_TrueProton;
  TH1D * hMuCL_TrueOthers;

  TH2D * PE_Lowest_CC0pi;
  TH2D * PE_Lowest_Other;
  
  TH2D * MCEfficiency;
  TH1D * MCEfficiency_Energy;
  TH1D * TotalCC0piEvent_Energy;
  TH1D * TotalCC1piEvent_Energy;

  TPie * MuonRec_TruePDG;
  TPie * MuonRec_TruePDG_switch;
  double MuonRec_TruePDG_val[4];  
  double MuonRec_TruePDG_switch_val[4];  

  TH2D* Pmu_vs_IronDist, *Pp_vs_IronDist, *Ppi_vs_IronDist;
  TH2D* MuCL_vs_IronDist;

  if(Plots){
    if(_isPM){
      hMuCL_TrueMuon = new TH1D("hMuCL_TrueMuon","Muon confidence level for the sand muons",100,-50,0);//ML tmp 500bins->100bins
      hMuCL_TruePion = new TH1D("hMuCL_TruePion","Muon confidence level for the true muons",100,-50,0);//ML tmp 500bins->100bins
      hMuCL_TrueProton = new TH1D("hMuCL_TrueProton","Muon confidence level for the true protons",100,-50,0);//ML tmp 500bins->100bins
      hMuCL_TrueOthers = new TH1D("hMuCL_TrueOthers","Muon confidence level for the true protons",100,-50,0);//ML tmp 500bins->100bins
      MuCL_vs_IronDist=new TH2D("MuCL_vs_IronDist","highest #mu_{CL} track",50,0,100,20,MuonCut,0);
    }
    else {
      hMuCL_TrueMuon = new TH1D("hMuCL_TrueMuon","Muon confidence level for the sand muons",50,0,1);
      hMuCL_TruePion = new TH1D("hMuCL_TruePion","Muon confidence level for the true muons",50,0,1);
      hMuCL_TrueProton = new TH1D("hMuCL_TrueProton","Muon confidence level for the true protons",50,0,1);
      hMuCL_TrueOthers = new TH1D("hMuCL_TrueOthers","Muon confidence level for the true protons",50,0,1);
      MuCL_vs_IronDist=new TH2D("MuCL_vs_IronDist","highest #mu_{CL} track",50,0,100,20,MuonCut,1);
    }    

    MuCL_vs_IronDist->GetXaxis()->SetTitle("Equivalent iron length (cm)");
    MuCL_vs_IronDist->GetYaxis()->SetTitle("#mu_{CL}");

    Pmu_vs_IronDist=new TH2D("Pmu_vs_IronDist","",50,0,100,50,0,6);
    Pmu_vs_IronDist->GetXaxis()->SetTitle("Equivalent iron length");
    Pmu_vs_IronDist->GetYaxis()->SetTitle("Muon Momentum");

    Pp_vs_IronDist=new TH2D("Pp_vs_IronDist","",50,0,100,50,0,6);
    Pp_vs_IronDist->GetXaxis()->SetTitle("Equivalent iron length");
    Pp_vs_IronDist->GetYaxis()->SetTitle("Proton Momentum");

    Ppi_vs_IronDist=new TH2D("Ppi_vs_IronDist","",50,0,100,50,0,6);
    Ppi_vs_IronDist->GetXaxis()->SetTitle("Equivalent iron length");
    Ppi_vs_IronDist->GetYaxis()->SetTitle("Pion Momentum");

    PE_Lowest_CC0pi = new TH2D("PE_Lowest_CC0pi","",50,0,1,200,0,1000);
    PE_Lowest_Other = new TH2D("PE_Lowest_Other","",50,0,1,200,0,1000);
    
    MCEfficiency = new TH2D("MCEfficiency","",NBinsTrueMom,0,NBinsTrueMom,NBinsTrueAngle,0,NBinsTrueAngle);
    MCEfficiency_Energy = new TH1D("MCEfficiency_Energy","",100,0,10);
    TotalCC0piEvent_Energy = new TH1D("TotalCC0piEvent_Energy","",100,0,10);
    TotalCC1piEvent_Energy = new TH1D("TotalCC1piEvent_Energy","",100,0,10);
    MCEfficiency_Energy->Sumw2();TotalCC0piEvent_Energy->Sumw2();TotalCC1piEvent_Energy->Sumw2();

    
      char Name[256];char Title[256];

      for(int i=0;i<NFSIs;i++){
	sprintf(Name,"hFCFVTrueEvents%d",i);
	sprintf(Title,"",i);
	hFCFVTrueEvents[i] = new TH1D(Name,Title,100,0,10);
	//hMuCL->Sumw2();
	hFCFVTrueEvents[i]->GetXaxis()->SetTitle("E_{#nu} (GeV)");    
	hFCFVTrueEvents[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hMuCL%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi",i);
	hMuCL[i] = new TH1D(Name,Title,500,-50,0);
	//hMuCL->Sumw2();
	hMuCL[i]->GetXaxis()->SetTitle("#mu_{CL}");    hMuCL[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hMuCL_Lowest%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi, only for the lowest MuCL tracks if ntracks>1",i);
	hMuCL_Lowest[i] = new TH1D(Name,Title,500,-50,0);
	//hMuCL_Lowest->Sumw2();
	hMuCL_Lowest[i]->GetXaxis()->SetTitle("#mu_{CL_Lowest}");    hMuCL_Lowest[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hMuCL_2tracks%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi",i);
	hMuCL_2tracks[i] = new TH2D(Name,Title,500,-50,0,500,-50,0);
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
	hRecMom[i] = new TH1D(Name,Title,20,0,100);
	//hRecMom->Sumw2();
	hRecMom[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
	hRecMom[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle->Sumw2();
	hRecAngle[i]->GetXaxis()->SetTitle("Angle (°)");
	hRecAngle[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecMom_CC0pi%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CC0pi[i] = new TH1D(Name,Title,20,0,100);
	//hRecMom_CC0pi->Sumw2();
	hRecMom_CC0pi[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
	hRecMom_CC0pi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CC0pi%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CC0pi[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle_CC0pi->Sumw2();
	hRecAngle_CC0pi[i]->GetXaxis()->SetTitle("Angle_CC0pi (°)");
	hRecAngle_CC0pi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecMom_CC1pi%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CC1pi[i] = new TH1D(Name,Title,20,0,100);
	//hRecMom_CC1pi->Sumw2();
	hRecMom_CC1pi[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
	hRecMom_CC1pi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CC1pi%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CC1pi[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle_CC1pi->Sumw2();
	hRecAngle_CC1pi[i]->GetXaxis()->SetTitle("Angle_CC1pi (°)");
	hRecAngle_CC1pi[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hRecMom_CC1pi_restr%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CC1pi_restr[i] = new TH1D(Name,Title,14,0,100);
	//hRecMom_CC1pi_restr->Sumw2();
	hRecMom_CC1pi_restr[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
	hRecMom_CC1pi_restr[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CC1pi_restr%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CC1pi_restr[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle_CC1pi_restr->Sumw2();
	hRecAngle_CC1pi_restr[i]->GetXaxis()->SetTitle("Angle_CC1pi (°)");
	hRecAngle_CC1pi_restr[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hRecMom_CC1pi_full%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CC1pi_full[i] = new TH1D(Name,Title,20,0,100);
	//hRecMom_CC1pi_full->Sumw2();
	hRecMom_CC1pi_full[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
	hRecMom_CC1pi_full[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CC1pi_full%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CC1pi_full[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle_CC1pi_full->Sumw2();
	hRecAngle_CC1pi_full[i]->GetXaxis()->SetTitle("Angle_CC1pi (°)");
	hRecAngle_CC1pi_full[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hRecMom_CCNpi%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom_CCNpi[i] = new TH1D(Name,Title,20,0,100);
	//hRecMom_CCNpi->Sumw2();
	hRecMom_CCNpi[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
	hRecMom_CCNpi[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle_CCNpi%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle_CCNpi[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle_CCNpi->Sumw2();
	hRecAngle_CCNpi[i]->GetXaxis()->SetTitle("Angle_CCNpi (°)");
	hRecAngle_CCNpi[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hSampleSecondTrack%d",i);
	sprintf(Title,"Number of tracks for the %dth-fsi",i);
	hSampleSecondTrack[i] = new TH1D(Name,Title,6,0,6);
	//hSampleSecondTrack->Sumw2();
	hSampleSecondTrack[i]->GetXaxis()->SetTitle("Track sample");
	hSampleSecondTrack[i]->GetYaxis()->SetTitle("Number of events");

      }

  }    
  double NEvents=0;double NEventsLost=0;double NEventsLostDetected=0;
  
  //cout<<"going to read events"<<endl;
  for(int ievt=0;ievt<nevt;ievt++){
    POTCount+=POT;
    //if(ievt%10000==0){cout<<ievt<<", "<<"POT processed="<<POTCount<<", xsec="<<IsXsec<<endl;
    if(ievt%100000==0)cout<<"Number of events lost due to detection (%)="<<100.*NEventsLostDetected/NEvents<<", due to detection & muon not found="<<100.*NEventsLost/NEvents<<endl;
    if(IsData && ievt==(nevt-1)) cout<<"Final POT="<<POTCount<<endl;
    wtree->GetEntry(ievt);
    if(IsData){FSIInt=0;Num_Int=0;}
    else weight*=ScalingMC;//Adjust the MC distribution to the amount of data we process
    if(IsSand) weight*=(1+SandReweight);
    ////////////////DETERMINE THE REWEIGHTING OF THE EVENT IF SYSTEMATICS ERROR FLUX
    if(Systematics_Flux){
      for(int i=0;i<NBinsEnergyFlux;i++){
	if(Enu<BinningEnergyFlux[i+1]) {
	  Errorweight=FluxVector[i];
	  break;
	}
      }
      weight*=(1+Errorweight);
    }
    else if(Systematics_Xsec){
      if(ievt%10000==0)cout<<"weight="<<weight<<", and after reweight="<<weight*ReWeight[dial]<<endl;
      weight=weight*ReWeight[dial];
    }
    //cout<<ReWeight[5]<<endl;

    //********************
    if(IsSand||IsBkgH || IsBkgV){
      //if(weight>1000) continue;
      continue;
    }
    //******************** tmp
    
    bool IsNuMu=((!IsSand) && (!IsAnti) && (!IsNuE) && (!IsBkgH) && (!IsBkgV));
    if(IsData && (!GoodSpill || !Spill)) continue;
    if((/*IsFV &&*/ !IsData) || IsData){ 
      // why that condition here ?? too easy to remove all external background this way

#ifdef DEBUG
      //cout<<"FV"<<", nbasrec="<<nIngBasRec<<endl;
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
      if(!IsData && IsFV && IsNuMu ){
	//efficiency is 100 <-> all FCFV true events
	hFCFVTrueEvents[FSIInt]->Fill(Enu,weight);

	if(FSIInt<3 && Selection==1){
	  TotalCC0piEvent[BinTrueMom][BinTrueAngle]+=weight;
	  MCTrueEvents->Fill(BinTrueMom+1,BinTrueAngle+1,weight);
	  if(Plots) TotalCC0piEvent_Energy->Fill(Enu,weight);
	}
	else if(FSIInt==3 && Selection==2){
	  TotalCC1pi+=weight;
	  TotalCC1piEvent[BinTrueMom][BinTrueAngle]+=weight;
	  if(Plots) TotalCC1piEvent_Energy->Fill(Enu,weight);
	}
      }

      nEvents[0]+=weight;
      nEventsInter[0][FSIInt]+=weight;

      for(int irec=0;irec<nIngBasRec;irec++){

	bool MuonFound=false;
	int MuonTrue;int MuonRec=0;
	int PionRec=0;
	int LowestMuCL=0;
	bool Trash=false;
	//if(IsData) FSIInt=0;
	//if(Num_Int==1) FSIInt=1;
	//else if(Num_Int==2) FSIInt=2;
	//else if(Num_Int>2 && Num_Int<20) FSIInt=3;
	//else if(Num_Int<30) FSIInt=4;
	//else if(Num_Int>=30) FSIInt=5;
	NEvents+=weight;
	if(!IsDetected[irec]) {NEventsLostDetected+=weight;continue;}
	nEvents[1]+=weight;
	nEventsInter[1][FSIInt]+=weight;

	if(Selection>=1 && !SelectionFV[irec]) continue; // 2 is CC1pi
	if(Selection==0 && !SelectionOV[irec]) continue;
	nEvents[2]+=weight;
	nEventsInter[2][FSIInt]+=weight;

	bool AllTracksWellReconstructed=true;
	int nBadRecTracks=0;

	// we use different variable for PM and WM muon-likelihood
	double *_mucl=(_isPM? CLMuon_Plan[irec] : CLMuon_KS[irec]);



	// ML tmp
	/*
	if(_isPM) {
	  _mucl=CLMuon_Likelihood[irec];
	  ProtonCut=0.3;
	  MuonCut=0.4;
	}
	*/
	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  if(!IsReconstructed[irec][itrk]) {
	    AllTracksWellReconstructed=false;
	    nBadRecTracks++;
	  }
	}



	//if(!AllTracksWellReconstructed){cout<<"one track not well reconstructed"<<endl;continue;}
	/*	
#ifdef DEBUG
	cout<<"Number of tracks="<<nTracks[irec]<<endl;
	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  cout<<"Sample="<<Sample[irec][itrk]<<", MuCL="<<TMath::Log(CLMuon_Plan[irec][itrk])<<", iron distance="<<IronDistance[irec][itrk]<<endl;
	}
	#endif*/
	////////////////////DETERMINE MUON TRACK/////////////////////////////////
	  for(int itrk=0;itrk<nTracks[irec];itrk++){
	    //if(CLMuon_Plan[irec][itrk]==0) CLMuon_Plan[irec][itrk]=1e-30;
	    if(IsReconstructed[irec][itrk]){
	      //cout<<CLMuon_Plan[irec][itrk]<<endl;


	      if(_isPM) _mucl[itrk]=TMath::Log(_mucl[itrk]);
  	      //CLMuon_Plan[irec][itrk]=TMath::Log(CLMuon_Plan[irec][itrk]);

	      //if(CLMuon_Plan[irec][itrk]<-50) CLMuon_Plan[irec][itrk]=-50;
	      //cout<<TypeOfTrack[irec][itrk]<<endl;
	      if(_mucl[itrk]!=_mucl[itrk] || (_isPM && _mucl[itrk]==-1)){ Trash=true; cout<<"problem, event trashed"<<endl;continue;}
	      if(TypeOfTrack[irec][itrk]==13) MuonTrue=itrk;
	      if(_mucl[itrk]>=_mucl[MuonRec]) {PionRec=MuonRec; MuonRec=itrk;}
	      else if(_mucl[itrk]>=_mucl[PionRec]) PionRec=itrk;
	      
	      MuonFound=true;
	      if(_mucl[itrk]<_mucl[LowestMuCL]) LowestMuCL=itrk;
	      
	      if(Plots){
		//if(SelectionOV[irec]){//MC case: muon at 99%. If data, one should only take long tracks since gamma contamination
		
		if(SelectionFV[irec] && !IsData){

		  if(IsFV && TMath::Abs(TypeOfTrack[irec][itrk])==13 && Sample[irec][itrk]==3) { //stopping muons
		    Pmu_vs_IronDist->Fill(IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio),Momentum[irec][itrk],weight);	  		   
		    double Leq=IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio);
		    if(Leq<6 && Momentum[irec][itrk]>0.5) cout<<Leq<<" plastic="<<PlasticDistance[irec][itrk]<<" angle="<<TrackAngle[irec][itrk]<<" mom="<<Momentum[irec][itrk]<<" last channels="<<LastChannelIX[irec][itrk]<<" "<<LastChannelIY[irec][itrk]<<endl;
		  }
		  
		  if(IsFV && TMath::Abs(TypeOfTrack[irec][itrk])==2212 && Sample[irec][itrk]==3)//stopping protons
		    Pp_vs_IronDist->Fill(IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio),Momentum[irec][itrk],weight);	  		   
		  
		  if(IsFV && TMath::Abs(TypeOfTrack[irec][itrk])==211 && Sample[irec][itrk]==3)//stopping pions
		    Ppi_vs_IronDist->Fill(IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio),Momentum[irec][itrk],weight);	  		   
		  
		  if(TMath::Abs(TypeOfTrack[irec][itrk])==13) hMuCL_TrueMuon->Fill(_mucl[itrk],weight);
		  else if(TMath::Abs(TypeOfTrack[irec][itrk])==211) hMuCL_TruePion->Fill(_mucl[itrk],weight);
		  else if(TMath::Abs(TypeOfTrack[irec][itrk])==2212) hMuCL_TrueProton->Fill(_mucl[itrk],weight);
		  else hMuCL_TrueOthers->Fill(_mucl[itrk],weight);
		  //cout<<"hello, particle is="<<TMath::Abs(TypeOfTrack[irec][itrk])<<endl;
		}
	      }
	    }
	  }


	  if(Selection==2){ // CC1pi
	    if(TMath::Abs(TypeOfTrack[irec][MuonRec])==13) MuonRec_TruePDG_val[0]+=weight;
	    else if(TMath::Abs(TypeOfTrack[irec][MuonRec])==211) MuonRec_TruePDG_val[1]+=weight;
	    else if(TMath::Abs(TypeOfTrack[irec][MuonRec])==2212) MuonRec_TruePDG_val[2]+=weight;
	    else MuonRec_TruePDG_val[3]+=weight;
	    
	    if(Sample[irec][MuonRec]<Sample[irec][PionRec]){
	      // switch PionRec && MuonRec
	      int pion_tmp=MuonRec;
	      int mu_tmp=PionRec;

	      if(TMath::Abs(TypeOfTrack[irec][mu_tmp])==13) MuonRec_TruePDG_switch_val[0]+=weight;
	      else if(TMath::Abs(TypeOfTrack[irec][mu_tmp])==211) MuonRec_TruePDG_switch_val[1]+=weight;
	      else if(TMath::Abs(TypeOfTrack[irec][mu_tmp])==2212) MuonRec_TruePDG_switch_val[2]+=weight;
	      else MuonRec_TruePDG_switch_val[3]+=weight;

	      if(false){
		// ML 2017-02-03 I remove it because it doesn't improve anything and add some phase space limitation
		PionRec=pion_tmp;
		MuonRec=mu_tmp;
	      }
	    }
	  }


	  if(!MuonFound){
	    /*
	    cout<<"Largest CL not found, check number of tracks="<<nTracks[irec]<<", mucl values="<<endl;
	    for(int itrk=0;itrk<nTracks[irec];itrk++){
	      cout<<"("<<_mucl[itrk]<<", is reconstructed="<<IsReconstructed[irec][itrk]<<") ";
	    }
	    cout<<endl;*/
	    NEventsLost+=weight;
	    continue;
	  }
	  int MuonLike=0;int ProtonLike=0;int Undetermined=0;
	    for(int itrk=0;itrk<nTracks[irec];itrk++){
	      if(!IsReconstructed[irec][itrk]) continue;
	      if(_mucl[itrk]>MuonCut) MuonLike++;
	      else if(_mucl[itrk]<ProtonCut && (_isPM? true: _mucl[itrk]>=0)) ProtonLike++;
		      // for WM, mucl=-1 corresponds to undetremined tracks (to few isohits)
	      else Undetermined++;
	    }
	    double rnd=rand->Uniform(0,1);

	    
	    //if((nTracks[irec]==1 && CLMuon_Plan[irec][0]>-3) || ((nTracks[irec]==2 && TMath::Max(CLMuon_Plan[irec][0],CLMuon_Plan[irec][1])>-3 && TMath::Min(CLMuon_Plan[irec][0],CLMuon_Plan[irec][1])<-3) && (TMath::Min(CLMuon_Plan[irec][0],CLMuon_Plan[irec][1])>=0)) || ((nTracks[irec]==3 && CLMuon_Plan[irec][MuonRec]>-3 && CLMuon_Plan[irec][(MuonRec+1)%3]<-3 && CLMuon_Plan[irec][(MuonRec+2)%3]<-3)) || ((nTracks[irec]==4 && CLMuon_Plan[irec][MuonRec]>-3 && CLMuon_Plan[irec][(MuonRec+1)%4]<-3 && CLMuon_Plan[irec][(MuonRec+2)%4]<-3 && CLMuon_Plan[irec][(MuonRec+3)%4]<-3))){
	      
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
      
	    // for CC1pi MuonSample2 is set to 5
	    if(Selection==1 && (Sample[irec][MuonRec]!=MuonSample1)&&(Sample[irec][MuonRec]!=MuonSample2)) continue; // ML done later for CC1pi

	    int nRecTracks=nTracks[irec]-nBadRecTracks;

	    if(Plots){
	      if(nRecTracks==1) hMuCL[FSIInt]->Fill(_mucl[MuonRec],weight);
	      if(nRecTracks==2) hMuCL_2tracks[FSIInt]->Fill(_mucl[MuonRec],_mucl[LowestMuCL],weight);
	      if(nRecTracks==2 && MuonLike==1 && Undetermined==0){
		hMuCL_Lowest[FSIInt]->Fill(_mucl[LowestMuCL],weight);
		if(_mucl[LowestMuCL]<-1){
		  if(FSIInt<3) PE_Lowest_CC0pi->Fill(ProportionHighPE[irec][LowestMuCL],MeanHighPE[irec][LowestMuCL],weight);
		  else PE_Lowest_Other->Fill(ProportionHighPE[irec][LowestMuCL],MeanHighPE[irec][LowestMuCL],weight);		    
		}
	      }
	    }
	    hNTracks[FSIInt]->Fill(nRecTracks,weight);
	    hRecMom[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	    hRecAngle[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	    //if(FSIInt<3) cout<<"CC0pi true, "<<MuonLike<<", "<<Undetermined<<", "<<ProtonLike<<endl;
	    //else cout<<MuonLike<<", "<<Undetermined<<", "<<ProtonLike<<endl;
      	    if(MuonLike==1 && Undetermined==0/*&& nRecTracks<=2*/){//CC0pi
	      //cout<<"Interaction value="<<FSIInt<<", Ion distance="<<IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio)<<endl;

	      //TEMPORARY
	      if(nRecTracks==2) hSampleSecondTrack[FSIInt]->Fill(Sample[irec][LowestMuCL]);
	      //
		
	      if(Selection==1){
		DataSelected[BinRecMom][BinRecAngle]+=weight;
		//cout<<"Bin="<<BinRecMom<<","<<BinRecAngle<<", data="<<DataSelected[BinRecMom][BinRecAngle]<<endl;
		if(!IsData){
		  if(FSIInt<3 && IsNuMu){
		    if(Plots) MCEfficiency_Energy->Fill(Enu,weight);
		    MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
		    Efficiency[BinTrueMom][BinTrueAngle]+=weight;
#ifdef DEBUG
		    if((IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio))<40 && TrueMomentumMuon>1){
		      cout<<"evt number="<<ievt<<", number of tracks="<<nRecTracks<<endl;
		      for(int itrk=0;itrk<nRecTracks;itrk++){
			cout<<"CL="<<_mucl[itrk]<<", distance="<<IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio)<<", sample="<<Sample[irec][itrk]<<endl;
		      }
		    }
#endif
		  }
		  else{
		    BkgSelected[BinRecMom][BinRecAngle]+=weight;
		  }
		}
	      }
	      hRecMom_CC0pi[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	      hRecAngle_CC0pi[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);

	    }
	    else if(MuonLike==2 &&(/*_isPM?true:*/ Undetermined==0)){//Side band CC1pi - for WM I reject Undertermined tracks (mucl=-1)
	      // for PM, no requirement on Undetermined==0 was done -> I make it now
	      
	      nEvents[3]+=weight;
	      nEventsInter[3][FSIInt]+=weight;
	      
	      if(nRecTracks>3) continue;
	      if(Sample[irec][MuonRec]>2){//only ingrid tracks
		hRecMom_CC1pi_full[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
		hRecAngle_CC1pi_full[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
		if(Sample[irec][MuonRec]==3) MuCL_vs_IronDist->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),_mucl[MuonRec],weight);
	      }
	      
	      if(Selection==2){
		nEvents[4]+=weight;
		nEventsInter[4][FSIInt]+=weight;
		DataSelected_full[BinRecMom][BinRecAngle]+=weight;
		if(!IsData){
		  if(FSIInt==3 && IsNuMu){
		    if(Plots) MCEfficiency_Energy->Fill(Enu,weight);
		    MCSelected_full[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
		    Efficiency_full[BinTrueMom][BinTrueAngle]+=weight;
		  }
		  else{
		    BkgSelected_full[BinRecMom][BinRecAngle]+=weight;
		  }
		}
		

		if((Sample[irec][MuonRec]!=MuonSample1)&&(Sample[irec][MuonRec]!=MuonSample2)) continue;
		nEvents[5]+=weight;
		nEventsInter[5][FSIInt]+=weight;
		hRecMom_CC1pi[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
		hRecAngle_CC1pi[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
				
		DataSelected[BinRecMom][BinRecAngle]+=weight;
		if(!IsData){
		  if(FSIInt==3 && IsNuMu){
		    if(Plots) MCEfficiency_Energy->Fill(Enu,weight);
		    MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
		    Efficiency[BinTrueMom][BinTrueAngle]+=weight;
		  }
		  else{
		    BkgSelected[BinRecMom][BinRecAngle]+=weight;
		  }
		}
		if(Sample[irec][MuonRec]==MuonSample1){
		  nEvents[6]+=weight;
		  nEventsInter[6][FSIInt]+=weight;
		  double leq=IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio);
		  hRecMom_CC1pi_restr[FSIInt]->Fill(leq,weight);
		  hRecAngle_CC1pi_restr[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
		}
	      }
	    }
	    else if(MuonLike>=3){//Side band CCNpi
	      hRecMom_CCNpi[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	      hRecAngle_CCNpi[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	    }	     
      }
    }
  }

  double TotalSelected=0,TotalTrue0=0,TotalTrueSelected=0;
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      if(!IsData ){
	  TotalTrueSelected+=Efficiency[c0][c1];

	if(Selection==1 && TotalCC0piEvent[c0][c1]!=0){
	  TotalTrue0+=TotalCC0piEvent[c0][c1];
	  Efficiency[c0][c1]/=TotalCC0piEvent[c0][c1];
	}
	else if(Selection==2 && TotalCC1piEvent[c0][c1]!=0){
	  TotalTrue0+=TotalCC1piEvent[c0][c1];
	  Efficiency[c0][c1]/=TotalCC1piEvent[c0][c1];
	}

	//	cout<<"Efficiency in true pmu bin "<<c0<<" and true thetamu bin "<<c1<<" is: "<<Efficiency[c0][c1]<<endl;
	if(Plots) MCEfficiency->SetBinContent(c0+1,c1+1,Efficiency[c0][c1]);
      }
    }
  }
  cout<<"OVERALL EFFICIENCY="<<TotalTrueSelected/TotalTrue0*100.<<"%"<<endl;
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	TotalSelected+=DataSelected[e0][e1];
	//	cout<<"Purity in rec pmu bin "<<e0<<" and rec thetamu bin "<<e1<<" is: "<<1-BkgSelected[e0][e1]/DataSelected[e0][e1]<<endl;
    }
  }
  cout<<"OVERALL PURITY="<<TotalTrueSelected/TotalSelected*100.<<"%"<<endl;
  cout<<"OVERALL SELECTED EVENTS="<<TotalSelected<<endl;
  cout<<"TOTAL GENERATED CC"<<(Selection==2)<<"PI IN FV="<<TotalTrue0<<endl;
  if(Selection==2){
    for(int i=0;i<nCuts;i++){
      cout<<myCuts[i]<<": \tselected="<<nEvents[i]<<" \tpurity="<<100.*nEventsInter[i][3]/nEvents[i]<<" \tefficiency="<<nEventsInter[i][3]/TotalTrue0*100<<endl;
      CutEfficiency->SetBinContent(i+1,nEventsInter[i][3]/TotalTrue0*100);
      CutPurity->SetBinContent(i+1,nEventsInter[i][3]/nEvents[i]*100);
    }
  }


  if(Plots){
    if(Selection==1) MCEfficiency_Energy->Divide(TotalCC0piEvent_Energy);
    else if(Selection==2) MCEfficiency_Energy->Divide(TotalCC1piEvent_Energy);
    int colors[4]={kBlue,kGreen+3,kRed,kGray};
    MuonRec_TruePDG=new TPie("MuonRecTruePDG","MuonRec true PDG",4,MuonRec_TruePDG_val,colors);
    MuonRec_TruePDG_switch=new TPie("MuonRecTruePDGswitch","MuonRec true PDG",4,MuonRec_TruePDG_switch_val,colors);
    MuonRec_TruePDG->SetLabelFormat("#splitline{%txt}{(%perc)}");
    MuonRec_TruePDG_switch->SetLabelFormat("#splitline{%txt}{(%perc)}");
  }
    
  char OutNameEvent[256];
  char OutNameEvent_root[256];
  ofstream fEvent;
  sprintf(OutNameEvent,"%s%s.txt",outnameevent,(_isPM?"PM":"WM"));
  sprintf(OutNameEvent_root,"%s%s.root",outnameevent,(_isPM?"PM":"WM"));

  if(Systematics_Flux){
    //sprintf(OutNameEvent,"files/MCSelected_Systematics%d_%d.txt",ErrorType,ErrorIteration);
    //sprintf(OutNameEvent_root,"files/MCSelected_Systematics%d_%d.root",ErrorType,ErrorIteration);
    for(int i=1;i<=NBinsEnergyFlux;i++){
      double Nneut=NeutrinoFlux->GetBinContent(i);
      NeutrinoFlux->SetBinContent(i,Nneut+Nneut*FluxVector[i-1]);
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

  //if(!IsData){
    TFile * wfile = new TFile(OutNameEvent_root,"recreate");
    MCTrueEvents->Write();
    CutEfficiency->Write();
    CutPurity->Write();
    wfile->Close();
    //}
  
  delete MCTrueEvents;
  delete CutPurity;
  delete CutEfficiency;
  
  if(Plots){
    TFile * wfile = new TFile(OutNameEvent_root,"update");
    if(!IsData){
      THStack * Stack_MuCL = new THStack("Stack_MuCL","");
      THStack * Stack_MuCL_Lowest = new THStack("Stack_MuCL_Lowest","");
      THStack * Stack_NTracks = new THStack("Stack_NTracks","");
      THStack * Stack_RecMom = new THStack("Stack_RecMom","");
      THStack * Stack_RecAngle = new THStack("Stack_RecAngle","");
      THStack * Stack_FCFVTrueEvents = new THStack("Stack_FCFVTrueEvents","");
      THStack * Stack_RecMom_CC0pi = new THStack("Stack_RecMom_CC0pi","");
      THStack * Stack_RecAngle_CC0pi = new THStack("Stack_RecAngle_CC0pi","");
      THStack * Stack_RecMom_CC1pi_restr = new THStack("Stack_RecMom_CC1pi_restr","");
      THStack * Stack_RecAngle_CC1pi_restr = new THStack("Stack_RecAngle_CC1pi_restr","");
      THStack * Stack_RecMom_CC1pi = new THStack("Stack_RecMom_CC1pi","");
      THStack * Stack_RecAngle_CC1pi = new THStack("Stack_RecAngle_CC1pi","");
      THStack * Stack_RecMom_CC1pi_full = new THStack("Stack_RecMom_CC1pi_full","");
      THStack * Stack_RecAngle_CC1pi_full = new THStack("Stack_RecAngle_CC1pi_full","");
      THStack * Stack_RecMom_CCNpi = new THStack("Stack_RecMom_CCNpi","");
      THStack * Stack_RecAngle_CCNpi = new THStack("Stack_RecAngle_CCNpi","");
      THStack * Stack_SampleSecondTrack = new THStack("Stack_SampleSecondTrack","");

      ProduceStack(hMuCL,Stack_MuCL);
      ProduceStack(hMuCL_Lowest,Stack_MuCL_Lowest);
      //No stack produced for 2tracks 2D plot, but should still normalise as done in the produce stack!
      //for(int fsi=0;fsi<NFSIs;fsi++) hMuCL_2tracks[fsi]->Scale((5.86e20/1e21)*100/1000);
      
      ProduceStack(hNTracks,Stack_NTracks);
      ProduceStack(hRecMom,Stack_RecMom);
      ProduceStack(hRecAngle,Stack_RecAngle);
      ProduceStack(hFCFVTrueEvents,Stack_FCFVTrueEvents);
      ProduceStack(hRecMom_CC0pi,Stack_RecMom_CC0pi);
      ProduceStack(hRecAngle_CC0pi,Stack_RecAngle_CC0pi);
      ProduceStack(hRecMom_CC1pi_restr,Stack_RecMom_CC1pi_restr);
      ProduceStack(hRecAngle_CC1pi_restr,Stack_RecAngle_CC1pi_restr);
      ProduceStack(hRecMom_CC1pi,Stack_RecMom_CC1pi);
      ProduceStack(hRecAngle_CC1pi,Stack_RecAngle_CC1pi);
      ProduceStack(hRecMom_CC1pi_full,Stack_RecMom_CC1pi_full);
      ProduceStack(hRecAngle_CC1pi_full,Stack_RecAngle_CC1pi_full);
      ProduceStack(hRecMom_CCNpi,Stack_RecMom_CCNpi);
      ProduceStack(hRecAngle_CCNpi,Stack_RecAngle_CCNpi);
      
      ProduceStack(hSampleSecondTrack,Stack_SampleSecondTrack);

      TLegend * leg = new TLegend(0.7,0.8,1,1);
      if(hRecMom[0]) leg->AddEntry(hRecMom[0],"CC0#pi");
      if(hRecMom[1]) leg->AddEntry(hRecMom[1],"MEC");
      if(hRecMom[3]) leg->AddEntry(hRecMom[3],"CC1#pi");
      if(hRecMom[4]) leg->AddEntry(hRecMom[4],"CC1#pi^{0}");
      if(hRecMom[5]) leg->AddEntry(hRecMom[5],"CCother");
      if(hRecMom[6]) leg->AddEntry(hRecMom[6],"NC");
            
      //TFile * fMC = new TFile("plots/MCPlots.root","RECREATE");

      Stack_FCFVTrueEvents->Write();
      Stack_MuCL->Write();
      Stack_MuCL_Lowest->Write();
      Stack_NTracks->Write();
      Stack_RecMom->Write();
      Stack_RecAngle->Write();
      Stack_RecMom_CC0pi->Write();
      Stack_RecAngle_CC0pi->Write();
      Stack_RecMom_CC1pi_restr->Write();
      Stack_RecAngle_CC1pi_restr->Write();
      Stack_RecMom_CC1pi->Write();
      Stack_RecAngle_CC1pi->Write();
      Stack_RecMom_CC1pi_full->Write();
      Stack_RecAngle_CC1pi_full->Write();
      Stack_RecMom_CCNpi->Write();
      Stack_RecAngle_CCNpi->Write();
      Stack_SampleSecondTrack->Write();
      Pmu_vs_IronDist->Write();
      Pp_vs_IronDist->Write();
      Ppi_vs_IronDist->Write();
      MuCL_vs_IronDist->Write();
      MCEfficiency->Write();
      MCEfficiency_Energy->Write();
      TotalCC0piEvent_Energy->Write();
      TotalCC1piEvent_Energy->Write();
      leg->Write();
      
      for(int fsi=0;fsi<NFSIs;fsi++){
	hMuCL[fsi]->Write();
	hMuCL_Lowest[fsi]->Write();
	hMuCL_2tracks[fsi]->Write();
	hNTracks[fsi]->Write();
      }
      hMuCL_TrueMuon->Write();
      hMuCL_TruePion->Write();
      hMuCL_TrueProton->Write();
      hMuCL_TrueOthers->Write();
      PE_Lowest_CC0pi->Write();
      PE_Lowest_Other->Write();
      MuonRec_TruePDG->Write();
      MuonRec_TruePDG_switch->Write();
      //fMC->Close();
    }
    else{
      //TFile * fData = new TFile("plots/DataPlots.root","RECREATE");
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

      //fData->Close();
    }
    wfile->Close();
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
    delete MCEfficiency;
    delete MCEfficiency_Energy;
    delete TotalCC0piEvent_Energy;
  }
  delete rand;
}

int main(int argc, char ** argv){
  cout<<"welcome, here is the CC0pi selection"<<endl;
  char * fName = new char[256];
  bool MC=false;
  bool Data=true;
  int Sample=0;
  bool IsPlots=true;
  bool Systematics=false;
  int SelectedError_Source=-1;
  double SelectedError_Variation=-1;
  char * InNameEvent = new char[256];
  char * OutNameEvent = new char[256];
  bool isPM=true;   
  bool Systematics_Flux=false;
  bool Systematics_Xsec=false;
  bool Systematics_Detector=false;
  int ErrorType=0;int n=0;
  int Xsec_dial=0;
  int File_Number=0;
  TVectorD * FluxVector=new TVectorD();
     
  int c=-1;
  while ((c = getopt(argc, argv, "w:ms:e:v:i:o:p:W")) != -1) {
    switch(c){
    case 'w':
      SandReweight=atof(optarg);
      break;
    case 'm':
      MC=true;
      Data=false;
      break;
    case 'W':
      isPM=false;
      break;
    case 's':
      Sample=atoi(optarg);//0 for sand muon, 1 for CC0pi, 2 for CC1pi
      break;
    case 'e':
      Systematics=true;
      SelectedError_Source=atoi(optarg);
    case 'v':
      SelectedError_Variation=atoi(optarg);
      break;
    case 'i':
      InNameEvent=optarg;
      break;
    case 'o':
      OutNameEvent=optarg;
      break;
    case 'p':
      DataEquivalent=atof(optarg);     //cout<<"Please enter the amount of data you'd like to mimic (I will adjust MC stat.) in unit of 1e21 POT:"<<endl;
      break;
    }
  }
    
  cout<<"welcome"<<endl;
  cout<<"Detector is "<<(isPM ? "PM" : "WM" )<<endl;
  //  XS->Xsec::Initialize();
  InitializeGlobal();

  if(!isPM){
    ProtonCut=ProtonCut_WM;
    MuonCut=MuonCut_WM;
  }
  if(Sample==2)    MuonSample2=5;
 
 
   TChain * chain = new TChain("wtree");
   if(Data){
     for(int i=StartRun;i<=EndRun;i++){
     //for(int i=14510;i<=14510;i++){
       for(int j=StartSubRun;j<=EndSubRun;j++){
	 //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_%08d_%04d.root",i,j);	   
	 //sprintf(fName,"root_input/CC0piTree%d.root",i);
	 //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_ReWeight.root",i);
	 //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Old_%d_Plan.root",i);
	 //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root",i);
	 //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_RandomPE.root",i);	
	 fName=Form(InNameEvent,i,j);
	 cout<<fName<<endl;
	 chain->Add(fName);
       }
     }

   }
   else{
     //cout<<"Please enter the amount of data you'd like to mimic (I will adjust MC stat.) in unit of 1e21 POT:"<<endl;
     //cin>>DataEquivalent;
     ScalingMC=DataEquivalent/NMCfiles;//one MC file is equivalent to 1e21 POT
       
     for(int i=0;i<NMCfiles;i++){
       //sprintf(fName,i);
       //sprintf(fName,"%s",InNameEvent);
       //sprintf(fName,"root_input/CC0piTree%d.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_ReWeight.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Old_%d_Plan.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_RandomPE.root",i);

       //       fName=Form(InNameEvent,i);
       sprintf(fName,"${INSTALLREPOSITORY}/XS/root_input/XSFormat_%s_Run1_%d_PlanDev.root",(isPM?"PM":"WM"),i);
       cout<<fName<<endl;
       chain->Add(fName);
     }
   }
   if(SelectedError_Source>=Systematics_Detector_Start && SelectedError_Source<=Systematics_Detector_End){
     Systematics_Detector=true;
     Distributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM);	 
   }
   else if(SelectedError_Source>=Systematics_Flux_Start && SelectedError_Source<=Systematics_Flux_End){
     Systematics_Flux=true;
     TFile * FluxError = new TFile("files/ErrorVarFinal.root");
     FluxVector = new TVectorD();
     cout<<"variation in flux #"<<SelectedError_Variation<<endl;
     sprintf(Name,"Var[%d]",( (int) SelectedError_Variation ));
     FluxVector = (TVectorD*) FluxError->Get(Name);
     File_Number= ( (int) SelectedError_Variation);
     //Distributions(chain,Data,Sample,IsPlots,SelectedError_Source,SelectedError_Variation,OutNameEvent);
     //Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent);
     Distributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM);
     FluxError->Close();
   }
   else if(SelectedError_Source>=Systematics_Xsec_Start && SelectedError_Source<=Systematics_Xsec_End){
     Systematics_Xsec=true;
     Xsec_dial=((int) SelectedError_Variation);
     File_Number=Xsec_dial;
     Distributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM); 
   }
   else{
     cout<<"Simple selection on the MC or data"<<endl;
     Distributions(chain,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM);
   }
   
   return 0;
}

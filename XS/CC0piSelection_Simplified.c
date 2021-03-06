//CAREFUL: ONLY STOPPED SAMPLE RIGHT NOW (SAMPLE==3)
//For no, we have a bug correction in this code, that set CL=1e-30 if this is 0. But try to avoid this if possible -> Check the CUmulative distributions
//For now, isreconstructed is not systematically applied
//For now, I removed the track witdth and matching parameters!
//CAREFUL: I REMOVED THE FIRST MC FILE!

// TO FOLLOW THE IMPORTANT STEPS OF THE CODE: SEARCH FOR THE KEYWORD "KEY" IN THE COMMENTS

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
#include <TProfile2D.h>
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
//MVA libs
//#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
//My libs
#include "setup.h"
//#include "Reconstruction.cc"
//#define DEBUG
//#define DEBUG2
//#define DEBUG3
//#define MVATRAINING
//#define MVAMUONTRAINING
//#define MVAPROTONTRAINING
#define MVAREADING
//#define TEMPORARY
//#define DEBUGMVA
//#define MOMENTUMSHIFT
bool BkgSub=false;

//If wish to cut part of the track for the MVA. This is the starting bin of the relative track length where to use the MVA.
int FirstHit=0;

char Type[32];
char Name[256];
char Name0[256];

double ScalingMC;
double DataEquivalent=1;
double SandReweight=0.6;
double MuonSample1=3;
double MuonSample2=3;//Careful:sample 1 and 2 are ALWAYS THE RESTRAINED SAMPLE (See plot names).
double MuonSample3=3;
double MuonCut=0.9;
double ProtonCut=0.9;
double MuonMVACut;double ProtonMVACut;
double MuonMVACut_CC0pi=0.06;
double ProtonMVACut_CC0pi=0.1;
double MuonMVACut_CC1pi=0.06;
double ProtonMVACut_CC1pi=0.0;

//double ProtonMVACut=0.1;

//Temp for 1D cut
//#define MVA2DCut
double MuonMVACut_slope=0.0048;
double MuonMVACut_origin=-0.10;
const double Limit2DCut=45;

//  
double MuonCut_WM=0.3;// to be tuned
double ProtonCut_WM=0.3;// to be tuned
//

//#define CCCOH_TUNING
//#define PI_LIKELIHOOD

//
// ML added 2017/10/16 -- input/output in degrees
double Get2TrackAngle(double ang1, double thetax1, double thetay1,double ang2, double thetax2, double thetay2){
  double * Track1 = new double[3];
  double * Track2 = new double[3];

  Track1[0]=TMath::Tan(TMath::Pi()/180.*thetax1)*TMath::Cos(TMath::Pi()/180.*ang1);
  Track1[1]=TMath::Tan(TMath::Pi()/180.*thetay1)*TMath::Cos(TMath::Pi()/180.*ang1);
  Track1[2]=TMath::Cos(TMath::Pi()/180.*ang1);
  double NormTrack1=TMath::Sqrt(Track1[0]*Track1[0]+Track1[1]*Track1[1]+Track1[2]*Track1[2]);

  Track2[0]=TMath::Tan(TMath::Pi()/180.*thetax2)*TMath::Cos(TMath::Pi()/180.*ang2);
  Track2[1]=TMath::Tan(TMath::Pi()/180.*thetay2)*TMath::Cos(TMath::Pi()/180.*ang2);
  Track2[2]=TMath::Cos(TMath::Pi()/180.*ang2);
  double NormTrack2=TMath::Sqrt(Track2[0]*Track2[0]+Track2[1]*Track2[1]+Track2[2]*Track2[2]);

  for(int i=0;i<3;i++){
    Track1[i]/=NormTrack1;
    Track2[i]/=NormTrack2;
  }

  double AngleInTracks=180/TMath::Pi()*TMath::ACos(Track1[0]*Track2[0]+Track1[1]*Track2[1]+Track1[2]*Track2[2]);
  return AngleInTracks;

}


const int NSimplifiedPDG=4;
int DetermineSimplifiedPDG(int track_pdg){
  if(TMath::Abs(track_pdg)==13) return 0;
  else if(TMath::Abs(track_pdg)==211) return 1;
  else if(TMath::Abs(track_pdg)==2212) return 2;
  else return 3;
}
//

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
      h[4]->SetFillColor(kAzure+7);//CCpi0
      h[5]->SetFillColor(kBlue+2);//CCnPi+/- (Cc other)
      h[6]->SetFillColor(kGray);//NC
      if(fsi>6 && fsi<=11){
	if(fsi==7) h[7]->SetFillColor(kYellow);//Sand
	if(fsi==8) h[8]->SetFillColor(kYellow+4);//Horizontal Module Bkg
	if(fsi==9) h[9]->SetFillColor(kYellow+7);//Vertical Module Bkg
	if(fsi==10) h[10]->SetFillColor(kMagenta);//AntiNuMu
	if(fsi==11) h[11]->SetFillColor(kGreen);//NuE
	if(fsi==12) h[12]->SetFillColor(kWhite);//NuE
	}
    }
    h[fsi]->SetLineColor(1);
    h[fsi]->SetLineWidth(2);
    hStack->Add(h[fsi]);
  }
  h[0]->SetTitle("CC0#pi (0p)");
  h[1]->SetTitle("CC0#pi (1p)");
  h[2]->SetTitle("CC0#pi (>2p)");
  h[3]->SetTitle("CC1#pi^{#pm}");
  h[4]->SetTitle("CC1#pi^{0}");
  h[5]->SetTitle("CCother");
  h[6]->SetTitle("NC");
  h[7]->SetTitle("Sand");
  h[8]->SetTitle("INGRID Bkg H");
  h[9]->SetTitle("INGRID Bkg V");
  h[10]->SetTitle("anti-#nu_{#mu}");
  h[11]->SetTitle("#nu_{e}");
  h[12]->SetTitle("Scinti. bkg");
}

void ProduceStackParticles(TH1D * hmu, TH1D * hpi, TH1D * hp, TH1D * hothers, THStack * hStack){

  hmu->GetYaxis()->SetTitleOffset(1.3);

  hmu->SetFillColor(kBlue+2);
  //hmu->SetFillStyle(3001);
  hpi->SetFillColor(kGreen+2);
  //hpi->SetFillStyle(3001);
  hp->SetFillColor(kRed);
  //hp->SetFillStyle(3001);
  hothers->SetFillColor(kGray);
  //hothers->SetFillStyle(3001);

  hmu->SetLineColor(1);
  hmu->SetLineWidth(2);
  hStack->Add(hmu);
  hpi->SetLineColor(1);
  hpi->SetLineWidth(2);
  hStack->Add(hpi);
  hp->SetLineColor(1);
  hp->SetLineWidth(2);
  hStack->Add(hp);
  hothers->SetLineColor(1);
  hothers->SetLineWidth(2);
  hStack->Add(hothers);

  hmu->SetTitle("#mu");
  hpi->SetTitle("#pi");
  hp->SetTitle("p");
  hothers->SetTitle("Others");

}

//void InitialiseTable(double DataSelected[NBinsRecMom][NBinsRecAngle],double MCSelected[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle],double BkgSelected[NBinsRecMom][NBinsRecAngle],double Efficiency[NBinsTrueMom][NBinsTrueAngle],double TotalCC0piEvent[NBinsTrueMom][NBinsTrueAngle],double TotalCC1piEvent[NBinsTrueMom][NBinsTrueAngle]){
void InitialiseTable(double ** DataSelected,double **** MCSelected,double ** BkgSelected,double ** Efficiency,double ** TotalCC0piEvent,double ** TotalCC1piEvent){
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


void CC0piDistributions(TChain * wtree,TChain * wtreeMVA, bool IsData,int Selection,bool Plots,bool Systematics_Flux,int File_Number,TVectorD FluxVector,bool Systematics_Xsec,int dial,bool Systematics_Detector,int ErrorType,char * outnameevent, bool _isPM){

  cout<<"hello"<<endl;
  int nevt=(int) wtree->GetEntries();
  cout<<"number of events="<<nevt<<endl;
  int err;
  int sig;
  TH1D * NeutrinoFlux;
  TRandom3 * rand = new TRandom3();


  ///////////////////////////////LOAD THE NEUTRINO FLUX
  //TFile * f = new TFile("files/nd34_tuned_11bv3.1_250ka.root");
  //NeutrinoFlux = (TH1D*) f->Get("ing3_tune_numu");
  TFile * f = new TFile("Flux/tune_ingbg.root");
  NeutrinoFlux = (TH1D*) f->Get("nd2_tune_numu");
  NeutrinoFlux->SetDirectory(0);//Detach the histogram from the TFile (one can close the file)
  f->Close();

  //TO DO
  const int nCuts=7;
  string myCuts[nCuts]={"Generated","Reconstruted + one INGRID track","FCFV","CC1pi Selection","2 or 3 tracks","MuTrk is INGRID stop/through","MuTrk is INGRID stop"};//For CC1pi
  if(Selection==1){
    myCuts[3]="CC0pi Selection";
    myCuts[4]="1 or 2 tracks";
  }
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
  TH1D * CutNEvents[NFSIs];
  for(int i=0;i<NFSIs;i++) CutNEvents[i] = new TH1D(Form("CutNEvents%d",i),Form("Number of events after each cut for the fsi %d",i),nCuts,0,nCuts);
  //END TO DO


  
  //##################################################################
  //##################################################################
  // 0. PREPARE THE TREE FOR READING, B. QUILAIN 2017/05/08
  //##################################################################
  //##################################################################

  
  //##################################################################
  // 0.0 PREPARE THE INPUT TREE FROM XS_CC0PI_PLAN

  int FSIInt;//Final state true information
  int Num_Int;//Interaction number
  int nTracks;//number of tracks reconstructed per vertex
  float weight;//weight of each event
  bool IsFV;//True vertex is in FV
  bool IsSand;
  bool IsAnti;
  bool IsBkgH;
  bool IsBkgV;
  bool IsNuE;
  bool IsDetected;//Event is reconstructed (means >=1 vertex reconstructed)
  bool SelectionFV;bool SelectionOV;//Vertex is reconstructed within/out of FV
  float OpeningAngle;//Opening angle between reconstructed tracks - only relevant for 2 track samples
  float CoplanarityAngle;//Coplanarity angle between reconstructed tracks - only relevant for 2 track samples
  float Enu;//True neutrino energy
  int nIngBasRec;//Number of vertices reconstructed per event
  bool NewEvent;//true if new event, false if it is the same event than previous but only another reconstructed one.
  float Errorweight;
  float TrueMomentumMuon, TrueAngleMuon;//True muon kinematic variables
  float TrueMomentumPion, TrueAnglePion;//True pion kinematic variables
  float POT;//Number of POT
  float TrackAngle[LimitTracks];//Reconstructed angle of the reconstructed track
  float TrackThetaY[LimitTracks];//Reconstructed 2D angle of the reconstructed track
  float TrackThetaX[LimitTracks];//Reconstructed 2D angle of the reconstructed track
  int TypeOfTrack[LimitTracks];//pdg of the reconstructed track. Careful that it is based on an algorithm defined in Reconstruction.cc
  float CLMuon[LimitTracks];// = new vector<float> [LimitTracks];
  float CLMuon_Plan[LimitTracks];// = new vector<float> [LimitTracks];
  float CLMuon_KS[LimitTracks];// = new vector<float> [LimitTracks];
  float CLMuon_Likelihood[LimitTracks];// = new vector<float> [LimitTracks];
  float TotalCharge[LimitTracks];//Charge of the track per unit distance
  float Momentum[LimitTracks];
  float IronDistance[LimitTracks];
  float PlasticDistance[LimitTracks];
  int Sample[LimitTracks];//Geometric properties of the track, defined in Reconstruction::SelectTrackSample
  bool IsReconstructed[LimitTracks];
  float ReWeight[NDials];
  float TrackWidth[LimitTracks];
  float EnergyDeposition[LimitTracks][LimitHits];// = new vector<float> [LimitTracks];
  int Spill;
  int GoodSpill;
  float ChosenCL[LimitTracks];// = new vector<float> [LimitTracks];
  
  /////////////////////////////////PREPARE FOR READING/////////////////////////////////////////
  
  wtree->SetBranchAddress("Spill",&Spill);
  wtree->SetBranchAddress("GoodSpill",&GoodSpill);
  wtree->SetBranchAddress("FSIInt",&FSIInt);
  wtree->SetBranchAddress("nIngBasRec",&nIngBasRec);
  wtree->SetBranchAddress("InteractionType",&Num_Int);
  wtree->SetBranchAddress("nTracks",&nTracks);
  wtree->SetBranchAddress("NewEvent",&NewEvent);
  wtree->SetBranchAddress("nIngBasRec",&nIngBasRec);
  wtree->SetBranchAddress("InteractionType",&Num_Int);  
  wtree->SetBranchAddress("nTracks",&nTracks);
  wtree->SetBranchAddress("OpeningAngle",&OpeningAngle);
  wtree->SetBranchAddress("CoplanarityAngle",&CoplanarityAngle);
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
  wtree->SetBranchAddress("SelectionFV",&SelectionFV);
  wtree->SetBranchAddress("SelectionOV",&SelectionOV);
  wtree->SetBranchAddress("IsDetected",&IsDetected);
  wtree->SetBranchAddress("OpeningAngle",&OpeningAngle);
  wtree->SetBranchAddress("CoplanarityAngle",&CoplanarityAngle);
  wtree->SetBranchAddress(Form("ReWeight[%d]",NDials),ReWeight);
  wtree->SetBranchAddress("POT",&POT);


  for(int itrk=0;itrk<LimitTracks;itrk++){
    wtree->SetBranchAddress(Form("TrackAngle_track%d",itrk),&(TrackAngle[itrk]));
    wtree->SetBranchAddress(Form("TrackThetaY_track%d",itrk),&(TrackThetaY[itrk]));
    wtree->SetBranchAddress(Form("TrackThetaX_track%d",itrk),&(TrackThetaX[itrk]));
    wtree->SetBranchAddress(Form("TrackWidth_track%d",itrk),&(TrackWidth[itrk]));
    wtree->SetBranchAddress(Form("TypeOfTrack_track%d",itrk),&(TypeOfTrack[itrk]));
    wtree->SetBranchAddress(Form("IsReconstructed_track%d",itrk),&(IsReconstructed[itrk]));
    wtree->SetBranchAddress(Form("Sample_track%d",itrk),&(Sample[itrk]));
    wtree->SetBranchAddress(Form("CLMuon_track%d",itrk),&(CLMuon[itrk]));
    wtree->SetBranchAddress(Form("CLMuon_Plan_track%d",itrk),&(CLMuon_Plan[itrk]));
    wtree->SetBranchAddress(Form("CLMuon_Likelihood_track%d",itrk),&(CLMuon_Likelihood[itrk]));
    wtree->SetBranchAddress(Form("CLMuon_KS_track%d",itrk),&(CLMuon_KS[itrk]));
    wtree->SetBranchAddress(Form("IronDistance_track%d",itrk),&(IronDistance[itrk]));
    wtree->SetBranchAddress(Form("PlasticDistance_track%d",itrk),&(PlasticDistance[itrk]));
    wtree->SetBranchAddress(Form("TotalCharge_track%d",itrk),&(TotalCharge[itrk]));
    wtree->SetBranchAddress(Form("Momentum_track%d",itrk),&(Momentum[itrk]));
  
    for(int ihit=0; ihit<LimitHits;ihit++){
      wtree->SetBranchAddress(Form("EnergyDeposition_track%d_hit%d",itrk,ihit),&(EnergyDeposition[itrk][ihit]));
    }      
  }
  // 0.0. END
  //##################################################################


  
  //##################################################################
  // 0.1 PREPARE THE INPUT TREE IF WE USE MVA PID (BDT CURRENTLY)

#ifdef MVAREADING
    cout<<"Start reading the BDT"<<endl;
    //I do not use a pointer because of memory leak in TMVA package for now, B. Quilain
    TMVA::Reader tmvareader;
    
    //Add the relevant variables
    float SampleTrackMVA;
    float EnergyDepositionMVA[LimitHits];
    float  EquivalentIronDistanceTrackMVA;
    tmvareader.AddVariable("Sample",&SampleTrackMVA);
    for(int ihit=FirstHit;ihit<LimitHits;ihit++){
      
      tmvareader.AddVariable(Form("EnergyDeposition_hit%d",ihit),&(EnergyDepositionMVA[ihit]));
    }
    //tmvareader.AddVariable(Form("(IronDistance+(PlasticDistance/(%3.3f)))",IronCarbonRatio),&EquivalentIronDistanceTrackMVA);
    tmvareader.AddSpectator(Form("(IronDistance+(PlasticDistance/(%3.3f)))",IronCarbonRatio),&EquivalentIronDistanceTrackMVA);
    
    //Define spectator
    float TrackAngleTrackMVA, TypeOfTrackMVA, CLMuon_LikelihoodTrackMVA, TotalChargeTrackMVA, MomentumTrackMVA, IsReconstructedTrackMVA;
    float FSIIntMVA, SpillMVA, GoodSpillMVA, NewEventMVA, nIngBasRecMVA, InteractionTypeMVA, EnuMVA, TrueMomentumMuonMVA, TrueAngleMuonMVA, IsFVMVA, IsSandMVA, IsAntiMVA, IsBkgHMVA, IsBkgVMVA, IsNuEMVA, SelectionFVMVA, SelectionOVMVA, IsDetectedMVA, POTMVA;
   
    tmvareader.AddSpectator("FSIInt",&FSIIntMVA);
    tmvareader.AddSpectator("Spill",&SpillMVA);
    tmvareader.AddSpectator("GoodSpill",&GoodSpillMVA);
    tmvareader.AddSpectator("NewEvent",&NewEventMVA);
    tmvareader.AddSpectator("nIngBasRec",&nIngBasRecMVA);
    tmvareader.AddSpectator("InteractionType",&InteractionTypeMVA);
    tmvareader.AddSpectator("Enu",&EnuMVA);
    tmvareader.AddSpectator("TrueMomentumMuon",&TrueMomentumMuonMVA);
    tmvareader.AddSpectator("TrueAngleMuon",&TrueAngleMuonMVA);
    tmvareader.AddSpectator("IsFV",&IsFVMVA);
    tmvareader.AddSpectator("IsSand",&IsSandMVA);
    tmvareader.AddSpectator("IsAnti",&IsAntiMVA);
    tmvareader.AddSpectator("IsBkgH",&IsBkgHMVA);
    tmvareader.AddSpectator("IsBkgV",&IsBkgVMVA);
    tmvareader.AddSpectator("IsNuE",&IsNuEMVA);
    tmvareader.AddSpectator("SelectionFV",&SelectionFVMVA);
    tmvareader.AddSpectator("SelectionOV",&SelectionOVMVA);
    tmvareader.AddSpectator("IsDetected",&IsDetectedMVA);
    tmvareader.AddSpectator("POT",&POTMVA);
    
    tmvareader.AddSpectator("TrackAngle",&TrackAngleTrackMVA);
    tmvareader.AddSpectator("TypeOfTrack",&TypeOfTrackMVA);
    tmvareader.AddSpectator("CLMuon_Likelihood",&CLMuon_LikelihoodTrackMVA);
    tmvareader.AddSpectator("TotalCharge",&TotalChargeTrackMVA);
    tmvareader.AddSpectator("Momentum",&MomentumTrackMVA);
    tmvareader.AddSpectator("IsReconstructed",&IsReconstructedTrackMVA);

  
    cout<<"End of assigning variables for the read BDT"<<endl;
    
    //tmvareader->BookMVA(TMVA::Types::kBDT,"weights/TMVAClassification_BDT.weights.xml");
    tmvareader.BookMVA("BDT method","weights/TMVAClassificationMuon_BDT.weights.xml");
    cout<<"End of initial reading the BDT"<<endl;


    cout<<"Start reading the Second BDT"<<endl;
    //I do not use a pointer because of memory leak in TMVA package for now
    TMVA::Reader tmvareader2;
    
    //Add the relevant variables
    tmvareader2.AddVariable("Sample",&SampleTrackMVA);
    for(int ihit=FirstHit;ihit<LimitHits;ihit++){   
      tmvareader2.AddVariable(Form("EnergyDeposition_hit%d",ihit),&(EnergyDepositionMVA[ihit]));
    }
    
    //Define spectator
    tmvareader2.AddSpectator("FSIInt",&FSIIntMVA);
    tmvareader2.AddSpectator("Spill",&SpillMVA);
    tmvareader2.AddSpectator("GoodSpill",&GoodSpillMVA);
    tmvareader2.AddSpectator("NewEvent",&NewEventMVA);
    tmvareader2.AddSpectator("nIngBasRec",&nIngBasRecMVA);
    tmvareader2.AddSpectator("InteractionType",&InteractionTypeMVA);
    tmvareader2.AddSpectator("Enu",&EnuMVA);
    tmvareader2.AddSpectator("TrueMomentumMuon",&TrueMomentumMuonMVA);
    tmvareader2.AddSpectator("TrueAngleMuon",&TrueAngleMuonMVA);
    tmvareader2.AddSpectator("IsFV",&IsFVMVA);
    tmvareader2.AddSpectator("IsSand",&IsSandMVA);
    tmvareader2.AddSpectator("IsAnti",&IsAntiMVA);
    tmvareader2.AddSpectator("IsBkgH",&IsBkgHMVA);
    tmvareader2.AddSpectator("IsBkgV",&IsBkgVMVA);
    tmvareader2.AddSpectator("IsNuE",&IsNuEMVA);
    tmvareader2.AddSpectator("SelectionFV",&SelectionFVMVA);
    tmvareader2.AddSpectator("SelectionOV",&SelectionOVMVA);
    tmvareader2.AddSpectator("IsDetected",&IsDetectedMVA);
    tmvareader2.AddSpectator("POT",&POTMVA);
    
    tmvareader2.AddSpectator("TrackAngle",&TrackAngleTrackMVA);
    tmvareader2.AddSpectator("TypeOfTrack",&TypeOfTrackMVA);
    tmvareader2.AddSpectator("CLMuon_Likelihood",&CLMuon_LikelihoodTrackMVA);
    tmvareader2.AddSpectator("TotalCharge",&TotalChargeTrackMVA);
    tmvareader2.AddSpectator(Form("(IronDistance+(PlasticDistance/(%3.3f)))",IronCarbonRatio),&EquivalentIronDistanceTrackMVA);
    tmvareader2.AddSpectator("Momentum",&MomentumTrackMVA);
    tmvareader2.AddSpectator("IsReconstructed",&IsReconstructedTrackMVA);

  
    cout<<"End of assigning variables for the read BDT"<<endl;
    
    //tmvareader2->BookMVA(TMVA::Types::kBDT,"weights/TMVAClassification_BDT.weights.xml");
    tmvareader2.BookMVA("BDT method","weights/TMVAClassificationProton_BDT.weights.xml");
    cout<<"End of initial reading the second BDT"<<endl;

#endif
  //Data will be sent from the wtree -> tmvareader! So we should prepare variables like "EquivalentIronDistance"
  //How to send data from the wtree -> tmvareader
  //double MVA = tmvareader->EvaluateMVA( "BDT",0);
  // 0.1. END
  //##################################################################

    
  // 0. END
  //##################################################################
  //##################################################################








     
  //##################################################################
  //##################################################################
  // 1. TRAINING OF THE MVA, B. QUILAIN 2017/05/08
  //   -> THIS TRAINING SHOULD BE DONE WITH PARTICLE GUN INPUT
  //##################################################################
  //##################################################################

#ifdef MVATRAINING
    
  //##############################TRAINING CASE######################################
  cout<<"Starting to train a BDT"<<endl;
  TCut preselection = "(IsFV) && (SelectionFV) && (IsDetected) && (IsReconstructed)";//Only MC in FV is used for training
  TMVA::Factory * factory;
  
#ifdef MVAMUONTRAINING  
  //Create the factory/prepare the forest
  //TFile * MVAoutputMuon = new TFile("src/MVAparticleMuon_1000trees_BugTrueParticleSolved.root","RECREATE");//Output file Name
  TFile * MVAoutputMuon = new TFile( (_isPM)?"src/MVAparticleMuon_1000trees.root":"src/MVAWMparticleMuon_1000trees.root","RECREATE");//Output file Name
  factory = new TMVA::Factory((_isPM)?"TMVAClassificationMuon":"TMVAWMClassificationMuon", MVAoutputMuon,"!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  //Add the relevant variables
  factory->AddVariable("Sample","Sample of the track", "", 'I' );
  for(int ihit=FirstHit;ihit<LimitHits;ihit++){
    factory->AddVariable(Form("EnergyDeposition_hit%d",ihit),Form("Energy depositon of hit %d",ihit), "", 'F' );
  }
  //factory->AddVariable(Form("EquivalentIronDistance := ( IronDistance + (PlasticDistance/(%3.3f)) )",IronCarbonRatio),'F');
  factory->AddSpectator(Form("EquivalentIronDistance := (IronDistance+(PlasticDistance/(%3.3f)))",IronCarbonRatio),"distance in iron of track",'F');

  //factory->AddVariable(Form("EquivalentIronDistance := IronDistance + PlasticDistance / 7.81"),'F');

  //Add spectator variables
  factory->AddSpectator("FSIInt","interaction ID after FSI",'I');
  factory->AddSpectator("Spill","spill number",'I');
  factory->AddSpectator("GoodSpill","good spill flag",'I');
  factory->AddSpectator("NewEvent","flag if it is a new simulated/data event",'I');
  factory->AddSpectator("nIngBasRec","number of reconstructed vertexes for the event",'I');
  factory->AddSpectator("InteractionType","interaction ID at the vertex",'I');
  factory->AddSpectator("Enu","E_{#nu} true",'F');
  factory->AddSpectator("TrueMomentumMuon","p_{#mu] true",'F');
  factory->AddSpectator("TrueAngleMuon","#theta_{#mu] true",'F');
  factory->AddSpectator("IsFV","flag for true vertex in FV",'I');
  factory->AddSpectator("IsSand","flag for sand #mu",'I');
  factory->AddSpectator("IsAnti","flag for #overline{#nu}_{#mu}",'I');
  factory->AddSpectator("IsBkgH","flag for bkg from INGRID Horiz.",'I');
  factory->AddSpectator("IsBkgV","flag for bkg from INGRID Vert.",'I');
  factory->AddSpectator("IsNuE","flag for #nu_{e}",'I');
  factory->AddSpectator("SelectionFV","flag for reconstructed in FV",'I');
  factory->AddSpectator("SelectionOV","flag for reconstructed out of FV",'I');
  factory->AddSpectator("IsDetected","flag for reconstructed vertex",'I');
  factory->AddSpectator("POT","number of corresponding POTs",'I');

  factory->AddSpectator("TrackAngle","angle of the tracke w.r.t z axis",'F');
  factory->AddSpectator("TypeOfTrack","pdf of track",'F');
  factory->AddSpectator("CLMuon_Likelihood","#mu_{CL} of track",'F');
  factory->AddSpectator("TotalCharge","Total dE/dx of track",'F');
  factory->AddSpectator("Momentum","true momentum of track",'F');
  factory->AddSpectator("IsReconstructed","reconstruction status of track",'I');
  
  //factory->AddSpectator("ReWeight[175]","reweight",'F');

  //Case 1: Classification
  //Define the signal and background trees.
  TCut signalMuon = "TypeOfTrack == 13 || TypeOfTrack == -13";
  factory->SetInputTrees(wtreeMVA,signalMuon && preselection,(!signalMuon) && preselection);

  //Add weights to the tree
  factory->SetSignalWeightExpression("weight");
  factory->SetBackgroundWeightExpression("weight");

  //Book the trees <=> prepare the forest
  //Define the method to use for MVA:
  factory->BookMethod( TMVA::Types::kBDT,"BDT","!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
  
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  
  //Keep in mind here that I only use 1/2 of trees for training and 1/2 for testing (default). This option can be changed!
  cout<<"End of BDT training"<<endl;
  
  MVAoutputMuon->Close();
  factory->Delete();
#endif

#ifdef MVAPROTONTRAINING  
  //Second BDT training
  cout<<"Start the training of the second BDT: 0pi is wished (pi0 or pi+), only protons -> protons vs pions (and other particles)"<<endl;
  //TFile * MVAoutputProton = new TFile("src/MVAparticleProton_1000trees_BugTrueParticleSolved_90PerCentEventsUsedForTraining.root","RECREATE");//Output file Name
  TFile * MVAoutputProton = new TFile((_isPM)?"src/MVAparticleProton_1000trees.root":"src/MVAWMparticleProton_1000trees.root","RECREATE");//Output file Name
  factory = new TMVA::Factory((_isPM)?"TMVAClassificationProton":"TMVAWMClassificationProton", MVAoutputProton,"!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  //Add the relevant variables
  factory->AddVariable("Sample","Sample of the track", "", 'I' );
  for(int ihit=FirstHit;ihit<LimitHits;ihit++){  
      factory->AddVariable(Form("EnergyDeposition_hit%d",ihit),Form("Energy deposition of hit %d",ihit), "", 'F' );
  }

    //Add spectator variables
    factory->AddSpectator("FSIInt","interaction ID after FSI",'I');
    factory->AddSpectator("Spill","spill number",'I');
    factory->AddSpectator("GoodSpill","good spill flag",'I');
    factory->AddSpectator("NewEvent","flag if it is a new simulated/data event",'I');
    factory->AddSpectator("nIngBasRec","number of reconstructed vertexes for the event",'I');
    factory->AddSpectator("InteractionType","interaction ID at the vertex",'I');
    factory->AddSpectator("Enu","E_{#nu} true",'F');
    factory->AddSpectator("TrueMomentumMuon","p_{#mu] true",'F');
    factory->AddSpectator("TrueAngleMuon","#theta_{#mu] true",'F');
    factory->AddSpectator("IsFV","flag for true vertex in FV",'I');
    factory->AddSpectator("IsSand","flag for sand #mu",'I');
    factory->AddSpectator("IsAnti","flag for #overline{#nu}_{#mu}",'I');
    factory->AddSpectator("IsBkgH","flag for bkg from INGRID Horiz.",'I');
    factory->AddSpectator("IsBkgV","flag for bkg from INGRID Vert.",'I');
    factory->AddSpectator("IsNuE","flag for #nu_{e}",'I');
    factory->AddSpectator("SelectionFV","flag for reconstructed in FV",'I');
    factory->AddSpectator("SelectionOV","flag for reconstructed out of FV",'I');
    factory->AddSpectator("IsDetected","flag for reconstructed vertex",'I');
    factory->AddSpectator("POT","number of corresponding POTs",'I');
    
    factory->AddSpectator("TrackAngle","angle of the tracke w.r.t z axis",'F');
    factory->AddSpectator("TypeOfTrack","pdf of track",'F');
    factory->AddSpectator("CLMuon_Likelihood","#mu_{CL} of track",'F');
    factory->AddSpectator("TotalCharge","Total dE/dx of track",'F');
    factory->AddSpectator(Form("EquivalentIronDistance := (IronDistance+(PlasticDistance/(%3.3f)))",IronCarbonRatio),"distance in iron of track",'F');
    factory->AddSpectator("Momentum","true momentum of track",'F');
    factory->AddSpectator("IsReconstructed","reconstruction status of track",'I');
    
    //Case 1: Classification
    //Define the signal and background trees.
    TCut signalProton = "TypeOfTrack == 2212 || TypeOfTrack == -2212";
    factory->SetInputTrees(wtreeMVA,signalProton && preselection,(!(signalProton)) && preselection);

    //Add weights to the tree
    factory->SetSignalWeightExpression("weight");
    factory->SetBackgroundWeightExpression("weight");

    //New 2018/02/16 B.Q
    //int nTrain = (int) 0.9*nevt;
    //int nTest = nevt - nTrain;
    //cout<<"Number of events used for training = "<<nTrain<<", and for testing = "<<nTest;
    //factory->PrepareTrainingAndTestTree("","nTrain_Signal=nTrain:nTest_Signal=nTest:nTrain_Background=nTrain:nTest_Background=nTest");
    //
    
    factory->BookMethod( TMVA::Types::kBDT,"BDT","!H:V:NTrees=10:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    
    //
    // Train MVAs using the set of training events
    factory->TrainAllMethods();
    
    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();
    
    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();
    
    
    //Keep in mind here that I only use 1/2 of trees for training and 1/2 for testing (default). This option can be changed!
    
    cout<<"End of second BDT training"<<endl;

    MVAoutputProton->Close();
    factory->Delete();
#endif
    //#################################################################################
    cout<<"Finished after training"<<endl;
#endif

    // 1. END
    //##################################################################
    //##################################################################



    //TO DO
    /*
  TH1D * FluxDistribution = new TH1D("FluxDistribution","",NBinsEnergyFlux,BinningEnergyFlux);
  for(int i=0;i<NBinsEnergyFlux;i++){
    if(Systematics_Flux) FluxDistribution->SetBinContent(i+1,hNormNuMu_Binned->GetBinContent(i+1)+hNormNuMu_Binned->GetBinContent(i+1)*FluxVector[i]);
    else FluxDistribution->SetBinContent(i+1,0.);
  }
*/
#ifdef DEBUG2
  cout<<"Right after MVA declaration"<<endl;
#endif
  //END TO DO







  
  
  //##################################################################
  //##################################################################
  // 2. PREPARE THE HISTOGRAMS BEFORE EVENT LOOPS, B. QUILAIN 2017/05/08
  //##################################################################
  //##################################################################


  
  // KEY
  //##################################################################
  // 2.0. VARIABLES BEING PASSED TO UNFOLDING (OR OTHER METHODS) THROUGH
  // ROOT AND TXT FILE -> THEY ARE COMPULSORY!

  //Difference between variables and variables_full: the first can have a track sample selection (example: INGRID stop) while the second has no track sample selection
  double POTCount=0;//Number of POT (for data)

  double ** DataSelected = new double*[NBinsRecMom];
  double ** BkgSelected = new double*[NBinsRecMom];
  double ** DataSelected_full = new double*[NBinsRecMom];
  double ** BkgSelected_full = new double*[NBinsRecMom];
  for(int i=0;i<NBinsRecMom;i++){
    DataSelected[i] = new double[NBinsRecAngle];
    BkgSelected[i] = new double[NBinsRecAngle];
    DataSelected_full[i] = new double[NBinsRecAngle];
    BkgSelected_full[i] = new double[NBinsRecAngle];
  }

  double ** Efficiency = new double*[NBinsTrueMom];
  double ** Efficiency_full = new double*[NBinsTrueMom];
  double ** TotalCC0piEvent = new double*[NBinsTrueMom];
  double ** TotalCC1piEvent = new double*[NBinsTrueMom];
  for(int i=0;i<NBinsTrueMom;i++){
    Efficiency[i] = new double[NBinsTrueAngle];
    Efficiency_full[i] = new double[NBinsTrueAngle];
    TotalCC0piEvent[i] = new double[NBinsTrueAngle];
    TotalCC1piEvent[i] = new double[NBinsTrueAngle];
  }

  double **** MCSelected = new double***[NBinsTrueMom];
  double **** MCSelected_full = new double***[NBinsTrueMom];
  for(int i=0;i<NBinsTrueMom;i++){
    MCSelected[i] = new double**[NBinsTrueAngle];
    MCSelected_full[i] = new double**[NBinsTrueAngle];
    for(int j=0;j<NBinsTrueAngle;j++){
      MCSelected[i][j] = new double*[NBinsRecMom];
      MCSelected_full[i][j] = new double*[NBinsRecMom];
      for(int k=0;k<NBinsRecMom;k++){
	MCSelected[i][j][k] = new double[NBinsRecAngle];
	MCSelected_full[i][j][k] = new double[NBinsRecAngle];
      }
    }
  }

  double TotalCC1pi=0;

  InitialiseTable(DataSelected,MCSelected,BkgSelected,Efficiency,TotalCC0piEvent,TotalCC0piEvent);//Set the variables to 0
  InitialiseTable(DataSelected_full,MCSelected_full,BkgSelected_full,Efficiency_full,TotalCC0piEvent,TotalCC0piEvent);//Set the variables to 0

  TH2D * MCTrueEvents = new TH2D("MCTrueEvents","Number of true simulated CC0pi (if selection==1) or CC1pi (if selection==2) in the FV",NBinsTrueMom,BinningTrueMom,NBinsTrueAngle,BinningTrueAngle);
  // 2.0. END
  //##################################################################

  

  
  // KEY
  //##################################################################
  // 2.1. NOT MANDATORY: ALL HISTOGRAMS TO DRAW PLOTS FOR CUT TUNING,
  // DRAWING PLOTS ETC...

  
  // 2.1.0. TUNING THE CUTS ON PID:
  TH1D * hMuCL[NFSIs];// 1 track sample without any selection on the track type
  TH2D * hMuCL_2tracks[NFSIs];// 2 track sample without any selection on the track type
  TH1D * hMuCL_Lowest[NFSIs];// lowest mucl of the 2 track sample, but only for the CC0pi reconstructed sample

  TH1D * hMVAMuondiscriminant_1track[NFSIs];
  TH1D * hMVAProtondiscriminant_2tracks[NFSIs];
  TH2D * hMVAMuonVSProtondiscriminant_2tracks[NFSIs];
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track[NFSIs];
  TH2D * hMVAMuondiscriminantVSPDG_1track[NFSIs];

  TH2D * hMVAMuondiscriminantVSDistance_1track[NFSIs];
  TH2D * hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[NFSIs];
  
  TH1D * TrueParticleType_1track[NFSIs];
  TH2D * TrueParticleType_2tracks[NFSIs];

  TH2D * hMVA_MuonVsProtonDiscriminant_TrueMuon;
  TH2D * hMVA_MuonVsProtonDiscriminant_TruePion;
  TH2D * hMVA_MuonVsProtonDiscriminant_TrueProton;
  TH1D * hMuCL_TrueMuon;
  TH1D * hMuCL_TruePion;
  TH1D * hMuCL_TrueProton;
  TH1D * hMuCL_TrueOthers;
  TH2D * hMVAMuondiscriminantVSDistance_TrueParticle[NSimplifiedPDG];

  TH1D * hLMVA_1track[NFSIs];
  TH2D * hLMVA_2tracks[NFSIs];
  TH1D * hLMVA_LargestLMVA[NFSIs];
  TH1D * hLMVA_SmallestLMVA[NFSIs];

  TH1D * hLMVA_TrueParticle[NSimplifiedPDG][NSamples];

  // 2.1.0. END
  
  TH1D * hNTracks[NFSIs];
  TH1D * hRecMom[NFSIs];
  TH1D * hRecAngle[NFSIs];
  TH1D * hFCFVTrueEvents[NFSIs];

  TH1D * hNTracks_CC0pi[NFSIs];
  TH1D * hRecMom_CC0pi[NFSIs];
  TH1D * hRecAngle_CC0pi[NFSIs];
  TH1D * hRecMom_CC0pi_restr[NFSIs];
  TH1D * hRecAngle_CC0pi_restr[NFSIs];
  TH1D * hRecMom_CC0pi_full[NFSIs];
  TH1D * hRecAngle_CC0pi_full[NFSIs];

  TH1D * hRecMom_CC1pi[NFSIs];
  TH1D * hRecAngle_CC1pi[NFSIs];
  TH1D * hRecMom_CC1pi_restr[NFSIs];
  TH1D * hRecAngle_CC1pi_restr[NFSIs];
  TH1D * hRecMom_CC1pi_full[NFSIs];
  TH1D * hRecAngle_CC1pi_full[NFSIs];

  TH1D * hRecMom_CCNpi[NFSIs];
  TH1D * hRecAngle_CCNpi[NFSIs];

  TH1D * hSampleSecondTrack[NFSIs];
  TH1D * hTotalChargeSecondTrack[NFSIs];
  TH1D * hTotalChargePerDistanceSecondTrack[NFSIs];
  TH1D * hTotalChargePerDistanceFirstTrack[NFSIs];
  TH1D * hOpeningAngle[NFSIs];
  TH1D * hCoplanarityAngle[NFSIs];

  TH2D * hRecCC0pi_MomentumResolution;
  TH2D * hRecCC0pi_AngleResolution;
  
  //Energy Deposition
  TH3D * EnergyDepositionLength_Muon[NSamples];
  TH2D * hEnergyDepositionLength_Muon[NSamples];
  TProfile2D * pEnergyDepositionLength_Muon[NSamples];
  TProfile * pEnergyDepositionLength_Muon_1D[NSamples];
  TH1D * hEnergyDepositionLength_Muon_1D[NSamples];

  TH3D * EnergyDepositionLength_Pion[NSamples];
  TH2D * hEnergyDepositionLength_Pion[NSamples];
  TProfile2D * pEnergyDepositionLength_Pion[NSamples];
  TProfile * pEnergyDepositionLength_Pion_1D[NSamples];
  TH1D * hEnergyDepositionLength_Pion_1D[NSamples];

  TH3D * EnergyDepositionLength_Proton[NSamples];
  TH2D * hEnergyDepositionLength_Proton[NSamples];
  TProfile2D * pEnergyDepositionLength_Proton[NSamples];
  TProfile * pEnergyDepositionLength_Proton_1D[NSamples];
  TH1D * hEnergyDepositionLength_Proton_1D[NSamples];



  

  TH2D * MCEfficiency_PTheta_ThinBinning;
  TH2D * TotalCC0piEvent_PTheta_ThinBinning;
  TH2D * TotalCC1piEvent_PTheta_ThinBinning;
  TH1D * MCEfficiency_Momentum_ThinBinning;
  TH1D * TotalEvent_Momentum_ThinBinning;
  TH1D * MCEfficiency_Angle_ThinBinning;
  TH1D * TotalEvent_Angle_ThinBinning;

  TH2D * MCReconstructionEfficiency_PTheta_ThinBinning;
  TH1D * MCReconstructionEfficiency_Momentum_ThinBinning;
  TH1D * MCReconstructionEfficiency_Angle_ThinBinning;

  TH2D * MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning;
  TH1D * MCReconstructionAndTopologyEfficiency_Momentum_ThinBinning;
  TH1D * MCReconstructionAndTopologyEfficiency_Angle_ThinBinning;

  /*
    TH1D * MCEfficiency_Momentum_ThinBinning_Relative;
    TH1D * MCEfficiency_Angle_ThinBinning_Relative;
    TH1D * MCReconstructionEfficiency_Momentum_ThinBinning_Relative;
    TH1D * MCReconstructionEfficiency_Angle_ThinBinning_Relative;
    TH1D * MCReconstructionAndTopologyEfficiency_Momentum_ThinBinning_Relative;
    TH1D * MCReconstructionAndTopologyEfficiency_Angle_ThinBinning_Relative;
  */
  
  TH2D * MCEfficiency;
  TH2D * MCEfficiency_PTheta;
  TH1D * MCEfficiency_PTheta_Momentum;
  TH1D * MCEfficiency_PTheta_Angle;

  TH1D * MCEfficiency_Energy;
  TH1D * TotalCC0piEvent_Energy;
  TH1D * TotalCC1piEvent_Energy;

  TPie * MuonRec_TruePDG;
  TPie * MuonRec_TruePDG_switch;
  double MuonRec_TruePDG_val[4];  
  double MuonRec_TruePDG_switch_val[4];  

  TH2D* Pmu_vs_IronDist, *Pp_vs_IronDist, *Ppi_vs_IronDist;
  TH2D* MuCL_vs_IronDist;

  TH2D * MuonID = new TH2D("MuonID","",50,0.,2.,30,0.,90.);
  TH2D * MuonIDMVA = new TH2D("MuonIDMVA","",50,0.,2.,30,0.,90.);
  TH2D * MuonIDTotal = new TH2D("MuonIDTotal","",50,0.,2.,30,0.,90.);
    
  if(Plots){
    if(_isPM){
      hMVA_MuonVsProtonDiscriminant_TrueMuon = new TH2D("hMVA_MuonVsProtonDiscriminant_TrueMuon","Muon vs Proton discriminant for the true muons",100,-.5,.5,100,-.5,.5);//ML tmp 500bins->100bins
      hMVA_MuonVsProtonDiscriminant_TruePion = new TH2D("hMVA_MuonVsProtonDiscriminant_TruePion","Muon vs Proton discriminant for the true pions",100,-.5,.5,100,-.5,.5);//ML tmp 500bins->100bins
      hMVA_MuonVsProtonDiscriminant_TrueProton = new TH2D("hMVA_MuonVsProtonDiscriminant_TrueProton","Muon vs Proton discriminant for the true protons",100,-.5,.5,100,-.5,.5);//ML tmp 500bins->100bins
      
      hMuCL_TrueMuon = new TH1D("hMuCL_TrueMuon","Muon confidence level for the sand muons",50,0,1);//ML tmp 500bins->100bins
      hMuCL_TruePion = new TH1D("hMuCL_TruePion","Muon confidence level for the true muons",50,0,1);//ML tmp 500bins->100bins
      hMuCL_TrueProton = new TH1D("hMuCL_TrueProton","Muon confidence level for the true protons",50,0,1);//ML tmp 500bins->100bins
      hMuCL_TrueOthers = new TH1D("hMuCL_TrueOthers","Muon confidence level for the true protons",50,0,1);//ML tmp 500bins->100bins
      MuCL_vs_IronDist=new TH2D("MuCL_vs_IronDist","highest #mu_{CL} track",50,0,100,20,MuonCut,0);
    }
    else {
      hMVA_MuonVsProtonDiscriminant_TrueMuon = new TH2D("hMVA_MuonVsProtonDiscriminant_TrueMuon","Muon vs Proton discriminant for the sand muons",100,-.5,.5,100,-.5,.5);//ML tmp 500bins->100bins
      hMVA_MuonVsProtonDiscriminant_TruePion = new TH2D("hMVA_MuonVsProtonDiscriminant_TruePion","Muon vs Proton discriminant for the sand muons",100,-.5,.5,100,-.5,.5);//ML tmp 500bins->100bins
      hMVA_MuonVsProtonDiscriminant_TrueProton = new TH2D("hMVA_MuonVsProtonDiscriminant_TrueProton","Muon vs Proton discriminant for the sand muons",100,-.5,.5,100,-.5,.5);//ML tmp 500bins->100bins

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

     
    for(int is=0;is<NSamples;is++){    
      EnergyDepositionLength_Muon[is] = new TH3D(Form("EnergyDepositionLength_Muon_%d",is),"",20,0,100,LimitHits,0,1.,100,0,200);
      EnergyDepositionLength_Muon[is]->GetXaxis()->SetTitle("Penetration in iron (cm)");
      EnergyDepositionLength_Muon[is]->GetYaxis()->SetTitle("Distance in track length unit");
      EnergyDepositionLength_Muon[is]->GetZaxis()->SetTitle("dE/dx (p.e/cm)");
      EnergyDepositionLength_Muon[is]->Sumw2();
	
      EnergyDepositionLength_Pion[is] = new TH3D(Form("EnergyDepositionLength_Pion_%d",is),"",20,0,100,LimitHits,0,1.,100,0,200);
      EnergyDepositionLength_Pion[is]->GetXaxis()->SetTitle("Penetration in iron (cm)");
      EnergyDepositionLength_Pion[is]->GetYaxis()->SetTitle("Distance in track length unit");
      EnergyDepositionLength_Pion[is]->GetZaxis()->SetTitle("dE/dx (p.e/cm)");
      EnergyDepositionLength_Pion[is]->Sumw2();
      
      EnergyDepositionLength_Proton[is] = new TH3D(Form("EnergyDepositionLength_Proton_%d",is),"",20,0,100,LimitHits,0,1.,100,0,200);
      EnergyDepositionLength_Proton[is]->GetXaxis()->SetTitle("Penetration in iron (cm)");
      EnergyDepositionLength_Proton[is]->GetYaxis()->SetTitle("Distance in track length unit");
      EnergyDepositionLength_Proton[is]->GetZaxis()->SetTitle("dE/dx (p.e/cm)");
      EnergyDepositionLength_Proton[is]->Sumw2();

    }
    
    MCEfficiency_PTheta_ThinBinning = new TH2D("MCEfficiency_PTheta_ThinBinning","",200,0,10,60,0,180);
    MCEfficiency_PTheta_ThinBinning->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
    MCEfficiency_PTheta_ThinBinning->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
    TotalCC0piEvent_PTheta_ThinBinning = new TH2D("TotalCC0piEvent_PTheta_ThinBinning","",200,0,10,60,0,180);
    TotalCC0piEvent_PTheta_ThinBinning->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
    TotalCC0piEvent_PTheta_ThinBinning->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
    TotalCC1piEvent_PTheta_ThinBinning = new TH2D("TotalCC1piEvent_PTheta_ThinBinning","",200,0,10,60,0,180);
    TotalCC1piEvent_PTheta_ThinBinning->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
    TotalCC1piEvent_PTheta_ThinBinning->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
    MCEfficiency_PTheta_ThinBinning->Sumw2();TotalCC0piEvent_PTheta_ThinBinning->Sumw2();TotalCC1piEvent_PTheta_ThinBinning->Sumw2();
    
    MCReconstructionEfficiency_PTheta_ThinBinning = new TH2D("MCReconstructionEfficiency_PTheta_ThinBinning","",200,0,10,60,0,180);
    MCReconstructionEfficiency_PTheta_ThinBinning->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
    MCReconstructionEfficiency_PTheta_ThinBinning->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
    MCReconstructionEfficiency_PTheta_ThinBinning->Sumw2();TotalCC0piEvent_PTheta_ThinBinning->Sumw2();TotalCC1piEvent_PTheta_ThinBinning->Sumw2();
    
    MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning = new TH2D("MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning","",200,0,10,60,0,180);
    MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
    MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
    MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Sumw2();TotalCC0piEvent_PTheta_ThinBinning->Sumw2();TotalCC1piEvent_PTheta_ThinBinning->Sumw2();
    
    MCEfficiency = new TH2D("MCEfficiency","",NBinsTrueMomSignal,0,NBinsTrueMomSignal,NBinsTrueAngleSignal,0,NBinsTrueAngleSignal);
    MCEfficiency_PTheta = new TH2D("MCEfficiency","",NBinsTrueMomSignal,BinningTrueMomSignal,NBinsTrueAngleSignal,BinningTrueAngleSignal);
    MCEfficiency_PTheta_Momentum = new TH1D("MCEfficiency_PTheta_Momentum","",NBinsTrueMomSignal,BinningTrueMomSignal);
    MCEfficiency_PTheta_Angle = new TH1D("MCEfficiency_PTheta_Angle","",NBinsTrueAngleSignal,BinningTrueAngleSignal);

    MCEfficiency_Energy = new TH1D("MCEfficiency_Energy","",100,0,10);
    MCEfficiency_Energy->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    TotalCC0piEvent_Energy = new TH1D("TotalCC0piEvent_Energy","",100,0,10);
    TotalCC0piEvent_Energy->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    TotalCC1piEvent_Energy = new TH1D("TotalCC1piEvent_Energy","",100,0,10);
    TotalCC1piEvent_Energy->GetXaxis()->SetTitle("E_{#nu} (GeV)");
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
      hMuCL[i] = new TH1D(Name,Title,500,0.,1.);
      //hMuCL->Sumw2();
      hMuCL[i]->GetXaxis()->SetTitle("#mu_{CL}");    hMuCL[i]->GetYaxis()->SetTitle("Number of events");
	
      sprintf(Name,"hMuCL_Lowest%d",i);
      sprintf(Title,"Muon confidence level for the %dth-fsi, only for the lowest MuCL tracks if ntracks>1",i);
      hMuCL_Lowest[i] = new TH1D(Name,Title,500,0.,1.);
      //hMuCL_Lowest->Sumw2();
      hMuCL_Lowest[i]->GetXaxis()->SetTitle("#mu_{CL_Lowest}");    hMuCL_Lowest[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hMuCL_2tracks%d",i);
      sprintf(Title,"Muon confidence level for the %dth-fsi",i);
      hMuCL_2tracks[i] = new TH2D(Name,Title,500,0.,1.,500,0.,1.);
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
	
      sprintf(Name,"hNTracks_CC0pi%d",i);
      sprintf(Title,"Number of tracks for the %dth-fsi",i);
      hNTracks_CC0pi[i] = new TH1D(Name,Title,5,0,5);
      //hNTracks_CC0pi->Sumw2();
      hNTracks_CC0pi[i]->GetXaxis()->SetTitle("Number of reconstructed tracks");
      hNTracks_CC0pi[i]->GetYaxis()->SetTitle("Number of events");

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
	
      sprintf(Name,"hRecMom_CC0pi_restr%d",i);
      sprintf(Title,"Distance in iron for the %dth-fsi",i);
      hRecMom_CC0pi_restr[i] = new TH1D(Name,Title,20,0,100);
      //hRecMom_CC0pi_restr->Sumw2();
      hRecMom_CC0pi_restr[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
      hRecMom_CC0pi_restr[i]->GetYaxis()->SetTitle("Number of events");
	
      sprintf(Name,"hRecAngle_CC0pi_restr%d",i);
      sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
      hRecAngle_CC0pi_restr[i] = new TH1D(Name,Title,30,0,90);
      //hRecAngle_CC0pi_restr->Sumw2();
      hRecAngle_CC0pi_restr[i]->GetXaxis()->SetTitle("Angle_CC0pi (°)");
      hRecAngle_CC0pi_restr[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hRecMom_CC0pi_full%d",i);
      sprintf(Title,"Distance in iron for the %dth-fsi",i);
      hRecMom_CC0pi_full[i] = new TH1D(Name,Title,20,0,100);
      //hRecMom_CC0pi_full->Sumw2();
      hRecMom_CC0pi_full[i]->GetXaxis()->SetTitle("Equivalent length in iron (cm)");
      hRecMom_CC0pi_full[i]->GetYaxis()->SetTitle("Number of events");
	
      sprintf(Name,"hRecAngle_CC0pi_full%d",i);
      sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
      hRecAngle_CC0pi_full[i] = new TH1D(Name,Title,30,0,90);
      //hRecAngle_CC0pi_full->Sumw2();
      hRecAngle_CC0pi_full[i]->GetXaxis()->SetTitle("Angle_CC0pi (°)");
      hRecAngle_CC0pi_full[i]->GetYaxis()->SetTitle("Number of events");
      //
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
      //hSampleSecondTrack[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hSampleSecondTrack[i] = new TH1D(Name,Title,6,0,6);
      //hSampleSecondTrack->Sumw2();
      hSampleSecondTrack[i]->GetXaxis()->SetTitle("Track sample");
      hSampleSecondTrack[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hTotalChargeSecondTrack%d",i);
      sprintf(Title,"Number of tracks for the %dth-fsi",i);
      //hTotalChargeSecondTrack[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hTotalChargeSecondTrack[i] = new TH1D(Name,Title,100,0,3000);
      //hTotalChargeSecondTrack->Sumw2();
      hTotalChargeSecondTrack[i]->GetXaxis()->SetTitle("Track charge / dx");
      hTotalChargeSecondTrack[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hTotalChargePerDistanceSecondTrack%d",i);
      sprintf(Title,"Number of tracks for the %dth-fsi",i);
      //hTotalChargePerDistanceSecondTrack[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hTotalChargePerDistanceSecondTrack[i] = new TH1D(Name,Title,120,0,200);
      //hTotalChargePerDistanceSecondTrack->Sumw2();
      hTotalChargePerDistanceSecondTrack[i]->GetXaxis()->SetTitle("Track Charge / Total Distance");
      hTotalChargePerDistanceSecondTrack[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hTotalChargePerDistanceFirstTrack%d",i);
      sprintf(Title,"Number of tracks for the %dth-fsi",i);
      //hTotalChargePerDistanceFirstTrack[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hTotalChargePerDistanceFirstTrack[i] = new TH1D(Name,Title,120,0,200);
      //hTotalChargePerDistanceFirstTrack->Sumw2();
      hTotalChargePerDistanceFirstTrack[i]->GetXaxis()->SetTitle("Track Charge / Total Distance");
      hTotalChargePerDistanceFirstTrack[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hOpeningAngle%d",i);
      sprintf(Title,"Opening angle for the %dth-fsi",i);
      //hOpeningAngle[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hOpeningAngle[i] = new TH1D(Name,Title,60,0,180);
      //hOpeningAngle->Sumw2();
      hOpeningAngle[i]->GetXaxis()->SetTitle("Opening angle (#circ)");
      hOpeningAngle[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hCoplanarityAngle%d",i);
      sprintf(Title,"Coplanarity angle for the %dth-fsi",i);
      //hCoplanarityAngle[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hCoplanarityAngle[i] = new TH1D(Name,Title,60,0,180);
      //hCoplanarityAngle->Sumw2();
      hCoplanarityAngle[i]->GetXaxis()->SetTitle("Coplanarity angle (#circ)");
      hCoplanarityAngle[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hMVAMuondiscriminant_1track%d",i);
      sprintf(Title,"MVAMuon discriminant of the 1 track sample for the %dth-fsi",i);
      //hMVAMuondiscriminant_1track[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAMuondiscriminant_1track[i] = new TH1D(Name,Title,100,-0.5,0.5);
      //hMVAMuondiscriminant_1track->Sumw2();
      hMVAMuondiscriminant_1track[i]->GetXaxis()->SetTitle("#mu_{MVA}");
      hMVAMuondiscriminant_1track[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hMVAMuondiscriminantVSPDG_1track%d",i);
      sprintf(Title,"MVAMuon discriminant of the 1 track sample for the %dth-fsi",i);
      //hMVAMuondiscriminant_1track[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAMuondiscriminantVSPDG_1track[i] = new TH2D(Name,Title,4,0,4,100,-0.5,0.5);
      //hMVAMuondiscriminant_1track->Sumw2();
      hMVAMuondiscriminantVSPDG_1track[i]->GetXaxis()->SetTitle("PDG");
      hMVAMuondiscriminantVSPDG_1track[i]->GetYaxis()->SetTitle("#mu_{MVA}");

      sprintf(Name,"hMVAProtondiscriminant_2tracks%d",i);
      sprintf(Title,"MVAProton discriminant of the 2 track sample (only for the track having the highest MVAProton) for the %dth-fsi",i);
      //hMVAProtondiscriminant_2tracks[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAProtondiscriminant_2tracks[i] = new TH1D(Name,Title,100,-0.5,0.5);
      //hMVAProtondiscriminant_2tracks->Sumw2();
      hMVAProtondiscriminant_2tracks[i]->GetXaxis()->SetTitle("p_{MVA} of the track having maximal p_{MVA}");
      hMVAProtondiscriminant_2tracks[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hMVAMuonVSProtondiscriminant_2tracks%d",i);
      sprintf(Title,"MVA discriminant for the 2 track sample: MVAMuon of most muon-like track vs MVAProton of most p-like track for the %dth-fsi",i);
      //hMVAMuonVSProtondiscriminant_2tracks[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAMuonVSProtondiscriminant_2tracks[i] = new TH2D(Name,Title,100,-0.5,0.5,100,-0.5,0.5);
      //hMVAMuonVSProtondiscriminant_2tracks->Sumw2();
      hMVAMuonVSProtondiscriminant_2tracks[i]->GetXaxis()->SetTitle("#mu_{MVA} of the track having maximal #mu_{MVA}");
      hMVAMuonVSProtondiscriminant_2tracks[i]->GetYaxis()->SetTitle("p_{MVA} of the track having maximal p_{MVA}");

      sprintf(Name,"hMVAMuondiscriminantVSMuonMomentum_1track%d",i);
      sprintf(Title,"MVAMuon discriminantVSMuonMomentum of the 1 track sample for the %dth-fsi",i);
      //hMVAMuondiscriminantVSMuonMomentum_1track[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAMuondiscriminantVSMuonMomentum_1track[i] = new TH2D(Name,Title,40,0,2,100,-0.5,0.5);
      //hMVAMuondiscriminantVSMuonMomentum_1track->Sumw2();
      hMVAMuondiscriminantVSMuonMomentum_1track[i]->GetXaxis()->SetTitle("True p_{#mu}");
      hMVAMuondiscriminantVSMuonMomentum_1track[i]->GetYaxis()->SetTitle("#mu_{MVA}");

      sprintf(Name,"hMVAMuondiscriminantVSDistance_1track%d",i);
      sprintf(Title,"MVAMuon discriminantVSDistance of the 1 track sample for the %dth-fsi",i);
      //hMVAMuondiscriminantVSDistance_1track[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAMuondiscriminantVSDistance_1track[i] = new TH2D(Name,Title,20,0,100,100,-0.5,0.5);
      //hMVAMuondiscriminantVSDistance_1track->Sumw2();
      hMVAMuondiscriminantVSDistance_1track[i]->GetXaxis()->SetTitle("d_{#mu}");
      hMVAMuondiscriminantVSDistance_1track[i]->GetYaxis()->SetTitle("#mu_{MVA}");

      sprintf(Name,"hMVAProtondiscriminantVSDistance_2tracks_LowestMVA%d",i);
      sprintf(Title,"MVAProton discriminantVSDistance of the 1 track sample for the %dth-fsi",i);
      //hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[i] = new TH2D(Name,Title,20,0,100,100,-0.5,0.5);
      //hMVAProtondiscriminantVSDistance_2tracks_LowestMVA->Sumw2();
      hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[i]->GetXaxis()->SetTitle("Plastic distance p_{p}");
      hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[i]->GetYaxis()->SetTitle("p_{MVA}");

  
      sprintf(Name,"TrueParticleType_1track%d",i);
      sprintf(Title,"Associated true particle to the 1 track sample for the %dth-fsi",i);
      TrueParticleType_1track[i] = new TH1D(Name,Title,4,0,4);
      //TrueParticleType_1track->Sumw2();
      TrueParticleType_1track[i]->GetXaxis()->SetTitle("Simplified PDG");
      TrueParticleType_1track[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"TrueParticleType_2tracks%d",i);
      sprintf(Title,"Associated true particle to the 2 tracks sample for the %dth-fsi",i);
      TrueParticleType_2tracks[i] = new TH2D(Name,Title,4,0,4,4,0,4);
      //TrueParticleType_2tracks->Sumw2();
      TrueParticleType_2tracks[i]->GetXaxis()->SetTitle("Simplified PDG of track 1");
      TrueParticleType_2tracks[i]->GetYaxis()->SetTitle("Simplified PDG of track 2");


      
      sprintf(Name,"hLMVA_1track%d",i);
      sprintf(Title,"MVAMuon discriminant of the 1 track sample for the %dth-fsi",i);
      //hLMVA_1track[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hLMVA_1track[i] = new TH1D(Name,Title,150,-1.5,1.5);
      //hLMVA_1track->Sumw2();
      hLMVA_1track[i]->GetXaxis()->SetTitle("Discriminator");
      hLMVA_1track[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hLMVA_LargestLMVA%d",i);
      sprintf(Title,"MVAMuon discriminant of the 1 track sample for the %dth-fsi",i);
      //hLMVA_LargestLMVA[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hLMVA_LargestLMVA[i] = new TH1D(Name,Title,150,-1.5,1.5);
      //hLMVA_LargestLMVA->Sumw2();
      hLMVA_LargestLMVA[i]->GetXaxis()->SetTitle("Discriminator");
      hLMVA_LargestLMVA[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hLMVA_SmallestLMVA%d",i);
      sprintf(Title,"MVAMuon discriminant of the 1 track sample for the %dth-fsi",i);
      //hLMVA_SmallestLMVA[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hLMVA_SmallestLMVA[i] = new TH1D(Name,Title,150,-1.5,1.5);
      //hLMVA_SmallestLMVA->Sumw2();
      hLMVA_SmallestLMVA[i]->GetXaxis()->SetTitle("Discriminator");
      hLMVA_SmallestLMVA[i]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"hLMVA_2tracks%d",i);
      sprintf(Title,"LMVA discriminant of the 2 track sample for the %dth-fsi",i);
      //hLMVA_2tracks[i] = new TH2D(Name,Title,6,0,6,40,0,2000);
      hLMVA_2tracks[i] = new TH2D(Name,Title,150,-1.5,1.5,150,-1.5,1.5);
      //hLMVA_2tracks->Sumw2();
      hLMVA_2tracks[i]->GetXaxis()->SetTitle("Discriminator 1st track");
      hLMVA_2tracks[i]->GetYaxis()->SetTitle("Discriminator 2nd track");
      
    }

    hRecCC0pi_MomentumResolution = new TH2D("hRecCC0pi_MomentumResolution","",50,0,100,150,0,5);
    hRecCC0pi_AngleResolution = new TH2D("hRecCC0pi_AngleResolution","",30,0,90,30,0,90);

    for(int i=0;i<NSimplifiedPDG;i++){
      sprintf(Name,"hMVAMuondiscriminantVSDistance_TrueParticle%d",i);
      sprintf(Title,"Associated true particle to the 1 track sample for the %dth-particle",i);
      hMVAMuondiscriminantVSDistance_TrueParticle[i] = new TH2D(Name,Title,20,0,100,100,-0.5,0.5);
      //hMVAMuondiscriminantVSDistance_TrueParticle->Sumw2();
      hMVAMuondiscriminantVSDistance_TrueParticle[i]->GetXaxis()->SetTitle("d_{#mu}");
      hMVAMuondiscriminantVSDistance_TrueParticle[i]->GetYaxis()->SetTitle("#mu_{MVA}");

      for(int j=0;j<NSamples;j++){
	sprintf(Name,"hLMVA_TrueParticle%d_Sample%d",i,j);
	sprintf(Title,"Associated true particle to the %dth-particle ans sample %d",i);
	hLMVA_TrueParticle[i][j] = new TH1D(Name,Title,150,-1.5,1.5);
	//hLMVA_TrueParticle->Sumw2();
	hLMVA_TrueParticle[i][j]->GetXaxis()->SetTitle("Discriminator");
	hLMVA_TrueParticle[i][j]->GetYaxis()->SetTitle("Number of events");
      }
    }
  }
  
  // 2.1. END
  //##################################################################

  // 2. END
  //##################################################################
  //##################################################################









  
  double NEvents=0;double NEventsFV=0;double NEventsLost=0;double NEventsLostDetected=0;

  int BinTrueInteraction=0;
  int BinTrueMom=0;
  int BinTrueAngle=0;
  //cout<<"going to read events"<<endl;






#ifdef MOMENTUMSHIFT
  TF1 * Shift = new TF1("Shift","[0]+[1]*x",TotalCC0piEvent_PTheta_ThinBinning->GetXaxis()->GetXmin(),TotalCC0piEvent_PTheta_ThinBinning->GetXaxis()->GetXmax());
  Shift->SetParameter(0,0.);
  Shift->SetParameter(1,1./2.);//So it reach a weight y=1 at x=2 GeV
#endif
  
  //##################################################################
  //##################################################################
  // START EVENT LOOP
  //##################################################################
  //##################################################################

  for(int ievt=0;ievt<nevt;ievt++){
    POTCount+=POT;

    //##################################################################
    //##################################################################
    // CASE OF A NEW SIMULATED EVENT, A NEW SIMVERTEX.
    //##################################################################
    //##################################################################

    wtree->GetEntry(ievt);
#ifdef MOMENTUMSHIFT
    if(TrueMomentumMuon<=2){
      weight *= Shift->Eval(TrueMomentumMuon);
      if(ievt%5000==0) cout<<"pmu="<<TrueMomentumMuon<<", weight is scaled by "<<Shift->Eval(TrueMomentumMuon)<<endl;
    }
#endif

    if(ievt%50000==0){
      cout<<"Event "<<ievt<<"/"<<nevt<<" processed"<<endl;
      cout<<"Percent of events in the true FV (%)="<<100.*NEventsFV/NEvents<<", Number of events lost due to detection (%)="<<100.*NEventsLostDetected/NEvents<<", due to detection & muon not found="<<100.*NEventsLost/NEvents<<endl;
      cout<<"POT="<<POTCount<<endl;
      if(ievt==(nevt-1)) cout<<"Final POT="<<POTCount<<endl;
    }
    //##################################################################
    //##################################################################



    
    //##################################################################
    //##################################################################
    // PREPARE THE REWEIGHTING FACTOR FOR DIFFERENT CASES
    //##################################################################
    //##################################################################
    
    if(IsData){FSIInt=0;Num_Int=0;}
    else {
      //********************
      if(IsSand||IsBkgH||IsBkgV){
	if(weight>1000) continue;
	// continue;
      }
      //******************** tmp
      weight*=ScalingMC;//Adjust the MC distribution to the amount of data we process
    }
    if(IsSand) weight*=(1+SandReweight);

    ////////////////DETERMINE THE REWEIGHTING OF THE EVENT IF SYSTEMATICS ERROR FLUX
    if(Systematics_Flux){
      for(int i=0;i<NBinsEnergyFlux;i++){
	if(Enu<BinningEnergyFlux[i+1]) {
	  Errorweight=FluxVector[i];
	  break;
	}
      }
      if(ievt%10000==0)cout<<"weight="<<weight<<", and after reweight="<<weight*(1+Errorweight)<<endl;
      weight*=(1+Errorweight);
    }
    else if(Systematics_Xsec){
      if(ievt%10000==0)cout<<"weight="<<weight<<", and after reweight="<<weight*ReWeight[dial]<<endl;
      weight=weight*ReWeight[dial];
    }

    //##################################################################
    //##################################################################



    

    //##################################################################
    //##################################################################
    // SELECT THE PID: Bayesian likelihood or Boosted Decision Tree
    //##################################################################
    //##################################################################

    //Select your CL
    for(int itrk=0;itrk<LimitTracks;itrk++){
      //ChosenCL[itrk]=(_isPM? CLMuon_Plan[itrk] : CLMuon_KS[itrk]);
      ChosenCL[itrk]=(_isPM? CLMuon_Likelihood[itrk] : CLMuon_KS[itrk]);
      
      // ML tmp
      /*
	if(_isPM) {
	ChosenCL[itrk]=CLMuon_Likelihood[itrk];
	ProtonCut=0.3;
	MuonCut=0.4;
	}
      */
    }
    //##################################################################
    //##################################################################

    
    bool IsNuMu=((!IsSand) && (!IsAnti) && (!IsNuE) && (!IsBkgH) && (!IsBkgV));
    //cout<<"Is data ="<<IsData<<", good spill="<<GoodSpill<<", spill="<<Spill<<endl;
    if(IsData && (!GoodSpill || !Spill)) continue;

    if((/*IsFV &&*/ !IsData) || IsData){ 
      // why that condition here ?? too easy to remove all external background this way

#ifdef DEBUG3
      cout<<"Entering within the FV, nbasrec="<<nIngBasRec<<endl;
#endif

      //##################################################################
      // If New Event (=new simulated vertex or new DAQ time window)
      // Determine the bin corresponding to the true neutrino energy, and
      // true momentum angle and momentum.

      // Fill also MC events
      
      if(NewEvent){
	if(FSIInt<3) BinTrueInteraction=1;
	else if(FSIInt==3) BinTrueInteraction=2;
	else BinTrueInteraction=3;
	
	BinTrueMom=0;
	BinTrueAngle=0;

	if(FSIInt<3){
	  for(int i=0;i<NBinsTrueMomSignal;i++){
	    if(TrueMomentumMuon<BinningTrueMomSignal[i+1]){BinTrueMom=i;break;}
	  }
	  for(int i=0;i<NBinsTrueAngleSignal;i++){
	    if(TrueAngleMuon<BinningTrueAngleSignal[i+1]){BinTrueAngle=i;break;}
	  }
	}//
	else if(FSIInt==3){
	  for(int i=0;i<NBinsTrueMomSB;i++){
	    if(TrueMomentumMuon<BinningTrueMomSB[i+1]){BinTrueMom=i+NBinsTrueMomSignal;break;}
	  }
	  for(int i=0;i<NBinsTrueAngleSB;i++){
	    if(TrueAngleMuon<BinningTrueAngleSB[i+1]){BinTrueAngle=i+NBinsTrueAngleSignal;break;}
	  }
	}
	else{
	  for(int i=0;i<NBinsTrueMomTrash;i++){
	    if(TrueMomentumMuon<BinningTrueMomTrash[i+1]){BinTrueMom=i+NBinsTrueMomSignal+NBinsTrueMomSB;break;}
	  }
	  for(int i=0;i<NBinsTrueAngleTrash;i++){
	    if(TrueAngleMuon<BinningTrueAngleTrash[i+1]){BinTrueAngle=i+NBinsTrueAngleSignal+NBinsTrueAngleSB;break;}
	  }
	}

	
	//cout<<FSIInt<<", "<<Num_Int<<", bin true mom="<<BinTrueMom<<", bin true angle="<<BinTrueAngle<<endl;

      
	if(!IsData){

	  if(IsNuMu /*&& IsFV*/){
	    if(IsFV) hFCFVTrueEvents[FSIInt]->Fill(Enu,weight);
	    nEvents[0]+=weight;
	    nEventsInter[0][FSIInt]+=weight;
	    
	    if(FSIInt<3){
	      if(Selection == 1) {
		TotalCC0piEvent[BinTrueMom][BinTrueAngle]+=weight;
		MCTrueEvents->Fill(BinTrueMom+1,BinTrueAngle+1,weight);
	      }
	      if(Plots){
		TotalCC0piEvent_Energy->Fill(Enu,weight);
		TotalCC0piEvent_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	      }
	    }
	    else if(FSIInt == 3){
	      if(Selection ==2){
		TotalCC1pi+=weight;
		TotalCC1piEvent[BinTrueMom][BinTrueAngle]+=weight;
		MCTrueEvents->Fill(BinTrueMom+1,BinTrueAngle+1,weight);
	      }
	      if(Plots){
		TotalCC1piEvent_Energy->Fill(Enu,weight);
		TotalCC1piEvent_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	      }
	    }
	  }
	}
      }//End of NewEvent

      
      bool MuonFound=false;
      int MuonTrue=0;int MuonRec=0;
      int PionTrue=0; int PionRec=0;
      int MuonTrueMVA=0;int MuonRecMVA=0;
      int LowestMuCL=0;
      bool Trash=false;
      
      double LMVA[nTracks];
      double MVAdiscriminant[nTracks];
      double MVAdiscriminant2[nTracks];
      double MVAdiscriminant3[nTracks];
      
      double LargestMVAdiscriminant_1track=-999;
      double LargestMVAdiscriminant2_1track=-999;
      double LargestMVAdiscriminant_2tracks=-999;
      double LargestMVAdiscriminant2_2tracks=-999;

      double LargestLMVA=-999;
      double SmallestLMVA=999;
      
      // if(IsData && IsDetected!=0) cout<<"Is detected = "<<IsDetected<<", Is FV rec = "<<SelectionFV<<endl;
      //First cuts!
      NEvents+=weight;
      if(IsFV) NEventsFV += weight;
      if(!IsDetected) {NEventsLostDetected+=weight; continue;}
      nEvents[1]+=weight;
      nEventsInter[1][FSIInt]+=weight;
      
      if(Selection>=1 && !SelectionFV) continue;
      if(Selection==0 && !SelectionOV) continue;
      nEvents[2]+=weight;
      nEventsInter[2][FSIInt]+=weight;
      

      bool AllTracksWellReconstructed=true;
      int nBadRecTracks=0;
      for(int itrk=0;itrk<nTracks;itrk++){
	if(!IsReconstructed[itrk]){
	  AllTracksWellReconstructed=false;
	  nBadRecTracks++;
	}
      }

#ifdef TEMPORARY
      bool OneINGRIDTrack=false;
      for(int itrk=0;itrk<nTracks;itrk++){
	//if(IsReconstructed[itrk]){
	  if(IronDistance[itrk]>0){OneINGRIDTrack=true;break;}
	  //}
      }
      //if(!OneINGRIDTrack){cout<<"No INGRID track"<<endl; continue;}
#endif
	     
      //if(!AllTracksWellReconstructed){cout<<"one track not well reconstructed"<<endl;continue;}
      	
#ifdef DEBUG3
      cout<<"Number of tracks="<<nTracks<<endl;
      for(int itrk=0;itrk<nTracks;itrk++){
      cout<<"Sample="<<Sample[itrk]<<", MuCL="<<TMath::Log(CLMuon_Plan[itrk])<<", iron distance="<<IronDistance[itrk]<<endl;
      }
#endif

      ////////////////////DETERMINE MUON TRACK/////////////////////////////////
      for(int itrk=0;itrk<nTracks;itrk++){
	//if(CLMuon_Plan[itrk]==0) CLMuon_Plan[itrk]=1e-30;
	if(IsReconstructed[itrk]){
	
#ifdef MVAREADING
	  //MVA	  	  
	  
	  SampleTrackMVA=(float) Sample[itrk];
	  for(int ihit=FirstHit;ihit<LimitHits;ihit++){
	    EnergyDepositionMVA[ihit]=EnergyDeposition[itrk][ihit];
	  }

	  MVAdiscriminant[itrk] = tmvareader.EvaluateMVA("BDT method");
	  MVAdiscriminant2[itrk] = tmvareader2.EvaluateMVA("BDT method");

	  LMVA[itrk] = TMath::Log((MVAdiscriminant[itrk]+0.5)/(MVAdiscriminant2[itrk]+0.5));
	  //double LRatio = (MVAdiscriminant[itrk]+0.5)/(MVAdiscriminant2[itrk]+0.5);//since both MVA spread from -0.5 to +0.5 -> Bring them back to 1
	  //double LMVA = LRatio / (1 + LRatio);//In order to tend to 1 when likelihood ratio LRatio -> inf, and to 0 when LRatio -> 0 
	  //double LMVA = TMath::Log(LRatio);//In order to tend to 1 when likelihood ratio LRatio -> inf, and to 0 when LRatio -> 0 
	  //cout<<"Particle type = "<<TypeOfTrack[itrk]<<", sample = "<<Sample[itrk]<<", Likelihood MVA ="<<LMVA<<", LRatio = "<<LRatio<<", mu MVA = "<<MVAdiscriminant[itrk]<<", p MVA = "<<MVAdiscriminant2[itrk]<<endl;
	    
	  if(TypeOfTrack[itrk]==13) MuonTrueMVA=itrk;
	  if(MVAdiscriminant[itrk]>=MVAdiscriminant[MuonRec]) MuonRec=itrk;
	  MuonFound=true;

	  if(LMVA[itrk]>LargestLMVA) LargestLMVA = LMVA[itrk];
	  if(LMVA[itrk]<SmallestLMVA) SmallestLMVA = LMVA[itrk];
	  
	  if(nTracks==1){
	    if(MVAdiscriminant[itrk]>LargestMVAdiscriminant_1track) LargestMVAdiscriminant_1track=MVAdiscriminant[itrk];
	    if(MVAdiscriminant2[itrk]>LargestMVAdiscriminant2_1track) LargestMVAdiscriminant2_1track=MVAdiscriminant2[itrk];
	  }
	  else if(nTracks==2){
	    if(MVAdiscriminant[itrk]>LargestMVAdiscriminant_2tracks) LargestMVAdiscriminant_2tracks=MVAdiscriminant[itrk];
	    if(MVAdiscriminant2[itrk]>LargestMVAdiscriminant2_2tracks) LargestMVAdiscriminant2_2tracks=MVAdiscriminant2[itrk];
	  }
	
	  double eqdist=IronDistance[itrk]+(PlasticDistance[itrk]/IronCarbonRatio);
	  int SimplifiedPDG=DetermineSimplifiedPDG(TypeOfTrack[itrk]);
	  hLMVA_TrueParticle[SimplifiedPDG][Sample[itrk]]->Fill(LMVA[itrk],weight);
	  if((Sample[itrk]==MuonSample1) || (Sample[itrk]==MuonSample2) ||(Sample[itrk]==MuonSample3)){
	    hMVAMuondiscriminantVSDistance_TrueParticle[SimplifiedPDG]->Fill(eqdist,MVAdiscriminant[itrk],weight);
	  }

	  if(Sample[itrk] == 2){
	    if(eqdist >= 10){
	      if(TMath::Abs(TypeOfTrack[itrk])==13) hMVA_MuonVsProtonDiscriminant_TrueMuon->Fill(MVAdiscriminant[itrk],MVAdiscriminant2[itrk],weight);
	      else if(TMath::Abs(TypeOfTrack[itrk])==211) hMVA_MuonVsProtonDiscriminant_TruePion->Fill(MVAdiscriminant[itrk],MVAdiscriminant2[itrk],weight);
	      else if(TMath::Abs(TypeOfTrack[itrk])==2212) hMVA_MuonVsProtonDiscriminant_TrueProton->Fill(MVAdiscriminant[itrk],MVAdiscriminant2[itrk],weight);
	    }
	  }
#else
	  
	  if(ChosenCL[itrk]!=ChosenCL[itrk] || (_isPM && ChosenCL[itrk]==-1)){ Trash=true; cout<<"problem, event trashed"<<endl;continue;}
	  if(TypeOfTrack[itrk]==13) MuonTrue=itrk;
	  if(ChosenCL[itrk]>=ChosenCL[MuonRec]) {PionRec=MuonRec; MuonRec=itrk;}
	  else if(ChosenCL[itrk]>=ChosenCL[PionRec]) PionRec=itrk;
	  MuonFound=true;
	  if(ChosenCL[itrk]<ChosenCL[LowestMuCL]) LowestMuCL=itrk;

#endif
	  if(Plots){
	    if(SelectionFV && !IsData){
	      if(TMath::Abs(TypeOfTrack[itrk])==13) hMuCL_TrueMuon->Fill(ChosenCL[itrk],weight);
	      else if(TMath::Abs(TypeOfTrack[itrk])==211) hMuCL_TruePion->Fill(ChosenCL[itrk],weight);
	      else if(TMath::Abs(TypeOfTrack[itrk])==2212) hMuCL_TrueProton->Fill(ChosenCL[itrk],weight);
	      else hMuCL_TrueOthers->Fill(ChosenCL[itrk],weight);
		//cout<<"hello, particle is="<<TMath::Abs(TypeOfTrack[irec][itrk])<<endl;
	    }
	  }

	}

	if(Plots){
	  //if(SelectionOV[irec]){//MC case: muon at 99%. If data, one should only take long tracks since gamma contamination
	    
	  if(SelectionFV && !IsData){

	    double eqdist=IronDistance[itrk]+(PlasticDistance[itrk]/IronCarbonRatio);
	      
	    if(IsFV && TMath::Abs(TypeOfTrack[itrk])==13 && Sample[itrk]==3) { //stopping muons
	      Pmu_vs_IronDist->Fill(eqdist,Momentum[itrk],weight);	  		   
	    }
		  
	    if(IsFV && TMath::Abs(TypeOfTrack[itrk])==2212 && Sample[itrk]==3)//stopping protons
	      Pp_vs_IronDist->Fill(eqdist,Momentum[itrk],weight);	  		   
		  
	    if(IsFV && TMath::Abs(TypeOfTrack[itrk])==211 && Sample[itrk]==3)//stopping pions
	      Ppi_vs_IronDist->Fill(eqdist,Momentum[itrk],weight);	  		   
		  
	      
	    double distratio = 0.01;//not 0, to avoid being just limits between two bins.
#ifdef DEBUGMVA
	    cout<<"Type of particle="<<TypeOfTrack[itrk]<<", sample="<<Sample[itrk]<<", equivalent iron distance="<<eqdist<<", distance in CH/Fe="<<PlasticDistance[itrk]<<"/"<<IronDistance[itrk]<<endl;
#endif
	    for(int ihit=0;ihit<LimitHits;ihit++){
#ifdef DEBUGMVA
	      if(ihit==0) cout<<setprecision(3)<<"dE/dx = [ ";
	      if(ihit<LimitHits-1) cout<<EnergyDeposition[itrk][ihit]<<", ";
	      else cout<<EnergyDeposition[itrk][ihit]<<" ]"<<endl;
#endif
	      if(TMath::Abs(TypeOfTrack[itrk])==13){
		if(EnergyDeposition[itrk][ihit]!=0) EnergyDepositionLength_Muon[Sample[itrk]]->Fill(eqdist,distratio,EnergyDeposition[itrk][ihit]);
	      }
	      else if(TMath::Abs(TypeOfTrack[itrk])==211){
		if(EnergyDeposition[itrk][ihit]!=0) EnergyDepositionLength_Pion[Sample[itrk]]->Fill(eqdist,distratio,EnergyDeposition[itrk][ihit]);
	      }
	      else if(TMath::Abs(TypeOfTrack[itrk])==2212){
		if(EnergyDeposition[itrk][ihit]!=0) EnergyDepositionLength_Proton[Sample[itrk]]->Fill(eqdist,distratio,EnergyDeposition[itrk][ihit]);
	      }
	      distratio += 1./LimitHits;
	      //cout<<"distratio="<<distratio<<", LmitHits="<<LimitHits<<endl;
	    }

	  }
	}
      }
  
	  
      if(Selection==2){ // CC1pi
	if(TMath::Abs(TypeOfTrack[MuonRec])==13) MuonRec_TruePDG_val[0]+=weight;
	else if(TMath::Abs(TypeOfTrack[MuonRec])==211) MuonRec_TruePDG_val[1]+=weight;
	else if(TMath::Abs(TypeOfTrack[MuonRec])==2212) MuonRec_TruePDG_val[2]+=weight;
	else MuonRec_TruePDG_val[3]+=weight;
	    
	if(Sample[MuonRec]<Sample[PionRec]){
	  // switch PionRec && MuonRec
	  int pion_tmp=MuonRec;
	  int mu_tmp=PionRec;

	  if(TMath::Abs(TypeOfTrack[mu_tmp])==13) MuonRec_TruePDG_switch_val[0]+=weight;
	  else if(TMath::Abs(TypeOfTrack[mu_tmp])==211) MuonRec_TruePDG_switch_val[1]+=weight;
	  else if(TMath::Abs(TypeOfTrack[mu_tmp])==2212) MuonRec_TruePDG_switch_val[2]+=weight;
	  else MuonRec_TruePDG_switch_val[3]+=weight;

	  if(false){
	    // ML 2017-02-03 I remove it because it doesn't improve anything and add some phase space limitation
	    PionRec=pion_tmp;
	    MuonRec=mu_tmp;
	  }
	}
      }
          

      if(!MuonFound){
	NEventsLost+=weight;
	continue;
      }

      
      int MuonLike=0;int ProtonLike=0;int Undetermined=0;
      //cout<<"new"<<endl;
      for(int itrk=0;itrk<nTracks;itrk++){
	if(!IsReconstructed[itrk]) continue;
#ifdef MVAREADING
	double eqdist=IronDistance[itrk]+(PlasticDistance[itrk]/IronCarbonRatio);
	double CutValue=MuonMVACut;
	//if(IsData) cout<<"Track #"<<itrk<<", eqdist="<<eqdist<<", sample="<<Sample[itrk]<<", muon mva="<<MVAdiscriminant[itrk]<<", p mva="<<MVAdiscriminant2[itrk]<<endl;
#ifdef MVA2DCut
	if(eqdist<Limit2DCut) CutValue=MuonMVACut_slope*eqdist+MuonMVACut_origin;
#endif
	//TEMP//
	if(Sample[itrk]==2) CutValue=-0.05;
	
	//if(MVAdiscriminant[itrk]>MuonMVACut) MuonLike++;
	if(MVAdiscriminant[itrk]>CutValue){
	  if(MVAdiscriminant2[itrk]<=ProtonMVACut) MuonLike++;
	  else Undetermined++;
	}
	else{
	  if(MVAdiscriminant2[itrk]>ProtonMVACut) ProtonLike++;
	  else Undetermined++;
	}
	// for WM, mucl=-1 corresponds to undetremined tracks (to few isohits)
	//cout<<"This track has mumva="<<MVAdiscriminant[itrk]<<", pmva="<<MVAdiscriminant2[itrk]<<", number of mu-like="<<MuonLike<<", "<<ProtonLike<<", "<<Undetermined<<endl;

#else	    
	if(ChosenCL[itrk]>MuonCut) MuonLike++;
	else if(ChosenCL[itrk]<ProtonCut && (_isPM? true: ChosenCL[itrk]>=0)) ProtonLike++;
	// for WM, mucl=-1 corresponds to undetremined tracks (to few isohits)
	else Undetermined++;
#endif
      }
      //double rnd=rand->Uniform(0,1);
            
      //if((nTracks==1 && ChosenCL[0]>-3) || ((nTracks==2 && TMath::Max(ChosenCL[0],ChosenCL[1])>-3 && TMath::Min(ChosenCL[0],ChosenCL[1])<-3) && (TMath::Min(ChosenCL[0],ChosenCL[1])>=0)) || ((nTracks==3 && ChosenCL[MuonRec]>-3 && ChosenCL[(MuonRec+1)%3]<-3 && ChosenCL[(MuonRec+2)%3]<-3)) || ((nTracks==4 && ChosenCL[MuonRec]>-3 && ChosenCL[(MuonRec+1)%4]<-3 && ChosenCL[(MuonRec+2)%4]<-3 && ChosenCL[(MuonRec+3)%4]<-3))){
      
      //if(TrackWidth[MuonRec]>1.1 && IronDistance[MuonRec]<15) continue;
      //else if(TrackWidth[MuonRec]>1.2 && IronDistance[MuonRec]<20) continue;
      //else if(TrackWidth[MuonRec]>1.3 && IronDistance[MuonRec]<30) continue;
      //else if(TrackWidth[MuonRec]>1.5) continue;
      //cout<<weight<<endl;
      int BinRecMom=0;
      int BinRecAngle=0;
      bool old=true;
      double MuonSample=Sample[MuonRec];//Used only in the line just below to separate through-going vs stopping Reconstructed Bin
      double EquivalentIronDistance=IronDistance[MuonRec]+(PlasticDistance[MuonRec]/IronCarbonRatio);
      if(MuonSample>=4) EquivalentIronDistance = BinningRecMom[NBinsRecMom-1] + (BinningRecMom[NBinsRecMom] - BinningRecMom[NBinsRecMom-1])/2.;
      else if (EquivalentIronDistance>=BinningRecMom[NBinsRecMom-1]){
	EquivalentIronDistance = BinningRecMom[NBinsRecMom-2];//To create a trash bins for stopped sampl;e
      cout<<"A trash bin for the stopped sample has been created at "<<BinningRecMom[NBinsRecMom-2]<<" cm"<<endl;
      } 
      for(int i=0;i<NBinsRecMom;i++){
	if((EquivalentIronDistance)<BinningRecMom[i+1]){BinRecMom=i;break;}
      }
      for(int i=0;i<NBinsRecAngle;i++){
	if(TrackAngle[MuonRec]<BinningRecAngle[i+1]){BinRecAngle=i;break;}
      }
      
      if(EquivalentIronDistance >= BinningRecMom[NBinsRecMom]) cout<<"You are loosing some events due to limits in your reconstruction binning. Please check. value = "<<EquivalentIronDistance<<" vs "<<BinningRecMom[NBinsRecMom]<<endl;
      
      
      ///////////////////THE ACTUAL SELECTION///////////////////////////////////////
      if(MuonTrue==MuonRec) MuonID->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
      if(MuonTrueMVA==MuonRec) MuonIDMVA->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
      MuonIDTotal->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
      //if(Sample[MuonTrueMVA] == 2){
      //if(MuonTrueMVA==MuonRec) pMisID->Fill(0,,weight);
      //}
      
      int nRecTracks=nTracks-nBadRecTracks;

    
      if((Sample[MuonRec]==MuonSample1) || (Sample[MuonRec]==MuonSample2) || (Sample[MuonRec]==MuonSample3)){
	if(Plots){
	  // for CC1pi MuonSample2 is set to 5
	  if(nRecTracks==1){
	    hMuCL[FSIInt]->Fill(ChosenCL[MuonRec],weight);
	    int SimplifiedPDG=DetermineSimplifiedPDG(TypeOfTrack[0]);
	    TrueParticleType_1track[FSIInt]->Fill(SimplifiedPDG);
	  }
	  if(nRecTracks==2){
	    if(MuonRec==0) LowestMuCL=1;
	    else LowestMuCL=0;
	    hMuCL_2tracks[FSIInt]->Fill(ChosenCL[MuonRec],ChosenCL[LowestMuCL],weight);
	    int SimplifiedPDG[2]={DetermineSimplifiedPDG(TypeOfTrack[0]),DetermineSimplifiedPDG(TypeOfTrack[1])};
	    TrueParticleType_2tracks[FSIInt]->Fill(SimplifiedPDG[0],SimplifiedPDG[1]);
	    if(MuonLike==1 && Undetermined==0){
	      hMuCL_Lowest[FSIInt]->Fill(ChosenCL[LowestMuCL],weight);
	    }
	  }

#ifdef MVAREADING
	  //if(FSIInt==3) cout<<"FSI type="<<FSIInt<<", Number of mu-like / p-like / undet = "<<MuonLike<<" / "<<ProtonLike<<" / "<<Undetermined<<endl;
	
	  //MVA
	  if(Plots){
	    hLMVA_LargestLMVA[FSIInt]->Fill(LargestLMVA,weight);
	    if(nTracks>1) hLMVA_SmallestLMVA[FSIInt]->Fill(SmallestLMVA,weight);
	    if(nTracks==1){
	      int SimplifiedPDG=DetermineSimplifiedPDG(TypeOfTrack[0]);
	      hMVAMuondiscriminantVSPDG_1track[FSIInt]->Fill(SimplifiedPDG,LargestMVAdiscriminant_1track,weight);
	      hMVAMuondiscriminant_1track[FSIInt]->Fill(LargestMVAdiscriminant_1track,weight);
	      hMVAMuondiscriminantVSMuonMomentum_1track[FSIInt]->Fill(TrueMomentumMuon,LargestMVAdiscriminant_1track,weight);
	      hMVAMuondiscriminantVSDistance_1track[FSIInt]->Fill(EquivalentIronDistance,LargestMVAdiscriminant_1track,weight);
	      hLMVA_1track[FSIInt]->Fill(LargestLMVA,weight);
	    }
	    if(nTracks==2){
	      hMVAProtondiscriminant_2tracks[FSIInt]->Fill(LargestMVAdiscriminant2_2tracks,weight);
	      hMVAMuonVSProtondiscriminant_2tracks[FSIInt]->Fill(LargestMVAdiscriminant_2tracks,LargestMVAdiscriminant2_2tracks,weight);
	      //Small search of the proton track: 2 track sample so proton track is the non muon track.
	      int ProtonRec=-1;//double EquivalentIronDistance_Proton=0;
	      if(MuonRec==0) ProtonRec=1;
	      else ProtonRec=0;
	      //EquivalentIronDistance_Proton=IronDistance[ProtonRec]+(PlasticDistance[ProtonRec]/IronCarbonRatio);
	      //
	      hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[FSIInt]->Fill(PlasticDistance[ProtonRec],LargestMVAdiscriminant2_2tracks,weight);//The proton is not required to stop, so checking is distance is less meaningful, since it can escape.
	      hLMVA_2tracks[FSIInt]->Fill(LargestLMVA,SmallestLMVA,weight);
	    }
	  }
#endif
	}
	hNTracks[FSIInt]->Fill(nRecTracks,weight);
	hRecMom[FSIInt]->Fill(EquivalentIronDistance,weight);
	hRecAngle[FSIInt]->Fill(TrackAngle[MuonRec],weight);
      }
    
      if(Selection == 1 && FSIInt<3 && IsNuMu){
	MCReconstructionEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	//if((Sample[MuonTrueMVA]==MuonSample1)||(Sample[MuonTrueMVA]==MuonSample2)) MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	if(Sample[MuonTrueMVA]>=2 && Sample[MuonTrueMVA]!=4) MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
      }
      if(Selection == 2 && FSIInt==3 && IsNuMu){
	MCReconstructionEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	//if((Sample[MuonTrueMVA]==MuonSample1)||(Sample[MuonTrueMVA]==MuonSample2)) MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	if(Sample[MuonTrueMVA]>=2 && Sample[MuonTrueMVA]!=4) MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	
      }






      
      ///////////////////////////////////////////
      if(MuonLike==1 && Undetermined==0 /*&& nRecTracks<=2*/){//CC0pi
	//cout<<"Interaction value="<<FSIInt<<", Ion distance="<<IronDistance[MuonRec]+(PlasticDistance[MuonRec]/IronCarbonRatio)<<endl;
	//cout<<"passed the cut"<<endl;
	//TEMPORARY
	nEvents[3]+=weight;
	nEventsInter[3][FSIInt]+=weight;

	if(nRecTracks==2) hSampleSecondTrack[FSIInt]->Fill(Sample[LowestMuCL]);
	//
	/*		
	if(Selection==1){

	  DataSelected[BinRecMom][BinRecAngle]+=weight;
	  //cout<<"Bin="<<BinRecMom<<","<<BinRecAngle<<", data="<<DataSelected[BinRecMom][BinRecAngle]<<endl;
	  if(!IsData){
	    if(FSIInt<3 && IsNuMu){
	      if(Plots) MCEfficiency_Energy->Fill(Enu,weight);
	      MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	      Efficiency[BinTrueMom][BinTrueAngle]+=weight;
#ifdef DEBUG
	      if(EquivalentIronDistance<40 && TrueMomentumMuon>1){
		cout<<"evt number="<<ievt<<", number of tracks="<<nRecTracks<<endl;
		for(int itrk=0;itrk<nRecTracks;itrk++){
		  cout<<"CL="<<ChosenCL[itrk]<<", distance="<<EquivalentIronDistance<<", sample="<<Sample[itrk]<<endl;
		}
	      }
#endif
	    }
	    else{
	      BkgSelected[BinRecMom][BinRecAngle]+=weight;
	    }
	  }
	  hRecMom_CC0pi[FSIInt]->Fill(EquivalentIronDistance,weight);
	  hRecAngle_CC0pi[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	}
	*/

	if(nRecTracks>2) continue;
	if(Sample[MuonRec]>2){//only ingrid tracks
	  hRecMom_CC0pi_full[FSIInt]->Fill(EquivalentIronDistance,weight);
	  hRecAngle_CC0pi_full[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	  if(Sample[MuonRec]==3) MuCL_vs_IronDist->Fill(EquivalentIronDistance,ChosenCL[MuonRec],weight);
	}
	      
	if(Selection==1){
	  nEvents[4]+=weight;
	  nEventsInter[4][FSIInt]+=weight;
	  DataSelected_full[BinRecMom][BinRecAngle]+=weight;
	  if(!IsData){
	    if(!BkgSub) MCSelected_full[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;//New: likelihood without background substraction case.
	    else if(FSIInt<3 && IsNuMu) MCSelected_full[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	    
	    if(FSIInt<3 && IsNuMu){
	      Efficiency_full[BinTrueMom][BinTrueAngle]+=weight;
	    }
	    else{
	      BkgSelected_full[BinRecMom][BinRecAngle]+=weight;
	    }
	  }
	  if((Sample[MuonRec]!=MuonSample1)&&(Sample[MuonRec]!=MuonSample2)&&(Sample[MuonRec]!=MuonSample3)) continue;
	  nEvents[5]+=weight;
	  nEventsInter[5][FSIInt]+=weight;
	  hRecMom_CC0pi[FSIInt]->Fill(EquivalentIronDistance,weight);
	  hRecAngle_CC0pi[FSIInt]->Fill(TrackAngle[MuonRec],weight);
				
	  DataSelected[BinRecMom][BinRecAngle]+=weight;
	  hNTracks_CC0pi[FSIInt]->Fill(nRecTracks,weight);

	  if(!IsData){
	    //cout<<FSIInt<<", "<<Num_Int<<", bin true mom="<<BinTrueMom<<", bin true angle="<<BinTrueAngle<<endl;
	    if(!BkgSub){MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;/*cout<<"filled"<<endl;*/}//New: likelihood without background substraction case.
	    else if(FSIInt<3 && IsNuMu) MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	    if(FSIInt<3 && IsNuMu){
	      if(Plots){
		MCEfficiency_Energy->Fill(Enu,weight);
		MCEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	      }
	      Efficiency[BinTrueMom][BinTrueAngle]+=weight;
	    }
	    else{
	      BkgSelected[BinRecMom][BinRecAngle]+=weight;
	    }

	    hRecCC0pi_MomentumResolution->Fill(EquivalentIronDistance,TrueMomentumMuon,weight);
	    hRecCC0pi_AngleResolution->Fill(TrackAngle[MuonRec],TrueAngleMuon,weight);

	  }
	  if(Sample[MuonRec]==MuonSample1 || Sample[MuonRec]==MuonSample2){
	    nEvents[6]+=weight;
	    nEventsInter[6][FSIInt]+=weight;
	    double leq=EquivalentIronDistance;
	    hRecMom_CC0pi_restr[FSIInt]->Fill(leq,weight);
	    hRecAngle_CC0pi_restr[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	  }
	}

      }
      else if((MuonLike==2 && Undetermined==0) || (MuonLike==1 && Undetermined==1)){//Side band CC1pi - for WM I reject Undertermined tracks (mucl=-1)
	// for PM, no requirement on Undetermined==0 was done -> I make it now
	
	nEvents[3]+=weight;
	nEventsInter[3][FSIInt]+=weight;
	      
	if(nRecTracks>3) continue;
	if(Sample[MuonRec]>2){//only ingrid tracks
	  hRecMom_CC1pi_full[FSIInt]->Fill(EquivalentIronDistance,weight);
	  hRecAngle_CC1pi_full[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	  if(Sample[MuonRec]==3) MuCL_vs_IronDist->Fill(EquivalentIronDistance,ChosenCL[MuonRec],weight);
	}
	      
	if(Selection==2){
	  nEvents[4]+=weight;
	  nEventsInter[4][FSIInt]+=weight;
	  DataSelected_full[BinRecMom][BinRecAngle]+=weight;
	  
	  if(!IsData){
	    if(!BkgSub) MCSelected_full[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;//New: likelihood without background substraction case.
	    else if(FSIInt==3 && IsNuMu) MCSelected_full[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	    if(FSIInt==3 && IsNuMu){
	      //if(Plots) MCEfficiency_Energy->Fill(Enu,weight);
	      //MCSelected_full[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	      Efficiency_full[BinTrueMom][BinTrueAngle]+=weight;
	    }
	    else{
	      BkgSelected_full[BinRecMom][BinRecAngle]+=weight;
	    }
	  }
		

	  if((Sample[MuonRec]!=MuonSample1)&&(Sample[MuonRec]!=MuonSample2)&&(Sample[MuonRec]!=MuonSample3)) continue;
	  nEvents[5]+=weight;
	  nEventsInter[5][FSIInt]+=weight;
	  hRecMom_CC1pi[FSIInt]->Fill(EquivalentIronDistance,weight);
	  hRecAngle_CC1pi[FSIInt]->Fill(TrackAngle[MuonRec],weight);
				
	  DataSelected[BinRecMom][BinRecAngle]+=weight;
	  if(!IsData){
	    if(!BkgSub) MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;//New: likelihood without background substraction case.
	    else if(FSIInt==3 && IsNuMu) MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	    if(FSIInt==3 && IsNuMu){
	      if(Plots){
		MCEfficiency_Energy->Fill(Enu,weight);
		MCEfficiency_PTheta_ThinBinning->Fill(TrueMomentumMuon,TrueAngleMuon,weight);
	      }
	      //MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
	      Efficiency[BinTrueMom][BinTrueAngle]+=weight;
	    }
	    else{
	      BkgSelected[BinRecMom][BinRecAngle]+=weight;
	    }
	  }
	  if(Sample[MuonRec]==MuonSample1 || Sample[MuonRec]==MuonSample2){
	    nEvents[6]+=weight;
	    nEventsInter[6][FSIInt]+=weight;
	    double leq=EquivalentIronDistance;
	    hRecMom_CC1pi_restr[FSIInt]->Fill(leq,weight);
	    hRecAngle_CC1pi_restr[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	  }
	}
	else if(Selection==1){//Side band case
	  /*
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
	  */	

	  if((Sample[MuonRec]!=MuonSample1)&&(Sample[MuonRec]!=MuonSample2)&&(Sample[MuonRec]!=MuonSample3)) continue;
	  //nEvents[5]+=weight;
	  //nEventsInter[5][FSIInt]+=weight;
	  hRecMom_CC1pi[FSIInt]->Fill(EquivalentIronDistance,weight);
	  hRecAngle_CC1pi[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	  /*			
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
	  if(Sample[MuonRec]==MuonSample1){
	    nEvents[6]+=weight;
	    nEventsInter[6][FSIInt]+=weight;
	    double leq=EquivalentIronDistance;
	    hRecMom_CC1pi_restr[FSIInt]->Fill(leq,weight);
	    hRecAngle_CC1pi_restr[FSIInt]->Fill(TrackAngle[MuonRec],weight);
	    }*/
	}

      }
      else if(MuonLike>=3){//Side band CCNpi
	hRecMom_CCNpi[FSIInt]->Fill(EquivalentIronDistance,weight);
	hRecAngle_CCNpi[FSIInt]->Fill(TrackAngle[MuonRec],weight);
      }	     
    }
  }
      	    
#ifdef DEBUG2
  cout<<"End of event loop"<<endl;
#endif

  if(Plots){//Should be done before Efficiency[c0][c1] becomes the real efficiency. Should be done when it is still the number of selected events.
    if(!IsData ){
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	double TotalSelectedMomentum=0,TotalTrue0Momentum=0;
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  if(Selection==1) TotalTrue0Momentum += TotalCC0piEvent[c0][c1];
	  else if(Selection==2) TotalTrue0Momentum += TotalCC1piEvent[c0][c1];
	  TotalSelectedMomentum += Efficiency[c0][c1];
	}
	double EfficiencyMomentum = TotalSelectedMomentum;
	if(TotalTrue0Momentum!=0) EfficiencyMomentum /= TotalTrue0Momentum;
	if(Plots) MCEfficiency_PTheta_Momentum->SetBinContent(c0+1,EfficiencyMomentum);
      }

      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	double TotalSelectedAngle=0,TotalTrue0Angle=0;
	for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	  if(Selection==1) TotalTrue0Angle += TotalCC0piEvent[c0][c1];
	  else if(Selection==2) TotalTrue0Angle += TotalCC1piEvent[c0][c1];
	  TotalSelectedAngle += Efficiency[c0][c1];
	}
	double EfficiencyAngle = TotalSelectedAngle;
	if(TotalTrue0Angle!=0) EfficiencyAngle /= TotalTrue0Angle;
	if(Plots) MCEfficiency_PTheta_Angle->SetBinContent(c1+1,EfficiencyAngle);
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
	if(Plots){
	  MCEfficiency->SetBinContent(c0+1,c1+1,Efficiency[c0][c1]);
	  MCEfficiency_PTheta->SetBinContent(c0+1,c1+1,Efficiency[c0][c1]);
	}
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

  for(int i=0;i<nCuts;i++){
    for(int j=0;j<NFSIs;j++) CutNEvents[j]->SetBinContent(i+1,nEventsInter[i][j]);
  }
  if(Selection==2){
    for(int i=0;i<nCuts;i++){
      cout<<myCuts[i]<<": \tselected="<<nEvents[i]<<" \tpurity="<<100.*nEventsInter[i][3]/nEvents[i]<<" \tefficiency="<<nEventsInter[i][3]/TotalTrue0*100<<endl;
      CutEfficiency->SetBinContent(i+1,nEventsInter[i][3]/TotalTrue0*100);
      CutPurity->SetBinContent(i+1,nEventsInter[i][3]/nEvents[i]*100);

      cout<<"\\hline"<<endl;
      cout << " & Total & Purity & Efficiency & CC0$\\pi$ & CC1$\\pi^{\\pm}$ & CC$\\pi^{0}$ & CCothers & NC \\\\"<<endl<<"\\hline"<<endl;
      cout<<setprecision(0)<<std::fixed<<myCuts[i]<<" & "<<nEvents[i]<<" & "<< 100*(nEventsInter[i][3])/nEvents[i]<<" $\\%$ & "<< 100*(nEventsInter[i][3])/TotalTrue0 <<" $\\%$ & "<<(nEventsInter[i][0]+nEventsInter[i][1]+nEventsInter[i][2])<<" & "<<nEventsInter[i][3]<<" & "<<nEventsInter[i][4]<<" & "<<nEventsInter[i][5]<<" & "<<nEventsInter[i][6]<<" \\\\"<<endl; 
    }
  }
  else if(Selection==1){
    cout<<"\\hline"<<endl;
    cout << " & Total & Purity & Efficiency & CC0$\\pi$ & CC1$\\pi^{\\pm}$ & CC$\\pi^{0}$ & CCothers & NC \\\\"<<endl<<"\\hline"<<endl;
    for(int i=0;i<nCuts;i++){
      cout<<setprecision(0)<<std::fixed<<myCuts[i]<<" & "<<nEvents[i]<<" & "<< 100*(nEventsInter[i][0]+nEventsInter[i][1]+nEventsInter[i][2])/nEvents[i]<<" $\\%$ & "<< 100*(nEventsInter[i][0]+nEventsInter[i][1]+nEventsInter[i][2])/TotalTrue0 <<" $\\%$ & "<<(nEventsInter[i][0]+nEventsInter[i][1]+nEventsInter[i][2])<<" & "<<nEventsInter[i][3]<<" & "<<nEventsInter[i][4]<<" & "<<nEventsInter[i][5]<<" & "<<nEventsInter[i][6]<<" \\\\"<<endl;
      CutEfficiency->SetBinContent(i+1,(nEventsInter[i][0]+nEventsInter[i][1]+nEventsInter[i][2])/TotalTrue0);
      CutPurity->SetBinContent(i+1,100*(nEventsInter[i][0]+nEventsInter[i][1]+nEventsInter[i][2])/nEvents[i]);
    }
    cout<<"\\hline"<<endl;
    cout<<endl<<endl;

    
    cout<<"Side bande CC1pi:"<<endl;
    double TotalSideBand=0;
    for(int ifsi=0;ifsi<NFSIs;ifsi++) TotalSideBand+=hRecMom_CC1pi[ifsi]->Integral();
    cout<<"\\hline"<<endl;
    cout << " & CC0$\\pi$ & CC1$\\pi^{\\pm}$ & CC$\\pi^{0}$ & CCothers & NC \\\\"<<endl<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;
    cout<<"Nevents & "<<(hRecMom_CC1pi[0]->Integral()+hRecMom_CC1pi[1]->Integral()+hRecMom_CC1pi[2]->Integral())<<" & "<<hRecMom_CC1pi[3]->Integral() << " & "<< hRecMom_CC1pi[4]->Integral() <<" & "<< hRecMom_CC1pi[5]->Integral()<<" & "<< hRecMom_CC1pi[6]->Integral() <<"\\\\"<<endl; 
    cout<<"Nevents & "<<100*(hRecMom_CC1pi[0]->Integral()+hRecMom_CC1pi[1]->Integral()+hRecMom_CC1pi[2]->Integral())/TotalSideBand<<" $\\%$ & "<<100*hRecMom_CC1pi[3]->Integral()/TotalSideBand << " $\\%$ & "<< 100*hRecMom_CC1pi[4]->Integral()/TotalSideBand <<" $\\%$ & "<< 100*hRecMom_CC1pi[5]->Integral()/TotalSideBand<<" $\\%$ & "<< 100*hRecMom_CC1pi[6]->Integral()/TotalSideBand << " $\\%$ \\\\"<<endl; 
    cout<<"\\hline"<<endl;

    //
    double SelectedCC0pi_1track=0;
    double SelectedCC0pi_2tracks=0;
    double SelectedAll_1track=0;
    double SelectedAll_2tracks=0;
    
    double TotalCC0pi_1track=0;
    double TotalCC0pi_2tracks=0;
    double TotalAll_1track=0;
    double TotalAll_2tracks=0;
    for(int ifsi=0;ifsi<NFSIs;ifsi++){
      if(ifsi<3){
	SelectedCC0pi_1track += hNTracks_CC0pi[ifsi]->GetBinContent(2);
	SelectedCC0pi_2tracks += hNTracks_CC0pi[ifsi]->GetBinContent(3);
	TotalCC0pi_1track += hNTracks[ifsi]->GetBinContent(2);
	TotalCC0pi_2tracks += hNTracks[ifsi]->GetBinContent(3);
      }
      SelectedAll_1track += hNTracks_CC0pi[ifsi]->GetBinContent(2);
      SelectedAll_2tracks += hNTracks_CC0pi[ifsi]->GetBinContent(3);    
      TotalAll_1track += hNTracks[ifsi]->GetBinContent(2);
      TotalAll_2tracks += hNTracks[ifsi]->GetBinContent(3);    
    }
    cout<<"Check of the CC0pi Selection:"<<endl;
    cout<<"\\hline"<<endl;
    cout<<"Original Purity & Efficiency & Purity \\\\"<<endl;
    cout<<100*TotalCC0pi_1track/TotalAll_1track<<" & "<<100*SelectedCC0pi_1track/TotalCC0pi_1track<<" & "<<100*SelectedCC0pi_1track/SelectedAll_1track<<" \\\\ "<<endl;
    cout<<100*TotalCC0pi_2tracks/TotalAll_2tracks<<" & "<<100*SelectedCC0pi_2tracks/TotalCC0pi_2tracks<<" & "<<100*SelectedCC0pi_2tracks/SelectedAll_2tracks<<" \\\\ "<<endl;
    cout<<"\\hline"<<endl;
  }
  
  if(Plots){
    
    //Efficiency
    MCEfficiency_Momentum_ThinBinning = (TH1D*) MCEfficiency_PTheta_ThinBinning->ProjectionX("MCEfficiency_Momentum_ThinBinning",0,MCEfficiency_PTheta_ThinBinning->GetNbinsX());
    MCEfficiency_Angle_ThinBinning = (TH1D*) MCEfficiency_PTheta_ThinBinning->ProjectionY("MCEfficiency_Angle_ThinBinning",0,MCEfficiency_PTheta_ThinBinning->GetNbinsY());

    MCReconstructionEfficiency_Momentum_ThinBinning = (TH1D*) MCReconstructionEfficiency_PTheta_ThinBinning->ProjectionX("MCReconstructionEfficiency_Momentum_ThinBinning",0,MCReconstructionEfficiency_PTheta_ThinBinning->GetNbinsX());
    MCReconstructionEfficiency_Angle_ThinBinning = (TH1D*) MCReconstructionEfficiency_PTheta_ThinBinning->ProjectionY("MCReconstructionEfficiency_Angle_ThinBinning",0,MCReconstructionEfficiency_PTheta_ThinBinning->GetNbinsY());

    MCReconstructionAndTopologyEfficiency_Momentum_ThinBinning = (TH1D*) MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->ProjectionX("MCReconstructionAndTopologyEfficiency_Momentum_ThinBinning",0,MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->GetNbinsX());
    MCReconstructionAndTopologyEfficiency_Angle_ThinBinning = (TH1D*) MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->ProjectionY("MCReconstructionAndTopologyEfficiency_Angle_ThinBinning",0,MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->GetNbinsY());

    
    if(Selection==1){
      MCEfficiency_Energy->Divide(TotalCC0piEvent_Energy);

      TotalEvent_Momentum_ThinBinning = (TH1D*) TotalCC0piEvent_PTheta_ThinBinning->ProjectionX("TotalEvent_Momentum_ThinBinning",0,TotalCC0piEvent_PTheta_ThinBinning->GetNbinsX());
      TotalEvent_Angle_ThinBinning = (TH1D*) TotalCC0piEvent_PTheta_ThinBinning->ProjectionY("TotalEvent_Angle_ThinBinning",0,TotalCC0piEvent_PTheta_ThinBinning->GetNbinsY());

      MCEfficiency_PTheta_ThinBinning->Divide(TotalCC0piEvent_PTheta_ThinBinning);

      MCReconstructionEfficiency_PTheta_ThinBinning->Divide(TotalCC0piEvent_PTheta_ThinBinning);

      MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Divide(TotalCC0piEvent_PTheta_ThinBinning);
    }
    else if(Selection==2){
      MCEfficiency_Energy->Divide(TotalCC1piEvent_Energy);

      TotalEvent_Momentum_ThinBinning = (TH1D*) TotalCC1piEvent_PTheta_ThinBinning->ProjectionX("TotalEvent_Momentum_ThinBinning",0,TotalCC1piEvent_PTheta_ThinBinning->GetNbinsX());
      TotalEvent_Angle_ThinBinning = (TH1D*) TotalCC1piEvent_PTheta_ThinBinning->ProjectionY("TotalEvent_Angle_ThinBinning",0,TotalCC1piEvent_PTheta_ThinBinning->GetNbinsY());

      MCEfficiency_PTheta_ThinBinning->Divide(TotalCC1piEvent_PTheta_ThinBinning);

      MCReconstructionEfficiency_PTheta_ThinBinning->Divide(TotalCC1piEvent_PTheta_ThinBinning);

      MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Divide(TotalCC1piEvent_PTheta_ThinBinning);
    }
    MCEfficiency_Momentum_ThinBinning->Divide(TotalEvent_Momentum_ThinBinning);
    MCEfficiency_Angle_ThinBinning->Divide(TotalEvent_Angle_ThinBinning);

    MCReconstructionEfficiency_Momentum_ThinBinning->Divide(TotalEvent_Momentum_ThinBinning);
    MCReconstructionEfficiency_Angle_ThinBinning->Divide(TotalEvent_Angle_ThinBinning);

    MCReconstructionAndTopologyEfficiency_Momentum_ThinBinning->Divide(TotalEvent_Momentum_ThinBinning);
    MCReconstructionAndTopologyEfficiency_Angle_ThinBinning->Divide(TotalEvent_Angle_ThinBinning);

    /*
    //This is to calculate the relative variation of the efficiency
    MCEfficiency_Momentum_ThinBinning_Relative = (TH1D*) MCEfficiency_Momentum_ThinBinning->Clone("MCEfficiency_Momentum_ThinBinning");
    MCEfficiency_Angle_ThinBinning_Relative = (TH1D*) MCEfficiency_Angle_ThinBinning->Clone("MCEfficiency_Angle_ThinBinning");
    for(int ibinx=1;ibinx<=MCEfficiency_Momentum_ThinBinning->GetNbinsX();ibinx++){
      double previousBin=MCEfficiency_Momentum_ThinBinning->GetBinContent(ibinx-1);
      double currentBin=MCEfficiency_Momentum_ThinBinning->GetBinContent(ibinx);
      double relative=currentBin-previousBin;
      if((currentBin+previousBin) !=0) relative /= (currentBin+previousBin);
      else relative = 0;
      MCEfficiency_Momentum_ThinBinning_Relative->SetBinContent(ibinx,relative);
    }
    for(int ibinx=1;ibinx<=MCEfficiency_Angle_ThinBinning->GetNbinsX();ibinx++){
      double previousBin=MCEfficiency_Angle_ThinBinning->GetBinContent(ibinx-1);
      double currentBin=MCEfficiency_Angle_ThinBinning->GetBinContent(ibinx);
      double relative=currentBin-previousBin;
      if((currentBin+previousBin) !=0) relative /= (currentBin+previousBin);
      else relative = 0;
      MCEfficiency_Angle_ThinBinning_Relative->SetBinContent(ibinx,relative);
    }*/
    //End efficiency

    int colors[4]={kBlue,kGreen+3,kRed,kGray};
    MuonRec_TruePDG=new TPie("MuonRecTruePDG","MuonRec true PDG",4,MuonRec_TruePDG_val,colors);
    MuonRec_TruePDG_switch=new TPie("MuonRecTruePDGswitch","MuonRec true PDG",4,MuonRec_TruePDG_switch_val,colors);
    MuonRec_TruePDG->SetLabelFormat("#splitline{%txt}{(%perc)}");
    MuonRec_TruePDG_switch->SetLabelFormat("#splitline{%txt}{(%perc)}");
  }
    
  char OutNameEvent[256];
  char OutNameEvent_root[256];
  ofstream fEvent;
  sprintf(OutNameEvent,"%s.txt",outnameevent);
  sprintf(OutNameEvent_root,"%s.root",outnameevent);


  if(Systematics_Flux){
    //sprintf(OutNameEvent,"files/MCSelected_Systematics%d_%d.txt",ErrorType,ErrorIteration);
    //sprintf(OutNameEvent_root,"files/MCSelected_Systematics%d_%d.root",ErrorType,ErrorIteration);
    for(int i=1;i<=NBinsEnergyFlux;i++){
      double Nneut=NeutrinoFlux->GetBinContent(i);
      NeutrinoFlux->SetBinContent(i,Nneut+Nneut*FluxVector[i-1]);
    }
  }
  /////
#ifdef DEBUG2
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    cout<<endl<<"Mom bin:"<<c0<<endl;
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      cout<<endl<<"Angle bin:"<<c1<<endl;
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  cout<<MCSelected[c0][c1][e0][e1]<<" ";
	}
      }
    }
  }
#endif
  /////
  
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
    CutEfficiency->Write();
    CutPurity->Write();
    wfile->Close();
#ifdef DEBUG2
  cout<<"Wrote the likelihood matrix in root file"<<endl;
#endif
  }
  delete MCTrueEvents;
  delete CutPurity;
  delete CutEfficiency;
  
  if(Plots){

#ifdef DEBUG2
  cout<<"Plot option activated. Will write several histo in the root file"<<endl;
#endif

    TFile * wfile = new TFile(OutNameEvent_root,"update");
    if(!IsData){
      THStack * Stack_MuCL_Particles = new THStack("Stack_MuCL_Particles","");

      THStack * Stack_MuCL = new THStack("Stack_MuCL","");
      THStack * Stack_MuCL_Lowest = new THStack("Stack_MuCL_Lowest","");
      THStack * Stack_NTracks = new THStack("Stack_NTracks","");
      THStack * Stack_RecMom = new THStack("Stack_RecMom","");
      THStack * Stack_RecAngle = new THStack("Stack_RecAngle","");
      THStack * Stack_FCFVTrueEvents = new THStack("Stack_FCFVTrueEvents","");
      THStack * Stack_NTracks_CC0pi = new THStack("Stack_NTracks_CC0pi","");
      THStack * Stack_RecMom_CC0pi_restr = new THStack("Stack_RecMom_CC0pi_restr","");
      THStack * Stack_RecAngle_CC0pi_restr = new THStack("Stack_RecAngle_CC0pi_restr","");
      THStack * Stack_RecMom_CC0pi = new THStack("Stack_RecMom_CC0pi","");
      THStack * Stack_RecAngle_CC0pi = new THStack("Stack_RecAngle_CC0pi","");
      THStack * Stack_RecMom_CC0pi_full = new THStack("Stack_RecMom_CC0pi_full","");
      THStack * Stack_RecAngle_CC0pi_full = new THStack("Stack_RecAngle_CC0pi_full","");
      THStack * Stack_RecMom_CC1pi_restr = new THStack("Stack_RecMom_CC1pi_restr","");
      THStack * Stack_RecAngle_CC1pi_restr = new THStack("Stack_RecAngle_CC1pi_restr","");
      THStack * Stack_RecMom_CC1pi = new THStack("Stack_RecMom_CC1pi","");
      THStack * Stack_RecAngle_CC1pi = new THStack("Stack_RecAngle_CC1pi","");
      THStack * Stack_RecMom_CC1pi_full = new THStack("Stack_RecMom_CC1pi_full","");
      THStack * Stack_RecAngle_CC1pi_full = new THStack("Stack_RecAngle_CC1pi_full","");
      THStack * Stack_RecMom_CCNpi = new THStack("Stack_RecMom_CCNpi","");
      THStack * Stack_RecAngle_CCNpi = new THStack("Stack_RecAngle_CCNpi","");
      THStack * Stack_SampleSecondTrack = new THStack("Stack_SampleSecondTrack","");
      THStack * Stack_TotalChargeSecondTrack = new THStack("Stack_TotalChargeSecondTrack","");
      THStack * Stack_TotalChargePerDistanceSecondTrack = new THStack("Stack_TotalChargePerDistanceSecondTrack","");
      THStack * Stack_TotalChargePerDistanceFirstTrack = new THStack("Stack_TotalChargePerDistanceFirstTrack","");
      THStack * Stack_OpeningAngle = new THStack("Stack_OpeningAngle","");
      THStack * Stack_CoplanarityAngle = new THStack("Stack_CoplanarityAngle","");
      THStack * Stack_MVAdiscriminant = new THStack("Stack_MVAdiscriminant","");
      THStack * Stack_MVAMuondiscriminant_1track = new THStack("Stack_MVAMuondiscriminant_1track","");
      THStack * Stack_MVAProtondiscriminant_2tracks = new THStack("Stack_MVAProtondiscriminant_2tracks","");
      THStack * Stack_CutNEvents = new THStack("Stack_CutNEvents","");
      
#ifdef DEBUG2
  cout<<"Stack plots declared"<<endl;
#endif
      ProduceStackParticles(hMuCL_TrueMuon,hMuCL_TruePion,hMuCL_TrueProton,hMuCL_TrueOthers,Stack_MuCL_Particles);
  
      ProduceStack(hMuCL,Stack_MuCL);
      ProduceStack(hMuCL_Lowest,Stack_MuCL_Lowest);
      //No stack produced for 2tracks 2D plot, but should still normalise as done in the produce stack!
      //for(int fsi=0;fsi<NFSIs;fsi++) hMuCL_2tracks[fsi]->Scale((5.86e20/1e21)*100/1000);
      
      ProduceStack(hNTracks,Stack_NTracks);
      ProduceStack(hRecMom,Stack_RecMom);
      ProduceStack(hRecAngle,Stack_RecAngle);
      ProduceStack(hFCFVTrueEvents,Stack_FCFVTrueEvents);
      ProduceStack(hRecMom_CC0pi_restr,Stack_RecMom_CC0pi_restr);
      ProduceStack(hRecAngle_CC0pi_restr,Stack_RecAngle_CC0pi_restr);
      ProduceStack(hRecMom_CC0pi,Stack_RecMom_CC0pi);
      ProduceStack(hRecAngle_CC0pi,Stack_RecAngle_CC0pi);
      ProduceStack(hRecMom_CC0pi_full,Stack_RecMom_CC0pi_full);
      ProduceStack(hRecAngle_CC0pi_full,Stack_RecAngle_CC0pi_full);
      ProduceStack(hNTracks_CC0pi,Stack_NTracks_CC0pi);

      ProduceStack(hRecMom_CC1pi_restr,Stack_RecMom_CC1pi_restr);
      ProduceStack(hRecAngle_CC1pi_restr,Stack_RecAngle_CC1pi_restr);
      ProduceStack(hRecMom_CC1pi,Stack_RecMom_CC1pi);
      ProduceStack(hRecAngle_CC1pi,Stack_RecAngle_CC1pi);
      ProduceStack(hRecMom_CC1pi_full,Stack_RecMom_CC1pi_full);
      ProduceStack(hRecAngle_CC1pi_full,Stack_RecAngle_CC1pi_full);
      ProduceStack(hRecMom_CCNpi,Stack_RecMom_CCNpi);
      ProduceStack(hRecAngle_CCNpi,Stack_RecAngle_CCNpi);
      
      ProduceStack(hSampleSecondTrack,Stack_SampleSecondTrack);
      ProduceStack(hTotalChargeSecondTrack,Stack_TotalChargeSecondTrack);
      ProduceStack(hTotalChargePerDistanceSecondTrack,Stack_TotalChargePerDistanceSecondTrack);
      ProduceStack(hTotalChargePerDistanceFirstTrack,Stack_TotalChargePerDistanceFirstTrack);
      ProduceStack(hOpeningAngle,Stack_OpeningAngle);
      ProduceStack(hCoplanarityAngle,Stack_CoplanarityAngle);
      ProduceStack(hMVAMuondiscriminant_1track,Stack_MVAMuondiscriminant_1track);
      ProduceStack(hMVAProtondiscriminant_2tracks,Stack_MVAProtondiscriminant_2tracks);

      ProduceStack(CutNEvents,Stack_CutNEvents);
      
#ifdef DEBUG2
  cout<<"Stack plots produced"<<endl;
#endif
     
      TLegend * leg = new TLegend(0.7,0.8,1,1);
      if(hRecMom[0]) leg->AddEntry(hRecMom[0],"CC0#pi");
      if(hRecMom[1]) leg->AddEntry(hRecMom[1],"MEC");
      if(hRecMom[3]) leg->AddEntry(hRecMom[3],"CC1#pi");
      if(hRecMom[4]) leg->AddEntry(hRecMom[4],"CC1#pi^{0}");
      if(hRecMom[5]) leg->AddEntry(hRecMom[5],"CCother");
      if(hRecMom[6]) leg->AddEntry(hRecMom[6],"NC");
            
      //TFile * fMC = new TFile("plots/MCPlots.root","RECREATE");

      Stack_FCFVTrueEvents->Write();

      Stack_MuCL_Particles->Write();
      Stack_MuCL->Write();
      Stack_MuCL_Lowest->Write();
      Stack_NTracks->Write();
      Stack_RecMom->Write();
      Stack_RecAngle->Write();
      Stack_RecMom_CC0pi_restr->Write();
      Stack_RecAngle_CC0pi_restr->Write();
      Stack_RecMom_CC0pi->Write();
      Stack_RecAngle_CC0pi->Write();
      Stack_RecMom_CC0pi_full->Write();
      Stack_RecAngle_CC0pi_full->Write();
      Stack_NTracks_CC0pi->Write();
      hRecCC0pi_MomentumResolution->Write();
      hRecCC0pi_AngleResolution->Write();
      
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
      Stack_TotalChargeSecondTrack->Write();
      Stack_TotalChargePerDistanceSecondTrack->Write();
      Stack_TotalChargePerDistanceFirstTrack->Write();
      Stack_OpeningAngle->Write();
      Stack_CoplanarityAngle->Write();
      Stack_MVAdiscriminant->Write();
      Stack_MVAMuondiscriminant_1track->Write();
      Stack_MVAProtondiscriminant_2tracks->Write();
      Stack_CutNEvents->Write();
      MCEfficiency->Write();
      MCEfficiency_PTheta->Write();
      MCEfficiency_PTheta_Momentum->Write();
      MCEfficiency_PTheta_Angle->Write();
      MCEfficiency_Energy->Write();

      MCEfficiency_PTheta_ThinBinning->Write();
      MCEfficiency_Momentum_ThinBinning->Write(); 
      MCEfficiency_Angle_ThinBinning->Write();

      MCReconstructionEfficiency_PTheta_ThinBinning->Write();
      MCReconstructionEfficiency_Momentum_ThinBinning->Write(); 
      MCReconstructionEfficiency_Angle_ThinBinning->Write();

      MCReconstructionAndTopologyEfficiency_PTheta_ThinBinning->Write();
      MCReconstructionAndTopologyEfficiency_Momentum_ThinBinning->Write(); 
      MCReconstructionAndTopologyEfficiency_Angle_ThinBinning->Write();

      //MCEfficiency_Momentum_ThinBinning_Relative->Write();
      //MCEfficiency_Angle_ThinBinning_Relative->Write();
    
      TotalCC0piEvent_Energy->Write();
      TotalCC1piEvent_Energy->Write();
      TotalCC0piEvent_PTheta_ThinBinning->Write();
      TotalCC1piEvent_PTheta_ThinBinning->Write();
      leg->Write();
      MuonID->Divide(MuonIDTotal);
      MuonID->Write();
      MuonIDMVA->Divide(MuonIDTotal);
      MuonIDMVA->Write();
      
      TDirectory * PIDCutTuning = gDirectory->mkdir("PIDCutTuning");
      PIDCutTuning->cd();
      for(int fsi=0;fsi<NFSIs;fsi++){
	hMuCL[fsi]->Write();
	hMuCL_Lowest[fsi]->Write();
	hMuCL_2tracks[fsi]->Write();
	hNTracks[fsi]->Write();
	hNTracks_CC0pi[fsi]->Write();
	hMVAMuondiscriminant_1track[fsi]->Write();
	hMVAMuondiscriminantVSPDG_1track[fsi]->Write();
	hMVAMuondiscriminantVSMuonMomentum_1track[fsi]->Write();
	hMVAMuondiscriminantVSDistance_1track[fsi]->Write();
	hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[fsi]->Write();
  	hMVAProtondiscriminant_2tracks[fsi]->Write();
	hMVAMuonVSProtondiscriminant_2tracks[fsi]->Write();
	TrueParticleType_2tracks[fsi]->Write();
	TrueParticleType_1track[fsi]->Write();

	hLMVA_1track[fsi]->Write();
	hLMVA_2tracks[fsi]->Write();
	hLMVA_LargestLMVA[fsi]->Write();
	hLMVA_SmallestLMVA[fsi]->Write();
	
	//hTotalChargePerDistanceFirstTrack[fsi]->Write();
	//hTotalChargePerDistanceSecondTrack[fsi]->Write();
 	//hSampleSecondTrack[fsi]->Write();
      }
      for(int i=0;i<NSimplifiedPDG;i++){
	hMVAMuondiscriminantVSDistance_TrueParticle[i]->Write();
	for(int j=0;j<NSamples;j++){
	  hLMVA_TrueParticle[i][j]->Write();
	}
      }
      hMVA_MuonVsProtonDiscriminant_TrueMuon->Write();
      hMVA_MuonVsProtonDiscriminant_TruePion->Write();
      hMVA_MuonVsProtonDiscriminant_TrueProton->Write();
      
      hMuCL_TrueMuon->Write();
      hMuCL_TruePion->Write();
      hMuCL_TrueProton->Write();

      hMuCL_TrueOthers->Write();
      MuonRec_TruePDG->Write();
      MuonRec_TruePDG_switch->Write();
      gDirectory->cd();
      
      TDirectory * MVAInputVariables = gDirectory->mkdir("MVAInputVariables");
      MVAInputVariables->cd();

      for(int is=0;is<NSamples;is++){

	//Energy deposition
	//Muon
	hEnergyDepositionLength_Muon[is] = (TH2D*) EnergyDepositionLength_Muon[is]->Project3D("zy");
	pEnergyDepositionLength_Muon[is] = (TProfile2D*) EnergyDepositionLength_Muon[is]->Project3DProfile("xy");
	pEnergyDepositionLength_Muon_1D[is] = (TProfile*) hEnergyDepositionLength_Muon[is]->ProfileX(Form("pEnergyDepositionLength_Muon_1D_%d",is));

	EnergyDepositionLength_Muon[is]->Write();
	pEnergyDepositionLength_Muon[is]->Write();
	hEnergyDepositionLength_Muon[is]->Write();
	
	hEnergyDepositionLength_Muon_1D[is] = (TH1D*) EnergyDepositionLength_Muon[is]->Project3D("y");
	for(int ibinx=1;ibinx<=pEnergyDepositionLength_Muon_1D[is]->GetNbinsX();ibinx++){
	  double RMS=((TH1D*) hEnergyDepositionLength_Muon[is]->ProjectionY("htemp",ibinx,ibinx))->GetRMS();
	  double Value=pEnergyDepositionLength_Muon_1D[is]->GetBinContent(ibinx);
	  hEnergyDepositionLength_Muon_1D[is]->SetBinContent(ibinx,Value);
	  hEnergyDepositionLength_Muon_1D[is]->SetBinError(ibinx,RMS/2.);
	}
	hEnergyDepositionLength_Muon_1D[is]->Write();

	//Pion
	hEnergyDepositionLength_Pion[is] = (TH2D*) EnergyDepositionLength_Pion[is]->Project3D("zy");
	pEnergyDepositionLength_Pion[is] = (TProfile2D*) EnergyDepositionLength_Pion[is]->Project3DProfile("xy");
	pEnergyDepositionLength_Pion_1D[is] = (TProfile*) hEnergyDepositionLength_Pion[is]->ProfileX(Form("pEnergyDepositionLength_Pion_1D_%d",is));

	EnergyDepositionLength_Pion[is]->Write();
	pEnergyDepositionLength_Pion[is]->Write();
	hEnergyDepositionLength_Pion[is]->Write();
	
	hEnergyDepositionLength_Pion_1D[is] = (TH1D*) EnergyDepositionLength_Pion[is]->Project3D("y");
	for(int ibinx=1;ibinx<=pEnergyDepositionLength_Pion_1D[is]->GetNbinsX();ibinx++){
	  double RMS=((TH1D*) hEnergyDepositionLength_Pion[is]->ProjectionY("htemp",ibinx,ibinx))->GetRMS();
	  double Value=pEnergyDepositionLength_Pion_1D[is]->GetBinContent(ibinx);
	  hEnergyDepositionLength_Pion_1D[is]->SetBinContent(ibinx,Value);
	  hEnergyDepositionLength_Pion_1D[is]->SetBinError(ibinx,RMS/2.);
	}
	hEnergyDepositionLength_Pion_1D[is]->Write();

	//Proton
	hEnergyDepositionLength_Proton[is] = (TH2D*) EnergyDepositionLength_Proton[is]->Project3D("zy");
	pEnergyDepositionLength_Proton[is] = (TProfile2D*) EnergyDepositionLength_Proton[is]->Project3DProfile("xy");
	pEnergyDepositionLength_Proton_1D[is] = (TProfile*) hEnergyDepositionLength_Proton[is]->ProfileX(Form("pEnergyDepositionLength_Proton_1D_%d",is));

	EnergyDepositionLength_Proton[is]->Write();
	pEnergyDepositionLength_Proton[is]->Write();
	hEnergyDepositionLength_Proton[is]->Write();
	
	hEnergyDepositionLength_Proton_1D[is] = (TH1D*) EnergyDepositionLength_Proton[is]->Project3D("y");
	for(int ibinx=1;ibinx<=pEnergyDepositionLength_Proton_1D[is]->GetNbinsX();ibinx++){
	  double RMS=((TH1D*) hEnergyDepositionLength_Proton[is]->ProjectionY("htemp",ibinx,ibinx))->GetRMS();
	  double Value=pEnergyDepositionLength_Proton_1D[is]->GetBinContent(ibinx);
	  hEnergyDepositionLength_Proton_1D[is]->SetBinContent(ibinx,Value);
	  hEnergyDepositionLength_Proton_1D[is]->SetBinError(ibinx,RMS/2.);
	}
	hEnergyDepositionLength_Proton_1D[is]->Write();


      }
      gDirectory->cd();

      
#ifdef DEBUG2
      cout<<"Plots wrote"<<endl;
#endif
      delete leg;
      
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
#ifdef DEBUG2
  cout<<"Writing output is closed. Now delete the pointers~."<<endl;
#endif

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
      delete hSampleSecondTrack[fsi];
      delete hTotalChargeSecondTrack[fsi];
      delete hTotalChargePerDistanceSecondTrack[fsi];
      delete hTotalChargePerDistanceFirstTrack[fsi];
      delete hOpeningAngle[fsi];
      delete hCoplanarityAngle[fsi];
      delete hMVAMuondiscriminant_1track[fsi];
      delete hMVAProtondiscriminant_2tracks[fsi];
      delete hMVAMuonVSProtondiscriminant_2tracks[fsi];
    }
    
    
#ifdef DEBUG2
  cout<<"End of deleting table of pointers."<<endl;
#endif

  delete hMVA_MuonVsProtonDiscriminant_TrueMuon;
  delete hMVA_MuonVsProtonDiscriminant_TruePion;
  delete hMVA_MuonVsProtonDiscriminant_TrueProton;
  
  delete hMuCL_TrueMuon;
  delete hMuCL_TruePion;
  delete hMuCL_TrueProton;
  
  delete MCEfficiency;
  delete MCEfficiency_PTheta;
  delete MCEfficiency_PTheta_Momentum;
  delete MCEfficiency_PTheta_Angle;
  delete MCEfficiency_Energy;
  delete TotalCC0piEvent_Energy;
  delete TotalCC1piEvent_Energy;
  delete MCEfficiency_PTheta_ThinBinning;
  delete TotalCC0piEvent_PTheta_ThinBinning;
  delete TotalCC1piEvent_PTheta_ThinBinning;
  delete hRecCC0pi_MomentumResolution;
  delete hRecCC0pi_AngleResolution;

#ifdef DEBUG2
  cout<<"End of deleting pointers."<<endl;
#endif

  }

  for(int i=0;i<NBinsRecMom;i++){
    delete DataSelected[i];
    delete BkgSelected[i];
    delete DataSelected_full[i];
    delete BkgSelected_full[i];
  }
  delete DataSelected;
  delete BkgSelected;
  delete DataSelected_full;
  delete BkgSelected_full;

  for(int i=0;i<NBinsTrueMom;i++){
    delete Efficiency[i];
    delete Efficiency_full[i];
    delete TotalCC0piEvent[i];
    delete TotalCC1piEvent[i];
  }
  delete Efficiency;
  delete Efficiency_full;
  delete TotalCC0piEvent;
  delete TotalCC1piEvent;

  
  for(int i=0;i<NBinsTrueMom;i++){
    for(int j=0;j<NBinsTrueAngle;j++){
      for(int k=0;k<NBinsRecMom;k++){
	delete MCSelected[i][j][k];
	delete MCSelected_full[i][j][k];
      }
      delete MCSelected[i][j];
      delete MCSelected_full[i][j];
    }
    delete MCSelected[i];
    delete MCSelected_full[i];
  }
  delete MCSelected;
  delete MCSelected_full;

  delete rand;
  
#ifdef MOMENTUMSHIFT
delete Shift;
#endif
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
  while ((c = getopt(argc, argv, "w:ms:e:v:i:o:p:WB")) != -1) {
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
    case 'B':
      BkgSub=true;     //cout<<"Please enter the amount of data you'd like to mimic (I will adjust MC stat.) in unit of 1e21 POT:"<<endl;
      break;
    }
  }
    
  cout<<"welcome"<<endl;
  cout<<"Detector is "<<(isPM ? "PM" : "WM" )<<endl;
  //  XS->Xsec::Initialize();
  InitializeGlobal(isPM,Sample);

  if(!isPM){
    ProtonCut=ProtonCut_WM;
    MuonCut=MuonCut_WM;
  }
  
  if(Sample==2){
    MuonMVACut = MuonMVACut_CC1pi;
    ProtonMVACut = ProtonMVACut_CC1pi;
  }
  else{
    MuonMVACut = MuonMVACut_CC0pi;
    ProtonMVACut = ProtonMVACut_CC0pi;
  }
  cout<<"Defined the MVA muon cut as:"<<MuonMVACut<<", and proton cut as:"<<ProtonMVACut<<endl;
   TChain * chain = new TChain("wtree");
   TChain * chainMVA = new TChain("wtreeMVA");
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
     cout<<NMCfiles<<endl;
     for(int i=1;i<NMCfiles;i++){
     	       //sprintf(fName,i);
	        //sprintf(fName,"%s",InNameEvent);
	        //sprintf(fName,"root_input/CC0piTree%d.root",i);
	        //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_ReWeight.root",i);
	        //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root",i);
	        //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Old_%d_Plan.root",i);
	        //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan.root",i);
	        //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d_Plan_RandomPE.root",i);
       fName=Form(InNameEvent,i);
       cout<<fName<<endl;
       chain->Add(fName);
#ifdef MVATRAINING
       cout<<"MVA training, add the tree:"<<fName<<endl;
       chainMVA->Add(fName);
#endif
     }
   }
   
   if(SelectedError_Source>=Systematics_Detector_Start && SelectedError_Source<=Systematics_Detector_End){
     Systematics_Detector=true;
     
     CC0piDistributions(chain,chainMVA,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM);	 
   }
   else if(SelectedError_Source>=Systematics_Flux_Start && SelectedError_Source<=Systematics_Flux_End){
     Systematics_Flux=true;
     TFile * FluxError = new TFile("Flux/ErrorVarFinal.root");
     FluxVector = new TVectorD();
     cout<<"variation in flux #"<<SelectedError_Variation<<endl;
     sprintf(Name,"Var[%d]",( (int) SelectedError_Variation ));
     FluxVector = (TVectorD*) FluxError->Get(Name);
     File_Number= ( (int) SelectedError_Variation);

     //Distributions(chain,Data,Sample,IsPlots,SelectedError_Source,SelectedError_Variation,OutNameEvent);
     //Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent);
     CC0piDistributions(chain,chainMVA,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM);
     FluxError->Close();
   }
   else if(SelectedError_Source>=Systematics_Xsec_Start && SelectedError_Source<=Systematics_Xsec_End){
     Systematics_Xsec=true;
     Xsec_dial=((int) SelectedError_Variation);
     File_Number=Xsec_dial;
     cout << "XS error is checked, using the dial #" << Xsec_dial << endl;
     CC0piDistributions(chain,chainMVA,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM); 
   }
   else{
     cout<<"Simple selection on the MC or data"<<endl;
     CC0piDistributions(chain,chainMVA,Data,Sample,IsPlots,Systematics_Flux,File_Number,*FluxVector,Systematics_Xsec,Xsec_dial,Systematics_Detector,ErrorType,OutNameEvent,isPM);
   }
   
   return 0;
}

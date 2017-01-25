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
#include <TApplication.h>
#include <TBox.h>
#include "setup.h"
//#include "Xsec.cc"
//#include "Reconstruction.cc"
//Xsec * XS = new Xsec();
//double IronCarbonRatio=7.87/1.03;
TLegend * leg_Sample4;
TLegend * leg_Sample4_Data;

void ProduceStack(TH1D * h[11], THStack * hStack){
  for(int fsi=1;fsi<11;fsi++){
    /*
      if(fsi==1) h[fsi]->SetFillColor(kRed);
      else if(fsi==3) h[fsi]->SetFillColor(kBlue);
      else if(fsi==4) h[fsi]->SetFillColor(kAzure+10);
      else if(fsi==5) h[fsi]->SetFillColor(kGreen-2);
      else if(fsi==6) h[fsi]->SetFillColor(kGray);
      else if(fsi==7) h[fsi]->SetFillColor(kYellow);
      else if(fsi==8) h[fsi]->SetFillColor(kYellow+2);
      else if(fsi==9) h[fsi]->SetFillColor(kYellow+4);
      else if(fsi==10) h[fsi]->SetFillColor(kMagenta);
    */
    h[fsi]->GetYaxis()->SetTitleOffset(1.3);
    if(fsi==1) h[fsi]->SetFillColor(kRed);
    //else if(fsi==2) h[fsi]->SetFillColor(kRed/*kOrange-3*/);
    else if(fsi==3) h[fsi]->SetFillColor(kBlue+2);
    else if(fsi==4) h[fsi]->SetFillColor(kAzure+7);
    else if(fsi==5) h[fsi]->SetFillColor(kAzure+10);
    else if(fsi==6) h[fsi]->SetFillColor(kGreen-2);
    else if(fsi==7) h[fsi]->SetFillColor(kYellow);
    else if(fsi==8) h[fsi]->SetFillColor(kYellow+2);
    else if(fsi==9) h[fsi]->SetFillColor(kYellow+4);
    else if(fsi==10) h[fsi]->SetFillColor(kMagenta);

    h[fsi]->SetLineColor(1);
    h[fsi]->SetLineWidth(2);
    
    /*
      if(fsi==1) h[fsi]->SetFillStyle(1001);
    
      else if(fsi==3) h[fsi]->SetFillStyle(3105);
      else if(fsi==4) h[fsi]->SetFillStyle(3154);
      else if(fsi==5) h[fsi]->SetFillStyle(3158);
    
      else if(fsi==3) h[fsi]->SetFillStyle(3002);
      else if(fsi==4) h[fsi]->SetFillStyle(3002);
      else if(fsi==5) h[fsi]->SetFillStyle(3002);
      else if(fsi==6) h[fsi]->SetFillStyle(3002);
      else if(fsi==7) h[fsi]->SetFillStyle(3002);
      else if(fsi==8) h[fsi]->SetFillStyle(3002);
      else if(fsi==9) h[fsi]->SetFillStyle(3002);
      else if(fsi==10) h[fsi]->SetFillStyle(3002);
    */
    hStack->Add(h[fsi]);
  }
}

void ProduceStack2(TH1D * h[11], THStack * hStack, bool Last){
  for(int fsi=1;fsi<11;fsi++){
    if(fsi==1) h[fsi]->SetFillColor(kRed);
   
    h[fsi]->GetYaxis()->SetTitleOffset(1.1);
    if(fsi==1) h[fsi]->SetFillColor(kRed);
    //else if(fsi==2) h[fsi]->SetFillColor(kRed/*kOrange-3*/);
    else if(fsi==3) h[fsi]->SetFillColor(kBlue+2);
    else if(fsi==4) h[fsi]->SetFillColor(kAzure+7);
    else if(fsi==5) h[fsi]->SetFillColor(kAzure+10);
    else if(fsi==6) h[fsi]->SetFillColor(kGreen-2);
    else if(fsi==7) h[fsi]->SetFillColor(kYellow);
    else if(fsi==8) h[fsi]->SetFillColor(kYellow+2);
    else if(fsi==9) h[fsi]->SetFillColor(kYellow+4);
    else if(fsi==10) h[fsi]->SetFillColor(kMagenta);
      
    //h[fsi]->SetLineColor(1);
    //h[fsi]->SetLineWidth(2);
    /*
      else if(fsi==3) h[fsi]->SetFillColor(kBlue);
      else if(fsi==4) h[fsi]->SetFillColor(kAzure+10);
      else if(fsi==5) h[fsi]->SetFillColor(kGreen-2);
      else if(fsi==6) h[fsi]->SetFillColor(kGray);
      else if(fsi==7) h[fsi]->SetFillColor(kYellow);
      else if(fsi==8) h[fsi]->SetFillColor(kYellow+2);
      else if(fsi==9) h[fsi]->SetFillColor(kYellow+4);
      else if(fsi==10) h[fsi]->SetFillColor(kMagenta);
    */
    if(fsi==1) h[fsi]->SetLineColor(kRed);
    //else if(fsi==2) h[fsi]->SetLineColor(kRed/*kOrange-3*/);
    else if(fsi==3) h[fsi]->SetLineColor(kBlue+2);
    else if(fsi==4) h[fsi]->SetLineColor(kAzure+7);
    else if(fsi==5) h[fsi]->SetLineColor(kAzure+10);
    else if(fsi==6) h[fsi]->SetLineColor(kGreen-2);
    else if(fsi==7) h[fsi]->SetLineColor(kYellow);
    else if(fsi==8) h[fsi]->SetLineColor(kYellow+2);
    else if(fsi==9) h[fsi]->SetLineColor(kYellow+4);
    else if(fsi==10) h[fsi]->SetLineColor(kMagenta);

    if(Last){
      if(fsi==1) h[fsi]->SetLineColor(1);
      else if(fsi==3) h[fsi]->SetLineColor(1);
      else if(fsi==4) h[fsi]->SetLineColor(1);
      else if(fsi==5) h[fsi]->SetLineColor(1);
      else if(fsi==6) h[fsi]->SetLineColor(1);
      else if(fsi==7) h[fsi]->SetLineColor(1);
      else if(fsi==8) h[fsi]->SetLineColor(1);
      else if(fsi==9) h[fsi]->SetLineColor(1);
      else if(fsi==10) h[fsi]->SetLineColor(1);
    }
    hStack->Add(h[fsi]);
  }
}

int main(int argc, char **argv)
{
  TApplication * app = new TApplication("app", &argc, argv);

  gROOT->SetStyle("Plain");
  //gStyle->SetFillColor(kWhite);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  //  XS->Xsec::Initialize();
  InitializeGlobal();

  char Name[256];
  //TFile * _file0 = new TFile("CC0pi_1000.root");
  //TTree * wtree = (TTree*) _file0->Get("wtree");  
  TChain * wtree = new TChain("wtree");
  for(int i=1;i<=1000;i++){
    if((i==165 || i==166 || i==167 || i==168) || (i==665 || i==666 || i==667 || i==668 || i==555)) continue;
    //if(i==17 || i==74) continue;
    sprintf(Name,"AllCC0pi5/CC0piTree%d.root",i);
    //sprintf(Name,"TempBNJ/CC0piTree%d.root",i);
    //sprintf(Name,"Temp0pi/CC0piTree%d.root",i);
    //sprintf(Name,"PEError/CC0piTree_%d_PE1.root",i);
    cout<<Name<<endl;
    wtree->Add(Name);
  }
  //if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
  //else continue;

  int nevt=(int) wtree->GetEntries();
 
  double POT;
  int FSIInt;
  int Num_Int;
  int nTracks[10];
  double weight;
  bool IsFV;
  bool IsDetected[10];
  double Enu;
  int nIngBasRec;
  double Errorweight;
  double TrueMomentumMuon, TrueAngleMuon;
  bool IsSand;
  bool IsAnti;
  bool IsBkgH;
  bool IsBkgV;
  bool IsNuE;
  int VertexZ[10];
  double TrackAngle[10][20];
  int TypeOfTrack[10][20];
  double CLMuon[10][20];// = new vector<double> [10][20];
  double TrackWidth[10][20];
  bool GeometricTrack[10][20];
  int LastChannelINGRIDX[10][20];
  int LastChannelINGRIDY[10][20];

  double Momentum[10][20];
  double IronDistance[10][20];
  double PlasticDistance[10][20];
  int Sample[10][20];
  bool IsReconstructed[10][20];
  double CriteriaAngleX[10][20];
  double CriteriaAngleY[10][20];
  double CriteriaHalfWayX[10][20];
  double CriteriaHalfWayY[10][20];
  double ReWeight[175];
  int GoodSpill;
  int Spill;

  TBranch* BrMC_FSIInt = wtree->GetBranch("FSIInt");
  BrMC_FSIInt->SetAddress(&FSIInt);
  wtree->SetBranchAddress("FSIInt",&FSIInt);

  TBranch* BrMC_nIngBasRec = wtree->GetBranch("nIngBasRec");
  BrMC_nIngBasRec->SetAddress(&nIngBasRec);
  wtree->SetBranchAddress("nIngBasRec",&nIngBasRec);

  TBranch * BrMC_Num_Int = wtree->GetBranch("InteractionType");
  BrMC_Num_Int->SetAddress(&Num_Int);
  wtree->SetBranchAddress("InteractionType",&Num_Int);
  
  TBranch * BrMC_nTracks = wtree->GetBranch("nTracks[10]");
  BrMC_nTracks->SetAddress(&nTracks);
  wtree->SetBranchAddress("nTracks[10]",&nTracks);
 
  TBranch * BrMC_weight = wtree->GetBranch("weight");
  BrMC_weight->SetAddress(&weight);
  wtree->SetBranchAddress("weight",&weight);

  TBranch * BrMC_IsSand = wtree->GetBranch("IsSand");
  BrMC_IsSand->SetAddress(&IsSand);
  wtree->SetBranchAddress("IsSand",&IsSand);
  
  TBranch * BrMC_IsAnti = wtree->GetBranch("IsAnti");
  BrMC_IsAnti->SetAddress(&IsAnti);
  wtree->SetBranchAddress("IsAnti",&IsAnti);
  
  TBranch * BrMC_IsBkgH = wtree->GetBranch("IsBkgH");
  BrMC_IsBkgH->SetAddress(&IsBkgH);
  wtree->SetBranchAddress("IsBkgH",&IsBkgH);

  TBranch * BrMC_IsBkgV = wtree->GetBranch("IsBkgV");
  BrMC_IsBkgV->SetAddress(&IsBkgV);
  wtree->SetBranchAddress("IsBkgV",&IsBkgV);

  TBranch * BrMC_IsNuE = wtree->GetBranch("IsNuE");
  BrMC_IsNuE->SetAddress(&IsNuE);
  wtree->SetBranchAddress("IsNuE",&IsNuE);

  TBranch * BrMC_Enu = wtree->GetBranch("Enu");
  BrMC_Enu->SetAddress(&Enu);
  wtree->SetBranchAddress("Enu",&Enu);

  TBranch * BrMC_TrueMomentumMuon = wtree->GetBranch("TrueMomentumMuon");
  BrMC_TrueMomentumMuon->SetAddress(&TrueMomentumMuon);
  wtree->SetBranchAddress("TrueMomentumMuon",&TrueMomentumMuon);

  TBranch * BrMC_TrueAngleMuon = wtree->GetBranch("TrueAngleMuon");
  BrMC_TrueAngleMuon->SetAddress(&TrueAngleMuon);
  wtree->SetBranchAddress("TrueAngleMuon",&TrueAngleMuon);

  TBranch * BrMC_IsFV = wtree->GetBranch("IsFV");
  BrMC_IsFV->SetAddress(&IsFV);
  wtree->SetBranchAddress("IsFV",&IsFV);
  
  TBranch * BrMC_IsDetected = wtree->GetBranch("IsDetected[10]");
  BrMC_IsDetected->SetAddress(&IsDetected);
  wtree->SetBranchAddress("IsDetected[10]",&IsDetected);
  /*
    TBranch * BrMC_VertexZ = wtree->GetBranch("VertexZ[10]");
    BrMC_VertexZ->SetAddress(&VertexZ);
    wtree->SetBranchAddress("VertexZ[10]",&VertexZ);
  */
  TBranch * BrMC_TrackAngle = wtree->GetBranch("TrackAngle[10][20]");
  BrMC_TrackAngle->SetAddress(TrackAngle);
  wtree->SetBranchAddress("TrackAngle[10][20]",TrackAngle);
 
  TBranch * BrMC_TypeOfTrack = wtree->GetBranch("TypeOfTrack[10][20]");
  BrMC_TypeOfTrack->SetAddress(TypeOfTrack);
  wtree->SetBranchAddress("TypeOfTrack[10][20]",TypeOfTrack);

  TBranch * BrMC_CLMuon = wtree->GetBranch("CLMuon[10][20]");
  BrMC_CLMuon->SetAddress(&CLMuon);
  wtree->SetBranchAddress("CLMuon[10][20]",&CLMuon);

  TBranch * BrMC_Momentum = wtree->GetBranch("Momentum[10][20]");
  BrMC_Momentum->SetAddress(Momentum);
  wtree->SetBranchAddress("Momentum[10][20]",Momentum);

  TBranch * BrMC_IronDistance = wtree->GetBranch("IronDistance[10][20]");
  BrMC_IronDistance->SetAddress(IronDistance);
  wtree->SetBranchAddress("IronDistance[10][20]",IronDistance);

  TBranch * BrMC_PlasticDistance = wtree->GetBranch("PlasticDistance[10][20]");
  BrMC_PlasticDistance->SetAddress(PlasticDistance);
  wtree->SetBranchAddress("PlasticDistance[10][20]",PlasticDistance);

  TBranch * BrMC_Sample = wtree->GetBranch("Sample[10][20]");
  BrMC_Sample->SetAddress(Sample);
  wtree->SetBranchAddress("Sample[10][20]",Sample);

  TBranch * BrMC_IsReconstructed = wtree->GetBranch("IsReconstructed[10][20]");
  BrMC_IsReconstructed->SetAddress(IsReconstructed);
  wtree->SetBranchAddress("IsReconstructed[10][20]",IsReconstructed);

  TBranch * BrMC_CriteriaAngleX = wtree->GetBranch("CriteriaAngleX[10][20]");
  BrMC_CriteriaAngleX->SetAddress(CriteriaAngleX);
  wtree->SetBranchAddress("CriteriaAngleX[10][20]",CriteriaAngleX);
  TBranch * BrMC_CriteriaAngleY = wtree->GetBranch("CriteriaAngleY[10][20]");
  BrMC_CriteriaAngleY->SetAddress(CriteriaAngleY);
  wtree->SetBranchAddress("CriteriaAngleY[10][20]",CriteriaAngleY);
  TBranch * BrMC_CriteriaHalfWayX = wtree->GetBranch("CriteriaHalfWayX[10][20]");
  BrMC_CriteriaHalfWayX->SetAddress(CriteriaHalfWayX);
  wtree->SetBranchAddress("CriteriaHalfWayX[10][20]",CriteriaHalfWayX);
  TBranch * BrMC_CriteriaHalfWayY = wtree->GetBranch("CriteriaHalfWayY[10][20]");
  BrMC_CriteriaHalfWayY->SetAddress(CriteriaHalfWayY);
  wtree->SetBranchAddress("CriteriaHalfWayY[10][20]",CriteriaHalfWayY);

  TBranch * BrMC_ReWeight = wtree->GetBranch("ReWeight[175]");
  BrMC_ReWeight->SetAddress(ReWeight);
  wtree->SetBranchAddress("ReWeight[175]",ReWeight);

  TBranch * BrMC_LastChannelINGRIDX = wtree->GetBranch("LastChannelINGRIDX[10][20]");
  BrMC_LastChannelINGRIDX->SetAddress(LastChannelINGRIDX);
  wtree->SetBranchAddress("LastChannelINGRIDX[10][20]",LastChannelINGRIDX);

  TBranch * BrMC_LastChannelINGRIDY = wtree->GetBranch("LastChannelINGRIDY[10][20]");
  BrMC_LastChannelINGRIDY->SetAddress(LastChannelINGRIDY);
  wtree->SetBranchAddress("LastChannelINGRIDY[10][20]",LastChannelINGRIDY);

  TBranch * BrMC_TrackWidth = wtree->GetBranch("TrackWidth[10][20]");
  BrMC_TrackWidth->SetAddress(TrackWidth);
  wtree->SetBranchAddress("TrackWidth[10][20]",TrackWidth);

  TBranch * BrMC_GeometricTrack = wtree->GetBranch("GeometricTrack[10][20]");
  BrMC_GeometricTrack->SetAddress(GeometricTrack);
  wtree->SetBranchAddress("GeometricTrack[10][20]",GeometricTrack);


    int NBinsIron=12;

    double DistIronBin[NBinsIron+1];
    for(int i=0;i<NBinsIron+1;i++){
    DistIronBin[0]=0;
    DistIronBin[1]=10;
    if(i>1 && i<=9) DistIronBin[i]=20+5*(i-2);
    DistIronBin[10]=60;
    DistIronBin[11]=70;
    DistIronBin[12]=100;
    }

    int NBinsMom=5;
    double DistMomBin[NBinsMom+1];
    for(int i=0;i<NBinsMom+1;i++){
    DistMomBin[0]=0;
    DistMomBin[1]=0.4;
    DistMomBin[2]=0.6;
    DistMomBin[3]=0.8;
    DistMomBin[4]=1.0;
    //DistMomBin[5]=2.0;
    DistMomBin[5]=10.0;
    }
    const int NBinsAngle=10;
    double DistAngleBin[NBinsAngle+1];
    for(int i=0;i<NBinsAngle+1;i++){
    if(i!=NBinsAngle) DistAngleBin[i]=i*5;
    else DistAngleBin[i]=90;//TMath::Cos(90*TMath::Pi()/180);
    cout<<"bin="<<i<<", bin value="<<DistAngleBin[i]<<endl;
    }

  TRandom3 * r1 = new TRandom3();
  const int NBinsHalf=50;
  double DistHalfBin[NBinsHalf+1];
  for(int i=0;i<NBinsHalf+1;i++){
    if(i!=NBinsHalf) DistHalfBin[i]=i;
    else DistHalfBin[i]=50;//TMath::Cos(90*TMath::Pi()/180);
  }

  TH1D * NTracksTotal = new TH1D("NTracksTotal","",6,0,6);
  TH1D * NTracksTotalSample4 = new TH1D("NTracksTotalSample4","",6,0,6);
  TH1D * CLMuonTotal = new TH1D("CLMuonTotal","",100,0,1);
  TH1D * SampleTotal = new TH1D("SampleTotal","",6,0,6);
  TH1D * SampleTotalSample4 = new TH1D("SampleTotalSample4","",6,0,6);
  TH1D * IronDistanceTotal = new TH1D("IronDistanceTotal","",20,0,100);
  TH1D * IronDistanceTotalSample4_1Track = new TH1D("IronDistanceTotalSample4_1Track","",20,0,100);
  TH1D * IronDistanceTotalSample4_2Tracks = new TH1D("IronDistanceTotalSample4_2Tracks","",20,0,100);

  TH1D * AngleTotal = new TH1D("AngleTotal","",30,0,90);

  TH1D * NTracksTotalSample4_Data = new TH1D("NTracksTotalSample4_Data","",6,0,6);
  TH1D * NTracksTotal_Data = new TH1D("NTracksTotal_Data","",6,0,6);
  TH1D * CLMuonTotal_Data = new TH1D("CLMuonTotal_Data","",100,0,1);
  TH1D * SampleTotal_Data = new TH1D("SampleTotal_Data","",6,0,6);
  TH1D * SampleTotalSample4_Data = new TH1D("SampleTotalSample4_Data","",6,0,6);
  TH1D * IronDistanceTotal_Data = new TH1D("IronDistanceTotal_Data","",20,0,100);
  TH1D * AngleTotal_Data = new TH1D("AngleTotal_Data","",30,0,90);
  TH1D * IronDistanceTotalSample4_1Track_Data = new TH1D("IronDistanceTotalSample4_1Track_Data","",20,0,100);
  TH1D * IronDistanceTotalSample4_2Tracks_Data = new TH1D("IronDistanceTotalSample4_2Tracks_Data","",20,0,100);


  TH1D * Selected_NTracks[6][11];//6=Sample of the reconstruction, 6=FSIInt
  TH1D * Reconstructed_NTracks[6][11];
  TH1D * Interacting_Enu[6][11];
  TH1D * Interacting_NTracks[6][11];
  TH1D * Selected_Momentum[6][11];
  TH1D * Reconstructed_Momentum[6][11];
  TH1D * Interacting_Momentum[6][11];
  TH1D * Selected_IronDistance[6][11];
  TH1D * Selected_IronDistance_MC[6];
  TH1D * Reconstructed_IronDistance[6][11];
  TH1D * Interacting_IronDistance[6][11];

  TH1D * Selected_Angle[6][11];
  TH1D * Reconstructed_Angle[6][11];
  TH1D * Interacting_Angle[6][11];

  TH1D * Selected_ReconstructedAngle[6][11];
  TH1D * Reconstructed_ReconstructedAngle[6][11];
  TH1D * Interacting_ReconstructedAngle[6][11];

  TH1D * Selected_Higher[4];
  TH1D * Interacting_Higher[4];
  TH1D * Selected_Lower[4];
  TH1D * Interacting_Lower[4];
  TH1D * InFV_Interacting[11];
  TH1D * InFV_Reconstructed[11];
  //TH1D * POTCount = new TH1D("POTCount","",2,0,2);

  TH2D * Selected2D = new TH2D("Selected2D","",NBinsIron,DistIronBin,NBinsAngle,DistAngleBin);
  TH2D * Selected2D_Data = new TH2D("Selected2D_Data","",NBinsIron,DistIronBin,NBinsAngle,DistAngleBin);
  TH1D * PID_Selected_IronDistance[6][11];
  TH1D * PIDAndMatching_Selected_IronDistance[6][11];
  TH1D * PID_Selected_ReconstructedAngle[6][11];
  TH1D * PIDAndMatching_Selected_ReconstructedAngle[6][11];
  TH1D * INGRIDTrackWidth[11];// = new TH1D("Problematic_TrackWidth","",100,0.,10.);
  TH1D * INGRIDAngleCriteria[11];// = new TH1D("Problematic_TrackWidth","",100,0.,10.);
  TH1D * INGRIDPositionCriteria[11];// = new TH1D("Problematic_TrackWidth","",100,0.,10.);

  TH1D * Reconstructed_NTracks_ATrackIsSample4[11];
  TH1D * Reconstructed_NTracks_ATrackIsSample5[11];
  TH1D * Selected_NTracks_ATrackIsSample4[11];
  TH1D * Selected_NTracks_ATrackIsSample5[11];

  TH1D * PID_Selected_IronDistance_Data[6];
  TH1D * PIDAndMatching_Selected_IronDistance_Data[6];
  TH1D * PID_Selected_ReconstructedAngle_Data[6];
  TH1D * PIDAndMatching_Selected_ReconstructedAngle_Data[6];
  TH1D * Selected_NTracks_ATrackIsSample4_Data;
  TH1D * Selected_NTracks_ATrackIsSample5_Data;
  TH1D * INGRIDTrackWidth_Data;// = new TH1D("Problematic_TrackWidth","",100,0.,10.);
  TH1D * INGRIDAngleCriteria_Data;// = new TH1D("Problematic_TrackWidth","",100,0.,10.);
  TH1D * INGRIDPositionCriteria_Data;// = new TH1D("Problematic_TrackWidth","",100,0.,10.);



  TH1D * Selected_NTracks_Data[6];//6=Sample of the reconstruction, 6=FSIInt
  TH1D * Reconstructed_NTracks_Data[6];
  TH1D * Interacting_NTracks_Data[6];
  TH1D * Selected_IronDistance_Data[6];
  TH1D * Reconstructed_IronDistance_Data[6];
  TH1D * Interacting_IronDistance_Data[6];
  TH1D * Selected_ReconstructedAngle_Data[6];
  TH1D * Reconstructed_ReconstructedAngle_Data[6];
  TH1D * Interacting_ReconstructedAngle_Data[6];
  TH1D * CLMuon_Distribution[6];
  TH1D * CLMuon_Distribution_Data[6];
  TH1D * CLMuon_Distribution_Particle[6][3];

  TH2D * hCriteriaAngleX[6];
  TH2D * hCriteriaAngleX_Data[6];
  TH2D * hCriteriaAngleY[6];
  TH2D * hCriteriaAngleY_Data[6];

  TH2D * hCriteriaHalfWayX[6];
  TH2D * hCriteriaHalfWayX_Data[6];
  TH2D * hCriteriaHalfWayY[6];
  TH2D * hCriteriaHalfWayY_Data[6];

  TH2D * hCriteriaHalfWayXY_Data[NBinsIron];
  TH2D * hCriteriaAngleXY_Data[NBinsIron];
  TH2D * hCriteriaHalfWayXY[NBinsIron];
  TH2D * hCriteriaAngleXY[NBinsIron];
  TH2D * MomentumXIron_ReconstructedForTrueMuon = new TH2D("MomentumXIron_ReconstructedForTrueMuon","",100,0,10,100,0,100);

  TH1D * Problematic_PlasticDistance = new TH1D("Problematic_PlasticDistance","",100,0.,100.);
  TH1D * Problematic_Angle = new TH1D("Problematic_Angle","",30,0.,90.);
  TH1D * Problematic_NTracks = new TH1D("Problematic_NTracks","",10,0.,10.);
  TH1D * Problematic_CLMuon = new TH1D("Problematic_CLMuon","",100,0.,1.);
  TH1D * Problematic_TrackWidth = new TH1D("Problematic_TrackWidth","",100,0.,10.);
  TH1D * Problematic_Geom = new TH1D("Problematic_Geom","",2,0.,2.);
  TH1D * Problematic_LastChannel = new TH1D("Problematic_LastChannel","",25,0.,25.);

  TH1D * Problematic_PlasticDistance_Data = new TH1D("Problematic_PlasticDistance_Data","",100,0.,100.);
  TH1D * Problematic_Angle_Data = new TH1D("Problematic_Angle_Data","",30,0.,90.);
  TH1D * Problematic_NTracks_Data = new TH1D("Problematic_NTracks_Data","",10,0.,10.);
  TH1D * Problematic_CLMuon_Data = new TH1D("Problematic_CLMuon_Data","",100,0.,1.);
  TH1D * Problematic_TrackWidth_Data = new TH1D("Problematic_TrackWidth_Data","",100,0.,10.);
  TH1D * Problematic_Geom_Data = new TH1D("Problematic_Geom_Data","",2,0.,2.);
  TH1D * Problematic_LastChannel_Data = new TH1D("Problematic_LastChannel_Data","",25,0.,25.);

  TH2D * TrackWidthXIronDistance = new TH2D("TrackWidthXIronDistance","",20,0,100,100,0.,10.);
  TH2D * TrackWidthXIronDistance_Data = new TH2D("TrackWidthXIronDistance_Data","",20,0,100,100,0.,10.);
  TH2D * TrackWidthXIronDistance_NuE = new TH2D("TrackWidthXIronDistance_NuE","",20,0,100,100,0.,10.);

  TH2D * AngleRecXTrue = new TH2D("AngleRecXTrue","",90,0.,90.,90,0.,90.);
  TH2D * AngleRecXTrue_Selected = new TH2D("AngleRecXTrue_Selected","",90,0.,90.,90,0.,90.);
  AngleRecXTrue->Sumw2();
  AngleRecXTrue_Selected->Sumw2();

  //D("
  for(int i=0;i<NBinsIron;i++){
    sprintf(Name,"hCriteriaAngleXY[%d]",i);
    hCriteriaAngleXY[i] = new TH2D(Name,"",30,0,90,30,0,90);//50,0.,50.);
    hCriteriaAngleXY[i]->Sumw2();
    sprintf(Name,"hCriteriaHalfWayXY[%d]",i);
    hCriteriaHalfWayXY[i] = new TH2D(Name,"",NBinsHalf,DistHalfBin,NBinsHalf,DistHalfBin);//50,0.,50.);
    hCriteriaHalfWayXY[i]->Sumw2();

    sprintf(Name,"hCriteriaAngleXY_Data[%d]",i);
    hCriteriaAngleXY_Data[i] = new TH2D(Name,"",30,0,90,30,0,90);//50,0.,50.);
    hCriteriaAngleXY_Data[i]->Sumw2();
    sprintf(Name,"hCriteriaHalfWayXY_Data[%d]",i);
    hCriteriaHalfWayXY_Data[i] = new TH2D(Name,"",NBinsHalf,DistHalfBin,NBinsHalf,DistHalfBin);//50,0.,50.);
    hCriteriaHalfWayXY_Data[i]->Sumw2();
  }

  for(int i=0;i<6;i++){
    sprintf(Name,"hCriteriaHalfWayX[%d]",i);
    hCriteriaHalfWayX[i] = new TH2D(Name,"",20,0,100,NBinsHalf,DistHalfBin);//50,0.,50.);
    hCriteriaHalfWayX[i]->Sumw2();
    sprintf(Name,"hCriteriaHalfWayY[%d]",i);
    hCriteriaHalfWayY[i] = new TH2D(Name,"",20,0,100,NBinsHalf,DistHalfBin);
    hCriteriaHalfWayY[i]->Sumw2();

    sprintf(Name,"hCriteriaAngleX[%d]",i);
    hCriteriaAngleX[i] = new TH2D(Name,"",20,0,100,30,0,90);
    hCriteriaAngleX[i]->Sumw2();
    sprintf(Name,"hCriteriaAngleY[%d]",i);
    hCriteriaAngleY[i] = new TH2D(Name,"",20,0,100,30,0,90);
    hCriteriaAngleY[i]->Sumw2();

    sprintf(Name,"hCriteriaHalfWayX_Data[%d]",i);
    hCriteriaHalfWayX_Data[i] = new TH2D(Name,"",20,0,100,NBinsHalf,DistHalfBin);
    hCriteriaHalfWayX_Data[i]->Sumw2();
    sprintf(Name,"hCriteriaHalfWayY_Data[%d]",i);
    hCriteriaHalfWayY_Data[i] = new TH2D(Name,"",20,0,100,NBinsHalf,DistHalfBin);
    hCriteriaHalfWayY_Data[i]->Sumw2();

    sprintf(Name,"hCriteriaAngleX_Data[%d]",i);
    hCriteriaAngleX_Data[i] = new TH2D(Name,"",20,0,100,30,0,90);
    hCriteriaAngleX_Data[i]->Sumw2();
    sprintf(Name,"hCriteriaAngleY_Data[%d]",i);
    hCriteriaAngleY_Data[i] = new TH2D(Name,"",20,0,100,30,0,90);
    hCriteriaAngleY_Data[i]->Sumw2();


    sprintf(Name,"CLMuon_Distribution[%d]",i);
    CLMuon_Distribution[i] = new TH1D(Name,"",200,-2,2);
    CLMuon_Distribution[i]->Sumw2();

    for(int a=0;a<3;a++){
      sprintf(Name,"CLMuon_Distribution_Particle[%d][%d]",i,a);
      CLMuon_Distribution_Particle[i][a] = new TH1D(Name,"",200,-2,2);
      CLMuon_Distribution_Particle[i][a]->Sumw2();
    }

    sprintf(Name,"CLMuon_Distribution_Data[%d]",i);
    CLMuon_Distribution_Data[i] = new TH1D(Name,"",200,-2,2);
    CLMuon_Distribution_Data[i]->Sumw2();

    sprintf(Name,"Selected_IronDistance_MC[%d]",i);
    Selected_IronDistance_MC[i] = new TH1D(Name,"",20,0,100);
    Selected_IronDistance_MC[i]->Sumw2();
  }
  for(int j=1;j<11;j++){
    sprintf(Name,"InFV_Interacting[%d]",j);
    InFV_Interacting[j] = new TH1D(Name,"",100,0,10);
    InFV_Interacting[j]->GetXaxis()->SetTitle("In FV interactions");

    sprintf(Name,"INGRIDTrackWidth[%d]",j);
    INGRIDTrackWidth[j] = new TH1D(Name,"Average track width distribution",40,1,5); 
    INGRIDTrackWidth[j]->GetXaxis()->SetTitle("Average transverse track width");
    INGRIDTrackWidth[j]->GetYaxis()->SetTitle("Number of tracks");

    sprintf(Name,"INGRIDAngleCriteria[%d]",j);
    INGRIDAngleCriteria[j] = new TH1D(Name,"PM and INGRID track angle difference",15,0,45); 
    INGRIDAngleCriteria[j]->GetXaxis()->SetTitle("Angle difference (#circ)");
    INGRIDAngleCriteria[j]->GetYaxis()->SetTitle("Number of tracks");

    sprintf(Name,"INGRIDPositionCriteria[%d]",j);
    INGRIDPositionCriteria[j] = new TH1D(Name,"PM and INGRID track transverse position difference",30,0,15); 
    INGRIDPositionCriteria[j]->GetXaxis()->SetTitle("Position difference (cm)");
    INGRIDPositionCriteria[j]->GetYaxis()->SetTitle("Number of tracks");


    sprintf(Name,"InFV_Reconstructed[%d]",j);
    InFV_Reconstructed[j] = new TH1D(Name,"",100,0,10);
    InFV_Reconstructed[j]->GetXaxis()->SetTitle("In FV interactions");

    sprintf(Name,"Reconstructed_NTracks_ATrackIsSample4[%d]",j);
    Reconstructed_NTracks_ATrackIsSample4[j] = new TH1D(Name,"",9,1,10);
    Reconstructed_NTracks_ATrackIsSample4[j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Reconstructed_NTracks_ATrackIsSample4[j]->GetYaxis()->SetTitle("Number of events");

    sprintf(Name,"Reconstructed_NTracks_ATrackIsSample5[%d]",j);
    Reconstructed_NTracks_ATrackIsSample5[j] = new TH1D(Name,"",9,1,10);
    Reconstructed_NTracks_ATrackIsSample5[j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Reconstructed_NTracks_ATrackIsSample5[j]->GetYaxis()->SetTitle("Number of events");

    sprintf(Name,"Selected_NTracks_ATrackIsSample4[%d]",j);
    Selected_NTracks_ATrackIsSample4[j] = new TH1D(Name,"",9,1,10);
    Selected_NTracks_ATrackIsSample4[j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Selected_NTracks_ATrackIsSample4[j]->GetYaxis()->SetTitle("Number of events");

    sprintf(Name,"Selected_NTracks_ATrackIsSample5[%d]",j);
    Selected_NTracks_ATrackIsSample5[j] = new TH1D(Name,"",9,1,10);
    Selected_NTracks_ATrackIsSample5[j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Selected_NTracks_ATrackIsSample5[j]->GetYaxis()->SetTitle("Number of events");


    for(int i=0;i<6;i++){
      sprintf(Name,"Selected_NTracks[%d][%d]",i,j);
      Selected_NTracks[i][j] = new TH1D(Name,"",9,1,10);
      Selected_NTracks[i][j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
      Selected_NTracks[i][j]->GetYaxis()->SetTitle("Number of events");

      //Selected_NTracks[i][j]->Sumw2();
      sprintf(Name,"Reconstructed_NTracks[%d][%d]",i,j);
      Reconstructed_NTracks[i][j] = new TH1D(Name,"",9,1,10);
      Reconstructed_NTracks[i][j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
      Reconstructed_NTracks[i][j]->GetYaxis()->SetTitle("Number of events");

      //Reconstructed_NTracks[i][j]->Sumw2();
      sprintf(Name,"Interacting_NTracks[%d][%d]",i,j);
      Interacting_NTracks[i][j] = new TH1D(Name,"",9,1,10);
      //Interacting_NTracks[i][j]->Sumw2();
      Interacting_NTracks[i][j]->GetXaxis()->SetTitle("Number of reconstructed tracks");
      Interacting_NTracks[i][j]->GetYaxis()->SetTitle("Number of events");

      sprintf(Name,"Selected_Momentum[%d][%d]",i,j);
      Selected_Momentum[i][j] = new TH1D(Name,"",200,0,10);

      sprintf(Name,"PID_Selected_IronDistance[%d][%d]",i,j);
      PID_Selected_IronDistance[i][j] = new TH1D(Name,"",20,0,100);
      sprintf(Name,"PIDAndMatching_Selected_IronDistance[%d][%d]",i,j);
      PIDAndMatching_Selected_IronDistance[i][j] = new TH1D(Name,"",20,0,100);

      sprintf(Name,"PID_Selected_ReconstructedAngle[%d][%d]",i,j);
      PID_Selected_ReconstructedAngle[i][j] = new TH1D(Name,"",30,0,90);
      sprintf(Name,"PIDAndMatching_Selected_ReconstructedAngle[%d][%d]",i,j);
      PIDAndMatching_Selected_ReconstructedAngle[i][j] = new TH1D(Name,"",30,0,90);


      //Selected_Momentum[i][j]->Sumw2();
      sprintf(Name,"Reconstructed_Momentum[%d][%d]",i,j);
      Reconstructed_Momentum[i][j] = new TH1D(Name,"",200,0,10);
      //Reconstructed_Momentum[i][j]->Sumw2();
      sprintf(Name,"Interacting_Momentum[%d][%d]",i,j);
      Interacting_Momentum[i][j] = new TH1D(Name,"",200,0,10);
      //Interacting_Momentum[i][j]->Sumw2();

      sprintf(Name,"Selected_IronDistance[%d][%d]",i,j);
      Selected_IronDistance[i][j] = new TH1D(Name,"",20,0,100);
   
      sprintf(Name,"Reconstructed_IronDistance[%d][%d]",i,j);
      Reconstructed_IronDistance[i][j] = new TH1D(Name,"",20,0,100);
      //Reconstructed_IronDistance[i][j]->Sumw2();
      sprintf(Name,"Interacting_IronDistance[%d][%d]",i,j);
      Interacting_IronDistance[i][j] = new TH1D(Name,"",20,0,100);
      //Interacting_IronDistance[i][j]->Sumw2();

      sprintf(Name,"Selected_Angle[%d][%d]",i,j);
      Selected_Angle[i][j] = new TH1D(Name,"",30,0,90);
      //Selected_Angle[i][j]->Sumw2();
      sprintf(Name,"Reconstructed_Angle[%d][%d]",i,j);
      Reconstructed_Angle[i][j] = new TH1D(Name,"",30,0,90);
      //Reconstructed_Angle[i][j]->Sumw2();
      sprintf(Name,"Interacting_Angle[%d][%d]",i,j);
      Interacting_Angle[i][j] = new TH1D(Name,"",30,0,90);
      //Interacting_Angle[i][j]->Sumw2();

      sprintf(Name,"Selected_ReconstructedAngle[%d][%d]",i,j);
      Selected_ReconstructedAngle[i][j] = new TH1D(Name,"",30,0,90);
      //Selected_ReconstructedAngle[i][j]->Sumw2();
      sprintf(Name,"Reconstructed_ReconstructedAngle[%d][%d]",i,j);
      Reconstructed_ReconstructedAngle[i][j] = new TH1D(Name,"",30,0,90);
      //Reconstructed_ReconstructedAngle[i][j]->Sumw2();
      sprintf(Name,"Interacting_ReconstructedAngle[%d][%d]",i,j);
      Interacting_ReconstructedAngle[i][j] = new TH1D(Name,"",30,0,90);
      //Interacting_ReconstructedAngle[i][j]->Sumw2();

    }
  }

  for(int type=0;type<4;type++){
    sprintf(Name,"Selected_Higher[%d]",type);
    Selected_Higher[type] = new TH1D(Name,"",11,0,1.1);

    sprintf(Name,"Interacting_Higher[%d]",type);
    Interacting_Higher[type] = new TH1D(Name,"",11,0,1.1);

    sprintf(Name,"Selected_Lower[%d]",type);
    Selected_Lower[type] = new TH1D(Name,"",11,0,1.1);

    sprintf(Name,"Interacting_Lower[%d]",type);
    Interacting_Lower[type] = new TH1D(Name,"",11,0,1.1);
  }

  double weightMC;
  double CountSand=0;
  cout<<"REMOVE THE FV CUT WHEN EFFICIENCY IS FINISHED, AND REMOVE ALSO THE FIRST DETECTED CUT"<<endl;

  for(int ievt=0;ievt<nevt;ievt++){
    if(ievt%10000==0) cout<<ievt<<endl;
    wtree->GetEntry(ievt);
    if(IsSand) FSIInt=7;
    else if(IsBkgH) FSIInt=8;
    else if(IsBkgH) FSIInt=9;
    else if(IsNuE || IsAnti) FSIInt=10;

    //if(FSIInt>7) continue;
    ///if(FSIInt==5) cout<<FSIInt<<endl;    
   
    if(weight>1e2) continue;
    weightMC+=weight;
    //cout<<FSIInt<<endl;
    if(FSIInt==2) FSIInt=1;
    //if(IsFV){
     
    if(IsFV) InFV_Interacting[FSIInt]->Fill(Enu,weight);
    bool FirstRec=true;
    for(int irec=0;irec<nIngBasRec;irec++){
      bool Trash=false;
      //cout<<"new rec"<<endl;
      bool MuonFound=false;
      int MuonTrue;int MuonRec=0;
      //if(weight>1e3) continue;//cout<<weight<<endl;
      CountSand++;//=weight; 
      //weight*=1.1;
      //if(IsSand) cout<<FSIInt<<", "<<Sample[irec][MuonRec]<<endl;
      //cout<<"ok"<<endl;
      if(IsDetected[irec]){
	if(IsFV && FirstRec){
	  InFV_Reconstructed[FSIInt]->Fill(Enu,weight);
	  //FirstRec=false;
	} 
	//if(FSIInt==5) cout<<"yeah="<<FSIInt<<endl;    

	//IronDistance[irec][MuonRec]=IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio);
	//if(VertexZ[irec]>10) continue;
	//if(IsAnti) cout<<"IsAnti, weight="<<weight<<endl;
	//else cout<<weight<<endl;
	//cout<<"IsFV, interaction="<<FSIInt<<endl;
	if(nTracks[irec]>20) nTracks[irec]=20;
	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  //cout<<"before test"<<endl;
	  if(IsReconstructed[irec][itrk]){
	    IronDistance[irec][itrk]=IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio);
	    double Test=r1->Uniform();
	    if(Sample[irec][itrk]==5 && IronDistance[irec][itrk]>58.5 && Test>0.95) Sample[irec][itrk]=4;
 
	    //if(nTracks[irec]==1) cout<<"Track sample="<<Sample[irec][itrk]<<endl;
	    if(CLMuon[irec][itrk]!=CLMuon[irec][itrk] || CLMuon[irec][itrk]==-1){ Trash=true;continue;}
	    if(TypeOfTrack[irec][itrk]==0){
	      MuonTrue=itrk;
	      if(Sample[irec][itrk]==4) MomentumXIron_ReconstructedForTrueMuon->Fill(TrueMomentumMuon,IronDistance[irec][itrk],weight);
	      if(Sample[irec][itrk]==5) AngleRecXTrue->Fill(TrueAngleMuon,TrackAngle[irec][itrk],weight);
	    }
	    //cout<<CLMuon[irec][itrk]<<", "<<CLMuon[irec][MuonRec]<<endl;
	    if(CLMuon[irec][itrk]>=CLMuon[irec][MuonRec]) MuonRec=itrk;
	    //cout<<Sample[irec][itrk]<<", "<<FSIInt<<", irec="<<irec<<", MuonRec="<<MuonRec<<endl;
	    MuonFound=true;
	  }
	}
	if(Trash) continue; 
	NTracksTotal->Fill(nTracks[irec],weight);
	if(!MuonFound) continue;
	//if(IsBkgH) cout<<"Is Bkg Horiz, weight="<<weight<<endl;
	//if(IsBkgV) cout<<"Is Bkg Verti, weight="<<weight<<endl;
	//if(IsNuE) cout<<"Is NuE, weight="<<weight<<endl;
	bool FirstSample4=true;
	bool FirstSample5=true;
	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  if(IsReconstructed[irec][itrk]==true){
	    if(Sample[irec][itrk]>=4){
	      if(IsNuE) TrackWidthXIronDistance_NuE->Fill(IronDistance[irec][itrk],TrackWidth[irec][itrk],weight);
	      else if(!IsNuE && !IsSand && !IsAnti && !IsBkgH && !IsBkgV) TrackWidthXIronDistance->Fill(IronDistance[irec][itrk],TrackWidth[irec][itrk],weight);
	    }
	    if(Sample[irec][itrk]==5) AngleTotal->Fill(TrackAngle[irec][itrk],weight);
	    CLMuonTotal->Fill(CLMuon[irec][itrk],weight);

	    if(TypeOfTrack[irec][itrk]<3) CLMuon_Distribution_Particle[Sample[irec][itrk]][TypeOfTrack[irec][itrk]]->Fill(CLMuon[irec][itrk],weight);
	    CLMuon_Distribution[Sample[irec][itrk]]->Fill(CLMuon[irec][itrk],weight);
	    SampleTotal->Fill(Sample[irec][itrk],weight);

	    if(Sample[irec][itrk]==5){
	      if(FirstSample5){
		if(IsFV && FirstRec) Reconstructed_NTracks_ATrackIsSample5[FSIInt]->Fill(nTracks[irec],weight);
		Selected_NTracks_ATrackIsSample5[FSIInt]->Fill(nTracks[irec],weight);
		FirstSample5=false;
	      }
	    }
	    if(Sample[irec][itrk]==4){
	      if(FirstSample4){
		if(IsFV && FirstRec){
		  Reconstructed_NTracks_ATrackIsSample4[FSIInt]->Fill(nTracks[irec],weight);
		  /*		    if(FSIInt==1){
		    cout<<"Is a track, weight="<<weight<<endl;
		    cout<<"Integral="<<Reconstructed_NTracks_ATrackIsSample4[FSIInt]->Integral();
		    }
		  */
		}
		Selected_NTracks_ATrackIsSample4[FSIInt]->Fill(nTracks[irec],weight);
		FirstSample4=false;
	      }
	      IronDistanceTotal->Fill(IronDistance[irec][itrk],weight);
	      NTracksTotalSample4->Fill(nTracks[irec],weight);
	      if(nTracks[irec]==1) IronDistanceTotalSample4_1Track->Fill(IronDistance[irec][itrk],weight);
	      if(nTracks[irec]==2){
		if(itrk==MuonRec) IronDistanceTotalSample4_2Tracks->Fill(IronDistance[irec][MuonRec],weight);
		for(int itrk2=0;itrk2<nTracks[irec];itrk2++){
		  if(itrk==itrk2) continue;
		  SampleTotalSample4->Fill(Sample[irec][itrk2],weight);
		}
	      }
	    }
	  }
	}

	//if(Sample[irec][MuonRec]==4) IronDistanceTotal->Fill(IronDistance[irec][MuonRec],weight);
	//if(Sample[irec][MuonRec]==5) AngleTotal->Fill(TrackAngle[irec][MuonRec],weight);
	//if(IronDistance[irec][MuonRec]>=58.5){Sample[irec][MuonRec]=4; cout<<"Careful, change this value for sample through going plots!!!!!!!"<<endl;}
	//cout<<"ntracks="<<nTracks[irec]<<", test="<<Sample[irec][MuonRec]<<", "<<FSIInt<<", irec="<<irec<<", MuonRec="<<MuonRec<<endl;

	if(IsFV && FirstRec) Reconstructed_NTracks[Sample[irec][MuonRec]][FSIInt]->Fill(nTracks[irec],weight);
	Reconstructed_Momentum[Sample[irec][MuonRec]][FSIInt]->Fill(TrueMomentumMuon,weight);
	Reconstructed_IronDistance[Sample[irec][MuonRec]][FSIInt]->Fill(IronDistance[irec][MuonRec],weight);
	Reconstructed_Angle[Sample[irec][MuonRec]][FSIInt]->Fill(TrueAngleMuon,weight);
	Reconstructed_ReconstructedAngle[Sample[irec][MuonRec]][FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	   
	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  if(IsReconstructed[irec][itrk]==true){
	   
	    for(double mucl=0;mucl<=1;mucl=mucl+0.1){
	      if(CLMuon[irec][itrk]>=mucl) Selected_Higher[TypeOfTrack[irec][itrk]]->Fill(mucl+1e-3,weight);
	      else if(CLMuon[irec][itrk]<=mucl){
		//IsReconstructed[irec][itrk]=0;
		//cout<<"type="<<TypeOfTrack[irec][itrk]<<", is rec="<<IsReconstructed[irec][itrk]<<endl;
		Selected_Lower[TypeOfTrack[irec][itrk]]->Fill(mucl+1e-3,weight);
		//if(mucl==.7) cout<<"passing"<<endl;
	      }
	      //if(mucl==1) cout<<"interacting"<<endl<<endl;
	      Interacting_Higher[TypeOfTrack[irec][itrk]]->Fill(mucl+1e-3,weight);
	      Interacting_Lower[TypeOfTrack[irec][itrk]]->Fill(mucl+1e-3,weight);
	    }

	    if(Sample[irec][itrk]==4){
	      if(IronDistance[irec][itrk]>10){
		int BinIron=(int) (IronDistance[irec][itrk]/5)-2;
		if(IronDistance[irec][itrk]>90) BinIron=NBinsIron-1;
		hCriteriaAngleXY[BinIron]->Fill(CriteriaAngleX[irec][itrk],CriteriaAngleY[irec][itrk],weight);
		hCriteriaHalfWayXY[BinIron]->Fill(CriteriaHalfWayX[irec][itrk],CriteriaHalfWayY[irec][itrk],weight);
	      }
	    }
	    hCriteriaHalfWayX[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaHalfWayX[irec][itrk],weight);
	    hCriteriaHalfWayY[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaHalfWayY[irec][itrk],weight);
	    hCriteriaAngleX[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaAngleX[irec][itrk],weight);
	    hCriteriaAngleY[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaAngleY[irec][itrk],weight);
	    

	  }
	}
	
	if((nTracks[irec]==1 && CLMuon[irec][0]>0.7) || ((nTracks[irec]==2 && TMath::Max(CLMuon[irec][0],CLMuon[irec][1])>0.7 && TMath::Min(CLMuon[irec][0],CLMuon[irec][1])<0.3) && (TMath::Min(CLMuon[irec][0],CLMuon[irec][1])>=0)) || ((nTracks[irec]==3 && CLMuon[irec][MuonRec]>0.7 && CLMuon[irec][(MuonRec+1)%3]<0.3 && CLMuon[irec][(MuonRec+2)%3]<0.3)) || ((nTracks[irec]==4 && CLMuon[irec][MuonRec]>0.7 && CLMuon[irec][(MuonRec+1)%4]<0.3 && CLMuon[irec][(MuonRec+2)%4]<0.3 && CLMuon[irec][(MuonRec+3)%4]<0.3))){
	  //if(FSIInt==1 && Sample[irec][MuonRec]==4) cout<<"The interaction is Selected, weight="<<weight<<endl;

	    
	  //if(Sample[irec][MuonRec]==4 && PlasticDistance[irec][MuonRec]<25) continue;
	  //if(Sample[irec][MuonRec]>=4) TrackWidthXIronDistance->Fill(IronDistance[irec][MuonRec],TrackWidth[irec][MuonRec],weight);

	  if(Sample[irec][MuonRec]==4 && IronDistance[irec][MuonRec]<17){
	    Problematic_PlasticDistance->Fill(PlasticDistance[irec][MuonRec],weight);
	    Problematic_Angle->Fill(TrackAngle[irec][MuonRec],weight);
	    Problematic_NTracks->Fill(nTracks[irec],weight);
	    Problematic_CLMuon->Fill(CLMuon[irec][MuonRec],weight);
	    Problematic_TrackWidth->Fill(TrackWidth[irec][MuonRec],weight);
	    Problematic_Geom->Fill(GeometricTrack[irec][MuonRec],weight);
	    int LS;
	    if(TMath::Abs(LastChannelINGRIDX[irec][MuonRec]-11.5)>TMath::Abs(LastChannelINGRIDY[irec][MuonRec]-11.5)) LS=LastChannelINGRIDX[irec][MuonRec];
	    else LS=LastChannelINGRIDY[irec][MuonRec];
	    Problematic_LastChannel->Fill(LS,weight);

	  }
	  //if(Sample[irec][MuonRec]==4 /*&& IronDistance[irec][MuonRec]<25*/ && (CriteriaAngleX[irec][MuonRec]>25 || CriteriaAngleY[irec][MuonRec]>25 || CriteriaHalfWayX[irec][MuonRec]>7 || CriteriaHalfWayY[irec][MuonRec]>7 || TrackWidth[irec][MuonRec]>1.1)) continue; /*
	  //if(Sample[irec][MuonRec]==4 && TrackWidth[irec][MuonRec]>1.1 && IronDistance[irec][MuonRec]<25) continue;	      
	  //if(Sample[irec][MuonRec]==4 && TrackWidth[irec][MuonRec]>1.5 && IronDistance[irec][MuonRec]>25) continue;	 
	  PID_Selected_IronDistance[Sample[irec][MuonRec]][FSIInt]->Fill(IronDistance[irec][MuonRec],weight);
	  PID_Selected_ReconstructedAngle[Sample[irec][MuonRec]][FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);

	  if(Sample[irec][MuonRec]==4){
	    INGRIDAngleCriteria[FSIInt]->Fill(max(CriteriaAngleX[irec][MuonRec],CriteriaAngleY[irec][MuonRec]),weight);
	    INGRIDPositionCriteria[FSIInt]->Fill(max(CriteriaHalfWayX[irec][MuonRec],CriteriaHalfWayY[irec][MuonRec]),weight);
	  }


	  if(Sample[irec][MuonRec]>=4 && (CriteriaAngleX[irec][MuonRec]>25 || CriteriaAngleY[irec][MuonRec]>25 || CriteriaHalfWayX[irec][MuonRec]>7 || CriteriaHalfWayY[irec][MuonRec]>7)) continue; 	      
	  PIDAndMatching_Selected_IronDistance[Sample[irec][MuonRec]][FSIInt]->Fill(IronDistance[irec][MuonRec],weight);
	  PIDAndMatching_Selected_ReconstructedAngle[Sample[irec][MuonRec]][FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);

	 
	  if(Sample[irec][MuonRec]>=4){
	    if(Sample[irec][MuonRec]==4) INGRIDTrackWidth[FSIInt]->Fill(TrackWidth[irec][MuonRec],weight);
	    if(TrackWidth[irec][MuonRec]>1.1 && IronDistance[irec][MuonRec]<15) continue;
	    else if(TrackWidth[irec][MuonRec]>1.2 && IronDistance[irec][MuonRec]<20) continue;
	    else if(TrackWidth[irec][MuonRec]>1.3 && IronDistance[irec][MuonRec]<30) continue;
	    else if(TrackWidth[irec][MuonRec]>1.5) continue;
	  }
	      
	  //else if(Sample[irec][MuonRec]==4 && TrackWidth[irec][MuonRec]>1.1 && IronDistance[irec][MuonRec]<15) continue;	      
	  //if(Sample[irec][MuonRec]==4 && TrackWidth[irec][MuonRec]>1.5 && IronDistance[irec][MuonRec]>25) continue;	  
	  if(Sample[irec][MuonRec]==5) AngleRecXTrue_Selected->Fill(TrueAngleMuon,TrackAngle[irec][MuonRec],weight);
	  //if(FSIInt==5) cout<<weight<<endl;
	  //if(FSIInt==1 && Sample[irec][MuonRec]==4) cout<<"The interaction is finally Selected, weight="<<weight<<endl;

	  Selected_NTracks[Sample[irec][MuonRec]][FSIInt]->Fill(nTracks[irec],weight);
	  /*	      if(FSIInt==1 && Sample[irec][MuonRec]==4){
	    cout<<"The interaction is finally Selected, weight="<<weight<<endl;
	    cout<<"Integral="<< Selected_NTracks[Sample[irec][MuonRec]][FSIInt]->Integral();
	    }
	  */
	  if(Sample[irec][MuonRec]==4) Selected2D->Fill(IronDistance[irec][MuonRec],TrackAngle[irec][MuonRec],weight);

	  Selected_Momentum[Sample[irec][MuonRec]][FSIInt]->Fill(TrueMomentumMuon,weight);
	  Selected_IronDistance[Sample[irec][MuonRec]][FSIInt]->Fill(IronDistance[irec][MuonRec],weight);
	  Selected_IronDistance_MC[Sample[irec][MuonRec]]->Fill(IronDistance[irec][MuonRec],weight);
	  Selected_Angle[Sample[irec][MuonRec]][FSIInt]->Fill(TrueAngleMuon,weight);
	  Selected_ReconstructedAngle[Sample[irec][MuonRec]][FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	}
	//if(IsFV){
	//Interacting_NuE[Sample[irec][MuonRec]][FSIInt]->Fill(Enu,weight);
	Interacting_NTracks[Sample[irec][MuonRec]][FSIInt]->Fill(nTracks[irec],weight);
	Interacting_Momentum[Sample[irec][MuonRec]][FSIInt]->Fill(TrueMomentumMuon,weight);
	Interacting_IronDistance[Sample[irec][MuonRec]][FSIInt]->Fill(IronDistance[irec][MuonRec],weight);
	Interacting_Angle[Sample[irec][MuonRec]][FSIInt]->Fill(TrueAngleMuon,weight);
	Interacting_ReconstructedAngle[Sample[irec][MuonRec]][FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	//}
	if(FirstRec) FirstRec=false;
      } 
	 
    }
  }
  //}
  
  cout<<"Number of sands="<<CountSand<<endl;
  TChain * wtreeData = new TChain("wtree");
  for(int i=13000;i<17265;i++){
    //for(int i=13000;i<17500;i++){
    //for(int i=14510;i<14511;i++){
    if(i>=14337 && i<=14429) continue;
    /*
      if(i<13048) continue;
      else{
      if(i>13094 && i<13107) continue;
      else if(i>13159 && i<13167) continue;
      else if(i>13175 && i<14510) continue;
      else if(i<16901 && i>14570) continue;
      }
    */
    //if(i>=16500 && i<=16900) continue;
    //if(i==14500) continue;
    sprintf(Name,"AllCC0pi5/CC0piTree_%d.root",i);
    wtreeData->Add(Name);
  }  
  
  /*
    for(int i=1;i<=100;i++){
    //sprintf(Name,"CC0piData6/CC0piTree%d.root",i);
    sprintf(Name,"PEError/CC0piTree_%d_PE1.root",i);
    cout<<Name<<endl;
    wtreeData->Add(Name);
    }
  */
  //wtreeData->Add("AllCC0piData.root");
  int nevtData=(int) wtreeData->GetEntries();

  TBranch* Br_Spill = wtreeData->GetBranch("Spill");
  Br_Spill->SetAddress(&Spill);
  wtreeData->SetBranchAddress("Spill",&Spill);

  TBranch* Br_GoodSpill = wtreeData->GetBranch("GoodSpill");
  Br_GoodSpill->SetAddress(&GoodSpill);
  wtreeData->SetBranchAddress("GoodSpill",&GoodSpill);

  TBranch* Br_FSIInt = wtreeData->GetBranch("FSIInt");
  Br_FSIInt->SetAddress(&FSIInt);
  wtreeData->SetBranchAddress("FSIInt",&FSIInt);

  TBranch* Br_nIngBasRec = wtreeData->GetBranch("nIngBasRec");
  Br_nIngBasRec->SetAddress(&nIngBasRec);
  wtreeData->SetBranchAddress("nIngBasRec",&nIngBasRec);

  TBranch * Br_Num_Int = wtreeData->GetBranch("InteractionType");
  Br_Num_Int->SetAddress(&Num_Int);
  wtreeData->SetBranchAddress("InteractionType",&Num_Int);
  
  TBranch * Br_nTracks = wtreeData->GetBranch("nTracks[10]");
  Br_nTracks->SetAddress(&nTracks);
  wtreeData->SetBranchAddress("nTracks[10]",&nTracks);
 
  TBranch * Br_weight = wtreeData->GetBranch("weight");
  Br_weight->SetAddress(&weight);
  wtreeData->SetBranchAddress("weight",&weight);

  TBranch * Br_Enu = wtreeData->GetBranch("Enu");
  Br_Enu->SetAddress(&Enu);
  wtreeData->SetBranchAddress("Enu",&Enu);

  TBranch * Br_TrueMomentumMuon = wtreeData->GetBranch("TrueMomentumMuon");
  Br_TrueMomentumMuon->SetAddress(&TrueMomentumMuon);
  wtreeData->SetBranchAddress("TrueMomentumMuon",&TrueMomentumMuon);

  TBranch * Br_TrueAngleMuon = wtreeData->GetBranch("TrueAngleMuon");
  Br_TrueAngleMuon->SetAddress(&TrueAngleMuon);
  wtreeData->SetBranchAddress("TrueAngleMuon",&TrueAngleMuon);

  TBranch * Br_IsFV = wtreeData->GetBranch("IsFV");
  Br_IsFV->SetAddress(&IsFV);
  wtreeData->SetBranchAddress("IsFV",&IsFV);

  TBranch * Br_POT = wtreeData->GetBranch("POT");
  Br_POT->SetAddress(&POT);
  wtreeData->SetBranchAddress("POT",&POT);
  /*
    TBranch * Br_VertexZ = wtreeData->GetBranch("VertexZ[10]");
    Br_VertexZ->SetAddress(&VertexZ);
    wtreeData->SetBranchAddress("VertexZ[10]",&VertexZ);
  */
  TBranch * Br_IsDetected = wtreeData->GetBranch("IsDetected[10]");
  Br_IsDetected->SetAddress(&IsDetected);
  wtreeData->SetBranchAddress("IsDetected[10]",&IsDetected);
  
  TBranch * Br_TrackAngle = wtreeData->GetBranch("TrackAngle[10][20]");
  Br_TrackAngle->SetAddress(TrackAngle);
  wtreeData->SetBranchAddress("TrackAngle[10][20]",TrackAngle);
 
  TBranch * Br_TypeOfTrack = wtreeData->GetBranch("TypeOfTrack[10][20]");
  Br_TypeOfTrack->SetAddress(TypeOfTrack);
  wtreeData->SetBranchAddress("TypeOfTrack[10][20]",TypeOfTrack);

  TBranch * Br_CLMuon = wtreeData->GetBranch("CLMuon[10][20]");
  Br_CLMuon->SetAddress(&CLMuon);
  wtreeData->SetBranchAddress("CLMuon[10][20]",&CLMuon);

  TBranch * Br_Momentum = wtreeData->GetBranch("Momentum[10][20]");
  Br_Momentum->SetAddress(Momentum);
  wtreeData->SetBranchAddress("Momentum[10][20]",Momentum);

  TBranch * Br_IronDistance = wtreeData->GetBranch("IronDistance[10][20]");
  Br_IronDistance->SetAddress(IronDistance);
  wtreeData->SetBranchAddress("IronDistance[10][20]",IronDistance);

  TBranch * Br_PlasticDistance = wtreeData->GetBranch("PlasticDistance[10][20]");
  Br_PlasticDistance->SetAddress(PlasticDistance);
  wtreeData->SetBranchAddress("PlasticDistance[10][20]",PlasticDistance);

  TBranch * Br_Sample = wtreeData->GetBranch("Sample[10][20]");
  Br_Sample->SetAddress(Sample);
  wtreeData->SetBranchAddress("Sample[10][20]",Sample);

  TBranch * Br_IsReconstructed = wtreeData->GetBranch("IsReconstructed[10][20]");
  Br_IsReconstructed->SetAddress(IsReconstructed);
  wtreeData->SetBranchAddress("IsReconstructed[10][20]",IsReconstructed);

  TBranch * Br_CriteriaAngleX = wtreeData->GetBranch("CriteriaAngleX[10][20]");
  Br_CriteriaAngleX->SetAddress(CriteriaAngleX);
  wtreeData->SetBranchAddress("CriteriaAngleX[10][20]",CriteriaAngleX);
  TBranch * Br_CriteriaAngleY = wtreeData->GetBranch("CriteriaAngleY[10][20]");
  Br_CriteriaAngleY->SetAddress(CriteriaAngleY);
  wtreeData->SetBranchAddress("CriteriaAngleY[10][20]",CriteriaAngleY);
  TBranch * Br_CriteriaHalfWayX = wtreeData->GetBranch("CriteriaHalfWayX[10][20]");
  Br_CriteriaHalfWayX->SetAddress(CriteriaHalfWayX);
  wtreeData->SetBranchAddress("CriteriaHalfWayX[10][20]",CriteriaHalfWayX);
  TBranch * Br_CriteriaHalfWayY = wtreeData->GetBranch("CriteriaHalfWayY[10][20]");
  Br_CriteriaHalfWayY->SetAddress(CriteriaHalfWayY);
  wtreeData->SetBranchAddress("CriteriaHalfWayY[10][20]",CriteriaHalfWayY);

  TBranch * Br_ReWeight = wtreeData->GetBranch("ReWeight[175]");
  Br_ReWeight->SetAddress(ReWeight);
  wtreeData->SetBranchAddress("ReWeight[175]",ReWeight);

  TBranch * Br_LastChannelINGRIDX = wtreeData->GetBranch("LastChannelINGRIDX[10][20]");
  Br_LastChannelINGRIDX->SetAddress(LastChannelINGRIDX);
  wtreeData->SetBranchAddress("LastChannelINGRIDX[10][20]",LastChannelINGRIDX);

  TBranch * Br_LastChannelINGRIDY = wtreeData->GetBranch("LastChannelINGRIDY[10][20]");
  Br_LastChannelINGRIDY->SetAddress(LastChannelINGRIDY);
  wtreeData->SetBranchAddress("LastChannelINGRIDY[10][20]",LastChannelINGRIDY);

  TBranch * Br_TrackWidth = wtreeData->GetBranch("TrackWidth[10][20]");
  Br_TrackWidth->SetAddress(TrackWidth);
  wtreeData->SetBranchAddress("TrackWidth[10][20]",TrackWidth);

  TBranch * Br_GeometricTrack = wtreeData->GetBranch("GeometricTrack[10][20]");
  Br_GeometricTrack->SetAddress(GeometricTrack);
  wtreeData->SetBranchAddress("GeometricTrack[10][20]",GeometricTrack);

  INGRIDTrackWidth_Data = new TH1D("INGRIDTrackWidth_Data","Average track width distribution",40,1,5); 
  INGRIDTrackWidth_Data->GetXaxis()->SetTitle("Average transverse track width");
  INGRIDTrackWidth_Data->GetYaxis()->SetTitle("Number of tracks");

  INGRIDAngleCriteria_Data = new TH1D(Name,"PM and INGRID track angle difference",15,0,45); 
  INGRIDAngleCriteria_Data->GetXaxis()->SetTitle("Angle difference (#circ)");
  INGRIDAngleCriteria_Data->GetYaxis()->SetTitle("Number of tracks");

  INGRIDPositionCriteria_Data = new TH1D(Name,"PM and INGRID track transverse position difference",30,0,15); 
  INGRIDPositionCriteria_Data->GetXaxis()->SetTitle("Position difference (cm)");
  INGRIDPositionCriteria_Data->GetYaxis()->SetTitle("Number of tracks");


  sprintf(Name,"Selected_NTracks_ATrackIsSample4_Data");
  Selected_NTracks_ATrackIsSample4_Data = new TH1D(Name,"",9,1,10);
  Selected_NTracks_ATrackIsSample4_Data->GetXaxis()->SetTitle("Number of reconstructed tracks");
  Selected_NTracks_ATrackIsSample4_Data->GetYaxis()->SetTitle("Number of events");

  sprintf(Name,"Selected_NTracks_ATrackIsSample5_Data");
  Selected_NTracks_ATrackIsSample5_Data = new TH1D(Name,"",9,1,10);
  Selected_NTracks_ATrackIsSample5_Data->GetXaxis()->SetTitle("Number of reconstructed tracks");
  Selected_NTracks_ATrackIsSample5_Data->GetYaxis()->SetTitle("Number of events");

  for(int i=0;i<6;i++){
    sprintf(Name,"Selected_NTracks_Data[%d]",i);
    Selected_NTracks_Data[i] = new TH1D(Name,"",9,1,10);
    Selected_NTracks_Data[i]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Selected_NTracks_Data[i]->GetYaxis()->SetTitle("Number of events");
    //Selected_NTracks_Data[i]->Sumw2();
    sprintf(Name,"Reconstructed_NTracks_Data[%d]",i);
    Reconstructed_NTracks_Data[i] = new TH1D(Name,"",9,1,10);
    Reconstructed_NTracks_Data[i]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Reconstructed_NTracks_Data[i]->GetYaxis()->SetTitle("Number of events");
    //Reconstructed_NTracks_Data[i]->Sumw2();
    sprintf(Name,"Interacting_NTracks_Data[%d]",i);
    Interacting_NTracks_Data[i] = new TH1D(Name,"",9,1,10);
    //Interacting_NTracks_Data[i]->Sumw2();
    Interacting_NTracks_Data[i]->GetXaxis()->SetTitle("Number of reconstructed tracks");
    Interacting_NTracks_Data[i]->GetYaxis()->SetTitle("Number of events");

    sprintf(Name,"Selected_IronDistance_Data[%d]",i);
    Selected_IronDistance_Data[i] = new TH1D(Name,"",20,0,100);
   
    sprintf(Name,"PID_Selected_IronDistance_Data[%d]",i);
    PID_Selected_IronDistance_Data[i] = new TH1D(Name,"",20,0,100);
    sprintf(Name,"PIDAndMatching_Selected_IronDistance_Data[%d]",i);
    PIDAndMatching_Selected_IronDistance_Data[i] = new TH1D(Name,"",20,0,100);

    sprintf(Name,"PID_Selected_ReconstructedAngle_Data[%d]",i);
    PID_Selected_ReconstructedAngle_Data[i] = new TH1D(Name,"",30,0,90);
    sprintf(Name,"PIDAndMatching_Selected_ReconstructedAngle_Data[%d]",i);
    PIDAndMatching_Selected_ReconstructedAngle_Data[i] = new TH1D(Name,"",30,0,90);

    //Selected_IronDistance_Data[i]->Sumw2();
    sprintf(Name,"Reconstructed_IronDistance_Data[%d]",i);
    Reconstructed_IronDistance_Data[i] = new TH1D(Name,"",20,0,100);
    //Reconstructed_IronDistance_Data[i]->Sumw2();
    sprintf(Name,"Interacting_IronDistance_Data[%d]",i);
    Interacting_IronDistance_Data[i] = new TH1D(Name,"",20,0,100);
    //Interacting_IronDistance_Data[i]->Sumw2();

    sprintf(Name,"Selected_ReconstructedAngle_Data_Data[%d]",i);
    Selected_ReconstructedAngle_Data[i] = new TH1D(Name,"",30,0,90);
    //Selected_ReconstructedAngle_Data[i]->Sumw2();
    sprintf(Name,"Reconstructed_ReconstructedAngle_Data[%d]",i);
    Reconstructed_ReconstructedAngle_Data[i] = new TH1D(Name,"",30,0,90);
    //Reconstructed_ReconstructedAngle_Data[i]->Sumw2();
    sprintf(Name,"Interacting_ReconstructedAngle_Data[%d]",i);
    Interacting_ReconstructedAngle_Data[i] = new TH1D(Name,"",30,0,90);
    //Interacting_ReconstructedAngle_Data[i]->Sumw2();
  }

  double POTCount=0;
  for(int ievt=0;ievt<nevtData;ievt++){
    if(ievt%10000==0) cout<<ievt<<endl;
    wtreeData->GetEntry(ievt);
    if(!GoodSpill || !Spill) continue;
    //if(IsSand) FSIInt=7;
    //else if(IsBkgH) FSIInt=8;
    //else if(IsBkgH) FSIInt=9;
    //else if(IsNuE || IsAnti) FSIInt=10;
    //if(FSIInt>7) continue;

    if(weight>1e2) continue;

    POTCount+=POT;
    if(ievt%10000==0) cout<<POTCount<<endl;
    if(FSIInt==2) FSIInt=1;
    
    //if(IsFV){
    //if(nIngBasRec!=0) cout<<nIngBasRec<<endl;
    for(int irec=0;irec<nIngBasRec;irec++){
      bool MuonFound=false;
      int MuonTrue;int MuonRec=0;
      bool Trash=false;
      if(IsDetected[irec]){
	//if(VertexZ[irec]>10) continue;

	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  //cout<<"before test"<<endl;
	  if(IsReconstructed[irec][itrk]){
	    IronDistance[irec][itrk]=IronDistance[irec][itrk]+(PlasticDistance[irec][itrk]/IronCarbonRatio);
	    double Test=r1->Uniform();
	    if(Sample[irec][itrk]==5 && IronDistance[irec][itrk]>58.5 && Test>0.95) Sample[irec][itrk]=4;

	    //if(nTracks[irec]==1) cout<<"Track sample="<<Sample[irec][itrk]<<endl;
	    //if(CLMuon[irec][itrk]!=CLMuon[irec][itrk]) continue;
	    if(CLMuon[irec][itrk]!=CLMuon[irec][itrk] || CLMuon[irec][itrk]==-1){ Trash=true; continue;}
	    if(TypeOfTrack[irec][itrk]==0) MuonTrue=itrk;
	    if(CLMuon[irec][itrk]>=CLMuon[irec][MuonRec]) MuonRec=itrk;
	    //cout<<"Muon Rec="<<MuonRec<<endl;
	    MuonFound=true;
	  }
	}

	if(Trash) continue;
	bool FirstSample4=true;
	bool FirstSample5=true;

	//cout<<CLMuon[irec][MuonRec]<<endl;
	NTracksTotal_Data->Fill(nTracks[irec],weight);
	if(!MuonFound) continue;
	for(int itrk=0;itrk<nTracks[irec];itrk++){
	  if(IsReconstructed[irec][itrk]==true){
	    if(Sample[irec][itrk]==5) AngleTotal_Data->Fill(TrackAngle[irec][itrk],weight);
	    CLMuonTotal_Data->Fill(CLMuon[irec][itrk],weight);
	    //if(Sample[irec][itrk]<4) cout<<CLMuon[irec][itrk]<<", weight="<<weight<<endl;
	    CLMuon_Distribution_Data[Sample[irec][itrk]]->Fill(CLMuon[irec][itrk],weight);
	    SampleTotal_Data->Fill(Sample[irec][itrk],weight);

	    if(Sample[irec][itrk]==4){
	      if(FirstSample4){
		Selected_NTracks_ATrackIsSample4_Data->Fill(nTracks[irec],weight);
		FirstSample4=false;
	      }

	      IronDistanceTotal_Data->Fill(IronDistance[irec][itrk],weight);
	      NTracksTotalSample4_Data->Fill(nTracks[irec],weight);
	      if(nTracks[irec]==1) IronDistanceTotalSample4_1Track_Data->Fill(IronDistance[irec][itrk],weight);
	      if(nTracks[irec]==2){
		if(MuonRec==itrk) IronDistanceTotalSample4_2Tracks_Data->Fill(IronDistance[irec][MuonRec],weight);
		for(int itrk2=0;itrk2<nTracks[irec];itrk2++){
		  if(itrk==itrk2) continue;
		  SampleTotalSample4_Data->Fill(Sample[irec][itrk2],weight);
		}
	      }
	    }

	    if(Sample[irec][itrk]==4){
	      if(IronDistance[irec][itrk]>10){
		int BinIron=(int) (IronDistance[irec][itrk]/5)-2;
		if(IronDistance[irec][itrk]>90) BinIron=NBinsIron-1;
		hCriteriaAngleXY_Data[BinIron]->Fill(CriteriaAngleX[irec][itrk],CriteriaAngleY[irec][itrk],weight);
		hCriteriaHalfWayXY_Data[BinIron]->Fill(CriteriaHalfWayX[irec][itrk],CriteriaHalfWayY[irec][itrk],weight);
	      }
	    }
	    hCriteriaHalfWayX_Data[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaHalfWayX[irec][itrk],weight);
	    hCriteriaHalfWayY_Data[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaHalfWayY[irec][itrk],weight);
	    hCriteriaAngleX_Data[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaAngleX[irec][itrk],weight);
	    hCriteriaAngleY_Data[Sample[irec][itrk]]->Fill(IronDistance[irec][itrk],CriteriaAngleY[irec][itrk],weight);


	  }
	}
	//if(Sample[irec][MuonRec]==4) IronDistanceTotal_Data->Fill(IronDistance[irec][MuonRec],weight);
	//if(Sample[irec][MuonRec]==5) AngleTotal_Data->Fill(TrackAngle[irec][MuonRec],weight);


	Reconstructed_NTracks_Data[Sample[irec][MuonRec]]->Fill(nTracks[irec],weight);
	Reconstructed_IronDistance_Data[Sample[irec][MuonRec]]->Fill(IronDistance[irec][MuonRec],weight);
	Reconstructed_ReconstructedAngle_Data[Sample[irec][MuonRec]]->Fill(TrackAngle[irec][MuonRec],weight);
  
	if((nTracks[irec]==1 && CLMuon[irec][0]>0.7) || ((nTracks[irec]==2 && TMath::Max(CLMuon[irec][0],CLMuon[irec][1])>0.7 && TMath::Min(CLMuon[irec][0],CLMuon[irec][1])<0.3) && (TMath::Min(CLMuon[irec][0],CLMuon[irec][1])>=0)) || ((nTracks[irec]==3 && CLMuon[irec][MuonRec]>0.7 && CLMuon[irec][(MuonRec+1)%3]<0.3 && CLMuon[irec][(MuonRec+2)%3]<0.3)) || ((nTracks[irec]==4 && CLMuon[irec][MuonRec]>0.7 && CLMuon[irec][(MuonRec+1)%4]<0.3 && CLMuon[irec][(MuonRec+2)%4]<0.3 && CLMuon[irec][(MuonRec+3)%4]<0.3))){
	  //if((nTracks[irec]==1 && CLMuon[irec][0]>0.7) || (nTracks[irec]==2 && TMath::Max(CLMuon[irec][0],CLMuon[irec][1])>0.7 && TMath::Min(CLMuon[irec][0],CLMuon[irec][1])<0.3) && (TMath::Min(CLMuon[irec][0],CLMuon[irec][1])>=0) || (nTracks[irec]==3 && CLMuon[irec][MuonRec]>0.7 && CLMuon[irec][(MuonRec+1)%3]<0.3 && CLMuon[irec][(MuonRec+2)%3]<0.3) || (nTracks[irec]==4 && CLMuon[irec][MuonRec]>0.7 && CLMuon[irec][(MuonRec+1)%4]<0.3 && CLMuon[irec][(MuonRec+2)%4]<0.3 && CLMuon[irec][(MuonRec+3)%4]<0.3)){
	  //cout<<"ok"<<endl;
	  //cout<<FSIInt<<", "<<Sample[irec][MuonRec]<<endl;
	  //if(Sample[irec][MuonRec]==4) cout<<"AngleX="<<CriteriaAngleX[irec][MuonRec]<<", AngleY="<<CriteriaAngleY[irec][MuonRec]<<", Haflway X="<<CriteriaHalfWayX[irec][MuonRec]<<", HalfwayY="<<CriteriaHalfWayY[irec][MuonRec]<<endl;
	  if(Sample[irec][MuonRec]>=4) TrackWidthXIronDistance_Data->Fill(IronDistance[irec][MuonRec],TrackWidth[irec][MuonRec],weight);

	  if(Sample[irec][MuonRec]==4){
	    if(IronDistance[irec][MuonRec]>10){
	      int BinIron=(int) (IronDistance[irec][MuonRec]/5)-2;
	      if(IronDistance[irec][MuonRec]>90) BinIron=NBinsIron-1;
	      hCriteriaAngleXY_Data[BinIron]->Fill(CriteriaAngleX[irec][MuonRec],CriteriaAngleY[irec][MuonRec],weight);
	      hCriteriaHalfWayXY_Data[BinIron]->Fill(CriteriaHalfWayX[irec][MuonRec],CriteriaHalfWayY[irec][MuonRec],weight);
	    }
	  }
	  //if(Sample[irec][MuonRec]==4 && PlasticDistance[irec][MuonRec]<25) continue;
	  if(Sample[irec][MuonRec]==4 && IronDistance[irec][MuonRec]<17){
	    Problematic_PlasticDistance_Data->Fill(PlasticDistance[irec][MuonRec],weight);
	    Problematic_Angle_Data->Fill(TrackAngle[irec][MuonRec],weight);
	    Problematic_NTracks_Data->Fill(nTracks[irec],weight);
	    Problematic_CLMuon_Data->Fill(CLMuon[irec][MuonRec],weight);
	    Problematic_TrackWidth_Data->Fill(TrackWidth[irec][MuonRec],weight);
	    Problematic_Geom_Data->Fill(GeometricTrack[irec][MuonRec],weight);
	    int LS;
	    if(TMath::Abs(LastChannelINGRIDX[irec][MuonRec]-11.5)>TMath::Abs(LastChannelINGRIDY[irec][MuonRec]-11.5)) LS=LastChannelINGRIDX[irec][MuonRec];
	    else LS=LastChannelINGRIDY[irec][MuonRec];
	    Problematic_LastChannel_Data->Fill(LS,weight);

	  }
	      
	  //if(Sample[irec][MuonRec]==4 /*&& IronDistance[irec][MuonRec]<25*/ && (CriteriaAngleX[irec][MuonRec]>25 || CriteriaAngleY[irec][MuonRec]>25 || CriteriaHalfWayX[irec][MuonRec]>7 || CriteriaHalfWayY[irec][MuonRec]>7 || TrackWidth[irec][MuonRec]>1.1)) continue; 
	  //if(Sample[irec][MuonRec]==4 && TrackWidth[irec][MuonRec]>1.1 && IronDistance[irec][MuonRec]<25) continue;	      
	  //if(Sample[irec][MuonRec]==4 && TrackWidth[irec][MuonRec]>1.5 && IronDistance[irec][MuonRec]>25) continue;	    
	  PID_Selected_IronDistance_Data[Sample[irec][MuonRec]]->Fill(IronDistance[irec][MuonRec],weight);
	  PID_Selected_ReconstructedAngle_Data[Sample[irec][MuonRec]]->Fill(TrackAngle[irec][MuonRec],weight);
	  if(Sample[irec][MuonRec]==4){
	    INGRIDAngleCriteria_Data->Fill(max(CriteriaAngleX[irec][MuonRec],CriteriaAngleY[irec][MuonRec]),weight);
	    INGRIDPositionCriteria_Data->Fill(max(CriteriaHalfWayX[irec][MuonRec],CriteriaHalfWayY[irec][MuonRec]),weight);
	  }

	  if(Sample[irec][MuonRec]>=4 && (CriteriaAngleX[irec][MuonRec]>25 || CriteriaAngleY[irec][MuonRec]>25 || CriteriaHalfWayX[irec][MuonRec]>7 || CriteriaHalfWayY[irec][MuonRec]>7)) continue;
	  PIDAndMatching_Selected_IronDistance_Data[Sample[irec][MuonRec]]->Fill(IronDistance[irec][MuonRec],weight);
	  PIDAndMatching_Selected_ReconstructedAngle_Data[Sample[irec][MuonRec]]->Fill(TrackAngle[irec][MuonRec],weight);


	  if(Sample[irec][MuonRec]>=4){
	    if(Sample[irec][MuonRec]==4) INGRIDTrackWidth_Data->Fill(TrackWidth[irec][MuonRec],weight);
	    if(TrackWidth[irec][MuonRec]>1.1 && IronDistance[irec][MuonRec]<15) continue;
	    else if(TrackWidth[irec][MuonRec]>1.2 && IronDistance[irec][MuonRec]<20) continue;
	    else if(TrackWidth[irec][MuonRec]>1.3 && IronDistance[irec][MuonRec]<30) continue;
	    else if(TrackWidth[irec][MuonRec]>1.5) continue;
	  }
	  if(Sample[irec][MuonRec]==4) Selected2D_Data->Fill(IronDistance[irec][MuonRec],TrackAngle[irec][MuonRec],weight);
	  Selected_NTracks_Data[Sample[irec][MuonRec]]->Fill(nTracks[irec],weight);
	  Selected_IronDistance_Data[Sample[irec][MuonRec]]->Fill(IronDistance[irec][MuonRec],weight);
	  Selected_ReconstructedAngle_Data[Sample[irec][MuonRec]]->Fill(TrackAngle[irec][MuonRec],weight);
	}
	Interacting_NTracks_Data[Sample[irec][MuonRec]]->Fill(nTracks[irec],weight);
	Interacting_IronDistance_Data[Sample[irec][MuonRec]]->Fill(IronDistance[irec][MuonRec],weight);
	Interacting_ReconstructedAngle_Data[Sample[irec][MuonRec]]->Fill(TrackAngle[irec][MuonRec],weight);
      }
    }
    //}
  }

  //cout<<endl<<endl<<endl<<POTCount<<endl;
  double ScalingFactor=POTCount/(9.91e21);
  //double ScalingFactor=1/(2.7);//(9.91e21);
  //double ScalingFactor=1;//POTCount/weightMC;
 
  for(int j=1;j<11;j++){
    InFV_Interacting[j]->Scale(ScalingFactor);
    InFV_Reconstructed[j]->Scale(ScalingFactor);

    Reconstructed_NTracks_ATrackIsSample4[j]->Scale(ScalingFactor);
    Reconstructed_NTracks_ATrackIsSample5[j]->Scale(ScalingFactor);

    Selected_NTracks_ATrackIsSample4[j]->Scale(ScalingFactor);
    Selected_NTracks_ATrackIsSample5[j]->Scale(ScalingFactor);

    INGRIDTrackWidth[j]->Scale(ScalingFactor);
    INGRIDAngleCriteria[j]->Scale(ScalingFactor);
    INGRIDPositionCriteria[j]->Scale(ScalingFactor);

    for(int i=0;i<6;i++){
      Selected_NTracks[i][j]->Scale(ScalingFactor);
      Reconstructed_NTracks[i][j]->Scale(ScalingFactor);
      Interacting_NTracks[i][j]->Scale(ScalingFactor);
      Selected_Momentum[i][j]->Scale(ScalingFactor);
      Reconstructed_Momentum[i][j]->Scale(ScalingFactor);
      Interacting_Momentum[i][j]->Scale(ScalingFactor);


      PID_Selected_IronDistance[i][j]->Scale(ScalingFactor);
      PIDAndMatching_Selected_IronDistance[i][j]->Scale(ScalingFactor);
      PID_Selected_ReconstructedAngle[i][j]->Scale(ScalingFactor);
      PIDAndMatching_Selected_ReconstructedAngle[i][j]->Scale(ScalingFactor);

      Selected_IronDistance[i][j]->Scale(ScalingFactor);
      Reconstructed_IronDistance[i][j]->Scale(ScalingFactor);
      Interacting_IronDistance[i][j]->Scale(ScalingFactor);
      Selected_Angle[i][j]->Scale(ScalingFactor);
      Reconstructed_Angle[i][j]->Scale(ScalingFactor);
      Interacting_Angle[i][j]->Scale(ScalingFactor);
      Selected_ReconstructedAngle[i][j]->Scale(ScalingFactor);
      Reconstructed_ReconstructedAngle[i][j]->Scale(ScalingFactor);
      Interacting_ReconstructedAngle[i][j]->Scale(ScalingFactor);
    }
  }

  for(int i=0;i<6;i++){
    Selected_IronDistance_MC[i]->Scale(ScalingFactor);
  }
  THStack * Stack_Reconstructed_NTracks = new THStack("Stack_Reconstructed_NTracks","");

  THStack * Stack_InFV_Interacting = new THStack("Stack_InFV_Interacting","");
  THStack * Stack_InFV_Reconstructed = new THStack("Stack_InFV_Reconstructed","");
  THStack * Stack_Selected_NTracks_Sample4 = new THStack("Stack_Selected_NTracks_Sample4","");
  THStack * Stack_Reconstructed_NTracks_Sample4 = new THStack("Stack_Reconstructed_NTracks_Sample4","");
  THStack * Stack_Selected_Momentum_Sample4 = new THStack("Stack_Selected_Momentum_Sample4","");
  THStack * Stack_Reconstructed_Momentum_Sample4 = new THStack("Stack_Reconstructed_Momentum_Sample4","");
  THStack * Stack_Selected_IronDistance_Sample4 = new THStack("Stack_Selected_IronDistance_Sample4","");
  THStack * Stack_PID_Selected_IronDistance_Sample4 = new THStack("Stack_PID_Selected_IronDistance_Sample4","");
  THStack * Stack_PIDAndMatching_Selected_IronDistance_Sample4 = new THStack("Stack_PIDAndMatching_Selected_IronDistance_Sample4","");
  THStack * Stack_PID_Selected_ReconstructedAngle_Sample4 = new THStack("Stack_PID_Selected_ReconstructedAngle_Sample4","");
  THStack * Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4 = new THStack("Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4","");

  THStack * Stack_Reconstructed_Angle_Sample4 = new THStack("Stack_Reconstructed_Angle_Sample4","");
  THStack * Stack_Selected_Angle_Sample4 = new THStack("Stack_Selected_Angle_Sample4","");
  THStack * Stack_Reconstructed_ReconstructedAngle_Sample4 = new THStack("Stack_Reconstructed_ReconstructedAngle_Sample4","");
  THStack * Stack_Selected_ReconstructedAngle_Sample4 = new THStack("Stack_Selected_ReconstructedAngle_Sample4","");

  THStack * Stack_Reconstructed_IronDistance_Sample4 = new THStack("Stack_Reconstructed_IronDistance_Sample4","");

  leg_Sample4= new TLegend(0.70,0.4,.99,.92);
  leg_Sample4->SetFillColor(0);
  for(int fsi=1;fsi<11;fsi++){
    if(fsi==1) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"CC-0Pi");
    //else if(fsi==2) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"CC-0Pi (Multi-protons)");
    else if(fsi==3) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"CC-1Pi");
    else if(fsi==4) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"CC-NPi");
    else if(fsi==5) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"CC-Other");
    else if(fsi==6) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"NC");
    else if(fsi==7) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"Wall Bkg");
    else if(fsi==8) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"INGRID H Bkg");
    else if(fsi==9) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"INGRID V Bkg");
    else if(fsi==10) leg_Sample4->AddEntry(Reconstructed_NTracks[4][fsi],"#nu_{#mu}+#nu_{e}");
  }

  leg_Sample4_Data= new TLegend(0.70,0.4,.99,.92);
  leg_Sample4_Data->SetFillColor(0);
  leg_Sample4_Data->AddEntry(Selected_IronDistance_Data[4],"Data");
  
  for(int fsi=1;fsi<11;fsi++){
    if(fsi==1) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"CC-0Pi");
    //else if(fsi==2) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"CC-0Pi (Multi-protons)");
    else if(fsi==3) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"CC-1Pi");
    else if(fsi==4) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"CC-NPi");
    else if(fsi==5) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"CC-Other");
    else if(fsi==6) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"NC");
    else if(fsi==7) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"Wall Bkg");
    else if(fsi==8) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"INGRID H Bkg");
    else if(fsi==9) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"INGRID V Bkg");
    else if(fsi==10) leg_Sample4_Data->AddEntry(Selected_IronDistance[4][fsi],"#nu_{#mu}+#nu_{e}");
  }

  THStack * Stack_Selected_NTracks_Sample5 = new THStack("Stack_Selected_NTracks_Sample5","");
  THStack * Stack_Reconstructed_NTracks_Sample5 = new THStack("Stack_Reconstructed_NTracks_Sample5","");
  THStack * Stack_Selected_Momentum_Sample5 = new THStack("Stack_Selected_Momentum_Sample5","");
  THStack * Stack_Reconstructed_Momentum_Sample5 = new THStack("Stack_Reconstructed_Momentum_Sample5","");
  THStack * Stack_Selected_Angle_Sample5 = new THStack("Stack_Selected_Angle_Sample5","");
  THStack * Stack_Reconstructed_Angle_Sample5 = new THStack("Stack_Reconstructed_Angle_Sample5","");
  THStack * Stack_Selected_ReconstructedAngle_Sample5 = new THStack("Stack_Selected_ReconstructedAngle_Sample5","");
  THStack * Stack_Reconstructed_ReconstructedAngle_Sample5 = new THStack("Stack_Reconstructed_ReconstructedAngle_Sample5","");
  THStack * Stack_INGRIDTrackWidth = new THStack("Stack_INGRIDTrackWidth","");
  THStack * Stack_INGRIDAngleCriteria = new THStack("Stack_INGRIDAngleCriteria","");
  THStack * Stack_INGRIDPositionCriteria = new THStack("Stack_INGRIDPositionCriteria","");

  TLegend * leg_Sample5= new TLegend(0.70,0.55,.99,.92);
  leg_Sample5->SetFillColor(0);

  TH1D * Selected_NTracks_Sample4_All;
  TH1D * Reconstructed_NTracks_Sample4_All;
  TH1D * Interacting_NTracks_Sample4_All;

  TH1D * Selected_NTracks_Sample5_All;
  TH1D * Reconstructed_NTracks_Sample5_All;
  TH1D * Interacting_NTracks_Sample5_All;

  TH1D * Selected_Momentum_Sample4_All;
  TH1D * Reconstructed_Momentum_Sample4_All;
  TH1D * Interacting_Momentum_Sample4_All;

  TH1D * Selected_IronDistance_Sample4_All;
  TH1D * Reconstructed_IronDistance_Sample4_All;
  TH1D * Interacting_IronDistance_Sample4_All;

  TH1D * Selected_Momentum_Sample5_All;
  TH1D * Reconstructed_Momentum_Sample5_All;
  TH1D * Interacting_Momentum_Sample5_All;

  TH1D * Selected_Angle_Sample5_All;
  TH1D * Reconstructed_Angle_Sample5_All;
  TH1D * Interacting_Angle_Sample5_All;
 
  TH1D * Selected_ReconstructedAngle_Sample5_All;
  TH1D * Reconstructed_ReconstructedAngle_Sample5_All;
  TH1D * Interacting_ReconstructedAngle_Sample5_All;

 
  for(int fsi=1;fsi<11;fsi++){
    if(fsi==1){
      Selected_Momentum_Sample4_All = (TH1D*) Selected_Momentum[4][fsi]->Clone("Selected_Momentum_Sample4_All");
      Reconstructed_Momentum_Sample4_All = (TH1D*) Reconstructed_Momentum[4][fsi]->Clone("Reconstructed_Momentum_Sample4_All");
      Interacting_Momentum_Sample4_All = (TH1D*) Interacting_Momentum[4][fsi]->Clone("Interacting_Momentum_Sample4_All");

      Selected_IronDistance_Sample4_All = (TH1D*) Selected_IronDistance[4][fsi]->Clone("Selected_IronDistance_Sample4_All");
      Reconstructed_IronDistance_Sample4_All = (TH1D*) Reconstructed_IronDistance[4][fsi]->Clone("Reconstructed_IronDistance_Sample4_All");
      Interacting_IronDistance_Sample4_All = (TH1D*) Interacting_IronDistance[4][fsi]->Clone("Interacting_IronDistance_Sample4_All");


      Selected_NTracks_Sample4_All = (TH1D*) Selected_NTracks[4][fsi]->Clone("Selected_NTracks_Sample4_All");
      Reconstructed_NTracks_Sample4_All = (TH1D*) Reconstructed_NTracks[4][fsi]->Clone("Reconstructed_NTracks_Sample4_All");
      Interacting_NTracks_Sample4_All = (TH1D*) Interacting_NTracks[4][fsi]->Clone("Interacting_NTracks_Sample4_All");


      Selected_Momentum_Sample5_All = (TH1D*) Selected_Momentum[5][fsi]->Clone("Selected_Momentum_Sample5_All");
      Reconstructed_Momentum_Sample5_All = (TH1D*) Reconstructed_Momentum[5][fsi]->Clone("Reconstructed_Momentum_Sample5_All");
      Interacting_Momentum_Sample5_All = (TH1D*) Interacting_Momentum[5][fsi]->Clone("Interacting_Momentum_Sample5_All");

      Selected_Angle_Sample5_All = (TH1D*) Selected_Angle[5][fsi]->Clone("Selected_Angle_Sample5_All");
      Reconstructed_Angle_Sample5_All = (TH1D*) Reconstructed_Angle[5][fsi]->Clone("Reconstructed_Angle_Sample5_All");
      Interacting_Angle_Sample5_All = (TH1D*) Interacting_Angle[5][fsi]->Clone("Interacting_Angle_Sample5_All");

      Selected_ReconstructedAngle_Sample5_All = (TH1D*) Selected_ReconstructedAngle[5][fsi]->Clone("Selected_ReconstructedAngle_Sample5_All");
      Reconstructed_ReconstructedAngle_Sample5_All = (TH1D*) Reconstructed_ReconstructedAngle[5][fsi]->Clone("Reconstructed_ReconstructedAngle_Sample5_All");
      Interacting_ReconstructedAngle_Sample5_All = (TH1D*) Interacting_ReconstructedAngle[5][fsi]->Clone("Interacting_ReconstructedAngle_Sample5_All");

      Selected_NTracks_Sample5_All = (TH1D*) Selected_NTracks[5][fsi]->Clone("Selected_NTracks_Sample5_All");
      Reconstructed_NTracks_Sample5_All = (TH1D*) Reconstructed_NTracks[5][fsi]->Clone("Reconstructed_NTracks_Sample5_All");
      Interacting_NTracks_Sample5_All = (TH1D*) Interacting_NTracks[5][fsi]->Clone("Interacting_NTracks_Sample5_All");

      //InFV_Interacting[fsi] = (TH1D*) InFV_Interacting[fsi]->Clone("InFV_Interacting");
    }
    else if(fsi!=2){
      Selected_NTracks_Sample4_All->Add(Selected_NTracks[4][fsi]);
      Reconstructed_NTracks_Sample4_All->Add(Reconstructed_NTracks[4][fsi]);
      Interacting_NTracks_Sample4_All->Add(Interacting_NTracks[4][fsi]);

      Selected_Momentum_Sample4_All->Add(Selected_Momentum[4][fsi]);
      Reconstructed_Momentum_Sample4_All->Add(Reconstructed_Momentum[4][fsi]);
      Interacting_Momentum_Sample4_All->Add(Interacting_Momentum[4][fsi]);

      Selected_NTracks_Sample5_All->Add(Selected_NTracks[5][fsi]);
      Reconstructed_NTracks_Sample5_All->Add(Reconstructed_NTracks[5][fsi]);
      Interacting_NTracks_Sample5_All->Add(Interacting_NTracks[5][fsi]);

      Selected_Momentum_Sample5_All->Add(Selected_Momentum[5][fsi]);
      Reconstructed_Momentum_Sample5_All->Add(Reconstructed_Momentum[5][fsi]);
      Interacting_Momentum_Sample5_All->Add(Interacting_Momentum[5][fsi]);
    }
  }



  ProduceStack(Selected_NTracks[4],Stack_Selected_NTracks_Sample4);
  ProduceStack(Reconstructed_NTracks[4],Stack_Reconstructed_NTracks_Sample4);
  for(int j=0;j<6;j++){
    bool Last=false;
    if(j==5) Last=true;
    ProduceStack2(Reconstructed_NTracks[j],Stack_Reconstructed_NTracks,Last);
  }
  ProduceStack(Selected_Momentum[4],Stack_Selected_Momentum_Sample4);
  ProduceStack(Reconstructed_Momentum[4],Stack_Reconstructed_Momentum_Sample4);
  ProduceStack(Selected_IronDistance[4],Stack_Selected_IronDistance_Sample4);
  ProduceStack(PID_Selected_IronDistance[4],Stack_PID_Selected_IronDistance_Sample4);
  ProduceStack(PIDAndMatching_Selected_IronDistance[4],Stack_PIDAndMatching_Selected_IronDistance_Sample4);
  ProduceStack(PID_Selected_ReconstructedAngle[4],Stack_PID_Selected_ReconstructedAngle_Sample4);
  ProduceStack(PIDAndMatching_Selected_ReconstructedAngle[4],Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4);
 
  ProduceStack(Reconstructed_IronDistance[4],Stack_Reconstructed_IronDistance_Sample4);
  ProduceStack(Selected_Angle[4],Stack_Selected_Angle_Sample4);
  ProduceStack(Reconstructed_Angle[4],Stack_Reconstructed_Angle_Sample4);
  ProduceStack(Selected_ReconstructedAngle[4],Stack_Selected_ReconstructedAngle_Sample4);
  ProduceStack(Reconstructed_ReconstructedAngle[4],Stack_Reconstructed_ReconstructedAngle_Sample4);


  ProduceStack(Selected_NTracks[5],Stack_Selected_NTracks_Sample5);
  ProduceStack(Reconstructed_NTracks[5],Stack_Reconstructed_NTracks_Sample5);
  ProduceStack(Selected_Momentum[5],Stack_Selected_Momentum_Sample5);
  ProduceStack(Reconstructed_Momentum[5],Stack_Reconstructed_Momentum_Sample5);
  ProduceStack(Selected_Angle[5],Stack_Selected_Angle_Sample5);
  ProduceStack(Reconstructed_Angle[5],Stack_Reconstructed_Angle_Sample5);
  ProduceStack(Selected_ReconstructedAngle[5],Stack_Selected_ReconstructedAngle_Sample5);
  ProduceStack(Reconstructed_ReconstructedAngle[5],Stack_Reconstructed_ReconstructedAngle_Sample5);

  ProduceStack(InFV_Interacting,Stack_InFV_Interacting);
  ProduceStack(InFV_Reconstructed,Stack_InFV_Reconstructed);

  ProduceStack(INGRIDTrackWidth,Stack_INGRIDTrackWidth);
  ProduceStack(INGRIDPositionCriteria,Stack_INGRIDPositionCriteria);
  ProduceStack(INGRIDAngleCriteria,Stack_INGRIDAngleCriteria);


  TCanvas * c0 = new TCanvas("c0","Reconstructed by TN-160 selection: All samples");
  Stack_Reconstructed_NTracks->SetTitle("Interaction type with reconstructed track");
  Stack_Reconstructed_NTracks->Draw();
  TH1D * hStack_Reconstructed_NTracks = (TH1D*) Stack_Reconstructed_NTracks->GetHistogram();
  hStack_Reconstructed_NTracks->GetXaxis()->SetTitle("Number of reconstructed tracks");
  hStack_Reconstructed_NTracks->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c0->SaveAs("NTracksReconstructed.pdf");
  for(int fsi=1;fsi<11;fsi++){
    int j=4;
    if(fsi==1) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==3) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==4) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==5) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==6) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==7) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==8) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==9) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
    else if(fsi==10) Reconstructed_NTracks[j][fsi]->SetLineColor(1);
  }

  TCanvas * c13 = new TCanvas();
  Stack_InFV_Interacting->SetTitle("Interaction type with E_{#nu}");
  Stack_InFV_Interacting->Draw();
  TH1D * hStack_InFV_Interacting = (TH1D*) Stack_InFV_Interacting->GetHistogram();
  hStack_InFV_Interacting->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hStack_InFV_Interacting->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c13->SaveAs("InFV_Interacting.pdf");

  TCanvas * c14 = new TCanvas();
  Stack_InFV_Reconstructed->SetTitle("Interaction type with E_{#nu}");
  Stack_InFV_Reconstructed->Draw();
  TH1D * hStack_InFV_Reconstructed = (TH1D*) Stack_InFV_Reconstructed->GetHistogram();
  hStack_InFV_Reconstructed->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hStack_InFV_Reconstructed->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c14->SaveAs("InFV_Reconstructed.pdf");

  TCanvas * c15 = new TCanvas();
  TH1D * InFV_Interacting_Clone = new TH1D("InFV_Interacting_Clone","",100,0,10);
  TH1D * InFV_Reconstructed_Clone = new TH1D("InFV_Reconstructed_Clone","",100,0,10);
  for(int fsi=1;fsi<11;fsi++){
    if(fsi==1){
      InFV_Interacting_Clone->Add(InFV_Interacting[fsi],1.);
      InFV_Reconstructed_Clone->Add(InFV_Reconstructed[fsi],1.);
    }
  }
  InFV_Interacting_Clone->Sumw2();
  InFV_Reconstructed_Clone->Sumw2();
  InFV_Reconstructed_Clone->Divide(InFV_Interacting_Clone);
  /*
    for(int ibinx=1;ibinx<=InFV_Interacting_Clone->GetNbinsX();ibinx++){
    double Value=InFV_Reconstructed_Clone->GetBinContent(ibinx);
    if(InFV_Interacting_Clone->GetBinContent(ibinx)!=0) Value/=InFV_Interacting_Clone->GetBinContent(ibinx);
    cout<<"Bin="<<ibinx<<", Value="<<Value<<endl;
    InFV_Reconstructed_Clone->SetBinContent(ibinx,Value);
    }
  */
  
  for(int ibinx=1;ibinx<=InFV_Interacting_Clone->GetNbinsX();ibinx++){
    //cout<<"Bin value="<<InFV_Reconstructed_Clone->GetBinContent(ibinx)<<", second="<<InFV_Interacting_Clone->GetBinContent(ibinx)<<endl;
    double Value;
    Value=(InFV_Reconstructed_Clone->GetBinContent(ibinx)*(1-InFV_Reconstructed_Clone->GetBinContent(ibinx)));
    //cout<<"Value first="<<Value<<", ";
    Value=TMath::Sqrt(Value);
    //cout<<"Value second="<<Value<<", ";
    if(InFV_Interacting_Clone->GetBinContent(ibinx)!=0) Value/=TMath::Sqrt(InFV_Interacting_Clone->GetBinContent(ibinx));
    //cout<<"Value third="<<Value<<endl;
    InFV_Reconstructed_Clone->SetBinError(ibinx,Value);
    //else InFV_Reconstructed_Clone->SetBinError(ibinx,1);//Temporary
  }

  InFV_Reconstructed_Clone->Draw();
  //leg_Sample4->Draw("same");
  c15->SaveAs("Efficiency.pdf");

  TCanvas * c1 = new TCanvas();
  Stack_Selected_NTracks_Sample4->SetTitle("INGRID stopped #mu: Interaction type with selected tracks");

  //(Stack_Selected_NTracks_Sample4->GetXaxis())->SetTitle("Number of reconstructed tracks");
  //(Stack_Selected_NTracks_Sample4->GetYaxis())->SetTitle("Number of events");
  Stack_Selected_NTracks_Sample4->Draw();
  TH1D * hStack_Selected_NTracks_Sample4 = (TH1D*) Stack_Selected_NTracks_Sample4->GetHistogram();
  hStack_Selected_NTracks_Sample4->GetXaxis()->SetTitle("Number of reconstructed tracks");
  hStack_Selected_NTracks_Sample4->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c1->SaveAs("Stopping_NTracksSelected.pdf");

  TCanvas * c2 = new TCanvas();
  Stack_Reconstructed_NTracks_Sample4->SetTitle("INGRID stopped #mu: Interaction type with reconstructed tracks");
  Stack_Reconstructed_NTracks_Sample4->Draw();
  TH1D * hStack_Reconstructed_NTracks_Sample4 = (TH1D*) Stack_Reconstructed_NTracks_Sample4->GetHistogram();
  hStack_Reconstructed_NTracks_Sample4->GetXaxis()->SetTitle("Number of reconstructed tracks");
  hStack_Reconstructed_NTracks_Sample4->GetYaxis()->SetTitle("Number of events");

  double TotalSample4_Reconstructed=Reconstructed_NTracks[4][1]->Integral()+Reconstructed_NTracks[4][2]->Integral()+Reconstructed_NTracks[4][3]->Integral()+Reconstructed_NTracks[4][4]->Integral()+Reconstructed_NTracks[4][5]->Integral()+Reconstructed_NTracks[4][6]->Integral()+Reconstructed_NTracks[4][7]->Integral()+Reconstructed_NTracks[4][8]->Integral()+Reconstructed_NTracks[4][9]->Integral()+Reconstructed_NTracks[4][10]->Integral();
  double TotalSample5_Reconstructed=Reconstructed_NTracks[5][1]->Integral()+Reconstructed_NTracks[5][2]->Integral()+Reconstructed_NTracks[5][3]->Integral()+Reconstructed_NTracks[5][4]->Integral()+Reconstructed_NTracks[5][5]->Integral()+Reconstructed_NTracks[5][6]->Integral()+Reconstructed_NTracks[5][7]->Integral()+Reconstructed_NTracks[5][8]->Integral()+Reconstructed_NTracks[5][9]->Integral()+Reconstructed_NTracks[5][10]->Integral();
  double TotalCC0pi_Reconstructed = Reconstructed_NTracks[0][1]->Integral()+Reconstructed_NTracks[1][1]->Integral()+Reconstructed_NTracks[2][1]->Integral()+Reconstructed_NTracks[3][1]->Integral()+Reconstructed_NTracks[4][1]->Integral()+Reconstructed_NTracks[5][1]->Integral();
  double TotalCC1pi_Reconstructed = Reconstructed_NTracks[0][3]->Integral()+Reconstructed_NTracks[1][3]->Integral()+Reconstructed_NTracks[2][3]->Integral()+Reconstructed_NTracks[3][3]->Integral()+Reconstructed_NTracks[4][3]->Integral()+Reconstructed_NTracks[5][3]->Integral();
  double TotalNC_Reconstructed = Reconstructed_NTracks[0][6]->Integral()+Reconstructed_NTracks[1][6]->Integral()+Reconstructed_NTracks[2][6]->Integral()+Reconstructed_NTracks[3][6]->Integral()+Reconstructed_NTracks[4][6]->Integral()+Reconstructed_NTracks[5][6]->Integral();



  double TotalSample4_Selected=Selected_NTracks[4][1]->Integral()+Selected_NTracks[4][2]->Integral()+Selected_NTracks[4][3]->Integral()+Selected_NTracks[4][4]->Integral()+Selected_NTracks[4][5]->Integral()+Selected_NTracks[4][6]->Integral()+Selected_NTracks[4][7]->Integral()+Selected_NTracks[4][8]->Integral()+Selected_NTracks[4][9]->Integral()+Selected_NTracks[4][10]->Integral();
  double TotalSample5_Selected=Selected_NTracks[5][1]->Integral()+Selected_NTracks[5][2]->Integral()+Selected_NTracks[5][3]->Integral()+Selected_NTracks[5][4]->Integral()+Selected_NTracks[5][5]->Integral()+Selected_NTracks[5][6]->Integral()+Selected_NTracks[5][7]->Integral()+Selected_NTracks[5][8]->Integral()+Selected_NTracks[5][9]->Integral()+Selected_NTracks[5][10]->Integral();
  double TotalCC0pi_Selected = Selected_NTracks[0][1]->Integral()+Selected_NTracks[1][1]->Integral()+Selected_NTracks[2][1]->Integral()+Selected_NTracks[3][1]->Integral()+Selected_NTracks[4][1]->Integral()+Selected_NTracks[5][1]->Integral();
  double TotalCC1pi_Selected = Selected_NTracks[0][3]->Integral()+Selected_NTracks[1][3]->Integral()+Selected_NTracks[2][3]->Integral()+Selected_NTracks[3][3]->Integral()+Selected_NTracks[4][3]->Integral()+Selected_NTracks[5][3]->Integral();
  double TotalNC_Selected = Selected_NTracks[0][6]->Integral()+Selected_NTracks[1][6]->Integral()+Selected_NTracks[2][6]->Integral()+Selected_NTracks[3][6]->Integral()+Selected_NTracks[4][6]->Integral()+Selected_NTracks[5][6]->Integral();

  double Total=0;
  double Total_Reconstructed;
  double Total_Selected;
  double TotalCC=0;
  double TotalCC1pi=0;
  double TotalNC=0;
  double TotalCC0pi=0;
  for(int fsi=1;fsi<11;fsi++){
    Total+=InFV_Interacting[fsi]->Integral();
    if(fsi<6) TotalCC+=InFV_Interacting[fsi]->Integral();
    if(fsi==6) TotalNC+=InFV_Interacting[fsi]->Integral();
    if(fsi==1) TotalCC0pi+=InFV_Interacting[fsi]->Integral();
    if(fsi==3) TotalCC1pi+=InFV_Interacting[fsi]->Integral();

    for(int j=0;j<6;j++) Total_Reconstructed+=Reconstructed_NTracks[j][fsi]->Integral();
    for(int j=0;j<6;j++) Total_Selected+=Reconstructed_NTracks[j][fsi]->Integral();
  }

  cout<<"Purity of the TN-160 reconstruction in CC-0pi after one track selection, for sample4:"<<Reconstructed_NTracks[4][1]->Integral()/TotalSample4_Reconstructed<<endl;

  cout<<"Purity of the TN-160 reconstruction in CC-0pi, for sample5:"<<Reconstructed_NTracks[5][1]->Integral()/TotalSample5_Reconstructed<<endl;
  cout<<"Purity of the TN-160 reconstruction in CC-0pi, All Samples:"<<TotalCC0pi_Reconstructed/Total_Reconstructed<<endl;

  cout<<"Purity of the TN-160 reconstruction in NC, for sample4:"<<Reconstructed_NTracks[4][6]->Integral()/TotalSample4_Reconstructed<<endl;
  cout<<"Purity of the TN-160 reconstruction in NC, for sample5:"<<Reconstructed_NTracks[5][6]->Integral()/TotalSample5_Reconstructed<<endl;
  cout<<"Purity of the TN-160 reconstruction in NC, All Samples:"<<TotalNC_Reconstructed/Total_Reconstructed<<endl;

  cout<<"FALSE: Efficiency of the TN-160 reconstruction in CC-0pi:"<< TotalCC0pi_Reconstructed/TotalCC0pi<<endl;
  cout<<"FALSE: Efficiency of the TN-160 reconstruction in CC-1pi:"<< TotalCC1pi_Reconstructed/TotalCC1pi<<endl;
  cout<<"FALSE: Efficiency of the TN-160 reconstruction in NC:"<< TotalNC_Reconstructed/TotalNC<<endl;
  cout<<"FALSE: Purity of the TN-160 reconstruction in CC-0pi:"<< InFV_Interacting[1]->Integral()/Total<<endl;
  cout<<"FALSE: Purity of the TN-160 reconstruction in CC:"<< TotalCC/Total<<endl;
  cout<<"FALSE: Purity of the TN-160 reconstruction in NC:"<< InFV_Interacting[6]->Integral()/Total<<endl;

  double TotalCC0pi_Sample4=0;
  double TotalCC0pi_Sample5=0;

  for(int fsi=1;fsi<11;fsi++){
 
    if(fsi==1) TotalCC0pi_Sample4+=Reconstructed_NTracks_ATrackIsSample4[fsi]->Integral();
    if(fsi==1) TotalCC0pi_Sample5+=Reconstructed_NTracks_ATrackIsSample5[fsi]->Integral();
  }

  double PriorPID_CC0pi=0;
  double PID_CC0pi=0;
  double PIDAndMatching_CC0pi=0;
  double PIDAndMatchingAndWidth_CC0pi=0;

  double PriorPID_CCOthers=0;
  double PID_CCOthers=0;
  double PIDAndMatching_CCOthers=0;
  double PIDAndMatchingAndWidth_CCOthers=0;

  double PriorPID_NC=0;
  double PID_NC=0;
  double PIDAndMatching_NC=0;
  double PIDAndMatchingAndWidth_NC=0;

  double PriorPID_Bkg=0;
  double PID_Bkg=0;
  double PIDAndMatching_Bkg=0;
  double PIDAndMatchingAndWidth_Bkg=0;

  double PriorPID_Total=0;
  double PID_Total=0;
  double PIDAndMatching_Total=0;
  double PIDAndMatchingAndWidth_Total=0;

  double PriorPID_CC0pi_Sample5=0;
  double PriorPID_Total_Sample5=0;
  double PID_CC0pi_Sample5=0;
  double PID_Total_Sample5=0;
  double PIDAndMatching_CC0pi_Sample5=0;
  double PIDAndMatching_Total_Sample5=0;
  double PIDAndMatchingAndWidth_CC0pi_Sample5=0;
  double PIDAndMatchingAndWidth_Total_Sample5=0;

  for(int fsi=1;fsi<11;fsi++){
    if(fsi==1){
      PriorPID_CC0pi+=Selected_NTracks_ATrackIsSample4[fsi]->Integral();
      PID_CC0pi+=PID_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatching_CC0pi+=PIDAndMatching_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatchingAndWidth_CC0pi+=Selected_IronDistance[4][fsi]->Integral();

      PriorPID_CC0pi_Sample5+=Selected_NTracks_ATrackIsSample5[fsi]->Integral();
      PID_CC0pi_Sample5+=PID_Selected_IronDistance[5][fsi]->Integral();
      PIDAndMatching_CC0pi_Sample5+=PIDAndMatching_Selected_IronDistance[5][fsi]->Integral();
      PIDAndMatchingAndWidth_CC0pi_Sample5+=Selected_IronDistance[5][fsi]->Integral();
    }
    else if(fsi<6){
      PriorPID_CCOthers+=Selected_NTracks_ATrackIsSample4[fsi]->Integral();
      PID_CCOthers+=PID_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatching_CCOthers+=PIDAndMatching_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatchingAndWidth_CCOthers+=Selected_IronDistance[4][fsi]->Integral();
    }
    else if(fsi==6){
      PriorPID_NC+=Selected_NTracks_ATrackIsSample4[fsi]->Integral();
      PID_NC+=PID_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatching_NC+=PIDAndMatching_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatchingAndWidth_NC+=Selected_IronDistance[4][fsi]->Integral();
    }
    else{
      PriorPID_Bkg+=Selected_NTracks_ATrackIsSample4[fsi]->Integral();
      PID_Bkg+=PID_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatching_Bkg+=PIDAndMatching_Selected_IronDistance[4][fsi]->Integral();
      PIDAndMatchingAndWidth_Bkg+=Selected_IronDistance[4][fsi]->Integral();
    }
    PriorPID_Total+=Selected_NTracks_ATrackIsSample4[fsi]->Integral();
    PID_Total+=PID_Selected_IronDistance[4][fsi]->Integral();
    PIDAndMatching_Total+=PIDAndMatching_Selected_IronDistance[4][fsi]->Integral();
    PIDAndMatchingAndWidth_Total+=Selected_IronDistance[4][fsi]->Integral();

    PriorPID_Total_Sample5+=Selected_NTracks_ATrackIsSample5[fsi]->Integral();
    PID_Total_Sample5+=PID_Selected_IronDistance[5][fsi]->Integral();
    PIDAndMatching_Total_Sample5+=PIDAndMatching_Selected_IronDistance[5][fsi]->Integral();
    PIDAndMatchingAndWidth_Total_Sample5+=Selected_IronDistance[5][fsi]->Integral();
  }

  cout<<"Selection:"<<endl;
  cout<<"Prior cut (only sample 4 selection): Number of events having at least one track sample 4 (CC0pi) = "<<PriorPID_CC0pi<<", Total="<<PriorPID_Total<<", Data="<<Selected_NTracks_ATrackIsSample4_Data->Integral()<<", Purity="<<PriorPID_CC0pi/PriorPID_Total<<endl;
  cout<<"Prior cut (only sample 4 selection): Number of events having a Mu like in sample 4 (CC0pi) = "<<PID_CC0pi<<", Total="<<PID_Total<<", Data="<<PID_Selected_IronDistance_Data[4]->Integral()<<", Purity="<<PID_CC0pi/PID_Total<<", Efficiency="<<PID_CC0pi/PriorPID_CC0pi<<endl;
  cout<<"Matching cut (only sample 4 selection): Number of events having a Mu like in sample 4 (CC0pi) = "<<PIDAndMatching_CC0pi<<", Total="<<PIDAndMatching_Total<<", Data="<<PIDAndMatching_Selected_IronDistance_Data[4]->Integral()<<", Purity="<<PIDAndMatching_CC0pi/PIDAndMatching_Total<<", Efficiency="<<PIDAndMatching_CC0pi/PID_CC0pi<<endl;
  cout<<"Width cut (only sample 4 selection): Number of events having a Mu like in sample 4 (CC0pi) = "<<PIDAndMatchingAndWidth_CC0pi<<", Total="<<PIDAndMatchingAndWidth_Total<<", Data="<<Selected_IronDistance_Data[4]->Integral()<<", Purity="<<PIDAndMatchingAndWidth_CC0pi/PIDAndMatchingAndWidth_Total<<", Efficiency="<<PIDAndMatchingAndWidth_CC0pi/PIDAndMatching_CC0pi<<endl;

  cout<<endl<<endl;
  cout<<" \\hline"<<endl;
  cout<<" & CC-0pi & CC-Others & NC & Other bkg & CC-0pi purity & CC-0pi efficiency & Total MC & Data"<<" \\"<<"\\"<<endl; 
  cout<<" \\hline"<<endl;
  cout<<"One track INGRID stopped & "<<PriorPID_CC0pi<<" & "<<PriorPID_CCOthers<<" & "<<PriorPID_NC<<" & "<<PriorPID_Bkg<<" & "<<PriorPID_CC0pi/PriorPID_Total<<" &  &"<<PriorPID_Total<<" & "<<Selected_NTracks_ATrackIsSample4_Data->Integral()<<" \\"<<"\\"<<endl;
  cout<<" \\hline"<<endl;
  cout<<"$\\mu_{CL}$ & "<<PID_CC0pi<<" & "<<PID_CCOthers<<" & "<<PID_NC<<" & "<<PID_Bkg<<" & "<<PID_CC0pi/PID_Total<<" & "<<PID_CC0pi/PriorPID_CC0pi<<" & "<<PID_Total<<" & "<<PID_Selected_IronDistance_Data[4]->Integral()<<" \\"<<"\\"<<endl;
  cout<<" \\hline"<<endl;
  cout<<"$\\mu$-like track matching & "<<PIDAndMatching_CC0pi<<" & "<<PIDAndMatching_CCOthers<<" & "<<PIDAndMatching_NC<<" & "<<PIDAndMatching_Bkg<<" & "<<PIDAndMatching_CC0pi/PIDAndMatching_Total<<" & "<<PIDAndMatching_CC0pi/PID_CC0pi<<" & "<<PIDAndMatching_Total<<" & "<<PIDAndMatching_Selected_IronDistance_Data[4]->Integral()<<" \\"<<"\\"<<endl;
  cout<<" \\hline"<<endl;
  cout<<"$\\mu$-like track with & "<<PIDAndMatchingAndWidth_CC0pi<<" & "<<PIDAndMatchingAndWidth_CCOthers<<" & "<<PIDAndMatchingAndWidth_NC<<" & "<<PIDAndMatchingAndWidth_Bkg<<" & "<<PIDAndMatchingAndWidth_CC0pi/PIDAndMatchingAndWidth_Total<<" & "<<PIDAndMatchingAndWidth_CC0pi/PIDAndMatching_CC0pi<<" & "<<PIDAndMatchingAndWidth_Total<<" & "<<Selected_IronDistance_Data[4]->Integral()<<" \\"<<"\\"<<endl;
  cout<<" \\hline"<<endl;




  cout<<"Purity of the Selected CC-0pi reconstruction in CC-0pi, for sample4:"<<Selected_NTracks[4][1]->Integral()/TotalSample4_Selected<<endl;
  cout<<"Purity of the Selected CC-0pi reconstruction in CC-0pi, for sample5:"<<Selected_NTracks[5][1]->Integral()/TotalSample5_Selected<<endl;
  cout<<"Purity of the Selected CC-0pi reconstruction in CC-0pi, All Samples:"<<TotalCC0pi_Selected/Total_Selected<<endl;

  cout<<"Purity of the Selected CC-0pi reconstruction in NC, for sample4:"<<Selected_NTracks[4][6]->Integral()/TotalSample4_Selected<<endl;
  cout<<"Purity of the Selected CC-0pi reconstruction in NC, for sample5:"<<Selected_NTracks[5][6]->Integral()/TotalSample5_Selected<<endl;
  cout<<"Purity of the Selected CC-0pi reconstruction in NC, All Samples:"<<TotalNC_Selected/Total_Selected<<endl;

  cout<<"FALSE: Efficiency of the Selected CC-0pi reconstruction in CC-0pi in sample 4:"<< Selected_NTracks[4][1]->Integral()/Reconstructed_NTracks_ATrackIsSample4[1]->Integral()<<endl;
  cout<<"FALSE: Efficiency of the Selected CC-0pi reconstruction in CC-0pi in sample 5:"<< Selected_NTracks[5][1]->Integral()/Reconstructed_NTracks_ATrackIsSample5[1]->Integral()<<endl;
  cout<<"FALSE: Efficiency of the Selected CC-0pi reconstruction in CC-0pi:"<< TotalCC0pi_Selected/TotalCC0pi<<endl;
  cout<<"FALSE: Efficiency of the Selected CC-0pi reconstruction in CC-1pi:"<< TotalCC1pi_Selected/TotalCC1pi<<endl;
  cout<<"FALSE: Efficiency of the Selected CC-0pi reconstruction in NC:"<< TotalNC_Selected/TotalNC<<endl;
 
  leg_Sample4->Draw("same");
  c2->SaveAs("Stopping_NTracksReconstructed.pdf");
  /*
    TCanvas * c20 = new TCanvas();
    Stack_Reconstructed_NTracks_Sample4->SetTitle("INGRID stopped #mu: Interaction type with reconstructed tracks");
    Stack_Reconstructed_NTracks_Sample4->Draw();
    TH1D * hStack_Reconstructed_NTracks_Sample4 = (TH1D*) Stack_Reconstructed_NTracks_Sample4->GetHistogram();
    hStack_Reconstructed_NTracks_Sample4->GetXaxis()->SetTitle("Number of reconstructed tracks");
    hStack_Reconstructed_NTracks_Sample4->GetYaxis()->SetTitle("Number of events");
    cout<<"Purity of the TN-160 reconstruction in CC-0pi:"<<Reconstructed_NTracks[4][1]->Integral()/Stack_Reconstructed_NTracks_Sample4->Integral()<<endl;
    leg_Sample4->Draw("same");
    c2->SaveAs("Stopping_NTracksReconstructed.pdf");
  */
  TCanvas * c3 = new TCanvas();
  Stack_Selected_Momentum_Sample4->SetTitle("CC-0pi like INGRID stopped #mu: Interaction type with p_{#mu}");
  Stack_Selected_Momentum_Sample4->Draw();
  TH1D * hStack_Selected_Momentum_Sample4 = (TH1D*) Stack_Selected_Momentum_Sample4->GetHistogram();
  hStack_Selected_Momentum_Sample4->GetXaxis()->SetTitle("p_{#mu} (GeV)");
  hStack_Selected_Momentum_Sample4->GetYaxis()->SetTitle("Number of events");
  hStack_Selected_Momentum_Sample4->GetXaxis()->SetRangeUser(0,2);
  leg_Sample4->Draw("same");
  c3->SaveAs("Stopping_MomentumSelected.pdf");

  TCanvas * c4 = new TCanvas();
  Stack_Reconstructed_Momentum_Sample4->SetTitle("INGRID stopped #mu: Interaction type with p_{#mu}");
  Stack_Reconstructed_Momentum_Sample4->Draw();
  TH1D * hStack_Reconstructed_Momentum_Sample4 = (TH1D*) Stack_Reconstructed_Momentum_Sample4->GetHistogram();
  hStack_Reconstructed_Momentum_Sample4->GetXaxis()->SetTitle("p_{#mu} (GeV)");
  hStack_Reconstructed_Momentum_Sample4->GetYaxis()->SetTitle("Number of events");
  hStack_Reconstructed_Momentum_Sample4->GetXaxis()->SetRangeUser(0,2);
  leg_Sample4->Draw("same");
  c4->SaveAs("Stopping_MomentumReconstructed.pdf");

  TCanvas * c30 = new TCanvas();
  THStack * Stack_Selected_Momentum_Sample4_Zoom = (THStack*) Stack_Selected_Momentum_Sample4->Clone("Stack_Selected_Momentum_Sample4_Zoom");
  Stack_Selected_Momentum_Sample4_Zoom->SetMinimum(0);
  Stack_Selected_Momentum_Sample4_Zoom->SetMaximum(2000);
  Stack_Selected_Momentum_Sample4_Zoom->Draw();
  TH1D * hStack_Selected_Momentum_Sample4_Zoom = (TH1D*) Stack_Selected_Momentum_Sample4_Zoom->GetHistogram();
  hStack_Selected_Momentum_Sample4_Zoom->GetXaxis()->SetRangeUser(0,1);
  //hStack_Selected_Momentum_Sample4_Zoom->GetYaxis()->SetRangeUser(0,2000);
  Stack_Selected_Momentum_Sample4_Zoom->Draw();
  leg_Sample4->Draw("same");
  c30->SaveAs("Stopping_MomentumSelected_Zoom.pdf");
  /*
    TCanvas * c40 = new TCanvas();
    THStack * Stack_Reconstructed_Momentum_Sample4_Zoom = (THStack*) Stack_Reconstructed_Momentum_Sample4->Clone("Stack_Reconstructed_Momentum_Sample4_Zoom");
    Stack_Reconstructed_Momentum_Sample4_Zoom->SetMinimum(0);
    Stack_Reconstructed_Momentum_Sample4_Zoom->SetMaximum(3000);
    Stack_Reconstructed_Momentum_Sample4_Zoom->Draw();
    TH1D * hStack_Reconstructed_Momentum_Sample4_Zoom = (TH1D*) Stack_Reconstructed_Momentum_Sample4_Zoom->GetHistogram();
    hStack_Reconstructed_Momentum_Sample4_Zoom->GetXaxis()->SetRangeUser(0,1);
    Stack_Reconstructed_Momentum_Sample4_Zoom->Draw();
    leg_Sample4->Draw("same");
    c40->SaveAs("Stopping_MomentumReconstructed_Zoom.pdf");
  */

  TCanvas * c040 = new TCanvas();
  Stack_INGRIDTrackWidth->SetTitle("CC-0pi like INGRID stopped #mu");
  Stack_INGRIDTrackWidth->Draw();
  TH1D * hStack_INGRIDTrackWidth = (TH1D*) Stack_INGRIDTrackWidth->GetHistogram();
  hStack_INGRIDTrackWidth->GetXaxis()->SetTitle("Average track width");
  hStack_INGRIDTrackWidth->GetYaxis()->SetTitle("Number of events");
  double MaxValue=std::max(hStack_INGRIDTrackWidth->GetBinContent(hStack_INGRIDTrackWidth->GetMaximumBin()),Selected_IronDistance_Data[4]->GetBinContent(Selected_IronDistance_Data[4]->GetMaximumBin()));
  hStack_INGRIDTrackWidth->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_INGRIDTrackWidth->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_INGRIDTrackWidth->Draw();

  INGRIDTrackWidth_Data->SetMarkerStyle(20);
  INGRIDTrackWidth_Data->Draw("E1same");
  //Selected_IronDistance_Data[4]->Draw("E1Psame");
  //Selected_IronDistance_Data[4]->Draw("E1same");

  leg_Sample4_Data->Draw("same");
  c040->SaveAs("INGRID_TrackWidth.pdf");



  TCanvas * c041 = new TCanvas();
  Stack_INGRIDAngleCriteria->SetTitle("CC-0pi like INGRID stopped #mu");
  Stack_INGRIDAngleCriteria->Draw();
  TH1D * hStack_INGRIDAngleCriteria = (TH1D*) Stack_INGRIDAngleCriteria->GetHistogram();
  hStack_INGRIDAngleCriteria->GetXaxis()->SetTitle("Average track width");
  hStack_INGRIDAngleCriteria->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_INGRIDAngleCriteria->GetBinContent(hStack_INGRIDAngleCriteria->GetMaximumBin()),Selected_IronDistance_Data[4]->GetBinContent(Selected_IronDistance_Data[4]->GetMaximumBin()));
  hStack_INGRIDAngleCriteria->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_INGRIDAngleCriteria->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_INGRIDAngleCriteria->Draw();
  INGRIDAngleCriteria_Data->SetMarkerStyle(20);
  INGRIDAngleCriteria_Data->Draw("E1same");
  leg_Sample4_Data->Draw("same");
  c041->SaveAs("INGRID_AngleCriteria.pdf");


  TCanvas * c042 = new TCanvas();
  Stack_INGRIDPositionCriteria->SetTitle("CC-0pi like INGRID stopped #mu");
  Stack_INGRIDPositionCriteria->Draw();
  TH1D * hStack_INGRIDPositionCriteria = (TH1D*) Stack_INGRIDPositionCriteria->GetHistogram();
  hStack_INGRIDPositionCriteria->GetXaxis()->SetTitle("Average track width");
  hStack_INGRIDPositionCriteria->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_INGRIDPositionCriteria->GetBinContent(hStack_INGRIDPositionCriteria->GetMaximumBin()),Selected_IronDistance_Data[4]->GetBinContent(Selected_IronDistance_Data[4]->GetMaximumBin()));
  hStack_INGRIDPositionCriteria->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_INGRIDPositionCriteria->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_INGRIDPositionCriteria->Draw();
  INGRIDPositionCriteria_Data->SetMarkerStyle(20);
  INGRIDPositionCriteria_Data->Draw("E1same");
  leg_Sample4_Data->Draw("same");
  c042->SaveAs("INGRID_PositionCriteria.pdf");






  TCanvas * c41 = new TCanvas();
  //c41->Divide(1,2);
  //c41->cd(1);
  Stack_Selected_IronDistance_Sample4->SetTitle("CC-0pi like INGRID stopped #mu: Interaction type with d_{#mu}");
  Stack_Selected_IronDistance_Sample4->Draw();
  TH1D * hStack_Selected_IronDistance_Sample4 = (TH1D*) Stack_Selected_IronDistance_Sample4->GetHistogram();
  hStack_Selected_IronDistance_Sample4->GetXaxis()->SetTitle("d_{#mu} (cm)");
  hStack_Selected_IronDistance_Sample4->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_Selected_IronDistance_Sample4->GetBinContent(hStack_Selected_IronDistance_Sample4->GetMaximumBin()),Selected_IronDistance_Data[4]->GetBinContent(Selected_IronDistance_Data[4]->GetMaximumBin()));
  hStack_Selected_IronDistance_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  //c41->cd(2);
  //Selected_IronDistance_Data[4]->Draw("same");
  //Selected_IronDistance_MC[4]->SetMarkerStyle(1);
  //Selected_IronDistance_MC[4]->Draw("HISTsame"); 
  Stack_Selected_IronDistance_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_Selected_IronDistance_Sample4->Draw();

  Selected_IronDistance_Data[4]->SetMarkerStyle(20);
  //Selected_IronDistance_Data[4]->Draw("E1Psame");
  Selected_IronDistance_Data[4]->Draw("E1same");

  leg_Sample4_Data->Draw("same");
  c41->SaveAs("Stopping_IronDistanceSelected.pdf");



  TCanvas * c410 = new TCanvas();
  //c41->Divide(1,2);
  //c41->cd(1);
  Stack_Selected_ReconstructedAngle_Sample4->SetTitle("CC-0pi like INGRID stopped #mu: Interaction type with #theta_{#mu}");
  Stack_Selected_ReconstructedAngle_Sample4->Draw();
  TH1D * hStack_Selected_ReconstructedAngle_Sample4 = (TH1D*) Stack_Selected_ReconstructedAngle_Sample4->GetHistogram();
  hStack_Selected_ReconstructedAngle_Sample4->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_Selected_ReconstructedAngle_Sample4->GetBinContent(hStack_Selected_ReconstructedAngle_Sample4->GetMaximumBin()),Selected_ReconstructedAngle_Data[4]->GetBinContent(Selected_ReconstructedAngle_Data[4]->GetMaximumBin()));
  hStack_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  //c41->cd(2);
  //Selected_ReconstructedAngle_Data[4]->Draw("same");
  Stack_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_Selected_ReconstructedAngle_Sample4->Draw();
  Selected_ReconstructedAngle_Data[4]->SetMarkerStyle(20);
  Selected_ReconstructedAngle_Data[4]->Draw("E1Psame");
  //Selected_ReconstructedAngle_MC[4]->SetMarkerStyle(1);
  //Stack_Selected_ReconstructedAngle_Sample4->Draw("same");
  leg_Sample4_Data->Draw("same");
  c410->SaveAs("Stopping_ReconstructedAngleSelected.pdf");



  TCanvas * c411 = new TCanvas();
  //c41->Divide(1,2);
  //c41->cd(1);
  Stack_PID_Selected_IronDistance_Sample4->SetTitle("CC-0pi like INGRID stopped #mu after PID");
  Stack_PID_Selected_IronDistance_Sample4->Draw();
  TH1D * hStack_PID_Selected_IronDistance_Sample4 = (TH1D*) Stack_PID_Selected_IronDistance_Sample4->GetHistogram();
  hStack_PID_Selected_IronDistance_Sample4->GetXaxis()->SetTitle("d_{#mu} (cm)");
  hStack_PID_Selected_IronDistance_Sample4->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_PID_Selected_IronDistance_Sample4->GetBinContent(hStack_PID_Selected_IronDistance_Sample4->GetMaximumBin()),PID_Selected_IronDistance_Data[4]->GetBinContent(PID_Selected_IronDistance_Data[4]->GetMaximumBin()));
  hStack_PID_Selected_IronDistance_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  //c41->cd(2);
  //PID_Selected_IronDistance_Data[4]->Draw("same");
  //PID_Selected_IronDistance_MC[4]->SetMarkerStyle(1);
  //PID_Selected_IronDistance_MC[4]->Draw("HISTsame"); 
  Stack_PID_Selected_IronDistance_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_PID_Selected_IronDistance_Sample4->Draw();

  PID_Selected_IronDistance_Data[4]->SetMarkerStyle(20);
  //PID_Selected_IronDistance_Data[4]->Draw("E1Psame");
  PID_Selected_IronDistance_Data[4]->Draw("E1same");

  leg_Sample4_Data->Draw("same");
  c411->SaveAs("Stopping_IronDistancePID_Selected.pdf");



  TCanvas * c412 = new TCanvas();
  //c41->Divide(1,2);
  //c41->cd(1);
  Stack_PID_Selected_ReconstructedAngle_Sample4->SetTitle("CC-0pi like INGRID stopped #mu after PID");
  Stack_PID_Selected_ReconstructedAngle_Sample4->Draw();
  TH1D * hStack_PID_Selected_ReconstructedAngle_Sample4 = (TH1D*) Stack_PID_Selected_ReconstructedAngle_Sample4->GetHistogram();
  hStack_PID_Selected_ReconstructedAngle_Sample4->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_PID_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_PID_Selected_ReconstructedAngle_Sample4->GetBinContent(hStack_PID_Selected_ReconstructedAngle_Sample4->GetMaximumBin()),PID_Selected_ReconstructedAngle_Data[4]->GetBinContent(PID_Selected_ReconstructedAngle_Data[4]->GetMaximumBin()));
  hStack_PID_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  //c41->cd(2);
  //PID_Selected_ReconstructedAngle_Data[4]->Draw("same");
  Stack_PID_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_PID_Selected_ReconstructedAngle_Sample4->Draw();
  PID_Selected_ReconstructedAngle_Data[4]->SetMarkerStyle(20);
  PID_Selected_ReconstructedAngle_Data[4]->Draw("E1Psame");
  //PID_Selected_ReconstructedAngle_MC[4]->SetMarkerStyle(1);
  //Stack_PID_Selected_ReconstructedAngle_Sample4->Draw("same");
  leg_Sample4_Data->Draw("same");
  c412->SaveAs("Stopping_ReconstructedAnglePID_Selected.pdf");




  TCanvas * c413 = new TCanvas();
  //c41->Divide(1,2);
  //c41->cd(1);
  Stack_PIDAndMatching_Selected_IronDistance_Sample4->SetTitle("CC-0pi like INGRID stopped #mu after PM/INGRID track matching");
  Stack_PIDAndMatching_Selected_IronDistance_Sample4->Draw();
  TH1D * hStack_PIDAndMatching_Selected_IronDistance_Sample4 = (TH1D*) Stack_PIDAndMatching_Selected_IronDistance_Sample4->GetHistogram();
  hStack_PIDAndMatching_Selected_IronDistance_Sample4->GetXaxis()->SetTitle("d_{#mu} (cm)");
  hStack_PIDAndMatching_Selected_IronDistance_Sample4->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_PIDAndMatching_Selected_IronDistance_Sample4->GetBinContent(hStack_PIDAndMatching_Selected_IronDistance_Sample4->GetMaximumBin()),PIDAndMatching_Selected_IronDistance_Data[4]->GetBinContent(PIDAndMatching_Selected_IronDistance_Data[4]->GetMaximumBin()));
  hStack_PIDAndMatching_Selected_IronDistance_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  //c41->cd(2);
  //PIDAndMatching_Selected_IronDistance_Data[4]->Draw("same");
  //PIDAndMatching_Selected_IronDistance_MC[4]->SetMarkerStyle(1);
  //PIDAndMatching_Selected_IronDistance_MC[4]->Draw("HISTsame"); 
  Stack_PIDAndMatching_Selected_IronDistance_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_PIDAndMatching_Selected_IronDistance_Sample4->Draw();

  PIDAndMatching_Selected_IronDistance_Data[4]->SetMarkerStyle(20);
  //PIDAndMatching_Selected_IronDistance_Data[4]->Draw("E1Psame");
  PIDAndMatching_Selected_IronDistance_Data[4]->Draw("E1same");

  leg_Sample4_Data->Draw("same");
  c413->SaveAs("Stopping_IronDistancePIDAndMatching_Selected.pdf");



  TCanvas * c414 = new TCanvas();
  //c41->Divide(1,2);
  //c41->cd(1);
  Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->SetTitle("CC-0pi like INGRID stopped #mu after PM/INGRID track matching");
  Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->Draw();
  TH1D * hStack_PIDAndMatching_Selected_ReconstructedAngle_Sample4 = (TH1D*) Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetHistogram();
  hStack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetTitle("Number of events");
  MaxValue=std::max(hStack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetBinContent(hStack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetMaximumBin()),PIDAndMatching_Selected_ReconstructedAngle_Data[4]->GetBinContent(PIDAndMatching_Selected_ReconstructedAngle_Data[4]->GetMaximumBin()));
  hStack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  //c41->cd(2);
  //PIDAndMatching_Selected_ReconstructedAngle_Data[4]->Draw("same");
  Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->GetYaxis()->SetRangeUser(0,MaxValue+0.1*MaxValue);
  Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->Draw();
  PIDAndMatching_Selected_ReconstructedAngle_Data[4]->SetMarkerStyle(20);
  PIDAndMatching_Selected_ReconstructedAngle_Data[4]->Draw("E1Psame");
  //PIDAndMatching_Selected_ReconstructedAngle_MC[4]->SetMarkerStyle(1);
  //Stack_PIDAndMatching_Selected_ReconstructedAngle_Sample4->Draw("same");
  leg_Sample4_Data->Draw("same");
  c414->SaveAs("Stopping_ReconstructedAnglePIDAndMatching_Selected.pdf");


  TCanvas * c415 = new TCanvas();
  Selected2D->Scale(ScalingFactor);
  Selected2D->SetTitle("MC: CC-0pi like INGRID stopped #mu #theta_{#mu}Xd_{#mu} distribution ");
  Selected2D->GetXaxis()->SetTitle("d_{#mu} (cm)");
  Selected2D->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
  Selected2D->GetZaxis()->SetTitle("Number of events");
  Selected2D->GetXaxis()->SetTitleOffset(1.5);
  Selected2D->GetYaxis()->SetTitleOffset(1.8);
  Selected2D->GetZaxis()->SetTitleOffset(1.2);
  Selected2D->GetXaxis()->SetRangeUser(0,80);
  Selected2D->GetYaxis()->SetRangeUser(0,30);
  Selected2D->Draw("LEGO1");
  Selected2D->SaveAs("SelectedMC.pdf");

  TCanvas * c416 = new TCanvas();
  Selected2D_Data->SetTitle("Data: CC-0pi like INGRID stopped #mu #theta_{#mu}Xd_{#mu} distribution ");
  Selected2D_Data->GetXaxis()->SetTitle("d_{#mu} (cm)");
  Selected2D_Data->GetYaxis()->SetTitle("#theta_{#mu} (#circ)");
  Selected2D_Data->GetZaxis()->SetTitle("Number of events");
  Selected2D_Data->GetXaxis()->SetTitleOffset(1.5);
  Selected2D_Data->GetYaxis()->SetTitleOffset(1.8);
  Selected2D_Data->GetZaxis()->SetTitleOffset(1.2);
  Selected2D_Data->GetXaxis()->SetRangeUser(0,80);
  Selected2D_Data->GetYaxis()->SetRangeUser(0,30);
  Selected2D_Data->Draw("LEGO1");
  Selected2D_Data->SaveAs("SelectedData.pdf");

  TCanvas * c417 = new TCanvas();
  TH2D * Selected2D_Data2 = (TH2D*) Selected2D_Data->Clone("Selected2D_Data2");
  cout.precision(2);
  for(int ibinx=1;ibinx<=Selected2D_Data->GetNbinsX();ibinx++){
    cout<<std::fixed<<"\\multirow{3}{*}{$d_{\\mu} \\in ["<<DistIronBin[ibinx-1]<<"~cm, "<<DistIronBin[ibinx]<<"~cm]$} & ";
    for(int ibiny=1;ibiny<=Selected2D_Data->GetNbinsY();ibiny++){
      double MC=Selected2D->GetBinContent(ibinx,ibiny);
      double Data=Selected2D_Data->GetBinContent(ibinx,ibiny);
      double ErrorData=TMath::Sqrt(Data);
      double Relative=Data-MC;
      if(MC!=0) Relative/=MC;
      Selected2D_Data2->SetBinContent(ibinx,ibiny,Relative);
      if(ibiny!=Selected2D_Data->GetNbinsY()) cout<<"Data: "<<Data<<" \\pm "<<ErrorData<<"& ";
      else cout<<"Data: "<<Data<<" \\pm "<<ErrorData<<" \\"<<"\\"<<endl;
    }
    cout<<" & ";
    for(int ibiny=1;ibiny<=Selected2D_Data->GetNbinsY();ibiny++){
      double MC=Selected2D->GetBinContent(ibinx,ibiny);
      if(ibiny!=Selected2D_Data->GetNbinsY()) cout<<"MC: "<<MC<<" & ";
      else cout<<"MC: "<<MC<<" \\"<<"\\"<<endl;
      double Data=Selected2D_Data->GetBinContent(ibinx,ibiny);
      double Relative=Data-MC;
      if(MC!=0) Relative/=MC;
    }
    cout<<" & ";
    for(int ibiny=1;ibiny<=Selected2D_Data->GetNbinsY();ibiny++){
      double MC=Selected2D->GetBinContent(ibinx,ibiny);
      double Data=Selected2D_Data->GetBinContent(ibinx,ibiny);
      double ErrorData=TMath::Sqrt(Data);
      double Error=ErrorData; if(MC!=0) Error/=MC;
      double Ratio=Data-MC;
      if(MC!=0) Ratio/=MC;
      if(ibiny!=Selected2D_Data->GetNbinsY()) cout<<"$\\frac{Data-MC}{MC}$: "<<Ratio*100<<" \\pm "<<Error*100<<" \\% & ";
      else cout<<"$\\frac{Data-MC}{MC}$: "<<Ratio*100<<" \\pm "<<Error*100<<"\\% \\"<<"\\"<<endl;
    }
    cout<<"\\hline"<<endl;
  }
  Selected2D_Data2->Draw("colz");


  TCanvas * c42 = new TCanvas();
  Stack_Reconstructed_IronDistance_Sample4->SetTitle("INGRID stopped #mu: Interaction type with p_{#mu}");
  Stack_Reconstructed_IronDistance_Sample4->Draw();
  TH1D * hStack_Reconstructed_IronDistance_Sample4 = (TH1D*) Stack_Reconstructed_IronDistance_Sample4->GetHistogram();
  hStack_Reconstructed_IronDistance_Sample4->GetXaxis()->SetTitle("d_{#mu} (cm)");
  hStack_Reconstructed_IronDistance_Sample4->GetYaxis()->SetTitle("Number of events");
  Reconstructed_IronDistance_Data[4]->Draw("same");
  leg_Sample4->Draw("same");
  c42->SaveAs("Stopping_IronDistanceReconstructed.pdf");


  TCanvas * c43 = new TCanvas();
  Stack_Reconstructed_Angle_Sample4->SetTitle("INGRID escaping #mu: Interaction type with #theta_{#mu}");
  Stack_Reconstructed_Angle_Sample4->Draw();
  
  TH1D * hStack_Reconstructed_Angle_Sample4 = (TH1D*) Stack_Reconstructed_Angle_Sample4->GetHistogram();
  hStack_Reconstructed_Angle_Sample4->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_Reconstructed_Angle_Sample4->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c43->SaveAs("Stopping_AngleReconstructed.pdf");

  TCanvas * c5 = new TCanvas();
  Stack_Selected_NTracks_Sample5->SetTitle("CC-0pi like INGRID escaping #mu: Interaction type with reconstructed tracks");
  //(Stack_Selected_NTracks_Sample5->GetXaxis())->SetTitle("p_{#mu} (GeV)");
  //(Stack_Selected_NTracks_Sample5->GetYaxis())->SetTitle("Number of events");
  Stack_Selected_NTracks_Sample5->Draw();
  TH1D * hStack_Selected_NTracks_Sample5 = (TH1D*) Stack_Selected_NTracks_Sample5->GetHistogram();
  hStack_Selected_NTracks_Sample5->GetXaxis()->SetTitle("Number of reconstructed tracks");
  hStack_Selected_NTracks_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c5->SaveAs("ThroughGoing_NTracksSelected.pdf");

  TCanvas * c6 = new TCanvas();
  Stack_Reconstructed_NTracks_Sample5->SetTitle("INGRID escaping #mu: Interaction type with reconstructed tracks");
  Stack_Reconstructed_NTracks_Sample5->Draw();
  TH1D * hStack_Reconstructed_NTracks_Sample5 = (TH1D*) Stack_Reconstructed_NTracks_Sample5->GetHistogram();
  hStack_Reconstructed_NTracks_Sample5->GetXaxis()->SetTitle("Number of reconstructed tracks");
  hStack_Reconstructed_NTracks_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c6->SaveAs("ThroughGoing_NTracksReconstructed.pdf");

  TCanvas * c7 = new TCanvas();
  Stack_Selected_Momentum_Sample5->SetTitle("CC-0pi like INGRID escaping #mu: Interaction type with p_{#mu}");
  Stack_Selected_Momentum_Sample5->Draw();
  TH1D * hStack_Selected_Momentum_Sample5 = (TH1D*) Stack_Selected_Momentum_Sample5->GetHistogram();
  hStack_Selected_Momentum_Sample5->GetXaxis()->SetTitle("p_{#mu} (GeV)");
  hStack_Selected_Momentum_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c7->SaveAs("ThroughGoing_MomentumSelected.pdf");

  TCanvas * c8 = new TCanvas();
  Stack_Reconstructed_Momentum_Sample5->SetTitle("INGRID escaping #mu: Interaction type with p_{#mu}");
  Stack_Reconstructed_Momentum_Sample5->Draw();
  TH1D * hStack_Reconstructed_Momentum_Sample5 = (TH1D*) Stack_Reconstructed_Momentum_Sample5->GetHistogram();
  hStack_Reconstructed_Momentum_Sample5->GetXaxis()->SetTitle("p_{#mu} (GeV)");
  hStack_Reconstructed_Momentum_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c8->SaveAs("ThroughGoing_MomentumReconstructed.pdf");

  TCanvas * c80 = new TCanvas();
  Stack_Reconstructed_Angle_Sample5->SetTitle("INGRID escaping #mu: Interaction type with #theta_{#mu}");
  Stack_Reconstructed_Angle_Sample5->Draw();
  TH1D * hStack_Reconstructed_Angle_Sample5 = (TH1D*) Stack_Reconstructed_Angle_Sample5->GetHistogram();
  hStack_Reconstructed_Angle_Sample5->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_Reconstructed_Angle_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c80->SaveAs("ThroughGoing_AngleReconstructed.pdf");

  TCanvas * c81 = new TCanvas();
  Stack_Reconstructed_ReconstructedAngle_Sample5->SetTitle("INGRID escaping #mu: Interaction type with reconstructed #theta_{#mu}");
  Stack_Reconstructed_ReconstructedAngle_Sample5->Draw();
  TH1D * hStack_Reconstructed_ReconstructedAngle_Sample5 = (TH1D*) Stack_Reconstructed_ReconstructedAngle_Sample5->GetHistogram();
  hStack_Reconstructed_ReconstructedAngle_Sample5->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_Reconstructed_ReconstructedAngle_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c81->SaveAs("ThroughGoing_ReconstructedAngleReconstructed.pdf");

  TCanvas * c82 = new TCanvas();
  //c82->Divide(1,2);
  //c82->cd(1);
  Stack_Selected_Angle_Sample5->SetTitle("INGRID escaping #mu: Interaction type with #theta_{#mu}");
  Stack_Selected_Angle_Sample5->Draw();
  TH1D * hStack_Selected_Angle_Sample5 = (TH1D*) Stack_Selected_Angle_Sample5->GetHistogram();
  hStack_Selected_Angle_Sample5->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_Selected_Angle_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  c82->SaveAs("ThroughGoing_AngleSelected.pdf");

  TCanvas * c83 = new TCanvas();
  //c82->Divide(1,2);
  //c83->cd(1);
  Stack_Selected_ReconstructedAngle_Sample5->SetTitle("INGRID escaping #mu: Interaction type with reconstructed #theta_{#mu}");
  Stack_Selected_ReconstructedAngle_Sample5->Draw();
  TH1D * hStack_Selected_ReconstructedAngle_Sample5 = (TH1D*) Stack_Selected_ReconstructedAngle_Sample5->GetHistogram();
  hStack_Selected_ReconstructedAngle_Sample5->GetXaxis()->SetTitle("#theta_{#mu} (#circ)");
  hStack_Selected_ReconstructedAngle_Sample5->GetYaxis()->SetTitle("Number of events");
  leg_Sample4->Draw("same");
  //c83->cd(2);
  Selected_ReconstructedAngle_Data[5]->SetMarkerStyle(20);
  Selected_ReconstructedAngle_Data[5]->Draw("Psame");
  c83->SaveAs("ThroughGoing_ReconstructedAngleSelected.pdf");

  /*
    TCanvas * c99 = new TCanvas();
    hStack_Selected_Momentum_Sample4->Draw();
    cout<<"Stopping sample:"<<endl;
    cout<<"Original Purity in CC0pi="<<Reconstructed_NTracks[4][1]->Integral()/Reconstructed_Momentum_Sample4_All->Integral()<<endl;
    cout<<"Purity after selection in CC0pi="<<Selected_NTracks[4][1]->Integral()/Selected_Momentum_Sample4_All->Integral()<<endl;
    cout<<"Cut efficiency="<<Selected_NTracks[4][1]->Integral()/Reconstructed_NTracks[4][1]->Integral()<<endl;
    cout<<"if removing NTracks>2:"<<endl;
    cout<<"Purity after selection in CC0pi="<<Selected_NTracks[4][1]->Integral(0,3)/Selected_NTracks_Sample4_All->Integral(0,3)<<endl;
    cout<<"Cut efficiency="<<Selected_NTracks[4][1]->Integral(0,3)/Reconstructed_NTracks[4][1]->Integral()<<endl;
    cout<<"test bin="<<Selected_NTracks[4][1]->GetBinCenter(3)<<", "<<Selected_NTracks[4][1]->Integral(2,3)<<", "<<Selected_NTracks[4][1]->Integral(0,10)<<endl;


    cout<<endl<<endl<<"INGRID escaping sample:"<<endl;
    cout<<"Original Purity in CC0pi="<<Reconstructed_NTracks[5][1]->Integral()/Reconstructed_Momentum_Sample5_All->Integral()<<endl;
    cout<<"Purity after selection in CC0pi="<<Selected_NTracks[5][1]->Integral()/Selected_Momentum_Sample5_All->Integral()<<endl;
    cout<<"Cut efficiency="<<Selected_NTracks[5][1]->Integral()/Reconstructed_NTracks[5][1]->Integral()<<endl;
    cout<<"if removing NTracks>2:"<<endl;
    cout<<"Purity after selection in CC0pi="<<Selected_NTracks[5][1]->Integral(0,3)/Selected_NTracks_Sample5_All->Integral(0,3)<<endl;
    cout<<"Cut efficiency="<<Selected_NTracks[5][1]->Integral(0,3)/Reconstructed_NTracks[5][1]->Integral()<<endl;
    cout<<"test bin="<<Selected_NTracks[5][1]->GetBinCenter(3)<<", "<<Selected_NTracks[5][1]->Integral(2,3)<<", "<<Selected_NTracks[5][1]->Integral(0,10)<<endl;

    TH1D * Selected_Higher_Purity = (TH1D*) Selected_Higher[0]->Clone("Selected_Higher_Purity");
    Selected_Higher_Purity->Sumw2();
    Selected_Higher_Purity->GetXaxis()->SetTitle("#mu_{CL} cut value");
    Selected_Higher_Purity->GetYaxis()->SetTitle("#mu purity");
    Selected_Higher_Purity->GetYaxis()->SetTitleOffset(1.3);
    Selected_Higher_Purity->GetXaxis()->SetRangeUser(0,1);

    TH1D * Selected_Higher_Purity_All = (TH1D*) Selected_Higher[0]->Clone("Selected_Higher_Purity_All");
    Selected_Higher_Purity_All->Sumw2();
    for(int type=1;type<4;type++) Selected_Higher_Purity_All->Add(Selected_Higher[type]);
    Selected_Higher_Purity->Divide(Selected_Higher_Purity_All);

    TH1D * Selected_Higher_Efficiency = (TH1D*) Selected_Higher[0]->Clone("Selected_Higher_Efficiency");
    Selected_Higher_Efficiency->Sumw2();
    Selected_Higher_Efficiency->GetXaxis()->SetTitle("#mu_{CL} cut value");
    Selected_Higher_Efficiency->GetYaxis()->SetTitle("#mu efficiency");
    Selected_Higher_Efficiency->GetYaxis()->SetTitleOffset(1.3);
    Selected_Higher_Efficiency->Divide(Interacting_Higher[0]);
    Selected_Higher_Efficiency->GetXaxis()->SetRangeUser(0,1);

    TCanvas * c9 = new TCanvas();
    Selected_Higher_Purity->Draw();
    c9->SaveAs("Muon_purity.pdf");
    TCanvas * c10 = new TCanvas();
    Selected_Higher_Efficiency->Draw();
    c10->SaveAs("Muon_efficiency.pdf");


    TH1D * Selected_Lower_Purity = (TH1D*) Selected_Lower[2]->Clone("Selected_Lower_Purity");
    Selected_Lower_Purity->Sumw2();
    Selected_Lower_Purity->GetXaxis()->SetTitle("#mu_{CL} cut value");
    Selected_Lower_Purity->GetYaxis()->SetTitle("p purity");
    Selected_Lower_Purity->GetYaxis()->SetTitleOffset(1.3);
    Selected_Lower_Purity->GetXaxis()->SetRangeUser(0,1);

    TH1D * Selected_Lower_Purity_All = (TH1D*) Selected_Lower[2]->Clone("Selected_Lower_Purity_All");
    Selected_Lower_Purity_All->Sumw2();
    for(int type=0;type<4;type++){ if(type!=2) Selected_Lower_Purity_All->Add(Selected_Lower[type]);}
    Selected_Lower_Purity->Divide(Selected_Lower_Purity_All);

    TH1D * Selected_Lower_Efficiency = (TH1D*) Selected_Lower[2]->Clone("Selected_Lower_Efficiency");
    Selected_Lower_Efficiency->Sumw2();
    Selected_Lower_Efficiency->GetXaxis()->SetTitle("#mu_{CL} cut value");
    Selected_Lower_Efficiency->GetYaxis()->SetTitle("p efficiency");
    Selected_Lower_Efficiency->GetYaxis()->SetTitleOffset(1.3);
    Selected_Lower_Efficiency->Divide(Interacting_Lower[2]);
    Selected_Lower_Efficiency->GetXaxis()->SetRangeUser(0,1);
  
    TCanvas * c11 = new TCanvas();
    Selected_Lower_Purity->Draw();
    c11->SaveAs("Proton_purity.pdf");
    TCanvas * c12 = new TCanvas();
    Selected_Lower_Efficiency->Draw();
    c12->SaveAs("Proton_efficiency.pdf");
  */
  /*
    TCanvas * c11 = new TCanvas();
    Selected_Muon_Purity_Sample5->Draw();
    c11->SaveAs("CC-0pi_purity_MuonGoing.pdf");
    TCanvas * c12 = new TCanvas();
    Selected_Muon_Efficiency_Sample5->Draw();
    c12->SaveAs("CC-0pi_efficiency_MuonGoing.pdf");

  */

  /*
    TFile * wfile = new TFile("UnfoldingInputs.root","RECREATE");
    wfile->Write();
    wfile->Close();*/
  /*
   */
  TCanvas * c00 = new TCanvas();
  NTracksTotal_Data->SetLineColor(2.);
  NTracksTotal_Data->Draw();
  NTracksTotal->SetLineColor(4.);
  NTracksTotal->Scale(ScalingFactor);
  NTracksTotal->Draw("same");
  TCanvas * c01 = new TCanvas();
  AngleTotal_Data->SetLineColor(2.);
  AngleTotal_Data->Draw();
  AngleTotal->SetLineColor(4.);
  AngleTotal->Scale(ScalingFactor);
  AngleTotal->Draw("same");
  TCanvas * c02 = new TCanvas();
  CLMuonTotal_Data->SetLineColor(2.);
  CLMuonTotal_Data->Draw();
  CLMuonTotal->SetLineColor(4.);
  CLMuonTotal->Scale(ScalingFactor);
  CLMuonTotal->Draw("same");
  TCanvas * c03 = new TCanvas();
  SampleTotal_Data->SetLineColor(2.);
  SampleTotal_Data->Draw();
  SampleTotal->SetLineColor(4.);
  SampleTotal->Scale(ScalingFactor);
  SampleTotal->Draw("same");
  TCanvas * c04 = new TCanvas();
  IronDistanceTotal_Data->SetLineColor(2.);
  IronDistanceTotal_Data->Draw();
  IronDistanceTotal->SetLineColor(4.);
  IronDistanceTotal->Scale(ScalingFactor);
  IronDistanceTotal->Draw("same");
  TCanvas * c05 = new TCanvas();
  NTracksTotalSample4_Data->SetLineColor(2.);
  NTracksTotalSample4_Data->Draw();
  NTracksTotalSample4->SetLineColor(4.);
  NTracksTotalSample4->Scale(ScalingFactor);
  NTracksTotalSample4->Draw("same");
  TCanvas * c06 = new TCanvas();
  SampleTotalSample4_Data->SetLineColor(2.);
  SampleTotalSample4_Data->Draw();
  SampleTotalSample4->SetLineColor(4.);
  SampleTotalSample4->Scale(ScalingFactor);
  SampleTotalSample4->Draw("same");
  TCanvas * c07 = new TCanvas();
  IronDistanceTotalSample4_1Track_Data->SetLineColor(2.);
  IronDistanceTotalSample4_1Track_Data->Draw();
  IronDistanceTotalSample4_1Track->SetLineColor(4.);
  IronDistanceTotalSample4_1Track->Scale(ScalingFactor);
  IronDistanceTotalSample4_1Track->Draw("same");
  TCanvas * c08 = new TCanvas();
  IronDistanceTotalSample4_2Tracks_Data->SetLineColor(2.);
  IronDistanceTotalSample4_2Tracks_Data->Draw();
  IronDistanceTotalSample4_2Tracks->SetLineColor(4.);
  IronDistanceTotalSample4_2Tracks->Scale(ScalingFactor);
  IronDistanceTotalSample4_2Tracks->Draw("same");
  TCanvas * c09[5];

  for(int i=0;i<6;i++){
    c09[i] = new TCanvas();
    if(i==1) continue;
    //c09->Divide(2,3);
    //c09->cd(i+1);
    CLMuon_Distribution[i]->Sumw2();
    CLMuon_Distribution[i]->Scale(ScalingFactor);
    CLMuon_Distribution[i]->SetLineColor(1);
    CLMuon_Distribution[i]->SetLineWidth(2);
    CLMuon_Distribution[i]->SetFillColor(4);
    CLMuon_Distribution[i]->SetFillStyle(3001);
    //CLMuon_Distribution_Data[i]->SetLineColor(2);
    //CLMuon_Distribution_Data[i]->Draw();
    CLMuon_Distribution[i]->Draw("HISTE1");
    if(i==0) CLMuon_Distribution[i]->SetTitle("PM stopped");
    //if(i==1) CLMuon_Distribution[i]->SetTitle("INGRID escaping");
    if(i==2) CLMuon_Distribution[i]->SetTitle("PM escaping");
    if(i==3) CLMuon_Distribution[i]->SetTitle("INGRID prematurely stopped");
    if(i==4) CLMuon_Distribution[i]->SetTitle("INGRID stopped");
    if(i==5) CLMuon_Distribution[i]->SetTitle("INGRID escaping");
    CLMuon_Distribution[i]->GetXaxis()->SetTitle("#mu_{CL}");
    CLMuon_Distribution[i]->GetYaxis()->SetTitle("Number of events");
    CLMuon_Distribution[i]->GetXaxis()->SetRangeUser(0,.99);
    sprintf(Name,"CLMuon_Topology%d.pdf",i);
    c09[i]->SaveAs(Name);
  }

  TCanvas * cCriteria[4];
  cCriteria[0] = new TCanvas();
  cCriteria[0]->Divide(1,2);
  cCriteria[0]->cd(1);
  hCriteriaAngleX_Data[4]->Draw("colz");
  cCriteria[0]->cd(2);
  hCriteriaAngleX[4]->Draw("colz");

  cCriteria[1] = new TCanvas();
  cCriteria[1]->Divide(1,2);
  cCriteria[1]->cd(1);
  hCriteriaAngleY_Data[4]->Draw("colz");
  cCriteria[1]->cd(2);
  hCriteriaAngleY[4]->Draw("colz");

  cCriteria[2] = new TCanvas();
  cCriteria[2]->Divide(1,2);
  cCriteria[2]->cd(1);
  hCriteriaHalfWayX_Data[4]->Draw("colz");
  cCriteria[2]->cd(2);
  hCriteriaHalfWayX[4]->Draw("colz");

  cCriteria[3] = new TCanvas();
  cCriteria[3]->Divide(1,2);
  cCriteria[3]->cd(1);
  hCriteriaHalfWayY_Data[4]->Draw("colz");
  cCriteria[3]->cd(2);
  hCriteriaHalfWayY[4]->Draw("colz");

  TCanvas * cCriteria_Angle = new TCanvas();
  cCriteria_Angle->Divide(3,3);
  for(int i=0;i<9;i++){
    cCriteria_Angle->cd(i+1);
    hCriteriaAngleXY[i]->Draw("colz");
  }

  TCanvas * cCriteria_Angle_Data = new TCanvas();
  cCriteria_Angle_Data->Divide(3,3);
  for(int i=0;i<9;i++){
    cCriteria_Angle_Data->cd(i+1);
    hCriteriaAngleXY_Data[i]->Draw("colz");
  }

  TCanvas * cCriteria_HalfWay = new TCanvas();
  cCriteria_HalfWay->Divide(3,3);
  for(int i=0;i<9;i++){
    cCriteria_HalfWay->cd(i+1);
    hCriteriaHalfWayXY[i]->Draw("colz");
  }

  TCanvas * cCriteria_HalfWay_Data = new TCanvas();
  cCriteria_HalfWay_Data->Divide(3,3);
  for(int i=0;i<9;i++){
    cCriteria_HalfWay_Data->cd(i+1);
    hCriteriaHalfWayXY_Data[i]->Draw("colz");
  }

  TCanvas * cPlasticDistance = new TCanvas();
  Problematic_PlasticDistance->Scale(ScalingFactor);
  Problematic_PlasticDistance->SetLineColor(4);
  Problematic_PlasticDistance_Data->SetLineColor(2);
  Problematic_PlasticDistance_Data->Draw();
  Problematic_PlasticDistance->Draw("same");

  TCanvas * cAngle = new TCanvas();
  Problematic_Angle->Scale(ScalingFactor);
  Problematic_Angle->SetLineColor(4);
  Problematic_Angle_Data->SetLineColor(2);
  Problematic_Angle_Data->Draw();
  Problematic_Angle->Draw("same");

  TCanvas * cNTracks = new TCanvas();
  Problematic_NTracks->Scale(ScalingFactor);
  Problematic_NTracks->SetLineColor(4);
  Problematic_NTracks_Data->SetLineColor(2);
  Problematic_NTracks_Data->Draw();
  Problematic_NTracks->Draw("same");

  TCanvas * cCLMuon = new TCanvas();
  Problematic_CLMuon->Scale(ScalingFactor);
  Problematic_CLMuon->SetLineColor(4);
  Problematic_CLMuon_Data->SetLineColor(2);
  Problematic_CLMuon_Data->Draw();
  Problematic_CLMuon->Draw("same");

  TCanvas * cTrackWidth = new TCanvas();
  Problematic_TrackWidth->Scale(ScalingFactor);
  Problematic_TrackWidth->SetLineColor(4);
  Problematic_TrackWidth_Data->SetLineColor(2);
  Problematic_TrackWidth_Data->Draw();
  Problematic_TrackWidth->Draw("same");

  TCanvas * cGeom = new TCanvas();
  Problematic_Geom->Scale(ScalingFactor);
  Problematic_Geom->SetLineColor(4);
  Problematic_Geom_Data->SetLineColor(2);
  Problematic_Geom_Data->Draw();
  Problematic_Geom->Draw("same");

  TCanvas * cLastChannel = new TCanvas();
  Problematic_LastChannel->Scale(ScalingFactor);
  Problematic_LastChannel->SetLineColor(4);
  Problematic_LastChannel_Data->SetLineColor(2);
  Problematic_LastChannel_Data->Draw();
  Problematic_LastChannel->Draw("same");


  TrackWidthXIronDistance->SetTitle("#mu: Average transverse track width with d_{#mu}");
  TrackWidthXIronDistance->GetXaxis()->SetTitle("d_{#mu} (cm)");
  TrackWidthXIronDistance->GetYaxis()->SetTitle("Average transverse track width (#channels)");

  TrackWidthXIronDistance_Data->SetTitle("Data: Average transverse track width with d_{#mu}");
  TrackWidthXIronDistance_Data->GetXaxis()->SetTitle("d_{#mu} (cm)");
  TrackWidthXIronDistance_Data->GetYaxis()->SetTitle("Average transverse track width (#channels)");

  TrackWidthXIronDistance_NuE->SetTitle("e^{-}: Average transverse track width with d_{#mu}");
  TrackWidthXIronDistance_NuE->GetXaxis()->SetTitle("d_{#mu} (cm)");
  TrackWidthXIronDistance_NuE->GetYaxis()->SetTitle("Average transverse track width (#channels)");

  TCanvas * cTrackWidthXIronDistance = new TCanvas();
  TProfile * pTrackWidthXIronDistance = (TProfile*) TrackWidthXIronDistance->ProfileX("pTrackWidthXIronDistance");
  TProfile * pTrackWidthXIronDistance_Data = (TProfile*) TrackWidthXIronDistance_Data->ProfileX("pTrackWidthXIronDistance_Data");

  pTrackWidthXIronDistance->SetLineColor(4);
  pTrackWidthXIronDistance->SetLineWidth(2);
  pTrackWidthXIronDistance_Data->SetLineColor(2);
  pTrackWidthXIronDistance_Data->SetLineWidth(2);
  pTrackWidthXIronDistance_Data->SetTitle("Average transverse track width with d_{#mu}");
  pTrackWidthXIronDistance_Data->GetYaxis()->SetRangeUser(0.5,2);
  pTrackWidthXIronDistance_Data->GetYaxis()->SetTitle("Average transverse track width (#channels)");
  pTrackWidthXIronDistance_Data->Draw();
  pTrackWidthXIronDistance->Draw("same");
  TLegend * legDataMC = new TLegend(0.7,0.7,0.94,0.94);
  legDataMC->SetFillColor(0);
  legDataMC->AddEntry(pTrackWidthXIronDistance_Data,"Data");
  legDataMC->AddEntry(pTrackWidthXIronDistance,"MC");
  legDataMC->Draw("same");

  TCanvas * cTrackWidthXIronDistance2D = new TCanvas();
  cTrackWidthXIronDistance2D->Divide(1,2);
  cTrackWidthXIronDistance2D->cd(1);
  TrackWidthXIronDistance_Data->Draw("colz");
  cTrackWidthXIronDistance2D->cd(2);
  TrackWidthXIronDistance->Draw("colz");

  TCanvas * cTrackWidthXIronDistance_NuE = new TCanvas();
  TProfile * pTrackWidthXIronDistance_NuE = (TProfile*) TrackWidthXIronDistance_NuE->ProfileX("pTrackWidthXIronDistance_NuE");
  pTrackWidthXIronDistance_NuE->SetLineColor(kGreen-2);
  pTrackWidthXIronDistance_NuE->SetLineWidth(2);
  pTrackWidthXIronDistance->SetLineColor(4);
  pTrackWidthXIronDistance->SetLineWidth(2);
  pTrackWidthXIronDistance->SetTitle("Average transverse track width with d_{#mu}");
  pTrackWidthXIronDistance->GetXaxis()->SetTitle("d_{#mu} (cm)");
  pTrackWidthXIronDistance->GetYaxis()->SetTitle("Average transverse track width (#channels)");
  TLegend * legNuE = new TLegend(0.7,0.7,0.94,0.94);
  legNuE->SetFillColor(0);
  pTrackWidthXIronDistance->GetYaxis()->SetRangeUser(0,3);
  legNuE->AddEntry(pTrackWidthXIronDistance,"#mu");
  legNuE->AddEntry(pTrackWidthXIronDistance_NuE,"e^{-}");
  pTrackWidthXIronDistance->Draw();
  pTrackWidthXIronDistance_NuE->Draw("same");
  legNuE->Draw("same");

  TH1D * MuWidth[3];
  MuWidth[0] = new TH1D("MuWidth[0]","",TrackWidthXIronDistance->GetNbinsY(),0,10);
  MuWidth[1] = new TH1D("MuWidth[1]","",TrackWidthXIronDistance->GetNbinsY(),0,10);
  MuWidth[2] = new TH1D("MuWidth[2]","",TrackWidthXIronDistance->GetNbinsY(),0,10);
  TH1D * EWidth[3];
  EWidth[0] = new TH1D("EWidth[0]","",TrackWidthXIronDistance->GetNbinsY(),0,10);
  EWidth[1] = new TH1D("EWidth[1]","",TrackWidthXIronDistance->GetNbinsY(),0,10);
  EWidth[2] = new TH1D("EWidth[2]","",TrackWidthXIronDistance->GetNbinsY(),0,10);

  TCanvas * canNuE[3];

  //TLine *line = new TLine(,0,8,0);

  for(int ibinx=1;ibinx<=TrackWidthXIronDistance->GetNbinsX();ibinx++){
    int Plot;
    if(TrackWidthXIronDistance_NuE->GetXaxis()->GetBinCenter(ibinx)<=20) Plot=0;
    else if(TrackWidthXIronDistance_NuE->GetXaxis()->GetBinCenter(ibinx)<=30) Plot=1;
    else Plot=2;
    //canNuE[Plot]->cd();
    MuWidth[Plot]->Add((TH1D*) TrackWidthXIronDistance->ProjectionY(Name,ibinx,ibinx),1);
    EWidth[Plot]->Add((TH1D*) TrackWidthXIronDistance_NuE->ProjectionY(Name,ibinx,ibinx),1);
  }
  for(int ip=0;ip<3;ip++){
    canNuE[ip] = new TCanvas();
    MuWidth[ip]->SetLineColor(4);
    EWidth[ip]->SetLineColor(kGreen-2);
    MuWidth[ip]->SetLineWidth(2);
    EWidth[ip]->SetLineWidth(2);
    MuWidth[ip]->SetTitle("Average track width distribution (#channels)");
    MuWidth[ip]->GetYaxis()->SetTitle("Event probability");
    MuWidth[ip]->GetXaxis()->SetTitle("Average transverse track width (#channels)");
    //MuWidth[ip]->GetYaxis()->SetRangeUser(0.,1.);
    MuWidth[ip]->GetXaxis()->SetRangeUser(1.,3.);
    MuWidth[ip]->DrawNormalized();
    EWidth[ip]->DrawNormalized("same");
    legNuE->Draw("same");
    sprintf(Name,"TrackWidth_IronBin%d.pdf",ip);
    canNuE[ip]->SaveAs(Name);
  }

  /*
    cTrackWidthXIronDistance_NuE->Divide(1,2);
    cTrackWidthXIronDistance_NuE->cd(1);
    pTrackWidthXIronDistance_NuE->Draw();
    cTrackWidthXIronDistance_NuE->cd(2);
    TrackWidthXIronDistance_NuE->Draw("colz");
  */
  TCanvas * cAngle2D = new TCanvas();
  cAngle2D->Divide(1,2);
  cAngle2D->cd(1);
  AngleRecXTrue->Draw("col");
  cAngle2D->cd(2);
  AngleRecXTrue_Selected->Draw("col");

  TCanvas * cIronXMom = new TCanvas();
  MomentumXIron_ReconstructedForTrueMuon->Draw("colz");
  /*
    TCanvas * cCriteria_Angle[NBinsIron];
    for(int i=0;i<NBinsIron;i++){
    cCriteria_Angle[i] = new TCanvas();
    cCriteria_Angle[i]->Divide(1,2);
    cCriteria_Angle[i]->cd(1);
    hCriteriaAngleXY_Data[i]->Draw("colz");
    cCriteria_Angle[i]->cd(2);
    hCriteriaAngleXY[i]->Draw("colz");
    }

    TCanvas * cCriteria_HalfWay[NBinsIron];
    for(int i=0;i<NBinsIron;i++){
    cCriteria_HalfWay[i] = new TCanvas();
    cCriteria_HalfWay[i]->Divide(1,2);
    cCriteria_HalfWay[i]->cd(1);
    hCriteriaHalfWayXY_Data[i]->Draw("colz");
    cCriteria_HalfWay[i]->cd(2);
    hCriteriaHalfWayXY[i]->Draw("colz");
    }
  */
  

  TCanvas * c010[6]; //= new TCanvas();
  TLegend * lMUCL = new TLegend(0.75,0.6,0.94,0.94);
  lMUCL->SetFillColor(0);
  //c010->Divide(2,3);
  for(int i=0;i<6;i++){
    c010[i] = new TCanvas();
    //c010->cd(i+1);
    for(int a=0;a<3;a++){
      CLMuon_Distribution_Particle[i][a]->Sumw2();
      CLMuon_Distribution_Particle[i][a]->Scale(ScalingFactor);
      if(i==0) CLMuon_Distribution_Particle[i][a]->SetTitle("PM stopped");
      else if(i==2) CLMuon_Distribution_Particle[i][a]->SetTitle("PM escaping");
      else if(i==3) CLMuon_Distribution_Particle[i][a]->SetTitle("INGRID prematurely stopped");
      else if(i==4) CLMuon_Distribution_Particle[i][a]->SetTitle("INGRID stopped");
      else if(i==5) CLMuon_Distribution_Particle[i][a]->SetTitle("INGRID escaping");

      CLMuon_Distribution_Particle[i][a]->GetXaxis()->SetRangeUser(0,0.99);
      
      CLMuon_Distribution_Particle[i][a]->GetXaxis()->SetTitle("#mu_{CL}");
      CLMuon_Distribution_Particle[i][a]->GetYaxis()->SetTitle("Number of events");
      CLMuon_Distribution_Particle[i][a]->GetYaxis()->SetTitleOffset(1.3);
      CLMuon_Distribution_Particle[i][a]->Scale(ScalingFactor);
      CLMuon_Distribution_Particle[i][a]->SetLineWidth(2);
      if(a==0) CLMuon_Distribution_Particle[i][a]->SetLineColor(4);
      else if(a==1) CLMuon_Distribution_Particle[i][a]->SetLineColor(kGreen+2);
      else if(a==2) CLMuon_Distribution_Particle[i][a]->SetLineColor(2);
      if(a==0) CLMuon_Distribution_Particle[i][a]->Draw("E1HIST");
      else CLMuon_Distribution_Particle[i][a]->Draw("E1HISTsame");
    }
    if(i==0){
      if(CLMuon_Distribution_Particle[i][0]) lMUCL->AddEntry(CLMuon_Distribution_Particle[i][0],"#mu");
      if(CLMuon_Distribution_Particle[i][1]) lMUCL->AddEntry(CLMuon_Distribution_Particle[i][1],"#pi");
      if(CLMuon_Distribution_Particle[i][2]) lMUCL->AddEntry(CLMuon_Distribution_Particle[i][2],"p");
    }
    lMUCL->Draw("same");
    sprintf(Name,"MUCL_Particles_Topology%d.pdf",i);
    c010[i]->SaveAs(Name);
  }

  app->Run();
  return 0;
}

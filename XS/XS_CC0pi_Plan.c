
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include <TError.h>
#include <TF1.h>
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
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <THStack.h>
int LimitTracks=20;
int LimitRecs=10;
int NDials=175;
#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "INGRID_Dimension.cc"
INGRID_Dimension * IngDim = new INGRID_Dimension();
#include "setup.h"
#include "Hit.h"
#include "PMdispRev.h"
#include "Corrections.cc"
#include "Reconstruction.cc"
#include "Xsec.cc"
Reconstruction * Rec = new Reconstruction();
Corrections * Cor = new Corrections();
Xsec * XS = new Xsec();
#define LIKELIHOODHERE
//#define ALLSTUDY
//#define DEBUG_PID
//#define DEBUG
//#define XSEC_ERROR
//#define DEBUG2
//#define GENERATEWIDTH
//double C[17]={-1.85,0.34,0.59,0.74,0.514,-0.37,1.25,-0.06,0.562,0.82,-0.47,0.6,-0.57,-0.45};

TSpline3 * s_PMIng_Plan; TSpline3 * s_PMSci_Plan;TSpline3 * s_Ing_Plan;
TSpline3 * sCL_PMIng_Plan; TSpline3 * sCL_PMSci_Plan;TSpline3 * sCL_Ing_Plan;
TF1 * CL_PMIng_Plan;TF1 * CL_PMSci_Plan;TF1 * CL_Ing_Plan;

TSpline3 * s_PMIng_Likelihood_Muon; TSpline3 * s_PMSci_Likelihood_Muon;TSpline3 * s_Ing_Likelihood_Muon;
TSpline3 * s_PMIng_Likelihood_NotMuon; TSpline3 * s_PMSci_Likelihood_NotMuon;TSpline3 * s_Ing_Likelihood_NotMuon;
TF1 * CL_PMIng_Muon;TF1 * CL_PMSci_Muon;TF1 * CL_Ing_Muon;
TF1 * CL_PMIng_NotMuon;TF1 * CL_PMSci_NotMuon;TF1 * CL_Ing_NotMuon;

TH1D * PDFParticle;double PMuon;double P_NotMuon;


Double_t Likelihood_Ing_Plan(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_Ing_Plan->Eval(pe[0]);
}
Double_t Likelihood_PMIng_Plan(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_PMIng_Plan->Eval(pe[0]);
}
Double_t Likelihood_PMSci_Plan(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_PMSci_Plan->Eval(pe[0]);
}

void LoadMuCLDistributions_Plan(){
  /********************************Load MuCL************************************/
  TFile * file_Plan = new TFile("/home/bquilain/CC0pi_XS/XS/src/PDFMuCL_Plan.root");
  s_PMIng_Plan = (TSpline3*) file_Plan->Get("s_PMIng");
  s_PMSci_Plan = (TSpline3*) file_Plan->Get("s_PMSci");
  s_Ing_Plan = (TSpline3*) file_Plan->Get("s_Ing");

  sCL_PMIng_Plan = (TSpline3*) file_Plan->Get("sCL_PMIng");
  sCL_PMSci_Plan = (TSpline3*) file_Plan->Get("sCL_PMSci");
  sCL_Ing_Plan = (TSpline3*) file_Plan->Get("sCL_Ing");

  double Start_PMIng_Plan=s_PMIng_Plan->GetXmin();double End_PMIng_Plan=s_PMIng_Plan->GetXmax();
  double Start_PMSci_Plan=s_PMSci_Plan->GetXmin();double End_PMSci_Plan=s_PMSci_Plan->GetXmax();
  double Start_Ing_Plan=s_Ing_Plan->GetXmin();double End_Ing_Plan=s_Ing_Plan->GetXmax();

  CL_PMIng_Plan = new TF1("CL_PMIng_Plan",Likelihood_PMIng_Plan,Start_PMIng_Plan,End_PMIng_Plan*10,2);
  CL_PMSci_Plan = new TF1("CL_PMSci_Plan",Likelihood_PMSci_Plan,Start_PMSci_Plan,End_PMSci_Plan*10,2);
  CL_Ing_Plan = new TF1("CL_Ing_Plan",Likelihood_Ing_Plan,Start_Ing_Plan,End_Ing_Plan*10,2);
  
  CL_PMIng_Plan->SetParameter(0,Start_PMIng_Plan);CL_PMIng_Plan->SetParameter(1,100/*End_PMIng_Plan*/);
  CL_PMSci_Plan->SetParameter(0,Start_PMSci_Plan);CL_PMSci_Plan->SetParameter(1,100/*End_PMSci_Plan*/);
  CL_Ing_Plan->SetParameter(0,Start_Ing_Plan);CL_Ing_Plan->SetParameter(1,100/*End_Ing_Plan*/);

}




Double_t Likelihood_Ing_Muon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_Ing_Likelihood_Muon->Eval(pe[0]);
}
Double_t Likelihood_PMIng_Muon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_PMIng_Likelihood_Muon->Eval(pe[0]);
}
Double_t Likelihood_PMSci_Muon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_PMSci_Likelihood_Muon->Eval(pe[0]);
}

Double_t Likelihood_Ing_NotMuon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_Ing_Likelihood_NotMuon->Eval(pe[0]);
}
Double_t Likelihood_PMIng_NotMuon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_PMIng_Likelihood_NotMuon->Eval(pe[0]);
}
Double_t Likelihood_PMSci_NotMuon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_PMSci_Likelihood_NotMuon->Eval(pe[0]);
}



void LoadMuCLDistributions_Likelihood(){
  /********************************Load MuCL************************************/
  TFile * file_Likelihood = new TFile("/home/bquilain/CC0pi_XS/XS/src/PDFMuCL_Likelihood.root");
  PDFParticle = (TH1D*) file_Likelihood->Get("PDFParticle");
  PDFParticle->Scale(1./PDFParticle->Integral());
  double Start_PDG=PDFParticle->GetXaxis()->GetXmin();
  PMuon=PDFParticle->GetBinContent(-Start_PDG+13+1);
  //PMuon=PDFParticle->GetBinContent(1e4+13+1);
  P_NotMuon=PDFParticle->Integral()-PMuon;
  cout<<"initialisation, pmuon="<<PMuon<<endl;
  
  cout<<"CAREFUL: if you change PDFPARTICLE histogram in GeneratePDF_Likelihood in a non-symmetric binning around 0, you should change the code here to obtain Pmu. It assumes for now that the PDFParticle histogram has a symmetric binning around 0"<<endl;
  
  s_PMIng_Likelihood_Muon = (TSpline3*) file_Likelihood->Get("s_PMIng_Muon");
  s_PMSci_Likelihood_Muon = (TSpline3*) file_Likelihood->Get("s_PMSci_Muon");
  s_Ing_Likelihood_Muon = (TSpline3*) file_Likelihood->Get("s_Ing_Muon");

  s_PMIng_Likelihood_NotMuon = (TSpline3*) file_Likelihood->Get("s_PMIng_NotMuon");
  s_PMSci_Likelihood_NotMuon = (TSpline3*) file_Likelihood->Get("s_PMSci_NotMuon");
  s_Ing_Likelihood_NotMuon = (TSpline3*) file_Likelihood->Get("s_Ing_NotMuon");
  

  double Start_PMIng_Muon=s_PMIng_Likelihood_Muon->GetXmin();double End_PMIng_Muon=s_PMIng_Likelihood_Muon->GetXmax();
  double Start_PMSci_Muon=s_PMSci_Likelihood_Muon->GetXmin();double End_PMSci_Muon=s_PMSci_Likelihood_Muon->GetXmax();
  double Start_Ing_Muon=s_Ing_Likelihood_Muon->GetXmin();double End_Ing_Muon=s_Ing_Likelihood_Muon->GetXmax();

  
  CL_PMIng_Muon = new TF1("CL_PMIng_Muon",Likelihood_PMIng_Muon,Start_PMIng_Muon,End_PMIng_Muon,2);
  CL_PMSci_Muon = new TF1("CL_PMSci_Muon",Likelihood_PMSci_Muon,Start_PMSci_Muon,End_PMSci_Muon,2);
  CL_Ing_Muon = new TF1("CL_Ing_Muon",Likelihood_Ing_Muon,Start_Ing_Muon,End_Ing_Muon,2);
  
  CL_PMIng_Muon->SetParameter(0,Start_PMIng_Muon);CL_PMIng_Muon->SetParameter(1,End_PMIng_Muon);
  CL_PMSci_Muon->SetParameter(0,Start_PMSci_Muon);CL_PMSci_Muon->SetParameter(1,End_PMSci_Muon);
  CL_Ing_Muon->SetParameter(0,Start_Ing_Muon);CL_Ing_Muon->SetParameter(1,End_Ing_Muon);


  double Start_PMIng_NotMuon=s_PMIng_Likelihood_NotMuon->GetXmin();double End_PMIng_NotMuon=s_PMIng_Likelihood_NotMuon->GetXmax();
  double Start_PMSci_NotMuon=s_PMSci_Likelihood_NotMuon->GetXmin();double End_PMSci_NotMuon=s_PMSci_Likelihood_NotMuon->GetXmax();
  double Start_Ing_NotMuon=s_Ing_Likelihood_NotMuon->GetXmin();double End_Ing_NotMuon=s_Ing_Likelihood_NotMuon->GetXmax();

  CL_PMIng_NotMuon = new TF1("CL_PMIng_NotMuon",Likelihood_PMIng_NotMuon,Start_PMIng_NotMuon,End_PMIng_NotMuon,2);
  CL_PMSci_NotMuon = new TF1("CL_PMSci_NotMuon",Likelihood_PMSci_NotMuon,Start_PMSci_NotMuon,End_PMSci_NotMuon,2);
  CL_Ing_NotMuon = new TF1("CL_Ing_NotMuon",Likelihood_Ing_NotMuon,Start_Ing_NotMuon,End_Ing_NotMuon,2);
  
  CL_PMIng_NotMuon->SetParameter(0,Start_PMIng_NotMuon);CL_PMIng_NotMuon->SetParameter(1,End_PMIng_NotMuon);
  CL_PMSci_NotMuon->SetParameter(0,Start_PMSci_NotMuon);CL_PMSci_NotMuon->SetParameter(1,End_PMSci_NotMuon);
  CL_Ing_NotMuon->SetParameter(0,Start_Ing_NotMuon);CL_Ing_NotMuon->SetParameter(1,End_Ing_NotMuon);

#ifdef DEBUG_PID
  cout<<"The mucl functions are defined from "<<Start_PMIng_Muon<<"p.e to "<<End_PMIng_Muon<<"p.e. After, it is assumed that it is equal a constant mucl value (the same as the last bin)"<<endl;
  for(double pe=0;pe<500;pe+=2){
    cout<<"pe="<<pe<<", mucl="<<CL_PMIng_Muon->Eval(pe)<<endl;
  }
  cout<<endl;
#endif
     //////////////////////////////////////////////////////////////////////////////////////////////////////////
/*  
  TFile * Temp = new TFile("Temp.root","RECREATE");
  CL_PMIng_Likelihood->Write();
  CL_PMSci_Likelihood->Write();
  CL_Ing_Likelihood->Write();
  Temp->Close(); 
*/
}


bool IsINGRID(int ch){
  bool Ing;
  if(ch<=7||ch>=24) Ing=true;
  else Ing=false;
  return Ing;
}

int GetMax(vector <int> V){
  int Max;
  if(V[0]>V[1] && V[0]>V[2]) Max=0;
  else if(V[0]>V[1] && V[0]<V[2])Max=2;
  else if(V[0]<V[1] && V[1]>V[2])Max=1;
  else if(V[0]<V[1] && V[1]<V[2])Max=2;
  return Max;
}

double DegRad(double angle){
  return angle*TMath::Pi()/180.;
}

double RadDeg(double angle){
  return angle*180./TMath::Pi();
}
//


int main(int argc, char **argv)
{
#ifdef LIKELIHOODHERE
  TH1D * hTest_PMIng = new TH1D("hTest_PMIng","",600,0,300);
  TH1D * hTest_PMSci = new TH1D("hTest_PMSci","",600,0,300);
  TH1D * hTest_Ing = new TH1D("hTest_Ing","",600,0,300);
  hTest_PMIng->Sumw2();  hTest_PMSci->Sumw2();  hTest_Ing->Sumw2();

  TH1D * hTestImmediate_PMIng = new TH1D("hTestImmediate_PMIng","",600,0,300);
  TH1D * hTestImmediate_PMSci = new TH1D("hTestImmediate_PMSci","",600,0,300);
  TH1D * hTestImmediate_Ing = new TH1D("hTestImmediate_Ing","",600,0,300);
  hTestImmediate_PMIng->Sumw2();  hTestImmediate_PMSci->Sumw2();  hTestImmediate_Ing->Sumw2();

  TH2D * CLTest_PMIng = new TH2D("CLTest_PMIng","",600,0,300,100,0,1);
  TH2D * CLTest_PMSci = new TH2D("CLTest_PMSci","",600,0,300,100,0,1);
  TH2D * CLTest_Ing = new TH2D("CLTest_Ing","",600,0,300,100,0,1);
  CLTest_PMIng->Sumw2();  CLTest_PMSci->Sumw2();  CLTest_Ing->Sumw2();
#endif
  double LowCL=0;double HighCL=0;
  int nLowCL=0;int nHighCL=0;
  bool Disp=false;
  int NeutrinoType=0;

  int c=-1;
  bool MC=false;
  int RandomIteration=-1;
  //char * RandomIteration = new char[256];
  //char * InputFileName = new char[256];
  //char * OutputFileName = new char[256];
  string InputFileName;string OutputFileName;
  
  int ErrorType;string ErrorValue;char cErrorValue[256];
  
   
  TFile * fPEAngle;
  TH1D * PEAngleData_PMIng;  
  TH1D * PEAngleMC_PMIng;  
  TH1D * PEAngleData_PMSci;  
  TH1D * PEAngleMC_PMSci;  
  double SystematicsPECorrected_PMIng;
  double SystematicsPECorrected_PMSci;

  int XSEC=false;string xsec_file;
  
  while ((c = getopt(argc, argv, "i:o:f:dmr:x:e:v:")) != -1) {
    switch(c){
    case 'i':
      InputFileName=optarg;
      break;
    case 'o':
      OutputFileName=optarg;
      break;
    case 'f':
      NeutrinoType=atoi(optarg);
      break;
    case 'd':
      Disp=true;
      break;
    case 'm':
      MC=true;
      break;
    case 'x':
      XSEC=true;
      xsec_file=optarg;
      break;
    case 'e':
      ErrorType=atoi(optarg);
      break;
    case 'v':
      ErrorValue=optarg;
      break;
    }
  }


  if(ErrorType==4){
    fPEAngle = new TFile((ErrorValue).c_str());
    PEAngleData_PMIng = (TH1D*) fPEAngle->Get("PEAngleData_PMIng");  
    PEAngleMC_PMIng = (TH1D*) fPEAngle->Get("PEAngleMC_PMIng");  
    PEAngleData_PMSci = (TH1D*) fPEAngle->Get("PEAngleData_PMSci");  
    PEAngleMC_PMSci = (TH1D*) fPEAngle->Get("PEAngleMC_PMSci");      
  }  

  double Nu_E;
  double TrueParticleNRJ=0;
  IngridEventSummary* evt = new IngridEventSummary();
  TBranch * Br;
  IngridSimVertexSummary * simver;//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
  IngridSimParticleSummary * SimPart;
  BeamInfoSummary * BeamSummary;
  PMAnaSummary * recon;	
  TTree * tree;

  //////////////////////////////////////
  LoadMuCLDistributions_Plan();
  LoadMuCLDistributions_Likelihood();
  
  /*
  TCanvas * alo = new TCanvas();
  f_PMIng->SetNpx(500);
  f_PMIng->SaveAs("testf.root");
  cout<<f_PMIng->Eval(16.8)<<endl;
  */
  int FSIInt=-1;
  int Num_Int=-1;
  int nTracks=-1;
  double weight=1;
  bool IsFV=false;
  bool IsDetected=false;
  bool IsSand=false;
  bool IsAnti=false;
  bool IsNuE=false;
  bool IsBkgH=false;
  bool IsBkgV=false;
  double POT;
  double Enu;
  double TrueAngleMuon;
  double TrueMomentumMuon;
  double TrueAnglePion;//ML
  double TrueMomentumPion;//ML
  int NIngBasRec;
  int GoodSpill;
  int Spill;
  int VnTracks[LimitRecs];
  bool VIsDetected[LimitRecs];
  bool VSelectionOV[LimitRecs];
  bool VSelectionFV[LimitRecs];
  double OpeningAngle[LimitRecs];
  double CoplanarityAngle[LimitRecs];
  double TrackAngle[LimitRecs][LimitTracks];
  double TrackThetaX[LimitRecs][LimitTracks];
  double TrackThetaY[LimitRecs][LimitTracks];
  int TypeOfTrack[LimitRecs][LimitTracks];
  double CLMuon[LimitRecs][LimitTracks];
  double CLMuon_Plan[LimitRecs][LimitTracks];
  double CLMuon_KS[LimitRecs][LimitTracks];
  double CLMuon_Likelihood[LimitRecs][LimitTracks];
  double ProportionHighPE[LimitRecs][LimitTracks];
  double MeanHighPE[LimitRecs][LimitTracks];
  double HighestPE[LimitRecs][LimitTracks];
  double TotalCharge[LimitRecs][LimitTracks];
  int NHits_PMIng[LimitRecs][LimitTracks];
  int NHits_PMSci[LimitRecs][LimitTracks];
  int NHits_Ing[LimitRecs][LimitTracks];
  int LastChannelINGRIDX[LimitRecs][LimitTracks];
  int LastChannelINGRIDY[LimitRecs][LimitTracks];
  double TrackWidth[LimitRecs][LimitTracks];
  double Momentum[LimitRecs][LimitTracks];
  double ID[LimitRecs][LimitTracks];
  double PD[LimitRecs][LimitTracks];
  int Sample[LimitRecs][LimitTracks];
  double CriteriaAngleX[LimitRecs][LimitTracks];
  double CriteriaAngleY[LimitRecs][LimitTracks];
  double CriteriaHalfWayX[LimitRecs][LimitTracks];
  double CriteriaHalfWayY[LimitRecs][LimitTracks];

  double ReWeight[NDials];
  bool IsReconstructed[LimitRecs][LimitTracks];
  bool GT[LimitRecs][LimitTracks];


  TFile * wfile = new TFile((OutputFileName).c_str(),"recreate");
  wfile->cd();
  TTree*              wtree    = new TTree("wtree","wtree");
  wtree->SetDirectory(wfile);
  wtree              -> Branch   ("InteractionType",&Num_Int,"Num_Int/I");
  wtree              -> Branch   ("nIngBasRec",&NIngBasRec,"nIngBasRec/I");
  wtree              -> Branch   ("FSIInt",&FSIInt,"FSIInt/I");
  wtree              -> Branch   ("nTracks[10]",VnTracks,"nTracks[10]/I");
  wtree              -> Branch   ("weight",&weight,"weight/D");
  wtree              -> Branch   ("ReWeight[175]",&ReWeight,"ReWeight[175]/D");
  wtree              -> Branch   ("IsFV",&IsFV,"IsFV/O");
  wtree              -> Branch   ("IsSand",&IsSand,"IsSand/O");
  wtree              -> Branch   ("IsAnti",&IsAnti,"IsAnti/O");
  wtree              -> Branch   ("IsNuE",&IsNuE,"IsNuE/O");
  wtree              -> Branch   ("IsBkgH",&IsBkgH,"IsBkgH/O");
  wtree              -> Branch   ("IsBkgV",&IsBkgV,"IsBkgV/O"); 
  wtree              -> Branch   ("IsDetected[10]",&VIsDetected,"IsDetected[10]/O");
  wtree              -> Branch   ("SelectionFV[10]",&VSelectionFV,"SelectionFV[10]/O");
  wtree              -> Branch   ("SelectionOV[10]",&VSelectionOV,"SelectionOV[10]/O");
  wtree              -> Branch   ("POT",&POT,"POT/D");
  wtree              -> Branch   ("GoodSpill",&GoodSpill,"GoodSpill/I");
  wtree              -> Branch   ("Spill",&Spill,"Spill/I");
  wtree              -> Branch   ("Enu",&Enu,"Enu/D");
  wtree              -> Branch   ("TrueAngleMuon",&TrueAngleMuon,"TrueAngleMuon/D");
  wtree              -> Branch   ("TrueMomentumMuon",&TrueMomentumMuon,"TrueMomentumMuon/D");
  wtree              -> Branch   ("TrueAnglePion",&TrueAnglePion,"TrueAnglePion/D");//ML
  wtree              -> Branch   ("TrueMomentumPion",&TrueMomentumPion,"TrueMomentumPion/D");//ML
  wtree              -> Branch   ("TrackWidth[10][20]",&TrackWidth,"TrackWidth[10][20]/D");
  wtree              -> Branch   ("TrackAngle[10][20]",&TrackAngle,"TrackAngle[10][20]/D");
  wtree              -> Branch   ("TrackThetaX[10][20]",&TrackThetaX,"TrackThetaX[10][20]/D");
  wtree              -> Branch   ("TrackThetaY[10][20]",&TrackThetaY,"TrackThetaY[10][20]/D");
  wtree              -> Branch   ("TypeOfTrack[10][20]",&TypeOfTrack,"TypeOfTrack[10][20]/I");
  wtree              -> Branch   ("Sample[10][20]",&Sample,"Sample[10][20]/I");
  wtree              -> Branch   ("TotalCharge[10][20]",&TotalCharge,"TotalCharge[10][20]/D");
  wtree              -> Branch   ("OpeningAngle[10]",&OpeningAngle,"OpeningAngle[10]/D");
  wtree              -> Branch   ("CoplanarityAngle[10]",&CoplanarityAngle,"CoplanarityAngle[10]/D");
  wtree              -> Branch   ("CLMuon[10][20]",&CLMuon,"CLMuon[10][20]/D"); 
  wtree              -> Branch   ("CLMuon_Plan[10][20]",&CLMuon_Plan,"CLMuon_Plan[10][20]/D"); 
  wtree              -> Branch   ("CLMuon_KS[10][20]",&CLMuon_KS,"CLMuon_KS[10][20]/D");
  wtree              -> Branch   ("CLMuon_Likelihood[10][20]",&CLMuon_Likelihood,"CLMuon_Likelihood[10][20]/D");
  wtree              -> Branch   ("ProportionHighPE[10][20]",&ProportionHighPE,"ProportionHighPE[10][20]/D");
  wtree              -> Branch   ("MeanHighPE[10][20]",&MeanHighPE,"MeanHighPE[10][20]/D");
  wtree              -> Branch   ("HighestPE[10][20]",&HighestPE,"HighestPE[10][20]/D");
  wtree              -> Branch   ("LastChannelINGRIDX[10][20]",&LastChannelINGRIDX,"LastChannelINGRIDX[10][20]/I");                 
  wtree              -> Branch   ("LastChannelINGRIDY[10][20]",&LastChannelINGRIDY,"LastChannelINGRIDY[10][20]/I");                 
  wtree              -> Branch   ("Momentum[10][20]",&Momentum,"Momentum[10][20]/D");
  wtree              -> Branch   ("IronDistance[10][20]",&ID,"IronDistance[10][20]/D");
  wtree              -> Branch   ("PlasticDistance[10][20]",&PD,"PlasticDistance[10][20]/D");
  wtree              -> Branch   ("IsReconstructed[10][20]",&IsReconstructed,"IsReconstructed[10][20]/O");
  wtree              -> Branch   ("GeometricTrack[10][20]",&GT,"GeometricTrack[10][20]/O");                    
  wtree              -> Branch   ("CriteriaAngleX[10][20]",&CriteriaAngleX,"CriteriaAngleX[10][20]/D");
  wtree              -> Branch   ("CriteriaAngleY[10][20]",&CriteriaAngleY,"CriteriaAngleY[10][20]/D");
  wtree              -> Branch   ("CriteriaHalfWayX[10][20]",&CriteriaHalfWayX,"CriteriaHalfWayX[10][20]/D");
  wtree              -> Branch   ("CriteriaHalfWayY[10][20]",&CriteriaHalfWayY,"CriteriaHalfWayY[10][20]/D");


  ////////////////////////XSEC error case//////////////////////
  TFile * _file1;
  TTree * weightstree;
  int nevt_weightstree;
  TArrayF * reweight = NULL;
  TBranch * Br_reweight;
  /////////////////////////////////////////////////////////////


  TH1D * Flux = new TH1D("Flux","",2,0,2);
  TH1D * POTCount = new TH1D("POTCount","",2,0,2);
  double FluxTot=0;


  cout<<"Welcome"<<endl;
  vector <HitTemp> HitV;
  TApplication theApp("App",0,0);
  int TrackSample;

  cout<<"Opening Events"<<endl;
  cout<<"For now, we are not trying to deal with shared hits (2 tracks passing in the same scintillators)"<<endl;
  cout<<"CAREFUL: change the weight to take into account regions where there is also INGRID type scintillators"<<endl;
  cout<<"Careful: Mom3 in Simpart is always Kinetic Energy. I calculated formula that link it to momentum. TruePartNRJ in the code is in fact the momentum here"<<endl;

  ////////////////////////////////////////PREPARE THE INPUT TREE AND BRANCHES////////////////////////////////////////////////
  TFile * _file0 = new TFile((InputFileName).c_str(),"read");
  if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
  else{cout<<"Not able to read the file"<<endl; return 0;}
  tree=(TTree*) _file0->Get("tree");
  tree->SetDirectory(_file0);
  if(tree!=tree){cout<<"Problem in reading the tree, seems not to be any"<<endl;return 0;}
  int nevt=(int) tree->GetEntries();
  cout<<"Total Number Of Events="<<nevt<<endl;

  //////////////////////////////////XSEC case/////////////////////////////////////
  if(XSEC){  
    _file1 = new TFile(xsec_file.c_str());
    if(_file1->IsOpen()) cout << _file1->GetName() <<" is open"<< endl ;
    else{cout<<"No file containing the cross section error effects"<<endl; return 0;}

    weightstree=(TTree*) _file1->Get("weightstree");
    if(weightstree!=weightstree) return 0;
    nevt_weightstree=(int) weightstree->GetEntries();
    if(nevt_weightstree!=nevt) cout<<"Problem in weights tree: number of events different from Input reconstructed files. Are you sure to have input the correct weightstree?"<<endl;
    
    reweight = NULL;
    Br_reweight = weightstree->GetBranch("weights");
    Br_reweight->SetAddress(&reweight);
    weightstree->SetBranchAddress("weights",&reweight);
  }


  
    Br=tree->GetBranch("fDefaultReco.");
    Br->SetAddress(&evt);
    tree->SetBranchAddress("fDefaultReco.",&evt);

 
    ////////////////////////////////////////////////START THE LOOP//////////////////////////////////////////////////

    for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
    if((ievt%100)==0) cout<<"Processing "<<ievt<<endl;
      evt->Clear("C");
      tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
    
      for(int ibas=0;ibas<LimitRecs;ibas++){
	for(int itrk=0;itrk<LimitTracks;itrk++){
	  TrackAngle[ibas][itrk]=-1;
	  TrackThetaX[ibas][itrk]=180;//ML
	  TrackThetaY[ibas][itrk]=180;//ML
	  TrackWidth[ibas][itrk]=-1;
	  TypeOfTrack[ibas][itrk]=-1;
	  TotalCharge[ibas][itrk]=-1.;
	  CLMuon[ibas][itrk]=-1;
	  CLMuon_Plan[ibas][itrk]=-1;
	  CLMuon_KS[ibas][itrk]=-1;
	  CLMuon_Likelihood[ibas][itrk]=-1;
	  ProportionHighPE[ibas][itrk]=-1;
	  MeanHighPE[ibas][itrk]=-1;
	  HighestPE[ibas][itrk]=-1;

	  NHits_PMIng[ibas][itrk]=0;
	  NHits_PMSci[ibas][itrk]=0;
	  NHits_Ing[ibas][itrk]=0;
	  LastChannelINGRIDX[ibas][itrk]=-1;
	  LastChannelINGRIDY[ibas][itrk]=-1;
	  Momentum[ibas][itrk]=-1;
	  ID[ibas][itrk]=-1;
	  PD[ibas][itrk]=-1;
	  Sample[ibas][itrk]=-1;
	  IsReconstructed[ibas][itrk]=false;
	  GT[ibas][itrk]=false;
	  CriteriaHalfWayX[ibas][itrk]=-1;
	  CriteriaHalfWayY[ibas][itrk]=-1;
	  CriteriaAngleX[ibas][itrk]=-1;
	  CriteriaAngleY[ibas][itrk]=-1;
	}
	VIsDetected[ibas]=false;
	VSelectionFV[ibas]=false;
	VSelectionOV[ibas]=false;
	VnTracks[ibas]=-1;
	OpeningAngle[ibas]=-1.;
	CoplanarityAngle[ibas]=-1.;
      }
    
      NIngBasRec=-1;
      TrueAngleMuon=-1;
      TrueMomentumMuon=-1;
      TrueAnglePion=-1;
      TrueMomentumPion=-1;
      IsFV=false;
      IsDetected=false;
      FSIInt=-1;
      Num_Int=-1;
      nTracks=0;
      weight=0;
      Enu=0;
      GoodSpill=0;
      Spill=0;

      
      if(XSEC){
	if(reweight->GetSize()!=NDials && ievt==0) cout<<"Problem: change NDials="<<NDials<<" into "<<reweight->GetSize()<<" in the CC0pi code"<<endl;
	weightstree->GetEntry(ievt);
	if(ievt%100==0) cout<<endl<<"*************************************************************************************************"<<endl;
	for(int i=0;i<reweight->GetSize();i++){
	  ReWeight[i]=reweight->GetAt(i);
	  int BinError= (int) (i/7);
	  if(ievt%100==0){
	    if(BinError==22) cout<<" The Bin ";
	    cout<<ReWeight[i]<<", ";
	  }
	  if(ReWeight[i]!=ReWeight[i]) cout<<"********************************************************************************************************************************************************************************************************************************************"<<endl;
	}
	if(ievt%100==0) cout<<endl<<endl;
      }


     
      /////////////////////////////////////DEFINE THE WEIGHT AND TRUE VERTEX PROPERTIES AND TRUE MUON TRUE PROPERTIES//////////////////////////////
      simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
      int mod;
      double TrueVertexPosition[3]={0,0,0};
      double norm;
      double totcrsne;

      if(MC){
	Nu_E=simver->nuE;
	Enu=Nu_E;
	Num_Int=simver->inttype;
	mod=simver->mod;
	TrueVertexPosition[0]=simver->xnu;
	TrueVertexPosition[1]=simver->ynu;
	TrueVertexPosition[2]=simver->znu;
#ifdef DEBUG
	//cout<<"TrueVertexPosition[0]="<<TrueVertexPosition[0]<<", TrueVertexPosition[1]="<<TrueVertexPosition[1]<<", TrueVertexPosition[2]="<<TrueVertexPosition[2]+120<<endl;
#endif
	 norm=simver->norm;
	 totcrsne=simver->totcrsne;
	 FluxTot+=norm;
	 Flux->Fill(1,norm);
	 vector <double> MuonTrue = Rec->Reconstruction::GetTrueMuonInformation(evt);

	 TrueAngleMuon=MuonTrue[1];
	 TrueMomentumMuon=MuonTrue[0];
	 IsFV=Rec->Reconstruction::IsFV(mod,TrueVertexPosition[0],TrueVertexPosition[1],TrueVertexPosition[2]);

	 XS->Xsec::DetermineNuType(IsSand,IsAnti,IsNuE,IsBkgH,IsBkgV,simver->nutype,simver->mod);
	 FSIInt=XS->Xsec::DetermineFSI(IsSand,IsAnti,IsNuE,IsBkgH,IsBkgV,evt);

	 if(FSIInt==3){
	   vector<double> PionTrue = Rec->Reconstruction::GetTruePionInformation(evt);
	   TrueAnglePion=PionTrue[1];
	   TrueMomentumPion=PionTrue[0];
	 }

	 weight = 1;
	 if(IsBkgH==1 || IsBkgV) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*7.87*58.5;
	 else if(IsSand==1) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*470);
	 else weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*46.2;
      }
      else weight = Cor->GetMCCorrections(1,mod);

      POT=0;
      int utime=0;
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





      /////////////////////////////////////COUNT THE NUMBER OF POT/////////////////////////////////////////
      for(int ib=0;ib<evt->NIngridBeamSummarys();ib++){
	BeamSummary = evt->GetBeamSummary(ib);
	utime=BeamSummary->trg_sec;
	for( int cyc=Scyc; cyc<Ncyc; cyc++ ){
	  POTCount->Fill(1,BeamSummary->ct_np[4][cyc-4+1]);
	  POT+=BeamSummary->ct_np[4][cyc-4+1];
	}
	if(BeamSummary->good_spill_flag==1) GoodSpill=1;
	if(BeamSummary->spill_flag==1) Spill=1;
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





      //////////////////////////////////START THE LOOP ON THE RECONSTRUCTION (RECONSTRUCTED VERTEXES////////////////
      //cout<<NIngBasRec<<endl;
      NIngBasRec= evt->NPMAnas();
      for(int irec=0;irec<NIngBasRec;irec++){//loop on reconstruction: considered as same event => same time cluster, vertex reconstructed...
	recon = (PMAnaSummary*) evt->GetPMAna(irec);	
	nTracks=recon->Ntrack;


	///////////////////////////////SELECTIONS ARE DEFINED AND CHECKED///////////////////////////

	Rec->Reconstruction::GetSelectionPM(&VSelectionFV[irec],&VSelectionOV[irec],recon,MC);
	/////////////////////////////////////////////////////////////////////////////////////





	//////////////////////////////////////PROPERTIES OF ALL THE RECONSTRUCTION AND ALL TRACKS OF IT///////////////////////////
	vector <Hit3D> VecTrk;
	VecTrk.clear();
	VecTrk=Rec->Reconstruction::Hit2DMatchingAllTracksPM(recon,MC);//Contains all hits in the reconstruction
	double RecVertexPosition[3]={0,0,0};//For now, the vertex is defined for each track separately!!!

	vector< vector<Hit3D> > VecDouble;
	for(int i=0;i<VecDouble.size();i++) VecDouble[i].clear();
	for(int itrk=0;itrk<nTracks;itrk++){
	  vector <Hit3D> VecT;
	  VecT.clear();
	  HitV.clear();
	  HitV=Rec->Reconstruction::EraseDoubleHitsPM(recon,itrk,HitV);
	  VecT=Rec->Reconstruction::Hit2DMatchingPM(evt,recon,HitV,VecT,MC);
	  VecDouble.push_back(VecT);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	/////////////////////////////////////////// INITIALISE THE VARIABLES ATTACHED TO INFORMATION ON THE TRACK///////////////
	  if(Disp) cout<<endl<<"FSIInt="<<FSIInt<<", ntracks="<<nTracks<<endl;

	  vector <int> InteractionTracks;
	  InteractionTracks.clear();
	  vector <int> InteractionRecTracks;
          InteractionRecTracks.clear();
	  vector <int> InteractionSample;
          InteractionSample.clear();
	  vector <int> InteractionAngle;
          InteractionAngle.clear();

	  bool Geom=false;

	  vector <RecTrack> VRecTrack;
	  VRecTrack.clear();	

	  IsDetected=true;
	  VIsDetected[irec]=true;
	  VnTracks[irec]=nTracks;
	////////////////////////////////////////////////////////////////////////////////////////////	








	for(int itrk=0;itrk<nTracks;itrk++){//loop on track
	  RecVertexPosition[0]=(recon->x)[itrk]/10;
	  RecVertexPosition[1]=(recon->y)[itrk]/10;
	  RecVertexPosition[2]=VecTrk.front().z;

	  HitV.clear(); 
	  vector <Hit3D> Vec;
	  vector <Hit3D> Vec2;
	  vector <Hit3D> VecAll;
	  
	  /*************************Vec contains 3D Hits without any double counted******************************/
	  HitV=Rec->Reconstruction::EraseDoubleHitsPM(recon,itrk,HitV);
	  Vec=Rec->Reconstruction::Hit2DMatchingPM(evt,recon,HitV,Vec,MC);
	  Vec=Rec->CountSharedHits(Vec,VecDouble,itrk);//fill used
	
	  if(Vec.size()==0){
	    cout<<"Stopped because no Hits are corresponding"<<endl;
	    continue;
	  }
	   

	  /*******************************************Geometric properties of the tracks are estimated******************************/
	    double AngleX,AngleY,dx;
	    vector <double> GeomVariables;
	    GeomVariables.clear();
            GeomVariables=Rec->Reconstruction::GetTrackAngle(Vec);
            dx=GeomVariables[2];
	    dx=1./TMath::Cos(DegRad((recon->angle)[itrk]));



	    /****************************************Determine the Criteria for matching tracks*************************************/
	    vector <double> CriteriaPMINGRID=Rec->Reconstruction::GetMatchingPMINGRID(Vec);
	    CriteriaAngleX[irec][itrk]=TMath::Abs(CriteriaPMINGRID[4]-CriteriaPMINGRID[0]);
	    CriteriaAngleY[irec][itrk]=TMath::Abs(CriteriaPMINGRID[5]-CriteriaPMINGRID[1]);
	    CriteriaHalfWayX[irec][itrk]=TMath::Abs(CriteriaPMINGRID[6]-CriteriaPMINGRID[2]);
	    CriteriaHalfWayY[irec][itrk]=TMath::Abs(CriteriaPMINGRID[7]-CriteriaPMINGRID[3]);




	    if(ErrorType==4){
	      int BinAnglePEWidth_PMIng=PEAngleData_PMIng->GetBinWidth(1);
	      int BinAngle=((int) ((recon->angle)[itrk]/BinAnglePEWidth_PMIng))+1;
	      SystematicsPECorrected_PMIng=PEAngleData_PMIng->GetBinContent(BinAngle)-PEAngleMC_PMIng->GetBinContent(BinAngle);
	      if(PEAngleMC_PMIng->GetBinContent(BinAngle)!=0) SystematicsPECorrected_PMIng/=PEAngleMC_PMIng->GetBinContent(BinAngle);
	      
	      int BinAnglePEWidth_PMSci=PEAngleData_PMSci->GetBinWidth(1);
	      BinAngle=((int) ((recon->angle)[itrk]/BinAnglePEWidth_PMSci))+1;
	      SystematicsPECorrected_PMSci=PEAngleData_PMSci->GetBinContent(BinAngle)-PEAngleMC_PMSci->GetBinContent(BinAngle);
	      if(PEAngleMC_PMSci->GetBinContent(BinAngle)!=0) SystematicsPECorrected_PMSci/=PEAngleMC_PMSci->GetBinContent(BinAngle);
	     
#ifdef DEBUG
	      cout<<"Track angle="<<(recon->angle)[itrk]<<", Systematics Ing="<<SystematicsPECorrected_PMIng<<", Sci="<<SystematicsPECorrected_PMSci<<endl;
#endif
	    }
  


	    ///////////////////////////////////////// DETERMINE THE TRACK TOPOLOGY AND DISTANCES CROSSED IN IRON AND CARBON/////////////////
	    //0. IMPORTANT: CORRECT A BUG IN THE PLANE DETERMINATION
	    sort(Vec.begin(),Vec.end());
	    //Vec=Rec->Reconstruction::DetermineINGRIDPlanes(Vec);//TO DO! AND ALSO FOR STARTING PLANE!
	    //(recon->ing_endpln)[itrk]=Vec.back().pln;
	    
	    //1. Determine if the track should have reached INGRID if linearly extrapolated
	    Geom=Rec->Reconstruction::HasGeomTrack(16,(recon->startxpln)[itrk],(recon->startxch)[itrk],DegRad((recon->thetax)[itrk]), (recon->startypln)[itrk],(recon->startych)[itrk],DegRad((recon->thetay)[itrk]));
	    //2. Determine the track sample
	    TrackSample=Rec->Reconstruction::SelectTrackSample((recon->pm_stop)[itrk],Geom,(recon->ing_trk)[itrk],(recon->ing_stop)[itrk],(recon->ing_endpln)[itrk]);
	    //3. Search for INGRID hits aligned with tracks in the PM, in the case the latter is not matched w/ any INGRID track
	                  /********************************Here 3D hits are reconstructed, but not only those from the track, all those of the reconstruction************************/
	    VecAll=Rec->Reconstruction::Hit2DMatchingClusterPM(evt,recon);	   
            VecAll=Rec->Reconstruction::IsInTrk(VecAll,VecTrk);
	                  /*******************************They are used to determine if tracks stops in INGRID in the first iron plane (before being reconstructed)*************/
	    //Vec=Rec->Reconstruction::SearchIngridHit(Vec,VecAll,DegRad((recon->thetax)[itrk]), DegRad((recon->thetay)[itrk]),TrackSample);


	    /****************************************Attenuation from the dx in the scintillator is applied*************************/
            Vec=Cor->Corrections::GetFiberAttenuation(Vec);
	    sort(Vec.begin(),Vec.end());
	    
	    //4. Determine the distance crossed in Plastic and Iron by the track
	    vector <double> Dist;//
	    Dist=Rec->Reconstruction::TrackPenetrationPM((recon->startxpln)[itrk],(recon->startxch)[itrk], DegRad((recon->thetax)[itrk]), (recon->startypln)[itrk],(recon->startych)[itrk], DegRad((recon->thetay)[itrk]),(recon->endxpln)[itrk],(recon->endxch)[itrk],(recon->endypln)[itrk],(recon->endych)[itrk],(recon->ing_startmod)[itrk],(recon->ing_startpln)[itrk], (recon->ing_endpln)[itrk],dx,(recon->pm_stop)[itrk]);
	    
	    if((recon->ing_endpln)[itrk]>=9) Dist[1]=58.5/TMath::Cos(DegRad((recon->angle)[itrk]));
	    else Dist[1]=6.5*((recon->ing_endpln)[itrk]-(recon->ing_startpln)[itrk])/TMath::Cos(DegRad((recon->angle)[itrk]));
	    if((recon->ing_endpln)[itrk]<2 && recon->ing_trk[itrk]) cout<<"ONLY 2 PLANES"<<endl;
	    if(dx!=dx || dx==0) {cout<<"Problem in dx evaluation"<<endl;continue;}   
	    if((recon->ing_startmod)[itrk]!=(recon->ing_endmod)[itrk]){cout<<"Stopped because Ingrid module is changing"<<endl;continue;}

	    if(Disp) cout<<"Track Number="<<itrk<<endl;
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	    
	    ///////////////////////////////////Determine the particle true associated to the reconstructed track/////////////////////
	    double TrkLength = Rec->Reconstruction::GetTrackLength(Vec);
	    //cout<<"Track Length="<<TrkLength<<endl;
	    int Particle =-1;
	    TrueParticleNRJ=-1;
	    if(MC){
	      int SimPartNumber=Rec->Reconstruction::GetTrackParticle(evt, recon, itrk, TrkLength);
	      SimPart=(IngridSimParticleSummary*) evt->GetSimParticle(SimPartNumber);
	      Particle =SimPart->pdg;
	      TrueParticleNRJ=TMath::Sqrt(SimPart->momentum[0]*SimPart->momentum[0]+SimPart->momentum[1]*SimPart->momentum[1]+SimPart->momentum[2]*SimPart->momentum[2]);
	    }

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	   

	    /***********************************************FILL THE VARIABLES*********************************************/
	    double MuCL;double CL_Likelihood;
#ifdef LIKELIHOODHERE
	    
	    double PEPlane[3][NPlnPM][2];
	    double PEPlaneTot[3][NPlnPM][2];
	    int Plane[3][NPlnPM][2];
	    double PECorrected; 

	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int i=0;i<3;i++){
		for(int iview=0;iview<2;iview++){
		  PEPlane[i][ipln][iview]=0;
		  PEPlaneTot[i][ipln][iview]=0;
		  Plane[i][ipln][iview]=0;
		}
	      }
	    }
	    
	    /*
	    for(int i=0;i<Vec.size();i++){
	       if(Vec[i].used>1) continue;

	      if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)) Plane[1][Vec[i].pln][Vec[i].view]++;
	      else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))) Plane[0][Vec[i].pln][Vec[i].view]++;
	      else Plane[2][Vec[i].pln][Vec[i].view]++;
	    }
	    
	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int iview=0;iview<2;iview++){
		
		if(Plane[0][ipln][iview]!=0 || Plane[1][ipln][iview]!=0){//case PM & active plane
		  if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){
		    PEPlaneTot[0][ipln][iview]=f_PMSci_Plan->GetRandom(); 
		    //cout<<"PM, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlaneTot[0][ipln][iview]<<endl;
		    hTestImmediate_PMSci->Fill(PEPlaneTot[0][ipln][iview],weight);
		  }
		  else{
		    PEPlaneTot[1][ipln][iview]=f_PMIng_Plan->GetRandom();
		    //cout<<"PM, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlaneTot[1][ipln][iview]<<endl;
		    hTestImmediate_PMIng->Fill(PEPlaneTot[1][ipln][iview],weight);
		  }
		}
		
		else if(Plane[2][ipln][iview]!=0){
		  PEPlaneTot[2][ipln][iview]=f_Ing_Plan->GetRandom();
		  //cout<<"INGRID, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlaneTot[2][ipln][iview]<<endl;
		  hTestImmediate_Ing->Fill(PEPlaneTot[2][ipln][iview],weight);
		}
	      }
	    }
	    
	    for(int i=0;i<Vec.size();i++){
	      if(Vec[i].used>1) continue;
	      
	      if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
		PECorrected=(PEPlaneTot[0][Vec[i].pln][Vec[i].view]+PEPlaneTot[1][Vec[i].pln][Vec[i].view])/(Plane[0][Vec[i].pln][Vec[i].view]+Plane[1][Vec[i].pln][Vec[i].view]);
		Vec[i].pecorr=PECorrected*dx;
	      }
	      else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
		PECorrected=(PEPlaneTot[0][Vec[i].pln][Vec[i].view]+PEPlaneTot[1][Vec[i].pln][Vec[i].view])/(Plane[0][Vec[i].pln][Vec[i].view]+Plane[1][Vec[i].pln][Vec[i].view]);
		Vec[i].pecorr=PECorrected*dx*1.3*INGRIDSCIBAR;
	      }
		
	      else if(!(Vec[i].mod==16)){
		PECorrected=PEPlaneTot[2][Vec[i].pln][Vec[i].view]/Plane[2][Vec[i].pln][Vec[i].view];
		Vec[i].pecorr=PECorrected*dx;
	      }
	    }
	*/

	    //Idea is to remove first hits: check the first plane, whether x or y. Remove the first fourth planes
	    //cout<<"Number of hits in the reco="<<VecTrk.size()<<", in the track="<<Vec.size()<<endl;

	    //cout<<"x="<<Vec[0].x<<","<<RecVertexPosition[0]<<endl;
	    //cout<<"y="<<Vec[0].y<<","<<RecVertexPosition[1]<<endl;
	    //cout<<"first="<<Vec[0].pln<<", z="<<Vec[0].z<<","<<RecVertexPosition[2]<<endl;


	    
	    ////////////////////////////Procedure to fix a bug between vertex position of the track and firt hits -> to solve later
	    /*
	    double StartPlnX=0;
	    double StartPlnY=0;
	    for(int i=0;i<Vec.size();i++){
	      if(Vec[i].view==0){
		StartPlnX=Vec[i].pln;
		break;
	      }
	    }
	    for(int i=0;i<Vec.size();i++){
	      if(Vec[i].view==1){
		StartPlnY=Vec[i].pln;
		break;
	      }
	    }
	    

	      (recon->startxpln)[itrk]=StartPlnX;
	      (recon->startypln)[itrk]=StartPlnY;
	      ///////////////////////////////////////////////////////////
	      */
	    
	    
	    ////////////////////////////APPLY DX CORRECTION/////////////////////////////////////////////
	    Vec=Cor->Corrections::GetDXCorrection(Vec,dx);
	    sort(Vec.begin(),Vec.end());

	    
	    double CLPlan=1;
	    double CLLikelihood_Muon=1;double CLLikelihood_NotMuon=1;
	    
	    int NCLHits=0;//Number of hits in the vertex plane or the plane just after
	    for(int i=0;i<Vec.size();i++){

	      if(Vec[i].used>1) continue;
	      if(Vec[i].mod==16){
		if(Vec[i].view==0){
		  if((Vec[i].pln==(recon->startxpln)[itrk]) || (Vec[i].pln==((recon->startxpln)[itrk]+1)) ) continue;
		  else if(Vec[i].pln<(recon->startxpln)[itrk]) cout<<"Problem, the vertex of the track is downstream its first hit..."<<endl;
		}
		else{
		  if((Vec[i].pln==(recon->startypln)[itrk]) || (Vec[i].pln==((recon->startypln)[itrk]+1)) ) continue;
		  else if(Vec[i].pln<(recon->startypln)[itrk]) cout<<"Problem, the vertex of the track is downstream its first hit..."<<endl;
		}
	      }
	      NCLHits++;
	    }
	    //cout<<"startx="<<(recon->startxpln)[itrk]<<", starty="<<(recon->startypln)[itrk]<<", hits essentials to this track only="<<NCLHits<<endl;
#ifdef DEBUG2
	    cout<<"Number of hits="<<Vec.size()<<endl;
#endif
	    for(int i=0;i<Vec.size();i++){
	      
	      if(Vec[i].used>1) continue;
	      if(NCLHits<2){
		if(Vec[i].mod==16){
		  if(Vec[i].view==0){ if(Vec[i].pln==(recon->startxpln)[itrk]) continue;}
		  else{ if(Vec[i].pln==(recon->startypln)[itrk]) continue;}
		}
	      }
	      else{
		if(Vec[i].mod==16){
		  if(Vec[i].view==0){ if((Vec[i].pln==(recon->startxpln)[itrk]) || (Vec[i].pln==((recon->startxpln)[itrk]+1)) ) continue;}
		  else{ if((Vec[i].pln==(recon->startypln)[itrk]) || (Vec[i].pln==((recon->startypln)[itrk]+1)) ) continue;}
		}
	      }
#ifdef DEBUG2	      
	      cout<<"Module="<<Vec[i].mod<<", plane="<<Vec[i].pln<<", channel="<<Vec[i].ch<<", view="<<Vec[i].view<<", pe="<<Vec[i].pe<<", corr="<<Vec[i].pecorr<<endl;
#endif
	      
	      if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
		PECorrected=Vec[i].pecorr;
		if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_PMIng;
		PEPlane[1][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[1][Vec[i].pln][Vec[i].view]++;
	      }
	      else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
		PECorrected=Vec[i].pecorr;
		if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_PMSci;
		PEPlane[0][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[0][Vec[i].pln][Vec[i].view]++;
	      }
	      else{
		PECorrected=Vec[i].pecorr;
		if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_PMIng;
		PEPlane[2][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[2][Vec[i].pln][Vec[i].view]++;
	      }
	    }//End of loop over hits
	    //cout<<"Now compare"<<endl;
	    
	    	    
	    double cl;
	    double cllikelihood_muon;
	    double cllikelihood_notmuon;
	    int NHits=0;
	    int HighPECount=0;double MaxPE=0;double AverageHighPE=0;
	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int iview=0;iview<2;iview++){
		//if(NHits>2) break;
		if(PEPlane[0][ipln][iview]!=0 || PEPlane[1][ipln][iview]!=0){//case PM & active plane
		  if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){
		    //cout<<"PM, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<endl;
		    cl = sCL_PMSci_Plan->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    hTest_PMSci->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);
		    CLTest_PMSci->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],cl,weight);
		    CLPlan *=cl;
		    
		    cllikelihood_muon=CL_PMSci_Muon->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    cllikelihood_notmuon=CL_PMSci_NotMuon->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    CLLikelihood_Muon *= cllikelihood_muon;
		    CLLikelihood_NotMuon *= cllikelihood_notmuon;
#ifdef DEBUG_PID
		    cout<<"PM, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<", cl="<<cl<<", cllikelihoodmuon="<<cllikelihood_muon<<", not muon="<<cllikelihood_notmuon<<", number of hits SciBar="<<Plane[0][ipln][iview]<<", number of hits INGRID="<<Plane[1][ipln][iview]<<endl;
#endif
		  }
		  else{
		    //cout<<"PM, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<endl;
		    hTest_PMIng->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);
		    cl = sCL_PMIng_Plan->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    CLTest_PMIng->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],cl,weight);
		    CLPlan *=cl;
		    
		    cllikelihood_muon=CL_PMIng_Muon->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    cllikelihood_notmuon=CL_PMIng_NotMuon->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    CLLikelihood_Muon *= cllikelihood_muon;
		    CLLikelihood_NotMuon *= cllikelihood_notmuon;
#ifdef DEBUG_PID		    
		    cout<<"PM, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<", cl="<<cl<<", cllikelihoodmuon="<<cllikelihood_muon<<", not muon="<<cllikelihood_notmuon<<", number of hits SciBar="<<Plane[0][ipln][iview]<<", number of hits INGRID="<<Plane[1][ipln][iview]<<endl;
#endif
		    //CLLikelihood_Muon *= CL_PMIng_Muon->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		    //CLLikelihood_NotMuon *= CL_PMIng_NotMuon->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		  }
		  NHits++;
		  if(cl<0.1){
		    HighPECount++;
		    AverageHighPE+=PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview];
		    if(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]>MaxPE) MaxPE=PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview];
		  }
		}
		
		else if(Plane[2][ipln][iview]!=0){
		  //cout<<"INGRID, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[2][ipln][iview]<<endl;
		
		  hTest_Ing->Fill(PEPlane[2][ipln][iview],weight);
		  cl = sCL_Ing_Plan->Eval(PEPlane[2][ipln][iview]);
		  CLTest_Ing->Fill(PEPlane[2][ipln][iview],cl,weight);	
			    
		  CLPlan *=cl;
		  
		  cllikelihood_muon=CL_Ing_Muon->Eval(PEPlane[2][ipln][iview]);
		  cllikelihood_notmuon=CL_Ing_NotMuon->Eval(PEPlane[2][ipln][iview]);
		  CLLikelihood_Muon *= cllikelihood_muon;
		  CLLikelihood_NotMuon *= cllikelihood_notmuon;
#ifdef DEBUG_PID		  
		  cout<<"INGRID, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[2][ipln][iview]<<", cl="<<cl<<", cllikelihoodmuon="<<cllikelihood_muon<<", not muon="<<cllikelihood_notmuon<<", number of hits="<<Plane[2][ipln][iview]<<endl;
#endif
		  //CLLikelihood_Muon *= CL_Ing_Muon->Eval(PEPlane[2][ipln][iview]);
		  //CLLikelihood_NotMuon *= CL_Ing_NotMuon->Eval(PEPlane[2][ipln][iview]);
		  NHits++;
		  if(cl<0.1){
		    HighPECount++;
		    AverageHighPE+=PEPlane[2][ipln][iview];
		    if(PEPlane[2][ipln][iview]>MaxPE) MaxPE=PEPlane[2][ipln][iview];
		  }
		}
	      }
	      //cout<<CLLikelihood_Muon<<", "<<CLLikelihood_NotMuon<<endl;
	    }
	
	    double Sum=0;
	    for(int i=0;i<NHits;i++){
	      //for(int i=0;i<NHits[1];i++){
	      Sum+=(pow(-TMath::Log(CLPlan),i)/TMath::Factorial(i));
	    }
 
	   //cout<<"LogCL="<<-TMath::Log(CLPlan2)<<", Sum="<<Sum<<", NHits="<<NHits<<endl;
	    MuCL=CLPlan*Sum;
#ifdef DEBUG_PID
	    if(TrueMomentumMuon>1){
	      cout<<isinf(TMath::Log(MuCL))<<endl;
	      cout<<"Event number="<<ievt<<", Particle type="<<Particle<<", Distance="<<Dist[1]<<", MuCL="<<MuCL<<", MuCL Likelihood="<<CLLikelihood_Muon*PMuon/(CLLikelihood_Muon*PMuon+CLLikelihood_NotMuon*P_NotMuon)<<endl<<endl;
		//", Sum part="<<Sum<<", Product part="<<CLPlan<<", number of hits="<<NHits<<", number of entries in Vec="<<Vec.size()<<", number of hits registered="<<recon->NhitTs(itrk)<<endl;
	    }
#endif
	    ProportionHighPE[irec][itrk]=((double) HighPECount/NHits);
	    MeanHighPE[irec][itrk]=AverageHighPE;
	    if(HighPECount!=0) MeanHighPE[irec][itrk]/=HighPECount;
	    HighestPE[irec][itrk]=MaxPE;
	    //cout<<setprecision(3)<<"PMuon="<<PMuon<<", Like to be muon="<<CLLikelihood_Muon<<", not muon="<<CLLikelihood_NotMuon<<endl;
	    
	    CL_Likelihood=CLLikelihood_Muon*PMuon/(CLLikelihood_Muon*PMuon+CLLikelihood_NotMuon*P_NotMuon);
	    //cout<<"Proportions of Max PE="<<ProportionHighPE[irec][itrk]<<", Highest PE="<<HighestPE[irec][itrk]<<", Mean High PE="<<MeanHighPE[irec][itrk]<<endl;
	    //cout<<"MuCL plan="<<MuCL<<", pdg="<<Particle<<endl<<", Likelihood Mucl="<<CL_Likelihood<<endl;

	    
#endif

	    double CL=1;//XS->Xsec::GetMuCL(CL_Ing,CL_PMIng,CL_PMSci,Vec, dx, TrackSample,SystematicPE, RandomIteration);
	    double CL_Plan=MuCL;//XS->Xsec::GetMuCL_Plan(CL_Ing_Plan,CL_PMIng_Plan,CL_PMSci_Plan,Vec, dx, TrackSample,SystematicPE, RandomIteration);
	    if(Disp){
	      //if(TrackSample>=3 && Particle!=13){
	      //cout<<"Final selection ="<<Particle<<endl;
	      VecAll=Rec->Reconstruction::ClusterPM(evt,recon, nTracks);	
	      double Slpe[2]={GeomVariables[3],GeomVariables[7]};
	      double b[2]={GeomVariables[4],GeomVariables[8]};
	      double Zi[2]={GeomVariables[5],GeomVariables[9]};
	      double Zf[2]={GeomVariables[6],GeomVariables[10]};

	      EvtDisp(Vec);
	      //}
	    }


	    
	    vector <double> LastChan;
	LastChan=Rec->Reconstruction::GetLastINGRIDChannel(Vec,TrackSample);
	LastChannelINGRIDY[irec][itrk]=LastChan[0];
	LastChannelINGRIDX[irec][itrk]=LastChan[1];
	TrackAngle[irec][itrk]=(recon->angle)[itrk];
	TrackThetaX[irec][itrk]=(recon->thetax)[itrk];
	TrackThetaY[irec][itrk]=(recon->thetay)[itrk];
	TrackWidth[irec][itrk]=Rec->Reconstruction::GetINGRIDTrackWidth(Vec);
	GT[irec][itrk]=Geom;
	IsReconstructed[irec][itrk]=true;
	TypeOfTrack[irec][itrk]=Particle;
	CLMuon[irec][itrk]=CL;
	CLMuon_Plan[irec][itrk]=CL_Plan;
	CLMuon_Likelihood[irec][itrk]=CL_Likelihood;
	if(CLMuon_Plan[irec][itrk]<0.05){
	  LowCL+=weight;
	  nLowCL++;
	}
	else if(CLMuon_Plan[irec][itrk]>0.95){
	  HighCL+=weight;
	  nHighCL++;
	}
	CLMuon_KS[irec][itrk]=(recon->mucl)[itrk];
	Momentum[irec][itrk]=TrueParticleNRJ;
	PD[irec][itrk]=Dist[0];
	if(TrackSample>=2) ID[irec][itrk]=Dist[1];
	else ID[irec][itrk]=0;
	Sample[irec][itrk]=TrackSample;
	TotalCharge[irec][itrk]=Rec->Reconstruction::GetTrackEnergyPerDistance(Vec,dx);
	//cout<<"Energy of the track="<<Rec->Reconstruction::GetTrackEnergy(Vec)<<endl;

	
	//5. Determine the unit vector, opening angle and coplanarity
	if(VnTracks[irec]==2){
	  vector <double> Kinematic;
	  Kinematic.clear();
	  for(int itrk2=0;itrk2<VnTracks[irec];itrk2++){
	    if(itrk2==itrk) continue;
	    else{
	      Kinematic=Rec->Reconstruction::GetKinematic(DegRad((recon->angle)[itrk]), DegRad((recon->thetax)[itrk]), DegRad((recon->thetay)[itrk]), DegRad((recon->angle)[itrk2]), DegRad((recon->thetax)[itrk2]), DegRad((recon->thetay)[itrk2]));
	      
	    }
	  }
	  OpeningAngle[irec]=RadDeg(Kinematic[0]);
	  CoplanarityAngle[irec]=RadDeg(Kinematic[1]);
	  //cout<<"Opening="<<OpeningAngle[irec]<<", coplanarity="<<CoplanarityAngle[irec]<<endl;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
#ifdef DEBUG2
	if(Sample[irec][itrk]>=3){
	  cout<<"New Track, having momentum="<<Momentum[irec][itrk]<<endl;
	  cout<<"pm stop="<<(recon->pm_stop)[itrk]<<", has ingrid track="<<(recon->ing_trk)[itrk]<<", is stopped in ingrid="<<(recon->ing_stop)[itrk]<<", ingrid plane starts="<<(recon->ing_startpln)[itrk]<<", ingrid plane stop="<<(recon->ing_endpln)[itrk]<<", "<<endl; 
	  /*  for(int ihit=0;ihit<Vec.size();ihit++){
	    cout<<"Hit, module="<<Vec[ihit].mod<<", plane="<<Vec[ihit].pln<<", view="<<Vec[ihit].view<<", channel="<<Vec[ihit].ch<<endl;
	    }*/
	}
#endif
	//cout<<"FSIInt test="<<FSIInt<<", typeoftrack="<<TypeOfTrack[irec][itrk]<<endl;
	Vec.clear();
	Vec2.clear();
	VecAll.clear();
	}//Tracks
	   
      }//Recons
      wtree->Fill(); 
    }//Evt
    if(XSEC) _file1->Close();

    Br->Delete();
    cout<<"writing"<<endl;
    cout<<"Low CL="<<LowCL<<", High="<<HighCL<<endl;
    cout<<"Nb Low CL="<<nLowCL<<", Nb High="<<nHighCL<<endl;
    wfile->cd();
    wtree  -> Write();
#ifdef DEBUG
    f_PMIng_Plan->SetNpx(300);
    f_PMSci_Plan->SetNpx(300);
    f_Ing_Plan->SetNpx(300);
    f_PMIng_Plan->Write();
    f_PMSci_Plan->Write();
    f_Ing_Plan->Write();

    CL_PMIng_Plan->Write();
    CL_PMSci_Plan->Write();
    CL_Ing_Plan->Write();

    CLTest_PMIng->Write();
    CLTest_PMSci->Write();
    CLTest_Ing->Write();

    //hTest_PMIng->Scale(1./hTest_PMIng->GetMaximum());
    //hTest_PMSci->Scale(1./hTest_PMSci->GetMaximum());
    //hTest_Ing->Scale(1./hTest_Ing->GetMaximum());
 
    hTest_PMIng->Write();
    hTest_PMSci->Write();
    hTest_Ing->Write();

    hTestImmediate_PMIng->Write();
    hTestImmediate_PMSci->Write();
    hTestImmediate_Ing->Write();

#endif
    wfile  -> Write();
    wfile  -> Close();
    if(Disp) theApp.Run();    
    return(0);
}
  

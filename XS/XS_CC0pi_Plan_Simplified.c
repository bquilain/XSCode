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
#include <TGraphErrors.h>
#include <TSpline.h>
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
#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "setup.h"
#include "Hit.h"
#include "PMdispRev.h"
#include "Corrections.cc"
#include "Reconstruction.cc"
#include "Xsec.cc"

#define PID
//#define LIKELIHOODHERE

// you may choose ONE of the following refinements:
// (without any, you build a muCL as muon vs others)
//   - PI_LIKELIHOOD is to build a pion vs proton CL for non-mu-like particles
//   - MUPI_LIKELIHOOD is to build the muCL as muon+pion vs proton
//       -> ML 2017/08/08 the 2nd turns out to be the most efficient for WM
//       -> the 2nd and the baseline (w/o any) give the same results for PM

//#define PI_LIKELIHOOD
//#define MUPI_LIKELIHOOD
 

#define MVA
//#define MVACHECKING
//#define INTERPOLATION
//#define ALLSTUDY
//#define DEBUG_PID
//#define DEBUG_PID_BDT
//#define DEBUG
//#define XSEC_ERROR
//#define DEBUG2
//#define GENERATEWIDTH
//double C[17]={-1.85,0.34,0.59,0.74,0.514,-0.37,1.25,-0.06,0.562,0.82,-0.47,0.6,-0.57,-0.45};
const double positionerror=0.5;
const double energyerror=0.5;

TSpline3 * s_PMIng_Plan; TSpline3 * s_PMSci_Plan;TSpline3 * s_Ing_Plan;
TSpline3 * sCL_PMIng_Plan; TSpline3 * sCL_PMSci_Plan;TSpline3 * sCL_Ing_Plan;
TF1 * CL_PMIng_Plan;TF1 * CL_PMSci_Plan;TF1 * CL_Ing_Plan;

TSpline3 * s_PMIng_Likelihood_Muon; TSpline3 * s_PMSci_Likelihood_Muon;TSpline3 * s_Ing_Likelihood_Muon;
TSpline3 * s_PMIng_Likelihood_NotMuon; TSpline3 * s_PMSci_Likelihood_NotMuon;TSpline3 * s_Ing_Likelihood_NotMuon;
TSpline3 * s_WM_Likelihood_Muon[3], *s_WM_Likelihood_NotMuon[3];
TSpline3 * s_WM_Likelihood_Pion[3], *s_WM_Likelihood_NotPion[3];
TSpline3 * s_Ing_Likelihood_Pion, *s_Ing_Likelihood_NotPion;

TF1 * CL_PMIng_Muon;TF1 * CL_PMSci_Muon;TF1 * CL_Ing_Muon;
TF1 * CL_PMIng_NotMuon;TF1 * CL_PMSci_NotMuon;TF1 * CL_Ing_NotMuon;
TF1 * CL_WM_Muon, * CL_WM_NotMuon;
TF1 * CL_WM_Pion, *CL_WM_NotPion,* CL_Ing_Pion, *CL_Ing_NotPion;

TH1D * PDFParticle;double PMuon;double P_NotMuon;double PPion, P_NotPion;


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
  TFile * file_Plan = new TFile("$INSTALLREPOSITORY/XS/src/PDFMuCL_Plan.root");
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

Double_t Likelihood_WM_Muon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_WM_Likelihood_Muon[(int)pe[1]]->Eval(pe[0]);
}
Double_t Likelihood_WM_NotMuon(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_WM_Likelihood_NotMuon[(int)pe[1]]->Eval(pe[0]);
}
#ifdef PI_LIKELIHOOD
Double_t Likelihood_WM_Pion(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_WM_Likelihood_Pion[(int)pe[1]]->Eval(pe[0]);
}
Double_t Likelihood_WM_NotPion(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_WM_Likelihood_NotPion[(int)pe[1]]->Eval(pe[0]);
}
Double_t Likelihood_Ing_Pion(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_Ing_Likelihood_Pion->Eval(pe[0]);
}
Double_t Likelihood_Ing_NotPion(Double_t *pe, Double_t *par){
  if(pe[0]<par[0]) pe[0]=par[0];
  else if(pe[0]>par[1]) pe[0]=par[1];
  double Likelihood=s_Ing_Likelihood_NotPion->Eval(pe[0]);
}
#endif

void LoadMuCLDistributions_Likelihood(bool PM=true,int ErrorType=0, int ErrorValue=0, bool BirksRescaled=false){
  /********************************Load MuCL************************************/
  TFile * file_Likelihood;
#ifdef MUPI_LIKELIHOOD
  if(ErrorType==0) file_Likelihood = new TFile(Form("$INSTALLREPOSITORY/XS/src/PDFMuPiCL_Likelihood%s.root",(PM?"":"_WM")));
  else if(ErrorType==5 && BirksRescaled) file_Likelihood = new TFile(Form("$INSTALLREPOSITORY/XS/src/PDFMuPiCL_Likelihood_Birks%d%s.root",ErrorValue,(PM?"":"_WM")));
  else file_Likelihood = new TFile(Form("$INSTALLREPOSITORY/XS/src/PDFMuPiCL_Likelihood%s_Systematics%d_%d.root",(PM?"":"_WM"),ErrorType,ErrorValue));
#else
  file_Likelihood = new TFile(Form("$INSTALLREPOSITORY/XS/src/PDFMuCL_Likelihood%s.root",(PM?"":"_WM")));
#endif
  PDFParticle = (TH1D*) file_Likelihood->Get("PDFParticle");
  PDFParticle->Scale(1./PDFParticle->Integral());
  double Start_PDG=PDFParticle->GetXaxis()->GetXmin();
  PMuon=PDFParticle->GetBinContent(-Start_PDG+13+1);
  PPion=PDFParticle->GetBinContent(-Start_PDG-211+1)+PDFParticle->GetBinContent(-Start_PDG+211+1);
#ifdef MUPI_LIKELIHOOD
  PMuon+=PPion;
#endif
  P_NotMuon=PDFParticle->Integral()-PMuon;
  PPion/=P_NotMuon;
  P_NotPion=1-PPion;
#ifdef MUPI_LIKELIHOOD
  cout<<"initialisation, p(mu+pi)="<<PMuon<<" ,pnot(mu+pi)="<<P_NotMuon<<endl;
#else
  cout<<"initialisation, pmuon="<<PMuon<<" ,pnotmuon="<<P_NotMuon<<endl;
  cout<<"among the remaining, ppion="<<PPion<<" ,pnotpion="<<P_NotPion<<endl;
#endif
  cout<<"CAREFUL: if you change PDFPARTICLE histogram in GeneratePDF_Likelihood in a non-symmetric binning around 0, you should change the code here to obtain Pmu. It assumes for now that the PDFParticle histogram has a symmetric binning around 0"<<endl;

  if(PM) {
    s_PMIng_Likelihood_Muon = (TSpline3*) file_Likelihood->Get("s_PMIng_Muon");
    s_PMSci_Likelihood_Muon = (TSpline3*) file_Likelihood->Get("s_PMSci_Muon");
    s_PMIng_Likelihood_NotMuon = (TSpline3*) file_Likelihood->Get("s_PMIng_NotMuon");
    s_PMSci_Likelihood_NotMuon = (TSpline3*) file_Likelihood->Get("s_PMSci_NotMuon");
  }
  else{
    for (int i=0;i<3;i++){
      s_WM_Likelihood_Muon[i] = (TSpline3*) file_Likelihood->Get(Form("s_WM_Muon_%d",i));
      s_WM_Likelihood_NotMuon[i] = (TSpline3*) file_Likelihood->Get(Form("s_WM_NotMuon_%d",i));
#ifdef PI_LIKELIHOOD  
      s_WM_Likelihood_Pion[i] = (TSpline3*) file_Likelihood->Get(Form("s_WM_Pion_%d",i));
      s_WM_Likelihood_NotPion[i] = (TSpline3*) file_Likelihood->Get(Form("s_WM_NotPion_%d",i));
#endif
    }
  }
#ifdef PI_LIKELIHOOD  
  s_Ing_Likelihood_Pion = (TSpline3*) file_Likelihood->Get("s_Ing_Pion");
  s_Ing_Likelihood_NotPion = (TSpline3*) file_Likelihood->Get("s_Ing_NotPion");
#endif
  s_Ing_Likelihood_Muon = (TSpline3*) file_Likelihood->Get("s_Ing_Muon");
  s_Ing_Likelihood_NotMuon = (TSpline3*) file_Likelihood->Get("s_Ing_NotMuon");

  if(PM){
    double Start_PMIng_Muon=s_PMIng_Likelihood_Muon->GetXmin();double End_PMIng_Muon=s_PMIng_Likelihood_Muon->GetXmax();
    double Start_PMSci_Muon=s_PMSci_Likelihood_Muon->GetXmin();double End_PMSci_Muon=s_PMSci_Likelihood_Muon->GetXmax();
    CL_PMIng_Muon = new TF1("CL_PMIng_Muon",Likelihood_PMIng_Muon,Start_PMIng_Muon,End_PMIng_Muon,2);
    CL_PMSci_Muon = new TF1("CL_PMSci_Muon",Likelihood_PMSci_Muon,Start_PMSci_Muon,End_PMSci_Muon,2);
    CL_PMIng_Muon->SetParameter(0,Start_PMIng_Muon);CL_PMIng_Muon->SetParameter(1,End_PMIng_Muon);
    CL_PMSci_Muon->SetParameter(0,Start_PMSci_Muon);CL_PMSci_Muon->SetParameter(1,End_PMSci_Muon);

    double Start_PMIng_NotMuon=s_PMIng_Likelihood_NotMuon->GetXmin();double End_PMIng_NotMuon=s_PMIng_Likelihood_NotMuon->GetXmax();
    double Start_PMSci_NotMuon=s_PMSci_Likelihood_NotMuon->GetXmin();double End_PMSci_NotMuon=s_PMSci_Likelihood_NotMuon->GetXmax();
    CL_PMIng_NotMuon = new TF1("CL_PMIng_NotMuon",Likelihood_PMIng_NotMuon,Start_PMIng_NotMuon,End_PMIng_NotMuon,2);
    CL_PMSci_NotMuon = new TF1("CL_PMSci_NotMuon",Likelihood_PMSci_NotMuon,Start_PMSci_NotMuon,End_PMSci_NotMuon,2);
    CL_PMIng_NotMuon->SetParameter(0,Start_PMIng_NotMuon);CL_PMIng_NotMuon->SetParameter(1,End_PMIng_NotMuon);
    CL_PMSci_NotMuon->SetParameter(0,Start_PMSci_NotMuon);CL_PMSci_NotMuon->SetParameter(1,End_PMSci_NotMuon);
  }
  else{
    double Start_WM_Muon=s_WM_Likelihood_Muon[0]->GetXmin();double End_WM_Muon=s_WM_Likelihood_Muon[0]->GetXmax();
    double Start_WM_NotMuon=s_WM_Likelihood_NotMuon[0]->GetXmin();double End_WM_NotMuon=s_WM_Likelihood_NotMuon[0]->GetXmax();
    CL_WM_Muon = new TF1("CL_WM_Muon",Likelihood_WM_Muon,Start_WM_Muon,End_WM_Muon,2);
    CL_WM_NotMuon = new TF1("CL_WM_NotMuon",Likelihood_WM_NotMuon,Start_WM_NotMuon,End_WM_NotMuon,2);
    CL_WM_Muon->SetParameter(0,Start_WM_Muon);CL_WM_Muon->SetParameter(1,End_WM_Muon);
    CL_WM_NotMuon->SetParameter(0,Start_WM_NotMuon);CL_WM_NotMuon->SetParameter(1,End_WM_NotMuon);
#ifdef PI_LIKELIHOOD
    double Start_WM_Pion=s_WM_Likelihood_Pion[0]->GetXmin();double End_WM_Pion=s_WM_Likelihood_Pion[0]->GetXmax();
    double Start_WM_NotPion=s_WM_Likelihood_NotPion[0]->GetXmin();double End_WM_NotPion=s_WM_Likelihood_NotPion[0]->GetXmax();
    CL_WM_Pion = new TF1("CL_WM_Pion",Likelihood_WM_Pion,Start_WM_Pion,End_WM_Pion,2);
    CL_WM_NotPion = new TF1("CL_WM_NotPion",Likelihood_WM_NotPion,Start_WM_NotPion,End_WM_NotPion,2);
    CL_WM_Pion->SetParameter(0,Start_WM_Pion);CL_WM_Pion->SetParameter(1,End_WM_Pion);
    CL_WM_NotPion->SetParameter(0,Start_WM_NotPion);CL_WM_NotPion->SetParameter(1,End_WM_NotPion);
#endif 
  }
  double Start_Ing_Muon=s_Ing_Likelihood_Muon->GetXmin();double End_Ing_Muon=s_Ing_Likelihood_Muon->GetXmax();
  CL_Ing_Muon = new TF1("CL_Ing_Muon",Likelihood_Ing_Muon,Start_Ing_Muon,End_Ing_Muon,2);
  CL_Ing_Muon->SetParameter(0,Start_Ing_Muon);CL_Ing_Muon->SetParameter(1,End_Ing_Muon);
  double Start_Ing_NotMuon=s_Ing_Likelihood_NotMuon->GetXmin();double End_Ing_NotMuon=s_Ing_Likelihood_NotMuon->GetXmax();
  CL_Ing_NotMuon = new TF1("CL_Ing_NotMuon",Likelihood_Ing_NotMuon,Start_Ing_NotMuon,End_Ing_NotMuon,2);
  CL_Ing_NotMuon->SetParameter(0,Start_Ing_NotMuon);CL_Ing_NotMuon->SetParameter(1,End_Ing_NotMuon);
#ifdef PI_LIKELIHOOD
  double Start_Ing_Pion=s_Ing_Likelihood_Pion->GetXmin();double End_Ing_Pion=s_Ing_Likelihood_Pion->GetXmax();
  CL_Ing_Pion = new TF1("CL_Ing_Pion",Likelihood_Ing_Pion,Start_Ing_Pion,End_Ing_Pion,2);
  CL_Ing_Pion->SetParameter(0,Start_Ing_Pion);CL_Ing_Pion->SetParameter(1,End_Ing_Pion);
  double Start_Ing_NotPion=s_Ing_Likelihood_NotPion->GetXmin();double End_Ing_NotPion=s_Ing_Likelihood_NotPion->GetXmax();
  CL_Ing_NotPion = new TF1("CL_Ing_NotPion",Likelihood_Ing_NotPion,Start_Ing_NotPion,End_Ing_NotPion,2);
  CL_Ing_NotPion->SetParameter(0,Start_Ing_NotPion);CL_Ing_NotPion->SetParameter(1,End_Ing_NotPion);
#endif

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


int FSIInt=-1;
int Num_Int=-1;
int nTracks=-1;
float weight=1;
bool IsFV=false;
bool IsSand=false;
bool IsAnti=false;
bool IsNuE=false;
bool IsBkgH=false;
bool IsBkgV=false;
bool IsSciBkg=false;
float POT;
float Enu;
float TrueAngleMuon;
float TrueMomentumMuon;
float TrueAnglePion;
float TrueMomentumPion;
bool NewEvent;
int NIngBasRec;
int GoodSpill;
int Spill;
bool VIsDetected;
bool VSelectionOV;
bool VSelectionFV;
float OpeningAngle;
float CoplanarityAngle;

float TrackAngle[LimitTracks];
float TrackThetaY[LimitTracks];
float TrackThetaX[LimitTracks];
int TypeOfTrack[LimitTracks];
float CLMuon[LimitTracks];
float CLMuon_Plan[LimitTracks];
float CLMuon_KS[LimitTracks];
float CLMuon_Likelihood[LimitTracks];
float TotalCharge[LimitTracks];

float TrackWidth[LimitTracks];
float Momentum[LimitTracks];
float ID[LimitTracks];
float PD[LimitTracks];
int Sample[LimitTracks];
float CriteriaAngleX[LimitTracks];
float CriteriaAngleY[LimitTracks];
float CriteriaHalfWayX[LimitTracks];
float CriteriaHalfWayY[LimitTracks];
float EnergyDeposition[LimitTracks][LimitHits];
float EnergyDepositionSpline[LimitTracks][LimitHits];
int NViewsPerPlaneEnergyDeposition[LimitTracks][LimitHits];
int NViewsPerPlaneEnergyDepositionNonIsolated[LimitTracks][LimitHits];

float ReWeight[NDials];
bool IsReconstructed[LimitTracks];
bool GT[LimitTracks];

void ResetInputVariables(bool newEventReset=true){
  for(int itrk=0;itrk<LimitTracks;itrk++){
    TrackAngle[itrk]=-1;
    TrackThetaX[itrk]=-1;
    TrackThetaY[itrk]=-1;
    TrackWidth[itrk]=-1;
    TypeOfTrack[itrk]=-1;
    TotalCharge[itrk]=-1.;
    CLMuon[itrk]=-1;
    CLMuon_Plan[itrk]=-1;
    CLMuon_KS[itrk]=-1;
    CLMuon_Likelihood[itrk]=-1;
    
    Momentum[itrk]=-1;
    ID[itrk]=0;
    PD[itrk]=0;
    Sample[itrk]=-1;
    IsReconstructed[itrk]=false;
    GT[itrk]=false;
    CriteriaHalfWayX[itrk]=-1;
    CriteriaHalfWayY[itrk]=-1;
    CriteriaAngleX[itrk]=-1;
    CriteriaAngleY[itrk]=-1;
    for(int ihit=0;ihit<LimitHits;ihit++){
      EnergyDeposition[itrk][ihit]=0;
      EnergyDepositionSpline[itrk][ihit]=0;
      NViewsPerPlaneEnergyDeposition[itrk][ihit]=0;
      NViewsPerPlaneEnergyDepositionNonIsolated[itrk][ihit]=0;
    }
  }
  
  VIsDetected=false;
  VSelectionFV=false;
  VSelectionOV=false;
  OpeningAngle=-1.;
  CoplanarityAngle=-1.;
  NewEvent=false;
  nTracks=-1;

  if(newEventReset){
    NIngBasRec=-1;
    TrueAngleMuon=-1;
    TrueMomentumMuon=-1;
    TrueAnglePion=-1;
    TrueMomentumPion=-1;
    IsFV=false;
    FSIInt=-1;
    Num_Int=-1;
    weight=0;
    Enu=0;
    GoodSpill=0;
    Spill=0;
    IsSand=false;
    IsAnti=false;
    IsNuE=false;
    IsBkgH=false;
    IsBkgV=false;
    IsSciBkg=false;
  }

}


#ifdef MVA
float TrackAngleMVA;
int TypeOfTrackMVA;
float CLMuonMVA;
float CLMuon_PlanMVA;
float CLMuon_KSMVA;
float CLMuon_LikelihoodMVA;
float TotalChargeMVA;
float TrackWidthMVA;
float MomentumMVA;
float IDMVA;
float PDMVA;
int SampleMVA;
bool IsReconstructedMVA;
bool GTMVA;
float EnergyDepositionMVA[LimitHits];

void ResetInputVariablesMVA(){//Careful, only reset variables that are track dependent. It assumes that ResetInputVariables will reset the variables will reset the variables that are not track-dependent but event-dependent. In the version of 2017/02/15, it is the case. Please try to maintain this when you modify the code, or else, the event-dependent variables of wtreeMVA will never be reset!
  TrackAngleMVA=-1;
  TrackWidthMVA=-1;
  TypeOfTrackMVA=-1;
  TotalChargeMVA=-1.;
  CLMuonMVA=-1;
  CLMuon_PlanMVA=-1;
  CLMuon_KSMVA=-1;
  CLMuon_LikelihoodMVA=-1;  
  MomentumMVA=-1;
  IDMVA=0;
  PDMVA=0;
  SampleMVA=-1;
  IsReconstructedMVA=false;
  GTMVA=false;
  for(int ihit=0;ihit<LimitHits;ihit++){
    EnergyDepositionMVA[ihit]=0;
  }

}
#endif


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
  bool WM=false;
  int RandomIteration=-1;
  //char * RandomIteration = new char[256];
  //char * InputFileName = new char[256];
  //char * OutputFileName = new char[256];
  string InputFileName;string OutputFileName;
  
  int ErrorType=0;string ErrorValue;char cErrorValue[256];

  int nbad=0;//for checking  
   
  TFile * fPEAngle;
  TH1D * PEAngleData_Ing;  
  TH1D * PEAngleMC_Ing;  
  TH1D * PEAngleData_PMIng;  
  TH1D * PEAngleMC_PMIng;  
  TH1D * PEAngleData_PMSci;  
  TH1D * PEAngleMC_PMSci;  
  TH1D * PEAngleData_WM;  
  TH1D * PEAngleMC_WM;  
  double SystematicsPECorrected_Ing;
  double SystematicsPECorrected_WM;
  double SystematicsPECorrected_PMIng;
  double SystematicsPECorrected_PMSci;

  int XSEC=false;string xsec_file;
  int retuned=false;
  bool BirksRescaled=false;
  bool PIDSyst=false;

  while ((c = getopt(argc, argv, "i:o:f:dmrPx:e:v:wt:")) != -1) {
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
    case 'w':
      WM=true;
      break;
    case 'x':
      XSEC=true;
      xsec_file=optarg;
      break;
    case 't':
      retuned=true;
      xsec_file=optarg;
      break;
    case 'e':
      ErrorType=atoi(optarg);
      break;
    case 'v':
      ErrorValue=optarg;
      break;
    case 'r':
      BirksRescaled=true;
      break;
    case 'P':
      PIDSyst=true;
      break;
    }
  }

  XSEC=XSEC || retuned;

  bool PM=!WM;
  cout<<"Detector is "<<(PM?"PM":"WM")<<endl;

  Reconstruction * Rec = new Reconstruction(PM);
  Corrections * Cor = new Corrections(PM);
  Xsec * XS = new Xsec(PM);
  //XS->Initialize();
  InitializeGlobal(PM);

  if(ErrorType==4){
    fPEAngle = new TFile((ErrorValue).c_str());
    PEAngleData_Ing = (TH1D*) fPEAngle->Get("PEAngleData_Ing");  
    PEAngleMC_Ing = (TH1D*) fPEAngle->Get("PEAngleMC_Ing");  
    if(PM){
      PEAngleData_PMIng = (TH1D*) fPEAngle->Get("PEAngleData_PMIng");  
      PEAngleMC_PMIng = (TH1D*) fPEAngle->Get("PEAngleMC_PMIng");  
      PEAngleData_PMSci = (TH1D*) fPEAngle->Get("PEAngleData_PMSci");  
      PEAngleMC_PMSci = (TH1D*) fPEAngle->Get("PEAngleMC_PMSci");      
    }
    else {
      PEAngleData_WM = (TH1D*) fPEAngle->Get("PEAngleData_WM");  
      PEAngleMC_WM = (TH1D*) fPEAngle->Get("PEAngleMC_WM");      

      //ML tmp -- I use the INGRID MC/data comparisons from PM runs
      TFile* fPEAngle_tmp = new TFile("$INSTALLREPOSITORY/XS/files/PEXAngle_PM.root");
      PEAngleData_Ing = (TH1D*) fPEAngle_tmp->Get("PEAngleData_Ing");  
      PEAngleMC_Ing = (TH1D*) fPEAngle_tmp->Get("PEAngleMC_Ing");  
      //
    }
  }  

  int birksindex=0;
  double CBIRKS_values[3]={0.0185,0.0208,0.0231};
  if(ErrorType==5) birksindex=atoi(ErrorValue.c_str())*2;


  double Nu_E;
  double TrueParticleNRJ=0;
  IngridEventSummary* evt = new IngridEventSummary();
  IngridSimVertexSummary * simver;//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
  IngridSimParticleSummary * SimPart;
  BeamInfoSummary * BeamSummary;
  PMAnaSummary * recon;	
  TTree * tree;

  //////////////////////////////////////
#ifdef LIKELIHOODHERE
  LoadMuCLDistributions_Plan();
  if(PIDSyst){  
    if(ErrorType==5)  LoadMuCLDistributions_Likelihood(PM,5,birksindex,BirksRescaled);
    else if(ErrorType==4)  LoadMuCLDistributions_Likelihood(PM,4);
  }
  else  LoadMuCLDistributions_Likelihood(PM);
#endif

  TFile * wfile = new TFile((OutputFileName).c_str(),"recreate");
  wfile->cd();
  TTree*              wtree    = new TTree("wtree","wtree");
  wtree->SetDirectory(wfile);
  wtree              -> Branch   ("InteractionType",&Num_Int,"Num_Int/I");
  wtree              -> Branch   ("nIngBasRec",&NIngBasRec,"nIngBasRec/I");
  wtree              -> Branch   ("FSIInt",&FSIInt,"FSIInt/I");
  wtree              -> Branch   ("weight",&weight,"weight/F");
  wtree              -> Branch   (Form("ReWeight[%d]",NDials),ReWeight,Form("ReWeight[%d]/F",NDials));
  wtree              -> Branch   ("IsFV",&IsFV,"IsFV/O");
  wtree              -> Branch   ("IsSand",&IsSand,"IsSand/O");
  wtree              -> Branch   ("IsAnti",&IsAnti,"IsAnti/O");
  wtree              -> Branch   ("IsNuE",&IsNuE,"IsNuE/O");
  wtree              -> Branch   ("IsBkgH",&IsBkgH,"IsBkgH/O");
  wtree              -> Branch   ("IsBkgV",&IsBkgV,"IsBkgV/O");
  wtree              -> Branch   ("IsSciBkg",&IsSciBkg,"IsSciBkg/O"); 
  wtree              -> Branch   ("POT",&POT,"POT/F");
  wtree              -> Branch   ("GoodSpill",&GoodSpill,"GoodSpill/I");
  wtree              -> Branch   ("Spill",&Spill,"Spill/I");
  wtree              -> Branch   ("Enu",&Enu,"Enu/F");
  wtree              -> Branch   ("TrueAngleMuon",&TrueAngleMuon,"TrueAngleMuon/F");
  wtree              -> Branch   ("TrueMomentumMuon",&TrueMomentumMuon,"TrueMomentumMuon/F");
  wtree              -> Branch   ("TrueAnglePion",&TrueAnglePion,"TrueAnglePion/F");
  wtree              -> Branch   ("TrueMomentumPion",&TrueMomentumPion,"TrueMomentumPion/F");
  wtree              -> Branch   ("NewEvent",&NewEvent,"NewEvent/O");
  wtree              -> Branch   ("nTracks",&nTracks,"nTracks/I");
  wtree              -> Branch   ("IsDetected",&VIsDetected,"IsDetected/O");
  wtree              -> Branch   ("SelectionFV",&VSelectionFV,"SelectionFV/O");
  wtree              -> Branch   ("SelectionOV",&VSelectionOV,"SelectionOV/O");
  wtree              -> Branch   ("OpeningAngle",&(OpeningAngle),"OpeningAngle/F");
  wtree              -> Branch   ("CoplanarityAngle",&(CoplanarityAngle),"CoplanarityAngle/F");
  //wtree              -> Branch   (Form("EnergyDepositionDistance[%d][%d]",LimitTracks,LimitHits),EnergyDepositionDistance,Form("EnergyDepositionDistance[%d][%d]/F",LimitTracks,LimitHits));

  for(int itrk=0;itrk<LimitTracks;itrk++){
    wtree              -> Branch   (Form("TrackWidth_track%d",itrk),&(TrackWidth[itrk]),Form("TrackWidth_track%d/F",itrk));
    wtree              -> Branch   (Form("TrackAngle_track%d",itrk),&(TrackAngle[itrk]),Form("TrackAngle_track%d/F",itrk));
    wtree              -> Branch   (Form("TrackThetaY_track%d",itrk),&(TrackThetaY[itrk]),Form("TrackThetaY_track%d/F",itrk));
    wtree              -> Branch   (Form("TrackThetaX_track%d",itrk),&(TrackThetaX[itrk]),Form("TrackThetaX_track%d/F",itrk));
    wtree              -> Branch   (Form("TypeOfTrack_track%d",itrk),&(TypeOfTrack[itrk]),Form("TypeOfTrack_track%d/I",itrk));
    //wtree              -> Branch   (Form("TrackAngle_track%d",itrk),&(TrackAngle[itrk]),Form("TrackAngle_track%d/F",itrk));
    wtree              -> Branch   (Form("IsReconstructed_track%d",itrk),&(IsReconstructed[itrk]),Form("IsReconstructed_track%d/O",itrk));
    wtree              -> Branch   (Form("Sample_track%d",itrk),&(Sample[itrk]),Form("Sample_track%d/I",itrk));
    
    wtree              -> Branch   (Form("CLMuon_track%d",itrk),&(CLMuon[itrk]),Form("CLMuon_track%d/F",itrk));
    wtree              -> Branch   (Form("CLMuon_Plan_track%d",itrk),&(CLMuon_Plan[itrk]),Form("CLMuon_Plan_track%d/F",itrk));
    wtree              -> Branch   (Form("CLMuon_KS_track%d",itrk),&(CLMuon_KS[itrk]),Form("CLMuon_KS_track%d/F",itrk));
    wtree              -> Branch   (Form("CLMuon_Likelihood_track%d",itrk),&(CLMuon_Likelihood[itrk]),Form("CLMuon_Likelihood_track%d/F",itrk));
    wtree              -> Branch   (Form("IronDistance_track%d",itrk),&(ID[itrk]),Form("IronDistance_track%d/F",itrk));
    wtree              -> Branch   (Form("PlasticDistance_track%d",itrk),&(PD[itrk]),Form("PlasticDistance_track%d/F",itrk));
    wtree              -> Branch   (Form("GeometricTrack_track%d",itrk),&(GT[itrk]),Form("GeometricTrack_track%d/O",itrk));
    wtree              -> Branch   (Form("TotalCharge_track%d",itrk),&(TotalCharge[itrk]),Form("TotalCharge_track%d/F",itrk));
    wtree              -> Branch   (Form("Momentum_track%d",itrk),&(Momentum[itrk]),Form("Momentum_track%d/F",itrk));

    for(int ihit=0; ihit<LimitHits;ihit++){
      wtree              -> Branch   (Form("EnergyDeposition_track%d_hit%d",itrk,ihit),&(EnergyDeposition[itrk][ihit]),Form("EnergyDeposition_track%d_hit%d/F",itrk,ihit));
    }
  }

#ifdef MVA
  //One entry per track. When loop over track, should reset the track-depedent variables only, not the 
  
  TTree*              wtreeMVA    = new TTree("wtreeMVA","wtreeMVA");
  wtreeMVA->SetDirectory(wfile);
  wtreeMVA              -> Branch   ("InteractionType",&Num_Int,"Num_Int/I");
  wtreeMVA              -> Branch   ("nIngBasRec",&NIngBasRec,"nIngBasRec/I");
  wtreeMVA              -> Branch   ("FSIInt",&FSIInt,"FSIInt/I");
  wtreeMVA              -> Branch   ("weight",&weight,"weight/F");
  wtreeMVA              -> Branch   (Form("ReWeight[%d]",NDials),ReWeight,Form("ReWeight[%d]/F",NDials));
  wtreeMVA              -> Branch   ("IsFV",&IsFV,"IsFV/O");
  wtreeMVA              -> Branch   ("IsSand",&IsSand,"IsSand/O");
  wtreeMVA              -> Branch   ("IsAnti",&IsAnti,"IsAnti/O");
  wtreeMVA              -> Branch   ("IsNuE",&IsNuE,"IsNuE/O");
  wtreeMVA              -> Branch   ("IsBkgH",&IsBkgH,"IsBkgH/O");
  wtreeMVA              -> Branch   ("IsBkgV",&IsBkgV,"IsBkgV/O"); 
  wtreeMVA              -> Branch   ("IsSciBkg",&IsSciBkg,"IsSciBkg/O"); 
  wtreeMVA              -> Branch   ("POT",&POT,"POT/F");
  wtreeMVA              -> Branch   ("GoodSpill",&GoodSpill,"GoodSpill/I");
  wtreeMVA              -> Branch   ("Spill",&Spill,"Spill/I");
  wtreeMVA              -> Branch   ("Enu",&Enu,"Enu/F");
  wtreeMVA              -> Branch   ("TrueAngleMuon",&TrueAngleMuon,"TrueAngleMuon/F");
  wtreeMVA              -> Branch   ("TrueMomentumMuon",&TrueMomentumMuon,"TrueMomentumMuon/F");
  wtreeMVA              -> Branch   ("TrueAnglePion",&TrueAnglePion,"TrueAnglePion/F");
  wtreeMVA              -> Branch   ("TrueMomentumPion",&TrueMomentumPion,"TrueMomentumPion/F");
  wtreeMVA              -> Branch   ("NewEvent",&NewEvent,"NewEvent/O");
  wtreeMVA              -> Branch   ("nTracks",&nTracks,"nTracks/I");
  wtreeMVA              -> Branch   ("IsDetected",&VIsDetected,"IsDetected/O");
  wtreeMVA              -> Branch   ("SelectionFV",&VSelectionFV,"SelectionFV/O");
  wtreeMVA              -> Branch   ("SelectionOV",&VSelectionOV,"SelectionOV/O");
  wtreeMVA              -> Branch   ("OpeningAngle",&(OpeningAngle),"OpeningAngle/F");
  wtreeMVA              -> Branch   ("CoplanarityAngle",&(CoplanarityAngle),"CoplanarityAngle/F");
  wtreeMVA              -> Branch   (Form("TrackWidth"),&(TrackWidthMVA),Form("TrackWidth/F"));
  wtreeMVA              -> Branch   (Form("TrackAngle"),&(TrackAngleMVA),Form("TrackAngle/F"));
  wtreeMVA              -> Branch   (Form("TypeOfTrack"),&(TypeOfTrackMVA),Form("TypeOfTrack/I"));
  wtreeMVA              -> Branch   (Form("IsReconstructed"),&(IsReconstructedMVA),Form("IsReconstructed/O"));
  wtreeMVA              -> Branch   (Form("Sample"),&(SampleMVA),Form("Sample/I"));
  wtreeMVA              -> Branch   (Form("CLMuon"),&(CLMuonMVA),Form("CLMuon/F"));
  wtreeMVA              -> Branch   (Form("CLMuon_Plan"),&(CLMuon_PlanMVA),Form("CLMuon_Plan/F"));
  wtreeMVA              -> Branch   (Form("CLMuon_KS"),&(CLMuon_KSMVA),Form("CLMuon_KS/F"));
  wtreeMVA              -> Branch   (Form("CLMuon_Likelihood"),&(CLMuon_LikelihoodMVA),Form("CLMuon_Likelihood/F"));
  wtreeMVA              -> Branch   (Form("IronDistance"),&(IDMVA),Form("IronDistance/F"));
  wtreeMVA              -> Branch   (Form("PlasticDistance"),&(PDMVA),Form("PlasticDistance/F"));
  wtreeMVA              -> Branch   (Form("GeometricTrack"),&(GTMVA),Form("GeometricTrack/O"));
  wtreeMVA              -> Branch   (Form("TotalCharge"),&(TotalChargeMVA),Form("TotalCharge/F"));
  wtreeMVA              -> Branch   (Form("Momentum"),&(MomentumMVA),Form("Momentum/F"));


  for(int ihit=0; ihit<LimitHits;ihit++){
    wtreeMVA              -> Branch   (Form("EnergyDeposition_hit%d",ihit),&(EnergyDepositionMVA[ihit]),Form("EnergyDeposition_hit%d/F",ihit));
  }

#endif
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
  int TrackSample;
#ifdef DISP
  TApplication theApp("App",0,0);
#endif
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
  tree->SetBranchAddress("fDefaultReco.",&evt);

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
    weightstree->SetBranchAddress("weights",&reweight);
  }
 
  ////////////////////////////////////////////////START THE LOOP//////////////////////////////////////////////////
 
  for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
    if((ievt%100)==0) cout<<"Processing "<<ievt<<endl;
    evt->Clear("C");
    ResetInputVariables();
    tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
    NewEvent=true;
      
    if(XSEC){
      weightstree->GetEntry(ievt);
      if(reweight->GetSize()!=NDials && ievt==0) cout<<"Problem: change NDials="<<NDials<<" into "<<reweight->GetSize()<<" in the CC0pi code"<<endl;
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

      XS->Xsec::DetermineNuType(IsSand,IsAnti,IsNuE,IsBkgH,IsBkgV,IsSciBkg,simver->nutype,simver->mod,TrueVertexPosition);
      FSIInt=XS->Xsec::DetermineFSI(IsSand,IsAnti,IsNuE,IsBkgH,IsBkgV,IsSciBkg,evt);

      // cout<<IsSciBkg<<" "<<TrueVertexPosition[0]<<" "<<TrueVertexPosition[1]<<" "<<TrueVertexPosition[2]+120<<endl;

      if(FSIInt==3 || FSIInt==12){
	vector<double> PionTrue = Rec->Reconstruction::GetTruePionInformation(evt);
	TrueAnglePion=PionTrue[1];
	TrueMomentumPion=PionTrue[0];
      }

      weight = 1;
      if(IsBkgH || IsBkgV) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(7.87*58.5+22);// ML small contribution from scinti added
      else if(IsSand) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*300); // ML 2017/11/27 in fact the vertex is drawn only in a 3 meters shell (not 4.7m)
      else  if(PM) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*46.2;
      else weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*50.;//for WM - ML 2017/05/05 the water tank is 50cm deep
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
    NIngBasRec= evt->NPMAnas();
      
    if(NIngBasRec==0){//Case of no reconstruction!
      wtree->Fill();
      ResetInputVariables();
    }
      
    for(int irec=0;irec<NIngBasRec;irec++){//loop on reconstruction: considered as same event => same time cluster, vertex reconstructed...
      recon = (PMAnaSummary*) evt->GetPMAna(irec);	
      nTracks=recon->Ntrack;
      if(nTracks>LimitTracks) nTracks=LimitTracks;

      ///////////////////////////////SELECTIONS ARE DEFINED AND CHECKED///////////////////////////

      Rec->Reconstruction::GetSelectionPM(&VSelectionFV,&VSelectionOV,recon,MC);
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
      // VecDouble[itrk] is the vector of all hits of track itrk, sorted by increasing z
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
      VIsDetected=true;
      //cout<<"In IngBasRec, num of bas rec="<<irec<<", Is Detected="<<VIsDetected<<endl;
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
	Vec=Rec->CountSharedHits(Vec,VecDouble,itrk);//fill the variables used
	
	if(Vec.size()==0){
	  cout<<itrk<<" Stopped because no Hits are corresponding "<<recon->NhitTs(itrk)<<endl;
	  continue;
	}
	   

	/*******************************************Geometric properties of the tracks are estimated******************************/
	// following are not used (already computed in PMAna, LoliAna)
	double AngleX,AngleY,dx;
	vector <double> GeomVariables;
	/*	    GeomVariables.clear();
	  GeomVariables=Rec->Reconstruction::GetTrackAngle(Vec);
	  dx=GeomVariables[2];*/
	dx=1./TMath::Cos(DegRad((recon->angle)[itrk]));
	// cout<<"********* 3d angle="<<recon->angle[itrk]<<" ****************"<<endl;



	/**************************************Determine the Criteria for matching tracks*************************************/
	// cout<<Vec.size()<<" "<<recon->NhitTs(itrk)<<endl;
	//for(int i=0;i<Vec.size();i++) cout<<Vec[i].view<<endl;

	vector <double> CriteriaPMINGRID=Rec->Reconstruction::GetMatchingPMINGRID(Vec);
	CriteriaAngleX[itrk]=TMath::Abs(CriteriaPMINGRID[4]-CriteriaPMINGRID[0]);
	CriteriaAngleY[itrk]=TMath::Abs(CriteriaPMINGRID[5]-CriteriaPMINGRID[1]);
	CriteriaHalfWayX[itrk]=TMath::Abs(CriteriaPMINGRID[6]-CriteriaPMINGRID[2]);
	CriteriaHalfWayY[itrk]=TMath::Abs(CriteriaPMINGRID[7]-CriteriaPMINGRID[3]);




	if(ErrorType==4){
	  if(PM){
	    int BinAnglePEWidth_PMIng=PEAngleData_PMIng->GetBinWidth(1);
	    int BinAngle=((int) ((recon->angle)[itrk]/BinAnglePEWidth_PMIng))+1;
	    if(BinAngle<10) {SystematicsPECorrected_PMIng=PEAngleData_PMIng->GetBinContent(BinAngle)-PEAngleMC_PMIng->GetBinContent(BinAngle);
	      if(PEAngleMC_PMIng->GetBinContent(BinAngle)!=0 && PEAngleData_PMIng->GetBinContent(BinAngle)!=0)
		SystematicsPECorrected_PMIng/=PEAngleMC_PMIng->GetBinContent(BinAngle);
	      else SystematicsPECorrected_PMIng=0.;
	    }
	    else SystematicsPECorrected_PMIng=0.;
	      
	    int BinAnglePEWidth_PMSci=PEAngleData_PMSci->GetBinWidth(1);
	    BinAngle=((int) ((recon->angle)[itrk]/BinAnglePEWidth_PMSci))+1;
	    SystematicsPECorrected_PMSci=PEAngleData_PMSci->GetBinContent(BinAngle)-PEAngleMC_PMSci->GetBinContent(BinAngle);
	    if(PEAngleMC_PMSci->GetBinContent(BinAngle)!=0 && PEAngleData_PMSci->GetBinContent(BinAngle)!=0)
	      SystematicsPECorrected_PMSci/=PEAngleMC_PMSci->GetBinContent(BinAngle);
	    else SystematicsPECorrected_PMSci=0.;
	  }

	  // ML added 2017/07/21
	  int BinAnglePEWidth_Ing=PEAngleData_Ing->GetBinWidth(1);
	  int BinAngle=((int) ((recon->angle)[itrk]/BinAnglePEWidth_Ing))+1;
	  SystematicsPECorrected_Ing=PEAngleData_Ing->GetBinContent(BinAngle)-PEAngleMC_Ing->GetBinContent(BinAngle);
	  if(BinAngle<9){
	    if(PEAngleMC_Ing->GetBinContent(BinAngle)!=0 && PEAngleData_Ing->GetBinContent(BinAngle)!=0)
	      SystematicsPECorrected_Ing/=PEAngleMC_Ing->GetBinContent(BinAngle);
	    else SystematicsPECorrected_Ing=0.;
	  }
	  else  SystematicsPECorrected_Ing=0.;
	     
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
	Geom=Rec->Reconstruction::HasGeomTrack((PM?16:15),(recon->startxpln)[itrk],(recon->startxch)[itrk],DegRad((recon->thetax)[itrk]), (recon->startypln)[itrk],(recon->startych)[itrk],DegRad((recon->thetay)[itrk]));
	//2. Determine the track sample
	TrackSample=Rec->Reconstruction::SelectTrackSample((recon->pm_stop)[itrk],Geom,(recon->ing_trk)[itrk],(recon->ing_stop)[itrk],(recon->ing_endpln)[itrk]);
	//3. Search for INGRID hits aligned with tracks in the PM, in the case the latter is not matched w/ any INGRID track
	/********************************Here 3D hits are reconstructed, but not only those from the track, all those of the reconstruction************************/
	VecAll=Rec->Reconstruction::Hit2DMatchingClusterPM(evt,recon,MC);	   
	VecAll=Rec->Reconstruction::IsInTrk(VecAll,VecTrk);
	/*******************************They are used to determine if tracks stops in INGRID in the first iron plane (before being reconstructed)*************/
	//Vec=Rec->Reconstruction::SearchIngridHit(Vec,VecAll,DegRad((recon->thetax)[itrk]), DegRad((recon->thetay)[itrk]),TrackSample);
	// ML 20170127 since it is also done is PMAna/LoliAna, I don't check it

	/****************************************Attenuation from the dx in the scintillator is applied*************************/
	Vec=Cor->Corrections::GetFiberAttenuation(Vec);
	sort(Vec.begin(),Vec.end());
      
	//4. Determine the distance crossed in Plastic and Iron by the track
	vector <double> Dist;//
	if(dx!=dx || dx==0) {cout<<"Problem in dx evaluation"<<endl;continue;}   
	Dist=Rec->Reconstruction::TrackPenetrationPM((recon->startxpln)[itrk],(recon->startxch)[itrk], DegRad((recon->thetax)[itrk]), (recon->startypln)[itrk],(recon->startych)[itrk], DegRad((recon->thetay)[itrk]),(recon->endxpln)[itrk],(recon->endxch)[itrk],(recon->endypln)[itrk],(recon->endych)[itrk],(recon->ing_startmod)[itrk],(recon->ing_startpln)[itrk], (recon->ing_endpln)[itrk], dx, TrackSample, Vec);
	         
	//if((recon->ing_endpln)[itrk]>=9) Dist[1]=58.5/TMath::Cos(DegRad((recon->angle)[itrk]));
	//else Dist[1]=6.5*((recon->ing_endpln)[itrk]-(recon->ing_startpln)[itrk])/TMath::Cos(DegRad((recon->angle)[itrk]));
	//if((recon->ing_endpln)[itrk]<2 && recon->ing_trk[itrk]) cout<<"ONLY 2 PLANES"<<endl;
	if(dx!=dx || dx==0) {cout<<"Problem in dx evaluation"<<endl;continue;}   
	//if((recon->ing_startmod)[itrk]!=(recon->ing_endmod)[itrk]){cout<<"Stopped because Ingrid module is changing"<<endl;continue;}

	if(Disp) cout<<"Track Number="<<itrk<<endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	    
	///////////////////////////////////Determine the particle true associated to the reconstructed track/////////////////////
	double TrkLength = Rec->Reconstruction::GetTrackLength(Vec);
	//cout<<"Track Length="<<TrkLength<<endl;
	int Particle =-1;
	TrueParticleNRJ=-1;
	if(MC){
	  //cout<<"Number of Sim particle in this event="<<evt->NIngridSimParticles()<<", type of the first ="<<evt->GetSimParticle(0)->pdg<<endl;;
	  //cout<<"pm stop="<<(recon->pm_stop)[itrk]<<", has ingrid track="<<(recon->ing_trk)[itrk]<<", is stopped in ingrid="<<(recon->ing_stop)[itrk]<<", ingrid plane starts="<<(recon->ing_startpln)[itrk]<<", ingrid plane stop="<<(recon->ing_endpln)[itrk]<<", "<<endl; 
	  int SimPartNumber=Rec->Reconstruction::GetTrackParticleBugged(evt, recon, itrk, TrkLength);
	  //SimPart=(IngridSimParticleSummary*) evt->GetSimParticle(SimPartNumberBugged);
	  //int ParticleBugged =SimPart->pdg;
	  //int SimPartNumber=Rec->Reconstruction::GetTrackParticle(evt, recon, itrk, TrkLength);
	  SimPart=(IngridSimParticleSummary*) evt->GetSimParticle(SimPartNumber);
	  Particle =SimPart->pdg;
	  //cout<<"Old bugged pdg = "<<ParticleBugged<<", new = "<<Particle<<endl<<endl;
	  //	  if(Particle==2112) cout<<"event "<<ievt<<" reco "<<irec<<" track="<<itrk<<" is associated to neutron"<<endl; 
	  TrueParticleNRJ=TMath::Sqrt(SimPart->momentum[0]*SimPart->momentum[0]+SimPart->momentum[1]*SimPart->momentum[1]+SimPart->momentum[2]*SimPart->momentum[2]);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////APPLY DX CORRECTION/////////////////////////////////////////////
	Vec=Cor->Corrections::GetDXCorrection(Vec,dx);
	// ML tmp
	Cor->GetDXCorrectionWM(Vec,DegRad(recon->angle[itrk]),DegRad(recon->thetax[itrk]),DegRad(recon->thetay[itrk]));
	//
	sort(Vec.begin(),Vec.end());
	   

	/***********************************************FILL THE VARIABLES*********************************************/
#ifdef PID
	double MuCL;double CL_Likelihood, CL_Likelihood_Pi;
      
	double PEPlane[3][NPlnPM][2];
	double PEPlaneTot[3][NPlnPM][2];
	int Plane[3][NPlnPM][2];
	int PlaneNonIsolated[3][NPlnPM][2];
	double PECorrected;
	double DistanceBarycenter[3][NPlnPM][2];

	for(int ipln=0;ipln<NPlnPM;ipln++){
	  for(int i=0;i<3;i++){
	    for(int iview=0;iview<2;iview++){
	      PEPlane[i][ipln][iview]=0;
	      PEPlaneTot[i][ipln][iview]=0;
	      Plane[i][ipln][iview]=0;
	      PlaneNonIsolated[i][ipln][iview]=0;
	      DistanceBarycenter[i][ipln][iview]=0;
	    }
	  }
	}

	if(PM){
	  double CLPlan=1;
	  double CLLikelihood_Muon=1;double CLLikelihood_NotMuon=1;

	  //1. Count the number of hits that are "useable" for the CL (not too close from the vertex). The goal (later used) is that for events where the number of hits useable < 2, we will use hits near the vertex
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


	  //2. Real loop on the event:
#ifdef DEBUG2
	  cout<<"Number of hits="<<Vec.size()<<endl;
#endif
	  sort(Vec.begin(),Vec.end());
	  double ControlTrkLength = TMath::Sqrt(pow(Vec.back().x-Vec.front().x,2)+pow(Vec.back().y-Vec.front().y,2)+pow(Vec.back().z-Vec.front().z,2));

	  for(int i=0;i<Vec.size();i++){

	    //2.a Remove hits near the vertex
	    if(NCLHits<=2){
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

	    //2. b Estimate the distance of the hit from the track vertex
	    double Distance = TMath::Sqrt(pow(Vec[i].x-Vec.front().x,2)+pow(Vec[i].y-Vec.front().y,2)+pow(Vec[i].z-Vec.front().z,2));
	      
#ifdef DEBUG2
	    //if(DistanceBarycenter/Dist[2] > 1.) cout<<endl<<endl<<endl<<endl<<"CAREFUL: Track is too long, Track Length="<<Dist[2]<<endl<<endl<<endl;
	    cout<<"Module="<<Vec[i].mod<<", plane="<<Vec[i].pln<<", channel="<<Vec[i].ch<<", view="<<Vec[i].view<<", pe="<<Vec[i].pe<<", corr="<<Vec[i].pecorr<<endl;
#endif

	    //2. c Count the number of hits per plane/view -> Check if we have neighbour hits on the same plane: used later for measuring transverse track width + Barycenter of the hit in one plane.
	    //All hits
	    if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)) PlaneNonIsolated[1][Vec[i].pln][Vec[i].view]++;
	    else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))) PlaneNonIsolated[0][Vec[i].pln][Vec[i].view]++;
	    else PlaneNonIsolated[2][Vec[i].pln][Vec[i].view]++;


	    double dedz,MEV2PE=46.,MEV2PE_PM=40.;
	    double BirksCorr=1.;
	    //Isolated hits

	    //2.d Remove hits used in 2 tracks if we have more than 2 hits independent from other tracks
	    if(NCLHits>2 && Vec[i].used>1) continue;
	    //	    cout<<recon->angle[itrk]<<" "<<Vec[i].pecorr<<" ";

	    //2.e Specific to PM: Classify the hits between INGRID PM, SciBar, and INGRID. Store the total PE per plane/view for each type
	    if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
	      PECorrected=Vec[i].pecorr;
	      dedz=Vec[i].pecorr/MEV2PE_PM/(1.+.09)/.275;// already dz-corrected
	      if(BirksRescaled) {
		BirksCorr=1-dedz*(CBIRKS_values[birksindex]-CBIRKS_values[1])/(1.+dedz*CBIRKS_values[birksindex]);
		PECorrected*=BirksCorr;
		cout<<"de/dz="<<dedz<<" (Mev/cm), Delta CB="<<CBIRKS_values[birksindex]-CBIRKS_values[1]<<" BirksCorr="<<BirksCorr<<endl;
	      }
	      if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_PMIng;
	      PEPlane[1][Vec[i].pln][Vec[i].view]+=PECorrected;
	      Plane[1][Vec[i].pln][Vec[i].view]++;
#ifdef MVA
	      DistanceBarycenter[1][Vec[i].pln][Vec[i].view]+=Vec[i].dist_plastic/IronCarbonRatio+Vec[i].dist_iron;
#endif
	    }
	    else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
	      PECorrected=Vec[i].pecorr;
	      dedz=Vec[i].pecorr/MEV2PE_PM/(1.+.09)/.275/1.88;
	      if(BirksRescaled) {
		BirksCorr=1-dedz*(CBIRKS_values[birksindex]-CBIRKS_values[1])/(1.+dedz*CBIRKS_values[birksindex]);
		PECorrected*=BirksCorr;
	      }
	      if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_PMSci;
	      PEPlane[0][Vec[i].pln][Vec[i].view]+=PECorrected;
	      Plane[0][Vec[i].pln][Vec[i].view]++;
#ifdef MVA
	      DistanceBarycenter[0][Vec[i].pln][Vec[i].view]+=Vec[i].dist_plastic/IronCarbonRatio+Vec[i].dist_iron;
#endif
	    }
	    else{
	      PECorrected=Vec[i].pecorr;
	      dedz=Vec[i].pecorr/MEV2PE/(1.+.09)/.275;
	      if(BirksRescaled) {
		BirksCorr=1-dedz*(CBIRKS_values[birksindex]-CBIRKS_values[1])/(1.+dedz*CBIRKS_values[birksindex]);
		PECorrected*=BirksCorr;
	      }
	      if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_Ing;
	      PEPlane[2][Vec[i].pln][Vec[i].view]+=PECorrected;
	      Plane[2][Vec[i].pln][Vec[i].view]++;
#ifdef MVA
	      DistanceBarycenter[2][Vec[i].pln][Vec[i].view]+=Vec[i].dist_plastic/IronCarbonRatio+Vec[i].dist_iron;
#endif
	    }
	    //cout<<PECorrected<<endl;
	  }//End of loop over hits
	  //cout<<"Now compare"<<endl;
	    

	  double cl;
	  double cllikelihood_muon;
	  double cllikelihood_notmuon;
	  int HighPECount=0;double MaxPE=0;double AverageHighPE=0;
	  int NHits=0;


	  //3. Estimate each plane barycenter for the MVA. It allows to store: position of the hit vs total charge in the plane
	  for(int ipln=0;ipln<NPlnPM;ipln++){
	    for(int iview=0;iview<2;iview++){
	      //if(NHits>2) break;
	      //cout<<"Pln="<<ipln<<", view="<<iview<<", PM is active="<<(PEPlane[0][ipln][iview]!=0 || PEPlane[1][ipln][iview]!=0)<<", INGRID active="<<(PEPlane[2][ipln][iview]!=0)<<endl;
#ifdef MVA	
	      //All hits
	      if(PlaneNonIsolated[0][ipln][iview]!=0 || PlaneNonIsolated[1][ipln][iview]!=0){//case PM & active plane
		double Barycenter = (DistanceBarycenter[0][ipln][iview]+DistanceBarycenter[1][ipln][iview]) / (PlaneNonIsolated[0][ipln][iview]+PlaneNonIsolated[1][ipln][iview]);
		double RelativeBarycenter = Barycenter / Dist[2];
		double bindistancesize = (1. / (LimitHits-1));//number of separation between intervals = number of intervals -1, as usual...
		int bindistance = ((int) (RelativeBarycenter / bindistancesize));
		if(bindistance >= LimitHits) bindistance=LimitHits-1;
		NViewsPerPlaneEnergyDepositionNonIsolated[itrk][bindistance]++;
	      }
	      if(PlaneNonIsolated[2][ipln][iview]!=0){
		double Barycenter = (DistanceBarycenter[2][ipln][iview]) / (PlaneNonIsolated[2][ipln][iview]);
		double RelativeBarycenter = Barycenter / Dist[2];
		double bindistancesize = (1. / (LimitHits-1));//number of separation between intervals = number of intervals -1, as usual...
		int bindistance = ((int) (RelativeBarycenter / bindistancesize));
		if(bindistance >= LimitHits) bindistance=LimitHits-1;
		//		    cout<<itrk<<"/"<<LimitTracks<<" "<<bindistance<<"/"<<LimitHits<<" "<<Barycenter<<" "<<Dist[2]<<" "<<Dist[0]<<" "<<Dist[1]<<" "<<dx<<" "<<recon->angle[itrk]<<" "<<TrackSample<<" "<<Vec.size()<<endl;
		//for(int i=0;i<Vec.size();i++) cout<<Vec[i].pln<<" "<<Vec[i].ch<<" "<<Vec[i].view<<" "<<Vec[i].z<<endl; 
		NViewsPerPlaneEnergyDepositionNonIsolated[itrk][bindistance]++;
	      }
#endif

	      //4. Estimate the MuCL by going plane by plane. Compare to the expectation using charge deposition
	      //Isolated
	      if(Plane[0][ipln][iview]!=0 || Plane[1][ipln][iview]!=0){//case PM & active plane
#ifdef LIKELIHOODHERE
		if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){//Should use SciBar or INGRID scintillator PDF
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
	      
#endif//End LIKELIHOODHERE
		
#ifdef MVA
		double Barycenter = (DistanceBarycenter[0][ipln][iview]+DistanceBarycenter[1][ipln][iview]) / (Plane[0][ipln][iview]+Plane[1][ipln][iview]);
		double RelativeBarycenter = Barycenter / Dist[2];
		double bindistancesize = (1. / (LimitHits-1));//number of separation between intervals = number of intervals -1, as usual...
		int bindistance = ((int) (RelativeBarycenter / bindistancesize));
		if(bindistance >= LimitHits) bindistance=LimitHits-1;

#ifdef DEBUG_PID_BDT
		cout<<"PM, pln="<<ipln<<", view="<<iview<<", Distance="<<RelativeBarycenter<<", bindistance="<<bindistance<<", Energy deposition="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<endl;
#endif
		EnergyDeposition[itrk][bindistance]+=(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
		NViewsPerPlaneEnergyDeposition[itrk][bindistance]++;
		    
#endif
	      }
	      if(Plane[2][ipln][iview]!=0){
		//cout<<"INGRID, plane="<<ipln<<", view="<<iview<<", PE="<<PEPlane[2][ipln][iview]<<endl;
#ifdef LIKELIHOODHERE		
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
		NHits++;

		if(cl<0.1){
		  HighPECount++;
		  AverageHighPE+=PEPlane[2][ipln][iview];
		  if(PEPlane[2][ipln][iview]>MaxPE) MaxPE=PEPlane[2][ipln][iview];
		}
#endif
		//CLLikelihood_Muon *= CL_Ing_Muon->Eval(PEPlane[2][ipln][iview]);
		//CLLikelihood_NotMuon *= CL_Ing_NotMuon->Eval(PEPlane[2][ipln][iview]);
#ifdef MVA
		double Barycenter = (DistanceBarycenter[2][ipln][iview]) / (Plane[2][ipln][iview]);
		double RelativeBarycenter = Barycenter / Dist[2];
		double bindistancesize = (1. / (LimitHits-1));//number of separation between intervals = number of intervals -1, as usual...
		int bindistance = ((int) (RelativeBarycenter / bindistancesize));
		if(bindistance >= LimitHits) bindistance=LimitHits-1;

#ifdef DEBUG_PID_BDT
		cout<<"INGRID, pln="<<ipln<<", view="<<iview<<", Distance="<<RelativeBarycenter<<", bindistance="<<bindistance<<", Energy deposition="<<PEPlane[2][ipln][iview]<<endl;
#endif
		EnergyDeposition[itrk][bindistance]+=PEPlane[2][ipln][iview];
		NViewsPerPlaneEnergyDeposition[itrk][bindistance]++;
#endif
	      }
	    }
	    //cout<<CLLikelihood_Muon<<", "<<CLLikelihood_NotMuon<<endl;
	  }
	
#ifdef LIKELIHOODHERE
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
	    
	  CL_Likelihood=CLLikelihood_Muon*PMuon/(CLLikelihood_Muon*PMuon+CLLikelihood_NotMuon*P_NotMuon);
	  //cout<<"MuCL plan="<<MuCL<<", pdg="<<Particle<<endl<<", Likelihood Mucl="<<CL_Likelihood<<endl;
#endif
#endif
	}//if PM

	else {
	  double cllikelihood_muon,cllikelihood_notmuon;
	  double CLLikelihood_Muon=1,CLLikelihood_NotMuon=1;
	  double cllikelihood_pion,cllikelihood_notpion;
	  double CLLikelihood_Pion=1,CLLikelihood_NotPion=1;

	  int BinAnglePEWidth_WM;
	  if(ErrorType==4) BinAnglePEWidth_WM=PEAngleData_WM->GetBinWidth(1);

	  for(int i=0;i<Vec.size();i++){
	    if(Vec[i].used>1) continue;

#ifdef MVA
	    double Distance = TMath::Sqrt(pow(Vec[i].x-Vec.front().x,2)+pow(Vec[i].y-Vec.front().y,2)+pow(Vec[i].z-Vec.front().z,2));
#endif
	    if(Vec[i].mod==15){
	      bool grid=(Vec[i].ch>=40);
	      double angle2D=(Vec[i].view==0?recon->thetax[itrk]:recon->thetay[itrk]);
	      double normalAngle=RadDeg(Rec->GetNormalAngleRad(grid,DegRad(recon->angle[itrk]),DegRad(angle2D)));
	      int normalAngleSlice=0;
	      if(normalAngle>80) normalAngleSlice=2;
	      else if(normalAngle>70) normalAngleSlice=1;

	      if(Vec[i].pecorr<5) continue; // avoid hits with de/dz < 5 pe/3mm
	      double PECorrected=Vec[i].pecorr;
	      if(ErrorType==4){
		int BinAngle=((int) normalAngle/BinAnglePEWidth_WM)+1;
		if(PEAngleMC_WM->GetBinContent(BinAngle)!=0 && PEAngleData_WM->GetBinContent(BinAngle)!=0) 
		  SystematicsPECorrected_WM=(PEAngleData_WM->GetBinContent(BinAngle)-PEAngleMC_WM->GetBinContent(BinAngle))/PEAngleMC_WM->GetBinContent(BinAngle);
		else SystematicsPECorrected_WM=0.;

		PECorrected=PECorrected+PECorrected*SystematicsPECorrected_WM;
	      }	      

	      double BirksCorr=1.;double MEV2PE_WM=28.; 
	      double dedz=Vec[i].pecorr/.3/MEV2PE_WM/(1.+.18);
	      if(BirksRescaled){
		BirksCorr=1-dedz*(CBIRKS_values[birksindex]-CBIRKS_values[1])/(1.+dedz*CBIRKS_values[birksindex]);
		PECorrected*=BirksCorr;
		cout<<"de/dz="<<dedz<<" (Mev/cm), Delta CB="<<CBIRKS_values[birksindex]-CBIRKS_values[1]<<" BirksCorr="<<BirksCorr<<endl;
	      }
#ifdef LIKELIHOODHERE
	      cllikelihood_muon=CL_WM_Muon->Eval(PECorrected,normalAngleSlice);
	      cllikelihood_notmuon=CL_WM_NotMuon->Eval(PECorrected,normalAngleSlice);
	      CLLikelihood_Muon*=cllikelihood_muon;
	      CLLikelihood_NotMuon*=cllikelihood_notmuon;
#ifdef PI_LIKELIHOOD
	      cllikelihood_pion=CL_WM_Pion->Eval(PECorrected,normalAngleSlice);
	      cllikelihood_notpion=CL_WM_NotPion->Eval(PECorrected,normalAngleSlice);
	      CLLikelihood_Pion*=cllikelihood_pion;
	      CLLikelihood_NotPion*=cllikelihood_notpion;
#endif
#endif

#ifdef MVA
	      double RelativePosition = Distance/TrkLength;
	      double bindistancesize = (1. / (LimitHits-1));//number of separation between intervals = number of intervals -1, as usual...
	      int bindistance = ((int) (RelativePosition / bindistancesize));
	      if(bindistance >= LimitHits) bindistance=LimitHits-1;
	      EnergyDeposition[itrk][bindistance]+=PECorrected;
	      NViewsPerPlaneEnergyDepositionNonIsolated[itrk][bindistance]++;
#ifdef DEBUG_PID
	      cout<<"WM, Track Length: "<< TrkLength <<", bin: "<< bindistance <<", Distance from the vertex: "<<Distance<<", Relative pos: "<<RelativePosition<<", Energy Dep: "<<PECorrected<<endl;
#endif
#endif
	    }
	    else{
	      if(Vec[i].pecorr<7) continue; // avoid hits with de/dz < 7 pe/cm

	      double PECorrected=Vec[i].pecorr;
	      double BirksCorr=1.;double MEV2PE=46.; 
	      double dedz=Vec[i].pecorr/MEV2PE/(1.+.09)/1.08;
	      if(ErrorType==4) PECorrected=PECorrected+PECorrected*SystematicsPECorrected_Ing;
	      if(BirksRescaled){
		BirksCorr=1-dedz*(CBIRKS_values[birksindex]-CBIRKS_values[1])/(1.+dedz*CBIRKS_values[birksindex]);
		PECorrected*=BirksCorr;
	      }
	      
#ifdef LIKELIHOODHERE
	      cllikelihood_muon=CL_Ing_Muon->Eval(PECorrected);
	      cllikelihood_notmuon=CL_Ing_NotMuon->Eval(PECorrected);
	      CLLikelihood_Muon*=cllikelihood_muon;
	      CLLikelihood_NotMuon*=cllikelihood_notmuon;
#ifdef PI_LIKELIHOOD
	      cllikelihood_pion=CL_Ing_Pion->Eval(PECorrected);
	      cllikelihood_notpion=CL_Ing_NotPion->Eval(PECorrected);
	      CLLikelihood_Pion*=cllikelihood_pion;
	      CLLikelihood_NotPion*=cllikelihood_notpion;	      
#endif
#endif
#ifdef MVA
	  double RelativePosition = Distance/TrkLength;
	  double bindistancesize = (1. / (LimitHits-1));//number of separation between intervals = number of intervals -1, as usual...
	  int bindistance = ((int) (RelativePosition / bindistancesize));
	  if(bindistance >= LimitHits) bindistance=LimitHits-1;
	  EnergyDeposition[itrk][bindistance]+=PECorrected;
	  NViewsPerPlaneEnergyDepositionNonIsolated[itrk][bindistance]++;
	  //cout<<"INGRID, Distance from the vertex: "<<Distance<<", Relative pos: "<<RelativePosition<<", Energy Dep: "<<PECorrected<<endl;
#ifdef DEBUG_PID
	  cout<<"INGRID, Track Length: "<< TrkLength <<", bin: "<< bindistance <<", Distance from the vertex: "<<Distance<<", Relative pos: "<<RelativePosition<<", Energy Dep: "<<PECorrected<<endl;
#endif
#endif
	    }
	  }
#ifdef LIKELIHOOD
	  CL_Likelihood=CLLikelihood_Muon*PMuon/(CLLikelihood_Muon*PMuon+CLLikelihood_NotMuon*P_NotMuon);
#ifdef PI_LIKELIHOOD	
	  CL_Likelihood_Pi=CLLikelihood_Pion*PPion/(CLLikelihood_Pion*PPion+CLLikelihood_NotPion*P_NotPion);
#endif
#endif
	  //  if(CL_Likelihood<0 || CL_Likelihood>1) cout<<CLLikelihood_Muon<<" "<<CLLikelihood_NotMuon<<endl;
	  //  cout<<CL_Likelihood<<endl;
	}


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


	//if(Particle==2112) cout<<"event "<<ievt<<" reco "<<irec<<" track="<<itrk<<" is associated to neutron -- mucl="<<CL_Likelihood<<endl; 
	    
	TrackAngle[itrk]=Rec->GetBeamAngle(DegRad(recon->angle[itrk]),DegRad(recon->thetay[itrk]),DegRad(recon->thetax[itrk]));
	TrackThetaX[itrk]=(recon->thetax)[itrk];
	TrackThetaY[itrk]=(recon->thetay)[itrk];
	// ML: TrackAngle is wrt beam direction; ThetaX,ThetaY are wrt z axis (horizontal)
	TrackWidth[itrk]=Rec->Reconstruction::GetINGRIDTrackWidth(Vec);
	GT[itrk]=Geom;
	IsReconstructed[itrk]=true;
	TypeOfTrack[itrk]=Particle;
#ifdef PI_LIKELIHOOD
	CLMuon[itrk]=CL_Likelihood_Pi;
#else
	CLMuon[itrk]=CL;
#endif
	//if(CL_Plan!=-1) CLMuon_Plan[itrk]=TMath::Log(CL_Plan);
	CLMuon_Plan[itrk]=CL_Plan;
	CLMuon_Likelihood[itrk]=CL_Likelihood;
	if(CLMuon_Plan[itrk]<0.05){ // ML 20170412
	  LowCL+=weight;
	  nLowCL++;
	}
	else if(CLMuon_Plan[itrk]>0.95){ // ML 20170412
	  HighCL+=weight;
	  nHighCL++;
	}
	CLMuon_KS[itrk]=(recon->mucl)[itrk];
#ifdef DEBUG_PID
	cout<<CLMuon_KS[itrk]<<endl;
	if(recon->mucl[itrk]!=recon->mucl[itrk]){ cout<<"*************"<<endl;nbad++;}
#endif
	Momentum[itrk]=TrueParticleNRJ;
	PD[itrk]=Dist[0];
	if(TrackSample>=2) ID[itrk]=Dist[1];
	else ID[itrk]=0;
	Sample[itrk]=TrackSample;
	TotalCharge[itrk]=Rec->Reconstruction::GetTrackEnergyPerDistance(Vec,dx);
	//cout<<"Energy of the track="<<Rec->Reconstruction::GetTrackEnergy(Vec)<<endl;


	//5. Determine the unit vector, opening angle and coplanarity
	if(nTracks==2){
	  vector <double> Kinematic;
	  Kinematic.clear();
	  for(int itrk2=0;itrk2<nTracks;itrk2++){
	    if(itrk2==itrk) continue;
	    else{
	      Kinematic=Rec->Reconstruction::GetKinematic(DegRad((recon->angle)[itrk]), DegRad((recon->thetax)[itrk]), DegRad((recon->thetay)[itrk]), DegRad((recon->angle)[itrk2]), DegRad((recon->thetax)[itrk2]), DegRad((recon->thetay)[itrk2]));
		  
	    }
	  }
	  OpeningAngle=RadDeg(Kinematic[0]);
	  CoplanarityAngle=RadDeg(Kinematic[1]);
	  //cout<<"Opening="<<OpeningAngle<<", coplanarity="<<CoplanarityAngle<<endl;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef DEBUG_PID_BDT
	//cout<<"Total number of hits="<<recon->NhitTs(itrk)<<", hits of vector="<<Vec.size()<<endl;
	//for(int i=0;i<Vec.size();i++){
	//cout<<"Module="<<Vec[i].mod<<", plane="<<Vec[i].pln<<", channel="<<Vec[i].ch<<", view="<<Vec[i].view<<", pe="<<Vec[i].pe<<", corr="<<Vec[i].pecorr<</*", Distance="<<DistanceBarycenter/Dist[2]<<*/endl;
	//}
	    
	cout<<"Type of particle="<<TypeOfTrack[itrk]<<", sample="<<Sample[itrk]<<", equivalent iron distance="<<ID[itrk]+(PD[itrk]/IronCarbonRatio)<<", distance in CH/Fe="<<PD[itrk]<<"/"<<ID[itrk]<<endl;
	    
	//Energydeposition with length
	for(int ihit=0;ihit<LimitHits;ihit++){
	  if(ihit==0) cout<<setprecision(3)<<"dE/dx = [ ";
	  if(ihit<LimitHits-1) cout<<EnergyDeposition[itrk][ihit]<<", ";
	  else cout<<EnergyDeposition[itrk][ihit]<<" ]"<<endl;
	}

	//NViews per plane
	for(int ihit=0;ihit<LimitHits;ihit++){
	  if(ihit==0) cout<<setprecision(3)<<"Number of views = [ ";
	  if(ihit<LimitHits-1) cout<<NViewsPerPlaneEnergyDeposition[itrk][ihit]<<", ";
	  else cout<<NViewsPerPlaneEnergyDeposition[itrk][ihit]<<" ]"<<endl;
	}
#endif


#ifdef MVA
	for(int ihit=0;ihit<LimitHits;ihit++){
	  //For EnergyDeposition
	  if(NViewsPerPlaneEnergyDeposition[itrk][ihit]!=0){
	    EnergyDeposition[itrk][ihit]/=NViewsPerPlaneEnergyDeposition[itrk][ihit];
	  }
	}

#endif

      
	Vec.clear();
	Vec2.clear();
	VecAll.clear();
      	    
#ifdef MVA
	TrackAngleMVA=TrackAngle[itrk];
	TrackWidthMVA=TrackWidth[itrk];
	TypeOfTrackMVA=TypeOfTrack[itrk];
	TotalChargeMVA=TotalCharge[itrk];
	CLMuonMVA=CLMuon[itrk];
	CLMuon_PlanMVA=CLMuon_Plan[itrk];
	CLMuon_KSMVA=CLMuon_KS[itrk];
	CLMuon_LikelihoodMVA=CLMuon_Likelihood[itrk];
	MomentumMVA=Momentum[itrk];
	IDMVA=ID[itrk];
	PDMVA=PD[itrk];
	SampleMVA=Sample[itrk];
	IsReconstructedMVA=IsReconstructed[itrk];
	GTMVA=GT[itrk];
	for(int ihit=0;ihit<LimitHits;ihit++){
	  // cout<<"hit: "<<ihit<<", Energy: "<<EnergyDeposition[itrk][ihit]<<endl;
	  EnergyDepositionMVA[ihit]=EnergyDeposition[itrk][ihit];
	}
	wtreeMVA->Fill();
	ResetInputVariablesMVA();
#endif
	
      }//Tracks
      //xcout<<"Is Detected="<<VIsDetected<<endl;
      wtree->Fill();
      ResetInputVariables(false);	
	
    }//Recons
  }//Evt
  if(XSEC) {
    _file1->Close();
  }

  cout<<"writing"<<endl;
  cout<<"Low CL="<<LowCL<<", High="<<HighCL<<endl;
  cout<<"Nb Low CL="<<nLowCL<<", Nb High="<<nHighCL<<endl;
  cout<<"Nb bad mucl="<<nbad<<endl;

#ifdef DEBUG
  /*
    f_PMIng_Plan->SetNpx(300);
    f_PMSci_Plan->SetNpx(300);
    f_Ing_Plan->SetNpx(300);
    f_PMIng_Plan->Write();
    f_PMSci_Plan->Write();
    f_Ing_Plan->Write();
  */
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

#ifdef DISP
  if(Disp) theApp.Run();
#endif
  return(0);
}
  

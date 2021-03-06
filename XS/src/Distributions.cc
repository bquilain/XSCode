
#ifndef Distributions_cc
#define Distributions_cc
/************CLibs********/
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
/********ROOT Libs*******/
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
#include <TH1.h>
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
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TSpline.h>
/**********INGRID Libs*************/
#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "INGRID_Dimension.cc"
/*********Hit Libs****************/
#include "Hit.h"
#include "Distributions.h"
#include "setup.h"
TSpline3 * s_PMIng;
TSpline3 * s_PMSci;
TSpline3 * s_Ing;
TF1 * f_PMIng;
TF1 * f_PMSci;
TF1 * f_Ing;
TF1 * CL_PMIng;
TF1 * CL_PMSci;
TF1 * CL_Ing;

Double_t Distributions::fPDF_PMIng(Double_t *pe, Double_t *par){
  return s_PMIng->Eval(pe[0])/par[0];
}
Double_t Distributions::fPDF_PMSci(Double_t *pe, Double_t *par){
  return s_PMSci->Eval(pe[0])/par[0];
}
Double_t Distributions::fPDF_Ing(Double_t *pe, Double_t *par){
  return s_Ing->Eval(pe[0])/par[0];
}

Double_t Distributions::Cumulative_PMIng(Double_t *pe, Double_t *par){
  if(pe[0]>=par[1]){ return (1.-f_PMIng->Integral(par[1],pe[0])/f_PMIng->Integral(par[1],par[2]));}
  else{ return (1.-f_PMIng->Integral(pe[0],par[1])/f_PMIng->Integral(par[0],par[1]));}
}
Double_t Distributions::Cumulative_PMSci(Double_t *pe, Double_t *par){
   if(pe[0]>=par[1]){ return (1.-f_PMSci->Integral(par[1],pe[0])/f_PMSci->Integral(par[1],par[2]));}
  else{ return (1.-f_PMSci->Integral(pe[0],par[1])/f_PMSci->Integral(par[0],par[1]));}
}
Double_t Distributions::Cumulative_Ing(Double_t *pe, Double_t *par){
  if(pe[0]>=par[1]){ return (1.-f_Ing->Integral(par[1],pe[0])/f_Ing->Integral(par[1],par[2]));}
  else{ return (1.-f_Ing->Integral(pe[0],par[1])/f_Ing->Integral(par[0],par[1]));}
}

#endif

#ifndef setup_h
#define setup_h
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
#include "Reconstruction.h"
double INGRIDSCIBAR=2.0;
Int_t Scyc =  4;
Int_t Ncyc = 12;
int NMod=17;
int NPln=11;
int NPlnPM=18;
int NView=2;
int NCh=24;
int NChPM=32;
double IronCarbonRatio=7.87/1.03;
double MeV2PE=46.0;
double MeV2PEPM=38.6;
double Beam[3]={0,-TMath::Sin(3.8*TMath::Pi()/180.),TMath::Cos(3.8*TMath::Pi()/180.)};
double C[17]={-1.85,0.34,0.59,0.74,0.514,-0.37,1.25,-0.06,0.562,0.82,-0.47,0.6,-0.57,-0.45};
const double Mmu=105.66/1000.;
const double Mp=938.272/1000.;
const double Mn=939.5659/1000.;
const int LimitTracks = 5;
const int LimitHits = 25;
const int NDials=175;

const int StartRun=14000;
const int EndRun=15999;
const int StartSubRun=0;
const int EndSubRun=300;
const int StartRunList=14;//Which list will be read (14=>14000.txt). Necessary to use only run processed
const int EndRunList=15;//Which list will be read (14=>14000.txt). Necessary to use only run processed
double NMCfiles=500;
double DataPOT=0.58;//In units of 10^21 POT
const int StartError=0;
const int EndError=0;//17;//34;//34;//41 for 2012 ;//17
int NFluxFiles;
int StartXsec=0;int EndXsec=24;int NXsecVariations=7;
int CenterXsecVariations=(int) (NXsecVariations-1-((double) (NXsecVariations-1)/2));

int NE[EndError+1];
double Step[EndError+1];
double Start[EndError+1];
double End[EndError+1];
double Nominal; double Err;
const int Systematics_Detector_Start=1;
const int Systematics_Detector_End=15;
const int Systematics_Flux_Start=16;
const int Systematics_Flux_End=16;
const int Systematics_Xsec_Start=17;
const int Systematics_Xsec_End=Systematics_Xsec_Start+24;
bool EStatistics=true;//if true, estimate stat. error after unfolding
const int NStatisticalVariations=1000;//number of stat. varied toy experiments to evaluate the stat. error.

//For particle gun
const int npdg=4;
int pdgValues[npdg]={13,211,-211,2212};

const int NBinsEnergyFlux=43;
double BinningEnergyFlux[NBinsEnergyFlux+1];
const int NBinsTrueEnergy=6;
const int NBinsRecEnergy=6;
double BinningTrueEnergy[NBinsTrueEnergy+1];
double BinningRecEnergy[NBinsRecEnergy+1];
const int NFSIs=12;//cc0pi+0p,cc0pi+1p,cc0pi+morep,cc1pi,ccnpi,ccpi0,nc,+all bkg
const int NBinsTrueMom=5;
const int NBinsTrueAngle=5;
const int NBinsRecMom=17;
const int NBinsRecAngle=30;
double BinningTrueMom[NBinsTrueMom+1];
double BinningTrueAngle[NBinsTrueAngle+1];
double BinningRecMom[NBinsRecMom+1];
double BinningRecAngle[NBinsRecAngle+1];


/*

const int NEnergyBins=43;
double EnergyBins[NBinsEnergyFlux+1];

const int NInteractionTypes=11;

const int NEnergyBins=43;

const int NBinsEnu=6;
const int NBinsIron=17;
const int NBinsMom=5;
const int NBinsAngle=4;
const int NBinsAngleRec=4;
double DistIronBin[NBinsIron+1];
double DistMomBin[NBinsMom+1];
double DistAngleBin[NBinsAngle+1];
double DistAngleRecBin[NBinsAngle+1];


const int NBinsTrueAngle=90;
double DistTrueAngleBin[NBinsAngle+1];



const int _NBinsMom=5;//2014/11/15
//const int _NBinsMom=10;
const int _NBinsIron=17;
double _DistIronBin[_NBinsIron+1];
double _DistMomBin[_NBinsMom+1];
const int _NBinsAngle=10;
double _DistAngleBin[_NBinsAngle+1];
const int _NBinsTrueAngle=10;
double _DistTrueAngleBin[_NBinsAngle+1];
*/

#endif

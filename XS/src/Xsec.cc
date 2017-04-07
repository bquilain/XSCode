#ifndef Xsec_cc
#define Xsec_cc
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
#include <TRandom3.h>
#include <TChain.h>
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
#include "Xsec.h"
#include "Reconstruction.cc"
#include "setup.h"
INGRID_Dimension * IngD = new INGRID_Dimension();
Reconstruction * Reco = new Reconstruction();
//#define DEBUG
void Xsec::Initialize(){

  for(int i=0;i<=NBinsEnergyFlux;i++){
    if(i==0) BinningEnergyFlux[i]=0;
    else if(i>=1 && i<36) BinningEnergyFlux[i]=0.5+(i-1)*0.1;
    else if(i>=36 && i<42) BinningEnergyFlux[i]=4+(i-36)*1;
    BinningEnergyFlux[42]=10; BinningEnergyFlux[NBinsEnergyFlux]=30;
  }


    for(int i=0;i<NBinsTrueMom+1;i++){
      BinningTrueMom[0]=0;
      BinningTrueMom[1]=0.5;
      BinningTrueMom[2]=0.7;
      BinningTrueMom[3]=1.0;
      BinningTrueMom[4]=5.0;
      BinningTrueMom[5]=30.0;
    }

    for(int i=0;i<NBinsRecMom+1;i++){
      BinningRecMom[0]=0;
      BinningRecMom[1]=10;    
      if(i>1 && i<NBinsRecMom) BinningRecMom[i]=10+5*(i-1);
      if(i==NBinsRecMom) BinningRecMom[i]=150;
    }

    BinningTrueAngle[0]=0;
    BinningTrueAngle[1]=10;
    BinningTrueAngle[2]=20;
    BinningTrueAngle[3]=30;
    BinningTrueAngle[4]=60;
    BinningTrueAngle[5]=180;

    for(int i=0;i<NBinsRecAngle+1;i++){
      BinningRecAngle[i]=3*i;
    }
    //BinningRecAngle[0]=0;
    //BinningRecAngle[1]=180;

    /*
    BinningTrueMom[0]=0;
    BinningTrueMom[1]=0.8;
    BinningTrueMom[2]=30;

    BinningRecMom[0]=0.;
    BinningRecMom[1]=50.;    
    BinningRecMom[2]=150;    

  for(int i=0;i<NBinsTrueMom+1;i++){
    BinningTrueMom[0]=0;
    BinningTrueMom[1]=0.5;
    BinningTrueMom[2]=0.7;
    BinningTrueMom[3]=1.0;
    BinningTrueMom[4]=5.0;
    BinningTrueMom[5]=30.0;
  }



  for(int i=0;i<NBinsRecAngle+1;i++){
    BinningRecAngle[0]=0;
    BinningRecAngle[1]=20;
    BinningRecAngle[2]=30;
    BinningRecAngle[3]=60;
    BinningRecAngle[4]=180;
  }
*/
    ////////////////////////////INITIALIZE ERROR NOW////////////////////////////////////
    for(int n=StartError;n<=EndError;n++){
      NE[n]=1;
      Start[n]=1;
      Step[n]=1;

      if(Systematics_Detector_End>=Systematics_Flux_Start || Systematics_Flux_Start>=Systematics_Xsec_Start) cout<<"################################################################################################################### STOP THE PROCESS, THERE IS A PROBLEM IN THE NUMBERING OF YOUR SYST. ERROR SOURCES. CHECK SYSTEMATICS_DETECTOR_START/END...###########################################################################"<<endl;
	 
      if(n>=Systematics_Detector_Start && n<=Systematics_Detector_End){//16: flux error
      if(n==1){//2: Stat
	Start[n]=0;
	NE[n]=1000;
	Step[n]=1;
      }
      else if(n==2){//2: DN
	//NE[n]=2;//Number of DN values tested
	Start[n]=0;
	NE[n]=10;
	Step[n]=1;
      }
      else if(n==3){//3: Hit efficiency
	NE[n]=1;//Only one probability to mask the channel, given by data - MC. So we try to recover data inefficiency from MC distribution. So only one variation, the one needed to recover data.
	//NE[n]=10;//TEMP
      }
      else if(n==4){//4: absolute difference MC/data light yield with track angle
	NE[n]=1;
      }
      else if(n==5) NE[n]=2;//Birks constant variations
      else if(n==6){//6: Beam related background (in fact, this mainly evaluate sand muons)
	Nominal=0.6;
	Err=TMath::Sqrt(Nominal*Nominal+0.2*0.2+0.2*0.2);
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==7){//7: 2D reconstruction error
	Nominal=3;
	Err=1;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=5;
      }
      else if(n==8){    //8: VetoUpstreamCriteria: nominal=0 planes, vary from 0->2 per 1 plane step
	Nominal=0;
	Err=1;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==9){//9: VetoEdgeCriteria: nominal=80cm, vary 70->100 per 10cm steps
	Nominal=80;
	Err=10;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=4;
      }
      else if(n==11){//11: Vertexing, plane tolerance: nominal=2, vary 2->4 per 1 plane steps
	Nominal=2;
	Err=1;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==12){//12: Vertexing, transverse tolerance: nominal=15cm, vary 15cm->20cm per 2.5cm steps
	Nominal=15;
	Err=2.5;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==13){//13: Track matching, plane tolerance: nominal=4, vary 3->5 per 1 plane steps
	Nominal=4;
	Err=1;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==14){//14: INGRID/PM tracks angle matching: nominal=35°, vary 30°->40° per 5° steps
	Nominal=35;
	Err=5;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      else if(n==15){//15: INGRID/PM tracks transverse position matching: nominal=8.5cm, vary 7.5cm->9.5 per 1cm steps
	Nominal=8.5;
	Err=1.;
	Start[n]=Nominal-Err;
	Step[n]=Err;
	NE[n]=3;
      }
      }
      else if(n>=Systematics_Flux_Start && n<=Systematics_Flux_End){//16: flux error
	Nominal=0;
	Err=1;
	Start[n]=Nominal;
	Step[n]=Err;
	NE[n]=500;
	NFluxFiles=NE[n];
      }
      else if(n>=Systematics_Xsec_Start && n<=Systematics_Xsec_End){//17: cross section error
	//one should separate the place in the reweight vector, and the real starts of the Xsection.
	//Real start=0. Place in the Xsection vector is (StartXsec+n)*NXsecVariations.
	Start[n]=(StartXsec+(n-Systematics_Xsec_Start))*NXsecVariations;
	Step[n]=1;
	NE[n]=NXsecVariations;
      }
      End[n]=Start[n]+(NE[n]+1)*Step[n];
    }
    
//0. No Error, nominal case
//1. TO DO
//2. Dark noise, variations
//3. Hit efficiency, variations within the difference data and MC
//4. Light yield, variation of PE with angle btw data and MC.
//5. Light yield, Birks quenching effect.
//6: Beam related background (in fact, this mainly evaluate sand muons)
//7: 2D reconstruction error
//8: VetoUpstreamCriteria: nominal=0 planes, vary from 0->2 per 1 plane step
//9: VetoEdgeCriteria: nominal=80cm, vary 70->100 per 10cm steps
//10: FVCriteria: nominal=100cm, vary 50->100 per 10 cms steps
//11: Vertexing, plane tolerance: nominal=2, vary 2->4 per 1 plane steps
//12: Vertexing, transverse tolerance: nominal=15cm, vary 15cm->20cm per 2.5cm steps
//13: Track matching, plane tolerance: nominal=4, vary 3->5 per 1 plane steps
//14: INGRID/PM tracks angle matching: nominal=35°, vary 30°->40° per 5° steps
//15: INGRID/PM tracks transverse position matching: nominal=8.5cm, vary 7.5cm->9.5 per 1cm steps
//16:Flux
//17: Xsec

}


void Xsec::DetermineNuType(int IsSand,int IsAnti,int IsNuE,int IsBkgH,int IsBkgV,int nutype, int intmode){
  IsSand=0;
  IsAnti=0;
  IsNuE=0;
  IsBkgH=0;
  IsBkgV=0;
  
  if(intmode==16){
    //if(nutype==1) return 0;
    if(nutype==2) IsAnti=1;
    else if(nutype==3) IsNuE=1;
  }
  else if(intmode<7) IsBkgH=1;
  else if(intmode<14) IsBkgV=1;
  else IsSand=1;
}


int Xsec::DetermineFSI(int IsSand,int IsAnti,int IsNuE,int IsBkgH,int IsBkgV,IngridEventSummary * evt){
  int FSIMuons=0;int FSIPions=0;int FSIProtons=0;int FSIPions0=0;int FSIOther=0;int FSINeutralPions=0;
  int FSIInt;
  IngridSimParticleSummary * SimPart;

  for(int is=0;is<evt->NIngridSimParticles();is++){
    SimPart=(IngridSimParticleSummary*) evt->GetSimParticle(is);
    //cout<<"Particle type="<<SimPart->pdg<<", NRJ="<<SimPart->momentum[0]<<", "<<SimPart->momentum[3]<<endl;

    if(TMath::Abs(SimPart->pdg)==211) FSIPions++;
    else if(TMath::Abs(SimPart->pdg)==111) FSINeutralPions++;
    else if(TMath::Abs(SimPart->pdg)==2212) FSIProtons++;
    else if(TMath::Abs(SimPart->pdg)==13) FSIMuons++;
    else if(TMath::Abs(SimPart->pdg)==111) FSIPions0++;
    else if(TMath::Abs(SimPart->pdg)!=2112) FSIOther++;
  }
  /*
  if(FSIMuons==1){
    if(FSINeutralPions==0){
      if(FSIPions==0) FSIInt=0+FSIProtons;
      if(FSIPions==1) FSIInt=10+FSIProtons;
      else FSIInt=20+FSIProtons;
    }
    else FSIInt=30+FSIProtons;
  }
  else{
    if(FSINeutralPions==0) FSIInt=50+FSIProtons;
    else FSIInt=60+FSIProtons;
  }
  */

  if(FSIMuons==1){
    if(FSINeutralPions==0){
      if(FSIPions==0){
	if(FSIProtons<2) FSIInt=FSIProtons;
	else FSIInt=2;
      }
      else if(FSIPions==1) FSIInt=3;
      else FSIInt=4;
    }
    else if(FSINeutralPions==1){
      if(FSIPions==0) FSIInt=5;
      else FSIInt=4;//If more pion that only one pi^0, consider that it is a Npions (example: CC 1pi+ + 1pi^0) 
    }
    else FSIInt=4;
  }
  else FSIInt=6;

  if(IsSand) FSIInt=7;
  else if(IsBkgH) FSIInt=8;
  else if(IsBkgV) FSIInt=9;
  else if(IsAnti) FSIInt=10;
  else if(IsNuE) FSIInt=11;
  
  return FSIInt;
}



double Xsec::GetMuCL(TF1 * CL_Ing, TF1 * CL_PMIng, TF1 * CL_PMSci, vector <Hit3D> Vec, double dx, int TrackSample, bool SystematicPE, int RandomPE){
  
  double CL=1;
  double PECorrected; 
  int NHits[3]={0,0,0};
  
  for(int i=0;i<Vec.size();i++){
    if(Vec[i].used>1) continue;
    //cout<<Vec[i].pecorr<<", "<<dx<<endl;
    
    if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
      PECorrected=Vec[i].pecorr;
      CL *= CL_PMIng->Eval(PECorrected);
      //cout<<"PMIng, MuCL="<<CL_PMIng->Eval(PECorrected)<<endl;
      NHits[1]++;
    }
   
    else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
      PECorrected=Vec[i].pecorr;
      CL *= CL_PMSci->Eval(PECorrected);      
      //cout<<"PMSci, MuCL="<<CL_PMIng->Eval(PECorrected)<<endl;
      NHits[2]++;
    }
    else if(!(Vec[i].mod==16)){
      PECorrected=Vec[i].pecorr;
      CL *= CL_Ing->Eval(PECorrected);      
      //cout<<"Ing, MuCL="<<CL_PMIng->Eval(PECorrected)<<endl;
      NHits[0]++;
    }
  }//End of loop over hits
  double Sum=0;
  for(int i=0;i<(NHits[0]+NHits[1]+NHits[2]);i++){
  //for(int i=0;i<NHits[1];i++){
    Sum+=(pow(-TMath::Log(CL),i)/TMath::Factorial(i));
  }
  //cout<<"CL Temp="<<CL<<", Sum="<<Sum<<", NHits="<<(NHits[0]+NHits[1]+NHits[2])<<endl;
  double MuCL=CL*Sum;
  //cout<<"CL final="<<MuCL<<endl<<endl;  
  return MuCL;
}

double Xsec::GetMuCL_Plan(TF1 * CL_Ing, TF1 * CL_PMIng, TF1 * CL_PMSci, vector <Hit3D> Vec, double dx, int TrackSample, bool SystematicPE, int RandomPE){
  double CL=1;
  double PECorrected; 
  int NHits=0;

  double PEPlane[3][NPlnPM][2];
  int Plane[3][NPlnPM][2];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int i=0;i<3;i++){
      for(int iview=0;iview<2;iview++){
      PEPlane[i][ipln][iview]=0;
      Plane[i][ipln][iview]=0;
      }
    }
  }
  
  //cout<<CL_PMIng->Eval(50)<<endl;
  for(int i=0;i<Vec.size();i++){

    if(Vec[i].used>1) continue;

    if(Vec[i].mod==16 && Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
      PECorrected=Vec[i].pecorr;
      PEPlane[1][Vec[i].pln][Vec[i].view]+=PECorrected;
      Plane[1][Vec[i].pln][Vec[i].view]++;
      //cout<<Vec[i].mod<<", plane="<<Vec[i].pln<<", view="<<Vec[i].view<<", chan="<<Vec[i].ch<<", x="<<Vec[i].x<<", y="<<Vec[i].y<<", pe="<<Vec[i].pe<<", corr="<<Vec[i].pecorr<<", PECorrected="<<PECorrected<<endl;
      //cout<<"hello"<<endl;
    }
    else if(Vec[i].mod==16 && !(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
      PECorrected=Vec[i].pecorr;
      PEPlane[0][Vec[i].pln][Vec[i].view]+=PECorrected;
      Plane[0][Vec[i].pln][Vec[i].view]++;
      //cout<<Vec[i].mod<<", plane="<<Vec[i].pln<<", view="<<Vec[i].view<<", chan="<<Vec[i].ch<<", x="<<Vec[i].x<<", y="<<Vec[i].y<<", pe="<<Vec[i].pe<<", corr="<<Vec[i].pecorr<<", PECorrected="<<PECorrected<<endl;
    }
    else{
      PECorrected=Vec[i].pecorr;
      PEPlane[2][Vec[i].pln][Vec[i].view]+=PECorrected;
      Plane[2][Vec[i].pln][Vec[i].view]++;
      //cout<<Vec[i].mod<<", plane="<<Vec[i].pln<<", view="<<Vec[i].view<<", chan="<<Vec[i].ch<<", x="<<Vec[i].x<<", y="<<Vec[i].y<<", pe="<<Vec[i].pe<<", corr="<<Vec[i].pecorr<<", PECorrected="<<PECorrected<<endl;
      }
  }//End of loop over hits
  for(int ipln=0;ipln<NPlnPM;ipln++){
    //cout<<"Pln="<<ipln;
    for(int iview=0;iview<2;iview++){
      //cout<<", View="<<iview<<", ";
      if(PEPlane[0][ipln][iview]!=0 || PEPlane[1][ipln][iview]!=0){//case PM & active plane
	NHits++;
	//cout<<"PM ,";
	if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){
	  CL *= CL_PMSci->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
	  //cout<<"SciBar type, PE="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<", CL="<<CL_PMSci->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview])<<endl;
	}
	else{
	  CL *= CL_PMIng->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]);
	  //cout<<"INGRID type, PE="<<PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview]<<", CL="<<CL_PMIng->Eval(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview])<<endl;
	}
      }
    
      else if(Plane[2][ipln][iview]!=0){
	NHits++;
	//cout<<"INGRID, ";
	CL *= CL_Ing->Eval(PEPlane[2][ipln][iview]);
	//cout<<"INGRID type, PE="<<PEPlane[2][ipln][iview]<<", CL="<<CL_PMIng->Eval(PEPlane[2][ipln][iview])<<endl;
	}
    }
  }

  double Sum=0;
  for(int i=0;i<NHits;i++){
  //for(int i=0;i<NHits[1];i++){
    Sum+=(pow(-TMath::Log(CL),i)/TMath::Factorial(i));
  }
  
  //cout<<"CL Temp="<<CL<<", Sum="<<Sum<<", NHits="<<NHits<<endl;
  double MuCL=CL*Sum;
  //cout<<"CL final="<<MuCL<<endl<<endl;
  
  return MuCL;
}

//Read the output file of CC0pi selection, but only the selected data
void Xsec::LoadInputFiles_UnfoldedData(char * InputName, double DataUnfoldedEvents[NBinsTrueMom][NBinsTrueAngle], double DataTrueEvents[NBinsTrueMom][NBinsTrueAngle]){
  
  TChain * wtree = new TChain("wtree");
  wtree->Add(InputName);

  int Iterations, Systematics, Statistics;
  double Events[NBinsTrueMom][NBinsTrueAngle];
  double EventsAll[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double TrueEvents[NBinsTrueMom][NBinsTrueAngle];

  TBranch* Br_Iterations = wtree->GetBranch("Iterations");
  Br_Iterations->SetAddress(&Iterations);
  wtree->SetBranchAddress("Iterations",&Iterations);

  TBranch* Br_Systematics = wtree->GetBranch("Systematics");
  Br_Systematics->SetAddress(&Systematics);
  wtree->SetBranchAddress("Systematics",&Systematics);

  TBranch* Br_Statistics = wtree->GetBranch("Statistics");
  Br_Statistics->SetAddress(&Statistics);
  wtree->SetBranchAddress("Statistics",&Statistics);

  TBranch* Br_Events = wtree->GetBranch("Events");
  Br_Events->SetAddress(Events);
  wtree->SetBranchAddress("Events",Events);

  TBranch* Br_EventsAll = wtree->GetBranch("EventsAll");
  Br_EventsAll->SetAddress(EventsAll);
  wtree->SetBranchAddress("EventsAll",EventsAll);

  TBranch* Br_TrueEvents = wtree->GetBranch("TrueEvents");
  Br_TrueEvents->SetAddress(TrueEvents);
  wtree->SetBranchAddress("TrueEvents",TrueEvents);

  int nevt=wtree->GetEntries();
  
  for(int ievt=0;ievt<nevt;ievt++){
    wtree->GetEntry(ievt);
    
    for(int c0=0;c0<NBinsTrueMom;c0++){
      for(int c1=0;c1<NBinsTrueAngle;c1++){
	DataUnfoldedEvents[c0][c1]=Events[c0][c1];

	double nEvts=0;
	for(int e0=0;e0<NBinsRecMom-1;e0++){//loop over effect 0
	  for(int e1=0;e1<NBinsRecAngle-1;e1++){//loop over effect 1
	    nEvts+=EventsAll[c0][c1][e0][e1];
	  }
	}
	DataTrueEvents[c0][c1]=nEvts;
	
      }
    }
  }
  delete wtree;
}


//Read the output file of CC0pi selection, but only the selected data
void Xsec::LoadInputFiles_OnlySelectedData(char * fDataName,double DataReconstructedEvents[NBinsRecMom][NBinsRecAngle]){

  double check=0;
  ifstream fEvent;

  cout<<"opening file "<<fDataName<<" for data input"<<endl;
  fEvent.open(fDataName,ios::in);
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      fEvent>>DataReconstructedEvents[e0][e1];
#ifdef DEBUG
      cout<<DataReconstructedEvents[e0][e1]<<" ";
#endif
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading data file"<<endl;
  fEvent.close();
}

//Read the output file of CC0pi selection.
void Xsec::LoadInputFiles(char * fDataName,char * fMCName, double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle], double DataReconstructedEvents[NBinsRecMom][NBinsRecAngle],double MCReconstructedEvents[NBinsRecMom][NBinsRecAngle], double MCReconstructedBkgEvents[NBinsRecMom][NBinsRecAngle], double Efficiency[NBinsTrueMom][NBinsTrueAngle],double *NumberOfPOT){

  double check=0;
  ifstream fEvent;

  cout<<"opening file "<<fDataName<<" for data input"<<endl;
  fEvent.open(fDataName,ios::in);
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      fEvent>>DataReconstructedEvents[e0][e1];
#ifdef DEBUG
      cout<<DataReconstructedEvents[e0][e1]<<" ";
#endif
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading data file"<<endl;
  fEvent.close();

  cout<<"opening file "<<fMCName<<" for mc input"<<endl;
  fEvent.open(fMCName,ios::in);
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      fEvent>>MCReconstructedEvents[e0][e1];
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading mc file, reconstructed events"<<endl;

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  fEvent>>MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	}
      }
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading mc file, true signal"<<endl;
  
    for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	fEvent>>MCReconstructedBkgEvents[e0][e1];
      }
    }
    fEvent>>check;
    if(check!=666) cout<<"problem in reading mc file, background"<<endl;
    
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	fEvent>>Efficiency[c0][c1];
      }
    } 
    fEvent>>check;
    if(check!=666) cout<<"problem in reading mc file, efficiency"<<endl;
    fEvent>>*NumberOfPOT;
    fEvent>>check;
    cout<<"number="<<*NumberOfPOT<<", check="<<check<<endl;
    if(check!=666) cout<<"problem in reading mc file, Number of POT"<<endl;
    fEvent.close();    
}

void Xsec::LoadNeutrinoFlux(TH1D * NeutrinoFlux){
  TFile * f = new TFile("files/nd34_tuned_11bv3.1_250ka.root");
  NeutrinoFlux = (TH1D*) f->Get("ing3_tune_numu");
  NeutrinoFlux->SetDirectory(0);//Detach the histogram from the TFile (one can close the file)
  f->Close();
}
			    
//First, build the likelihood matrix (tensor) based on the MC distribution
//Second, build the prior coming from the MC distribution. Whether to use this prior or not is defined in SetPrior function. Give also the true distribution (ptiot mc not normalised)
void Xsec::BuildLikelihood(double vLikelihood[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle], double vPriorMC[NBinsTrueMom][NBinsTrueAngle], double TrueEventsDistribution[NBinsTrueMom][NBinsTrueAngle], double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle]){
  double MCEvents[NBinsTrueMom][NBinsTrueAngle];
  double MCEventsTotal=0;

  //Prepare the likelihood distribution
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      MCEvents[c0][c1]=0;
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  vLikelihood[c0][c1][e0][e1]=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	  MCEvents[c0][c1]+=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	  MCEventsTotal+=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	}
      }

      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  if(MCEvents[c0][c1]!=0) vLikelihood[c0][c1][e0][e1]/=MCEvents[c0][c1];//P(Ej/Ci)=P(Ej&Ci)/P(Ci)
	}
      }
    }
  }

  //Prepare MC prior
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      vPriorMC[c0][c1]=MCEvents[c0][c1];
      if(MCEventsTotal!=0) vPriorMC[c0][c1]/=MCEventsTotal;
      TrueEventsDistribution[c0][c1]=MCEvents[c0][c1];
    }
  }
}

//First, choose the prior according to the user will: 
// If the number of iteration is 0, the user can choose whether to give it an original prior (by vInitialPrior and PriorMC=false) or to take the MC distribution as a prior (vInitialPriorMC and PriorMC=true).
// else, prior(iteration=N)=Posterior(iteration=N-1)
//Second, the prior is normalised.
void Xsec::SetPrior(double vPriorNormalised[NBinsTrueMom][NBinsTrueAngle], double vPrior[NBinsTrueMom][NBinsTrueAngle], double vInitialPriorMC[NBinsTrueMom][NBinsTrueAngle], double vInitialPrior[NBinsTrueMom][NBinsTrueAngle],double vPosterior[NBinsTrueMom][NBinsTrueAngle], bool PriorMC, double IterationStep){

  if(IterationStep==0){
    if(PriorMC){
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  vPrior[c0][c1]=vInitialPriorMC[c0][c1];
	}
      }
    }
    else{
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  vPrior[c0][c1]=vInitialPrior[c0][c1];
	}
      }
    }
  }
  else{
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	vPrior[c0][c1]=vPosterior[c0][c1];
      }
    }
  }
  
  double PriorEventsTotal=0;
  
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      PriorEventsTotal+=vPrior[c0][c1];
    }
  }
  
  if(PriorEventsTotal!=0){
    for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	vPriorNormalised[c0][c1]=vPrior[c0][c1]/PriorEventsTotal;
      }
    }
  }
}

//Build the unfolding "matrix" (tensor) using the likelihood and the prior chosen
void Xsec::BuildUnfolding(double vUnfolding[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle], double vLikelihood[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle], double vPrior[NBinsTrueMom][NBinsTrueAngle]){

  double MCEventsProbability[NBinsRecMom][NBinsRecAngle];

  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      MCEventsProbability[e0][e1]=0;
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  MCEventsProbability[e0][e1]+=vLikelihood[c0][c1][e0][e1]*vPrior[c0][c1];//Sum_cause(P(Ej|Ci)P(Ci))
	  vUnfolding[c0][c1][e0][e1]=vLikelihood[c0][c1][e0][e1]*vPrior[c0][c1];
       	}
      }
      for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
	for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
	  if(MCEventsProbability[e0][e1]!=0) vUnfolding[c0][c1][e0][e1]/=MCEventsProbability[e0][e1];//P(Ci/Ej)=P(Ej|Ci)P(Ci)/Sum_cause(P(Ej|Ci)P(Ci))
	}
      }
    }
  }

}
//Unfold the signal reconstructed distribution: data-mc
void Xsec::ApplyUnfolding(double vPosteriorEvents[NBinsTrueMom][NBinsTrueAngle], double vUnfolding[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle], double DataReconstructedEvents[NBinsRecMom][NBinsRecAngle], double MCReconstructedBkgEvents[NBinsRecMom][NBinsRecAngle]){

  //Build signal =data - background
  double DataReconstructedEvents_Signal[NBinsRecMom][NBinsRecAngle];
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      DataReconstructedEvents_Signal[e0][e1]=DataReconstructedEvents[e0][e1]-MCReconstructedBkgEvents[e0][e1];
      if(DataReconstructedEvents_Signal[e0][e1]<0) DataReconstructedEvents_Signal[e0][e1]=0;
      //DataReconstructedEvents_Signal[e0][e1]=DataReconstructedEvents[e0][e1];
    }
  }
#ifdef DEBUG

  cout<<"Test of the signal after bkg corection:"<<endl;
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      cout<<DataReconstructedEvents_Signal[e0][e1]<<" ";
    }
  }
  cout<<endl;

    cout<<"Re-Check the unfolding:"<<endl;
   for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
  cout<<vUnfolding[c0][c1][e0][e1]<<" ";
	}
      }
      cout<<endl;
    }
   }
#endif
  
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      vPosteriorEvents[c0][c1]=0;
	}
      }

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  vPosteriorEvents[c0][c1]+=vUnfolding[c0][c1][e0][e1]*DataReconstructedEvents_Signal[e0][e1];
	}
      }
      if(vPosteriorEvents[c0][c1]<0) vPosteriorEvents[c0][c1]=0;
    }
  }
}

//Generate statistical fluctuation in the DataReconstructedDistribution assuming Poissonian distribution in each reconstructed bin. If needed,
void Xsec::GenerateStatisticalFluctuations(double DataReconstructedEvents[NBinsRecMom][NBinsRecAngle]){

  TRandom3 * r = new TRandom3();
  r->SetSeed(0);

  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      double Mean=DataReconstructedEvents[e0][e1];
      DataReconstructedEvents[e0][e1]=r->PoissonD(Mean);//return a double and not an integer (number of neutrino detected is reweighted, and so not integer)
      //DataReconstructedEvents[e0][e1]=Mean+0.1*(Mean-r->PoissonD(Mean));//return a double and not an integer (number of neutrino detected is reweighted, and so not integer)
      //cout<<"TEMPORARY, bin="<<e0<<","<<e1<<", Mean="<<Mean<<", fluctuated="<<DataReconstructedEvents[e0][e1]<<endl; 
   }
  }

  delete r;
}

//Generate fluctuations of the MC in each bin that will affect the prior mc, the likelihood (and though the unfolding) and the background estimation (+the efficiency). The user should give a tensor of relative sigma errors one wants.
void Xsec::GenerateMCFluctuations(double MCReconstructedEvents_TrueSignal[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle],double RelativeSigma[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle]){

 TRandom3 * r = new TRandom3();
  r->SetSeed(0);

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1

	  double Mean=MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	  double Sigma=Mean*RelativeSigma[c0][c1][e0][e1];
	  MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]=r->Gaus(Mean,Sigma);
	  //if(Mean!=0) cout<<"TEMPORARY, bin="<<e0<<","<<e1<<", Mean="<<Mean<<", fluctuated="<<MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]<<endl;
	}
      }
    }
  }

  delete r;
}
			       
#endif

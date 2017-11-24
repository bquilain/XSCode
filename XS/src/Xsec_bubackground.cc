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

Xsec::Xsec(bool PM){
  _isPM=PM;
  Reco->SetDetector(PM);
}
Xsec::~Xsec(){}

void Xsec::SetDetector(bool PM){
  _isPM=PM;
  Reco->SetDetector(PM);
}
bool Xsec::GetDetector(){
  return _isPM;
}

void Xsec::Initialize(){

  InitializeGlobal(_isPM);
}

void Xsec::DetermineNuType(bool&IsSand,bool&IsAnti,bool&IsNuE,bool&IsBkgH,bool&IsBkgV,bool &IsSciBkg,int nutype, int intmode,double* VertexPosition){
  IsSand=0;
  IsAnti=0;
  IsNuE=0;
  IsBkgH=0;
  IsBkgV=0;
  IsSciBkg=0;
  if((_isPM && intmode==16) || (!_isPM && intmode==17)){
    //if(nutype==1) return 0;
    if(nutype==2) IsAnti=1;
    else if(nutype==3) IsNuE=1;
    else if(nutype==1){
      // let's determine if the vertex is in scintillator bars

      if(!_isPM && intmode==17){
	double x=VertexPosition[0];//cm
	double y=VertexPosition[1];//cm
	double z=VertexPosition[2]+120;//cm
	
	double xch,ych,zch;
	bool foundXY=false, foundZ=false;
	
	for(int ipln=0;ipln<PLNMAX;ipln++){
	  
	  foundZ=false;
	  int grid;
	  
	  for(int iview=0;iview<VIEWMAX;iview++){
	    IngD->get_pos_loli(15,iview,ipln,0,&xch,&ych,&zch); // scinti planes
	    if(fabs(z-zch)<loli_scinti_thick/2) {foundZ=true;grid=0;}
	    
	    IngD->get_pos_loli(15,iview,ipln,40,&xch,&ych,&zch); // grid 1
	    if(fabs(z-zch)<loli_scinti_width/2) {foundZ=true;grid=1;}
	    
	    IngD->get_pos_loli(15,iview,ipln,60,&xch,&ych,&zch); // grid 2
	    if(fabs(z-zch)<loli_scinti_width/2) {foundZ=true;grid=2;}
	  }
	  
	  if(foundZ){
	    if(grid==0){ // plane
	      bool foundX=false,foundY=false;
	      for (int ich=0;ich<40;ich++){
		IngD->get_pos_loli(15,0,ipln,ich,&xch,&ych,&zch);
		if(fabs(y-ych)<loli_scinti_width/2 && fabs(z-zch)<loli_scinti_thick/2) {foundY=true;}
		
		IngD->get_pos_loli(15,1,ipln,ich,&xch,&ych,&zch);
		if(fabs(x-xch)<loli_scinti_width/2 && fabs(z-zch)<loli_scinti_thick/2) {foundX=true;}	
	      }
	    if(foundX || foundY) IsSciBkg=true; // both can't be simultaneously satisfied
	    }
	    else {
	      bool foundX=false,foundY=false;
	      for (int ich=0;ich<40;ich++){
		IngD->get_pos_loli(15,0,ipln,ich,&xch,&ych,&zch);
		if(fabs(y-ych)<loli_scinti_thick/2 && fabs(z-zch)<loli_scinti_width/2) {foundY=true;}
	      }
	      for (int ich=0;ich<CHMAX;ich++){
		IngD->get_pos_loli(15,1,ipln,ich,&xch,&ych,&zch);
		if(fabs(x-xch)<loli_scinti_thick/2 && fabs(z-zch)<loli_scinti_width/2) {foundX=true;}
	      }
	      if(foundX || foundY) IsSciBkg=true; 
	    }
	    
	    break;
	  }
	}
      
      if (fabs(x)>50 || fabs(y)>50) IsSciBkg=false;
      }
    }
  }
  else if(intmode>=0 && intmode<7) IsBkgH=1;
  else if(intmode>=7 && intmode<14) IsBkgV=1;
  else IsSand=1;
}


int Xsec::DetermineFSI(int IsSand,int IsAnti,int IsNuE,int IsBkgH,int IsBkgV,int IsSciBkg,IngridEventSummary * evt){
  int FSIMuons=0;int FSIPions=0;int FSIProtons=0;int FSIOther=0;int FSINeutralPions=0;
  int FSIInt;
  IngridSimParticleSummary * SimPart;

  for(int is=0;is<evt->NIngridSimParticles();is++){
    SimPart=(IngridSimParticleSummary*) evt->GetSimParticle(is);
    //cout<<"Particle type="<<SimPart->pdg<<", NRJ="<<SimPart->momentum[0]<<", "<<SimPart->momentum[3]<<endl;

    if(TMath::Abs(SimPart->pdg)==211) FSIPions++;
    else if(TMath::Abs(SimPart->pdg)==111) FSINeutralPions++;
    else if(TMath::Abs(SimPart->pdg)==2212) FSIProtons++;
    else if(TMath::Abs(SimPart->pdg)==13) FSIMuons++;
    else if(TMath::Abs(SimPart->pdg)!=2112 && TMath::Abs(SimPart->pdg)!=22) FSIOther++;
  }

  if(FSIMuons==1){
    if(FSIOther==0){
      if(FSINeutralPions==0){
	if(FSIPions==0){
	  if(FSIProtons<2) FSIInt=FSIProtons;
	  else FSIInt=2;
	}
	else if(FSIPions==1) FSIInt=3;
	else FSIInt=5; // CCNpi+/-, N>1 are inside CCother
	//else FSIInt=4;
      }
      else if(FSINeutralPions==1 && FSIPions==0) FSIInt=4; // CC1pi0
      else FSIInt=5; // CCNpiMpi0, N,M>1 are inside CCother
    }
    else FSIInt=5;
  }
  else FSIInt=6;

  if(IsSand) FSIInt=7;
  else if(IsBkgH) FSIInt=8;
  else if(IsBkgV) FSIInt=9;
  else if(IsAnti) FSIInt=10;
  else if(IsNuE) FSIInt=11;
  else if(IsSciBkg) FSIInt=12;

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
void Xsec::LoadInputFiles_OnlyUnfoldedData(char * InputName, double ** DataUnfoldedEvents, double ** DataTrueEvents){
  
  TChain * wtree = new TChain("wtree");
  wtree->Add(InputName);
  cout<<"opening file "<<InputName<<" for data input"<<endl;
  
  int Iterations, Systematics, Statistics;
  double Events[NBinsTrueMom][NBinsTrueAngle];
  double EventsAll[NBinsTrueMom][NBinsTrueAngle][NBinsRecMom][NBinsRecAngle];
  double TrueEvents[NBinsTrueMom][NBinsTrueAngle];
  double XSection[NBinsTrueMom][NBinsTrueAngle];

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

  TBranch* Br_XSection = wtree->GetBranch("XSection");
  Br_XSection->SetAddress(XSection);
  wtree->SetBranchAddress("XSection",XSection);
    
  int nevt=wtree->GetEntries();
  
  for(int ievt=0;ievt<nevt;ievt++){
    wtree->GetEntry(ievt);
    
    for(int c0=0;c0<NBinsTrueMom;c0++){
      for(int c1=0;c1<NBinsTrueAngle;c1++){
	DataUnfoldedEvents[c0][c1]=Events[c0][c1];
	//DataUnfoldedEvents[c0][c1]=XSection[c0][c1];
#ifdef DEBUG
	cout<<"True bin ("<<c0<<","<<c1<<")     XS="<<DataUnfoldedEvents[c0][c1]<<", Nev="<<Events[c0][c1]<<endl;
#endif
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
void Xsec::LoadInputFiles_OnlySelectedData(char * fDataName,double ** DataReconstructedEvents){

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

//For side band: Read the output file of CC0pi selection, but only the selected data
void Xsec::LoadInputFiles_OnlySelectedDataSB(char * fDataName,char * fDataNameSB,double ** DataReconstructedEvents){

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


  //Side band
  cout<<"opening file "<<fDataNameSB<<" for data input"<<endl;
  fEvent.open(fDataNameSB,ios::in);
  for(int e0=0;e0<NBinsRecMomSB;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngleSB;e1++){//loop over effect 1
      fEvent>>DataReconstructedEvents[e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal];
#ifdef DEBUG
      cout<<DataReconstructedEvents[e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal]<<" ";
#endif
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading data file"<<endl;
  fEvent.close();
}

//Read the output file of CC0pi selection.
void Xsec::LoadInputFiles(char * fDataName,char * fMCName, double **** MCReconstructedEvents_TrueSignal, double ** DataReconstructedEvents,double ** MCReconstructedEvents, double ** MCReconstructedBkgEvents, double ** Efficiency,double *NumberOfPOT){
    //////////
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

//Read the output file of CC0pi selection, with side band.
void Xsec::LoadInputFilesSB(char * fDataName,char * fMCName, char * fDataNameSB,char * fMCNameSB, double **** MCReconstructedEvents_TrueSignal, double ** DataReconstructedEvents,double ** MCReconstructedEvents, double ** MCReconstructedBkgEvents, double ** Efficiency,double *NumberOfPOT){

  //First, fill the matrix with 0 since it will not be filled everywhere (matrix is diagonal by block):
  for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
      DataReconstructedEvents[e0][e1]=0;
      MCReconstructedEvents[e0][e1]=0;
      MCReconstructedBkgEvents[e0][e1]=0;
    }
  }    
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  MCReconstructedEvents_TrueSignal[c0][c1][e0][e1]=0;
	}
      }
    }
  }
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      Efficiency[c0][c1]=0;
    }
  }

  
  //////////////////FILL WITH SIGNAL
  double check=0;
  ifstream fEvent;

  cout<<"opening file "<<fDataName<<" for data input"<<endl;
  fEvent.open(fDataName,ios::in);
  for(int e0=0;e0<NBinsRecMomSignal;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngleSignal;e1++){//loop over effect 1
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
  for(int e0=0;e0<NBinsRecMomSignal;e0++){//loop over effect 0
    for(int e1=0;e1<NBinsRecAngleSignal;e1++){//loop over effect 1
      fEvent>>MCReconstructedEvents[e0][e1];
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading mc file, reconstructed events"<<endl;

  for(int c0=0;c0<NBinsTrueMomSignal;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngleSignal;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMomSignal;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngleSignal;e1++){//loop over effect 1
	  fEvent>>MCReconstructedEvents_TrueSignal[c0][c1][e0][e1];
	}
      }
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading mc file, true signal"<<endl;
  
    for(int e0=0;e0<NBinsRecMomSignal;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngleSignal;e1++){//loop over effect 1
	fEvent>>MCReconstructedBkgEvents[e0][e1];
      }
    }
    fEvent>>check;
    if(check!=666) cout<<"problem in reading mc file, background"<<endl;
    
    for(int c0=0;c0<NBinsTrueMomSignal;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngleSignal;c1++){//loop over cause 1
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

    
    ////////////////////////SIDE BAND////////////////////////////////////////
    cout<<"Side Band: opening file "<<fDataNameSB<<" for data input"<<endl;
    fEvent.open(fDataNameSB,ios::in);
        
    for(int e0=0;e0<NBinsRecMomSB;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngleSB;e1++){//loop over effect 1
	fEvent>>DataReconstructedEvents[e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal];
#ifdef DEBUG
	cout<<DataReconstructedEvents[e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal]<<" ";
#endif
      }
    }
    fEvent>>check;
    if(check!=666) cout<<"problem in reading data file"<<endl;
    fEvent.close();
    
    cout<<"Side Band: opening file "<<fMCNameSB<<" for mc input"<<endl;
    fEvent.open(fMCNameSB,ios::in);
    for(int e0=0;e0<NBinsRecMomSB;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngleSB;e1++){//loop over effect 1
	fEvent>>MCReconstructedEvents[e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal];
      }
    }
    fEvent>>check;
    if(check!=666) cout<<"problem in reading mc file, reconstructed events"<<endl;
    
    for(int c0=0;c0<NBinsTrueMomSB;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngleSB;c1++){//loop over cause 1
      for(int e0=0;e0<NBinsRecMomSB;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngleSB;e1++){//loop over effect 1
	  fEvent>>MCReconstructedEvents_TrueSignal[c0+NBinsTrueMomSignal][c1+NBinsTrueAngleSignal][e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal];
	}
      }
    }
  }
  fEvent>>check;
  if(check!=666) cout<<"problem in reading mc file, true signal"<<endl;
  
    for(int e0=0;e0<NBinsRecMomSB;e0++){//loop over effect 0
      for(int e1=0;e1<NBinsRecAngleSB;e1++){//loop over effect 1
	fEvent>>MCReconstructedBkgEvents[e0+NBinsRecMomSignal][e1+NBinsRecAngleSignal];
      }
    }
    fEvent>>check;
    if(check!=666) cout<<"problem in reading mc file, background"<<endl;
    
    for(int c0=0;c0<NBinsTrueMomSB;c0++){//loop over cause 0
      for(int c1=0;c1<NBinsTrueAngleSB;c1++){//loop over cause 1
	fEvent>>Efficiency[c0+NBinsTrueMomSignal][c1+NBinsTrueAngleSignal];
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
void Xsec::BuildLikelihood(double **** vLikelihood, double ** vPriorMC, double ** TrueEventsDistribution, double **** MCReconstructedEvents_TrueSignal){
  double ** MCEvents = new double*[NBinsTrueMom];
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    MCEvents[c0] = new double [NBinsTrueAngle];
  }
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

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    delete MCEvents[c0];
  }
  delete MCEvents;
  
}

//First, choose the prior according to the user will: 
// If the number of iteration is 0, the user can choose whether to give it an original prior (by vInitialPrior and PriorMC=false) or to take the MC distribution as a prior (vInitialPriorMC and PriorMC=true).
// else, prior(iteration=N)=Posterior(iteration=N-1)
//Second, the prior is normalised.
void Xsec::SetPrior(double ** vPriorNormalised, double ** vPrior, double ** vInitialPriorMC, double ** vInitialPrior,double ** vPosterior, bool PriorMC, double IterationStep){

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
void Xsec::BuildUnfolding(double **** vUnfolding, double **** vLikelihood, double ** vPrior){

  double ** MCEventsProbability = new double*[NBinsRecMom];
  for(int c0=0;c0<NBinsRecMom;c0++){//loop over cause 0
    MCEventsProbability[c0] = new double [NBinsRecAngle];
  }
  
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

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    delete MCEventsProbability[c0];
  }
  delete MCEventsProbability;

}
//Unfold the signal reconstructed distribution: data-mc
void Xsec::ApplyUnfolding(double ** vPosteriorEvents, double **** vUnfolding, double ** DataReconstructedEvents, double ** MCReconstructedBkgEvents){

  //Build signal =data - background
  double ** DataReconstructedEvents_Signal = new double*[NBinsRecMom];
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    DataReconstructedEvents_Signal[c0] = new double [NBinsRecAngle];
  }

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

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    delete DataReconstructedEvents_Signal[c0];
  }
  delete DataReconstructedEvents_Signal;

}

//Generate statistical fluctuation in the DataReconstructedDistribution assuming Poissonian distribution in each reconstructed bin. If needed,
void Xsec::GenerateStatisticalFluctuations(double ** DataReconstructedEvents){

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
void Xsec::GenerateMCFluctuations(double **** MCReconstructedEvents_TrueSignal,double **** RelativeSigma){

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

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
#include <TSpline.h>
//int LimitTracks=20;
int LimitRecs=5;
//int NDials=175;
#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "INGRID_Dimension.cc"

#include "setup.h"
#include "Hit.h"
#include "PMdispRev.h"
#include "Corrections.cc"
#include "Reconstruction.cc"
#include "Xsec.cc"
int NBins=600;
int StartPE=0;
int EndPE=300;

//#define TRAINGOODTRACKS
#define DEBUG
// you may choose ONE of the following refinements:
// (without any, you build a muCL as muon vs others)
//   - PI_LIKELIHOOD is to build a pion vs proton CL for non-mu-like particles
//   - MUPI_LIKELIHOOD is to build the muCL as muon+pion vs proton

#define PI_LIKELIHOOD
#define VSPROTON
//#define MUPI_LIKELIHOOD

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

/* already done in Reconstruction.cc
double DegRad(double angle){
  return angle*TMath::Pi()/180.;
}

double RadDeg(double angle){
  return angle*180./TMath::Pi();
}
*/

int main(int argc, char **argv)
{
  InitializeGlobal();

  bool Disp=false;
  int FileNumber=1;
  int c=-1;
  bool MC=false;
  bool SystematicPE=false;
  int RandomIteration;
  string RandomIteration_string;
  gErrorIgnoreLevel = kError;
  bool PM=true;

  while ((c = getopt(argc, argv, "f:d:mwr:")) != -1) {
    switch(c){
    case 'f':
      FileNumber=atoi(optarg);
      break;
    case 'd':
      Disp=true;
      break;
    case 'm':
      MC=true;
      break;
    case 'w':
      PM=false;
      break;
    case 'r':
      SystematicPE=true;
      RandomIteration=atoi(optarg);
      RandomIteration_string=optarg;
      break;
    }
  }

  cout<<"Welcome"<<endl;
  vector <HitTemp> HitV;
  TApplication theApp("App",0,0);


  INGRID_Dimension * IngDim = new INGRID_Dimension();
  Reconstruction * Rec = new Reconstruction(PM);
  Corrections * Cor = new Corrections(PM);  
  Xsec * xs = new Xsec(PM);
  xs->Xsec::Initialize();

  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));  
  char suffix[3];sprintf(suffix,(PM?"":"_WM"));  
  char * cMCOUT = getenv((PM?"MCOUTPUTSTORAGE":"MCOUTPUTSTORAGE_WM"));

  int type;//0=CC, 1=NC
  double Nu_E;
  char File[256];
  int nTracks=-1;
  float weight=1;
  int NIngBasRec;


  char *NameDist = new char[256];
  char * FileDist = new char[256];

  
#ifdef MUPI_LIKELIHOOD
  TFile * wfile = new TFile(Form("PID/PDFMuPiCL_Likelihood%d%s.root",FileNumber,suffix),"recreate");
#else
  //TFile * wfile = new TFile(Form("src/PDFMuCL_Likelihood%s.root",suffix),"recreate");
  TFile * wfile = new TFile(Form("PID/PDFMuCL_Likelihood%d%s.root",FileNumber,suffix),"recreate");
#endif
  TH1D * PDFMuCL_PMIng_Muon = new TH1D("PDFMuCL_PMIng_Muon","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_PMSci_Muon = new TH1D("PDFMuCL_PMSci_Muon","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_Ing_Muon = new TH1D("PDFMuCL_Ing_Muon","",NBins,StartPE,EndPE);
  PDFMuCL_PMIng_Muon->Sumw2();PDFMuCL_PMSci_Muon->Sumw2();PDFMuCL_Ing_Muon->Sumw2();

  TH1D * PDFMuCL_PMIng_NotMuon = new TH1D("PDFMuCL_PMIng_NotMuon","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_PMSci_NotMuon = new TH1D("PDFMuCL_PMSci_NotMuon","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_Ing_NotMuon = new TH1D("PDFMuCL_Ing_NotMuon","",NBins,StartPE,EndPE);
  PDFMuCL_PMIng_NotMuon->Sumw2();PDFMuCL_PMSci_NotMuon->Sumw2();PDFMuCL_Ing_NotMuon->Sumw2();

#ifdef PI_LIKELIHOOD
  TH1D * PDFPiCL_PMIng_Pion = new TH1D("PDFPiCL_PMIng_Pion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_PMSci_Pion = new TH1D("PDFPiCL_PMSci_Pion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_Ing_Pion = new TH1D("PDFPiCL_Ing_Pion","",NBins,StartPE,EndPE);

  TH1D * PDFPiCL_PMIng_NotPion = new TH1D("PDFPiCL_PMIng_NotPion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_PMSci_NotPion = new TH1D("PDFPiCL_PMSci_NotPion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_Ing_NotPion = new TH1D("PDFPiCL_Ing_NotPion","",NBins,StartPE,EndPE);
  PDFPiCL_PMIng_Pion->Sumw2();PDFPiCL_PMSci_Pion->Sumw2();PDFPiCL_Ing_Pion->Sumw2();
  PDFPiCL_PMIng_NotPion->Sumw2();PDFPiCL_PMSci_NotPion->Sumw2();PDFPiCL_Ing_NotPion->Sumw2();

  TH1D * PDFPiCL_WM_Pion[3], *PDFPiCL_WM_NotPion[3];
#endif
  TH1D * PDFMuCL_WM_Muon[3], *PDFMuCL_WM_NotMuon[3];

  for (int i =0;i<3;i++){ // normal angle slices
    PDFMuCL_WM_Muon[i]=new TH1D(Form("PDFMuCL_WM_Muon_%d",i),"",100,0,100);
    PDFMuCL_WM_NotMuon[i]=new TH1D(Form("PDFMuCL_WM_NotMuon_%d",i),"",100,0,100);
    PDFMuCL_WM_Muon[i]->Sumw2(); PDFMuCL_WM_NotMuon[i]->Sumw2();
#ifdef PI_LIKELIHOOD
    PDFPiCL_WM_Pion[i]=new TH1D(Form("PDFPiCL_WM_Pion_%d",i),"",100,0,100);
    PDFPiCL_WM_NotPion[i]=new TH1D(Form("PDFPiCL_WM_NotPion_%d",i),"",100,0,100);
    PDFPiCL_WM_Pion[i]->Sumw2(); PDFPiCL_WM_NotPion[i]->Sumw2();
#endif
  }

  TH1D * PDFParticle = new TH1D("PDFParticle","",2e4,-1e4,1e4);

  double PECorrected=0;
#ifdef TRAINGOODTRACKS
  const int NParticles = 4;
  TH1D * ProportionGoodHits[NParticles];
  TH1D * AngleDifferenceX[NParticles];
  TH1D * AngleDifferenceY[NParticles];
  TH2D * AngleDifferenceXProportionGoodHits[NParticles];
  for(int ip=0;ip<NParticles;ip++){
    ProportionGoodHits[ip] = new TH1D(Form("ProportionGoodHits%d",ip),"",100,0,1);
    AngleDifferenceX[ip] = new TH1D(Form("AngleDifferenceX%d",ip),"",180,-90,90);
    AngleDifferenceY[ip] = new TH1D(Form("AngleDifferenceY%d",ip),"",180,-90,90);
    AngleDifferenceXProportionGoodHits[ip] = new TH2D(Form("AngleDifferenceXProportionGoodHits%d",ip),"",180,-90,90,100,0,1);
  }
#endif
  cout<<"Opening Events"<<endl;

  IngridEventSummary* evt;
  IngridSimVertexSummary * simver = new IngridSimVertexSummary();//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
  IngridHitSummary * Hit = new IngridHitSummary();
  IngridSimHitSummary * SimHit= new IngridSimHitSummary();
  IngridSimHitSummary * SimHit2= new IngridSimHitSummary();
  IngridSimParticleSummary * SimPart= new IngridSimParticleSummary();
  BeamInfoSummary * BeamSummary = new BeamInfoSummary();
  PMAnaSummary * recon = new PMAnaSummary();	
  IngridHitSummary * hit = new IngridHitSummary();
  IngridHitSummary * hit2 =new IngridHitSummary();

  bool IsFV=false;
  bool IsSand=false;
  bool IsAnti=false;
  bool IsNuE=false;
  bool IsBkgH=false;
  bool IsBkgV=false;
  bool IsSciBkg=false;
  
  //TFile * _file0;
  TTree * tree;
  
  cout<<"opening"<<endl;
  sprintf(File,"%s/%sMC_Run1_%d_wNoise_ana_shrinked.root",cMCOUT,DetName,FileNumber);
 
  if(!MC) return 0;
    TFile * _file0 = new TFile(File);
    if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
    else return 0;
    tree=(TTree*) _file0->Get("tree");
    if(tree!=tree) return 0;
    int nevt=(int) tree->GetEntries();
    cout<<"Total Number Of Events="<<nevt<<endl;

    evt = new IngridEventSummary();
    tree->SetBranchAddress("fDefaultReco.",&evt);
 
    for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
      if((ievt%1000)==0) cout<<"Processing "<<ievt<<endl;
      evt->Clear();
      tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
      simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
      int mod;
      double norm;
      double totcrsne;

      if(MC){
	mod=simver->mod;
	double TrueVertexPosition[3] = {simver->xnu,simver->ynu,simver->znu};
	norm=simver->norm;
	totcrsne=simver->totcrsne;
	if(IsSand || IsAnti || IsBkgH || IsBkgV || IsNuE) continue;
	xs->Xsec::DetermineNuType(IsSand,IsAnti,IsNuE,IsBkgH,IsBkgV,IsSciBkg,simver->nutype,simver->mod,TrueVertexPosition);
	
	weight = 1;
	if(IsBkgH || IsBkgV) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(7.87*58.5+22);// ML small contribution from scinti added                           
	else if(IsSand) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*300); // ML 2017/11/27 in fact the vertex is drawn only in a 3 meters shell (not 4.7m)
	else  if(PM) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*46.2;
	else weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*50.;//for WM - ML 2017/05/05 the water tank is
	// ML 2017/08/04 the 1/100 corresponds to the file# 
	// it's a normalization factor so idc
      }

      NIngBasRec= evt->NPMAnas();

      for(int irec=0;irec<NIngBasRec;irec++){
	recon = (PMAnaSummary*) evt->GetPMAna(irec);	
	nTracks=recon->Ntrack;
      
	if(!(recon->vetowtracking) && !(recon->edgewtracking)) {

	  vector< vector<Hit3D> > VecDouble;
	  for(int i=0;i<VecDouble.size();i++){
	    VecDouble[i].clear();
	  }
	  for(int itrk=0;itrk<nTracks;itrk++){
	    vector <Hit3D> VecT;
	    VecT.clear();
	    HitV.clear();
	    HitV=Rec->Reconstruction::EraseDoubleHitsPM(recon,itrk,HitV);
	    VecT=Rec->Reconstruction::Hit2DMatchingPM(evt,recon,HitV,VecT,MC);
	    VecDouble.push_back(VecT);
	  }

	  
	  for(int itrk=0;itrk<nTracks;itrk++){//loop on track
	    bool SeveralHits=false;
	    HitV.clear(); 
	    vector <Hit3D> Vec;
	    HitV=Rec->Reconstruction::EraseDoubleHitsPM(recon,itrk,HitV);
	    Vec=Rec->Reconstruction::Hit2DMatchingPM(evt,recon,HitV,Vec,MC);
	    Vec=Rec->CountSharedHits(Vec,VecDouble,itrk);//fill used
	    //VecDoubleNeeded
	  
	    if(Vec.size()==0){
	      cout<<"Stopped because no Hits are corresponding"<<endl;
	      continue;
	    }
	    
	    double dx=1./TMath::Cos(DegRad((recon->angle)[itrk]));
	    Vec=Cor->Corrections::GetFiberAttenuation(Vec);
	    Vec=Cor->Corrections::GetDXCorrection(Vec,dx);
	    
	    //Track length:
	    double TrkLength= Rec->Reconstruction::GetTrackLength(Vec);
	    int ParticleNumber=TMath::Abs(Rec->Reconstruction::GetTrackParticleBugged(evt, recon, itrk,TrkLength));
	    int PDG=((IngridSimParticleSummary*) evt->GetSimParticle(ParticleNumber))->pdg;
#ifdef TRAINGOODTRACKS
	    SimPart = (IngridSimParticleSummary*) evt->GetSimParticle(ParticleNumber);
	    double thetaX=(TMath::ATan((SimPart->fpos[0]-SimPart->ipos[0])/(SimPart->fpos[2]-SimPart->ipos[2])))*180./TMath::Pi();
	    double thetaY=(TMath::ATan((SimPart->fpos[1]-SimPart->ipos[1])/(SimPart->fpos[2]-SimPart->ipos[2])))*180./TMath::Pi();
	    
	    int PSimplified = 0;	    
	    if(TMath::Abs(PDG) == 13) PSimplified = 0;
	    else if(TMath::Abs(PDG) == 211) PSimplified = 1;
	    else if(TMath::Abs(PDG) == 2212) PSimplified = 2;
	    else PSimplified = 3;
	    //cout<<recon->thetay[itrk]<<", "<<thetaX<<", "<<GetHitProportion(PDG,Vec)<<endl;
	    AngleDifferenceX[PSimplified]->Fill(recon->thetay[itrk] - thetaX,weight);
	    AngleDifferenceY[PSimplified]->Fill(recon->thetax[itrk] - thetaY,weight);
	    vector <int> PDGHits; PDGHits.clear();
	    PDGHits = Rec->TrackComposition(Vec);
	    double Proportions=0; double Total=0;
	    for(int p=0;p<NParticles;p++){
	      if(p == PSimplified) Proportions = PDGHits[p];
	      Total += PDGHits[p];
	    }
	    if(Total!=0) Proportions /= Total;
	    //cout<<"PDG: "<<PDG<<", Proportions: "<<Proportions<<", Type mu: "<<PDGHits[0]<<", Type pi: "<<PDGHits[1]<<", Type p: "<<PDGHits[2]<<endl;
	    ProportionGoodHits[PSimplified]->Fill(Proportions,weight);
	    AngleDifferenceXProportionGoodHits[PSimplified]->Fill(TMath::Abs(recon->thetay[itrk] - thetaX)+TMath::Abs(recon->thetax[itrk] - thetaY), Proportions,weight);
	    //bool TrackGoodForPID = false;
	    if(TMath::Abs(recon->thetay[itrk] - thetaX)+TMath::Abs(recon->thetax[itrk] - thetaY) > GoodTracksCut) continue;
#endif
	    PDFParticle->Fill(PDG,weight);
	    //cout<<PDG<<endl;
	    double PEPlane[3][NPlnPM][2];//one per type of scinti, plane and view
	    int Plane[3][NPlnPM][2];
	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int i=0;i<3;i++){
		for(int iv=0;iv<2;iv++){
		  PEPlane[i][ipln][iv]=0;
		  Plane[i][ipln][iv]=0;
		}
	      }
	    }
	    
	  //1. Count the number of hits that are "useable" for the CL (not too close from the vertex). The goal (later used) is that for events where the number of hits useable < 2 per view, we will use hits near the vertex
	  int NCLHits[NView]={0,0};//Number of hits in the vertex plane or the plane just after
	  int NCLHitsNonIsolated[NView]={0,0};//Number of hits in the vertex plane or the plane just after
	  for(int i=0;i<Vec.size();i++){

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
	    NCLHitsNonIsolated[Vec[i].view]++;
	    if(Vec[i].used<=1) NCLHits[Vec[i].view]++;
	  }
	  //cout<<"NHits:"<<Vec.size()<<", Non-isolated hits: "<<NCLHitsNonIsolated[0]+NCLHitsNonIsolated[1]<<", Isolated hits: "<<NCLHits[0]+NCLHits[1]<<endl;
	  int HitCut;
	  //if((NCLHitsNonIsolated[0] < CutNHitsPID) || (NCLHitsNonIsolated[1] < CutNHitsPID)) HitCut=0;
	  //else if((NCLHits[0] < CutNHitsPID) || (NCLHits[1] < CutNHitsPID)) HitCut=1;
	  //else HitCut=2;
	  if((NCLHitsNonIsolated[0]+NCLHitsNonIsolated[1]) < CutNHitsPID) HitCut=0;
	  else if((NCLHits[0]+NCLHits[1]) < CutNHitsPID) HitCut=1;
	  else HitCut=2;
	  //cout<<"PDG: "<<PDG<<", "<<HitCut<<endl;
	  //2. Real loop on the event:
#ifdef DEBUG2
	  cout<<"Number of hits="<<Vec.size()<<endl;
#endif
	  sort(Vec.begin(),Vec.end());
	  double ControlTrkLength = TMath::Sqrt(pow(Vec.back().x-Vec.front().x,2)+pow(Vec.back().y-Vec.front().y,2)+pow(Vec.back().z-Vec.front().z,2));

	  for(int i=0;i<Vec.size();i++){

	    //2.a Remove hits near the vertex
	    if(HitCut == 0){//Case where we have too few hits non isolated, we keep the second layer near vertex
	      if(Vec[i].mod==16){
		if(Vec[i].view==0){ if(Vec[i].pln==(recon->startxpln)[itrk]) continue;}
		else{ if(Vec[i].pln==(recon->startypln)[itrk]) continue;}
	      }
	      //cout<<"Not enough non-isolated hits: "<<endl;
	    }
	    else if(HitCut == 1){//Case where we have enough non isolated hits but not enough isolated: we use the non-isolated, but we cut the vertex that is too dangerous.
	      if(Vec[i].mod==16){
		if(Vec[i].view==0){ if((Vec[i].pln==(recon->startxpln)[itrk]) || (Vec[i].pln==((recon->startxpln)[itrk]+1)) ) continue;}
		else{ if((Vec[i].pln==(recon->startypln)[itrk]) || (Vec[i].pln==((recon->startypln)[itrk]+1)) ) continue;}
	      }
	      //cout<<"Not enough isolated hits"<<endl;
	    }
	    else{//Case where we have enough hits
	      if(Vec[i].used>1) continue;
	    }
	    
	    if(Vec[i].mod==16 && Rec->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
	      PECorrected=Vec[i].pecorr;
	      PEPlane[1][Vec[i].pln][Vec[i].view]+=PECorrected;
	      Plane[1][Vec[i].pln][Vec[i].view]++;
	    }
	    else if(Vec[i].mod==16 && !(Rec->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
	      PECorrected=Vec[i].pecorr;
	      PEPlane[0][Vec[i].pln][Vec[i].view]+=PECorrected;
	      Plane[0][Vec[i].pln][Vec[i].view]++;
	    }
	    else if(Vec[i].mod<14){
	      PECorrected=Vec[i].pecorr;
	      PEPlane[2][Vec[i].pln][Vec[i].view]+=PECorrected;
	      Plane[2][Vec[i].pln][Vec[i].view]++;
	    }
	  }//EndPE of loop over hits


	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int iview=0;iview<2;iview++){
#ifdef MUPI_LIKELIHOOOD
		if(PDG==13 || abs(PDG)==211){
#else 
		if(PDG==13){
#endif
		  if(Plane[0][ipln][iview]!=0 || Plane[1][ipln][iview]!=0){//case PM & active plane
		    if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){
		      PDFMuCL_PMSci_Muon->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
		    }
		    else PDFMuCL_PMIng_Muon->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more INGRID
		  }
		  else if(Plane[2][ipln][iview]!=0) PDFMuCL_Ing_Muon->Fill(PEPlane[2][ipln][iview],weight);
		}
		else {
		  if(Plane[0][ipln][iview]!=0 || Plane[1][ipln][iview]!=0){//case PM & active plane
		    if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){
		      PDFMuCL_PMSci_NotMuon->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#ifdef PI_LIKELIHOOD
		      if(abs(PDG)==211)
			PDFPiCL_PMSci_Pion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#ifdef VSPROTON
		      else if(abs(PDG) != 13)
			PDFPiCL_PMSci_NotPion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#else
		      else
			PDFPiCL_PMSci_NotPion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#endif
		      
#endif
		    }
		    else {
		      PDFMuCL_PMIng_NotMuon->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more INGRID
#ifdef PI_LIKELIHOOD
		      if(abs(PDG)==211)
			PDFPiCL_PMIng_Pion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#ifdef VSPROTON
		      else if(abs(PDG) != 13)
			PDFPiCL_PMIng_NotPion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#else
		      else
			PDFPiCL_PMIng_NotPion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#endif

#endif
		    }
		  }
		  else if(Plane[2][ipln][iview]!=0) {
		    PDFMuCL_Ing_NotMuon->Fill(PEPlane[2][ipln][iview],weight);		  
#ifdef PI_LIKELIHOOD	
		    if(abs(PDG)==211)
		      PDFPiCL_Ing_Pion->Fill(PEPlane[2][ipln][iview],weight);		  
#ifdef VSPROTON
		      else if(abs(PDG) != 13)
			PDFPiCL_Ing_NotPion->Fill(PEPlane[2][ipln][iview],weight);
#else
		    else
		      PDFPiCL_Ing_NotPion->Fill(PEPlane[2][ipln][iview],weight);
#endif
		    
#endif
		  }
		}
	      }
	    }
	    // end of PM+INGRID treatment

	    // ------------ now treat WM (ML 2017/08/04) ------
	    if(!PM){
	      Cor->GetDXCorrectionWM(Vec,DegRad(recon->angle[itrk]),DegRad(recon->thetax[itrk]),DegRad(recon->thetay[itrk]));
	      for(int i=0;i<Vec.size();i++){
		if(Vec[i].used>1) continue;

		bool grid=(Vec[i].ch>=40);
		double angle2D=(Vec[i].view==0?recon->thetax[itrk]:recon->thetay[itrk]);
		double normalAngle=RadDeg(Rec->GetNormalAngleRad(grid,DegRad(recon->angle[itrk]),DegRad(angle2D)));
		int normalAngleSlice=0;
		if(normalAngle>80) normalAngleSlice=2;
		else if(normalAngle>70) normalAngleSlice=1;
#ifdef MUPI_LIKELIHOOD
		if(PDG==13 || abs(PDG)==211)
#else
		if(PDG==13)
#endif
		  PDFMuCL_WM_Muon[normalAngleSlice]->Fill(Vec[i].pecorr,weight);
		else {
		  PDFMuCL_WM_NotMuon[normalAngleSlice]->Fill(Vec[i].pecorr,weight);
#ifdef PI_LIKELIHOOD
		  if(abs(PDG)==211)
		    PDFPiCL_WM_Pion[normalAngleSlice]->Fill(Vec[i].pecorr,weight);
#ifdef VSPROTON
		  else if(abs(PDG) != 13)
		    PDFPiCL_WM_NotPion[normalAngleSlice]->Fill(Vec[i].pecorr,weight);
#else
		  else
		    PDFPiCL_WM_NotPion[normalAngleSlice]->Fill(Vec[i].pecorr,weight);
#endif
		  
#endif
		}
	      }
	    }
	    // ---------- end of WM treatment ------------

	    Vec.clear();
	  }//Tracks
	}//cuts
      }//Evt
    }
    evt->Delete();
    _file0->Close();
    delete _file0;
  
  wfile->cd();

  PDFParticle->Write();

  if(PM){ 
    PDFMuCL_PMIng_Muon->Write();
    PDFMuCL_PMSci_Muon->Write();
    PDFMuCL_PMIng_NotMuon->Write();
    PDFMuCL_PMSci_NotMuon->Write();

#ifdef PI_LIKELIHOOD
    PDFPiCL_PMIng_Pion->Write();
    PDFPiCL_PMSci_Pion->Write();
    PDFPiCL_PMIng_NotPion->Write();
    PDFPiCL_PMSci_NotPion->Write();
#endif
  }
  else{

    for(int i=0;i<3;i++){
      PDFMuCL_WM_Muon[i]->Write();
      PDFMuCL_WM_NotMuon[i]->Write();

#ifdef PI_LIKELIHOOD
      PDFPiCL_WM_Pion[i]->Write();
      PDFPiCL_WM_NotPion[i]->Write();
#endif
    }
  }

  PDFMuCL_Ing_Muon->Write();
  PDFMuCL_Ing_NotMuon->Write();
#ifdef PI_LIKELIHOOD
  PDFPiCL_Ing_Pion->Write();
  PDFPiCL_Ing_NotPion->Write();
#endif

  
#ifdef TRAINGOODTRACKS
  for(int ip=0;ip<NParticles;ip++){
    ProportionGoodHits[ip]->Write();
    AngleDifferenceX[ip]->Write();
    AngleDifferenceY[ip]->Write();
    AngleDifferenceXProportionGoodHits[ip]->Write();
  }
#endif
  
  wfile  -> Close();

  return(0);
}
  

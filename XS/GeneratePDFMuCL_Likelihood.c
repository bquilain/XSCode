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


// you may choose ONE of the following refinements:
// (without any, you build a muCL as muon vs others)
//   - PI_LIKELIHOOD is to build a pion vs proton CL for non-mu-like particles
//   - MUPI_LIKELIHOOD is to build the muCL as muon+pion vs proton

//#define PI_LIKELIHOOD
#define MUPI_LIKELIHOOD

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
  int NFiles=100;
  int IFiles=1;
  int c=-1;
  bool MC=false;
  bool SystematicPE=false;
  int RandomIteration;
  string RandomIteration_string;
  gErrorIgnoreLevel = kError;
  bool PM=true;

  while ((c = getopt(argc, argv, "i:f:d:mwr:")) != -1) {
    switch(c){
    case 'i':
      IFiles=atoi(optarg);
      break;
    case 'f':
      NFiles=atoi(optarg);
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
  char * NameFolder = new char[256];
  char *intStr = new char[256];
  sprintf(intStr,"/%d/",NFiles);
  string NFiles_string = string(intStr);
  char *intStr2 = new char[256];
  sprintf(intStr2,"%d",NFiles);
  string NFiles_string2 = string(intStr2);

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
  TFile * wfile = new TFile(Form("src/PDFMuPiCL_Likelihood%s.root",suffix),"recreate");
#else
  TFile * wfile = new TFile(Form("src/PDFMuCL_Likelihood%s.root",suffix),"recreate");
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
  TH1D * PDFPiCL_PMIng_NotPion = new TH1D("PDFPiCL_PMIng_NotPion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_PMSci_NotPion = new TH1D("PDFPiCL_PMSci_NotPion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_Ing_NotPion = new TH1D("PDFPiCL_Ing_NotPion","",NBins,StartPE,EndPE);
  TH1D * PDFPiCL_Ing_Pion = new TH1D("PDFPiCL_Ing_Pion","",NBins,StartPE,EndPE);
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
  //TFile * _file0;
  TTree * tree;
  
  for(int n=IFiles;n<=NFiles;n++){//Loop over different files
    cout<<"opening"<<endl;
    sprintf(File,"%s/%sMC_Run1_%d_wNoise_ana.root",cMCOUT,DetName,n);

    if(!MC && n>=415 && n<419) continue;
    TFile * _file0 = new TFile(File);
    if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
    else continue;
    tree=(TTree*) _file0->Get("tree");
    if(tree!=tree) continue;
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
      double posX;
      double posY;
      double posZ;
      double norm;
      double totcrsne;

      if(MC){
	mod=simver->mod;
	posX=simver->xnu;
	posY=simver->ynu;
	posZ=simver->znu;
	norm=simver->norm;
	totcrsne=simver->totcrsne;
	weight = 1;//MCCorrections(1,mod);
	if(PM)	weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(1.03*46)/100;	
	else weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*50/100;	
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
	    //Track length:
	    double TrkLength= Rec->Reconstruction::GetTrackLength(Vec);
	    int ParticleNumber=TMath::Abs(Rec->Reconstruction::GetTrackParticle(evt, recon, itrk,TrkLength));
	    int PDG=((IngridSimParticleSummary*) evt->GetSimParticle(ParticleNumber))->pdg;
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
	    
	    // here PM and INGRID are treated
	    for(int i=0;i<Vec.size();i++){
	      if(Vec[i].used>1) continue;

	      if(Vec[i].mod==16 && Rec->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
		PECorrected=Vec[i].pecorr/dx;
		PEPlane[1][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[1][Vec[i].pln][Vec[i].view]++;
	      }
	      else if(Vec[i].mod==16 && !(Rec->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
		PECorrected=(Vec[i].pecorr/(1.3*dx))/INGRIDSCIBAR;
		PEPlane[0][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[0][Vec[i].pln][Vec[i].view]++;
	      }
	      else if(Vec[i].mod<14){
		PECorrected=Vec[i].pecorr/(dx);
		PEPlane[2][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[2][Vec[i].pln][Vec[i].view]++;
	      }
	    }//EndPE of loop over hits


	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int iview=0;iview<2;iview++){
#ifdef MUPI_LIKELIHOOOD
		if(PDG==13 || abs(PDG)==211){
#else 
		if(PDG==13 ){
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
		      else
			PDFPiCL_PMSci_NotPion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#endif
		    }
		    else {
		      PDFMuCL_PMIng_NotMuon->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more INGRID
#ifdef PI_LIKELIHOOD
		      if(abs(PDG)==211)
			PDFPiCL_PMIng_Pion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
		      else
			PDFPiCL_PMIng_NotPion->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#endif
		    }
		  }
		  else if(Plane[2][ipln][iview]!=0) {
		    PDFMuCL_Ing_NotMuon->Fill(PEPlane[2][ipln][iview],weight);		  
#ifdef PI_LIKELIHOOD	
		    if(abs(PDG)==211)
		      PDFPiCL_Ing_Pion->Fill(PEPlane[2][ipln][iview],weight);		  
		    else
		      PDFPiCL_Ing_NotPion->Fill(PEPlane[2][ipln][iview],weight);		  
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
		  else
		    PDFPiCL_WM_NotPion[normalAngleSlice]->Fill(Vec[i].pecorr,weight);
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
  }//Files
  
  wfile->cd();

  PDFParticle->Scale(1./PDFParticle->Integral());
  PDFParticle->Write();
  
  /*double Max_PMIng_Muon=PDFMuCL_PMIng_Muon->GetBinContent(PDFMuCL_PMIng_Muon->GetMaximumBin());
  double Max_PMSci_Muon=PDFMuCL_PMSci_Muon->GetBinContent(PDFMuCL_PMSci_Muon->GetMaximumBin());
  double Max_Ing_Muon=PDFMuCL_Ing_Muon->GetBinContent(PDFMuCL_Ing_Muon->GetMaximumBin());*/

  if(PM){
    PDFMuCL_PMIng_Muon->Scale(1./PDFMuCL_PMIng_Muon->Integral());
    PDFMuCL_PMSci_Muon->Scale(1./PDFMuCL_PMSci_Muon->Integral());
    PDFMuCL_PMIng_Muon->Write();
    PDFMuCL_PMSci_Muon->Write();

    double x_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];double ex_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];double y_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];double ey_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFMuCL_PMIng_Muon->GetNbinsX();ibinx++){
      x_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinCenter(ibinx);
      ex_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinWidth(ibinx)/2;
      y_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinContent(ibinx);
      ey_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinError(ibinx);
    }
    TGraphErrors * g_PMIng_Muon = new TGraphErrors(PDFMuCL_PMIng_Muon->GetNbinsX(),x_PMIng_Muon,y_PMIng_Muon,ex_PMIng_Muon,ey_PMIng_Muon);
    g_PMIng_Muon->Write("g_PMIng_Muon");
    TSpline3* s_PMIng_Muon = new TSpline3("s_PMIng_Muon",g_PMIng_Muon);
    s_PMIng_Muon->Write("s_PMIng_Muon");

    double x_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];double ex_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];double y_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];double ey_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFMuCL_PMSci_Muon->GetNbinsX();ibinx++){
      x_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinCenter(ibinx);
      ex_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinWidth(ibinx)/2;
      y_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinContent(ibinx);
      ey_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinError(ibinx);
    }
    TGraphErrors * g_PMSci_Muon = new TGraphErrors(PDFMuCL_PMSci_Muon->GetNbinsX(),x_PMSci_Muon,y_PMSci_Muon,ex_PMSci_Muon,ey_PMSci_Muon);
    g_PMSci_Muon->Write("g_PMSci_Muon");
    TSpline3 * s_PMSci_Muon = new TSpline3("s_PMSci_Muon",g_PMSci_Muon);
    s_PMSci_Muon->Write("s_PMSci_Muon");

#ifdef PI_LIKELIHOOD
    PDFPiCL_PMIng_Pion->Scale(1./PDFPiCL_PMIng_Pion->Integral());
    PDFPiCL_PMSci_Pion->Scale(1./PDFPiCL_PMSci_Pion->Integral());
    PDFPiCL_PMIng_Pion->Write();
    PDFPiCL_PMSci_Pion->Write();

    double x_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];double ex_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];double y_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];double ey_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMIng_Pion->GetNbinsX();ibinx++){
      x_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinCenter(ibinx);
      ex_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinWidth(ibinx)/2;
      y_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinContent(ibinx);
      ey_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMIng_Pion = new TGraphErrors(PDFPiCL_PMIng_Pion->GetNbinsX(),x_PMIng_Pion,y_PMIng_Pion,ex_PMIng_Pion,ey_PMIng_Pion);
    g_PMIng_Pion->Write("g_PMIng_Pion");
    TSpline3* s_PMIng_Pion = new TSpline3("s_PMIng_Pion",g_PMIng_Pion);
    s_PMIng_Pion->Write("s_PMIng_Pion");

    double x_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];double ex_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];double y_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];double ey_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMSci_Pion->GetNbinsX();ibinx++){
      x_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinCenter(ibinx);
      ex_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinWidth(ibinx)/2;
      y_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinContent(ibinx);
      ey_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMSci_Pion = new TGraphErrors(PDFPiCL_PMSci_Pion->GetNbinsX(),x_PMSci_Pion,y_PMSci_Pion,ex_PMSci_Pion,ey_PMSci_Pion);
    g_PMSci_Pion->Write("g_PMSci_Pion");
    TSpline3 * s_PMSci_Pion = new TSpline3("s_PMSci_Pion",g_PMSci_Pion);
    s_PMSci_Pion->Write("s_PMSci_Pion");
#endif
  }
  else{
    TGraphErrors * g_WM_Muon[3];
    TSpline3 * s_WM_Muon[3];
    TGraphErrors * g_WM_Pion[3];
    TSpline3 * s_WM_Pion[3];

    for(int i=0;i<3;i++){
      PDFMuCL_WM_Muon[i]->Scale(1./PDFMuCL_WM_Muon[i]->Integral());
      PDFMuCL_WM_Muon[i]->Write();

      double x_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];double ex_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];double y_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];double ey_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFMuCL_WM_Muon[i]->GetNbinsX();ibinx++){
	x_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinCenter(ibinx);
	ex_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinWidth(ibinx)/2;
	y_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinContent(ibinx);
	ey_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinError(ibinx);
      }
       g_WM_Muon[i] = new TGraphErrors(PDFMuCL_WM_Muon[i]->GetNbinsX(),x_WM_Muon,y_WM_Muon,ex_WM_Muon,ey_WM_Muon);
       g_WM_Muon[i]->Write(Form("g_WM_Muon_%d",i));
       s_WM_Muon[i] = new TSpline3(Form("s_WM_Muon_%d",i),g_WM_Muon[i]);
       s_WM_Muon[i]->Write(Form("s_WM_Muon_%d",i));

#ifdef PI_LIKELIHOOD
      PDFPiCL_WM_Pion[i]->Scale(1./PDFPiCL_WM_Pion[i]->Integral());
      PDFPiCL_WM_Pion[i]->Write();

      double x_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];double ex_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];double y_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];double ey_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFPiCL_WM_Pion[i]->GetNbinsX();ibinx++){
	x_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinCenter(ibinx);
	ex_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinWidth(ibinx)/2;
	y_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinContent(ibinx);
	ey_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinError(ibinx);
      }
      g_WM_Pion[i] = new TGraphErrors(PDFPiCL_WM_Pion[i]->GetNbinsX(),x_WM_Pion,y_WM_Pion,ex_WM_Pion,ey_WM_Pion);
      g_WM_Pion[i]->Write(Form("g_WM_Pion_%d",i));
      s_WM_Pion[i] = new TSpline3(Form("s_WM_Pion_%d",i),g_WM_Pion[i]);
      s_WM_Pion[i]->Write(Form("s_WM_Pion_%d",i));
#endif
    }
  }

  PDFMuCL_Ing_Muon->Scale(1./PDFMuCL_Ing_Muon->Integral()); 
  PDFMuCL_Ing_Muon->Write();

  double x_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];double ex_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];double y_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];double ey_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_Ing_Muon->GetNbinsX();ibinx++){
    x_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinCenter(ibinx);
    ex_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinWidth(ibinx)/2;
    y_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinContent(ibinx);
    ey_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_Muon = new TGraphErrors(PDFMuCL_Ing_Muon->GetNbinsX(),x_Ing_Muon,y_Ing_Muon,ex_Ing_Muon,ey_Ing_Muon);
  g_Ing_Muon->Write("g_Ing_Muon");
  TSpline3 * s_Ing_Muon = new TSpline3("s_Ing_Muon",g_Ing_Muon);
  s_Ing_Muon->Write("s_Ing_Muon");

#ifdef PI_LIKELIHOOD
  PDFPiCL_Ing_Pion->Scale(1./PDFPiCL_Ing_Pion->Integral()); 
  PDFPiCL_Ing_Pion->Write();

  double x_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];double ex_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];double y_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];double ey_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFPiCL_Ing_Pion->GetNbinsX();ibinx++){
    x_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinCenter(ibinx);
    ex_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinWidth(ibinx)/2;
    y_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinContent(ibinx);
    ey_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_Pion = new TGraphErrors(PDFPiCL_Ing_Pion->GetNbinsX(),x_Ing_Pion,y_Ing_Pion,ex_Ing_Pion,ey_Ing_Pion);
  g_Ing_Pion->Write("g_Ing_Pion");
  TSpline3 * s_Ing_Pion = new TSpline3("s_Ing_Pion",g_Ing_Pion);
  s_Ing_Pion->Write("s_Ing_Pion");
#endif

  // not muon

  /*  double Max_PMIng_NotMuon=PDFMuCL_PMIng_NotMuon->GetBinContent(PDFMuCL_PMIng_NotMuon->GetMaximumBin());
  double Max_PMSci_NotMuon=PDFMuCL_PMSci_NotMuon->GetBinContent(PDFMuCL_PMSci_NotMuon->GetMaximumBin());
  double Max_Ing_NotMuon=PDFMuCL_Ing_NotMuon->GetBinContent(PDFMuCL_Ing_NotMuon->GetMaximumBin());*/
  if(PM){
    PDFMuCL_PMIng_NotMuon->Scale(1./PDFMuCL_PMIng_NotMuon->Integral());
    PDFMuCL_PMSci_NotMuon->Scale(1./PDFMuCL_PMSci_NotMuon->Integral());
    PDFMuCL_PMIng_NotMuon->Write();
    PDFMuCL_PMSci_NotMuon->Write();

    double x_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];double ex_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];double y_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];double ey_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFMuCL_PMIng_NotMuon->GetNbinsX();ibinx++){
      x_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinCenter(ibinx);
      ex_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinWidth(ibinx)/2;
      y_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinContent(ibinx);
      ey_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinError(ibinx);
    }
    TGraphErrors * g_PMIng_NotMuon = new TGraphErrors(PDFMuCL_PMIng_NotMuon->GetNbinsX(),x_PMIng_NotMuon,y_PMIng_NotMuon,ex_PMIng_NotMuon,ey_PMIng_NotMuon);
    g_PMIng_NotMuon->Write("g_PMIng_NotMuon");
    TSpline3* s_PMIng_NotMuon = new TSpline3("s_PMIng_NotMuon",g_PMIng_NotMuon);
    s_PMIng_NotMuon->Write("s_PMIng_NotMuon");

    double x_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];double ex_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];double y_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];double ey_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFMuCL_PMSci_NotMuon->GetNbinsX();ibinx++){
      x_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinCenter(ibinx);
      ex_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinWidth(ibinx)/2;
      y_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinContent(ibinx);
      ey_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinError(ibinx);
    }
    TGraphErrors * g_PMSci_NotMuon = new TGraphErrors(PDFMuCL_PMSci_NotMuon->GetNbinsX(),x_PMSci_NotMuon,y_PMSci_NotMuon,ex_PMSci_NotMuon,ey_PMSci_NotMuon);
    g_PMSci_NotMuon->Write("g_PMSci_NotMuon");
    TSpline3 * s_PMSci_NotMuon = new TSpline3("s_PMSci_NotMuon",g_PMSci_NotMuon);
    s_PMSci_NotMuon->Write("s_PMSci_NotMuon");

#ifdef PI_LIKELIHOOD
    PDFPiCL_PMIng_NotPion->Scale(1./PDFPiCL_PMIng_NotPion->Integral());
    PDFPiCL_PMSci_NotPion->Scale(1./PDFPiCL_PMSci_NotPion->Integral());
    PDFPiCL_PMIng_NotPion->Write();
    PDFPiCL_PMSci_NotPion->Write();

    double x_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];double ex_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];double y_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];double ey_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMIng_NotPion->GetNbinsX();ibinx++){
      x_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinCenter(ibinx);
      ex_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinWidth(ibinx)/2;
      y_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinContent(ibinx);
      ey_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMIng_NotPion = new TGraphErrors(PDFPiCL_PMIng_NotPion->GetNbinsX(),x_PMIng_NotPion,y_PMIng_NotPion,ex_PMIng_NotPion,ey_PMIng_NotPion);
    g_PMIng_NotPion->Write("g_PMIng_NotPion");
    TSpline3* s_PMIng_NotPion = new TSpline3("s_PMIng_NotPion",g_PMIng_NotPion);
    s_PMIng_NotPion->Write("s_PMIng_NotPion");

    double x_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];double ex_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];double y_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];double ey_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMSci_NotPion->GetNbinsX();ibinx++){
      x_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinCenter(ibinx);
      ex_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinWidth(ibinx)/2;
      y_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinContent(ibinx);
      ey_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMSci_NotPion = new TGraphErrors(PDFPiCL_PMSci_NotPion->GetNbinsX(),x_PMSci_NotPion,y_PMSci_NotPion,ex_PMSci_NotPion,ey_PMSci_NotPion);
    g_PMSci_NotPion->Write("g_PMSci_NotPion");
    TSpline3 * s_PMSci_NotPion = new TSpline3("s_PMSci_NotPion",g_PMSci_NotPion);
    s_PMSci_NotPion->Write("s_PMSci_NotPion");
#endif
  }
  else{
    TGraphErrors * g_WM_NotMuon[3];
    TSpline3 * s_WM_NotMuon[3];
    TGraphErrors * g_WM_NotPion[3];
    TSpline3 * s_WM_NotPion[3];

    for(int i=0;i<3;i++){
      PDFMuCL_WM_NotMuon[i]->Scale(1./PDFMuCL_WM_NotMuon[i]->Integral());
      PDFMuCL_WM_NotMuon[i]->Write();

      double x_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];double ex_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];double y_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];double ey_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFMuCL_WM_NotMuon[i]->GetNbinsX();ibinx++){
	x_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinCenter(ibinx);
	ex_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinWidth(ibinx)/2;
	y_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinContent(ibinx);
	ey_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinError(ibinx);
      }
       g_WM_NotMuon[i] = new TGraphErrors(PDFMuCL_WM_NotMuon[i]->GetNbinsX(),x_WM_NotMuon,y_WM_NotMuon,ex_WM_NotMuon,ey_WM_NotMuon);
       g_WM_NotMuon[i]->Write(Form("g_WM_NotMuon_%d",i));
       s_WM_NotMuon[i] = new TSpline3(Form("s_WM_NotMuon_%d",i),g_WM_NotMuon[i]);
       s_WM_NotMuon[i]->Write(Form("s_WM_NotMuon_%d",i));

#ifdef PI_LIKELIHOOD
      PDFPiCL_WM_NotPion[i]->Scale(1./PDFPiCL_WM_NotPion[i]->Integral());
      PDFPiCL_WM_NotPion[i]->Write();

      double x_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];double ex_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];double y_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];double ey_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFPiCL_WM_NotPion[i]->GetNbinsX();ibinx++){
	x_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinCenter(ibinx);
	ex_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinWidth(ibinx)/2;
	y_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinContent(ibinx);
	ey_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinError(ibinx);
      }
       g_WM_NotPion[i] = new TGraphErrors(PDFPiCL_WM_NotPion[i]->GetNbinsX(),x_WM_NotPion,y_WM_NotPion,ex_WM_NotPion,ey_WM_NotPion);
       g_WM_NotPion[i]->Write(Form("g_WM_NotPion_%d",i));
       s_WM_NotPion[i] = new TSpline3(Form("s_WM_NotPion_%d",i),g_WM_NotPion[i]);
       s_WM_NotPion[i]->Write(Form("s_WM_NotPion_%d",i));
#endif
    }
  }


  PDFMuCL_Ing_NotMuon->Scale(1./PDFMuCL_Ing_NotMuon->Integral());
  PDFMuCL_Ing_NotMuon->Write();
  
  double x_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];double ex_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];double y_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];double ey_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_Ing_NotMuon->GetNbinsX();ibinx++){
    x_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinCenter(ibinx);
    ex_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinWidth(ibinx)/2;
    y_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinContent(ibinx);
    ey_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_NotMuon = new TGraphErrors(PDFMuCL_Ing_NotMuon->GetNbinsX(),x_Ing_NotMuon,y_Ing_NotMuon,ex_Ing_NotMuon,ey_Ing_NotMuon);
  g_Ing_NotMuon->Write("g_Ing_NotMuon");
  TSpline3 * s_Ing_NotMuon = new TSpline3("s_Ing_NotMuon",g_Ing_NotMuon);
  s_Ing_NotMuon->Write("s_Ing_NotMuon");

#ifdef LI_LIKELIHOOD
  PDFPiCL_Ing_NotPion->Scale(1./PDFPiCL_Ing_NotPion->Integral());
  PDFPiCL_Ing_NotPion->Write();
  
  double x_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];double ex_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];double y_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];double ey_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFPiCL_Ing_NotPion->GetNbinsX();ibinx++){
    x_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinCenter(ibinx);
    ex_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinWidth(ibinx)/2;
    y_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinContent(ibinx);
    ey_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_NotPion = new TGraphErrors(PDFPiCL_Ing_NotPion->GetNbinsX(),x_Ing_NotPion,y_Ing_NotPion,ex_Ing_NotPion,ey_Ing_NotPion);
  g_Ing_NotPion->Write("g_Ing_NotPion");
  TSpline3 * s_Ing_NotPion = new TSpline3("s_Ing_NotPion",g_Ing_NotPion);
  s_Ing_NotPion->Write("s_Ing_NotPion");
#endif




  
  wfile  -> Close();

  return(0);
}
  

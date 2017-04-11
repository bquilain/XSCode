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
INGRID_Dimension * IngDim = new INGRID_Dimension();
#include "setup.h"
#include "Hit.h"
#include "PMdispRev.h"
#include "Corrections.cc"
#include "Reconstruction.cc"
#include "Xsec.cc"
Reconstruction * Rec = new Reconstruction();
Corrections * Cor = new Corrections();
Xsec * xs = new Xsec();
int NBins=600;
int StartPE=0;
int EndPE=300;
TSpline3 * s_PMIng,* s_PMSci,* s_Ing;
TF1 * f_PMIng,* f_PMSci,* f_Ing;
#define DEBUG

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
/* already defined in Reconstruction.cc
double DegRad(double angle){
  return angle*TMath::Pi()/180.;
}

double RadDeg(double angle){
  return angle*180./TMath::Pi();
}
*/

//TSPLINE CONVERSION TO A TF1
Double_t fPDF_Ing(Double_t *pe, Double_t *par){
  if(s_Ing->Eval(pe[0])>=0){ return s_Ing->Eval(pe[0])/par[0];}
  else {return 0;}
}
Double_t fPDF_PMIng(Double_t *pe, Double_t *par){
  if(s_PMIng->Eval(pe[0])>=0){ return s_PMIng->Eval(pe[0])/par[0];}
  else {return 0;}
}
Double_t fPDF_PMSci(Double_t *pe, Double_t *par){
  if(s_PMSci->Eval(pe[0])>=0){ return s_PMSci->Eval(pe[0])/par[0];}
  else {return 0;}
}


//TF1 PDF CONVERSION TO CUMULATIVE FUNCTION
Double_t Cumulative_Ing_KS(Double_t *pe, Double_t *par){
  if(pe[0]<(par[1]-50.)){return (1.-f_Ing->Integral(par[0],pe[0])/f_Ing->Integral(par[0],par[1]));}
  else{return (1.-f_Ing->Integral(par[0],par[1]-50)/f_Ing->Integral(par[0],par[1]));}
}
Double_t Cumulative_PMIng_KS(Double_t *pe, Double_t *par){
  if(pe[0]<(par[1]-50.)){return (1.-f_PMIng->Integral(par[0],pe[0])/f_PMIng->Integral(par[0],par[1]));}
  else{return (1.-f_PMIng->Integral(par[0],par[1]-50)/f_PMIng->Integral(par[0],par[1]));}
}
Double_t Cumulative_PMSci_KS(Double_t *pe, Double_t *par){
  if(pe[0]<(par[1]-50.)){return (1.-f_PMSci->Integral(par[0],pe[0])/f_PMSci->Integral(par[0],par[1]));}
  else{return (1.-f_PMSci->Integral(par[0],par[1]-50)/f_PMSci->Integral(par[0],par[1]));}
}


Double_t Cumulative_Ing(Double_t *pe, Double_t *par){
  if(pe[0]<(par[2]-50.)){
    if(pe[0]>=par[1]){
      return (1.-f_Ing->Integral(par[1],pe[0])/f_Ing->Integral(par[1],par[2]));
    }
    else{ return (1.-f_Ing->Integral(pe[0],par[1])/f_Ing->Integral(par[0],par[1]));}
  }
  else{return (1.-f_Ing->Integral(par[1],par[2]-50)/f_Ing->Integral(par[1],par[2]));}
}

Double_t Cumulative_PMIng(Double_t *pe, Double_t *par){
  if(pe[0]<(par[2]-50.)){
    if(pe[0]>=par[1]){
      return (1.-f_PMIng->Integral(par[1],pe[0])/f_PMIng->Integral(par[1],par[2]));
    }
    else{ return (1.-f_PMIng->Integral(pe[0],par[1])/f_PMIng->Integral(par[0],par[1]));}
  }
  else{return (1.-f_PMIng->Integral(par[1],par[2]-50)/f_PMIng->Integral(par[1],par[2]));}
}

Double_t Cumulative_PMSci(Double_t *pe, Double_t *par){
  if(pe[0]<(par[2]-50.)){
    if(pe[0]>=par[1]){
      return (1.-f_PMSci->Integral(par[1],pe[0])/f_PMSci->Integral(par[1],par[2]));
    }
    else{ return (1.-f_PMSci->Integral(pe[0],par[1])/f_PMSci->Integral(par[0],par[1]));}
  }
  else{return (1.-f_PMSci->Integral(par[1],par[2]-50)/f_PMSci->Integral(par[1],par[2]));}
}

/*
Double_t Cumulative_Ing(Double_t *pe, Double_t *par){
  if(pe[0]<(par[2]-50.)){
    if(pe[0]>=par[1]){
      return (1.-f_Ing->Integral(par[1],pe[0])/f_Ing->Integral(par[1],par[2]));
    }
    else{ return (1.-f_Ing->Integral(pe[0],par[1])/f_Ing->Integral(par[0],par[1]));}
  }
  else{return (1.-f_Ing->Integral(par[1],par[2]-50)/f_Ing->Integral(par[1],par[2]));}
}

Double_t Cumulative_PMIng(Double_t *pe, Double_t *par){
  if(pe[0]<(par[2]-50.)){
    if(pe[0]>=par[1]){
      return (1.-f_PMIng->Integral(par[1],pe[0])/f_PMIng->Integral(par[1],par[2]));
    }
    else{ return (1.-f_PMIng->Integral(pe[0],par[1])/f_PMIng->Integral(par[0],par[1]));}
  }
  else{return (1.-f_PMIng->Integral(par[1],par[2]-50)/f_PMIng->Integral(par[1],par[2]));}
}

Double_t Cumulative_PMSci(Double_t *pe, Double_t *par){
  if(pe[0]<(par[2]-50.)){
    if(pe[0]>=par[1]){
      return (1.-f_PMSci->Integral(par[1],pe[0])/f_PMSci->Integral(par[1],par[2]));
    }
    else{ return (1.-f_PMSci->Integral(pe[0],par[1])/f_PMSci->Integral(par[0],par[1]));}
  }
  else{return (1.-f_PMSci->Integral(par[1],par[2]-50)/f_PMSci->Integral(par[1],par[2]));}
}
*/

int main(int argc, char **argv)
{
  bool Disp=false;
  int NFiles=100;
  int IFiles=1;
  int c=-1;
  bool MC=false;
  bool SystematicPE=false;
  int RandomIteration;
  string RandomIteration_string;
  gErrorIgnoreLevel = kError;

  while ((c = getopt(argc, argv, "i:f:d:mr:")) != -1) {
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

  int type;//0=CC, 1=NC
  double Nu_E;
  char File[256];
  int nTracks=-1;
  float weight=1;
  int NIngBasRec;


  char *NameDist = new char[256];
  char * FileDist = new char[256];



  TFile * wfile;
  if(MC) wfile    = new TFile("src/PDFMuCL_Plan_test.root","recreate");
  else wfile    = new TFile("src/PDFMuCL_Plan_Data_test.root","recreate");
  TH1D * PDFMuCL_PMIng = new TH1D("PDFMuCL_PMIng","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_PMSci = new TH1D("PDFMuCL_PMSci","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_Ing = new TH1D("PDFMuCL_Ing","",NBins,StartPE,EndPE);
  PDFMuCL_PMIng->Sumw2();PDFMuCL_PMSci->Sumw2();PDFMuCL_Ing->Sumw2();

#ifdef DEBUG
  TH1D * PDFMuCL_PMIng_raw = new TH1D("PDFMuCL_PMIng_raw","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_PMSci_raw = new TH1D("PDFMuCL_PMSci_raw","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_Ing_raw = new TH1D("PDFMuCL_Ing_raw","",NBins,StartPE,EndPE);

  TH1D * PDFMuCL_PMIng_fiberattenuationcorrected = new TH1D("PDFMuCL_PMIng_fiberattenuationcorrected","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_PMSci_fiberattenuationcorrected = new TH1D("PDFMuCL_PMSci_fiberattenuationcorrected","",NBins,StartPE,EndPE);
  TH1D * PDFMuCL_Ing_fiberattenuationcorrected = new TH1D("PDFMuCL_Ing_fiberattenuationcorrected","",NBins,StartPE,EndPE);
#endif
  
  double PECorrected=0;

  cout<<"Opening Events"<<endl;

  IngridEventSummary* evt;
  TBranch * Br;
  IngridSimVertexSummary * simver = new IngridSimVertexSummary();//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
  IngridHitSummary * Hit = new IngridHitSummary();
  IngridSimHitSummary * SimHit= new IngridSimHitSummary();
  IngridSimHitSummary * SimHit2= new IngridSimHitSummary();
  IngridSimParticleSummary * SimPart= new IngridSimParticleSummary();
  BeamInfoSummary * BeamSummary = new BeamInfoSummary();
  PMAnaSummary * recon = new PMAnaSummary();	
  IngridHitSummary * hit = new IngridHitSummary();
  IngridHitSummary * hit2 =new IngridHitSummary();

  TChain * tree = new TChain("tree");
  for(int n=IFiles;n<=NFiles;n++){//Loop over different files
    if(MC)tree->Add(Form("/export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root",n));
    else tree->Add(Form("/export/scraid2/data/bquilain/DataNew/ingrid_00014510_%04d_pmmergedKSPManabsd_woXTalk.root",n));
    cout<<"opening"<<endl;
  }
    //if(n>=415 && n<419) continue;
    //TFile * _file0 = new TFile(File);
    int nevt=(int) tree->GetEntries();
    cout<<"Total Number Of Events="<<nevt<<endl;

    evt = new IngridEventSummary();
    Br=tree->GetBranch("fDefaultReco.");
    Br->SetAddress(&evt);
    tree->SetBranchAddress("fDefaultReco.",&evt);
 
    for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
      if((ievt%100)==0) cout<<"Processing "<<ievt<<endl;
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
	weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*470)*1.61/100;	
      }

      NIngBasRec= evt->NPMAnas();

      for(int irec=0;irec<NIngBasRec;irec++){
	recon = (PMAnaSummary*) evt->GetPMAna(irec);	
	nTracks=recon->Ntrack;
      
	if((recon->vetowtracking) && !(recon->edgewtracking)) {


	  
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


#ifdef DEBUG2
	  cout<<endl<<"New vertex"<<endl;
#endif
	  
	  for(int itrk=0;itrk<nTracks;itrk++){//loop on track
#ifdef DEBUG2
	    cout<<"New track:"<<endl;
	    cout<<"Startxpln="<<(recon->startxpln)[itrk]/10.<<", Startypln="<<(recon->startypln)[itrk]/10.<<", Startxch="<<(recon->startxch)[itrk]/10.<<", Startych="<<(recon->startych)[itrk]/10.<<endl;
	    cout<<"Endxpln="<<(recon->endxpln)[itrk]/10.<<", Endypln="<<(recon->endypln)[itrk]/10.<<", Endxch="<<(recon->endxch)[itrk]/10.<<", Endych="<<(recon->endych)[itrk]/10.<<endl;
	    //cout<<"Gradient x="<<(recon->gradx)[itrk]<<", intcpt="<<(recon->intcptx)[itrk]<<", thetax="<<(recon->thetax)[itrk]<<endl;
	    //cout<<"Gradient y="<<(recon->grady][itrk]<<", intcpt="<<(recon->intcpty)[itrk]<<", thetay="<<(recon->thetay)[itrk]<<endl;
#endif
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

	    double PEPlane[3][NPlnPM][2];//one per type of scinti, plane and view
	    int Plane[3][NPlnPM][2];
#ifdef DEBUG
	    
	    double PEPlane_raw[3][NPlnPM][2];//one per type of scinti, plane and view
	    double PEPlane_fiberattenuationcorrected[3][NPlnPM][2];//one per type of scinti, plane and view
#endif
	    for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int i=0;i<3;i++){
		for(int iv=0;iv<2;iv++){
		  PEPlane[i][ipln][iv]=0;
		  Plane[i][ipln][iv]=0;
#ifdef DEBUG
		  PEPlane_raw[i][ipln][iv]=0;
		  PEPlane_fiberattenuationcorrected[i][ipln][iv]=0;
#endif
		}
	      }
	    }
	      
	    for(int i=0;i<Vec.size();i++){
	      if(Vec[i].used>1) continue;

	      if(Vec[i].mod==16 && Rec->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
		PECorrected=Vec[i].pecorr/dx;
		PEPlane[1][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[1][Vec[i].pln][Vec[i].view]++;
#ifdef DEBUG
		PEPlane_raw[1][Vec[i].pln][Vec[i].view]+=Vec[i].pe;
		PEPlane_fiberattenuationcorrected[1][Vec[i].pln][Vec[i].view]+=Vec[i].pecorr;
#endif
	      }
	      else if(Vec[i].mod==16 && !(Rec->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch))){
		PECorrected=(Vec[i].pecorr/(1.3*dx))/INGRIDSCIBAR;
		PEPlane[0][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[0][Vec[i].pln][Vec[i].view]++;

#ifdef DEBUG
		PEPlane_raw[0][Vec[i].pln][Vec[i].view]+=Vec[i].pe/(1.3*INGRIDSCIBAR);
		PEPlane_fiberattenuationcorrected[0][Vec[i].pln][Vec[i].view]+=Vec[i].pecorr/(1.3*INGRIDSCIBAR);
#endif
	      }
	      else{
		PECorrected=Vec[i].pecorr/(dx);
		PEPlane[2][Vec[i].pln][Vec[i].view]+=PECorrected;
		Plane[2][Vec[i].pln][Vec[i].view]++;
#ifdef DEBUG
		PEPlane_raw[2][Vec[i].pln][Vec[i].view]+=Vec[i].pe;
		PEPlane_fiberattenuationcorrected[2][Vec[i].pln][Vec[i].view]+=Vec[i].pecorr;
#endif
	      }
	    }//EndPE of loop over hits


	  for(int ipln=0;ipln<NPlnPM;ipln++){
	      for(int iview=0;iview<2;iview++){
		if(Plane[0][ipln][iview]!=0 || Plane[1][ipln][iview]!=0){//case PM & active plane
		  if(Plane[0][ipln][iview]>=Plane[1][ipln][iview]){
		    PDFMuCL_PMSci->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more SciBar
#ifdef DEBUG
		    PDFMuCL_PMSci_raw->Fill(PEPlane_raw[0][ipln][iview]+PEPlane_raw[1][ipln][iview],weight);//more SciBar
		    PDFMuCL_PMSci_fiberattenuationcorrected->Fill(PEPlane_fiberattenuationcorrected[0][ipln][iview]+PEPlane_fiberattenuationcorrected[1][ipln][iview],weight);//more SciBar
#endif
		  }
		  else{
		    PDFMuCL_PMIng->Fill(PEPlane[0][ipln][iview]+PEPlane[1][ipln][iview],weight);//more INGRID
#ifdef DEBUG
		    PDFMuCL_PMIng_raw->Fill(PEPlane_raw[0][ipln][iview]+PEPlane_raw[1][ipln][iview],weight);//more SciBar
		    PDFMuCL_PMIng_fiberattenuationcorrected->Fill(PEPlane_fiberattenuationcorrected[0][ipln][iview]+PEPlane_fiberattenuationcorrected[1][ipln][iview],weight);//more SciBar
#endif
	  
		  }
		}
		else if(Plane[2][ipln][iview]!=0){
		  PDFMuCL_Ing->Fill(PEPlane[2][ipln][iview],weight);
#ifdef DEBUG
		  PDFMuCL_PMIng_raw->Fill(PEPlane_raw[2][ipln][iview],weight);//more IngBar
		  PDFMuCL_PMIng_fiberattenuationcorrected->Fill(PEPlane_fiberattenuationcorrected[2][ipln][iview],weight);//more SciBar
#endif
		}
	      }
	    }

	    
	    Vec.clear();
	  }//Tracks
	}//cuts
      }//Evt
    }
    evt->Delete();
  
  double Max_PMIng=PDFMuCL_PMIng->GetBinContent(PDFMuCL_PMIng->GetMaximumBin());
  double Max_PMSci=PDFMuCL_PMSci->GetBinContent(PDFMuCL_PMSci->GetMaximumBin());
  double Max_Ing=PDFMuCL_Ing->GetBinContent(PDFMuCL_Ing->GetMaximumBin());
  PDFMuCL_PMIng->Scale(1./Max_PMIng);
  PDFMuCL_PMSci->Scale(1./Max_PMSci);
  PDFMuCL_Ing->Scale(1./Max_Ing);
#ifdef DEBUG
  double Max_PMIng_raw=PDFMuCL_PMIng_raw->GetBinContent(PDFMuCL_PMIng_raw->GetMaximumBin());
  double Max_PMSci_raw=PDFMuCL_PMSci_raw->GetBinContent(PDFMuCL_PMSci_raw->GetMaximumBin());
  double Max_Ing_raw=PDFMuCL_Ing_raw->GetBinContent(PDFMuCL_Ing_raw->GetMaximumBin());
  PDFMuCL_PMIng_raw->Scale(1./Max_PMIng_raw);
  PDFMuCL_PMSci_raw->Scale(1./Max_PMSci_raw);
  PDFMuCL_Ing_raw->Scale(1./Max_Ing_raw);

  double Max_PMIng_fiberattenuationcorrected=PDFMuCL_PMIng_fiberattenuationcorrected->GetBinContent(PDFMuCL_PMIng_fiberattenuationcorrected->GetMaximumBin());
  double Max_PMSci_fiberattenuationcorrected=PDFMuCL_PMSci_fiberattenuationcorrected->GetBinContent(PDFMuCL_PMSci_fiberattenuationcorrected->GetMaximumBin());
  double Max_Ing_fiberattenuationcorrected=PDFMuCL_Ing_fiberattenuationcorrected->GetBinContent(PDFMuCL_Ing_fiberattenuationcorrected->GetMaximumBin());
  PDFMuCL_PMIng_fiberattenuationcorrected->Scale(1./Max_PMIng_fiberattenuationcorrected);
  PDFMuCL_PMSci_fiberattenuationcorrected->Scale(1./Max_PMSci_fiberattenuationcorrected);
  PDFMuCL_Ing_fiberattenuationcorrected->Scale(1./Max_Ing_fiberattenuationcorrected);
#endif
  double x_PMIng[PDFMuCL_PMIng->GetNbinsX()];double ex_PMIng[PDFMuCL_PMIng->GetNbinsX()];double y_PMIng[PDFMuCL_PMIng->GetNbinsX()];double ey_PMIng[PDFMuCL_PMIng->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMIng->GetNbinsX();ibinx++){
    x_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinCenter(ibinx);
    ex_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinWidth(ibinx)/2;
    y_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinContent(ibinx);
    ey_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinError(ibinx);
  }
  TGraphErrors * g_PMIng = new TGraphErrors(PDFMuCL_PMIng->GetNbinsX(),x_PMIng,y_PMIng,ex_PMIng,ey_PMIng);
  //g_PMIng->Write("g_PMIng");
  s_PMIng = new TSpline3("s_PMIng",g_PMIng);

  double x_PMSci[PDFMuCL_PMSci->GetNbinsX()];double ex_PMSci[PDFMuCL_PMSci->GetNbinsX()];double y_PMSci[PDFMuCL_PMSci->GetNbinsX()];double ey_PMSci[PDFMuCL_PMSci->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMSci->GetNbinsX();ibinx++){
    x_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinCenter(ibinx);
    ex_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinWidth(ibinx)/2;
    y_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinContent(ibinx);
    ey_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinError(ibinx);
  }
  TGraphErrors * g_PMSci = new TGraphErrors(PDFMuCL_PMSci->GetNbinsX(),x_PMSci,y_PMSci,ex_PMSci,ey_PMSci);
  //g_PMSci->Write("g_PMSci");
  s_PMSci = new TSpline3("s_PMSci",g_PMSci);

    double x_Ing[PDFMuCL_Ing->GetNbinsX()];double ex_Ing[PDFMuCL_Ing->GetNbinsX()];double y_Ing[PDFMuCL_Ing->GetNbinsX()];double ey_Ing[PDFMuCL_Ing->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_Ing->GetNbinsX();ibinx++){
    x_Ing[ibinx-1]=PDFMuCL_Ing->GetBinCenter(ibinx);
    ex_Ing[ibinx-1]=PDFMuCL_Ing->GetBinWidth(ibinx)/2;
    y_Ing[ibinx-1]=PDFMuCL_Ing->GetBinContent(ibinx);
    ey_Ing[ibinx-1]=PDFMuCL_Ing->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing = new TGraphErrors(PDFMuCL_Ing->GetNbinsX(),x_Ing,y_Ing,ex_Ing,ey_Ing);
  //g_Ing->Write("g_Ing");
  s_Ing = new TSpline3("s_Ing",g_Ing);

  //################################CONVERT SPLINES TO TF1 (TO USE INTEGRAL FUNCTION) & build cumulative distribution
  
  //create PDF functions
  double Start_PMIng=s_PMIng->GetXmin();double End_PMIng=s_PMIng->GetXmax();
  double Start_PMSci=s_PMSci->GetXmin();double End_PMSci=s_PMSci->GetXmax();
  double Start_Ing=s_Ing->GetXmin();double End_Ing=s_Ing->GetXmax();
  //cout<<"Start="<<s_PMIng->GetXmin()<<", end="<<s_PMIng->GetXmax()<<endl;
  f_PMIng = new TF1("f_PMIng",fPDF_PMIng,Start_PMIng,End_PMIng,1);
  f_PMSci = new TF1("f_PMSci",fPDF_PMSci,Start_PMSci,End_PMSci,1);
  f_Ing = new TF1("f_Ing",fPDF_Ing,Start_Ing,End_Ing,1);
  
  f_PMIng->SetParameter(0,1.);
  f_PMSci->SetParameter(0,1.);
  f_Ing->SetParameter(0,1.);
  
  //cout<<f_PMIng->Eval(10)<<endl;
  double MaxX_PMIng=f_PMIng->GetMaximumX(Start_PMIng,End_PMIng);
  double MaxX_PMSci=f_PMSci->GetMaximumX(Start_PMSci,End_PMSci);
  double MaxX_Ing=f_Ing->GetMaximumX(Start_Ing,End_Ing);

  //create Cumulative
  TF1 * CL_PMIng = new TF1("CL_PMIng",Cumulative_PMIng_KS,Start_PMIng,End_PMIng,2);
  TF1 * CL_PMSci = new TF1("CL_PMSci",Cumulative_PMSci_KS,Start_PMSci,End_PMSci,2);
  TF1 * CL_Ing = new TF1("CL_Ing",Cumulative_Ing_KS,Start_Ing,End_Ing,2);

  CL_PMIng->SetParameter(0,Start_PMIng);CL_PMIng->SetParameter(1,End_PMIng);
  CL_PMSci->SetParameter(0,Start_PMSci);CL_PMSci->SetParameter(1,End_PMSci);
  CL_Ing->SetParameter(0,Start_Ing);CL_Ing->SetParameter(1,End_Ing);

  //##############################################################################################################
  
  int NBins_PMIng=g_PMIng->GetN();
  double * cx_PMIng=(double*) g_PMIng->GetX();
  double * cy_PMIng = new double[NBins_PMIng];
  for(int i=0;i<NBins_PMIng;i++){
    cy_PMIng[i]=CL_PMIng->Eval(cx_PMIng[i]);
  }
  TGraph * g_Cumulative_PMIng = new TGraphErrors(NBins_PMIng,cx_PMIng,cy_PMIng);
  TSpline3 * sCL_PMIng = new TSpline3("sCL_PMIng",g_Cumulative_PMIng);
		   
  int NBins_PMSci=g_PMSci->GetN();
  double * cx_PMSci=(double*) g_PMSci->GetX();
  double * cy_PMSci = new double[NBins_PMSci];
  for(int i=0;i<NBins_PMSci;i++){
    cy_PMSci[i]=CL_PMSci->Eval(cx_PMSci[i]);
  }
  TGraph * g_Cumulative_PMSci = new TGraphErrors(NBins_PMSci,cx_PMSci,cy_PMSci);
  TSpline3 * sCL_PMSci = new TSpline3("sCL_PMSci",g_Cumulative_PMSci);

  int NBins_Ing=g_Ing->GetN();
  double * cx_Ing=(double*) g_Ing->GetX();
  double * cy_Ing = new double[NBins_Ing];
  for(int i=0;i<NBins_Ing;i++){
    cy_Ing[i]=CL_Ing->Eval(cx_Ing[i]);
  }
  TGraph * g_Cumulative_Ing = new TGraphErrors(NBins_Ing,cx_Ing,cy_Ing);
  TSpline3 * sCL_Ing = new TSpline3("sCL_Ing",g_Cumulative_Ing);

  wfile->cd();
  PDFMuCL_PMIng->Write();
  PDFMuCL_PMSci->Write();
  PDFMuCL_Ing->Write();
#ifdef DEBUG
  PDFMuCL_PMIng_raw->Write();
  PDFMuCL_PMSci_raw->Write();
  PDFMuCL_Ing_raw->Write();
  PDFMuCL_PMIng_fiberattenuationcorrected->Write();
  PDFMuCL_PMSci_fiberattenuationcorrected->Write();
  PDFMuCL_Ing_fiberattenuationcorrected->Write();
#endif
  s_PMIng->Write("s_PMIng");
  s_PMSci->Write("s_PMSci");
  s_Ing->Write("s_Ing");
  //g_Cumulative_PMIng->Write("Cumulative_PMIng");
  //g_Cumulative_PMSci->Write("Cumulative_PMSci");
  //g_Cumulative_Ing->Write("Cumulative_Ing");
  sCL_PMIng->Write("sCL_PMIng");
  sCL_PMSci->Write("sCL_PMSci");
  sCL_Ing->Write("sCL_Ing");

  
  //cout<<"CL TF1 in three points="<<CL_PMIng->Eval(4)<<", "<<CL_PMIng->Eval(10)<<", "<<CL_PMIng->Eval(30)<<", and Spline="<<sCL_PMIng->Eval(4)<<", "<<sCL_PMIng->Eval(10)<<", "<<sCL_PMIng->Eval(30)<<endl;
  
  wfile  -> Close();

  return(0);
}
  

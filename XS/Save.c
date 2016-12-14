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
int LimitTracks=20;
int LimitRecs=5;
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
Xsec * xs = new Xsec();

int NBins=600;
int Start=0;
int End=300;

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
  wfile    = new TFile("src/PDFMuCL.root","recreate");
  TH1D * PDFMuCL_PMIng = new TH1D("PDFMuCL_PMIng","",NBins,Start,End);
  TH1D * PDFMuCL_PMSci = new TH1D("PDFMuCL_PMSci","",NBins,Start,End);
  TH1D * PDFMuCL_Ing = new TH1D("PDFMuCL_Ing","",NBins,Start,End);
  PDFMuCL_PMIng->Sumw2();PDFMuCL_PMSci->Sumw2();PDFMuCL_Ing->Sumw2();
  
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
  //TFile * _file0;
  TTree * tree;
  
  for(int n=IFiles;n<=NFiles;n++){//Loop over different files
    cout<<"opening"<<endl;
    sprintf(File,"/export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root",n);
    if(n>=415 && n<419) continue;
    TFile * _file0 = new TFile(File);
    if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
    else continue;
    tree=(TTree*) _file0->Get("tree");
    if(tree!=tree) continue;
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
	    
	    for(int i=0;i<Vec.size();i++){
	      //cout<<"Module="<<Vec[i].mod<<", plane="<<Vec[i].pln<<", view="<<Vec[i].view<<", x="<<Vec[i].x<<", y="<<Vec[i].y<<endl;
	      if(Vec[i].used>1) continue;
	      
	      if(Vec[i].mod==16 && IsINGRID(Vec[i].ch)){
		PECorrected=Vec[i].pecorr/dx;
		PDFMuCL_PMIng->Fill(PECorrected,weight);
		//cout<<"hello"<<endl;
	      }
	      else if(Vec[i].mod==16 && !(IsINGRID(Vec[i].ch))){
		PECorrected=Vec[i].pecorr/(1.3*dx);
		PDFMuCL_PMSci->Fill(PECorrected,weight);
	      }
	      else{
		PECorrected=Vec[i].pecorr/(dx);
		PDFMuCL_Ing->Fill(PECorrected,weight);
	      }
	    }//End of loop over hits
	    Vec.clear();
	  }//Tracks
	}//cuts
      }//Evt
    }
    evt->Delete();
    Br->Delete();
  }//Files
  
  wfile->cd();
  double Max_PMIng=PDFMuCL_PMIng->GetBinContent(PDFMuCL_PMIng->GetMaximumBin());
  double Max_PMSci=PDFMuCL_PMSci->GetBinContent(PDFMuCL_PMSci->GetMaximumBin());
  double Max_Ing=PDFMuCL_Ing->GetBinContent(PDFMuCL_Ing->GetMaximumBin());
  PDFMuCL_PMIng->Scale(1./Max_PMIng);
  PDFMuCL_PMSci->Scale(1./Max_PMSci);
  PDFMuCL_Ing->Scale(1./Max_Ing);
  
  PDFMuCL_PMIng->Write();
  PDFMuCL_PMSci->Write();
  PDFMuCL_Ing->Write();

  double x_PMIng[PDFMuCL_PMIng->GetNbinsX()];double ex_PMIng[PDFMuCL_PMIng->GetNbinsX()];double y_PMIng[PDFMuCL_PMIng->GetNbinsX()];double ey_PMIng[PDFMuCL_PMIng->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMIng->GetNbinsX();ibinx++){
    x_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinCenter(ibinx);
    ex_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinWidth(ibinx)/2;
    y_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinContent(ibinx);
    ey_PMIng[ibinx-1]=PDFMuCL_PMIng->GetBinError(ibinx);
  }
  TGraphErrors * g_PMIng = new TGraphErrors(PDFMuCL_PMIng->GetNbinsX(),x_PMIng,y_PMIng,ex_PMIng,ey_PMIng);
  g_PMIng->Write("g_PMIng");
  s_PMIng = new TSpline5("s_PMIng",g_PMIng);
  s_PMIng->Write("s_PMIng");

  double x_PMSci[PDFMuCL_PMSci->GetNbinsX()];double ex_PMSci[PDFMuCL_PMSci->GetNbinsX()];double y_PMSci[PDFMuCL_PMSci->GetNbinsX()];double ey_PMSci[PDFMuCL_PMSci->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMSci->GetNbinsX();ibinx++){
    x_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinCenter(ibinx);
    ex_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinWidth(ibinx)/2;
    y_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinContent(ibinx);
    ey_PMSci[ibinx-1]=PDFMuCL_PMSci->GetBinError(ibinx);
  }
  TGraphErrors * g_PMSci = new TGraphErrors(PDFMuCL_PMSci->GetNbinsX(),x_PMSci,y_PMSci,ex_PMSci,ey_PMSci);
  g_PMSci->Write("g_PMSci");
  s_PMSci = new TSpline5("s_PMSci",g_PMSci);
  s_PMSci->Write("s_PMSci");

    double x_Ing[PDFMuCL_Ing->GetNbinsX()];double ex_Ing[PDFMuCL_Ing->GetNbinsX()];double y_Ing[PDFMuCL_Ing->GetNbinsX()];double ey_Ing[PDFMuCL_Ing->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_Ing->GetNbinsX();ibinx++){
    x_Ing[ibinx-1]=PDFMuCL_Ing->GetBinCenter(ibinx);
    ex_Ing[ibinx-1]=PDFMuCL_Ing->GetBinWidth(ibinx)/2;
    y_Ing[ibinx-1]=PDFMuCL_Ing->GetBinContent(ibinx);
    ey_Ing[ibinx-1]=PDFMuCL_Ing->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing = new TGraphErrors(PDFMuCL_Ing->GetNbinsX(),x_Ing,y_Ing,ex_Ing,ey_Ing);
  g_Ing->Write("g_Ing");
  s_Ing = new TSpline5("s_Ing",g_Ing);
  s_Ing->Write("s_Ing");

  
  f_PMIng = new TF1("f_PMIng",xs->Xsec::fPDF_PMIng,Start,End,1);
  f_PMSci = new TF1("f_PMSci",xs->Xsec::fPDF_PMSci,Start,End,1);
  f_Ing = new TF1("f_Ing",xs->Xsec::fPDF_Ing,Start,End,1);
  f_PMIng->SetParameter(0,1.);
  f_PMSci->SetParameter(0,1.);
  f_Ing->SetParameter(0,1.);

  f_PMIng->Write();
  f_PMSci->Write();
  f_Ing->Write();

  TF1 * CL_PMIng = new TF1("CL_PMIng",xs->Xsec::Cumulative_PMIng,Start,End,1);
  TF1 * CL_PMSci = new TF1("CL_PMSci",xs->Xsec::Cumulative_PMSci,Start,End,1);
  TF1 * CL_Ing = new TF1("CL_Ing",xs->Xsec::Cumulative_Ing,Start,End,1);

  double MaxX_PMIng=f_PMIng->GetMaximumX(Start,End);
  double MaxX_PMSci=f_PMSci->GetMaximumX(Start,End);
  double MaxX_Ing=f_Ing->GetMaximumX(Start,End);
  CL_PMIng->SetParameter(0,MaxX_PMIng);
  CL_PMSci->SetParameter(0,MaxX_PMSci);
  CL_Ing->SetParameter(0,MaxX_Ing);

  cout<<MaxX_PMIng<<", "<<MaxX_PMSci<<", "<<MaxX_Ing<<endl;
  cout<<CL_PMIng->Eval(MaxX_PMIng)<<endl;
  
  CL_PMIng->Write();
  CL_PMSci->Write();
  CL_Ing->Write();

  
  wfile  -> Close();

  return(0);
}
  

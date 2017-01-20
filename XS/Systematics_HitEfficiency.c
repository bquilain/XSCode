
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
int LimitTracks=100;
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
Reconstruction * Rec = new Reconstruction();
Corrections * Cor = new Corrections();

//#define DEBUG
//#define XSEC_ERROR
//#define GENERATEWIDTH
//double C[17]={-1.85,0.34,0.59,0.74,0.514,-0.37,1.25,-0.06,0.562,0.82,-0.47,0.6,-0.57,-0.45};
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

  while ((c = getopt(argc, argv, "i:f:d:m")) != -1) {
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
    }
  }
  char * NameFolder = new char[256];
  char *intStr = new char[256];
  sprintf(intStr,"/%d/",NFiles);
  string NFiles_string = string(intStr);

  cout<<"Welcome"<<endl;
  vector <HitTemp> HitV;
  TApplication theApp("App",0,0);

  int type;//0=CC, 1=NC
  double Nu_E;
  char File[256];
  float MeV2PE=46.0; //MIP9
  double TrueParticleNRJ=0;
  int IntNumber=0;

  int FSIInt=-1;
  int Num_Int=-1;
  int nTracks=-1;
  float weight=1;
  bool IsFV=false;
  bool IsDetected=false;
  bool IsSand=false;
  float POT;
  float Enu;
  float TrueAngleMuon;
  float TrueMomentumMuon;
  int NIngBasRec;
  TFile * _InputEff_Ingrid;
  TH1D * InputEff_Ingrid;
  TFile * _InputEffMC_Ingrid;
  TH1D * InputEffMC_Ingrid;

  TFile * _InputEff_SciBar;
  TH1D * InputEff_SciBar;
  TFile * _InputEffMC_SciBar;
  TH1D * InputEffMC_SciBar;

  TRandom3 * Rand = new TRandom3();

  _InputEff_Ingrid = new TFile("HitEfficiency/PMEfficiencyX_Ing.root");
  InputEff_Ingrid = (TH1D*) _InputEff_Ingrid->Get("hEffX_Ing");
  _InputEff_SciBar = new TFile("HitEfficiency/PMEfficiencyX_Sci.root");
  InputEff_SciBar = (TH1D*) _InputEff_SciBar->Get("hEffX_Sci");
 
 _InputEffMC_Ingrid = new TFile("HitEfficiency/MC_PMEfficiencyX_Ing.root");
  InputEffMC_Ingrid = (TH1D*) _InputEffMC_Ingrid->Get("hEffX_Ing");
  _InputEffMC_SciBar = new TFile("HitEfficiency/MC_PMEfficiencyX_Sci.root");
  InputEffMC_SciBar = (TH1D*) _InputEffMC_SciBar->Get("hEffX_Sci");
  
  TFile*            wfile = new TFile("test.root", "recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* evt2 = new IngridEventSummary();
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary",
                                  &evt2,  64000,  99);
  
  cout<<"Opening Events"<<endl;

  for(int a=1;a<2;a++){
    if(!MC && a!=1) continue;
  
    //gErrorIgnoreLevel = Error;
  for(int n=IFiles;n<=NFiles;n++){//Loop over different files
    if(MC && a==0) sprintf(File,"/export/scraid2/data/bquilain/MCfiles/PM_MC_Sand%d_BirksCorrectedMIP40_ReWeight_wNoise_KSana_woXTalk.root",n); 
    else if(MC && a==1) sprintf(File,"/export/scraid2/data/bquilain/MCfiles/PM_MC_Beam%d_BirksCorrectedMIP40_ReWeight_wNoise_KSana_woXTalk.root",n);
    else sprintf(File,"/export/scraid2/data/bquilain/ingrid_00015769_%04d_pmmergedKSPMabsd_woXTalk.root",n);

    TFile * _file0=new TFile(File);
    if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
    else continue;
    TTree * tree=(TTree*) _file0->Get("tree");
    if(tree!=tree) continue;
    int nevt=(int) tree->GetEntries();
    cout<<"Total Number Of Events="<<nevt<<endl;

    IngridEventSummary* evt = new IngridEventSummary();
    TBranch * Br=tree->GetBranch("fDefaultReco.");
    Br->SetAddress(&evt);
    tree->SetBranchAddress("fDefaultReco.",&evt);

    for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
      if((ievt%100)==0) cout<<"Processing "<<ievt<<endl;
      evt->Clear();
      evt2->Clear();
      tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
      if(evt->NIngridSimVertexes()>1) cout<<"********************************"<<endl;
      IngridSimVertexSummary * simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
      evt2->AddSimVertex(simver);
      IngridHitSummary * Hit=new IngridHitSummary();
      IngridHitSummary * Hit2=new IngridHitSummary();
      IngridSimHitSummary * SimHit=new IngridSimHitSummary();
      IngridSimHitSummary * SimHit2=new IngridSimHitSummary();
      IngridSimParticleSummary * SimPart=new IngridSimParticleSummary();
      //int mod;
      double posX;
      double posY;
      double posZ;
      double norm;
      double totcrsne;

      if(MC){
	Nu_E=simver->nuE;
	Enu=Nu_E;
	Num_Int=simver->inttype;
	//cout<<endl<<endl<<"Num Int="<<Num_Int<<endl;
	IntNumber++;
	//mod=simver->mod;
	 posX=simver->xnu;
	 posY=simver->ynu;
	 posZ=simver->znu;

	 norm=simver->norm;
	 totcrsne=simver->totcrsne;
	 weight = 1;//MCCorrections(1,mod);
	 if(a==1) weight*=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(1.03*46.2)/100;
	 else weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*470)*1.61/100;
      }
      evt2->run=evt->run;
      evt2->event=evt->run;
      evt2->runmode=evt->runmode;
      evt2->trgid=evt->trgid;
      evt2->version=evt->version;
      evt2->date=evt->date;
      evt2->time=evt->time;
      evt2->trgtime=evt->trgtime;
      evt2->nd280nspill=evt->nd280nspill;

      int nsimhit = evt -> NIngridSimHits();
      for(int i=0;i<nsimhit;i++){
	SimHit=evt->GetIngridSimHit(i);
	evt2->AddIngridSimHit(SimHit);
      }
      
      int nsimpart = evt -> NIngridSimParticles();
      for(int i=0;i<nsimpart;i++){
	SimPart=evt->GetSimParticle(i);
	evt2->AddSimParticle(SimPart);
      }
      
      for( int mod=0; mod<NMod; mod++ ){   //### Module Loop                                                                      
	for( int cyc=0; cyc<24; cyc++ ){  //### Cycle Loop                                                       
	  
	  int ninghit = evt -> NIngridModHits(mod, cyc);
	  NIngBasRec= evt->NPMAnas();
	  for(int i=0;i<ninghit;i++){
	    Hit2->Clear();
	    Hit->Clear();
	    Hit2=evt->GetIngridModHit(i,mod,cyc);
	    Hit->mod=Hit2->mod;
	    Hit->cyc=Hit2->cyc;
	    Hit->view=Hit2->view;
	    Hit->pln=Hit2->pln;                  // plane number
	    Hit->ch=Hit2->ch;                   // strip number(0~23)
	    Hit->adc=Hit2->adc;                  // high gain ADC value
	    Hit->loadc=Hit2->loadc;                // low  gain ADC value
	    Hit->pe=Hit2->pe;                 // number of photoelectrons, without correction
	    Hit->lope=Hit2->lope;               // number of photoelectrons, without correction
	    Hit->pecorr=Hit2->pecorr;             // number of photoelectrons, with correction
	    Hit->tdc=Hit2->tdc;                // raw TDC value
	    Hit->time=Hit2->time;               // hit time (ns).
	    Hit->tnearhit=Hit2->tnearhit;           // minumum value of hit time difference (ns).
	    Hit->timecorr=Hit2->timecorr;           // hit time (ns). Only the first hit inside 
	    Hit->xy=Hit2->xy;                 // transverse position (cm), x or y
	    Hit->z=Hit2->z;                  // position (cm) along beam direction
	    Hit->addbasicrecon=Hit2->addbasicrecon;      // This hit is member of basic recon or not
	    Hit->dummy=Hit2->dummy;              // this is dummy(study for MPPC noise)
	    Hit->gocosmic=Hit2->gocosmic;           // For efficiency study with cosmic, it is denominator for efficiency
	    Hit->hitcosmic=Hit2->hitcosmic;          // For efficiency study with cosmic, it is numerator for efficiency
	    Hit->isohit=Hit2->isohit;
	    
	    for(int j=0;j<Hit2->NSimHits();j++){
	      SimHit=Hit2->GetIngridSimHit(j);
	      nsimhit=evt2->NIngridSimHits();
	      for(int k=0;k<nsimhit;k++){
		SimHit2=evt2->GetIngridSimHit(k);
		if((SimHit->pdg==SimHit2->pdg) && (SimHit->trackid==SimHit2->trackid) && (SimHit->edeposit==SimHit2->edeposit)){Hit->AddIngridSimHit(SimHit2);}
	      }
	    }
	 	    
	    evt2->AddIngridModHit(Hit,mod,cyc);
	    /*
	    bool TrackHit=false;
	    for(int irec=0;irec<NIngBasRec;irec++){
	      PMAnaSummary * recon = (PMAnaSummary*) evt->GetPMAna(irec);
	      if(recon->hitcyc==cyc){
		nTracks=recon->Ntrack;
		for(int itrk=0;itrk<nTracks;itrk++){
		  for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
		    Hit=recon->GetIngridHitTrk(ihit,itrk);
		    if(Hit->mod == Hit2->mod && Hit->tdc==Hit2->tdc && Hit->pln==Hit2->pln && Hit->ch==Hit2->ch && Hit->view==Hit2->view){
		      //TrackHit=true;
		      int ibin= (int) ((recon->angle)[itrk]/3);
		      if(Hit->mod==16 && Hit->pln!=0 && Hit->ch>7 && Hit->ch<24){
			double Inef = TMath::Abs(InputEff_SciBar->GetBinContent(ibin)-InputEffMC_SciBar->GetBinContent(ibin));
			double IRand=Rand->Uniform();
			//cout<<Inef<<", Irand="<<IRand<<endl;
			//if(IRand<Inef) continue;
			//else evt2->AddIngridModHit(Hit2,mod,cyc);
		      }
		      else{
			double Inef = TMath::Abs(InputEff_Ingrid->GetBinContent(ibin)-InputEffMC_Ingrid->GetBinContent(ibin));
                        double IRand=Rand->Uniform();
                        //if(IRand<Inef) continue;
                        //else evt2->AddIngridModHit(Hit2,mod,cyc);
                      }

		    }
		  }
		}
	      }
	    }
	    if(TrackHit==false) evt2->AddIngridModHit(Hit2,mod,cyc);*/
	  }
	}
      }
      wtree->Fill();
    }
  }
  }
  wfile->cd();
  wtree->Write();
  wfile->Close();
  return(0);
}
  


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
//int LimitTracks=100;
int LimitRecs=10;
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
/* already defined in Reconstruction.cc
double DegRad(double angle){
  return angle*TMath::Pi()/180.;
}

double RadDeg(double angle){
  return angle*180./TMath::Pi();
}
*/

int main(int argc, char **argv)
{

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  bool PM=true;
  int c=-1;
  string InputFileName;string OutputFileName;

  while ((c = getopt(argc, argv, "i:o:w")) != -1) {
    switch(c){
    case 'i':
      InputFileName=optarg;
      break;
    case 'o':
      OutputFileName=optarg;
      break;
    case 'w':
      PM=false;
      break;
    }
  }
  
  Xsec * xs = new Xsec(PM);
  xs->Xsec::Initialize();
  Reconstruction * Rec = new Reconstruction(PM);
  Corrections * Cor = new Corrections(PM);
  INGRID_Dimension * IngDim = new INGRID_Dimension();
  char DetName[2]; sprintf(DetName,(PM?"PM":"WM"));

  cout<<"Welcome"<<endl;
  int NIngBasRec;

  TFile * DataHitEffiency = new TFile(Form("%s/XS/files/DataHitEfficiency_%s.root",cINSTALLREPOSITORY,DetName));
  TFile * MCHitEffiency = new TFile(Form("%s/XS/files/MCHitEfficiency_%s.root",cINSTALLREPOSITORY,DetName));

  TH1D * DataEfficiency_PMIng,*DataEfficiency_PMSci,*DataEfficiency_WMPlan,*DataEfficiency_WMGrid;
  TH1D * MCEfficiency_PMIng,*MCEfficiency_PMSci,*MCEfficiency_WMPlan,*MCEfficiency_WMGrid;

  if(PM){
    DataEfficiency_PMIng = (TH1D*) DataHitEffiency->Get("Efficiency_PMIng");
    DataEfficiency_PMSci = (TH1D*) DataHitEffiency->Get("Efficiency_PMSci");
    MCEfficiency_PMIng = (TH1D*) MCHitEffiency->Get("Efficiency_PMIng");
    MCEfficiency_PMSci = (TH1D*) MCHitEffiency->Get("Efficiency_PMSci");
  }
  else {
    DataEfficiency_WMPlan = (TH1D*) DataHitEffiency->Get("Efficiency_WMPlan");
    DataEfficiency_WMGrid = (TH1D*) DataHitEffiency->Get("Efficiency_WMGrid");
    MCEfficiency_WMPlan = (TH1D*) MCHitEffiency->Get("Efficiency_WMPlan");
    MCEfficiency_WMGrid = (TH1D*) MCHitEffiency->Get("Efficiency_WMGrid");
  }

  TRandom3 * Rand = new TRandom3();
  Rand->SetSeed(0);
  
  TFile*            wfile = new TFile((OutputFileName).c_str(), "recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* evt2 = new IngridEventSummary();
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary",
                                  &evt2,  64000,  99);
  
  cout<<"Opening Events"<<endl;


  TFile * _file0=new TFile((InputFileName).c_str());
  if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
  else{return 0;}
  TTree * tree=(TTree*) _file0->Get("tree");
  if(tree!=tree){return 0;}
  int nevt=(int) tree->GetEntries();
  cout<<"Total Number Of Events="<<nevt<<endl;

  IngridEventSummary* evt = new IngridEventSummary();
  tree->SetBranchAddress("fDefaultReco.",&evt);
  
  IngridSimVertexSummary * simver;
  IngridHitSummary * Hit=new IngridHitSummary();
  IngridHitSummary * Hit2;
  IngridHitSummary * HitTrk;
  IngridSimHitSummary * SimHit;
  IngridSimHitSummary * SimHit2;
  IngridSimParticleSummary * SimPart;
  PMAnaSummary * recon;
    
  for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
    if((ievt%100)==0) cout<<"Processing "<<ievt<<endl;
    evt->Clear();
    evt2->Clear();
    tree->GetEntry(ievt);
    if(evt->NIngridSimVertexes()>1) cout<<"********************************"<<endl;
    simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?
    evt2->AddSimVertex(simver);
        
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
    //Goal is to:
    //1. Take each hit of the event: Hit2
    //2. Copy it into another hit: Hit
    //3. If this hit is within the reconstruction (loop on reconstruction hits, and then compare their properties): apply the masking randomly. Depending on the random result, add it to the total event hit of a new event.
    //4. If the hit is not within track, add it directly to the new event
    //-> we should keep all hits not reconstructed, and only some of them that belongs to a reconstructed track. 
    for( int mod=0; mod<NMod; mod++ ){   //### Module Loop                                                                      
      for( int cyc=0; cyc<23; cyc++ ){  //### Cycle Loop                                                       
	
	int ninghit = evt -> NIngridModHits(mod, cyc);
	NIngBasRec= evt->NPMAnas();
	for(int i=0;i<ninghit;i++){
	  //Hit2->Clear();
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
	  // pecorr contains pe_cross in case of WM MC hits
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
	    //SimHit->Clear();
	    SimHit=Hit2->GetIngridSimHit(j);
	    nsimhit=evt2->NIngridSimHits();
	    for(int k=0;k<nsimhit;k++){
	      //SimHit2->Clear();
	      SimHit2=evt2->GetIngridSimHit(k);
	      if((SimHit->pdg==SimHit2->pdg) && (SimHit->trackid==SimHit2->trackid) && (SimHit->edeposit==SimHit2->edeposit)){Hit->AddIngridSimHit(SimHit2);}
	    }
	  }
	  
	  
	  bool TrackHit=false;
	  //Before adding the hit, one will test if it is in a track. 
	  //If yes, do the test, if not add it.

	  for(int irec=0;irec<NIngBasRec;irec++){
	    //recon->Clear();
	    recon = (PMAnaSummary*) evt->GetPMAna(irec);
	    if(recon->hitcyc!=cyc) continue;
	    int nTracks=recon->Ntrack;
	      
	    for(int itrk=0;itrk<nTracks;itrk++){
	      //cout<<recon->Ntrack<<", hitcyc="<<recon->hitcyc<<", nhit="<<recon->NhitTs(itrk)<<endl;
	      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
		//HitTrk->Clear();
		if(TrackHit) continue;
		HitTrk=recon->GetIngridHitTrk(ihit,itrk);
		//cout<<HitTrk->mod<<", "<<Hit2->mod<<endl;
		if(HitTrk->mod == Hit2->mod && HitTrk->tdc==Hit2->tdc && HitTrk->pln==Hit2->pln && HitTrk->ch==Hit2->ch && HitTrk->view==Hit2->view){
		  TrackHit=true;
		  int BinRecAngle;
		  for(int i=0;i<=NBinsRecAngle;i++){
		    if((recon->angle)[itrk]<BinningRecAngle[i+1]){BinRecAngle=i;break;}
		  }

		  // ML 2017/07/03 I use alternative computation of the hit inefficiency
		  // --- if d/m > 1 ie d > m, then nothing to mask in MC and Inef<0
		  double d=0,m=0;
		  if(PM){
		    if(!(Rec->Reconstruction::IsINGRID(Hit2->mod,Hit2->pln,Hit2->ch))){
		      d=DataEfficiency_PMSci->GetBinContent(BinRecAngle);
		      m=MCEfficiency_PMSci->GetBinContent(BinRecAngle);
		      //double Inef =TMath::Abs(DataEfficiency_PMSci->GetBinContent(BinRecAngle)-MCEfficiency_PMSci->GetBinContent(BinRecAngle));
		    }
		    else{
		      d=DataEfficiency_PMIng->GetBinContent(BinRecAngle);
		      m=MCEfficiency_PMIng->GetBinContent(BinRecAngle);
		      //  double Inef = TMath::Abs(DataEfficiency_PMIng->GetBinContent(BinRecAngle)-MCEfficiency_PMIng->GetBinContent(BinRecAngle));
		    }
		  }
		  else{
		    bool IsGrid=(Hit->ch>=40);
		    if(!IsGrid){
		      d=DataEfficiency_WMPlan->GetBinContent(BinRecAngle);
		      m=MCEfficiency_WMPlan->GetBinContent(BinRecAngle);
		      //  double Inef = TMath::Abs(DataEfficiency_WMPlan->GetBinContent(BinRecAngle)-MCEfficiency_WMPlan->GetBinContent(BinRecAngle));
		    }
		    else{
		      d=DataEfficiency_WMPlan->GetBinContent(BinRecAngle);
		      m=MCEfficiency_WMPlan->GetBinContent(BinRecAngle);
		      //double Inef = TMath::Abs(DataEfficiency_WMGrid->GetBinContent(BinRecAngle)-MCEfficiency_WMGrid->GetBinContent(BinRecAngle));
		    }
		  }//WM
		  
		  double Inef=1.-(d/m);
		  double IRand=Rand->Uniform();
		  if(IRand<Inef){cout<<"Hit removed"<<endl;continue;}
		  else evt2->AddIngridModHit(Hit,mod,cyc);

		}//same hit
	      }//hit
	    }//trk
	  }//rec   
	  if(TrackHit==false) evt2->AddIngridModHit(Hit,mod,cyc);

	}
      }
    }
    wtree->Fill();

  }
  wfile->cd();
  wtree->Write();
  wfile->Close();
  return(0);
}
  

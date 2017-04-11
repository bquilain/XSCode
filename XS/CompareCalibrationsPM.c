#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>

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
#include <TH2.h>
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
#include <TRandom3.h>
#include "PMdispRev.h"

#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
//#include "PMAnaSummary.h"
#include "IngridConstants.h"
#include "setup.h"
#include "Reconstruction.cc"
Reconstruction *Rec=new Reconstruction();
#include "Corrections.cc"
Corrections *Cor=new Corrections();

bool IsINGRID(int ch){
  bool Ing;
  if(ch<=7||ch>=24) Ing=true;
  else Ing=false;
  return Ing;
}
/* already defined in Reconstruction.cc
double DegRad(double angle){
  return angle*TMath::Pi()/180.;
}

double RadDeg(double angle){
  return angle*180./TMath::Pi();
}
*/

int main(int argc, char **argv){
  cout<<"hello"<<endl;
  char type;
  //cout<<"Data or MC?(d/m)"<<endl;
  //cin>>type;
  bool Disp=false;
  bool MC=false;
  int c=-1;
  int NFiles=30;
  int IFiles=0;
  int IRuns=14510;
 int FRuns=14510;
 char Name[256];

  while ((c = getopt(argc, argv, "mdi:f:r:t:")) != -1) {
    switch(c){
    case 'm':
      MC=true;
      break;
    case 'd':
      Disp=true;
      break;
    case 'i':
      IFiles=atoi(optarg);
      break;
    case 'f':
      NFiles=atoi(optarg);
      break;
    case 'r':
      IRuns=atoi(optarg);
      break;
    case 't':
      FRuns=atoi(optarg);
      break;
    }
  }
  cout<<NFiles<<endl;
  double Icyc=4;
  double Ncyc=12;
  double Nmod=17;
  TFile* _file0 = NULL;
  char FileMC[500], FileData[256];
  cout<<"Hello Matt :)"<<endl;
  if(MC) cout<<"Analyzing MC"<<endl;
  else cout<<"Analyzing Data"<<endl;
  double Nhits, NhitsMC;
  IngridHitSummary * Hit = new IngridHitSummary();
  //if(!MC) sprintf(Name,"/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/Data_CalibrationPM%d.root",FRuns);
  //else sprintf(Name,"/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/MC_CalibrationPM%d.root",NFiles);
  if(!MC) sprintf(Name,"/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/Data_CalibrationPM.root");
  else sprintf(Name,"/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/MC_CalibrationPM.root");
  TFile * wfile = new TFile(Name,"recreate");
  TTree * wtree = new TTree("wtree","wtree");

  int utime=0;
  double Charge=0;
  double ChargeCorrected=0;
  double ChargeDXCorrected=0;
  double dx=0;
  int Module=0;
  int Plane=0;
  int View=0;
  int Chan=0;
  int Channel=0;
  double X=0;double Y=0;
  double AngleX=0;double AngleY=0;
  double Angle=0;
  double weight=0;
  bool VSelectionFV=false;
  bool VSelectionOV=false;
  int TrackSample=0;
  int Used=0;
  int MCSample=-1;

  wtree              -> Branch   ("utime",&utime,"utime/I");
  wtree              -> Branch   ("SelectionFV",&VSelectionFV,"SelectionFV/O");
  wtree              -> Branch   ("SelectionOV",&VSelectionOV,"SelectionOV/O");
  wtree              -> Branch   ("Charge",&Charge,"Charge/D");
  wtree              -> Branch   ("ChargeCorrected",&ChargeCorrected,"ChargeCorrected/D");
  wtree              -> Branch   ("ChargeDXCorrected",&ChargeDXCorrected,"ChargeDXCorrected/D");
  wtree              -> Branch   ("dx",&dx,"dx/D");
  wtree              -> Branch   ("Channel",&Channel,"Channel/I");
  wtree              -> Branch   ("Module",&Module,"Module/I");
  wtree              -> Branch   ("Plane",&Plane,"Plane/I");
  wtree              -> Branch   ("View",&View,"View/I");
  wtree              -> Branch   ("Chan",&Chan,"Chan/I");
  wtree              -> Branch   ("X",&X,"X/D");
  wtree              -> Branch   ("Y",&Y,"Y/D");
  wtree              -> Branch   ("AngleX",&AngleX,"AngleX/D");
  wtree              -> Branch   ("AngleY",&AngleY,"AngleY/D");
  wtree              -> Branch   ("Angle",&Angle,"Angle/D");
  wtree              -> Branch   ("weight",&weight,"weight/D");
  wtree              -> Branch   ("TrackSample",&TrackSample,"TrackSample/I");
  wtree              -> Branch   ("Used",&Used,"Used/I");
  wtree              -> Branch   ("MCSample",&MCSample,"MCSample/I");
 

  TH1F * POTCount = new TH1F("POTCount","",2,0,2);

  TTree * tree;
  IngridEventSummary* evt;
  TBranch * Br;
  PMAnaSummary * pmana;
  BeamInfoSummary * BeamSummary = new BeamInfoSummary();
  IngridHitSummary * inghitsum = new IngridHitSummary();
  IngridSimVertexSummary * simver;

  for(int a=0;a<1;a++){
    cout<<"a="<<a<<endl;
    if(!MC && a!=0) continue;
    
    for(int R=IRuns;R<=FRuns;R++){
      for(int f=IFiles;f<=NFiles;f++){
	if(MC && (f>=415 && f<=419)) continue;
	
	if(MC) {
	  if(a==0){
	    sprintf(FileMC,"/export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root",f);
	  }
	  else if(a==1) sprintf(FileMC,"/export/scraid2/data/bquilain/MCfiles/PMMC_Run1_%d_wNoise_ana.root",f);
	  MCSample=a;
	  POTCount->Fill(1.,1e21);
	}
	else sprintf(FileMC,"/export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root",R,f);    
	
	_file0=new TFile(FileMC);
	
      if(_file0->IsOpen()) cout << _file0->GetName() << "is open"<<endl ;
      else continue;

      tree=(TTree*) _file0->Get("tree");
      int nevt=(int) tree->GetEntries();
      cout<<nevt<<endl;
      evt = new IngridEventSummary();
      Br=tree->GetBranch("fDefaultReco.");
      Br->SetAddress(&evt);
      tree->SetBranchAddress("fDefaultReco.",&evt);
          
      for(int ievt=0;ievt<nevt;ievt++){
	if((ievt%100)==0) cout<<ievt<<endl;
	evt->Clear();//vire l'evt précédent
	tree->GetEntry(ievt);//charge l'evt grace au link avec la branche

	utime=0;
	weight=1;
	int NPMAnas;
	NPMAnas = evt->NPMAnas();
	if(MC){
	  simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));//il y a un numéro. On peut donc bien avoir plusieurs simvert/periode d'integ ;-)?    
	  double Nu_E=simver->nuE;
	  int Num_Int=simver->inttype;
	  int mod=simver->mod;
	  double posX=simver->xnu;
	  double posY=simver->ynu;
	  double posZ=simver->znu;
	  double norm=simver->norm;
	  double totcrsne=simver->totcrsne;
	  if(a==0) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*470);
	  else if(a==1) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*46.2;
	}
	else{
	  for(int ib=0;ib<evt->NIngridBeamSummarys();ib++){
	    BeamSummary = evt->GetBeamSummary(ib);
	    for( int cyc=Scyc; cyc<Ncyc; cyc++ ) POTCount->Fill(1,BeamSummary->ct_np[4][cyc-4+1]);
	    utime  =  BeamSummary->trg_sec;
	  }
	}

	double PE;

	for(int irec=0;irec<NPMAnas;irec++){
	  pmana = (PMAnaSummary*) evt->GetPMAna(irec);
	  Rec->Reconstruction::GetSelectionPM(&VSelectionFV,&VSelectionOV,pmana,MC);
	 
	  int cyc=pmana->hitcyc;
	    int inihit=0;
	    double PeTrack=0;
	    vector< vector<Hit3D> > VecDouble;
	    for(int i=0;i<VecDouble.size();i++){
	      VecDouble[i].clear();
	    }
	    VecDouble.clear();
	    
	    for(int itrk=0;itrk<pmana->Ntrack;itrk++){
	      vector <Hit3D> VecT;
	      vector <HitTemp> HitV;
	      VecT.clear();
	      HitV.clear();
	      HitV=Rec->Reconstruction::EraseDoubleHitsPM(pmana,itrk,HitV);
	      VecT=Rec->Reconstruction::Hit2DMatchingPM(evt,pmana,HitV,VecT,MC);
	      VecDouble.push_back(VecT);
	    }

	    for(int itrk=0;itrk<pmana->Ntrack;itrk++){
	      bool Geom=Rec->Reconstruction::HasGeomTrack(16,(pmana->startxpln)[itrk],(pmana->startxch)[itrk],DegRad((pmana->thetax)[itrk]), (pmana->startypln)[itrk],(pmana->startych)[itrk],DegRad((pmana->thetay)[itrk]));
	      TrackSample=Rec->Reconstruction::SelectTrackSample((pmana->pm_stop)[itrk],Geom,(pmana->ing_trk)[itrk],(pmana->ing_stop)[itrk],(pmana->ing_endpln)[itrk]);
	      vector <HitTemp> HitV;
	      vector <Hit3D> Vec;
	      Vec.clear();
	      vector <double> PECorrectedFiberOnly;
	      PECorrectedFiberOnly.clear();
	      
	      double SlopeX=TMath::Tan((pmana->thetax)[itrk]*TMath::Pi()/180);
	      double SlopeY=TMath::Tan((pmana->thetay)[itrk]*TMath::Pi()/180);
	      //dx=TMath::Sqrt(SlopeY*SlopeY+SlopeX*SlopeX+1);
	      dx=1./TMath::Cos(DegRad((pmana->angle)[itrk]));
	      
	      HitV=Rec->Reconstruction::EraseDoubleHitsPM(pmana,itrk,HitV);
	      Vec=Rec->Reconstruction::Hit2DMatchingPM(evt,pmana,HitV,Vec,MC);
	      Vec=Rec->CountSharedHits(Vec,VecDouble,itrk);

	      if(Vec.size()==0){
		cout<<"No 3D hits"<<", number of hits="<<pmana->NhitTs(itrk)<<", hit in the event="<<endl;
		continue;
	      }
	      Vec=Cor->Corrections::GetFiberAttenuation(Vec);
       
	      for(int ihit=0;ihit<Vec.size();ihit++) PECorrectedFiberOnly.push_back(Vec[ihit].pecorr);
	      Vec=Cor->Corrections::GetDXCorrection(Vec,dx);

	      for(int ihit=0;ihit<Vec.size();ihit++){
		Charge=Vec[ihit].pe;
                ChargeCorrected=PECorrectedFiberOnly[ihit];
                ChargeDXCorrected=Vec[ihit].pecorr;
		SlopeX=TMath::Tan((pmana->thetax)[itrk]*TMath::Pi()/180);
		SlopeY=TMath::Tan((pmana->thetay)[itrk]*TMath::Pi()/180);
		//dx=TMath::Sqrt(SlopeY*SlopeY+SlopeX*SlopeX+1);
		Module=Vec[ihit].mod;
		Plane=Vec[ihit].pln;
		View=Vec[ihit].view;
		Chan=Vec[ihit].ch;
		X=Vec[ihit].x;
		Y=Vec[ihit].y;
		AngleX=(pmana->thetax)[itrk];
		Angle=(pmana->angle)[itrk];
		AngleY=(pmana->thetay)[itrk];
		Used=Vec[ihit].used;
		if(Vec[ihit].mod==16 && Rec->Reconstruction::IsINGRID(Vec[ihit].mod,Vec[ihit].pln,Vec[ihit].ch)) Channel=1;
		else if(Vec[ihit].mod==16 && !(Rec->Reconstruction::IsINGRID(Vec[ihit].mod,Vec[ihit].pln,Vec[ihit].ch))) Channel=2;
		else Channel=0;
		if(dx<0.5) cout<<"*******************************dx2="<<dx<<endl;
		
		wtree->Fill();
		//cout<<"channel="<<Channel<<", pe="<<ChargeCorrected<<endl;
	      }

	      //cout<<"end track"<<endl; 
	    }//trk
	}//PMAnas
      }//Events
      _file0->Close();
      }//Number of files
    }
  }//Sand or Beam

  wfile->cd();
  POTCount->Write("POTCount");
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

  //File.Close();

  return(0);
}
  

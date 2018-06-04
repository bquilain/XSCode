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


//#define DEBUG

int main(int argc, char **argv){

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");


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
  bool PM=true;

  while ((c = getopt(argc, argv, "wmdi:f:r:t:")) != -1) {
    switch(c){
    case 'm':
      MC=true;
      break;
    case 'w':
      PM=false; // WaterModule
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
  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));

  InitializeGlobal(PM);
  Rec->SetDetector(PM);
  Cor->SetDetector(PM);

  char * cMCOUT = getenv((PM?"MCOUTPUTSTORAGE":"MCOUTPUTSTORAGE_WM"));
  char * cDATAOUT = getenv((PM?"DATAOUTPUTSTORAGE":"DATAOUTPUTSTORAGE_WM"));

  cout<<"Hello Matt :)"<<endl;
  if(MC) cout<<"Analyzing MC in the "<<DetName<<endl;
  else cout<<"Analyzing Data in the "<<DetName<<endl;
  double Nhits, NhitsMC;
  IngridHitSummary * Hit = new IngridHitSummary();
  if(!MC) sprintf(Name,"%s/XS/files/Data_Calibration%s.root",cINSTALLREPOSITORY,DetName);
  else sprintf(Name,"%s/XS/files/MC_Calibration%s.root",cINSTALLREPOSITORY,DetName);
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
  PMAnaSummary * pmana;
  BeamInfoSummary * BeamSummary = new BeamInfoSummary();
  IngridHitSummary * inghitsum = new IngridHitSummary();
  IngridSimVertexSummary * simver;

    
  for(int R=IRuns;R<=FRuns;R++){
    for(int f=IFiles;f<=NFiles;f++){
      // ---- choice of the input files -----
      //if(MC && PM && (f>=415 && f<=419)) continue; ML 2017/06/23
      
      if(MC) {
	//sprintf(FileMC,"%s/%sMC_Wall_Run1_%d_wNoise_ana.root",cMCOUT,DetName,f);
	sprintf(FileMC,"%s/%sMC_Run1_%d_wNoise_ana.root",cMCOUT,DetName,f);
	MCSample=0;//0 is for Sand Muons
	POTCount->Fill(1.,1e21); // per file
      }
      else {
	if(PM) sprintf(FileMC,"%s/DataNew/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root",cDATAOUT,R,f);    
	else sprintf(FileMC,"%s/ingrid_%08d_%04d_anadev.root",cDATAOUT,R,f);    
      }
	
      _file0=new TFile(FileMC);
      
      if(_file0->IsOpen()) cout << _file0->GetName() << " is open"<<endl ;
      else continue;


      // ---- loading the content from the chosen file ----
      tree=(TTree*) _file0->Get("tree");
      int nevt=(int) tree->GetEntries();
      cout<<nevt<<endl;
      evt = new IngridEventSummary();
      tree->SetBranchAddress("fDefaultReco.",&evt);

      int n0=0;
          
      // ---- loop over the events ----
      for(int ievt=0;ievt<nevt;ievt++){
	if((ievt%1000)==0) cout<<ievt<<endl;
	evt->Clear();
	tree->GetEntry(ievt);

	utime=0;
	weight=1;
	int NPMAnas = evt->NPMAnas();
	if(MC){
	  simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));
	  double Nu_E=simver->nuE;
	  int Num_Int=simver->inttype;
	  int mod=simver->mod;
	  double posX=simver->xnu;
	  double posY=simver->ynu;
	  double posZ=simver->znu;
	  double norm=simver->norm;
	  double totcrsne=simver->totcrsne;
	  weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(2.2*470);
	  //else if(a==1) weight=norm*totcrsne*pow(10.,-38.)*6.02*pow(10.,23.)*(PM?46.2:50);//ML 2017/06/13
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

	  bool sandMuon=pmana->Ntrack==1 && (min(pmana->startxpln[0],pmana->startypln[0])<(PM?1:2)) && (max(pmana->endxpln[0],pmana->endypln[0])>(PM?15:20));
	  if(!sandMuon) continue;// new definition of sand muons
	  //	  if(pmana->Ntrack!=1 || !VSelectionOV) continue; // ML 2017/06/13 - old def of the sand muon sample
	  n0++;

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
	    bool Geom=Rec->Reconstruction::HasGeomTrack((PM?16:15),(pmana->startxpln)[itrk],(pmana->startxch)[itrk],DegRad((pmana->thetax)[itrk]), (pmana->startypln)[itrk],(pmana->startych)[itrk],DegRad((pmana->thetay)[itrk]));
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
	    Vec=Cor->Corrections::GetDXCorrection(Vec,dx); //only PM+INGRID here
	    if(!PM) Cor->Corrections::GetDXCorrectionWM(Vec,DegRad(pmana->angle[itrk]),DegRad(pmana->thetax[itrk]),DegRad(pmana->thetay[itrk]));
	    // output is pe/3mm

	    // now fill the variables and write the tree
#ifdef DEBUG
	    cout<<ievt<<" "<<n0<<" "<<Vec.size()<<" "<<pmana->angle[itrk]<<" "<<pmana->thetax[itrk]<<" "<<pmana->thetay[itrk]<<endl;
#endif
	    for(int ihit=0;ihit<Vec.size();ihit++){
	      Charge=Vec[ihit].pe;
#ifdef DEBUG
	      if(n0==12) cout<<ihit<<" "<<Vec[ihit].mod<<" "<<Vec[ihit].view<<" "<<Vec[ihit].pln<<" "<<Vec[ihit].ch<<" "<<Vec[ihit].pe;
#endif
	      if(Charge<4.5) continue; // ML 2017/06/20
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
	      AngleY=(pmana->thetay)[itrk];
	      
	      // generalized here to `normal angle` (unchanged for scinti from planes)
	      bool grid=(Module==15 && Chan>=40);
	      Angle=Rec->NormalAngle(pmana->angle[itrk],AngleX,AngleY,View,grid) ;
#ifdef DEBUG 
	      if(n0==12) cout<<" normal angle="<<Angle<<" dz="<<3*ChargeCorrected/ChargeDXCorrected<<","<<Cor->dzWM(DegRad(pmana->angle[itrk]),DegRad(View?AngleY:AngleX),grid)<<" pe="<<Charge<<endl;
#endif
	      Used=Vec[ihit].used;
	      if(Vec[ihit].mod==16){//PM
 		if(Rec->Reconstruction::IsINGRID(Vec[ihit].mod,Vec[ihit].pln,Vec[ihit].ch)) Channel=1;
		else Channel=2;
	      }
	      else if(Vec[ihit].mod==15){//WM
		if(!grid) Channel=1;
		else Channel=2+View; // Channel 2 is GridX and channel 3 is GridY
	      }
	      else Channel=0;//INGRID
	      if(PM && dx<0.5) cout<<"*******************************dx2="<<dx<<endl;
		
	      wtree->Fill();
	      //cout<<"channel="<<Channel<<", pe="<<ChargeCorrected<<endl;
	    }

	    //cout<<"end track"<<endl; 
	  }//trk
	}//PMAnas
      }//Events
      _file0->Close();
      cout<<"Number of sand selected = "<<n0<<endl;
    }//Number of files
  }


  wfile->cd();
  POTCount->Write("POTCount");
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

  //File.Close();

  return(0);
}
  

//
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>

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

#include "TApplication.h"
//#include "analysis_beam.hxx"
#include "PMrecon_HitInfo.hxx"
#include "IngAna_HitInfo.hxx"
#include "IngAnaVertex_HitInfo.hxx"
//#include "PMdisp.h"

#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "PMReconSummary.h"
#include "NeutInfoSummary.h"

const static int   GapbwBunch     = 581 ;    //[nsec]
const static int   TDCOffset      = 300;     //[nsec]
const static int   fExpBeamTime   = 200;     //[nsec]
const static int   beamontime     = 100;     //ontime


//#define  STUDY_BADCHANNEL
vector<int>  bad_mod;
vector<int>  bad_pln;
vector<int>  bad_view;
vector<int>  bad_ch;
void Add_Bad_Channel();
bool Is_Bad_Channel(IngridHitSummary* thit);


const static char *dst_data          = "/home/t2kingrid/data/dst";
const static char *cosmic_data       = "/home/t2kingrid/data/cosmic";



long fcTime;
int main(int argc,char *argv[]){

#ifdef STUDY_BADCHANNEL
  cout << "It is BAD channel study mode. Is it O.K.?" << endl;
  //cin.get();
  Add_Bad_Channel();
#endif

  vector <int> HitV;
  char buff1[200];//,buff2[200];
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);
  Int_t run_number;
  Int_t sub_run_number=0;
  Int_t c=-1;
  char Output[300];
  bool rename = false;
  bool renamef= false;
  Int_t Scyc =  4;
  Int_t Ncyc = 12;
  Int_t Nmod = 17;
  Int_t Nini = 0;
  bool disp = false; 
  bool cosmic = false;

  Int_t Nactpln=0;
  Int_t Layerpe=0;

  while ((c = getopt(argc, argv, "r:s:f:cado:i:")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 'o':
      sprintf(Output,"%s",optarg);
      rename  = true;
      break;
    case 's':
      sub_run_number=atoi(optarg);
      break;
    case 'f':
      sprintf(buff1,"%s",optarg);
      run_number=0;
      renamef = true;
      break;
    case 'a':
      Scyc = 0;
      Ncyc = 23;
      break;
    case 'c':
      Scyc = 14;
      Ncyc = 16;
      cosmic = true;
      break;
    case 'd':
      disp = true;
      break;
    case 'i':
      Nini=atoi(optarg);
      break;

    }
  }
  FileStat_t fs;
  // ifstream timing;

  //sprintf(buff1,"%s/ingrid_%08d_%04d_processed.root",
  //  dst_data, run_number,sub_run_number);
  if(!renamef){
    sprintf(buff1,"%s/ingrid_%08d_%04d_recon.root",
	    dst_data, run_number,sub_run_number);
  }
  if(cosmic){
    sprintf(buff1,"%s/ingrid_%08d_%04d_recon.root",
	    cosmic_data, run_number,sub_run_number);

  }
  
  //######## read root file #############################
  //#####################################################
  cout<<"reading "<<buff1<<"....."<<endl;
  if(gSystem->GetPathInfo(buff1,fs)){
    cout<<"Cannot open file "<<buff1<<endl;
    exit(1);
  }
  cout << buff1 << endl;
  IngridEventSummary* evt = new IngridEventSummary();
  //TFile*            rfile = new TFile(Form("/home/daq/data/dst/ingrid_%08d_%04d_processed.root", 
  TFile*            rfile = new TFile(buff1,
				      "read");
  TTree*             tree = (TTree*)rfile -> Get("tree");
  TBranch*          EvtBr = tree->GetBranch("fDefaultReco.");
  EvtBr                   ->SetAddress(&evt);
  tree                    ->SetBranchAddress("fDefaultReco.", &evt);
  int                nevt = (int)tree -> GetEntries();
  cout << "Total # of events = " << nevt <<endl;

  //#### make rootfile after analysis #####
  //#######################################
  //sprintf(buff1,"/home/daq/data/dst/ingrid_%08d_0000_tclster.root",run_number);
  if( !rename ){
    sprintf(Output, "%s/ingrid_%08d_%04d_track.root", 
	    dst_data, run_number, sub_run_number); 
  }
  if(cosmic){
    sprintf(Output, "%s/ingrid_%08d_%04d_track.root", 
	    cosmic_data, run_number, sub_run_number); 
  }

  TFile*            wfile = new TFile(Output, "recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* wsummary = new IngridEventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary", 
				  &wsummary,  64000,  99);
  IngridHitSummary*    inghitsum;
  PMReconSummary* recon;
  //PMReconSummary* recon = new PMReconSummary();
  IngridBasicReconSummary* ingrecon = new IngridBasicReconSummary();
  //Hit                        hit;
  Track                        track;
  for(int ievt=Nini; ievt<nevt; ievt++){

    if(ievt%100==0)cout << "analyze event# " << ievt<<endl;
    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);

    for( int mod=0; mod<Nmod; mod++ ){   //### Module Loop
      for( int cyc=Scyc; cyc<Ncyc; cyc++ ){  //### Cycle Loop
	for(int i=0;i<htrack.size();i++){
	  htrack[i].clear();
	}
	for(int i=0;i<vtrack.size();i++){
          vtrack[i].clear();
	  //cout<<vtrack[i].hitid.size()<<endl;
        }
	htrack.clear();
	vtrack.clear();
	//cout<<"Module = "<<mod<<" , Cycle = "<<cyc<<endl;
	for( int view=0; view<2; view++ ){  //### View
	  //allhit.clear();
	  //int ninghit = evt -> NIngridModHits(mod, cyc);
	  int ningtrack = evt -> NPMModRecons(mod, cyc, view);
	  //if(ningtrack!=0) cout<<ningtrack<<endl;
	  if(ningtrack==0)continue;

	  for(int i=0; i<ningtrack; i++){
	    track.clear();

	    recon   = (PMReconSummary*) (evt -> GetPMModRecon(i, mod, cyc, view) );

	    if(i==0){
	      Nactpln = recon -> nactpln;
	      Layerpe = recon -> layerpe;
	    }

	    track.view = recon -> view;
	    track.ipln = recon -> startpln;
            track.fpln = recon -> endpln;
            track.ixy  = recon -> startxy;
            track.fxy  = recon -> endxy;
            track.iz   = recon -> startz;
	    //cout<<recon -> startz<<endl;
            track.fz   = recon -> endz;
            track.ang  = recon -> angle;
            track.slope= recon -> slope;
            track.intcpt = recon -> intcpt;
            track.veto = recon -> vetowtracking;
            track.edge = recon -> edgewtracking;
            track.vetodist = recon -> vetodist;
            track.clstime = recon -> clstime;
	    //cout<<recon->HitId.size()<<endl;
	    //cout<<"In the first 2D tracks from PMrecon, hits = ";
	    //cout<<"Hit ID size="<<recon->HitId.size()<<endl;
	    for(int ihit=0;ihit<recon->HitId.size();ihit++){
	      track.hitid.push_back(recon -> HitId[ihit]); 
	      //cout<<track.hitid[ihit]<<"  ";
	    }
	    //cout<<endl;
	    //cout<<"Track 2D size = "<< track.hitid.size()<<endl;
	    //cout<<recon->hello<<endl;
	    //track.nhits
	    if(view==0)
	      htrack.push_back(track);
	    else
	      vtrack.push_back(track);
	    //cout<<"Track = "<<htrack[i].hitid.size()<<endl;
	  }
	}
	//cout<<"before 3D matching"<<endl;
	fIngAna(mod);
	//cout<<"before Vertex Recons"<<endl;
	fIngVertex(mod);
	//cout<<"Vertex Recons done"<<endl;
	for(int i=0;i<ingtrack.size();i++){
	  //cout<<"before track "<<ingtrack.size()<<"  "<<i<<endl;
	  if(!(ingtrack[i].vtxtrk))continue;
	  //cout<<i<<endl;
	  ingrecon->Clear();
	  //cout<<i<<endl;
	  //HitV.clear();
	  ingrecon -> startxpln = ingtrack[i].startxpln;
	  ingrecon -> startypln = ingtrack[i].startypln;
	  ingrecon -> startxch  = ingtrack[i].startxch;
	  ingrecon -> startych  = ingtrack[i].startych;
	  ingrecon -> endxpln   = ingtrack[i].endxpln;
	  ingrecon -> endypln   = ingtrack[i].endypln;
	  ingrecon -> endxch    = ingtrack[i].endxch;
	  ingrecon -> endych    = ingtrack[i].endych;
	  ingrecon -> vertexz   = ingtrack[i].vertexz;
	  ingrecon -> vertexxz  = ingtrack[i].vertexxz;
	  ingrecon -> vertexyz  = ingtrack[i].vertexyz;
	  ingrecon -> thetax    = ingtrack[i].thetax;
	  ingrecon -> thetay    = ingtrack[i].thetay;
	  ingrecon -> angle     = ingtrack[i].angle;
	  ingrecon -> x         = ingtrack[i].x;
	  ingrecon -> y         = ingtrack[i].y;
	  ingrecon -> z         = ingtrack[i].z;
	  ingrecon -> zx        = ingtrack[i].zx;
	  ingrecon -> zy        = ingtrack[i].zy;
	  ingrecon -> vetowtracking = ingtrack[i].vetowtracking;
	  ingrecon -> edgewtracking = ingtrack[i].edgewtracking;
	  ingrecon -> vetodist = ingtrack[i].vetodist;
	  ingrecon -> Ntrack = ingtrack[i].Ntrack;
	  ingrecon -> hastrk = true;
	  ingrecon -> matchtrk = true;
	  ingrecon -> hitcyc = cyc;
	  ingrecon -> hitmod = mod;
	  ingrecon -> clstime = ingtrack[i].clstime;
	  ingrecon -> exptime = ((int)ingrecon -> clstime - TDCOffset )% GapbwBunch - fExpBeamTime;
	  //cout<<ingrecon -> Ntrack<<endl;
	  for(int itrk=0;itrk<ingrecon->Ntrack;itrk++){
	    HitV.clear();
	    //cout<<"Track Num="<<itrk<<endl;
	    //cout<<"number of hits="<<ingtrack[i].hitnum[itrk].size()<<endl;
	    for(int ihit=0;ihit<ingtrack[i].hitnum[itrk].size();ihit++){
	      inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(ingtrack[i].hitnum[itrk][ihit], mod, cyc) );
	      //if(ingrecon->Ntrack>1) cout<<"Hit identifier = "<<ingtrack[i].hitnum[itrk][ihit]<<endl;
	      //cout<<"before adding the hit"<<endl;
	      HitV.push_back(ingtrack[i].hitnum[itrk][ihit]);
	      
	      int twice=0;
	      //cout<<"Number of hits="<<HitV.size()<<endl;
	      for(int a=0;a<HitV.size()-1;a++){
		if(HitV[a]==ingtrack[i].hitnum[itrk][ihit]) {
		  //cout<<"twice"<<endl;
		  twice++;
		  break;
		}
	      }
	    
	    if(twice!=0) {
	      HitV.pop_back();
	      continue;
	    }
	    else ingrecon->AddIngridHitTrk(inghitsum,itrk); 
	    /*
	    for(int a=0;a<HitV.size();a++){
	      cout<<HitV[a]<<endl;
	    }
	    cout<<endl;*/
	    //cout<<"Track i ="<<i<<"Size of the first= "<<ingtrack[i].hitnum[0].size()<<"  "<<ingrecon ->NhitTs(itrk)<<endl;
	    }
	  }



	  if(fabs(ingrecon -> exptime) < beamontime)
	    ingrecon -> ontime = true;
	  else
	    ingrecon -> ontime = false;

	  ingrecon -> nactpln = Nactpln;
	  ingrecon -> layerpe = Layerpe;

	  evt   -> AddBasicRecon( ingrecon );
	  ingtrack[i].clear();
	}
	//cout<<"event is added"<<endl;

	
	
      }//cyc
    }//mod
    
    //fAnalysis(evt,Scyc,Ncyc);
    
    
    wsummary = evt;
    wtree -> Fill();
  }//Event Loop
  
  //######## Write and Close ####################
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();
  
  }

const int bad_num = 20;
int       bad_id[bad_num][4] =
  {// mod, view,  pln,  hh, 
    {    0,    1,    8,    2}, //cable damage channel 
    {    0,    1,    8,   14}, //cable damage channel
    {    0,    0,   13,   12}, //cable damage channel
    {    1,    1,    1,   14}, //cable damage channel
    {    5,    1,    2,    0}, //cable damage channel
    {    5,    1,    4,    2}, //high gain channel  
    {    5,    1,    4,   18}, //cable damage channel
    {    5,    1,    4,   20}, //cable damage channel
    {    5,    1,    7,    9}, //cable damage channel
    {    7,    0,    5,    0}, //pedestal channel    
    {    7,    0,    6,    9}, //pedestal channel   
    {    9,    1,    0,   12}, //cable damage chnnale
    {    9,    0,    4,   11}, //cable damage chnnale
    {   11,    0,    5,   13}, //pedestal shift channel
    {   15,    1,    8,   14}, //added 2010/11/1 cable damaged
    {   16,    1,    9,   21}, //added 2010/12/4 ??                
    {   16,    0,   14,    0}, //added 2010/12/4 ??               
    {   14,    1,   10,    4}, //added 2010/12/4 ??               
    {   14,    0,   10,   15}, //added 2010/12/4 ??               
    {   15,    0,   10,   15}  //added 2010/12/4 ??               
    //{    3,    1,    1,   11}  //added 2010/12/4 ??               
  };

void Add_Bad_Channel(){
  bad_mod. clear();
  bad_view.clear();
  bad_pln. clear();
  bad_ch.  clear();
  for(int ibad=0; ibad < bad_num; ibad++){
    bad_mod.  push_back ( bad_id[ibad][0] );
    bad_view. push_back ( bad_id[ibad][1] );
    bad_pln.  push_back ( bad_id[ibad][2] );
    bad_ch.   push_back ( bad_id[ibad][3] );
  }
}

bool Is_Bad_Channel(IngridHitSummary* thit){
  for(int ibad=0; ibad < bad_num; ibad++){
    if( thit -> mod  == bad_id[ibad][0]  &&
	thit -> view == bad_id[ibad][1]  &&
	thit -> pln  == bad_id[ibad][2]  &&
	thit -> ch   == bad_id[ibad][3]
	)
      return true;
  }
  return false;
}


//#define DEBUG_INGHIT
//#define INGRIDPID
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
#include "PMreconRevOfficial.hxx"
#include "PMAnaRevOfficial.hxx"

#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "NeutInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "PMReconSummary.h"
#include "PMAnaSummary.h"

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
void GetIngNonRecHits(IngridEventSummary* evt, int cyc); // for fIngHitPMJoint

const static char *dst_data          = "/home/kikawa/scraid1/data/pm_ingrid";
const static char *cosmic_data       = "/home/kikawa/scraid1/data/pm_ingrid";
//static const char *data_file_folder  = "/home/kikawa/scraid1/data/pm_ingrid";
//static const char *calib_file_folder = "/home/kikawa/scraid1/data/pm_ingrid";
//static const char *dailycheck_folder = "/home/kikawa/scraid1/data/pm_ingrid";



//#include "IngridConstants.h"
long fcTime;
int main(int argc,char *argv[]){

#ifdef STUDY_BADCHANNEL
  cout << "It is BAD channel study mode. Is it O.K.?" << endl;
  //cin.get();
  Add_Bad_Channel();
#endif

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
  int VertexingPlane  =pln_th;
  double VertexingChannel =ch_th;//mm
  int TrackMatching=diff_th;
  double AngleCut=ang_th;//degrees
  double TransverseCut=pos_th;//mm
  int nINGRIDPlanes=2; // number of planes fIngridHitPMJoint is looking at
  bool Error=false;
  int ErrorType=0;
  float ErrorValue=0.;
  requireIngridTrack=true;

  while ((c = getopt(argc, argv, "r:s:f:cado:i:e:v:N")) != -1) {
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
      requireIngridTrack=false;
      break;
    case 'd':
      disp = true;
      break;
    case 'i':
      Nini=atoi(optarg);
      break;
    case 'e':
      Error=true;
      ErrorType=atoi(optarg);
      break;
    case 'v':
      ErrorValue=atof(optarg);
      break;
    case 'N':
      requireIngridTrack=false;
      break;
    }
  }
  if(Error){
    if(ErrorType==11) VertexingPlane = (int) ErrorValue;
    else if(ErrorType==12) VertexingChannel = (double) ErrorValue;
    else if(ErrorType==13) TrackMatching = (int) ErrorValue;
    else if(ErrorType==14) AngleCut = (double) ErrorValue;
    else if(ErrorType==15) TransverseCut = (double) ErrorValue;
  }
  cout<<"Matching plane="<<TrackMatching<<", Vertexing plane="<<VertexingPlane<<", Channel="<<VertexingChannel<<"cm, Angle cut="<<AngleCut<<"°, Transverse cut="<<TransverseCut<<"cm"<<endl;
  cout<<"Fonction fIngHitPMJoint() looking for hits in first "<<nINGRIDPlanes<<" of INGRID"<<endl;
  cout<<(requireIngridTrack? "1 Ingrid track required" : "no Ingrid track required")<<endl;
  FileStat_t fs;
  // ifstream timing;

  //sprintf(buff1,"%s/ingrid_%08d_%04d_processed.root",
  //  dst_data, run_number,sub_run_number);
  if(!renamef){
    sprintf(buff1,"%s/ingrid_%08d_%04d_recon.root",
	    dst_data, run_number,sub_run_number);
  }
  /*  if(cosmic){
    sprintf(buff1,"%s/ingrid_%08d_%04d_recon.root",
	    cosmic_data, run_number,sub_run_number);
  }
  */
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
    sprintf(Output, "%s/ingrid_%08d_%04d_anas.root", 
	    dst_data, run_number, sub_run_number); 
  }
  /*
  if(cosmic){
    sprintf(Output, "%s/ingrid_%08d_%04d_anas.root", 
	    cosmic_data, run_number, sub_run_number); 
  }
  */
  TFile*            wfile = new TFile(Output, "recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* wsummary = new IngridEventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary", 
				  &wsummary,  64000,  99);
  IngridHitSummary*    inghitsum;
  IngridHitSummary* pmhit = new IngridHitSummary();
  PMReconSummary* recon;
  //PMReconSummary* recon = new PMReconSummary();
  PMAnaSummary* pmanasum = new PMAnaSummary();
  //Hit                        hit;
  TrackPM                      track;
  TrackIng                     ingtrack;

  Hits                        hit;
  
  

  for(int ievt=Nini; ievt<(Nini==0?nevt:Nini+1); ievt++){

    if(ievt%100==0)cout << "analyze event# " << ievt<<endl;
    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);

    
    for( int cyc=Scyc; cyc<Ncyc; cyc++ ){  //### Cycle Loop

      hingtrack.clear();
      vingtrack.clear();
      for( int mod=0; mod<7; mod++ ){  //### Module
	for( int view=0; view<2; view++ ){  //### View
	  int npmtrack = evt -> NPMModRecons(mod, cyc, view);
	  for(int i=0; i<npmtrack; i++){
	    recon   = (PMReconSummary*) (evt -> GetPMModRecon(i, mod, cyc, view) );
	    ingtrack.clear();

	    ingtrack.mod  = mod;
	    ingtrack.view = recon -> view;
	    ingtrack.ipln = recon -> startpln;
	    ingtrack.fpln = recon -> endpln;
	    ingtrack.iz   = recon -> startz;
	    ingtrack.fz   = recon -> endz;
	    ingtrack.ang  = recon -> angle;
	    ingtrack.slope= recon -> slope;
	    ingtrack.veto = recon -> vetowtracking;
	    ingtrack.edge = recon -> edgewtracking;
	    ingtrack.stop = recon -> modfc;
	    ingtrack.clstime = recon -> clstime;

	    if(view==0){
	      ingtrack.ixy  = recon -> startxy;
	      ingtrack.fxy  = recon -> endxy;
	      ingtrack.intcpt = recon -> intcpt;
	    }
	    else{
	      ingtrack.ixy  = recon -> startxy  + (mod-3)*1500;
	      ingtrack.fxy  = recon -> endxy    + (mod-3)*1500;
	      ingtrack.intcpt = recon -> intcpt + (mod-3)*1500;
	    }

	    for(int NHIT=0;NHIT<(recon->Nhits());NHIT++){
	      inghitsum = recon -> GetIngridHit(NHIT);
	      //cout<<"mod="<<inghitsum->mod<<endl;
	      hit.mod   = inghitsum->mod;
	      hit.view  = inghitsum->view;
	      hit.pln   = inghitsum->pln;
	      hit.ch    = inghitsum->ch;
	      hit.pe    = inghitsum->pe;
	      hit.lope  = inghitsum->lope;
	      hit.isohit= inghitsum->isohit;
	      hit.xy    = inghitsum->xy;
	      hit.z    = inghitsum->z;
	      hit.time    = inghitsum->time;
	      //if((inghitsum -> NSimHits()) > 0 && inghitsum -> cyc >= 0)
	      if((inghitsum -> NSimHits()) > 0)
     		hit.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	      else
		hit.pdg   = 0;
	      //cout<<hit.pdg<<endl;
	      ingtrack.hit.push_back(hit);

	    }

	    if(view==0)
	      hingtrack.push_back(ingtrack);
	    else
	      vingtrack.push_back(ingtrack);

	  }
	}
      }


      htrack.clear();
      vtrack.clear();
      for( int view=0; view<2; view++ ){  //### View
	int npmtrack = evt -> NPMModRecons(16, cyc, view);
	for(int i=0; i<npmtrack; i++){
	  recon   = (PMReconSummary*) (evt -> GetPMModRecon(i, 16, cyc, view) );
	  track.clear();

	  track.view = recon -> view;
	  track.ipln = recon -> startpln;
	  track.fpln = recon -> endpln;
	  track.ixy  = recon -> startxy;
	  track.fxy  = recon -> endxy;
	  track.iz   = recon -> startz;
	  track.fz   = recon -> endz;
	  track.ang  = recon -> angle;
	  track.slope= recon -> slope;
	  track.intcpt = recon -> intcpt;
	  track.veto = recon -> vetowtracking;
	  track.edge = recon -> edgewtracking;
	  track.stop = recon -> modfc;
	  track.clstime = recon -> clstime;

	  track.ing_trk = false;

	  for(int NHIT=0;NHIT<(recon->Nhits());NHIT++){
	    inghitsum = recon -> GetIngridHit(NHIT);
	    
	    hit.mod   = inghitsum->mod;
	    hit.view  = inghitsum->view;
	    hit.pln   = inghitsum->pln;
	    hit.ch    = inghitsum->ch;
	    hit.pe    = inghitsum->pe;
	    hit.lope  = inghitsum->lope;
	    hit.isohit= inghitsum->isohit;
	    hit.xy    = inghitsum->xy;
	    hit.z    = inghitsum->z;
	    hit.time    = inghitsum->time;

	    //if((inghitsum -> NSimHits()) > 0 && inghitsum -> cyc >= 0)
	    if((inghitsum -> NSimHits()) > 0)
	      hit.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	    else
	      hit.pdg   = 0;
	    
	    track.hit.push_back(hit);
	    
	  }

	  if(view==0)
	    htrack.push_back(track);
	  else
	    vtrack.push_back(track);
	}
      }

      GetIngNonRecHits(evt,cyc);//ML new


      //cout<<endl;
      if(!fPMAna(TrackMatching,VertexingPlane,VertexingChannel,AngleCut,TransverseCut, nINGRIDPlanes))continue;

      // isoHitCut is the minimal number of isolated it I require
      int isoHitCut=3;

      for(int i=0;i<pmtrack.size();i++){

	if(pmtrack[i].Ntrack == 0)continue;
	pmanasum -> Clear();

	reCalcIsohit(pmtrack[i],isoHitCut);

	//Vertex activity
	int ninghit = evt -> NIngridModHits(16, cyc);
	int xpln  = pmtrack[i].trk[0].startxpln;
	int ypln  = pmtrack[i].trk[0].startypln;
	float xch = pmtrack[i].trk[0].startxch;
	float ych = pmtrack[i].trk[0].startych;
	float vact=0;
	for(int ihit=0; ihit<ninghit; ihit++){
	  inghitsum  = (IngridHitSummary*) (evt -> GetIngridModHit(ihit,16,cyc));
	  if(Is_Bad_Channel( inghitsum ))continue;
	  vact += veract(inghitsum,xpln,ypln,xch,ych);
	}	

	pmanasum -> Ntrack = pmtrack[i].Ntrack;
	pmanasum -> Ningtrack = pmtrack[i].Ningtrack;
	pmanasum -> vetowtracking = pmtrack[i].vetowtracking;
	pmanasum -> edgewtracking = pmtrack[i].edgewtracking;
	pmanasum -> hitcyc = cyc;
	pmanasum -> clstime = pmtrack[i].clstime;
	pmanasum -> exptime = ((int)pmanasum -> clstime - TDCOffset )% GapbwBunch - fExpBeamTime;
	if(fabs(pmanasum -> exptime) < beamontime)
	  pmanasum -> ontime = true;
	else
	  pmanasum -> ontime = false;
	//pmanasum -> veract = vact;
	pmanasum -> clstimecorr = vact;//temporary clstimecorr is used to fill the vertex activity information

	pmanasum -> x           .clear();
	pmanasum -> y           .clear();
	pmanasum -> z           .clear();
	pmanasum -> zx          .clear();
	pmanasum -> zy          .clear();
	
	pmanasum -> startxpln   .clear();
	pmanasum -> startypln   .clear();
	pmanasum -> startxch    .clear();
	pmanasum -> startych    .clear();
	pmanasum -> endxpln     .clear();
	pmanasum -> endypln     .clear();
	pmanasum -> endxch      .clear();
	pmanasum -> endych      .clear();
	pmanasum -> thetax      .clear();
	pmanasum -> thetay      .clear();
	pmanasum -> angle       .clear();
	
	pmanasum -> ing_startmod.clear();
	pmanasum -> ing_endmod  .clear();
	pmanasum -> ing_startpln.clear();
	pmanasum -> ing_endpln  .clear();
	pmanasum -> ing_trk     .clear();
	pmanasum -> pm_stop     .clear();
	pmanasum -> ing_stop    .clear();
	pmanasum -> sci_range   .clear();
	pmanasum -> iron_range  .clear();
	pmanasum -> iron_pene   .clear();
	
	pmanasum -> veto        .clear();
	pmanasum -> edge        .clear();
	pmanasum -> pdg         .clear();
	pmanasum -> trkpe       .clear();
	pmanasum -> mucl        .clear();
	

	for(int t=0;t<pmtrack[i].trk.size();t++){//run over the 3d tracks found
	  for(int ihit=0;ihit<pmtrack[i].trk[t].hit.size();ihit++){//look at the hit registered
	    bool HitFound=false;
	    for(int imod=0;imod<Nmod;imod++){//vary the module number
	      int npmtrack = evt -> NPMModRecons(imod, cyc, pmtrack[i].trk[t].hit[ihit].view);//get the number of reconstruction in this module
	      for(int l=0; l<npmtrack; l++){//run over the reconstructions
		if(HitFound) break;
		recon   = (PMReconSummary*) (evt -> GetPMModRecon(l, imod, cyc, pmtrack[i].trk[t].hit[ihit].view) );
		for(int NHIT=0;NHIT<(recon->Nhits());NHIT++){//run over all hits of the recon
		  if(HitFound) break;
		  pmhit = recon -> GetIngridHit(NHIT);
		  if(pmhit->mod==pmtrack[i].trk[t].hit[ihit].mod && pmhit->view==pmtrack[i].trk[t].hit[ihit].view && pmhit->pln==pmtrack[i].trk[t].hit[ihit].pln && pmhit->ch==pmtrack[i].trk[t].hit[ihit].ch && pmhit->pe==pmtrack[i].trk[t].hit[ihit].pe && pmhit->lope==pmtrack[i].trk[t].hit[ihit].lope && pmhit->xy==pmtrack[i].trk[t].hit[ihit].xy && pmhit->z==pmtrack[i].trk[t].hit[ihit].z && pmhit->time==pmtrack[i].trk[t].hit[ihit].time){
		    pmhit->isohit=pmtrack[i].trk[t].hit[ihit].isohit; // correct the isohit variable
		    pmanasum -> AddIngridHitTrk(pmhit,t);//for each hit of the track, we check if it corresponds to a hit of the reconstruction. If yes, we add it ones, and then, we change the hit of the track we are looking for.
		    HitFound=true;
		  }
		}
	      } 
	    }
	    if(!HitFound){ // ie was not in any reconstruction -> from fIngHitPMJoint
	      pmhit = evt->GetIngridModHit(pmtrack[i].trk[t].hit[ihit].hit_id,pmtrack[i].trk[t].hit[ihit].mod,cyc);
#ifdef DEBUG_INGHIT
	      cout<<"trk="<<t<<" hit= mod "<<pmtrack[i].trk[t].hit[ihit].mod<<" view "<<pmtrack[i].trk[t].hit[ihit].view<<" pln "<<pmtrack[i].trk[t].hit[ihit].pln<<" ch "<<pmtrack[i].trk[t].hit[ihit].ch<<" added!"<<endl;
#endif
	      pmanasum->AddIngridHitTrk(pmhit,t);
	    }

	  }

	  pmanasum -> x           .push_back(pmtrack[i].trk[t].x);
	  pmanasum -> y           .push_back(pmtrack[i].trk[t].y);
	  pmanasum -> z           .push_back(pmtrack[i].trk[t].z);
	  pmanasum -> zx          .push_back(pmtrack[i].trk[t].intcptx);
	  pmanasum -> zy          .push_back(pmtrack[i].trk[t].intcpty);

	  pmanasum -> startxpln   .push_back(pmtrack[i].trk[t].startxpln);
	  pmanasum -> startypln   .push_back(pmtrack[i].trk[t].startypln);
	  pmanasum -> startxch    .push_back(pmtrack[i].trk[t].startxch);
	  pmanasum -> startych    .push_back(pmtrack[i].trk[t].startych);
	  pmanasum -> endxpln     .push_back(pmtrack[i].trk[t].endxpln);
	  pmanasum -> endypln     .push_back(pmtrack[i].trk[t].endypln);
	  pmanasum -> endxch      .push_back(pmtrack[i].trk[t].endxch);
	  pmanasum -> endych      .push_back(pmtrack[i].trk[t].endych);
	  pmanasum -> thetax      .push_back(pmtrack[i].trk[t].thetax);
	  pmanasum -> thetay      .push_back(pmtrack[i].trk[t].thetay);
	  pmanasum -> angle       .push_back(pmtrack[i].trk[t].angle);

	  pmanasum -> ing_startmod.push_back(pmtrack[i].trk[t].ing_startmod);
	  pmanasum -> ing_endmod  .push_back(pmtrack[i].trk[t].ing_endmod);
	  pmanasum -> ing_startpln.push_back(pmtrack[i].trk[t].ing_startpln);
	  pmanasum -> ing_endpln  .push_back(pmtrack[i].trk[t].ing_endpln);
	  pmanasum -> ing_trk     .push_back(pmtrack[i].trk[t].ing_trk);
	  pmanasum -> pm_stop     .push_back(pmtrack[i].trk[t].pm_stop);
	  pmanasum -> ing_stop    .push_back(pmtrack[i].trk[t].ing_stop);
	  pmanasum -> sci_range   .push_back(pmtrack[i].trk[t].sci_range);
	  pmanasum -> iron_range  .push_back(pmtrack[i].trk[t].iron_range);
	  pmanasum -> iron_pene   .push_back(pmtrack[i].trk[t].iron_pene);

	  pmanasum -> veto        .push_back(pmtrack[i].trk[t].vetowtracking);
	  pmanasum -> edge        .push_back(pmtrack[i].trk[t].edgewtracking);
	  pmanasum -> pdg         .push_back(pmtrack[i].trk[t].pdg);
	  pmanasum -> trkpe       .push_back(pmtrack[i].trk[t].trkpe);
	  pmanasum -> mucl        .push_back(pmtrack[i].trk[t].mucl);

	}
	//cout<<"In the End="<<(pmanasum -> GetIngridHitTrk(0,0))->pe<<endl;                                                                            
	//if(pmanasum->Ntrack>1) cout<<"In the End 222222="<<(pmanasum -> GetIngridHitTrk(0,1))->pe<<endl<<endl;
	//for(int t=0;t<pmtrack[i].trk.size();t++) cout<<"number of hits="<<pmanasum->NhitTs(t)<<endl;

	evt   -> AddPMAna(pmanasum);
      }
      
      
    }//cyc

    
    
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

void GetIngNonRecHits(IngridEventSummary* evt, int cyc){
  IngridHitSummary*    inghitsum;
  int imod,view,pln,ch;
  float pe,lope;
  int ninghit,pdg;
  memset(ingnonrechits_pe,  0,sizeof(ingnonrechits_pe));
  memset(ingnonrechits_lope,0,sizeof(ingnonrechits_lope));
  memset(ingnonrechits_pdg ,0,sizeof(ingnonrechits_pdg ));
  memset(ingnonrechits_id  ,0,sizeof(ingnonrechits_id ));
  

  // store info of all IngridHits
  for(int mod=0; mod<7; mod++){
    ninghit = evt -> NIngridModHits(mod, cyc);
    for(int ihit=0; ihit<ninghit; ihit++){
      inghitsum  = (IngridHitSummary*) (evt -> GetIngridModHit(ihit,mod,cyc));
      if(Is_Bad_Channel( inghitsum ))continue;
      view = inghitsum->view;
      ch   = inghitsum->ch;
      pln  = inghitsum->pln;
      pe   = inghitsum->pe;
      lope = inghitsum->lope;
      if((inghitsum -> NSimHits()) > 0)pdg = inghitsum -> GetIngridSimHit(0)->pdg;
      else pdg = 0;
      if(pln>=Cpln) continue; // pln>=11 are veto planes
      ingnonrechits_pe[mod][view][pln][ch]=pe;
      ingnonrechits_lope[mod][view][pln][ch]=lope;
      ingnonrechits_pdg [mod][view][pln][ch]=pdg;
      ingnonrechits_id [mod][view][pln][ch]=ihit;
    }
  }

  // now delete all reconstructed hits
  int hn=hingtrack.size();
  int vn=vingtrack.size();
  for(int i=0;i<hn;i++){
    int nhit=hingtrack[i].hit.size();
    for(int j=0;j<nhit;j++){
      imod=hingtrack[i].hit[j].mod;
      view=hingtrack[i].hit[j].view;
      pln=hingtrack[i].hit[j].pln;
      ch=hingtrack[i].hit[j].ch;
      ingnonrechits_pe[imod][view][pln][ch]=0;
      ingnonrechits_lope[imod][view][pln][ch]=0;
      ingnonrechits_pdg[imod][view][pln][ch]=0;
      ingnonrechits_id[imod][view][pln][ch]=0;
    }
  }

  for(int i=0;i<vn;i++){
    int nhit=vingtrack[i].hit.size();
    for(int j=0;j<nhit;j++){
      imod=vingtrack[i].hit[j].mod;
      view=vingtrack[i].hit[j].view;
      pln=vingtrack[i].hit[j].pln;
      ch=vingtrack[i].hit[j].ch;
      ingnonrechits_pe[imod][view][pln][ch]=0;
      ingnonrechits_lope[imod][view][pln][ch]=0;
      ingnonrechits_pdg[imod][view][pln][ch]=0;
      ingnonrechits_id[imod][view][pln][ch]=0;
    }
  }
}

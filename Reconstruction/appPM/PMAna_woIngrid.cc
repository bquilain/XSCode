
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
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
#include "PMreconRev.hxx"
#include "PMAna_woIngrid.hxx"

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
bool MC=false;

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

  TH1D * PEFV = new TH1D("PEFV","",200,0,100);
  TH1D * PEDist = new TH1D("PEDist","",200,0,100);
  TH1D * PEDist2 = new TH1D("PEDist2","",200,0,100);
  TH1D * PEDist3 = new TH1D("PEDist3","",200,0,100);
  TH1D * PE1h0v = new TH1D("PE1h0v","Only 1 horizontal track",200,0,100);
  TH1D * PE0h1v = new TH1D("PE0h1v","Only 1 vertical track",200,0,100);
  TH1D * PE1h1v = new TH1D("PE1h1v","1 horizontal& 1 vertical track",200,0,100);
  TH1D * PE2h1v = new TH1D("PE2h1v","2 horizontal tracks & 1 vertical track",200,0,100);
  TH1D * PE1h2v = new TH1D("PE1h2v","1 horizontal track & 2 vertical tracks",200,0,100);
  TH1D * PE2h2v = new TH1D("PE2h2v","2 horizontal & 2 vertical tracks",200,0,100);
  TH1D * PEMore = new TH1D("PEMore","More than 2 horizontal & 2 vertical tracks",200,0,100);

  TH1D * PE1h0v2 = new TH1D("PE1h0v2","Only 1 horizontal track",200,0,100);
  TH1D * PE0h1v2 = new TH1D("PE0h1v2","Only 1 vertical track",200,0,100);
  TH1D * PE1h1v2 = new TH1D("PE1h1v2","1 horizontal& 1 vertical track",200,0,100);
  TH1D * PE2h1v2 = new TH1D("PE2h1v2","2 horizontal tracks & 1 vertical track",200,0,100);
  TH1D * PE1h2v2 = new TH1D("PE1h2v2","1 horizontal track & 2 vertical tracks",200,0,100);
  TH1D * PE2h2v2 = new TH1D("PE2h2v2","2 horizontal & 2 vertical tracks",200,0,100);
  TH1D * PEMore2 = new TH1D("PEMore2","More than 2 horizontal & 2 vertical tracks",200,0,100);

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
  while ((c = getopt(argc, argv, "m:r:s:f:cado:i:")) != -1) {
    switch(c){
    case 'm':
      MC=true;
      break;
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
  char OutNameFV[256];
  sprintf(OutNameFV,"PEDist_%d_FVPureRead.root",sub_run_number);
  char OutName[256];
  sprintf(OutName,"PEDist_%d_FV_ReadByPMAna.root",sub_run_number);
  char OutName1h0v[256];
  sprintf(OutName1h0v,"PEDist_%d_FV_ReadByPMAna_1h0v.root",sub_run_number);
  char OutName0h1v[256];
  sprintf(OutName0h1v,"PEDist_%d_FV_ReadByPMAna_0h1v.root",sub_run_number);
  char OutName1h1v[256];
  sprintf(OutName1h1v,"PEDist_%d_FV_ReadByPMAna_1h1v.root",sub_run_number);
  char OutName2h1v[256];
  sprintf(OutName2h1v,"PEDist_%d_FV_ReadByPMAna_2h1v.root",sub_run_number);
  char OutName1h2v[256];
  sprintf(OutName1h2v,"PEDist_%d_FV_ReadByPMAna_1h2v.root",sub_run_number);
  char OutName2h2v[256];
  sprintf(OutName2h2v,"PEDist_%d_FV_ReadByPMAna_2h2v.root",sub_run_number);
  char OutName3h3v[256];
  sprintf(OutName3h3v,"PEDist_%d_FV_ReadByPMAna_3h3v.root",sub_run_number);

  char OutName1h0v2[256];
  sprintf(OutName1h0v2,"PEDist_%d_FV_AfterPMAna_1h0v.root",sub_run_number);
  char OutName0h1v2[256];
  sprintf(OutName0h1v2,"PEDist_%d_FV_AfterPMAna_0h1v.root",sub_run_number);
  char OutName1h1v2[256];
  sprintf(OutName1h1v2,"PEDist_%d_FV_AfterPMAna_1h1v.root",sub_run_number);
  char OutName2h1v2[256];
  sprintf(OutName2h1v2,"PEDist_%d_FV_AfterPMAna_2h1v.root",sub_run_number);
  char OutName1h2v2[256];
  sprintf(OutName1h2v2,"PEDist_%d_FV_AfterPMAna_1h2v.root",sub_run_number);
  char OutName2h2v2[256];
  sprintf(OutName2h2v2,"PEDist_%d_FV_AfterPMAna_2h2v.root",sub_run_number);
  char OutName3h3v2[256];
  sprintf(OutName3h3v2,"PEDist_%d_FV_AfterPMAna_3h3v.root",sub_run_number);

  char OutName2[256];
  sprintf(OutName2,"PEDist_%d_FV_AfterPMAna.root",sub_run_number);
  char OutName3[256];
  sprintf(OutName3,"PEDist_%d_FV_DirectAfterPMAna.root",sub_run_number);
 double Tot1Track=0;
      double Final1Track=0;
     
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
    sprintf(Output, "%s/ingrid_%08d_%04d_anas.root", 
	    dst_data, run_number, sub_run_number); 
  }
  if(cosmic){
    sprintf(Output, "%s/ingrid_%08d_%04d_anas.root", 
	    cosmic_data, run_number, sub_run_number); 
  }

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

  for(int ievt=Nini; ievt<nevt; ievt++){

    if(ievt%100==0)cout << "analyze event# " << ievt<<endl;
    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);

    
    
    for( int cyc=Scyc; cyc<Ncyc; cyc++ ){  //### Cycle Loop

      hingtrack.clear();
      vingtrack.clear();
      for( int mod=0; mod<7; mod++ ){  //### Module
	
	for( int view=0; view<2; view++ ){  //### View
	  int npmtrack = evt -> NPMModRecons(mod, cyc, view);//parcourt les traces 2D
	  //cout<<"npmtrack="<<npmtrack<<endl;
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

	      if((inghitsum -> NSimHits()) > 0)
     		hit.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	      else
		hit.pdg   = 0;

	      ingtrack.hit.push_back(hit);

	    }

	    if(view==0)
	      hingtrack.push_back(ingtrack);
	    else
	      vingtrack.push_back(ingtrack);



	  }
	}
      }//On a reécupéré les hits et on les a mis dans ingtrack.hit. ingtrack est une TrackIng => 1 trace 2D contenant l'info des hits. vingtrack et hingtrack sont des vecteurs de ces traces 2D => l'info est conservée dans ces classes TrackIng. ça c'est pour Ingrid!!! C'est grace a mes modifs que tout ça est conservé ;-)

      htrack.clear();
      vtrack.clear();
      //cout<<endl;
      double PE;double TrackNum=0;
      bool TrackFV=false;

      if(evt -> NPMModRecons(16, cyc, 0)==1 && evt -> NPMModRecons(16, cyc, 1)==0) TrackNum=0.5;
      else if(evt -> NPMModRecons(16, cyc, 0)==0 && evt -> NPMModRecons(16, cyc, 1)==1) TrackNum=1;
      else if(evt -> NPMModRecons(16, cyc, 0)==1 && evt -> NPMModRecons(16, cyc, 1)==1) TrackNum=1.5;
      else if(evt -> NPMModRecons(16, cyc, 0)==2 && evt -> NPMModRecons(16, cyc, 1)==1) TrackNum=2;
      else if(evt -> NPMModRecons(16, cyc, 0)==1 && evt -> NPMModRecons(16, cyc, 1)==2) TrackNum=2.5;
      else if(evt -> NPMModRecons(16, cyc, 0)==2 && evt -> NPMModRecons(16, cyc, 1)==2) TrackNum=3;
      else TrackNum=3.5;

      if(TrackNum==1.5 && (!((evt -> GetPMModRecon(0, 16, cyc, 0))->vetowtracking) && !((evt -> GetPMModRecon(0, 16, cyc, 1))->vetowtracking)) && (!((evt -> GetPMModRecon(0, 16, cyc, 0))->edgewtracking) && !((evt -> GetPMModRecon(0, 16, cyc, 1))->edgewtracking))){
	Tot1Track++;
	TrackFV=true;
      }
      //      if((evt -> NPMModRecons(16, cyc, 0) == evt -> NPMModRecons(16, cyc, 1)) && (!((evt -> GetPMModRecon(0, 16, cyc, 0))->vetowtracking) && !((evt -> GetPMModRecon(0, 16, cyc, 1))->vetowtracking)) && (!((evt ->GetPMModRecon(0, 16, cyc, 0))->edgewtracking) && !((evt -> GetPMModRecon(0, 16, cyc, 1))->edgewtracking))){

      for( int view=0; view<2; view++ ){  //### View
	int npmtrack = evt -> NPMModRecons(16, cyc, view);
	//if(npmtrack>0) cout<<"2D PM Tracks="<<npmtrack<<endl;
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

	  //cout<<"Track is in View="<<view<<", NHits="<<recon->Nhits()<<", Vetow="<<recon -> vetowtracking<<", Edge="<<recon -> edgewtracking<<", Angle="<<recon ->angle<<endl;
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
	   
	    if((inghitsum->pe+inghitsum->lope)/2<40) PE=inghitsum->pe;
	    else PE=inghitsum->lope;
	    if(!(recon -> vetowtracking) && !(recon -> edgewtracking) && inghitsum->pln!=0 && inghitsum->ch>7 && inghitsum->ch<24){
	      //cout<<"Hits high pe="<<inghitsum->pe<<", plane="<<inghitsum->pln<<", channel="<<inghitsum->ch<<endl;
	   
	      if(TrackNum==0.5) PE1h0v->Fill(PE);
	      else if(TrackNum==1) PE0h1v->Fill(PE);
	      else if(TrackNum==1.5) PE1h1v->Fill(PE);
	      else if(TrackNum==2) PE2h1v->Fill(PE);  
              else if(TrackNum==2.5) PE1h2v->Fill(PE);
              else if(TrackNum==3) PE2h2v->Fill(PE);
	      else if(TrackNum==3.5) PEMore->Fill(PE);

	      if(TrackFV) PEFV->Fill(PE);
	    }
	    if((inghitsum -> NSimHits()) > 0)
	      hit.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	    else
	      hit.pdg   = 0;
	    
	    track.hit.push_back(hit);
	    
	  }

	  //cout<<"2D Track Size="<<track.hit.size()<<endl;
	  if(view==0)
	    htrack.push_back(track);
	  else
	    vtrack.push_back(track);
	}
      }
      //de meme pour le PM, normalement htrack et vtrack sont des vecteurs contenant l'indo sur les traces elles meme contenant l'info sur les hits => check. Oui track est une TrackPM => oui on a bien un vecteur de hits 
      
      if(!fPMAna())continue;//applique la reconstruction. C'est la ou on peut perdre notre information sur les hits. Il faut donc aller regarder cette fonction
      //
      //cout<<pmtrack.size()<<endl;     
      for(int i=0;i<pmtrack.size();i++){

	if(pmtrack[i].Ntrack == 0)continue;
	pmanasum -> Clear();
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


	pmanasum -> x           .clear();
	pmanasum -> y           .clear();
	pmanasum -> z           .clear();
	pmanasum -> zx          .clear();
	pmanasum -> zx          .clear();
	
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
	

	//cout<<"Amount of 3D tracks:"<<endl;
	//cout<<pmanasum->Ntrack<<"   size="<<pmtrack[i].trk.size()<<endl;

	for(int t=0;t<pmtrack[i].trk.size();t++){
	  //cout<<"Track Number="<<t<<endl;
	  //cout<<"Size of Track in Ana="<<pmtrack[i].trk[t].hit.size()<<endl;
	  //cout<<"AMOUNT OF HITS IN EVT="<<evt->NIngridHits()<<" and in the other="<<evt->NIngridModHits(16,pmanasum -> hitcyc)<<endl;
	  //cout<<"New 3D Track, vetow="<<pmtrack[i].trk[t].vetowtracking<<", NHits="<<pmtrack[i].trk[t].hit.size()<<endl;
	 
	  if(TrackNum==1.5 && !(pmtrack[i].trk[t].vetowtracking) /*&& !(pmtrack[i].trk[t].edgewtracking)*/) Final1Track++;
	     
	  for(int ihit=0;ihit<pmtrack[i].trk[t].hit.size();ihit++){

	    //	    cout<<"Hit highpe="<<pmtrack[i].trk[t].hit[ihit].pe<<", Plane="<<pmtrack[i].trk[t].hit[ihit].pln<<" ,Channel="<<pmtrack[i].trk[t].hit[ihit].ch<<endl;
	    /*    if(!(pmtrack[i].trk[t].vetowtracking) && !(pmtrack[i].trk[t].edgewtracking) && pmtrack[i].trk[t].hit[ihit].pln!=0 && pmtrack[i].trk[t].hit[ihit].ch>7 && pmtrack[i].trk[t].hit[ihit].ch<24 && pmtrack[i].trk[t].hit[ihit].mod==16){
	      if((pmtrack[i].trk[t].hit[ihit].pe+pmtrack[i].trk[t].hit[ihit].lope)/2<40) PEDist3->Fill(pmtrack[i].trk[t].hit[ihit].pe);
	      else PEDist3->Fill(pmtrack[i].trk[t].hit[ihit].lope);
	      //cout<<"Hit highpe="<<pmtrack[i].trk[t].hit[ihit].pe<<", Plane="<<pmtrack[i].trk[t].hit[ihit].pln<<" ,Channel="<<pmtrack[i].trk[t].hit[ihit].ch<<endl;
	      if((pmtrack[i].trk[t].hit[ihit].pe+pmtrack[i].trk[t].hit[ihit].lope)/2<40) PE=pmtrack[i].trk[t].hit[ihit].pe;
	      else PE=pmtrack[i].trk[t].hit[ihit].lope;
	      if(TrackNum==0.5) PE1h0v2->Fill(PE);
	      else if(TrackNum==1) PE0h1v2->Fill(PE);
	      else if(TrackNum==1.5) PE1h1v2->Fill(PE);
	      else if(TrackNum==2) PE2h1v2->Fill(PE);
	      else if(TrackNum==2.5) PE1h2v2->Fill(PE);
	      else if(TrackNum==3) PE2h2v2->Fill(PE);
	      else if(TrackNum==3.5) PEMore2->Fill(PE);
	   
	    }
*/
	    bool HitFound=false;
	    for(int imod=0;imod<Nmod;imod++){
	      int npmtrack = evt -> NPMModRecons(imod, cyc, pmtrack[i].trk[t].hit[ihit].view);
	      for(int l=0; l<npmtrack; l++){
		if(HitFound) break;
		recon   = (PMReconSummary*) (evt -> GetPMModRecon(l, imod, cyc, pmtrack[i].trk[t].hit[ihit].view) );
		for(int NHIT=0;NHIT<(recon->Nhits());NHIT++){
		  if(HitFound) break;
		  pmhit = recon -> GetIngridHit(NHIT);
		  //cout<<pmhit->pe<<", mod="<<pmhit->mod<<", "<<pmtrack[i].trk[t].hit[ihit].mod<<endl;
		  if(pmhit->mod==pmtrack[i].trk[t].hit[ihit].mod && pmhit->view==pmtrack[i].trk[t].hit[ihit].view && pmhit->pln==pmtrack[i].trk[t].hit[ihit].pln && pmhit->ch==pmtrack[i].trk[t].hit[ihit].ch && pmhit->pe==pmtrack[i].trk[t].hit[ihit].pe && pmhit->lope==pmtrack[i].trk[t].hit[ihit].lope && pmhit->xy==pmtrack[i].trk[t].hit[ihit].xy && pmhit->z==pmtrack[i].trk[t].hit[ihit].z && pmhit->time==pmtrack[i].trk[t].hit[ihit].time){
		  pmanasum -> AddIngridHitTrk(pmhit,t);
		  //cout<<"Hit added with mod="<<pmhit->mod<<endl;
		  if(!(pmtrack[i].trk[t].vetowtracking) && !(pmtrack[i].trk[t].edgewtracking) && pmhit->pln!=0 && pmhit->ch>7 && pmhit->ch<24){
		    if((pmhit->pe+pmhit->lope)/2<40) PEDist2->Fill(pmhit->pe);
		    else PEDist2->Fill(pmhit->lope);
		    }
		    //cout<<pmhit->pe<<endl;
		  HitFound=true;
		  }
		}
              } 
            }
	    //cout<<pmhit->pe<<endl;
	    //pmanasum -> AddIngridHitTrk(pmhit,t);
	    //cout<<"Track Number="<<t<<"  , Energy="<<(pmanasum -> GetIngridHitTrk(ihit,t))->pe<<endl;
	    //pmanasum -> hit[pmtrack[i].trk[t]].pushback(pmtrack[i].trk[t].hit[ihit]);
	  }
	  if(pmanasum->Ntrack>1 && pmanasum->NhitTs(t)>99){
	    //cout<<"***********************************************************************************************************************************"<<endl;
	    //cout<<pmanasum->NhitTs(t)<<"  "<<pmanasum->nhitTs[t]<<endl;
	    //cout<<"Just After Hits were Filled="<<(pmanasum -> GetIngridHitTrk(0,0))->pe<<"   "<<(pmanasum -> GetIngridHitTrk(t,0))->pe<<endl;
	    //if(t==1) cout<<"Just After 2="<<(pmanasum -> GetIngridHitTrk(0,1))->pe<<endl;
	  }
	  pmanasum -> x           .push_back(pmtrack[i].trk[t].x);
	  pmanasum -> y           .push_back(pmtrack[i].trk[t].y);
	  pmanasum -> z           .push_back(pmtrack[i].trk[t].z);
	  pmanasum -> zx          .push_back(pmtrack[i].trk[t].zx);
	  pmanasum -> zx          .push_back(pmtrack[i].trk[t].zx);

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
	  //cout<<"INGRID Trk="<<pmanasum -> ing_trk.back()<<endl;//"  ";
	  pmanasum -> pm_stop     .push_back(pmtrack[i].trk[t].pm_stop);
	  pmanasum -> ing_stop    .push_back(pmtrack[i].trk[t].ing_stop);
	  pmanasum -> sci_range   .push_back(pmtrack[i].trk[t].sci_range);
	  pmanasum -> iron_range  .push_back(pmtrack[i].trk[t].iron_range);
	  pmanasum -> iron_pene   .push_back(pmtrack[i].trk[t].iron_pene);

	  pmanasum -> veto        .push_back(pmtrack[i].trk[t].vetowtracking);
	  pmanasum -> edge        .push_back(pmtrack[i].trk[t].edgewtracking);
	  pmanasum -> pdg         .push_back(pmtrack[i].trk[t].pdg);
	  pmanasum -> trkpe       .push_back(pmtrack[i].trk[t].trkpe);
       
	}

	//cout<<"In the End="<<(pmanasum -> GetIngridHitTrk(0,0))->pe<<endl;
	//if(pmanasum->Ntrack>1) cout<<"In the End 222222="<<(pmanasum -> GetIngridHitTrk(0,1))->pe<<endl<<endl;
	evt   -> AddPMAna(pmanasum);
	//cout<<endl;
      }
      
      
    }//cyc

    
    
    wsummary = evt;
    wtree -> Fill();
  }//Event Loop
  PEFV->SaveAs(OutNameFV);
  PEDist->SaveAs(OutName);
  PEDist2->SaveAs(OutName2);
  PEDist3->SaveAs(OutName3);
  
  PE1h0v->SaveAs(OutName1h0v);
  PE0h1v->SaveAs(OutName0h1v);
  PE1h1v->SaveAs(OutName1h1v);
  PE2h1v->SaveAs(OutName2h1v);
  PE1h2v->SaveAs(OutName1h2v);
  PE2h2v->SaveAs(OutName2h2v);
  PEMore->SaveAs(OutName3h3v);


  PE1h0v2->SaveAs(OutName1h0v2);
  PE0h1v2->SaveAs(OutName0h1v2);
  PE1h1v2->SaveAs(OutName1h1v2);
  PE2h1v2->SaveAs(OutName2h1v2);
  PE1h2v2->SaveAs(OutName1h2v2);
  PE2h2v2->SaveAs(OutName2h2v2);
  PEMore2->SaveAs(OutName3h3v2);
  //cout<<"Final Ratio="<<Final1Track/Tot1Track<<", Tot="<<Tot1Track<<endl;
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


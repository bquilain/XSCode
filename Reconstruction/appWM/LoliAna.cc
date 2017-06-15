#define USEINGTRK
#define USEINGHIT
//#define USEINGPID
//#define USEBACKTRK
//#define USEPARTRK
//#define USESHORTTRK
 
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
//#include "PMrecon.hxx"
#include "INGRID_Dimension.hxx"
#include "Lolirecon.hxx"
#include "LoliAna.hxx" 

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

#include "Lolidisp.h"

// already defined in IngridConstants
//const static int   GapbwBunch     = 581 ;    //[nsec]
//const static int   TDCOffset      = 300;     //[nsec]
//const static int   fExpBeamTime   = 200;     //[nsec]
//const static int   beamontime     = 100;     //ontime

//#define  STUDY_BADCHANNEL
vector<int>  bad_mod;
vector<int>  bad_pln;
vector<int>  bad_view;
vector<int>  bad_ch;
void Add_Bad_Channel();
void fMode(int mode);
bool Is_Bad_Channel(IngridHitSummary* thit);
void GetNonRecHits(IngridEventSummary* evt, int cyc);//newly added

//const static char *dst_data          = "/home/kikawa/scraid1/data/pm_ingrid";
//const static char *cosmic_data       = "/home/kikawa/scraid1/data/pm_ingrid";

//static const char *data_file_folder  = "/home/kikawa/scraid1/data/pm_ingrid";
//static const char *calib_file_folder = "/home/kikawa/scraid1/data/pm_ingrid";
//static const char *dailycheck_folder = "/home/kikawa/scraid1/data/pm_ingrid";



//#include "IngridConstants.h"
long fcTime;

const double fiberIndex=1.91; // correction; see talk of 2017/01/13
const double c_vac=30; //speed of light in vacuum, cm/ns

int main(int argc,char *argv[]){

//#ifdef STUDY_BADCHANNEL
//  cout << "It is BAD channel study mode. Is it O.K.?" << endl;
//  //cin.get();
//  Add_Bad_Channel();
//#endif

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
  bool useINGRID_PID=true;

  // to evaluate systematics
  int VertexingPlane  =pln_th;//planes
  double VertexingChannel =ch_th;//mm
  int TrackMatching=diff_th; // planes
  double AngleCut=ang_th;//degrees
  double TransverseCut=pos_th;//mm
  int nINGRIDPlanes=2; // number of planes fIngridHitPMJoint is looking at
  bool Error=false;
  int ErrorType=0;
  float ErrorValue=0.;

  bool requireIngridTrack=true; //ML 2016/11/24 for new option -N

  while ((c = getopt(argc, argv, "r:s:f:cado:i:e:v:NI")) != -1) {
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
    case 'I':
      useINGRID_PID=false;
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

  // ML 2016/11/24
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
  if(!renamef && cosmic){
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
  if(!rename && cosmic){
    sprintf(Output, "%s/ingrid_%08d_%04d_anas.root", 
	    cosmic_data, run_number, sub_run_number); 
  }

  TFile*            wfile = new TFile(Output, "recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* wsummary = new IngridEventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary", 
				  &wsummary,  64000,  99);
  IngridHitSummary*    inghitsum, *inghitsum2;
  IngridSimVertexSummary*   simver;
  PMReconSummary* recon;
  //PMReconSummary* recon = new PMReconSummary();
  PMAnaSummary* pmanasum = new PMAnaSummary();
  PMAnaSummary* pmanasum2 = new PMAnaSummary();
  Hit                        hit;
  TrackPM                      track;
  TrackIng                     ingtrack;

  Hits                        hits;
 
  LoadMuCL("$(INSTALLREPOSITORY)/Reconstruction/inc/sandmuon_WM_distributions_47k_cut4.5pe_normalAngleSlices.root",false);
  LoadMuCL_I("$(INSTALLREPOSITORY)/Reconstruction/inc/sandmuon_distributions_mc_cut4.5_MIP31.5_wINGRID.root");
  int biasedMuCL=0;

  cout<<(useINGRID_PID?"USE INGRID HITS FOR PID":"ONLY WM HITS FOR PID")<<endl;

  Initialize_INGRID_Dimension();

  //  nevt=Nini+1;

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
	  int npmtrack = evt -> NPMModRecons(mod, cyc, view);
	  //	  if(mod==3) cout<<" in mod 3 recons="<<npmtrack<<endl;
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
	      hits.clear(); //ML tmp
	      hits.mod   = inghitsum->mod;
	      hits.view  = inghitsum->view;
	      hits.pln   = inghitsum->pln;
	      hits.ch    = inghitsum->ch;
	      hits.pe    = inghitsum->pe+inghitsum->pe_cross; //2017-01-17 to read calib info and 2017/03/21 to read crosstalk
	      hits.lope  = inghitsum->lope;
	      hits.isohit= inghitsum->isohit;

	      hits.recon_id=i;
	      hits.hit_id=NHIT;

	      //if((inghitsum -> NSimHits()) > 0 && inghitsum -> cyc >= 0)
	      if((inghitsum -> NSimHits()) > 0)
     		hits.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	      else
		hits.pdg   = 0;

	      ingtrack.hit.push_back(hits);

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
	int npmtrack = evt -> NPMModRecons(15, cyc, view);
	  
	for(int i=0; i<npmtrack; i++){
	  recon   = (PMReconSummary*) (evt -> GetPMModRecon(i, 15, cyc, view) );
	  //recon   = (PMReconSummary*) (evt -> GetPMModRecon(i, 16, cyc, view) );
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
	    hits.clear(); //ML tmp
	    hits.mod   = inghitsum->mod;
	    hits.view  = inghitsum->view;
	    hits.pln   = inghitsum->pln;
	    hits.ch    = inghitsum->ch;
	    // since I did not rerun calibration with the new switching point I use pe instead of pecorr. 
	    hits.pe    = inghitsum->pe+inghitsum->pe_cross;
	    hits.lope  = inghitsum->lope;
	    // now I compute the switch
	    if(evt->NIngridSimVertexes()==0) //ie data
	      hits.pe=(0.5*(hits.pe+hits.lope)<43? hits.pe:hits.lope);

	    hits.isohit= inghitsum->isohit;

	    hits.recon_id=i;
	    hits.hit_id=NHIT;

	    //axis=0
	    INGRID_Dimension *fdim_temp = new INGRID_Dimension();
	    int reconpln,reconch;
            fdim_temp -> get_reconplnch_loli( hits.mod, hits.view, hits.pln, hits.ch, 0, &reconpln, &reconch ); //mod view pln axis
	    hits.pln= reconpln;
	    hits.ch = reconch;
	    delete fdim_temp;
	    
	    //if((inghitsum -> NSimHits()) > 0 && inghitsum -> cyc >= 0)
	    if((inghitsum -> NSimHits()) > 0)
	      hits.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	    else
	      hits.pdg   = 0;
	    
	    track.hit.push_back(hits);
	    
	  }

	  if(view==0)
	    htrack.push_back(track);
	  else
	    vtrack.push_back(track);
	}
      }
      
      GetNonRecHits(evt,cyc);//newly added
      
      // ML 2016/11/24 requireIngridTrack for new option -N
      if(!fPMAna(15,TrackMatching,VertexingPlane,VertexingChannel,AngleCut,TransverseCut,nINGRIDPlanes,requireIngridTrack))continue;
     
      // isoHitCut is the minimal number of isolated it I require
      int isoHitCut=2;
      //      cout<<"\nMatching done!"<<endl;
      for(int i=0;i<pmtrack.size();i++){
	//cout<<"reco="<<i<<" nb of tracks="<<pmtrack[i].trk.size()<<endl;	
	pmanasum->Clear();

	reCalcIsohit(pmtrack[i],isoHitCut);
	if(pmtrack[i].Ntrack == 0)continue;

	//Vertex activity
	int ninghit = evt -> NIngridModHits(15, cyc);
	int xpln  = pmtrack[i].trk[0].startxpln;
	int ypln  = pmtrack[i].trk[0].startypln;
	float xch = pmtrack[i].trk[0].startxch;
	float ych = pmtrack[i].trk[0].startych;
	//float vact=0;
	float vact[100];
	for(int j=0;j<100;j++){
		vact[j]=0;
	}
	for(int ihit=0; ihit<ninghit; ihit++){
	  inghitsum  = (IngridHitSummary*) (evt -> GetIngridModHit(ihit,15,cyc));
	  if(Is_Bad_Channel( inghitsum ))continue;
	  //vact += veract(inghitsum,xpln,ypln,xch,ych);
	  vact[0] += veract(inghitsum,xpln,ypln,xch,ych,3,80 ,0);
	  vact[1] += veract(inghitsum,xpln,ypln,xch,ych,6,160,0);
	  vact[2] += veract(inghitsum,xpln,ypln,xch,ych,3,80 ,1);
	  vact[3] += veract(inghitsum,xpln,ypln,xch,ych,6,160,1);
	}

	//use non rec hit
	vact[4] += veract(xpln,ypln,xch,ych,3,80 ,0);
	vact[5] += veract(xpln,ypln,xch,ych,6,160,0);
	vact[6] += veract(xpln,ypln,xch,ych,3,80 ,1);
	vact[7] += veract(xpln,ypln,xch,ych,6,160,1);

	pmanasum -> Ntrack = min(pmtrack[i].Ntrack,INGRIDRECON_MAXTRACKS);
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
	//pmanasum -> clstimecorr = vact;//temporary clstimecorr is used to fill the vertex activity information
	for(int j=0;j<100;j++){
		pmanasum -> vact[j]=0;
		pmanasum -> vact[j]=vact[j];
	}

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
	
	pmanasum->nhits=0;
	for(int itrk=0;itrk<INGRIDRECON_MAXTRACKS;itrk++) pmanasum->nhitTs[itrk]=0;

	for(int t=0;t<min((int)pmtrack[i].trk.size(),INGRIDRECON_MAXTRACKS);t++){

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

	  if(!calcMuCL(pmtrack[i].trk[t],isoHitCut,useINGRID_PID))  {
	    biasedMuCL++;
	    pmanasum->mucl .push_back( -1);
	  }
	  else 
	    pmanasum -> mucl        .push_back(pmtrack[i].trk[t].mucl);

	  // if(pmtrack[i].trk[t].mucl!=pmtrack[i].trk[t].mucl) cout<<"******************"<<endl;


	  //	  if(pmtrack[i].trk[t].pdg==13 && pmtrack[i].trk[t].ing_trk && pmtrack[i].trk[t].mucl<0.05) cout<<ievt<<" "<<i<<" "<<t<<" "<<pmtrack[i].trk[t].mucl<<endl;
	  
	  vector<Hits>::iterator iter;
	  for(iter=pmtrack[i].trk[t].hit.begin();iter<pmtrack[i].trk[t].hit.end();iter++){
	    Hits hit=(*iter);
	    //  cout<<hit.recon_id<<" "<<hit.hit_id<<" "<<hit.view<<" "<<hit.mod<<" "<<hit.pln<<endl;
	    inghitsum->Clear();

	    if(hit.recon_id!=-1)
	      inghitsum=((PMReconSummary*)evt->GetPMModRecon(hit.recon_id,hit.mod,cyc,hit.view))->GetIngridHit(hit.hit_id);
	    else
	      inghitsum=evt->GetIngridModHit(hit.hit_id,hit.mod,cyc);
	    
	    inghitsum->isohit=hit.isohit;
	    pmanasum->AddIngridHitTrk(inghitsum,t);
	  }//hit
	}
	evt   -> AddPMAna(pmanasum);

      }
     
      //here I correct for propagation in the fiber - ML 2016/12/16
      int nHits=evt->NIngridModHits(15,cyc);
      int nAna=evt->NPMAnas();
      //      if(nAna>0) cout << nAna <<endl;
      for(int ihit=0;ihit<nHits;ihit++){
	inghitsum=evt->GetIngridModHit(ihit, 15,cyc);
	bool found=false;
	for(int iana=0;iana<nAna;iana++){	
	  pmanasum2=evt->GetPMAna(iana);
	  int nTracks=pmanasum2->Ntrack;
	  int itrk;
	  for(itrk=0;itrk<nTracks;itrk++){
	    int nHits2=pmanasum2->NhitTs(itrk);
	    for(int ihit2=0;ihit2<nHits2;ihit2++){
	      inghitsum2=pmanasum2->GetIngridHitTrk(ihit2,itrk);
	      if(inghitsum->mod==inghitsum2->mod && inghitsum->view==inghitsum2->view && inghitsum->pln==inghitsum2->pln && inghitsum->ch==inghitsum2->ch){
		found=true;
		break;
	      }
	    }//ihit2
	    if(found) break;
	  }

	  if(found){
	    //now iana and itrk are set to the first reco and track inghitsum belongs to
	    double tx=pmanasum2->thetax[0];
	    double ty=pmanasum2->thetay[0];
	    double startx=(pmanasum2->startxch[0]-100)/10;//cm
	    double starty=(pmanasum2->startych[0]-100)/10;//cm
	    
	    double startz_x, startz_y, posixy; // I'm only interested in startz_x/y
	    INGRID_Dimension *a_temp = new INGRID_Dimension();
	    a_temp -> get_posi_lolirecon( 15, 0,pmanasum2->startxpln[0],0,0, &posixy, &startz_x);
	    a_temp -> get_posi_lolirecon( 15, 1,pmanasum2->startypln[0],0,0, &posixy, &startz_y);
	    delete a_temp;
	    
	    double fiberLength;
	    if(inghitsum->view==0) //missing coordinate is x
	      fiberLength=tan(ty*TMath::Pi()/180)*(inghitsum->z-startz_y)+starty;
	    else  // missing coordinate is y
	      fiberLength=100-(tan(tx*TMath::Pi()/180)*(inghitsum->z-startz_x)+startx);
	    
	    if(fiberLength>100) fiberLength=100;
	    if(fiberLength<0) fiberLength=0;
	    
	    inghitsum->timecorr -= fiberLength*fiberIndex/c_vac; // in ns
	    //  cout<<"found! "<<ihit<<" "<<iana<<","<<pmanasum2->ing_trk[itrk]<<" "<<nTracks<<","<<itrk<<" "<<nHits<<" "<<fiberLength<<" "<<fiberLength*fiberIndex/c_vac<<endl;
	    break;
	  }
	}//iana     
	//	if(!found) cout<<"not found..."<<inghitsum->mod<<" "<<inghitsum->pecorr<<" "<<endl;
      }//ihit

      //now I correct hitcls by the mean delay (fiberLength=50)
      for(int iana=0;iana<nAna;iana++){	
	pmanasum2=evt->GetPMAna(iana);
	pmanasum2->clstime -= 50*fiberIndex/c_vac; // in ns
      }


	  
	  if(disp && cyc==4){
	    for( int mod=0; mod<Nmod; mod++ ){   //### Module Loop
	     
	     int ninghit = evt -> NIngridModHits(mod, cyc);
	   
	     //if(ninghit<6)continue;
	   
	     for(int i=0; i<ninghit; i++){
	       inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(i, mod, cyc) );
	       if(inghitsum->cyc==-2)continue; //This is killed hit 
	       if(Is_Bad_Channel( inghitsum ))continue;
	       if(mod==15){
	     	if(inghitsum -> pe < 1.5) continue;
	       }
	       else{
	     	if(inghitsum -> pe < 2.5) continue;
	       }
	   
	       inghitsum -> addbasicrecon = false;
	   
	       hit.mod   = mod;
	       hit.id    = i;
	       hit.pe    = inghitsum -> pe;
	       hit.lope  = inghitsum -> lope;
	       hit.time  = inghitsum -> time;
	       hit.view  = inghitsum -> view;
	       hit.pln   = inghitsum -> pln;
	       hit.ch    = inghitsum -> ch;
	   
	       /*
	       if((inghitsum -> NSimHits()) > 0)
	         hit.pdg   = inghitsum -> GetIngridSimHit(0)->pdg;
	       else
	         hit.pdg   = 0;
	       */
	   
	       //hit.posxy = xyposi(mod, hit.view, hit.pln, hit.ch);
	       //hit.posz  = zposi (mod, hit.pln,  hit.ch,  hit.ch);
	   
	       allhit_for_disp.push_back(hit);
	     }
	    }

	    for(int i=0;i<pmtrack.size();i++){
	      cout<<(i+1);
	      if((i+1)==1)     cout<<"st track:";
	      else if((i+1)==2)cout<<"nd track:";
	      else if((i+1)==3)cout<<"rd track:";
	      else             cout<<"th track:";
	      cout<<" Veto:"<<pmtrack[i].vetowtracking
		  <<" Edge:"<<pmtrack[i].edgewtracking
		  <<" Endpln:"<<pmtrack[i].trk[0].ing_endpln
		  <<" pdg:"<<pmtrack[i].trk[0].pdg
		  <<" INGSTOP:"<<pmtrack[i].trk[0].ing_stop<<endl;
	    }

	    cout << "displayed event# " << ievt<<endl;
	    cout << "displayed cycle  " << cyc<<endl;

	    if(evt -> NIngridSimVertexes() > 0){
	      cout << "Nsimver " << evt -> NIngridSimVertexes() << endl;
	      simver   = (IngridSimVertexSummary*) (evt -> GetSimVertex(0));
	      simver_for_disp  = (IngridSimVertexSummary*) (evt -> GetSimVertex(0));
	      for(int j=0;j<evt->NIngridSimParticles();j++){
		cout << "secondary particle " << j << " pdg " << evt -> GetSimParticle(j)->pdg << 
			" momentum " << evt -> GetSimParticle(j)->momentum[0] << " " << evt -> GetSimParticle(j)->momentum[1] << " " << evt -> GetSimParticle(j)->momentum[2] << endl;
	        if(evt -> GetSimParticle(j)->pdg==13){
	          simpar_for_disp  = (IngridSimParticleSummary*) (evt -> GetSimParticle(j));
	        }
	      }
	      fMode(simver->inttype);
	    }


	    EvtDisp_Analoli();

	    allhit_for_disp.clear();
	  }

 
      
    }//cyc

    
    
    wsummary = evt;
    wtree -> Fill();
  }//Event Loop

//######## Write and Close ####################
  wfile->cd();
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

  cout<<wfile->GetName()<<" completed!"<<endl;
  cout<<"biased MuCL="<<biasedMuCL<<" (Nisohits<3)"<<endl;
}

//update from Koga, 2017/03/21, copy from Loli_addcrosstalk_slit.cc
const int bad_num = 6+13+22+2;
int       bad_id[bad_num][4] =
  {// mod, view,  pln,  ch, 
    //INGRID 6
    {    2,    0,   13,    2},
    {    3,    0,   10,   14},
    {    4,    1,   10,    9},
    {    4,    1,    0,   15},
    {    4,    1,    7,   21},
    {    4,    1,    4,   14},
    //WM MPPC 13
    {   15,    1,    2,   75},
    {   15,    1,    4,   21},
    {   15,    1,    4,   20},
    {   15,    1,    4,   49},
    {   15,    1,    4,   71},
    {   15,    1,    5,   19},
    {   15,    1,    5,   20},
    {   15,    1,    5,   71},
    {   15,    1,    5,   70},
    {   15,    1,    4,   48},
    {   15,    1,    5,   48},
    {   15,    1,    0,   74},
    {   15,    1,    0,    0},
    //WM scinti 22
    {   15,    0,    0,   20},
    {   15,    0,    0,   21},
    {   15,    0,    0,   26},
    {   15,    0,    0,   27},
    {   15,    0,    0,   28},
    {   15,    0,    0,   30},
    {   15,    0,    0,   54},
    {   15,    0,    0,   78},
    {   15,    0,    4,   40},
    {   15,    0,    4,   60},
    {   15,    0,    5,   56},
    {   15,    0,    7,   61},
    {   15,    0,    7,   79},
    {   15,    0,    0,   39},
    {   15,    1,    2,   40},
    {   15,    1,    3,   67},
    {   15,    1,    7,   24},
    {   15,    1,    7,   32},
    {   15,    1,    7,   63},
    {   15,    1,    7,   64},
    {   15,    1,    7,   71},
    {   15,    1,    7,   75},
    //PM 2
    {   16,    0,   14,    0},
    {   16,    1,    9,   21}
  };
/*
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
*/
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


//newly added
void GetNonRecHits(IngridEventSummary* evt, int cyc){
  IngridHitSummary*    inghitsum;
  int imod,view,pln,ch;
  float pe,lope;
  int ninghit,pdg;
  memset(nonrechits,        0,sizeof(nonrechits));
  memset(ingnonrechits,     0,sizeof(ingnonrechits));
  memset(ingnonrechits_lope,0,sizeof(ingnonrechits_lope));
  memset(ingnonrechits_pdg ,0,sizeof(ingnonrechits_pdg ));
  memset(ingnonrechits_id  ,0,sizeof(ingnonrechits_id ));
  
 
  ninghit = evt -> NIngridModHits(15, cyc);
  INGRID_Dimension *fdim_temp = new INGRID_Dimension();
  for(int ihit=0; ihit<ninghit; ihit++){
    inghitsum  = (IngridHitSummary*) (evt -> GetIngridModHit(ihit,15,cyc));
    if(Is_Bad_Channel( inghitsum ))continue;
    view = inghitsum->view;
    ch   = inghitsum->ch;
    pln  = inghitsum->pln;
    //pe   = PE(inghitsum->pe,inghitsum->lope,15,view,pln,ch);
    pe   = inghitsum->pe;

    //axis=0
    int reconpln,reconch;
    fdim_temp -> get_reconplnch_loli( 15, view, pln, ch, 0, &reconpln, &reconch ); //mod view pln axis
    pln= reconpln;
    ch = reconch;

    nonrechits[view][pln][ch]=pe;
  }
  delete fdim_temp;

  for(int mod=0; mod<7; mod++){
    ninghit = evt -> NIngridModHits(mod, cyc);
    for(int ihit=0; ihit<ninghit; ihit++){
      inghitsum  = (IngridHitSummary*) (evt -> GetIngridModHit(ihit,mod,cyc));
      if(Is_Bad_Channel( inghitsum ))continue;
      view = inghitsum->view;
      ch   = inghitsum->ch;
      pln  = inghitsum->pln;
      //pe   = PE(inghitsum->pe,inghitsum->lope,mod,view,pln,ch);
      pe   = inghitsum->pe;
      lope = inghitsum->lope;
      if((inghitsum -> NSimHits()) > 0)pdg = inghitsum -> GetIngridSimHit(0)->pdg;
      else pdg = 0;
      if(pln>=11) continue; // veto planes of INGRID modules
      ingnonrechits[mod][view][pln][ch]=pe;
      ingnonrechits_lope[mod][view][pln][ch]=lope;
      ingnonrechits_pdg [mod][view][pln][ch]=pdg;
      ingnonrechits_id [mod][view][pln][ch]=ihit;
    }
  }

  int hn,vn;
  hn=htrack.size();
  vn=vtrack.size();

  for(int i=0;i<hn;i++){
    int nhit=htrack[i].hit.size();
    for(int j=0;j<nhit;j++){
      view=htrack[i].hit[j].view;
      pln=htrack[i].hit[j].pln;
      ch=htrack[i].hit[j].ch;
      nonrechits[view][pln][ch]=0;
    }
  }

  for(int i=0;i<vn;i++){
    int nhit=vtrack[i].hit.size();
    for(int j=0;j<nhit;j++){
      view=vtrack[i].hit[j].view;
      pln=vtrack[i].hit[j].pln;
      ch=vtrack[i].hit[j].ch;
      nonrechits[view][pln][ch]=0;
    }
  }

  
  hn=hingtrack.size();
  vn=vingtrack.size();
  for(int i=0;i<hn;i++){
    int nhit=hingtrack[i].hit.size();
    for(int j=0;j<nhit;j++){
      imod=hingtrack[i].hit[j].mod;
      view=hingtrack[i].hit[j].view;
      pln=hingtrack[i].hit[j].pln;
      ch=hingtrack[i].hit[j].ch;
      ingnonrechits[imod][view][pln][ch]=0;
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
      ingnonrechits[imod][view][pln][ch]=0;
      ingnonrechits_lope[imod][view][pln][ch]=0;
      ingnonrechits_pdg[imod][view][pln][ch]=0;
      ingnonrechits_id[imod][view][pln][ch]=0;
    }
  }


}


void fMode(int mode){
  if(mode==1)
    cout<<"CCQE"<<endl;
  else if(mode==11)
    cout<<"CC1pi+"<<endl;
  else if(mode==12)
    cout<<"CC1pi0"<<endl;
  else if(mode==13)
    cout<<"CC1pi-"<<endl;
  else if(mode==16)
    cout<<"CCcohpi"<<endl;
  else if(mode==17)
    cout<<"CC1gamma"<<endl;
  else if(mode==21)
    cout<<"CCmultipi"<<endl;
  else if(mode==22)
    cout<<"CC1eta"<<endl;
  else if(mode==23)
    cout<<"CC1K"<<endl;
  else if(mode==26)
    cout<<"CCDIS"<<endl;

  else if(mode==31)
    cout<<"NC1pi0"<<endl;
  else if(mode==32)
    cout<<"NC1pi0"<<endl;
  else if(mode==33)
    cout<<"NC1pi-"<<endl;
  else if(mode==34)
    cout<<"NC1pi+"<<endl;
  else if(mode==36)
    cout<<"NCcohpi"<<endl;
  else if(mode==37)
    cout<<"NC1gamma"<<endl;
  else if(mode==38)
    cout<<"NC1gamma"<<endl;
  else if(mode==41)
    cout<<"NCmultipi"<<endl;
  else if(mode==42)
    cout<<"NC1eta"<<endl;
  else if(mode==43)
    cout<<"NC1eta"<<endl;
  else if(mode==44)
    cout<<"NC1K0"<<endl;
  else if(mode==45)
    cout<<"NC1K+"<<endl;
  else if(mode==46)
    cout<<"NCDIS"<<endl;
  else if(mode==51)
    cout<<"NCE"<<endl;
  else if(mode==52)
    cout<<"NCE"<<endl;
  else 
    cout<<mode<<endl;

}

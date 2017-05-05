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
#include <TRandom3.h>

#include "TApplication.h"
//#include "PMrecon.hxx"
//#include "Lolirecon.hxx"
#include "INGRID_Dimension.cxx"
//#include "PMdisp.h"
//#include "Lolidisp.h"

#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "PMReconSummary.h"
#include "NeutInfoSummary.h"

//#define  STUDY_BADCHANNEL
vector<int>  bad_mod;
vector<int>  bad_pln;
vector<int>  bad_view;
vector<int>  bad_ch;
void Add_Bad_Channel();
bool Is_Bad_Channel(IngridHitSummary* thit);
void fMode(int mode);

//define in IngridConstants.h
//const static char *dst_data          = "/home/t2kingrid/data/dst";
const static char *out_data          = "/home/t2kingrid/data/dst";
//const static char *cosmic_data       = "/home/t2kingrid/data/cosmic";

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
  //Int_t Nmod = 17;
  Int_t Nmod = 16;
  //Int_t Nmod = 14;
  Int_t Nini = 0;
  bool disp = false; 
  bool cosmic = false;
  Int_t seed =0;
  double crossrate_submod0  =0.0154; // from TN 7.3.3
  double crossrate_submod123=0.0067; // from TN 7.3.3
  while ((c = getopt(argc, argv, "r:s:f:cado:i:t:y:z:")) != -1) {
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
    case 't':
      seed=atoi(optarg);
      break;
    case 'y':
      crossrate_submod0=atof(optarg);
      break;
    case 'z':
      crossrate_submod123=atof(optarg);
      break;

    }
  }
  FileStat_t fs;
  // ifstream timing;

  if(!renamef){
    sprintf(buff1,"%s/ingrid_%08d_%04d_form.root",
	    dst_data, run_number,sub_run_number);
  }
  if(cosmic){
    sprintf(buff1,"%s/ingrid_%08d_%04d_form.root",
	    cosmic_data, run_number,sub_run_number);
  }
  
  //######## read root file #############################
  //#####################################################
  cout<<"reading "<<buff1<<"....."<<endl;
  if(gSystem->GetPathInfo(buff1,fs)){
    cout<<"Cannot open file "<<buff1<<endl;
    exit(1);
  }

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
  if( !rename ){
    sprintf(Output, "%s/ingrid_%08d_%04d_recon.root", 
	    out_data, run_number, sub_run_number); 
  }
  if(cosmic){
    sprintf(Output, "%s/ingrid_%08d_%04d_recon.root", 
	    cosmic_data, run_number, sub_run_number); 
  }

  std::cout << buff1 << std::endl;
  std::cout << Output << std::endl;


  TFile*            wfile = new TFile(Output, "recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* wsummary = new IngridEventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary", 
				  &wsummary,  64000,  99);
  IngridHitSummary*    inghitsum;
  IngridSimVertexSummary*   simver;
  //IngridBasicReconSummary*    ingbasicrecon;
  //PMReconSummary*    ingbasicrecon;
  //IngridHitSummary*  inghitsum_t;
  //IngridBasicReconSummary* recon = new IngridBasicReconSummary();
  PMReconSummary* recon = new PMReconSummary();
  //Hit                        hit;

  Initialize_INGRID_Dimension();

  
  const int VIEW=2;
  const int PLN=8;
  const int CH=80;
  INGRID_Dimension *fdim_temp = new INGRID_Dimension();
  TRandom3 r;
  r.SetSeed(seed);

  for(int ievt=Nini; ievt<nevt; ievt++){

    if(ievt%100==0)cout << "analyze event# " << ievt<<endl;
    //    cout << "analyze event# " << ievt<<endl;
    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);

    for( int cyc=Scyc; cyc<Ncyc; cyc++ ){  //### Cycle Loop
      for( int mod=15; mod<16; mod++ ){   //### Module Loop
	
	int ninghit = evt -> NIngridModHits(mod, cyc);
        double pe_cross[VIEW][PLN][CH];
        for(int view=0;view<VIEW;view++){
	  for(int pln =0;pln <PLN ;pln++){
	    for(int ch  =0;ch<CH    ;ch++){
	      pe_cross[view][pln][ch]=0;
	    }
	  }
        }

	for(int i=0; i<ninghit; i++){
	  inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(i, mod, cyc) );
	  if(inghitsum->cyc==-2)continue; //This is killed hit 

	  int view_s = inghitsum -> view;
	  int pln_s  = inghitsum -> pln;
	  int ch_s   = inghitsum -> ch;
	  double pe_s   = inghitsum -> pe;
          double posxy_s[2]; // in mm
	  //ML modif (corr) 2017/03/21
	  if(view_s==0){
	    if(inghitsum->gridcell_id_y1>20)posxy_s[0] = (inghitsum -> gridcell_id_y1-21)*5;
	    else posxy_s[0] = inghitsum -> gridcell_id_y1*5;
	    posxy_s[1] = inghitsum -> xy;
	  }
	  else {
	    posxy_s[0] = inghitsum -> xy;
	    if(inghitsum->gridcell_id_x1>20) posxy_s[1] = (inghitsum -> gridcell_id_x1-21)*5;
	    else posxy_s[1] = inghitsum -> gridcell_id_x1*5;
	  }
	  posxy_s[view_s]-=50; // to have center at 0cm
	    //double posz_s = inghitsum -> gridcell_id_y1;


          int view;        
          if(view_s==0)view=1;
          else if(view_s==1)view=0;
          int ch_min, ch_max;
          if(ch_s>=40 && ch_s<60)     {ch_min=40;ch_max=60;}
          else if(ch_s>=60 && ch_s<80){ch_min=60;ch_max=80;}
          else if(ch_s<40)continue;
	  //cout<<"reference hit (view,pln,ch,pe)="<<view_s<<" "<<pln_s<<" "<<ch_s<<" pe="<<pe_s<<endl;
	  //cout<<"position x="<<posxy_s[0]<<" y="<<posxy_s[1]<<" grid_id: x="<<inghitsum->gridcell_id_y1<<" y="<<inghitsum->gridcell_id_x1<<endl;
          for(int ch=ch_min;ch<ch_max;ch++){
            double posxy[3];
            fdim_temp->get_pos_loli(mod,view,pln_s,ch,&posxy[0],&posxy[1],&posxy[2]);//cm
	    //	    cout<<ch<<" view_s="<<view_s<<" posxy_s="<<posxy_s[view_s]<<" posxy="<<posxy[view_s]<<endl;
	    cout<<posxy_s[view_s]<<" "<<posxy[view_s]<<" "<<fabs(posxy_s[view_s]-posxy[view_s])<<endl;
            if( fabs(posxy_s[view_s]-posxy[view_s]) <5. && pln_s<2){ pe_cross[view][pln_s][ch]+= r.Poisson(pe_s * crossrate_submod0);cout<<pe_cross[view][pln_s][ch]<<endl;}//poisson
            if( fabs(posxy_s[view_s]-posxy[view_s]) <5. && pln_s>=2)pe_cross[view][pln_s][ch]+= r.Poisson(pe_s * crossrate_submod123);//poisson
          }
        }

	for(int i=0; i<ninghit; i++){
	  inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(i, mod, cyc) );
	  if(inghitsum->cyc==-2)continue; //This is killed hit 
          for(int view=0;view<VIEW;view++){
	    for(int pln =0;pln <PLN ;pln++){
	      for(int ch  =0;ch<CH    ;ch++){
		if( inghitsum->mod == mod && inghitsum->view == view && inghitsum->pln == pln && inghitsum->ch == ch && pe_cross[view][pln][ch]>0){
		  inghitsum->pe_cross += pe_cross[view][pln][ch];
		  //cout<<"\ncross-talk added to hit (view,pln,ch)="<<view<<" "<<pln<<" "<<ch<<", pe_cross="<<pe_cross[view][pln][ch]<<endl;
		  pe_cross[view][pln][ch]=0;
		}
	      }
	    }
          }
        }
        for(int view=0;view<VIEW;view++){
	  for(int pln =0;pln <PLN ;pln++){
	    for(int ch  =0;ch<CH    ;ch++){
	      if(pe_cross[view][pln][ch]>0){
		IngridHitSummary* hit = new IngridHitSummary();
		IngridSimHitSummary* simhit = new IngridSimHitSummary();
		hit->mod  = mod;	
		hit->view = view;	
		hit->pln  = pln;	
		hit->ch   = ch;
		hit->time = inghitsum->time;	
		hit->timecorr=inghitsum->timecorr;
		hit->pe = 0;
		hit->pecorr=0;
		hit->pe_cross = pe_cross[view][pln][ch];
		double xy, z;
		fdim_temp-> get_pos_loli_xy( mod, view, pln, ch, &xy, &z);
		inghitsum -> xy    = xy;
		inghitsum ->  z    =  z;

		if(Is_Bad_Channel( hit ))continue;
		evt->AddIngridModHit(hit, mod, cyc);
		evt->AddIngridSimHit(simhit);
		//cout<<"\ncross-talk added, new hit (view,pln,ch)="<<view<<" "<<pln<<" "<<ch<<", pe_cross="<<pe_cross[view][pln][ch]<<endl;
	      }
	    }
	  }
        }

      }//cyc
    }//mod

    wsummary = evt;
    wtree -> Fill();
  }//Event Loop

  if(fdim_temp) delete fdim_temp;

  //######## Write and Close ####################
  wfile->cd();
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

}


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



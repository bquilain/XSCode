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
  for(int ievt=Nini; ievt<nevt; ievt++){

    if(ievt%100==0)cout << "analyze event# " << ievt<<endl;
    //    cout << "analyze event# " << ievt<<endl;
    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);

    for( int cyc=Scyc; cyc<Ncyc; cyc++ ){  //### Cycle Loop
      for( int mod=15; mod<16; mod++ ){   //### Module Loop
	
	int ninghit = evt -> NIngridModHits(mod, cyc);
	cout<<"ninghit="<<ninghit<<endl;
	double pe_gridcell_x[42][8];
	double pe_gridcell_y[42][8];
	for(int i=0;i<42;i++){
		for(int j=0;j<8;j++){
			pe_gridcell_x[i][j]=0;
			pe_gridcell_y[i][j]=0;
		}
	}
	for(int i=0; i<ninghit; i++){
	  inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(i, mod, cyc) );
	  if(inghitsum->cyc==-2)continue; //This is killed hit 

	  if(inghitsum -> gridcell_id_x1!=-1){
	  	pe_gridcell_x[inghitsum -> gridcell_id_x1][inghitsum -> pln]+=inghitsum -> pe;
	  }
	  if(inghitsum -> gridcell_id_x2!=-1){
	  	pe_gridcell_x[inghitsum -> gridcell_id_x2][inghitsum -> pln]+=inghitsum -> pe;
	  }
	  if(inghitsum -> gridcell_id_y1!=-1){
	  	pe_gridcell_y[inghitsum -> gridcell_id_y1][inghitsum -> pln]+=inghitsum -> pe;
	  }
	  if(inghitsum -> gridcell_id_y2!=-1){
	  	pe_gridcell_x[inghitsum -> gridcell_id_y2][inghitsum -> pln]+=inghitsum -> pe;
	  }
	}

	TRandom3 r;
	for(int i=0;i<42;i++){
		for(int j=0;j<8;j++){
			if(pe_gridcell_x[i][j]!=0){
	    			INGRID_Dimension *fdim_temp = new INGRID_Dimension();
				int cross_n=0;
				int cross_ch[10];
          		        fdim_temp -> get_loli_scintiid_from_cellid( mod, 0, j, i, &cross_n, cross_ch); //mod view pln gridcellid crossn crossch
				for(int k=0;k<cross_n;k++){
				        ninghit = evt -> NIngridModHits(mod, cyc);
					for(int l=0;l<ninghit;l++){
	  					inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(l, mod, cyc) );
						if(inghitsum->mod == mod && inghitsum->view == 0 && inghitsum->pln == j && inghitsum->ch == cross_ch[k]  ){
							inghitsum->pe_cross += r.Poisson(pe_gridcell_x[i][j] * 0.02);//poisson
						}
						else if(l==ninghit-1){
							IngridHitSummary* hit = new IngridHitSummary();
							IngridSimHitSummary* simhit = new IngridSimHitSummary();
							hit->mod = mod;	
							hit->view = 0;	
							hit->pln = j;	
							hit->ch = cross_ch[k];
							if(cross_ch[k]<-1 || cross_ch[k]>79)std::cout << cross_ch[k] << std::endl;
							hit->time = inghitsum->time;	
							hit->pe = 0;
							hit->pe_cross = r.Poisson(pe_gridcell_x[i][j] * 0.02);//poisson;
							if(hit->pe_cross==0)continue;
							cout<<"xtalk hit added"<<endl;
							evt->AddIngridModHit(hit, mod, cyc);
							evt->AddIngridSimHit(simhit);
						}
					}
				}
			}
			if(pe_gridcell_y[i][j]!=0){
	    			INGRID_Dimension *fdim_temp = new INGRID_Dimension();
				int cross_n=0;
				int cross_ch[10];
          		        fdim_temp -> get_loli_scintiid_from_cellid( mod, 1, j, i, &cross_n, cross_ch); //mod view pln gridcellid crossn crossch
				for(int k=0;k<cross_n;k++){
				        ninghit = evt -> NIngridModHits(mod, cyc);
					for(int l=0;l<ninghit;l++){
	  					inghitsum   = (IngridHitSummary*) (evt -> GetIngridModHit(l, mod, cyc) );
						if(inghitsum->mod == mod && inghitsum->view == 1 && inghitsum->pln == j && inghitsum->ch == cross_ch[k]  ){
							inghitsum->pe_cross += r.Poisson(pe_gridcell_y[i][j] * 0.02);//poisson
						}
						else if(l==ninghit-1){
							IngridHitSummary* hit = new IngridHitSummary();
							IngridSimHitSummary* simhit = new IngridSimHitSummary();
							hit->mod = mod;	
							hit->view = 1;	
							hit->pln = j;	
							hit->ch = cross_ch[k];	
							hit->time = inghitsum->time;	
							hit->pe = 0;
							hit->pe_cross = r.Poisson(pe_gridcell_y[i][j] * 0.02);//poisson;
							if(hit->pe_cross==0)continue;
							evt->AddIngridModHit(hit, mod, cyc);
							evt->AddIngridSimHit(simhit);
						}
					}
				}
			}
		}
	}

      }//cyc
     }//mod

    wsummary = evt;
    wtree -> Fill();
  }//Event Loop

  //######## Write and Close ####################
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

}


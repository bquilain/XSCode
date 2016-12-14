//##### Standard C++ lib. ######
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>
#include <math.h> 
//##### Root Library ###########
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
#include <TSystem.h>
//##### INGRID Library #########
#include "INGRIDEVENTSUMMARY.h"
#include "BeamInfoSummary.h"
//##### BSD acc. Library #######
//#include "./bsd/BeamData.h"
#include "./bsd_new/BeamData/BeamData.h"

BEAM_DATA bsd;
BeamData& beam = BeamData::instance();
bool GoodSpillv01();
bool GoodSpill();
bool Horn320();
bool HornOff();


int main(int argc, char *argv[]){
  TROOT        root ("GUI","GUI");
  TApplication app  ("App",0,0);

  //beam.init_db("","");
  beam.init_db();

  int    run_number;
  int             c  =     -1;
  int    sub_run_number=0;
  bool   hornoff     = false;
  bool   horn320     = false;
  bool   v01         = false; //### using beam summry not quick summary
  while ((c = getopt(argc, argv, "r:s:voa")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      sub_run_number=atoi(optarg);
      break;
    case 'v':
      v01 = true;
      break;
    case 'o':
      hornoff = true;
      break;
    case 'a':
      horn320 = true;
      break;
    }
  }
  //###### Read INGRID Data #######
  //###############################
  IngridEventSummary* summary  = new IngridEventSummary();

  FileStat_t fs;
  char buff[300];
  sprintf(buff, "/export/scraid2/data/bquilain/Datafiles/ingrid_%08d_%04d_calib.root", run_number, sub_run_number);
  if(gSystem->GetPathInfo(buff, fs ) ){
    cout<<"Cannot open file "<<Form("/export/scraid2/data/bquilain/Datafiles/ingrid_%08d_%04d_calib.root", run_number, sub_run_number) << endl;
    exit(1);
  }

  cout<<"read..."<<endl;
  TFile*            rfile      = new TFile( Form("/export/scraid2/data/bquilain/Datafiles/ingrid_%08d_%04d_calib.root", run_number, sub_run_number),
					"read");

  TTree*             tree      = (TTree*)rfile -> Get("tree");
  TBranch*          EvtBr      = tree->GetBranch("fDefaultReco.");
  EvtBr                        ->SetAddress(&summary);
  tree                         ->SetBranchAddress("fDefaultReco.", &summary);
  int                nevt      = (int)tree -> GetEntries();
  cout << "Total # of events = " << nevt <<endl;
 

  //###### New INGRID Data ########
  //###############################
  BeamInfoSummary*        bsum  = new BeamInfoSummary();

  TFile*                 wfile  = new TFile( Form("/export/scraid2/data/bquilain/Datafiles/ingrid_%08d_%04d_wbsd.root", run_number, sub_run_number),
					"recreate");
  TTree*                 wtree  = new TTree("tree", "tree");
  IngridEventSummary* wsummary  = new IngridEventSummary(); 
  wtree                         -> Branch   ("fDefaultReco.",
					     "IngridEventSummary", 
					     &wsummary,  64000,   99);
  wtree -> SetMaxTreeSize(10000000000);

  ofstream fIngGood("INGRID_Good_Spill_List.txt");
  for(int ievt = 0; ievt < nevt; ++ievt){  
    if(ievt%1000==0) 
	cout<<"Run#"<<run_number<<" subrun#"<<sub_run_number<<"\tprocess event:"<<ievt<<endl;
    wsummary -> Clear("C");
    summary  -> Clear("C");
    bsum     -> Clear("C");
    tree     -> GetEntry(ievt);

    int nd280nspill = summary -> nd280nspill;
    int       utime = summary -> time;

    if( ! beam.get_spill( nd280nspill, utime, 10 ) ){
      cout << "Error in looking for spill " << nd280nspill 
	   << " "                           << utime 
	   << endl;
      continue;
    }
    int dt = abs ( beam.trigger_time_sec(0) - utime );

    if( dt > 10 ){
      cout << "Error in looking for spill " << nd280nspill 
 	   << " "                           << utime 
	   << " Time difference is too big "<< dt
	   << endl;
      continue; 
    }

    bool ok = false;
    if( v01 && GoodSpillv01() )
      ok = true;
    else if( !v01 && GoodSpill() && !horn320 && !hornoff )
      ok = true;
    else if( horn320 && Horn320() )
      ok = true;
    else if( hornoff && HornOff() )
      ok = true;

    if(!ok)
      continue;

    cout << nd280nspill << endl;
    bsd = beam.data();
    //### Set bunch_flag ####
    for(int cyc=0; cyc<23; cyc++){
      if(4<=cyc&&cyc<=9){
	if( bsd.ct_np [4][ cyc - 4 + 1] > 1e11 )  //Use CT05, all bunches
	  {
	    summary -> bunch_flag[cyc] = true;
	  }
	else summary -> bunch_flag[cyc] = false;
      }
      else summary -> bunch_flag[cyc]=false;
    }
    //#### Fill ########
    //##################
    bsum -> nurun           = bsd.nurun;
    bsum -> spill_flag      = bsd.spill_flag;
    bsum -> good_spill_flag = bsd.good_spill_flag;
    bsum -> trg_sec         = bsd.trg_sec[0];
    bsum -> spillnum        = bsd.spillnum;
    bsum -> nd280spillnum   = ( (bsd.spillnum) & 0xffff ) + 1;
    bsum -> good_spill_flag = bsd.good_spill_flag;
    bsum -> run_type        = bsd.run_type;
    for(int l=0; l<3; l++){
      bsum -> target_eff[l] = bsd.target_eff[l];
    }
    for(int l=0; l<5; l++){
      for(int k=0; k<9; k++){
	bsum -> beam_time[l][k] = bsd.beam_time[l][k];
      }
    }
    for(int l=0; l<5; l++){
      for(int k=0; k<9; k++){
	bsum -> ct_np    [l][k] = bsd.ct_np    [l][k];
      }
    }
    for(int l=0; l<3; l++){
      for(int k=0; k<5; k++){
	bsum -> hct[l][k]       = bsd.hct[l][k];
      }
    }
    for(int l=0; l<12; l++){
      bsum -> mumon[l]          = bsd.mumon[l];
    }
    for(int l=0; l<13; l++){
      bsum -> otr[l]            = bsd.otr[l];
    }
    summary -> AddBeamSummary(bsum);
    wsummary = summary;
    wtree   -> Fill();
  }
  wtree   -> Write();
  wfile   -> Write();
  wfile   -> Close();
}


bool GoodSpillv01(){
  if( beam.good_spill() )
    return true;
  else
    return false;
}

bool GoodSpill(){

  if( ! beam.spill_flag() )
    return false;

  //#### Bad Spill at Run31 #####   
  if( 1269071040  <= beam.trigger_time_sec(0) && beam.trigger_time_sec(0) <= 1269071460  )
    return false;
  if( 1269097200  <= beam.trigger_time_sec(0) && beam.trigger_time_sec(0) <= 1269124200  )
    return false;
  if( 1269143340  <= beam.trigger_time_sec(0) && beam.trigger_time_sec(0) <= 1269143460  )
    return false;


  //### bad spill (temporary method) ####
  if( 598833  <= beam.spill_number() && beam.spill_number() <= 598837  &&
      beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
    return false;
  if( 598937  <= beam.spill_number() && beam.spill_number() <= 598941  &&
      beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
    return false;
  if( 1201242 <= beam.spill_number() && beam.spill_number() <= 1201270 &&
      beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
    return false;
  if( 1202461 <= beam.spill_number() && beam.spill_number() <= 1202533 &&
      beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
    return false;
  if( 1234283 <= beam.spill_number() && beam.spill_number() <= 1234304 &&
      beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
    return false;
  if( 1235436 <= beam.spill_number() && beam.spill_number() <= 1235450 &&
      beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
    return false;
  bool normal_cond = ( beam.run_type() == 1 ) &&
    ( fabs ( beam.horn_current(0) - 250 ) < 5 ) && 
    ( fabs ( beam.horn_current(1) - 250 ) < 5 ) &&
    ( fabs ( beam.horn_current(2) - 250 ) < 5 );
  
  bool beamtrg     =
    ( beam.ct_np(4) > 1e11 );
  
  bool mucenter    = 
    ( sqrt( beam.mumon_si_x() * beam.mumon_si_x() + beam.mumon_si_y() * beam.mumon_si_y() ) < 10 );
  
  bool muct        =
    ( 31.7 < ( beam.mumon_si_totq() / beam.ct_np(4) * 1e9 ) ) &&
    ( ( beam.mumon_si_totq() / beam.ct_np(4) * 1e9 ) < 35.1 );
  
  if( normal_cond && beamtrg && mucenter && muct )
    return true;
  else
    return true;
}

bool Horn320(){
  cout << beam.horn_current(0) << endl;
  if( beam.ct_np(4) > 1e11 ){
    if( ( fabs ( beam.horn_current(0) - 320 ) < 5 ) &&
	( fabs ( beam.horn_current(1) - 320 ) < 5 ) &&
	( fabs ( beam.horn_current(2) - 320 ) < 5 )
	){
      return true;
    }
  }
  return false;
}
bool HornOff(){
  if( beam.ct_np(4) > 1e11 ){
    if( ( fabs ( beam.horn_current(0) - 0 ) < 5 ) &&
	( fabs ( beam.horn_current(1) - 0 ) < 5 ) &&
	( fabs ( beam.horn_current(2) - 0 ) < 5 )
	){
      return true;
    }
  }
  return false;
}

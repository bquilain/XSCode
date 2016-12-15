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
#include "./bsd/BeamData/BeamData.h"
//#include "./bsd/event.h"

int main(int argc, char *argv[]){
  TROOT        root ("GUI","GUI");
  TApplication app  ("App",0,0);
  BeamData& beam = BeamData::instance();
  //beam.init_db("","");
  beam.init_db();

  int    run_number;
  int             c  =     -1;
  int       no_beam  =  false;
  int    sub_run_number=0;
  while ((c = getopt(argc, argv, "r:s:n")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      sub_run_number=atoi(optarg);
      break;
    case 'n':
      no_beam = true;
      break;
    }
  }
  //###### Read INGRID Data #######
  //###############################
  IngridEventSummary* summary  = new IngridEventSummary();

  FileStat_t fs;
  char buff[300];
  bool ok = false;
  sprintf(buff, "/home/daq/data/dst/ingrid_%08d_%04d.root", run_number, sub_run_number);
  if(!gSystem->GetPathInfo(buff, fs ) ){
    cout<<"open file "<< buff << endl;
    goto Find;
  }
  sprintf(buff, "/media/disk/Run29_root_file/ingrid_%08d_%04d.root", run_number, sub_run_number);
  if(!gSystem->GetPathInfo(buff, fs ) ){
    cout<<"open file "<< buff << endl;
    goto Find;
  }
  sprintf(buff, "/media/disk/Run30_root_file/ingrid_%08d_%04d.root", run_number, sub_run_number);
  if(!gSystem->GetPathInfo(buff, fs ) ){
    cout<<"open file "<< buff << endl;
    goto Find;
  }
  
  exit(1);

 Find:


  //TFile*            rfile      = new TFile( Form("/home/daq/data/dst/ingrid_%08d_%04d.root", run_number, sub_run_number),
  TFile*            rfile      = new TFile( buff, "read" );


  TTree*             tree      = (TTree*)rfile -> Get("tree");
  TBranch*          EvtBr      = tree->GetBranch("fDefaultReco.");
  EvtBr                        ->SetAddress(&summary);
  tree                         ->SetBranchAddress("fDefaultReco.", &summary);
  int                nevt      = (int)tree -> GetEntries();
 
 

  BEAM_DATA bsd;
  ofstream fIngGood("INGRID_Good_Spill_List.txt");
  for(int ievt = 0; ievt < nevt; ++ievt){  
    if(ievt%1000==0)cout<<"Run#"<<run_number<<"\tprocess event:"<<ievt<<endl;
    tree     -> GetEntry(ievt);

    int nd280nspill = summary -> nd280nspill;
    int       utime = summary -> time; 

    if( ! beam.get_spill( nd280nspill, utime ) ){
      /*
      cout << "Error in looking for spill " << nd280nspill 
	   << " "                           << utime 
	   << endl;
      */
      continue;
    }
    int dt = abs ( beam.trigger_time_sec(0) - utime );

    if( dt > 10 ){
      /*
      cout << "Error in looking for spill " << nd280nspill 
 	   << " "                           << utime 
	   << " Time difference is too big "<< dt
	   << endl;
      */
      continue; 
   }

    if( ! beam.spill_flag() )
      continue;

    //### In order to define good spill@INGRID it self, ###
    //### dump the nurun, spill# and trg_sec[0]         ###
    fIngGood << beam.beam_runnumber()    << "\t"
	     << beam.spill_number()      << "\t"
	     << beam.trigger_time_sec(0) << "\t"
	     << endl;



    //#### Bad Spill at Run31 #####   
    if( 1269071040  <= beam.trigger_time_sec(0) && beam.trigger_time_sec(0) <= 1269071460  )
      continue;
    if( 1269097200  <= beam.trigger_time_sec(0) && beam.trigger_time_sec(0) <= 1269124200  )
      continue;
    if( 1269143340  <= beam.trigger_time_sec(0) && beam.trigger_time_sec(0) <= 1269143460  )
      continue;


    //### bad spill (temporary method) ####
    if( 598833  <= beam.spill_number() && beam.spill_number() <= 598837  )
      continue;
    if( 598937  <= beam.spill_number() && beam.spill_number() <= 598941  )
      continue;
    if( 1201242 <= beam.spill_number() && beam.spill_number() <= 1201270 )
      continue;
    if( 1202461 <= beam.spill_number() && beam.spill_number() <= 1202533 )
      continue;
    if( 1234283 <= beam.spill_number() && beam.spill_number() <= 1234304 )
      continue;
    if( 1235436 <= beam.spill_number() && beam.spill_number() <= 1235450 )
      continue;
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

    if( !(normal_cond && beamtrg && mucenter && muct) )
      continue;


  }
  fIngGood.close();
}


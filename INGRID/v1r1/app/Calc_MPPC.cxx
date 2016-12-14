//under construction
//////////////////////////////////////////////////
//make calibration constant at offline
//Made by Masashi Otani
//////////////////////////////////////////////////

/*
|||||||||||||||||||||||||||||||
consists 
1. main
2. ProcessFile
3. Read
4. 
main{
ProcessFile
->Read->Analysis(each event)
|||||||||||||||||||||||||||||||
*/

// ND280 software includes
#include "TMidasBank.hxx"
#include "TMidasFile.hxx"
#include "TMidasBankProxy.hxx"
#include "TMidasBankProxyRegistry.hxx"
#include "TND280RawEvent.hxx"
#include "TRawDataHeader.hxx"
#include "TRunInfoBank.hxx"

// oaRawEvent includes
#include "TTripTDigitBank.hxx"
#include "TMidasTripTDigitItr.hxx"
#include "TMidasTripTDigit.hxx"

#include "TMCMBank.hxx"
#include "TTriggerBank.hxx"

// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include <stdlib.h>

#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>
#include <sys/stat.h>
using namespace std;
// INGRID includes
#include "setup.hxx"
#include "ana_MPPC.cxx"
#include "INGRID_Ch_config.cxx"
#include "INGRID_BadCh_mapping.cxx"
#include <unistd.h> // using getopt      
FileStat_t            fs;

int                   run_number;
int                   sub_run_number;
Int_t                 NumEvt;
ana_MPPC*             fana_MPPC;
INGRID_Ch_config*     fINGRID_Ch_config;
INGRID_BadCh_mapping* fINGRID_BadCh_mapping;

//
const char* logain_dir="/home/ingrid/data_process/MPPC_calib/mppc_calib_table/logain";
const char* dir0="/home/ingrid/data_process/MPPC_calib/mppc_calib_table/0000";
const char* dir1="/home/ingrid/data_process/MPPC_calib/mppc_calib_table/0001";
const char* dir2="/home/ingrid/data_process/MPPC_calib/mppc_calib_table/0002";
const char* dir3="/home/ingrid/data_process/MPPC_calib/mppc_calib_table/0003";
const char* table_dir="/home/ingrid/data_process/MPPC_calib/mppc_calib_table/table";

bool     cond_flag;     //### if # of events < anabreak, false
ofstream calibfile;     //### New(2010/3/3) format for MPPC calibration
                        //### constant table.  
ifstream read_calibfile;
ofstream real_calibfile;
ofstream MPPCfile;      //### High Gain pedestal and gain
ofstream noisefile;     //### noise rate(High Gain)
ofstream sigmafile;
ofstream Tdcthreshold;  //### Tdc Threshold(High Gain)
ofstream LowGain;       //### Low Gain pedestal

TH1F*  fH_HighAdc[NumMod][2][NumTFB][NumCh];
TH1F*  fH_LowAdc[NumMod][2][NumTFB][NumCh];
Int_t  cAnaTrg;
Int_t  fMinAdcwTdcCut[NumMod][2][NumTFB][NumCh];
Int_t  anacounter; //
Int_t  anabreak;   //number of event for analysis

int StartTime;
int EndTime;


void Analysis(){
  cout<<"analysis"<<endl;
  //TCanvas *c1 = new TCanvas("c1","c1",10,10,500,500);
  char buff1[300];
  double pedestal,gain,noise,caa,lopedestal;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    cout<<"nummod:"<<nummod<<endl;
    for(Int_t view=0;view<2;view++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      	for(Int_t numch=0;numch<24;numch++){
	  int global_ch;

	  bool tpl,veto;
	  int tmp;int tmpch;
	  if(numtfb>10){tpl=false;veto=true;tmp=numtfb-11;}
	  else {tpl=true;veto=false;tmp=numtfb;}
	  if(view==0&&tpl){tmpch=numch+24;}
	  else if(view==1&&tpl){tmpch=numch;}
	  else if(veto){tmpch=numch;}



	  if(fINGRID_Ch_config->get_global_ch(&nummod,&view,&numtfb,&numch,&global_ch)&&!fINGRID_BadCh_mapping->badchannel(&nummod, &tmp, &tmpch, &tpl, &veto)){
	    if(!fana_MPPC->analysis_old_version(fH_HighAdc[nummod][view][numtfb][numch])){
	    

	    };
	    fana_MPPC -> analysis_logain(fH_LowAdc[nummod][view][numtfb][numch]);
	    //fH_LowAdc[nummod][view][numtfb][numch]->Draw();
	    //c1->Update();
	    //cout<<fana_MPPC->get_lowpedestal()<<endl;
	    //cin.get();


	    pedestal  = fana_MPPC->get_pedestal();
	    lopedestal= fana_MPPC->get_lowpedestal();
	    gain      = fana_MPPC->get_gain();
	    Int_t Nentry = fH_HighAdc[nummod][view][numtfb][numch]->GetEntries();
	    noise     = fana_MPPC->get_noise(Nentry);
	    caa       = fana_MPPC->get_crosstalk_and_afterpulse(Nentry);
	    double pedsigma = fana_MPPC->get_pedestal_sigma();

	    if(fabs(gain - 10)>6){
	      cout<<"Mod :"   << nummod	     
		  <<"view :"  << view
		  <<"pln :"   << numtfb
		  <<"ch :"    << numch
		  <<"gain :"  << gain
		  <<endl;
	      //fH_HighAdc[nummod][view][numtfb][numch]->Draw();
	      //c1->Update();
	      //cin.get();

	    }

	    if(noise > 1.5){
	      //fH_HighAdc[nummod][view][numtfb][numch]->Draw();
	      //c1->Update();
	      //cin.get();

	    }

	    bool tpl,veto;
	    int tmp;int tmpch;
	    if(numtfb>10){tpl=false;veto=true;tmp=numtfb-11;}
	    else {tpl=true;veto=false;tmp=numtfb;}
	    if(view==0){tmpch=numch+24;}
	    else{tmpch=numch;};
	    

	   
	    if(!fINGRID_BadCh_mapping->badchannel(&nummod, &tmp, &tmpch, &tpl, &veto)){
	    MPPCfile<<nummod   <<"\t"
                    <<view     <<"\t"
                    <<numtfb   <<"\t"
                    <<numch    <<"\t"
                    <<pedestal <<"\t"
                    <<gain     <<endl;
	    noisefile<<nummod   <<"\t"
                     <<view     <<"\t"
                     <<numtfb   <<"\t"
                     <<numch    <<"\t"
                     <<noise    <<"\t"
                     <<caa      <<endl;
	    sigmafile<<nummod   <<"\t"
                     <<view     <<"\t"
                     <<numtfb   <<"\t"
                     <<numch    <<"\t"
                     <<noise    <<"\t"
                     <<pedsigma <<endl;
	    LowGain  <<nummod   <<"\t"
                     <<view     <<"\t"
                     <<numtfb   <<"\t"
                     <<numch    <<"\t"
                     <<lopedestal <<"\t"
                     <<gain/10   <<endl;

	    calibfile <<nummod     << "\t"
                      <<view       << "\t"
                      <<numtfb     << "\t"
                      <<numch      << "\t"
                      <<pedestal   << "\t"
                      <<lopedestal << "\t"
                      <<gain       << "\n";


	    Tdcthreshold
	             <<nummod   <<"\t"
		     <<view     <<"\t"
                     <<numtfb   <<"\t"
                     <<numch    <<"\t"
 		     <<fMinAdcwTdcCut[nummod][view][numtfb][numch]  <<"\t"
		     <<1.0*(fMinAdcwTdcCut[nummod][view][numtfb][numch]-pedestal)/gain  <<"\t"
		     <<endl;

	    }
	    else{
	    MPPCfile<<nummod   <<"\t"
                    <<view     <<"\t"
                    <<numtfb   <<"\t"
                    <<numch    <<"\t"
                    <<"100"    <<"\t"
                    <<"100"    <<endl;
	    noisefile<<nummod   <<"\t"
                     <<view     <<"\t"
                     <<numtfb   <<"\t"
                     <<numch    <<"\t"
                     <<"0.1"    <<"\t"
                     <<"0.01"   <<endl;
	    Tdcthreshold
	             <<nummod   <<"\t"
		     <<view     <<"\t"
                     <<numtfb   <<"\t"
                     <<numch    <<"\t"
                     <<noise    <<"\t"
		     <<-777     <<"\t"
		     <<-777     <<"\t"
		     <<endl;


	    }
	    
	    

	  
	  }//active channel
	}//numch
      }//numtfb
    }//view
    cout<<"gain:"<<gain<<endl;
  }//nummod


  //############# write tdc threshold ######################
  //########################################################

}


void Event(ND::TND280RawEvent* re) {

  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    NumEvt = runinfo.GetSeqNumber();
  }

  Int_t TrgId;
  ND::THandle<ND::TTriggerBank> triggerBank;
  while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
    ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
    TrgId = (trigger.GetTriggerWord()>>48)&0xffff;
  }

  if(TrgId==cAnaTrg){//Beam Event
    anacounter++;
    // Loop over all banks of type TTripTHitBank
    ND::THandle<ND::TTripTDigitBank> triptBank;
    while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {
      // Create an iterator over digits
      ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
      while ( ! itr.EOD() ) {
	ND::TMidasTripTDigit digit(itr.Get());
	Int_t mod, view, plane, ch;
	Int_t rmm      = digit.GetRMMNum();
	Int_t trip     = digit.GetTripTNum();
	Int_t trip_ch  = digit.GetChannelNum();
	Int_t tfb      = digit.GetTFBNum();
	if(fINGRID_Ch_config->channel_configuration(&rmm,&tfb,&trip,&trip_ch,&mod,&view,&plane, &ch)){  //### INGRID using channel
	  Int_t highadc  =  digit.GetHighGainADC();
	  fH_HighAdc[mod][view][plane][ch] -> Fill(highadc);
	  Int_t loadc    =  digit.GetLowGainADC();
	  fH_LowAdc[mod][view][plane][ch]  -> Fill(loadc);

	  long  tdc      =  digit.GetTimeOffset();
	  if(NumEvt>0 && highadc> 100&& tdc<16777201 && highadc < fMinAdcwTdcCut[mod][view][plane][ch]){
	    fMinAdcwTdcCut[mod][view][plane][ch] = highadc;
	  }//Tdc threshold calibration

	 
	}
      }   // End of loop over digits in this bank
    }     // End of loop over banks of digits in this event
  }//TrgId==1
}
//____________________________________________________________________________


//____________________________________________________________________________
void ProcessFile(const char *FileName) {
  
  cout<<"processing file..."<<endl;
  ND::TMidasBankProxyRegistry::Instance().Print();
  
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    return;
  }

  ND::TMidasFile mf;
  mf.Open(FileName);

  // Loop over events in this file


  while ( ND::TND280RawEvent* re = mf.ReadRawEvent() ) {
    if(NumEvt%10==0)cout<<NumEvt<<endl;


    re->PromoteMidasBanks(false);
   
    Event(re);
    if( NumEvt == 0 ){           //### Start of RUN
      ND::TRawDataHeader header = re->GetHeader();
      StartTime = header.GetTimeStamp();
    }
    if( anacounter > anabreak ){ //### End of Making CalibTable
      ND::TRawDataHeader header = re->GetHeader();
      EndTime = header.GetTimeStamp();
      cond_flag=true;
      delete re;
      break;
    }
    delete re;

  }  // End loop over events
}

//_____________________________________________________________________________

void Book(){
  char buff1[300];
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t view=0;view<2;view++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
	for(Int_t numch=0;numch<24;numch++){
	  if(view==0)sprintf(buff1,"Mod%02dXPln%02dCh%02d",nummod,numtfb,numch);
	  if(view==1)sprintf(buff1,"Mod%02dYPln%02dCh%02d",nummod,numtfb,numch);
	  fH_HighAdc[nummod][view][numtfb][numch] = new TH1F(buff1,buff1,200,100,300);
	  if(view==0)sprintf(buff1,"Mod%02dXPln%02dCh%02dLow",nummod,numtfb,numch);
	  if(view==1)sprintf(buff1,"Mod%02dYPln%02dCh%02dLow",nummod,numtfb,numch);
	  fH_LowAdc[nummod][view][numtfb][numch] = new TH1F(buff1,buff1,300,0,300);
	}//numch
      }//numtfb
    }//View
  }//nummod
}//Book
void Initialize(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t view=0;view<2;view++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
	for(Int_t numch=0;numch<24;numch++){
	  fMinAdcwTdcCut[nummod][view][numtfb][numch]=1000;
	}//numch
      }//numtfb
    }//View
  }//nummod
}//Initialize


//_____________________________________________________________________________
int main(int argc,char *argv[])
{
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);
  cond_flag=false;
  anacounter=0;
  cout<<"start"<<endl;
  int c=-1;
  char FileName[300],txt_file_name[300];
  Int_t pedestal_and_gain_for_cosmic;
  cAnaTrg=1;
  anabreak=300;
  NumEvt=0;
  sub_run_number = 0;
  while ((c = getopt(argc, argv, "hr:s:t:b:")) != -1) {
    switch(c){
    case 'b':
      anabreak = atoi(optarg);
      break;
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      sub_run_number=atoi(optarg);
      break;
    case 'h':
      cout<<"-r [run number]\t midas file to root file"<<endl;
      break;
    case 't':
      cAnaTrg=atoi(optarg);
      break;
    case '?':
      cout<<"Unknown option"<<endl;
      cout<<"-r [run number]\t midas file to root file"<<endl;
      exit(1);
      break;
    }
  }//option end

  //
  sprintf(FileName,"%s/ingrid_%08d_%04d.daq.mid.gz",
	  data_file_folder,
	  run_number,
	  sub_run_number);
  cout<<"file name is :"<<FileName<<endl;
  fana_MPPC = new ana_MPPC();
  fINGRID_Ch_config = new INGRID_Ch_config();
  fINGRID_BadCh_mapping = new INGRID_BadCh_mapping();

/*
  sprintf(txt_file_name,
	  "%s/ingrid_%08d_%04d_MPPCcalib.txt",
	  table_dir,
	  run_number,
	  sub_run_number);  

  cout<<"calib file name is :"<<txt_file_name<<endl;
*/

  //
 // calibfile.open(txt_file_name);

  calibfile.open("temp.txt");

  //
  sprintf(txt_file_name,"%s/ingrid_%08d_0000.txt",dir0,run_number);
  calibfile<<"\t\t\t\n";
  MPPCfile.open(txt_file_name);

  //
  sprintf(txt_file_name,"%s/ingrid_%08d_0000_logain.txt",logain_dir,run_number);
  LowGain.open(txt_file_name);

  //
  sprintf(txt_file_name,"%s/ingrid_%08d_0001.txt",dir1,run_number);
  noisefile.open(txt_file_name);

  //
  sprintf(txt_file_name,"%s/ingrid_%08d_0003.txt",dir3,run_number);
  sigmafile.open(txt_file_name);

  //
  sprintf(txt_file_name,"%s/ingrid_%08d_0002.txt",dir2,run_number);
  Tdcthreshold.open(txt_file_name);

  Initialize();
  Book();
  cout<<"read to midas file..."<<endl;
  ProcessFile(FileName);
  Analysis();

  /*
  calibfile.close();
  sprintf(txt_file_name,
	  "/home/daq/data/calib_table/ingrid_%08d_%04d_MPPCcalib.txt",
	  run_number,
	  sub_run_number);  
  
  calibfile.open(txt_file_name, ios::app);
  */

  calibfile.close();

  read_calibfile.open("temp.txt");

  sprintf(txt_file_name,
	  "%s/ingrid_%08d_%04d_MPPCcalib.txt",
	  table_dir,
	  run_number,
	  sub_run_number);  

  std::cout << "open : " << txt_file_name << std::endl;
  real_calibfile.open(txt_file_name);

  real_calibfile << cond_flag << "\t"
		 << StartTime << "\t"
		 << EndTime   << "\t"
		 << std::endl;

  int a,b,cc,d,f; 
  double e,g;

  while( read_calibfile >> a >> b >> cc >> d >> e >> f >> g ){
    //
    real_calibfile << a << "\t"
		   << b << "\t"
		   << cc << "\t"
		   << d << "\t"
		   << e << "\t"
		   << f << "\t"
		   << g << std::endl;
  }



  MPPCfile.    close();
  noisefile.   close();
  sigmafile.   close();
  Tdcthreshold.close();
  real_calibfile.close();

  return 0;
}



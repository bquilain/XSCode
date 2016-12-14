//////////////////////////////////////////////////
//////////////////////////////////////////////////
//   dump INGRID MIDAS file to txt file         //
//                                              //
//           Made by Masashi Otani              //
//                                              //
//                ver1.0 2009/9/14              //
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//

//|||||||||||||||||||||||||||||||
// consists 
//   1. main
//   2. ProcessFile
//   3. Read
//   4. Print
//   main{
//     ProcessFile{
//        Read{
//            Print{
//            }
//        }
//     }
//   }
//|||||||||||||||||||||||||||||||


// ND280 software includes
#include "TMidasBank.hxx"
#include "TMidasFile.hxx"
#include "TMidasBankProxy.hxx"
#include "TMidasBankProxyRegistry.hxx"
#include "TND280RawEvent.hxx"
#include "TRawDataHeader.hxx"

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
#include <stdlib.h>

//standard C++ libraly includes
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>
#include <sys/stat.h>
#include <unistd.h>
using namespace std;
FileStat_t fs;

// INGRID includes
#include "INGRID_Ch_config.cxx"
int run_number;
INGRID_Ch_config *fINGRID_Ch_config;
Int_t rmm, trip, trip_ch, tfb;
Int_t mod, plane, ch;
bool veto_flag, tpl_flag;
Int_t highadc, cycle;
Long_t Time;

void Print(){
  cout<<" RMM:"<<rmm<<" TRIP:"<<trip<<" Ch:"<<trip_ch;
  cout<<" Mod.:"<<mod;
  if(tpl_flag)cout<<" TPL:"<<plane;
  if(veto_flag)cout<<" VETO:"<<plane;
  cout<<" Cyc:"<<cycle<<" ADC:"<<highadc<<" TDC:"<<Time<<endl;
}

void Event(ND::TND280RawEvent* re) {
  // Loop over all banks of type TTripTHitBank
  ND::THandle<ND::TTripTDigitBank> triptBank;
  while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {
    // Create an iterator over digits
    ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
    while ( ! itr.EOD() ) {
      ND::TMidasTripTDigit digit(itr.Get());
      rmm = digit.GetRMMNum();
      trip = digit.GetTripTNum();
      trip_ch = digit.GetChannelNum();
      tfb     =  digit.GetTFBNum();
      if(fINGRID_Ch_config->channel_configuration(&rmm,&tfb,&trip,&trip_ch,&mod,&plane,&ch,&tpl_flag,&veto_flag)){
	highadc  =  digit.GetHighGainADC();
	Time    =  digit.GetTimeOffset();
	const ND::TRawDataHeader& h = re->GetHeader();
	int evno = h.GetSeqNo();    // Event Sequence Number
	int iint = digit.GetIntegrationNum();    // = Capacitor number
	int icoff = triptBank->GetTFBStartCycle();
	int cycle = iint - icoff;
	if (cycle<0) cycle += 23;
	Print();
      }
    }   // End of loop over digits in this bank
  }     // End of loop over banks of digits in this event
}

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
    re->PromoteMidasBanks(false);
    Event(re);
    delete re;
  }  // End loop over events
}

//_____________________________________________________________________________

int main(int argc,char *argv[]){
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);
  int c=-1;
  char FileName[300];
  while ((c = getopt(argc, argv, "hr:")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      sprintf(FileName,"/online/daqdata/ingrid_%08d_0000.daq.mid",run_number);
      cout<<"file name is :"<<FileName<<endl;
      break;
    case 'h':
      cout<<"-r [run number]\t midas file to root file"<<endl;
      break;
    case '?':
      cout<<"Unknown option"<<endl;
      cout<<"-r [run number]\t midas file to root file"<<endl;
      exit(1);
      break;
    }
  }//option end
  fINGRID_Ch_config = new INGRID_Ch_config();
  cout<<"start to midas to root..."<<endl;
  ProcessFile(FileName);
  return 0;
}


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

bool     cond_flag;     //### if # of events < anabreak, false
ofstream calibfile;     //### New(2010/3/3) format for MPPC calibration
                        //### constant table.  



TH1F*  fH_HighAdc[NumMod][2][NumTFB][NumCh];
TH1F*  fH_LowAdc[NumMod][2][NumTFB][NumCh];
Int_t  cAnaTrg;
Int_t  fMinAdcwTdcCut[NumMod][2][NumTFB][NumCh];
Int_t  anacounter; //
Int_t  anabreak;   //number of event for analysis

int StartTime;
int EndTime;

bool  New(){
  TFile* fFile = new TFile();
  
}

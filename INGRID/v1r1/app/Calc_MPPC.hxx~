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
#include "PM_Ch_config.cxx" 
#include "INGRID_BadCh_mapping.cxx"
#include <unistd.h> // using getopt      

FileStat_t            fs;
ana_MPPC*             fana_MPPC;
INGRID_Ch_config*     fINGRID_Ch_config;
PM_Ch_config*            fPM_Ch_config; 
INGRID_BadCh_mapping* fINGRID_BadCh_mapping;
Int_t                 fMinAdcwTdcCut[NumMod][2][NumTFB][LayerNumCh];
Int_t                 cAnaTrg;
Int_t                 nCalib; //### In this code, calib. root file is made 
                              //### every ANABREAK events. this number 
                              //### shows current file number
char                  cCalibFile[300];


const static int   ANABREAK     = 300;
const static float ErrorMinGain =   5;
const static float ErrorMaxGain =  25;

  TH1F* fH_sHighAdc[NumMod][2][NumTFB][LayerNumCh];
  TH1F* fH_sLowAdc[NumMod][2][NumTFB][LayerNumCh];


TFile*   fTFile;
TTree*   tree;
TH1F*    fH_HighAdc[NumMod][2][NumTFB][LayerNumCh]; //### High ADC dist. 
TH1F*    fH_LowAdc [NumMod][2][NumTFB][LayerNumCh]; //### Low  ADC dist.
Int_t    IngRun;       //### INGRID run number
Int_t    IngSubRun;    //### INGRID sub run number
Int_t    StartEvt;     //### Start event number of this calib. file
Int_t    StartTime;    //### Start time         of this calib. file
Int_t    EndEvt;       //### End   event number of this calib. file
Int_t    EndTime;      //### End   time         of this calib. file
Int_t    AnaEvt;      //###
Int_t    AnaTime;      //###
Int_t    Version;      //### version of anaMPPC
Int_t    MakeTime;     //### time this calib. file is made
Int_t    TrgMode;      //### trigger mode

/*#####################Added by B.Quilain on 2013/03/05#####################*/
Float_t  HighADC [NumMod][2][NumTFB][LayerNumCh];
Float_t  LowADC [NumMod][2][NumTFB][LayerNumCh];
Int_t    NDChan [NumMod][2][NumTFB][LayerNumCh];

Float_t  HighGain [NumMod][2][NumTFB][LayerNumCh];
Float_t  HighPed  [NumMod][2][NumTFB][LayerNumCh];
Int_t    LowPed   [NumMod][2][NumTFB][LayerNumCh];
Float_t  MPPCNoise[NumMod][2][NumTFB][LayerNumCh];
Float_t  MPPCCandA[NumMod][2][NumTFB][LayerNumCh];
Float_t  TdcThre  [NumMod][2][NumTFB][LayerNumCh];
Bool_t   GoodCh   [NumMod][2][NumTFB][LayerNumCh];
Bool_t   GainCheck;
Bool_t   EvtCheck;


//____________________ Start of Book ________________________________________ 

void Book(string fname){
  cout << "make calib file " << fname << endl;
  fTFile    = new TFile  (fname.c_str(), 
			  "recreate");
  tree      = new TTree  ("calibtree",
			  "tree for calibration");
  tree     -> Branch     ("IngRun", 
			  &IngRun,
			  "IngRun/I");
  tree     -> Branch     ("IngSubRun", 
			  &IngSubRun,
			  "IngSubRun/I");
  tree     -> Branch     ("StartTime", 
			  &StartTime,
			  "StartTime/I");
  tree     -> Branch     ("StartTime", 
			  &StartEvt,
			  "StartEvt/I");
  tree     -> Branch     ("EndTime", 
			  &EndTime,
			  "EndTime/I");
  tree     -> Branch     ("EndEvt", 
			  &EndEvt,
			  "EndEvt/I");
  tree     -> Branch     ("AnaTime", 
			  &AnaTime,
			  "AnaTime/I");
  tree     -> Branch     ("AnaEvt", 
			  &AnaEvt,
			  "AnaEvt/I");
  tree     -> Branch     ("HighGain", 
			  HighGain,
			  Form("HighGain[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));
  tree     -> Branch     ("HighPed", 
			  HighPed,
			  Form("HighPed[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));

  /*#####################Added by B.Quilain on 2013/03/05#####################*/ 
 tree     -> Branch     ("HighADC",
                          HighADC,
                          Form("HighGain[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));
 tree     -> Branch     ("LowADC",
			 LowADC,
			 Form("LowGain[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));
 tree     -> Branch     ("NDChan",
                         NDChan,
                         Form("NDChan[%d][%d][%d][%d]/I",NumMod,2,NumTFB,LayerNumCh));


  tree     -> Branch     ("LowPed", 
			  LowPed,
			  Form("LowPed[%d][%d][%d][%d]/I",NumMod,2,NumTFB,LayerNumCh));

  tree     -> Branch     ("TdcThre", 
			  TdcThre,
			  Form("TdcThre[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));

  tree     -> Branch     ("MPPCNoise", 
			  MPPCNoise,
			  Form("MPPCNoise[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));

  tree     -> Branch     ("MPPCCandA", 
			  MPPCCandA,
			  Form("MPPCCandA[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));

  tree     -> Branch     ("GoodCh", 
			  GoodCh,
			  Form("GoodCh[%d][%d][%d][%d]/O",NumMod,2,NumTFB,LayerNumCh));
  tree     -> Branch     ("GainCheck", 
			  &GainCheck,
			  "GainCheck/O");
  tree     -> Branch     ("EvtCheck", 
			  &EvtCheck,
			  "EvtCheck/O");
  GainCheck = false;
  EvtCheck  = false;

  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
	for(int ich = 0; ich < LayerNumCh; ich++){
          HighADC        [imod][iview][itfb][ich] = -1.e-5;
          LowADC         [imod][iview][itfb][ich] = -1.e-5;
	  NDChan         [imod][iview][itfb][ich] = -1.e-5;
	  HighGain       [imod][iview][itfb][ich] = -1.e-5;
	  HighPed        [imod][iview][itfb][ich] = -1.e-5;
	  LowPed         [imod][iview][itfb][ich] = -1;
	  TdcThre        [imod][iview][itfb][ich] = -1.e-5;
	  fMinAdcwTdcCut [imod][iview][itfb][ich] = 1000; //because this is used for minimum value
	  GoodCh         [imod][iview][itfb][ich] = false;

	  fH_sHighAdc     [imod][iview][itfb][ich]
	    = new TH1F( Form("fH_sHighAdc%02d%d%02d%02d",imod,iview,itfb,ich),
			Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),

			//1024, 0, 1024);
			150, 100, 250);

	  fH_sLowAdc      [imod][iview][itfb][ich]
	    = new TH1F( Form("fH_sLowAdc%02d%d%02d%02d",imod,iview,itfb,ich),
			Form("Low Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),

			1024, 0, 1024);

	}//ich
      }//itfb
    }//iview
  }//imod
}
//___________________________________________________________________
//___________________________________________________________________
void MakeHist(){
  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
	for(int ich = 0; ich < LayerNumCh; ich++){

	  fH_HighAdc     [imod][iview][itfb][ich]
	    = new TH1F( Form("fH_HighAdc%02d%d%02d%02d",imod,iview,itfb,ich),
			Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),

			//1024, 0, 1024);
			150, 100, 250);

	  fH_LowAdc      [imod][iview][itfb][ich]
	    = new TH1F( Form("fH_LowAdc%02d%d%02d%02d",imod,iview,itfb,ich),
			Form("Low Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),

			1024, 0, 1024);
	}//ich
      }//itfb
    }//iview
  }//imod

}

//____________________ End of Book ___________________________________________ 
//____________________________________________________________________________

void fGainCheck(){
  GainCheck  = true;
  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
	for(int ich = 0; ich < LayerNumCh; ich++){
	  if( GoodCh[imod][iview][itfb][ich] ){
	    if( HighGain[imod][iview][itfb][ich] < ErrorMinGain ||
		HighGain[imod][iview][itfb][ich] > ErrorMaxGain 
		)
	      GainCheck = false;
	  }
	}
      }
    }
  }
}

//____________________ Start of Close ________________________________________ 
void Close(){
  cout << "close current calib. root file..." << endl;

  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
	for(int ich = 0; ich < LayerNumCh; ich++){
	  fH_sHighAdc[imod][iview][itfb][ich] -> Add(fH_HighAdc[imod][iview][itfb][ich] );
	  fH_sLowAdc[imod][iview][itfb][ich] -> Add(fH_LowAdc[imod][iview][itfb][ich] );

	}//ich
      }//itfb
    }//iview
  }//imod

  tree      -> Fill();
  tree      -> Write();
  fTFile    -> Write();
  fTFile    -> Close();
}
//____________________ End of Close __________________________________________ 
//____________________________________________________________________________

//____________________ Start of Clear ________________________________________ 
void HistClear(){
  cout << "clear the variables and histograms...." << endl;
  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
	for(int ich = 0; ich < LayerNumCh; ich++){
	  //fH_HighAdc[imod][iview][itfb][ich] -> Clear();
	  fH_HighAdc[imod][iview][itfb][ich] -> Reset();
	  fH_HighAdc[imod][iview][itfb][ich] -> SetName ( Form("fH_HighAdc%02d%d%02d%02d",imod,iview,itfb,ich) );
	  fH_HighAdc[imod][iview][itfb][ich] -> SetTitle( Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich) );



	  //fH_LowAdc [imod][iview][itfb][ich] -> Clear();
	  fH_LowAdc [imod][iview][itfb][ich] -> Reset();
	  fH_LowAdc [imod][iview][itfb][ich] -> SetName ( Form("fH_LowAdc%02d%d%02d%02d",imod,iview,itfb,ich) );
	  fH_LowAdc[imod][iview][itfb][ich] -> SetTitle( Form("Low Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich) );

	  //delete    fH_HighAdc[imod][iview][itfb][ich];
	  //delete    fH_LowAdc [imod][iview][itfb][ich];
	}//ich
      }//itfb
    }//iview
  }//imod
  //tree   -> Clear();
  //fTFile -> Clear();
}
//____________________ End of Book ___________________________________________ 
//____________________________________________________________________________

void StartCalib(ND::TND280RawEvent* re){

  ND::TRawDataHeader header = re->GetHeader();
  StartTime  = header.GetTimeStamp();
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    StartEvt                  = runinfo.GetSeqNumber();
  }
  AnaEvt     =    0;

  MakeTime   = time(NULL);
  Version    = fana_MPPC -> get_version();

}


void EndCalib(ND::TND280RawEvent* re){
  ND::TRawDataHeader header = re->GetHeader();
  EndTime  = header.GetTimeStamp();
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    EndEvt                    = runinfo.GetSeqNumber();
  }

}

int evtnum(ND::TND280RawEvent* re){
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    int num = runinfo.GetSeqNumber();
    return num;
  }

}
int trgid(ND::TND280RawEvent* re){

  ND::THandle<ND::TTriggerBank> triggerBank;
  while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
    ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
    return (trigger.GetTriggerWord()>>48)&0xffff;
  }


}

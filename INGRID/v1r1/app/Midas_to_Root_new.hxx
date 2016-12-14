#ifndef __MIDAS_TO_ROOT_NEW_HXX__
#define __MIDAS_TO_ROOT_NEW_HXX__

// ND280 software includes
#include "TMidasBank.hxx"
#include "TMidasFile.hxx"
#include "TMidasBankProxy.hxx"
#include "TMidasBankProxyRegistry.hxx"
#include "TND280RawEvent.hxx"
#include "TRawDataHeader.hxx"
// oaRawEvent includes
#include "TTripTDigitBank.hxx"
#include "TRunInfoBank.hxx"
#include "TMidasTripTDigitItr.hxx"
#include "TMidasTripTDigit.hxx"
#include "TMCMBank.hxx"
#include "TTriggerBank.hxx"
// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TTree.h"
//C++ libraly includes
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>
#include <algorithm>
#include <deque>
#include <vector>
#include <sys/stat.h>
#include <unistd.h> // using getopt      
using namespace std;
// INGRID includes
#include "setup.hxx"
#include "INGRID_Ch_config.cxx"
#include "INGRID_BadCh_mapping.cxx"
#include "INGRID_Dimension.cxx"

const static Int_t   cNumMod    = 14;
const static Int_t   cNumTPL    = 11;
const static Int_t   cNumVeto   =  4;
const static Int_t   cNumCh     = 24;
const static Int_t   cNumLyr    =  2;
const static Int_t   cNumCyc    = 23;
const static Int_t   cTrgBeam   =  1;
const static Int_t   cTrgCalib  =  2;
const static Int_t   cTrgCosmic =128;


//__________________________________________________________
Int_t                 NumEvt;
Int_t                 TrgId;
Int_t                 cAnaEvt;
Int_t                 cAnaTrg;
Int_t                 cAnaMod;
Int_t                 cAnaCyc;

class Hit{
public:
  Int_t      pln;
  Int_t     view;
  Int_t       ch;
  Int_t      adc;
  Float_t     pe;   //p.e.
  Long_t     tdc;  
  Long_t    time;   //nsec
  Long_t      t0;//
  Long_t   rawtdc;//
};
Hit fhit;
vector<Hit> hit[cNumMod][cNumCyc];
bool withtdc(const Hit& left, const Hit& right){
  return left.tdc < right.tdc;
};
bool withview(const Hit& left, const Hit& right){
  return left.view < right.view;
};

bool withpln(const Hit& left, const Hit& right){
  return left.pln < right.pln;
};
bool withch(const Hit& left, const Hit& right){
  return left.ch < right.ch;
};

bool fSortTdc(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withtdc);
};
bool fSortPln(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withpln);
};

const static Double_t cTimeRsr = 30; //nsec

void Print(vector<Hit> hit[][cNumCyc], Int_t nummod, Int_t numcyc){
  for(int i=0;i<hit[nummod][numcyc].size();i++){
    cout<<"pln:"   <<hit[nummod][numcyc][i].pln
	<<"\tview:"<<hit[nummod][numcyc][i].view
	<<"\tch:"  <<hit[nummod][numcyc][i].ch
	<<"\tadc:" <<hit[nummod][numcyc][i].adc
	<<"\tpe:"  <<hit[nummod][numcyc][i].pe
	<<"\ttdc:" <<hit[nummod][numcyc][i].tdc
	<<"\ttime:" <<hit[nummod][numcyc][i].time
	<<endl;

  }
  cin.get();
}
void Print(vector<Hit> hit){
  for(int i=0;i<hit.size();i++){
    cout<<"pln:"   <<hit[i].pln
	<<"\t view:"<<hit[i].view
	<<"\t ch:"  <<hit[i].ch
	<<"\t adc:" <<hit[i].adc
	<<"\t pe:"  <<hit[i].pe
	<<"\t tdc:" <<hit[i].tdc
	<<"\t time:" <<hit[i].time
	<<endl;
  }

}

Int_t cTdcRsr = 15; // 30/2


bool fFindTimeClster(vector<Hit> hit[][cNumCyc], vector<Hit> &hitclster, Int_t nummod, Int_t numcyc, Long_t &time){
  if(hit[nummod][numcyc].size()>=2){
    Int_t first = (Int_t)hit[nummod][numcyc][0].time;
    Int_t end = (Int_t)hit[nummod][numcyc][hit[nummod][numcyc].size()-1].time;
    Int_t nbin = (end-first)/5;
    if(nbin<1)nbin=1;
    TH1F *h = new TH1F("h","h",nbin, first,end);
    for(int i=0;i<hit[nummod][numcyc].size();i++){
      h->Fill(hit[nummod][numcyc][i].time);
    }
    Int_t max=h->GetMaximumBin();
    //    TH1F *h2 = new TH1F("h2","h2",nbin, first,end);
    time = max*(end-first)/nbin+first;// - 0.5*(end-first)/nbin;


    Int_t count=0;
    vector<Hit>::iterator it;
    for(it = hit[nummod][numcyc].begin();it != hit[nummod][numcyc].end();it++){
      if(fabs(it->time-time)<cTdcRsr){
	count++;
	//	h2->Fill(it->time);
	hitclster.push_back(*it);
	it = hit[nummod][numcyc].erase(it);
	it--;
      }
    }

    if(count>=4){
      //      TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
      //      h->Draw();
      //	    c1->Update();
      //	    cin.get();
      //      h2->SetLineColor(2);
      //      h2->Draw("same");
      //     c1->Update();
      //    cin.get();
      h->Reset();
      delete h;
      //    h2->Reset();
      //    delete h2;
      return true;

    }
    h->Reset();
    delete h;
    //h2->Reset();
    //delete h2;
    return false;
  }
  else return false;
}

//__________________________________________________________
const static Double_t     cTdcBin  = 2.5;
const static Int_t        cTrgBin  =  10;
Int_t                     fTrgTime ;
Double_t
   pedestal[cNumMod][cNumLyr][cNumTPL+cNumVeto][cNumCh];
Double_t
       gain[cNumMod][cNumLyr][cNumTPL+cNumVeto][cNumCh];

void   read_calib(int crun){
  for(int i=0; i<14; i++){
    for(int j=0; j<2; j++){
      for(int k=0; k<14; k++){
	for(int l=0; l<24; l++){
	  pedestal[i][j][k][l]=100;
	  gain[i][j][k][l]=100;
	}
      }
    }
  }
  char buff[300];
  sprintf(buff,"/home/daq/data/MPPC_calib/ingrid_%08d_0000.txt",crun);
  ifstream file(buff);
  int mod, view, pln, ch;
  while(file>>mod>>view>>pln>>ch){
    file>>pedestal[mod][view][pln][ch]>>gain[mod][view][pln][ch];
  }
  file.close();
}

void   adc_calib(vector<Hit> hit[][cNumCyc], Int_t nummod, Int_t numcyc){
  for(Int_t i = 0 ; i < hit[nummod][numcyc].size() ; i++ ){
    hit[nummod][numcyc][i].pe = 1.0 * 
            ( hit[nummod][numcyc][i].adc - pedestal [ nummod ]
	                                            [ hit[nummod][numcyc][i].view ] 
                                                    [ hit[nummod][numcyc][i].pln]
                                                    [ hit[nummod][numcyc][i].ch ]  
	     )/ gain   [ nummod ]
                       [ hit[nummod][numcyc][i].view] 
                       [ hit[nummod][numcyc][i].pln ]
                       [ hit[nummod][numcyc][i].ch  ];
    /*
    if(nummod==0&&hit[nummod][numcyc][i].view==0&&hit[nummod][numcyc][i].pln==0&&hit[nummod][numcyc][i].ch==0){
      cout<<hit[nummod][numcyc][i].adc<<" "<<pedestal[nummod][hit[nummod][numcyc][i].view][hit[nummod][numcyc][i].pln][hit[nummod][numcyc][i].ch]<<" "<<gain[nummod][hit[nummod][numcyc][i].view][hit[nummod][numcyc][i].pln][hit[nummod][numcyc][i].ch]<<" "<<hit[nummod][numcyc][i].pe<<endl;
      cin.get();
    }
    */
  }                    
};
void   adc_calibcls(vector<Hit> &hit, Int_t nummod, Int_t numcyc){
  for(Int_t i = 0 ; i < hit.size() ; i++ ){
    hit[i].pe = 1.0 *
            ( hit[i].adc - pedestal [ nummod ]
	                            [ hit[i].view ] 
	                            [ hit[i].pln]
	                            [ hit[i].ch ]  
	     )/ gain   [ nummod ]
                       [ hit[i].view] 
                       [ hit[i].pln ]
                       [ hit[i].ch  ];

  }                    
};
void   adc_calibclsd(vector<Hit> &hit, Int_t nummod, Int_t numcyc){
  for(Int_t i = 0 ; i < hit.size() ; i++ ){
    hit[i].pe = 1.0 *
            ( hit[i].adc - pedestal [ nummod ]
	                            [ hit[i].view ] 
	                            [ hit[i].pln]
	                            [ hit[i].ch ]  
	     )/ gain   [ nummod ]
                       [ hit[i].view] 
                       [ hit[i].pln ]
                       [ hit[i].ch  ];


    cout<<nummod<<" "<<hit[i].view<<" "<<hit[i].pln<<" "<<hit[i].ch;
    cout<<" "<<hit[i].adc<<" "<<pedestal [ nummod ]
                                    [ hit[i].view ] 
                                    [ hit[i].pln]
                                    [ hit[i].ch ]
	<<" "<<gain   [ nummod ]
                      [ hit[i].view] 
                      [ hit[i].pln ]
                      [ hit[i].ch ]
        <<" "<<hit[i].pe<<endl;
    cin.get();

  }                    
};

void   tdc_calib(vector<Hit> hit[][cNumCyc], Int_t nummod, Int_t numcyc){
  for(int i=0;i<hit[nummod][numcyc].size();i++){
    
    hit[nummod][numcyc][i].time = cTdcBin * hit[nummod][numcyc][i].tdc 
                                - cTrgBin * fTrgTime;

    /*
    cout<<"tdc :"       <<cTdcBin * hit[nummod][numcyc][i].rawtdc
	<<"\tt0 :"       <<cTdcBin * hit[nummod][numcyc][i].t0
	<<"\t time :"    <<hit[nummod][numcyc][i].time	
	<<"\t trgtime :" <<fTrgTime*cTrgBin
	<<endl;
    */
  }
  
};
//__________________________________________________________


TFile* fTFile;
TTree* tree;
TH1F*  fHPeMean[cNumMod];
TH1F*  fHActMap[cNumMod];
TH2F*  fHChActMap[cNumMod][cNumLyr];
Int_t  fActPln;
Int_t  fMissPln;
Int_t  Cycle;
Float_t fPeMean;
Int_t  fMod;
Long_t fTime;
Int_t  nHit;
Int_t  nSpill;
Long_t  UTime;
vector<int>       adc;
vector<double>     pe;
vector<int>      nsec;
vector<long>       tdc;
vector<int>       pln;
vector<int>        ch;      
vector<int>      view;      
//### for study MCM time problem ######
long MCMTime_a_now;
long MCMTime_a_before;
long MCMTime_b_now;
long MCMTime_b_before;
long MCMTime_c;
long MCMTime_diff;
const long constant = 100000000;
void Print(){
  cout<<"a now:"        <<MCMTime_a_now
      <<"\ta before:"   <<MCMTime_a_before
      <<"\tb now:"      <<MCMTime_b_now
      <<"\tb before:"   <<MCMTime_b_before
      <<"\tnow(a+b)-before(a+b)"<< ( MCMTime_a_now + MCMTime_b_now ) - (MCMTime_a_before + MCMTime_b_before )
      <<"\tMCMTime c :"         <<MCMTime_c
      <<endl;
}

Int_t  fActNum[11]; //threshold
void Book(Int_t run){
  char buff[300];
  sprintf(buff,"/home/daq/data/root_new/ingrid_%08d_0000.root",run);
  fTFile  = new TFile(buff,"recreate");
  tree    = new TTree("tree","tree");
  tree->Branch("nHit"   ,&nHit, "nHit/I");
  tree->Branch("adc"    ,&adc           );
  tree->Branch("pe"     ,&pe            );
  tree->Branch("nsec"   ,&nsec          );
  tree->Branch("tdc"    ,&tdc           );
  tree->Branch("view"   ,&view          );
  tree->Branch("pln"    ,&pln           );
  tree->Branch("ch"     ,&ch            );
  tree->Branch("nSpill" ,&nSpill, "nSpill/I");
  tree->Branch("UTime"  ,&UTime , "UTime/L" );
  tree->Branch("TrdId"  ,&TrgId , "TrgId/I" );
  //### for study MCM time problem ######
  tree   ->Branch("MCMTime_a_now",    &MCMTime_a_now,     "MCMTime_a_now/L");
  tree   ->Branch("MCMTime_b_now",    &MCMTime_b_now,     "MCMTime_b_now/L");
  tree   ->Branch("MCMTime_a_before", &MCMTime_a_before,  "MCMTime_a_before/L");
  tree   ->Branch("MCMTime_b_before", &MCMTime_b_before,  "MCMTime_b_before/L");
  tree   ->Branch("MCMTime_diff",     &MCMTime_diff,      "MCMTime_diff/L");
  tree   ->Branch("MCMTime_c",        &MCMTime_c,         "MCMTime_c/L");

  char buff1[300],buff2[300];
  tree->Branch("NumEvt",&NumEvt,"NumEvt/I");
  tree->Branch("Cycle" ,&Cycle, "Cycle/I");
  tree->Branch("fMod"  ,&fMod,  "fMod/I");

}
void Save(){
  tree->AutoSave();
}
void Write(){
  tree->Write();
  fTFile->Write();
  fTFile->Close();
}

#endif

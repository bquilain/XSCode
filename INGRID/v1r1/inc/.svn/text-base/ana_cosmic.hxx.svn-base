#ifndef _ANA_COSMIC_H
#define _ANA_COSMIC_H

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
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "TApplication.h"
#include "setup.hxx"


class ana_cosmic{
private:
public:
  Long_t event_number;
  Int_t NSeq;
  Int_t NCyc;
  Double_t pedestal[NumMod][NumTFB][NumCh];
  Double_t gain[NumMod][NumTFB][NumCh];
  Bool_t tpl_trigger_flag[NumTFB];
  Int_t trigger_type;
  Bool_t pattern[NumTFB];
  Long_t UTime;
  Int_t Adc[NumMod][NumTFB][NumCh][NumCyc];
  Long_t Tdc[NumMod][NumTFB][NumCh][NumCyc];
  Double_t pe[NumMod][NumTFB][NumCh][NumCyc];
  Bool_t Hit[NumMod][NumTFB][NumCh][NumCyc];
  Bool_t Tdc_Hit[NumMod][NumTFB][NumCh][NumCyc];
  Bool_t Sync_Hit[NumMod][NumTFB][NumCh][NumCyc];
  void init_Sync_Hit(int numcyc);
  Bool_t Act_Hit[NumMod][NumTFB][NumCh][NumCyc];
  Int_t Expect_Hit_Channel[NumMod][NumTFB][NumCyc][2];
  Double_t cos_x[NumCyc],cos_y[NumCyc],path_length[NumCyc];

  ana_cosmic();
  void set_pedestal_and_gain(Int_t file_number);
  void initialize_tript_data();

  void ana_trigger_tpl_ver1(ULong64_t twl);  
  Int_t get_trigger_type(){return trigger_type;};
  void ana_event_number(ULong64_t twl);  
  Long_t get_event_number(){return event_number;};
  void ana_trigger_type(ULong64_t twl);
  Bool_t pattern_match(Int_t alg);
  void print_pattern();
  Bool_t cosmic_or_not(){
    if(this->trigger_type==128)return true;
    return false;
  }
  void convert_adc_to_pe();
  void cut_with_pe(Double_t cut);
  void tdc_synchronize(Long_t cut);//090523 version
  void tdc_synchronize_old(Long_t cut);//090523 version
  void cut_with_tdc();
  void only_one_activity();//only one activity Sync_Hit
  Int_t number_of_activity_layer(Int_t numcyc);
  Int_t number_of_activity_plane(Int_t numcyc);
  void fit_track(Int_t numcyc);



  ///////////////////////
  void debug();
  ///////////////////////
};

#endif

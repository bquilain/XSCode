#ifndef _ANA_MPPC_H
#define _ANA_MPPC_H

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
#include "setup.hxx"

#include "TApplication.h"
#include "setup.hxx"
const Int_t MAXONEPE=250;
const Int_t MINONEPE=100;
//const Int_t MAXPEDESTAL=200;
const Int_t MAXPEDESTAL=4000;
//const Int_t MINPEDESTAL=80;
const Int_t MINPEDESTAL=3000;
//const Int_t FITRANGE=3; 
const Int_t FITRANGEp=8;  
const Int_t FITRANGEm=10;  
const Int_t FITRANGE=8;  
//const Int_t HIST_MIN=2400;
const Int_t HIST_MIN=2550;
//const Int_t HIST_MAX=2800;
const Int_t HIST_MAX=2650;
//const Int_t EXPECT_GAIN=9;
const Int_t EXPECT_GAIN=30;
const Int_t EXPECT_PEDESTAL=160;

class ana_MPPC{
private:
  Double_t pedestal_peak_pos,pedestal_peak_sigma,pedestal_peak_height;
  Double_t onepe_peak_pos,onepe_peak_sigma,onepe_peak_height;
  Double_t gain;
  Double_t noise;
  Double_t crosstalk_and_afterpulse;
  //Double_t number_of_pedestal_events;
  double number_of_pedestal_events;
  //Double_t number_of_onepe_events;
  double  number_of_onepe_events;
  Double_t mean_pe;

  Int_t fBinMax,fBinEdge;
  Double_t fBinWidth,fXMax;
  TH1F *noise_hist;
  TH1F *noise_hist_wo_pedestal; 
  Int_t fMinX_noise_hist_wo_pedestal;
  Double_t fChisquare,fNDF;
  Double_t without_sigma;

public:
  ana_MPPC();

  Bool_t analysis_old_version(TH1F *noise_hist);//made by Masashi Otani 2009/05/22
  //Bool_t analysis_old_version(TH1F *noise_hist);//made by Masashi Otani 2009/05/22
  Bool_t analysis(TH1F *noise_hist);//refine version of analysis made at 2009/04/20

  void set_without_sigma(Double_t a){without_sigma=a;}
  void analysis_no_mppc(TH1F *pedestal_hist);

  Bool_t set_gain_and_pedestal(Int_t temp,Int_t MPPC_HV,Int_t nmod,Int_t ntfb,Int_t nch);//set MPPC gain and pedestal when temperature = temp, HV = MPPC_HV. if the gain and pedestal file doesn't exist, return false:: if HV=70.9,MPPC_HV=709,tempereture=13.2,temp=13
  Bool_t set_gain_and_pedestal(Int_t run_number,Int_t nmod,Int_t ntfb,Int_t nch);

  Double_t get_gain();
  //Double_t get_num_pedestal(){return number_of_pedestal_events;};
  double get_num_pedestal(){return number_of_pedestal_events;};
  //Double_t get_num_onepe(){return number_of_onepe_events;};
  double get_num_onepe(){return number_of_onepe_events;};
  Double_t get_pedestal();
  Double_t get_pedestal_sigma();
  Double_t get_pedestal_height();
  Double_t get_onepe();
  Double_t get_noise(Int_t entry);//entry including number of cycles
  Double_t get_crosstalk_and_afterpulse(Int_t entry);//entry including number of cycles
  Double_t get_mean_pe();
  Double_t get_chi2(){return fChisquare;}
  Double_t get_ndf(){return fNDF;}
};

#endif

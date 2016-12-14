#ifndef _ANA_BEAM_H
#define _ANA_BEAM_H

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
#include "setup.hxx"

#include "TApplication.h"

class ana_beam{
private:
  Int_t nmod,nch,ncyc,ntfb;
  Double_t pedestal[NumMod][NumTFB][NumCh];
  Double_t pedestal_sigma[NumMod][NumTFB][NumCh];
  Double_t gain[NumMod][NumTFB][NumCh];
  Double_t pe[NumMod][NumTFB][NumCh][NumCyc];
  Long_t Tdc[NumMod][NumTFB][NumCh][NumCyc];
public:
  Long_t Tdc_x[NumMod][NumTFB][NumCyc];//active plane's TDC
  Double_t Pe_x[NumMod][NumTFB][NumCyc];//active plane's Pe
  Long_t Tdc_y[NumMod][NumTFB][NumCyc];//active plane's TDC
  Double_t Pe_y[NumMod][NumTFB][NumCyc];//active plane's Pe

  ana_beam();
  Bool_t set_pedestal_and_gain(int file_number);
  void adc_to_pe(Int_t *Adc,Int_t numcyc);
  void set_Tdc(Long_t *tdc,Int_t numcyc);
  void set_pe(Double_t *pe,Int_t numcyc);

  Double_t get_pe(Int_t nummod,Int_t numtfb,Int_t numch,Int_t numcyc){return pe[nummod][numtfb][numch][numcyc];};
  

  Bool_t ana_plane_activity(int numpl,int numcyc,double PE_CUT,Long_t TDC_CUT);
  Bool_t top_veto_activity(int numcyc,double PE_CUT);//made by masashi otani at 090518 top_veto=TPL[11]Ch[0]~Ch[21]
  Bool_t bottom_veto_activity(int numcyc,double PE_CUT);//made by masashi otani at 090518 bottom_veto=TPL[12]Ch[0]~Ch[21]
  Bool_t left_veto_activity(int numcyc,double PE_CUT);//made by masashi otani at 090518 left_veto=TPL[12]Ch[23]~Ch[45]
 

  //static const Double_t CUT_PE=2.5;
  //static const Double_t DIFF_TDC=10;
  Int_t hit_chx,hit_chy;
  Int_t get_hit_chx(){return hit_chx;};
  Int_t get_hit_chy(){return hit_chy;};

};
#endif

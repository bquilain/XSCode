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

class hit{
private:
  Int_t hit_flag;
  Int_t hit_flag_x,hit_flag_y;
  Int_t hit_layer[NumTFB][NumLayer];

public:
  hit();
  Int_t one_hit_every_one_layer(Double_t *pe,Int_t numcyc,Double_t cut,Int_t *hit_channel_x,Int_t *hit_channel_y);//pe[nummod][numtfb][numch][numcyc]
  Int_t one_hit_every_one_layer_x(Double_t *pe,Int_t numcyc,Double_t cut,Int_t *hit_channel_x);//pe[nummod][numtfb][numch][numcyc]
  Int_t one_hit_every_one_layer_y(Double_t *pe,Int_t numcyc,Double_t cut,Int_t *hit_channel_y);//pe[nummod][numtfb][numch][numcyc]
  Int_t three_hit_AND_side_layer(Double_t *pe,Int_t nmod,Int_t ntfb,Int_t nch,Int_t ncyc,Double_t cut);
  Int_t exist_Tdc(Long_t *Tdc,Int_t nmod,Int_t ntfb,Int_t nch,Int_t ncyc);

};

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

class hit_ver2{
private:
  Bool_t hit_channel[NumMod][NumTFB][LayerNumCh][2];
  Bool_t hit_layer[NumMod][NumTFB][2];//0:x layer 1:y layer
  Int_t hit_channel_id[NumMod][NumTFB][LayerNumCh][2];//hit channel at each tfb
  Int_t hit_channel_number[NumMod][NumTFB][2];//number of hit channels at each tfb
  Int_t hit_layer_id[NumMod][NumTFB][2];//hit channel at each tfb
  Int_t hit_layer_number[NumMod][2];//number of hit channels at each tfb


public:
  hit_ver2(Double_t *pe,Double_t cut,Int_t ncyc);

  Bool_t get_hit_channel(Int_t nmod,Int_t ntfb,Int_t nch,Int_t nlayer);

};

///////////////////////////////////////////////////////////
////this is the class to monitor
//// gain
//// pedestal
//// sigma of pedestal
////at offline
////made by Masashi Otani
////last update 2009/04/17
/////////////////////////////////////////////////////////// 

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

static const double MAX_RATIO=0.3;
static const double MIN_RATIO=-0.3;
static const double ERROR_RATIO=0.2;
class plot{
public:
  Bool_t flag_gain_error[NumMod][NumTFB][UseNumCh];
private:
  int file_number_1,file_number_2;
  char buff1[300];
  double temp;

  TCanvas *c1;
  Double_t pe[NumMod][NumTFB][NumCh][NumCyc];
  Double_t pedestal_1[NumMod][NumTFB][UseNumCh],pedestal_2[NumMod][NumTFB][UseNumCh];
  Double_t pedestal_sigma_1[NumMod][NumTFB][UseNumCh],pedestal_sigma_2[NumMod][NumTFB][UseNumCh],pedestal_sigma_ratio[NumMod][NumTFB][UseNumCh];
  Double_t gain_1[NumMod][NumTFB][UseNumCh],gain_2[NumMod][NumTFB][UseNumCh],gain_ratio[NumMod][NumTFB][UseNumCh];
  Double_t noise_1[NumMod][NumTFB][UseNumCh],noise_2[NumMod][NumTFB][UseNumCh],noise_ratio[NumMod][NumTFB][UseNumCh];
  Double_t Ch[UseNumCh];
  Double_t Mod[NumMod];
  Double_t TFB[NumTFB];
  
  TGraph *graph_gain_1[NumMod][NumTFB],*graph_gain_2[NumMod][NumTFB],*graph_gain_ratio[NumMod][NumTFB];
  TGraph *graph_noise_1[NumMod][NumTFB],*graph_noise_2[NumMod][NumTFB],*graph_noise_ratio[NumMod][NumTFB];
  TGraph *graph_pedestal_sigma_1[NumMod][NumTFB],*graph_pedestal_sigma_2[NumMod][NumTFB],*graph_pedestal_sigma_ratio[NumMod][NumTFB];
  
  TH1F *hist_Adc[NumMod][NumTFB][UseNumCh];
  TH1F *hist_cyc[NumCyc];
  TH1F *h_gain[NumMod][NumTFB],*h_gain_1[NumMod][NumTFB],*h_gain_2[NumMod][NumTFB],*h_gain_all;
  
  TH1F *h_pedestal_sigma[NumMod][NumTFB];
  TGraph *graph_pedestal_1[NumMod][NumTFB],*graph_pedestal_2[NumMod][NumTFB];

public:
  bool flag;
  plot(int ffile_number_1,int ffile_number_2);//read pedestal,gain,noise


  void plot_gain();//plot each channel`s gain per each planes
  void plot_gain_ratio();//plot each channel`s gain ratio per each planes
  void plot_pedestal();//plot each channel`s pedestal per each planes
  void plot_pedestal_sigma();//plot each channel`s pedestal sigma per each planes


  void plot_noise();
  void plot_noise_ratio();
  
  void hist_gain();
  void hist_gain_all();
  void hist_pedestal_sigma();

  void plot_adc();

  void evdisp(Double_t *npe);//easy event display


  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  ////put 0 to no mppc channel data of veto planes /////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  void veto_correction(int nummod,int numtfb){
    if(numtfb==11){
      for(int numch=22;numch<UseNumCh;numch++){
	pedestal_1[nummod][numtfb][numch]=0;gain_1[nummod][numtfb][numch]=0;noise_1[nummod][numtfb][numch]=0;
	pedestal_2[nummod][numtfb][numch]=0;gain_2[nummod][numtfb][numch]=0;noise_2[nummod][numtfb][numch]=0;
	pedestal_sigma_1[nummod][numtfb][numch]=0;gain_ratio[nummod][numtfb][numch]=0;noise_1[nummod][numtfb][numch]=0;
	pedestal_sigma_2[nummod][numtfb][numch]=0;noise_ratio[nummod][numtfb][numch]=0;noise_2[nummod][numtfb][numch]=0;

      }
    }
    if(numtfb==12){
      pedestal_1[nummod][numtfb][22]=0;noise_1[nummod][numtfb][22]=0;gain_1[nummod][numtfb][22]=0;
      pedestal_1[nummod][numtfb][23]=0;noise_1[nummod][numtfb][23]=0;gain_1[nummod][numtfb][23]=0;
      pedestal_1[nummod][numtfb][46]=0;noise_1[nummod][numtfb][46]=0;gain_1[nummod][numtfb][46]=0;
      pedestal_1[nummod][numtfb][47]=0;noise_1[nummod][numtfb][47]=0;gain_1[nummod][numtfb][47]=0;
      pedestal_2[nummod][numtfb][22]=0;noise_2[nummod][numtfb][22]=0;gain_2[nummod][numtfb][22]=0;
      pedestal_2[nummod][numtfb][23]=0;noise_2[nummod][numtfb][23]=0;gain_2[nummod][numtfb][23]=0;
      pedestal_2[nummod][numtfb][46]=0;noise_2[nummod][numtfb][46]=0;gain_2[nummod][numtfb][46]=0;
      pedestal_2[nummod][numtfb][47]=0;noise_2[nummod][numtfb][47]=0;gain_2[nummod][numtfb][47]=0;
      pedestal_sigma_1[nummod][numtfb][22]=0;gain_ratio[nummod][numtfb][22]=0;gain_1[nummod][numtfb][22]=0;
      pedestal_sigma_1[nummod][numtfb][23]=0;gain_ratio[nummod][numtfb][23]=0;gain_1[nummod][numtfb][23]=0;
      pedestal_sigma_1[nummod][numtfb][46]=0;gain_ratio[nummod][numtfb][46]=0;gain_1[nummod][numtfb][46]=0;
      pedestal_sigma_1[nummod][numtfb][47]=0;gain_ratio[nummod][numtfb][47]=0;gain_1[nummod][numtfb][47]=0;
      pedestal_sigma_2[nummod][numtfb][22]=0;noise_ratio[nummod][numtfb][22]=0;gain_2[nummod][numtfb][22]=0;
      pedestal_sigma_2[nummod][numtfb][23]=0;noise_ratio[nummod][numtfb][23]=0;gain_2[nummod][numtfb][23]=0;
      pedestal_sigma_2[nummod][numtfb][46]=0;noise_ratio[nummod][numtfb][46]=0;gain_2[nummod][numtfb][46]=0;
      pedestal_sigma_2[nummod][numtfb][47]=0;noise_ratio[nummod][numtfb][47]=0;gain_2[nummod][numtfb][47]=0;
      }

  }
};

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
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "TApplication.h"

#include "hit_ver2.hxx"


hit_ver2::hit_ver2(Double_t *pe,Double_t cut,Int_t ncyc){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    hit_layer_number[nummod][0]=0;
    hit_layer_number[nummod][1]=0;
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      //X layer
      hit_channel_number[nummod][numtfb][0]=0;
      for(Int_t numch=0;numch<LayerNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+ncyc)>cut){
	  hit_channel[nummod][numtfb][numch][0]=true;
	  hit_channel_id[nummod][numtfb][hit_channel_number[nummod][numtfb][0]][0]=numch;
	  hit_channel_number[nummod][numtfb][0]++;
	}//if
	else{
	  hit_channel[nummod][numtfb][numch][0]=false;
	}
      }//numch
      if(hit_channel_number[nummod][numtfb][0]>0){
	hit_layer[nummod][numtfb][0]=true;
	hit_layer_id[nummod][hit_layer_number[nummod][0]][0]=numtfb;
	hit_layer_number[nummod][0]++;
      }//if
      //Y layer
      for(Int_t numch=0;numch<LayerNumCh;numch++){
	Int_t real_numch=numch+UseNumCh;
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+real_numch*NumCyc+ncyc)>cut){
	  hit_channel[nummod][numtfb][numch][1]=true;
	  hit_channel_id[nummod][numtfb][hit_channel_number[nummod][numtfb][1]][0]=numch;
	  hit_channel_number[nummod][numtfb][1]++;
	}//if
	else{
	  hit_channel[nummod][numtfb][numch][1]=false;
	}
      }//numch
      if(hit_channel_number[nummod][numtfb][1]>0){
	hit_layer[nummod][numtfb][0]=true;
	hit_layer_id[nummod][hit_layer_number[nummod][1]][1]=numtfb;
	hit_layer_number[nummod][1]++;
      }//if
    }//numtfb
  }//nummod
}

Bool_t hit_ver2::get_hit_channel(Int_t nmod,Int_t ntfb,Int_t nch,Int_t nlayer){
 
  if(hit_channel[nmod][ntfb][nch][nlayer])return true;
  else return false;

}

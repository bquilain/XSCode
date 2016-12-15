#ifndef __ANA_DRYRUN_HXX__
#define __ANA_DRYRUN_HXX__

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





//__________________________________________________________


class Hit{
public:
  Int_t       id;   //b/w Hit class and IngridHitSummary
  Int_t      pln;
  Int_t     view;
  Int_t       ch;
  Int_t      adc;
  Float_t     pe;   //p.e.
  Long_t     tdc;  
  Long_t    nsec; 
  Long_t    time;   //nsec
  Long_t    tnearhit;   //nsec
  double   posxy;
  double    posz;
  void clear(){
    id   =  -1;
    pln  =  -1;
    view =  -1;
    ch   =  -1;
    adc  =  -1;
    pe   =  -1.e-5;
    tdc  =  -1;
    nsec =  -1;
    time =  -1;
    tnearhit = -1;
    posxy=  -1.e-5;
    posz =  -1.e-5;
  }
};

vector<Hit> allhit;

bool withtime(const Hit& left, const Hit& right){
  return left.time < right.time;
};

bool fSortTime(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withtime);
};



void Print(vector<Hit> hit){
  for(int i=0;i<hit.size();i++){
    cout<<"pln:"     <<hit[i].pln
	<<"\t view:" <<hit[i].view
	<<"\t ch:"   <<hit[i].ch
	<<"\t adc:"  <<hit[i].adc
	<<"\t pe:"   <<hit[i].pe
	<<"\t tdc:"  <<hit[i].tdc
	<<"\t time:" <<hit[i].time
	<<"\t posxy:"<<hit[i].posxy
	<<"\t posz:" <<hit[i].posz
	<<endl;
  }

}


Int_t cTdcRsr = 50; //nsec
float cPeCut  = 18;
//float cPeCut  = 0;


bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, 
		     Long_t &ctime, int maxhit=3){
  //fSortTime( hit );
  Int_t nhit = hit.size();

  if(nhit <= maxhit)return false;
  for(Int_t i=0 ; i<nhit-maxhit; i++){
    //if( fabs( hit[i].time - hit[i+maxhit-1].time) < cTdcRsr ){
    if( fabs( hit[i].time - hit[i+maxhit].time) < cTdcRsr ){

      /*
      long basetime = 0;
      //for( int j=0; j<maxhit; j++){
      for( int j=0; j<maxhit+1; j++){
	basetime += hit[i+j].time;
      }
      //basetime = basetime/(maxhit);
      basetime = basetime/(maxhit+1);
      */
     
      //#### change the definition of basetime   ###
      //#### ~2010/4/29 mean time of first hits  ###
      //#### 2010/4/29~ time of highest hit      ###
      long basetime;
      float highpe = 0;
      float sumpe  = 0;
      for(int j=0; j<maxhit+1; j++){
	if( hit[i+j].pe > highpe ){
	  basetime = hit[i+j].time;
	  highpe   = hit[i+j].pe;
	  sumpe    = sumpe + hit[i+j].pe;
	}
      }
      if( sumpe < cPeCut )
	continue;

      vector<Hit>::iterator it;
      Int_t ncount=0;
      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cTdcRsr ){
	  ncount++;
	}
      }
      if(ncount<=maxhit)continue;

      //###### clstering ######
      //#######################
      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cTdcRsr ){
	  hitclster.push_back(*it);
	  it = hit.erase(it);
	  it--;
	}
      }

      //###### caluculate time of clster with max p.e. ######
      //#####################################################
      highpe = 0;
      for(int j=0; j<hitclster.size(); j++){
	if( highpe < hitclster[j].pe ){
	  highpe  = hitclster[j].pe;
	  ctime   = hitclster[j].time;
	}
      }

      return true;
    }
  }//i
  return false;
}

bool fGetTNearHit(vector<Hit> &hit){
  int nhit = hit.size();
  if( nhit == 0 ){
    return true;
  }
  if( nhit == 1 ){
    hit[0].tnearhit = 0;
    return true;
  }

  for(int ihit=0; ihit<nhit; ihit++){
    if( ihit == nhit-1 )
      hit[ihit].tnearhit = abs( hit[ ihit ].time - hit[ ihit-1 ].time );
    if( ihit == 0 )
      hit[ihit].tnearhit = abs( hit[ ihit + 1 ].time - hit[ ihit ].time );

    else{
      long bnearhit = abs(hit[ihit+1].time - hit[ihit].  time);
      long anearhit = abs(hit[ihit]  .time - hit[ihit-1].time);
      hit[ihit].tnearhit = min( bnearhit, anearhit );
    }
  }
}

bool fMergeClster(vector<Hit> &hit1, vector<Hit> &hit2, 
		  vector<Hit> &mergehit ){


      //###### clstering ######
      //#######################
      vector<Hit>::iterator it;
      for(it = hit1.begin() ; it != hit1.end(); it++){
	mergehit.push_back(*it);
	it = hit1.erase(it);
	it--;
      }

      for(it = hit2.begin() ; it != hit2.end(); it++){
	mergehit.push_back(*it);
	it = hit2.erase(it);
	it--;
      }
}

#endif

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
#include "setup.hxx"


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

class Hit{
public:
  Int_t      pln;
  Int_t     view;
  Int_t       ch;
  Int_t      adc;
  Float_t     pe;   //p.e.
  Long_t     tdc;  
  Long_t    nsec;  
  Long_t    time;   //nsec
};
vector<Hit> allhit;
bool withtdc(const Hit& left, const Hit& right){
  return left.tdc > right.tdc;
};
bool withtime(const Hit& left, const Hit& right){
  return left.time < right.time;
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
bool fSortTime(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withtime);
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
                    // changed 09/11/17 when Rsr=20,
                    // ncount <=3 to ncount <=4


bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time){
  //before using this function please fSortTdc
  Int_t nhit = hit.size();

  if(nhit < 2)return false;
  for(Int_t i=1 ; i<nhit-1; i++){
    if( (fabs( hit[i].time - hit[i-1].time) < cTdcRsr &&
	 fabs( hit[i].time - hit[i+1].time) < cTdcRsr )   ){
      Long_t basetime = hit[i].time;
      time = basetime;
      //degug 1
      /*
      Float_t start = hit[0].time;
      Float_t end   = hit[nhit-1].time;
      cout<<hit[i].time<<"\t"<<start<<"\t"<<end<<endl;
      Int_t range = end - start;
      TH1F *h1 = new TH1F("h1","h1",range/10, start, end);
      for(Int_t j=0 ; j<nhit; j++){
	h1->Fill(hit[j].time);
      }
      */


      vector<Hit>::iterator it;
      Int_t ncount=0;
      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cTdcRsr ){
	  ncount++;
	}
      }

      if(ncount<=3)continue;

      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cTdcRsr ){
	  hitclster.push_back(*it);
	  it = hit.erase(it);
	  it--;
	}
      }


      //degug 1
      /*
      TH1F *h2 = new TH1F("h2","h2",range/10, start, end);
      for(Int_t j=0 ; j<hitclster.size(); j++){
	h2->Fill(hitclster[j].time);
      }
      TCanvas *c1 = new TCanvas("c1","c1",10,10,500,500);
      h1->SetLineWidth(3);
      h1->Draw();
      h2->SetLineColor(2);
      h2->Draw("same");
      c1->Update();
      h1->Reset();
      h2->Reset();
      delete h1;
      delete h2;
      */

      return true;
    }
  }//i
  return false;
}

/*
bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time){

  if(hit.size()>=2){
    Int_t first = (Int_t)hit[0].time;
    Int_t end   = (Int_t)hit[hit.size()-1].time;
    Int_t nbin = (end-first)/5;
    if(nbin<1)nbin=1;
    TH1F *h = new TH1F("h","h",nbin, first,end);
    for(int i=0;i<hit.size();i++){
      h->Fill(hit[i].time);
    }
    Int_t max=h->GetMaximumBin();
         TH1F *h2 = new TH1F("h2","h2",nbin, first,end);
    time = max*(end-first)/nbin+first;// - 0.5*(end-first)/nbin;


    Int_t count=0;
    vector<Hit>::iterator it;

    for(it = hit.begin();it != hit.end();it++){
      if(fabs(it->time-time)<cTdcRsr){
	count++;
		h2->Fill(it->time);
	hitclster.push_back(*it);
	it = hit.erase(it);
	it--;
      }
    }

    if(count>=4){
            TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
            h->Draw();
      	    c1->Update();
      	    cin.get();
            h2->SetLineColor(2);
            h2->Draw("same");
           c1->Update();
          cin.get();
      h->Reset();
      delete h;
          h2->Reset();
          delete h2;
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
*/
//__________________________________________________________
const static Double_t     cTdcBin  = 2.5;
const static Int_t        cTrgBin  =  10;
Int_t                     fTrgTime ;
Double_t
   pedestal[cNumMod][cNumLyr][cNumTPL+cNumVeto][cNumCh];
Double_t
       gain[cNumMod][cNumLyr][cNumTPL+cNumVeto][cNumCh];

void   read_calib(int crun){
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
    
    hit[nummod][numcyc][i].time = cTdcBin * hit[nummod][numcyc][i].tdc - cTrgBin * fTrgTime;

  }
  
};
//__________________________________________________________



#endif

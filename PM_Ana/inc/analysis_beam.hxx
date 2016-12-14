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
#include "INGRID_Dimension.cxx"
INGRID_Dimension *fINGRID_Dimension;

const static Int_t   cNumMod    = 14;
const static Int_t   cNumTPL    = 11;
const static Int_t   cNumVeto   =  4;
const static Int_t   cNumCh     = 24;
const static Int_t   cNumLyr    =  2;
const static Int_t   cNumCyc    = 23;
const static Int_t   cTrgBeam   =  1;
const static Int_t   cTrgCalib  =  2;
const static Int_t   cTrgCosmic =128;
const static Int_t   NumPln     =11;
const static Int_t   NIron     =10;
const static Int_t   Nscinti     =24;
const static Int_t   XLyr     =0;
const static Int_t   YLyr     =1;


//__________________________________________________________
Int_t                 TrgId;
Int_t                 cAnaEvt;
Int_t                 cAnaTrg;
Int_t                 cAnaMod;

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
    posxy=  -1.e-5;
    posz =  -1.e-5;
  }
};
vector<Hit> allhit;
bool withtdc(const Hit& left, const Hit& right){
  return left.tdc > right.tdc;
};
bool withposz(const Hit& left, const Hit& right){
  return left.posz < right.posz;
};

bool withview(const Hit& left, const Hit& right){
  return left.view < right.view;
};
bool withtime(const Hit& left, const Hit& right){
  return left.time < right.time;
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
bool fSortPosz(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withposz);
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


//old version modified 09/11/25
/*
bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time){
  //before using this function please fSortTdc
  Int_t nhit = hit.size();

  if(nhit < 2)return false;
  for(Int_t i=1 ; i<nhit-1; i++){
    if( (fabs( hit[i].time - hit[i-1].time) < cTdcRsr &&
	 fabs( hit[i].time - hit[i+1].time) < cTdcRsr )   ){
      Long_t basetime = hit[i].time;
      time = basetime;
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
      return true;
    }
  }//i
  return false;
}
*/

Int_t cTdcRsr = 50; //nsec
//modified version 09/12/01
//when caluculating cluster time, weight by # of p.e.
bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time, Long_t &ctime, int &dimch, int maxhit=3){
  //before using this function please fSortTime
  Int_t nhit = hit.size();

  if(nhit <= maxhit)return false;
  for(Int_t i=0 ; i<nhit-maxhit; i++){
    if( fabs( hit[i].time - hit[i+maxhit-1].time) < cTdcRsr ){

      long basetime = 0;
      for( int j=0; j<maxhit; j++){
	basetime += hit[i+j].time;
      }
      basetime = basetime/maxhit;
      /*
      long basetime = ( hit[i].time + 
			hit[i+1].time +
			hit[i+2].time +
			hit[i+3].time) / 4;
      */
      vector<Hit>::iterator it;
      Int_t ncount=0;
      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cTdcRsr ){
	  ncount++;
	}
      }
      if(ncount<=maxhit)continue;


      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cTdcRsr ){
	  hitclster.push_back(*it);
	  it = hit.erase(it);
	  it--;
	}
      }
      //###### caluculate fTime ##################
      //##########################################
      time = 0 ;
      for(int j=0; j<hitclster.size(); j++){
	time += hitclster[j].time;
      }
      time = 1.0 * time / hitclster.size();

      //###### caluculate fTime ##################
      //###### with weighting p.e. ###############
      //##########################################
      /*
      double totalpe = 0;
      for(int i=0; i<hitclster.size(); i++){
	totalpe += hitclster[i].pe;
      }

      double ctimetemp = 0 ;
      for(int i=0; i<hitclster.size(); i++){
	ctimetemp += 1.0 * (hitclster[i].time) * hitclster[i].pe  ;
      }
      ctimetemp = 1.0 * ctimetemp / (  totalpe);
      ctime = ctimetemp;
      */

      //###### caluculate fTime ##################
      //###### with max p.e. #####################
      //##########################################
      double highpe = 0;
      Int_t  dimview = 0;
      Int_t  dimpln  = 0;
      for(int j=0; j<hitclster.size(); j++){
	if( highpe < hitclster[j].pe ){
	  highpe  = hitclster[j].pe;
	  ctime   = hitclster[j].time;
	  dimview = hitclster[j].view;
	  dimpln  = hitclster[j].pln;
	}
      }
      highpe = 0;
      dimch  = -1;
      for(int j=0; j<hitclster.size(); j++){
	if( hitclster[j].view != dimview &&
	    hitclster[j].pln  == dimpln  &&
	    hitclster[j].pe   >  highpe){
	  highpe = hitclster[j].pe;
	  dimch  = hitclster[j].ch;
	}
      }
      if( dimch != -1 && dimview == 1 ) dimch = 23 - dimch;

      return true;
    }
  }//i
  return false;
}
/*
bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time){
  //before using this function please fSortTdc
  Int_t nhit = hit.size();

  if(nhit < 3)return false;
  for(Int_t i=0 ; i<nhit-3; i++){
    if( fabs( hit[i].time - hit[i+3].time) < cTdcRsr ){

      long basetime = ( hit[i].time + 
			hit[i+1].time +
			hit[i+2].time +
			hit[i+3].time) / 4;

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
      //###### caluculate fTime ##################
      time = 0 ;
      for(int i=0; i<hitclster.size(); i++){
	time += hitclster[i].time;
      }
      time = 1.0 * time / hitclster.size();

      return true;
    }
  }//i
  return false;
}
*/
Int_t cCoTime = 10;
bool fFindTimeCoincidence(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time){
  //before using this function please fSortTdc
  Int_t nhit = hit.size();

  if(nhit < 3)return false;
  for(Int_t i=0 ; i<nhit-2; i++){
    if( fabs( hit[i].time - hit[i+2].time) <= cCoTime ){

      long basetime = ( hit[i].time + 
			hit[i+1].time + 
			hit[i+2].time )/3;
      vector<Hit>::iterator it;

      for(it = hit.begin() ; it != hit.end(); it++){
	if( fabs( basetime - it->time) < cCoTime ){
	  hitclster.push_back(*it);
	  it = hit.erase(it);
	  it--;
	}
      }
      return true;
    }
  }//i
  return false;
}


bool fFindTimeClsterUsePe(vector<Hit> &hit, vector<Hit> &hitclster, Long_t &time, double threshold){
  //before using this function please fSortTdc
  Int_t nhit = hit.size();

  //delete lower than threshold
  vector<Hit>::iterator it;
  for(it = hit.begin() ; it != hit.end(); it++){
    if( it->pe < threshold ){
      it = hit.erase(it);
      it--;
    }
  }

  nhit       = hit.size();  //redefine nhit
  if(nhit < 3){
    return false;
  }
  for(Int_t i=0 ; i<nhit-3; i++){
    if( fabs( hit[i].time - hit[i+3].time) < cTdcRsr ){

      long basetime = ( hit[i].time + 
			hit[i+1].time +
			hit[i+2].time +
			hit[i+3].time) / 4;

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
      hit.clear();
      return true;
    }
  }//i
  hit.clear();
  return false;
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


//######### for analysis cosmic ######################
//####################################################
TGraph* TGtrackx;
TGraph* TGtracky;
TF1*    pol1x;
TF1*    pol1y;
class pos{
public:
  double posxy;
  double  posz;
};

/*
bool withposz(const pos& left, const pos& right){
  return left.posz < right.posz;
};
bool fSortposz(vector<pos> &a){
  std::stable_sort(a.begin(), a.end(), withposz);
};
*/
/*
bool LinFit(vector<Hit> &hitclster, 
	    float *ax, float *bx, 
	    float *ay, float *by){
  vector<pos>     xv,yv;
  double tempxy, tempz;
  pos temppos;
  int n=0;
  for(int i=0; i<hitclster.size(); i++){
    if(hitclster[i].pe>4.5)n++;
  }
  if(n<=3)return false;
  for(Int_t i=0; i<hitclster.size(); i++){
    //##### X Layer ###########################
    //############# ###########################
    //if(hitclster[i].view==0){
    if(hitclster[i].pe>4.5&&hitclster[i].view==0){
      if(hitclster[i].pln<=10){//Trackin plane
	if(fINGRID_Dimension->get_posXY(0 ,
					XLyr ,
					hitclster[i].pln ,
					hitclster[i].ch ,
					&tempxy , 
					&tempz
					)
	   )
	  {
	    temppos.posxy = tempxy;
	    temppos.posz = tempz;
	    xv.push_back(temppos);
	  }
      }//TPL
      if(hitclster[i].pln==13||hitclster[i].pln==14){//VETO plane
	if(fINGRID_Dimension->get_posVeto(0 ,
					XLyr ,
					hitclster[i].pln ,
					hitclster[i].ch ,
					&tempxy , 
					&tempz
					)
	   )
	  {
	    temppos.posxy = tempxy;
	    temppos.posz = tempz;
	    xv.push_back(temppos);
	  }
      }//Veto
    }//X
    //##### Y Layer ###########################
    //############# ###########################
    //if(hitclster[i].view==1){
    if(hitclster[i].pe>4.5&&hitclster[i].view==1){
      if(hitclster[i].pln<=10){//Trackin plane
	if(fINGRID_Dimension->get_posXY(0 ,
					YLyr ,
					hitclster[i].pln ,
					hitclster[i].ch ,
					&tempxy , 
					&tempz
					)
	   )
	  {
	    temppos.posxy = tempxy;
	    temppos.posz = tempz;
	    yv.push_back(temppos);
	  }
      }//TPL
      if(hitclster[i].pln==11||hitclster[i].pln==12){//VETO plane
	if(fINGRID_Dimension->get_posVeto(0 ,
					YLyr ,
					hitclster[i].pln ,
					hitclster[i].ch ,
					&tempxy , 
					&tempz
					)
	   )
	  {
	    temppos.posxy = tempxy;
	    temppos.posz = tempz;
	    yv.push_back(temppos);
	  }
      }//Veto
    }//Y
  }//hit loop


  //#######  Fit X Layer #######################
  //############################################
  fSortposz(xv);
  Int_t    nhitx = xv.size();
  Double_t     x[nhitx];
  Double_t    zx[nhitx];
  for(int i = 0 ; i < nhitx ; i++){
    x [i] = xv[i].posxy;
    zx[i] = xv[i].posz;
  }

  TGtrackx  =  new TGraph(nhitx, zx, x);

  if(DRAW){
    fCMcanvas -> cd();
    fPMviewX  -> cd();
    TGtrackx  -> Draw("same");
    pol1x     -> Draw("same");
    fCMcanvas -> Update();
  }

  pol1x     =  new    TF1("pol1x","pol1",-10,150);
  TGtrackx  -> Fit("pol1x","qn","",-10,150);
  Double_t    chi2x   = pol1x -> GetChisquare();
  Double_t    ndfx    = pol1x -> GetNDF();
  Double_t    chindfx = chi2x/ndfx;

  //#######  Fit X Layer #######################
  //############################################
  fSortposz(yv);
  Int_t    nhity = yv.size();
  Double_t     y[nhity];
  Double_t    zy[nhity];
  for(int i = 0 ; i < nhity ; i++){
    y [i] = yv[i].posxy;
    zy[i] = yv[i].posz;
  }
  TGtracky  =  new TGraph(nhity, zy, y);

  if(DRAW){
    fCMcanvas -> cd();
    fPMviewY  -> cd();
    TGtracky  -> Draw("same");
    pol1y     -> Draw("same");
    fCMcanvas -> Update();
  }

  pol1y     =  new    TF1("pol1y","pol1",-10,150);
  TGtracky  -> Fit("pol1y","qn","",-10,150);
  Double_t    chi2y   = pol1y -> GetChisquare();
  Double_t    ndfy    = pol1y -> GetNDF();
  Double_t    chindfy = chi2y/ndfy;
  //####### Clear and delete function ###########


  if(chindfx < 20 && chindfy < 20 ){
    *ax = pol1x -> GetParameter(1);
    *bx = pol1x -> GetParameter(0);
    *ay = pol1y -> GetParameter(1);
    *by = pol1y -> GetParameter(0);

    pol1x    -> Clear();
    pol1y    -> Clear();
    TGtrackx -> Clear();
    TGtracky -> Clear();
    delete     pol1x;
    delete     pol1y;
    delete  TGtrackx;
    delete  TGtracky;
    
    return true;
  }
  else{
    pol1x    -> Clear();
    pol1y    -> Clear();
    TGtrackx -> Clear();
    TGtracky -> Clear();
    delete     pol1x;
    delete     pol1y;
    delete  TGtrackx;
    delete  TGtracky;

    return false;
  }
}
*/

#endif

#ifndef __LOLIRECON_HXX__
#define __LOLIRECON_HXX__

//#define PMEDGE

// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
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

#include "Lolirecon.hxx"
//#include "LoliAna.hxx"
#include "INGRID_Dimension.cxx"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"


const int Cview = 2;
const int Cpln = 200;
const int Cch = 200;
const int Ccls = 200;
const int Ccell = 5000;
const int Chit = 5000;
const int Ctra = 5000;
//const int Cdist = 4;
const int Cdist = 6;
const int Cmod = 10;
const int Cdir = 3;
const int MAX_CELLVALUE=80;


//__________________________________________________________


class Hit{
public:
  Int_t       id;   //b/w Hit class and IngridHitSummary
  Int_t      mod;
  Int_t      pln;
  Int_t     view;
  Int_t       ch;
  Int_t     used;
  //Int_t      adc;
  Float_t     pe;   //p.e.
  Float_t     pe_cross;   //p.e.
  Float_t   lope;   //p.e.
  Long_t     tdc;  
  //Long_t    nsec; 
  Long_t    time;   //nsec
  //Long_t    tnearhit;   //nsec
  double   posxy;
  double    posz;
  //Int_t      pdg;

  void clear(){
    id   =  -1;
    mod  =  -1;
    pln  =  -1;
    view =  -1;
    ch   =  -1;
    //adc  =  -1;
    pe   =  -1.e-5;
    lope =  -1.e-5;
    tdc  =  -1;
    //nsec =  -1;
    time =  -1;
    //tnearhit = -1;
    posxy=  -1.e-5;
    posz =  -1.e-5;
    //pdg  =  -1;
  }
};
vector<Hit> allhit;
vector<Hit> allhit_for_disp;
vector<Hit> hitcls;
vector<Hit> hitcls_for_joint;
IngridSimVertexSummary* simver_for_disp;
IngridSimParticleSummary* simpar_for_disp;

class Track{
public:
  Int_t     mod;
  Int_t     view;
  Int_t     ipln;
  Int_t     fpln;
  Float_t    ixy;
  Float_t    fxy;
  Float_t     iz;
  Float_t     fz;
  Float_t  slope;
  Float_t intcpt;
  Float_t    ang;
  Float_t  clstime;
  Bool_t    veto;
  Bool_t    edge;
  Bool_t    stop;
  Float_t vetodist;
  //Float_t totpe;
  //Int_t isohit;
  vector<Int_t> hitid;
  vector<Bool_t> isohit;
  //Int_t pdg[4];

  void clear(){
    mod    =  -1;
    view   =  -1;
    ipln   =  -1;
    fpln   =  -1;
    ixy    =  -1.e-5;
    fxy    =  -1.e-5;
    iz     =  -1.e-5;
    fz     =  -1.e-5;
    slope  =  -1.e-5;
    intcpt =  -1.e-5;
    ang    =  -1.e-5;
    clstime=  -1.e-5;
    veto   = false;
    edge   = false;
    stop   = false;
    vetodist   =  -1.e-5;
    hitid.clear();
    isohit.clear();
    //totpe   =  -1.e-5;
    //isohit   =  -1;
    //memset(pdg,-1,sizeof(pdg));
  }
};

vector<Track> alltrack;


int chmax(int mod, int view, int pln, int axis=0){
  int maxch=0;
  if(mod<15)maxch=24;
  else if(mod==15){
        INGRID_Dimension *fdim = new INGRID_Dimension();
        maxch = fdim -> get_chmax_loli( mod, view, pln, axis); //mod view pln axis
	delete fdim;
  }
  else if(mod==16){
	maxch=32;
  }
  return maxch;
}

int plnmax(int mod){
  int maxpln=0;
  if(mod<15)maxpln=11;
  else if(mod==15){
        maxpln = 60;
  }
  else if(mod==16) maxpln=18;
  return maxpln;
}

int plnmax(int mod, int view, int pln, int axis=0){
  int maxpln=0;
  if(mod<15)maxpln=11;
  else if(mod==15){
        INGRID_Dimension *fdim = new INGRID_Dimension();
        maxpln = fdim -> get_plnmax_loli( mod, view, pln, axis); //mod view pln axis
	delete fdim;
  }
  else if(mod==16) maxpln=18;
  return maxpln;
}


float zposi(int mod,int view,int pln){
  double posiz;
  if(mod==16){
    if(view==0){
      if(pln==0)posiz=5;
      else posiz=46*pln+9;
    }
    else{
      if(pln==0)posiz=28;
      else posiz=46*pln+32;
    }
  }
  else if(mod==15){
    return 0;
  }
  else{
    if(view==0){
      posiz=5+105*pln+1074.5;
    }
    else{
      posiz=5+105*pln+1074.5+10;
    }
  }
  return posiz;
}


float zposi(int mod,int view,int pln, int ch, int axis=0){
  double posiz, posixy;
  if(mod==16){
    if(view==0){
      if(pln==0)posiz=5;
      else posiz=46*pln+9;
    }
    else{
      if(pln==0)posiz=28;
      else posiz=46*pln+32;
    }
  }
  else if(mod==15){
        INGRID_Dimension *fdim = new INGRID_Dimension();
        fdim -> get_posi_lolirecon( mod, view, pln, ch, axis, &posixy, &posiz); //mod view pln ch axis in recon
	//posiz = 10. * posiz; //cm -> mm
	if(axis==0)posiz = posiz*10  + 409.5; //shift
	if(axis==1)posiz = posixy*10 + 600.; //shift
	delete fdim;
  }
  else{
    if(view==0){
      //posiz=5+105*pln+1074.5;
      posiz=5+105*pln+1074.5+10;
    }
    else{
      //posiz=5+105*pln+1074.5+10;
      posiz=5+105*pln+1074.5;
    }
  }  
  return posiz;
}

float xyposi(int mod, int view, int pln,int ch, int axis=0){
  double posixy, posiz;
  
  if(mod==16){
    if(pln==0){
      posixy=ch*50+25;
    }
    else{
      if(ch<8)posixy=ch*50+25;
      else if(ch<24)posixy=412.5+(ch-8)*25;
      else posixy=(ch-8)*50+25;
    }
  }
  else if(mod==15){
      INGRID_Dimension *fdim = new INGRID_Dimension();
      fdim -> get_posi_lolirecon( mod, view, pln, ch, axis, &posixy, &posiz); //mod view pln ch axis in recon
      //posixy = 10. * posixy; //cm -> mm
      if(axis==0)posixy = posixy*10 + 600.; //shift
      if(axis==1)posixy = posiz*10  + 409.5; //shift
      delete fdim;
  }
  else{
    posixy=ch*50+25;
  }
  return posixy;
}

float vposiz(int mod,int ch){
  float X;
  //if(mod!=16)
  if(mod!=15)
    X=5+9.5+50*ch+1075;
  else
    X=5+9.5+50*ch;
  return X;
}


int vposixy(int mod,int pln){
  int Y=0;
  if(mod==0){
    if(pln==11)Y=(int)-105.75;  //right
    else if(pln==12)Y=1309;//left
    else if(pln==13)Y=-59; //botom
    else Y=1284;           //top
  }
  else if(mod<7){
    if(pln==11)Y=-191;  //right
    else if(pln==12)Y=1309;//left
    else if(pln==13)Y=-59; //botom
    else Y=1284;           //top
  }
  else if(mod==7){
    if(pln==11)Y=(int)-105.75;  //right
    else if(pln==12)Y=1309;//left
    else if(pln==13)Y=-59; //botom
    else Y=1284;           //top
  }
  else if(mod<14){
    if(pln==11)Y=(int)-105.75;  //right
    else if(pln==12)Y=1309;//left
    else if(pln==13)Y=-216; //botom
    else Y=1284;           //top
  }


  else{//proton module
    if(pln==18||pln==21)Y=-55;
    else Y=1255;
  }
  return Y;
}

float sciwidth(int mod, int view, int pln, int ch, int axis, int flag=0){
  float sciw;
  if(mod==15){
        INGRID_Dimension *fdim = new INGRID_Dimension();
	sciw = fdim -> get_sciwidth(mod,view,pln,ch,axis);
	sciw = sciw*10.;
	delete fdim;
	//added by koga 2016/7/11
	if(flag==1){sciw=0.5;}
  }
  else if(mod!=16||pln==0||ch<8||ch>=24)sciw=25;
  else sciw=12.5;

  return sciw;
};


float scithick(int mod, int view, int pln, float xy, int axis, int flag=0){
  float scith;
  if(mod==15){
        INGRID_Dimension *fdim = new INGRID_Dimension();
	scith = fdim -> get_scithick(mod,view,pln,0,axis); //mod view pln ch axis
	scith = scith*10.;
	delete fdim;
	//added by koga 2016/7/11
	if(flag==1){scith=0.5;}
  }
  else if(mod!=16||pln==0||xy<400||xy>800)scith=5;
  else scith=6.5;

  return scith;
};


float Yerr(int mod,int view, int pln,int ch, int axis, int edge){
  float yerr;
  if(edge==0)yerr=xyposi(mod,view,pln,ch,axis)-sciwidth(mod,view,pln,ch,axis);
  else       yerr=xyposi(mod,view,pln,ch,axis)+sciwidth(mod,view,pln,ch,axis);
  //if(edge==0)yerr=(xyposi(mod,view,pln,ch,axis)-sciwidth(mod,view,pln,ch,axis))/sqrt(12);
  //else       yerr=(xyposi(mod,view,pln,ch,axis)+sciwidth(mod,view,pln,ch,axis))/sqrt(12);
  //if(edge==0)yerr=xyposi(mod,view,pln,ch,axis)-sciwidth(mod,view,pln,ch,axis,1);
  //else       yerr=xyposi(mod,view,pln,ch,axis)+sciwidth(mod,view,pln,ch,axis,1);
  return yerr;
};






bool actpln[Cview][Cpln];
int nactpln;
int fNactpln(int mod){
  nactpln=0;
  memset(actpln,false,sizeof(actpln));
  for(int i=0; i<(int)hitcls.size(); i++){
    if(hitcls[i].pln<plnmax(mod, 0, 0, 0))actpln[hitcls[i].view][hitcls[i].pln]=true;
  }
  for(int i=0; i<plnmax(mod, 0, 0, 0); i++){
    if(mod!=16){
      if(actpln[0][i]&&actpln[1][i])nactpln++;
      else{
	actpln[0][i]=false;
	actpln[1][i]=false;
      }
    }
    else{
      if(actpln[0][i]&&actpln[1][i])nactpln++;
      else{
	actpln[0][i]=false;
	actpln[1][i]=false;
      }
    }

  }
  return nactpln;
};



float fLayerpe(int mod){
  float layerpe=0;
  for(int i=0; i<(int)hitcls.size(); i++){
    if(hitcls[i].pln<plnmax(mod, 0, 0, 0)){
      if(actpln[hitcls[i].view][hitcls[i].pln])layerpe+=hitcls[i].pe;
    }
  }
  if(nactpln==0)layerpe=0;
  else layerpe=layerpe/nactpln;
  return layerpe;
};


bool withtime(const Hit& left, const Hit& right){
  return left.time < right.time;
};

void fSortTime(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withtime);
};



//Int_t cTdcRsr = 50; //nsec
Int_t cTdcRsr = 50; //nsec
float cPeCut  = 18;


bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, 
		     Long_t &ctime){
  //fSortTime( hit );
  Int_t nhit = hit.size();
  int maxhit=5;
  if(hit[0].mod==MOD_WAGASCI){
	maxhit=1;
  };

  if(nhit <= maxhit)return false;

  
//comment out by koga for cosmic study 
//  for(Int_t i=0 ; i<nhit-maxhit; i++){
//    //if( fabs( hit[i].time - hit[i+maxhit-1].time) < cTdcRsr ){
//    if( fabs( hit[i].time - hit[i+maxhit].time) < cTdcRsr ){
//
//      //long basetime = 0;
//      //for( int j=0; j<maxhit; j++){
//      //for( int j=0; j<maxhit+1; j++){
//	//basetime += hit[i+j].time;
//      //}
//      //basetime = basetime/(maxhit);
//      //basetime = basetime/(maxhit+1);
//     
//      //#### change the definition of basetime   ###
//      //#### ~2010/4/29 mean time of first hits  ###
//      //#### 2010/4/29~ time of highest hit      ###
//      long basetime;    
//      float highpe = 0;
//      float sumpe  = 0;
//      for(int j=0; j<maxhit+1; j++){
//	if( hit[i+j].pe > highpe ){
//	  basetime = hit[i+j].time;
//	  highpe   = hit[i+j].pe;
//	  sumpe    = sumpe + hit[i+j].pe;
//	}
//      }
//      if(hit[0].mod!=MOD_WAGASCI){
//	      if( sumpe < cPeCut )continue;
//      }
//
//      vector<Hit>::iterator it;
//      Int_t ncount=0;
//      for(it = hit.begin() ; it != hit.end(); it++){
//	if( fabs( basetime - it->time) < cTdcRsr ){
//	  ncount++;
//	}
//      }
//      if(ncount<=maxhit)continue;
//
//      //###### clstering ######
//      //#######################
//      hitcls.clear(); //Reset hit clster
//      for(it = hit.begin() ; it != hit.end(); it++){
//	if( fabs( basetime - it->time) < cTdcRsr ){
//	  hitclster.push_back(*it);
//	  it = hit.erase(it);
//	  it--;
//	}
//      }
//
//      //###### caluculate time of clster with max p.e. ######
//      //#####################################################
//      highpe = 0;
//      for(int j=0; j<hitclster.size(); j++){
//	if( highpe < hitclster[j].pe ){
//	  highpe  = hitclster[j].pe;
//	  ctime   = hitclster[j].time;
//	}
//      }
//
//      return true;
//    }
//  }//i
//  return false;


      //#### change the definition of basetime   ###
      //#### ~2010/4/29 mean time of first hits  ###
      //#### 2010/4/29~ time of highest hit      ###
      //###### clstering ######
      //#######################
      vector<Hit>::iterator it;
      hitcls.clear(); //Reset hit clster
      for(it = hit.begin() ; it != hit.end(); it++){
	  hitclster.push_back(*it);
	  it = hit.erase(it);
	  it--;
      }

      //###### caluculate time of clster with max p.e. ######
      //#####################################################
      float highpe = 0;
      for(int j=0; j<(int)hitclster.size(); j++){
	if( highpe < hitclster[j].pe ){
	  highpe  = hitclster[j].pe;
	  ctime   = hitclster[j].time;
	}
      }
      return true;

}

/*
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
*/





int view,id[Cview][Chit],pln[Cview][Chit],ch[Cview][Chit],hit[Cview],hitnum[Cview][Cpln][Cch],used[Cview][Chit];//B2

float pe[Cview][Chit];
int PLN,PLN2,CH,CL,CL2,CLHIT,CLHIT2,CLHIT3,CELL,CELL2,NEI,TRA,TRA2,TRACELL,TRACL,TRACL2,HIT,DIST,DIST2,DIST3,dummy,TMP;

bool hitcl[Cview][Cpln][Ccls];//B2

int clchi[Cview][Cpln][Ccls],clchf[Cview][Cpln][Ccls],ncl[Cview][Cpln],numcl[Cview][Cpln][Ccls],clhit[Cview][Cpln][Ccls][Cch];//B2
float clcenter[Cview][Cpln][Ccls],clpe[Cview][Cpln][Ccls];//B2
int clused[Cview][Cpln][Ccls];//B2

int cellu[Cview][Cpln][Ccell][Cdist],celld[Cview][Cpln][Ccell][Cdist],ncell[Cview][Cpln][Cdist],value[Cview][Cpln][Ccell][Cdist],nvalue[Cview][Cpln][Ccell][Cdist];//B2

bool neibor[Cview][Cpln][Ccell][Cdist][Cdist];//B2

int neiu[Cview][Cpln][Ccell][Cdist][Cdist],neid[Cview][Cpln][Ccell][Cdist][Cdist],nnei[Cview][Cpln][Cdist][Cdist];//B2

int track_cell[Cview][Ccell][Cpln],track_pln[Cview][Ccell][Cpln],track_dist[Cview][Ccell][Cpln],ntracell[Cview][Ccell],ntracl[Cview][Ccell],ntracl2[Cview][Ccell],ntrack[Cview],ntrack2[Cview],ntrack3[Cview],trank[Cview][Cpln][Cpln],ncltra[Cview][Cpln],rank[Cview][Ccell];//B2

bool ttrack[Cview][Ccell],ttrack2[Cview][Ccell];//B2

int plane[Cview][Ctra][Cpln],clus[Cview][Ctra][Cpln];//B2

bool vetowtracking[Cview][Ctra],edgewtracking[Cview][Ctra],track_stop[Cview][Ctra];//B2

float vetodist[Cview][Ctra];

bool ovcl[Cview][Cpln][Ccls];//add//B2

float dis;

float XX[Cview][Ctra][Cpln],YY[Cview][Ctra][Cpln],Xe[Cview][Ctra][Cpln],Yel[Cview][Ctra][Cpln],Yeh[Cview][Ctra][Cpln];//B2

float x[Cdir],y[Cdir],xe[Cdir],yel[Cdir],yeh[Cdir];
float chi2,tmp;

float par[Cview][Ctra][2];//2 is number of fitting parameter(p0,p1)

float trape[Cview][Ctra];//B2

float ulimit;
bool upst,dwst;
//int select;
int nhit;
int ue,shita;

int hit_n[Cview][Chit],hit_id[Cview][Chit][Chit];//B2

bool isohit[Cview][Ctra][Chit];//B2

//double neibor_th = 100.;
//double neibor_th = 10.;
//double neibor_th = 5.;
double neibor_th = 3.;

double joint_th = 25.;



bool fTracking(int mod, int axis=0, int N_long=1){
  memset(hit,0,sizeof(hit));
  memset(hitnum,-1,sizeof(hitnum));
  if(mod!=15){
  	alltrack.clear();
  }
  for(int i=0; i<(int)hitcls.size(); i++){
    view=hitcls[i].view;
    id[view][hit[view]]=hitcls[i].id;
    pln[view][hit[view]]=hitcls[i].pln;
    ch[view][hit[view]]=hitcls[i].ch;
    pe[view][hit[view]]=hitcls[i].pe;
    used[view][hit[view]]=hitcls[i].used;
    //pdg[view][hit[view]]=hitcls[i].pdg;
    if(pln[view][hit[view]]<plnmax(mod,view,pln[view][hit[view]],axis))hitnum[view][pln[view][hit[view]]][ch[view][hit[view]]]=hit[view];
    hit[view]++;
  }

  for(int VIEW=0;VIEW<2;VIEW++){

    //*****Define clusters*****
    memset(ncl[VIEW],0,sizeof(ncl[VIEW]));
    memset(numcl[VIEW],0,sizeof(numcl[VIEW]));
    memset(clhit[VIEW],0,sizeof(clhit[VIEW]));
    for(PLN=0;PLN<plnmax(mod,VIEW,PLN,axis);PLN++){
      CH=0;
      while(1){
	if(hitnum[VIEW][PLN][CH]>=0){
	  clchi[VIEW][PLN][ncl[VIEW][PLN]]=CH;
	  clhit[VIEW][PLN][ncl[VIEW][PLN]][numcl[VIEW][PLN][ncl[VIEW][PLN]]]=hitnum[VIEW][PLN][CH];
	  while(1){
	    numcl[VIEW][PLN][ncl[VIEW][PLN]]++;
	    //don't make cluster WAGASCI
	    //if(mod==MOD_WAGASCI)break;
	    if( (mod==MOD_WAGASCI && !( (axis==0 && PLN%3==0 ) || (axis==1 && PLN%3==1) ) ) )break;
	    if(CH==chmax(mod,VIEW,PLN,axis)-1)break;
	    if(hitnum[VIEW][PLN][CH+1]<0)break;
	    CH++;
	    clhit[VIEW][PLN][ncl[VIEW][PLN]][numcl[VIEW][PLN][ncl[VIEW][PLN]]]=hitnum[VIEW][PLN][CH];
	  }
	  clchf[VIEW][PLN][ncl[VIEW][PLN]]=CH;
	  ncl[VIEW][PLN]++;
	}//if
	if(CH==chmax(mod,VIEW,PLN,axis)-1)break;
	CH++;
      }//while
    }//for PLN

    //*****Fill pe and center of clusters*****
    memset(clpe[VIEW],0,sizeof(clpe[VIEW]));
    memset(clcenter[VIEW],0,sizeof(clcenter[VIEW]));
    for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis);PLN++){
      for(CL=0;CL<ncl[VIEW][PLN];CL++){
        for(CLHIT=0;CLHIT<numcl[VIEW][PLN][CL];CLHIT++){
          clpe[VIEW][PLN][CL]+=pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]];
          clcenter[VIEW][PLN][CL]+=pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]]*xyposi(mod,VIEW,PLN,clchi[VIEW][PLN][CL]+CLHIT,axis);
        }
        clcenter[VIEW][PLN][CL]=clcenter[VIEW][PLN][CL]/clpe[VIEW][PLN][CL];
      }
    }
    //*****Used clusters*****
    memset(clused[VIEW],0,sizeof(clused[VIEW]));
    for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis);PLN++){
      for(CL=0;CL<ncl[VIEW][PLN];CL++){
	for(CLHIT=0;CLHIT<numcl[VIEW][PLN][CL];CLHIT++){
	  if(used[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==1){
	    clused[VIEW][PLN][CL]=1;
	  }
	}
      }
    }
    //*****Define cells*****
    memset(ncell[VIEW],0,sizeof(ncell[VIEW]));
    for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis)-1;PLN++){
      for(CL=0;CL<ncl[VIEW][PLN];CL++){
	//for(DIST=0;DIST<2;DIST++){
	for(DIST=0;DIST<6;DIST++){
	  //if(PLN==plnmax(mod, VIEW, PLN, axis)-2&&DIST==1)continue;
	  if(mod!=15 && DIST>1)continue;
	  if(PLN + DIST > plnmax(mod, VIEW, PLN, axis)-2)continue;
	  for(CL2=0;CL2<ncl[VIEW][PLN+DIST+1];CL2++){
	    //if(DIST==1&&fabs(clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+2][CL2])>92)continue;
	    //if(DIST==1&&fabs(clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+2][CL2])>150)continue;
	    if( mod==15 && fabs(clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+DIST+1][CL2])>120 )continue;
	    else if(DIST==1&&fabs(clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+2][CL2])>150)continue;
	    //added by koga 2016/8/4
	    if( clused[VIEW][PLN][CL]==1 || clused[VIEW][PLN+DIST+1][CL2]==1 )continue;
	    cellu[VIEW][PLN][ncell[VIEW][PLN][DIST]][DIST]=CL;
	    celld[VIEW][PLN][ncell[VIEW][PLN][DIST]][DIST]=CL2;
	    ncell[VIEW][PLN][DIST]++;
	  }
	}
      }
    }

    //*****Define neiborhoods*****
    memset(nnei[VIEW],0,sizeof(nnei[VIEW]));
    for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis)-2;PLN++){
      //for(DIST=0;DIST<2;DIST++){
      for(DIST=0;DIST<6;DIST++){
	//if(PLN==0&&DIST==1)continue;
	if(PLN-DIST<0)continue;
	for(CELL=0;CELL<ncell[VIEW][PLN-DIST][DIST];CELL++){
	  //for(DIST2=0;DIST2<2;DIST2++){
	  for(DIST2=0;DIST2<6;DIST2++){
	    //if(PLN==plnmax(mod, VIEW, PLN, axis)-3&&DIST2==1)continue;
	    if(PLN+DIST2>plnmax(mod, VIEW, PLN, axis)-3)continue;

	    if(mod!=15 && DIST==1&&DIST2==1)continue;
	    if(mod!=15 && DIST>1)continue;
	    if(mod!=15 && DIST2>1)continue;
	    //if(DIST>2&&DIST2>2)continue;//comment out by koga 2016/7/11

	    for(CELL2=0;CELL2<ncell[VIEW][PLN+1][DIST2];CELL2++){
	      if(celld[VIEW][PLN-DIST][CELL][DIST]==cellu[VIEW][PLN+1][CELL2][DIST2]){
		
		x[0]=zposi(mod,VIEW,PLN-DIST,0,axis);//mod view pln ch axis
		x[1]=zposi(mod,VIEW,PLN+1,0,axis);
		x[2]=zposi(mod,VIEW,PLN+2+DIST2,0,axis);
		
		y[0]=clcenter[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]];
		y[1]=clcenter[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]];
		y[2]=clcenter[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]];
	        if(DIST>2&&DIST2>2&&(fabs(y[2]-y[0])>60))continue;//added by koga 2016/7/11

		xe[0]=scithick(mod,VIEW,PLN-DIST,y[0],axis,0);
		xe[1]=scithick(mod,VIEW,PLN+1,   y[1],axis,0);
		xe[2]=scithick(mod,VIEW,PLN+2+DIST2,y[2],axis,0);

		yel[0]=y[0]-Yerr(mod,VIEW,PLN-DIST,clchi[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]],axis,0);
		yel[1]=y[1]-Yerr(mod,VIEW,PLN+1,clchi[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]],axis,0);
		yel[2]=y[2]-Yerr(mod,VIEW,PLN+2+DIST2,clchi[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]],axis,0);
		
		yeh[0]=-y[0]+Yerr(mod,VIEW,PLN-DIST,clchf[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]],axis,1);
		yeh[1]=-y[1]+Yerr(mod,VIEW,PLN+1,clchf[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]],axis,1);
		yeh[2]=-y[2]+Yerr(mod,VIEW,PLN+2+DIST2,clchf[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]],axis,1);
		
		TGraphAsymmErrors *graph=new TGraphAsymmErrors(3,x,y,xe,xe,yel,yeh);
		TF1 *f=new TF1("f","[0]+[1]*x");
		f->SetParameters(y[0]-x[0]*(y[2]-y[0])/(x[2]-x[0]),(y[2]-y[0])/(x[2]-x[0]));
		
		graph->Fit("f","Q");
		chi2=f->GetChisquare();
		//std::cout << chi2 << std::endl;
		graph->Delete();
		f->Delete();
		
		/*		
		if(mod==16){
		  if(DIST==0&&DIST2==0)ulimit=500;
		  else ulimit==280;
		}
		else{
		  if(DIST==0&&DIST2==0)ulimit=2200;
		  else ulimit==280;
		}
		*/
		
		//if(DIST==0&&DIST2==0)ulimit=3;
		if(mod==15)ulimit=neibor_th;
		else ulimit=3;

		if(chi2<ulimit){
		  neiu[VIEW][PLN][nnei[VIEW][PLN][DIST][DIST2]][DIST][DIST2]=CELL;
		  neid[VIEW][PLN][nnei[VIEW][PLN][DIST][DIST2]][DIST][DIST2]=CELL2;
		  nnei[VIEW][PLN][DIST][DIST2]++;
		  //break;
		}
		
		
	      }
	    }
	  }
	}
      }
    }



    //*****Define value of cells*****
    memset(value[VIEW],0,sizeof(value[VIEW]));
    memset(nvalue[VIEW],0,sizeof(nvalue[VIEW]));
    for(int jndex=0;jndex<MAX_CELLVALUE;jndex++){
      for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis)-2;PLN++){
	//for(DIST=0;DIST<2;DIST++){
	for(DIST=0;DIST<6;DIST++){
	  //for(DIST2=0;DIST2<2;DIST2++){
	  for(DIST2=0;DIST2<6;DIST2++){
	    //if(PLN==0&&DIST==1)continue;
	    if(PLN-DIST<0)continue;
	    //if(PLN==plnmax(mod, VIEW, PLN, axis)-3&&DIST2==1)continue;
	    if(PLN+DIST2>plnmax(mod, VIEW, PLN, axis)-3)continue;

	    //if(DIST==1&&DIST2==1)continue;	    
	    //if(DIST>2&&DIST2>2)continue;//comment out by koga 2016/7/11    

	    for(NEI=0;NEI<nnei[VIEW][PLN][DIST][DIST2];NEI++){
	      if(value[VIEW][PLN-DIST][neiu[VIEW][PLN][NEI][DIST][DIST2]][DIST]==value[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]){
		nvalue[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]=value[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]+1;
	      }
	    }
	  }
	}
      }
      for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis)-1;PLN++){
        for(DIST=0;DIST<6;DIST++){
          //if(PLN==plnmax(mod, VIEW, PLN, axis)-2&&DIST==1)continue;
          if(PLN+DIST>plnmax(mod, VIEW, PLN, axis)-2)continue;
          for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
            value[VIEW][PLN][CELL][DIST]=nvalue[VIEW][PLN][CELL][DIST];
          }
        }
      }
    }

    //*****Define tracks*****
    ntrack[VIEW]=0;
    memset(ntracell[VIEW],0,sizeof(ntracell[VIEW]));
    memset(trape[VIEW],0,sizeof(trape[VIEW]));
    memset(neibor[VIEW],true,sizeof(neibor[VIEW]));
    for(PLN=1;PLN<plnmax(mod, VIEW, PLN, axis)-1;PLN++){
      //for(DIST=0;DIST<2;DIST++){
      for(DIST=0;DIST<6;DIST++){
	//if(PLN==plnmax(mod, VIEW, PLN, axis)-2&&DIST==1)continue;
	if(PLN+DIST>plnmax(mod, VIEW, PLN, axis)-2)continue;
	for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
	  if(value[VIEW][PLN][CELL][DIST]>0){
	    if(PLN+DIST==plnmax(mod, VIEW, PLN, axis)-2){upst=true;}
	    else{
	      upst=true;
	      //for(DIST2=0;DIST2<2;DIST2++){
	      for(DIST2=0;DIST2<6;DIST2++){
		if(mod!=15 && DIST==1&&DIST2==1)continue;
		if(mod!=15 && DIST>1)continue;
		if(mod!=15 && DIST2>1)continue;
		//if(DIST>2&&DIST2>2)continue;//comment out by koga 2016/7/11
		for(NEI=0;NEI<nnei[VIEW][PLN+DIST][DIST][DIST2];NEI++){
		  //if(cellu[PLN][CELL]==cellu[PLN][neiu[PLN][NEI]])upst=false;
		  if(CELL==neiu[VIEW][PLN+DIST][NEI][DIST][DIST2])upst=false;
		  if(CELL==neiu[VIEW][PLN+DIST][NEI][DIST][DIST2]&&neibor[VIEW][PLN+DIST][NEI][DIST][DIST2]){upst=true;break;}
		}
	      }
	    }
	    if(upst){

	      track_pln[VIEW][ntrack[VIEW]][0]=PLN;
	      track_cell[VIEW][ntrack[VIEW]][0]=CELL;
	      track_dist[VIEW][ntrack[VIEW]][0]=DIST;

	      plane[VIEW][ntrack[VIEW]][0]=track_pln[VIEW][ntrack[VIEW]][0]+track_dist[VIEW][ntrack[VIEW]][0]+1;
	      clus[VIEW][ntrack[VIEW]][0]=celld[VIEW][track_pln[VIEW][ntrack[VIEW]][0]][track_cell[VIEW][ntrack[VIEW]][0]][track_dist[VIEW][ntrack[VIEW]][0]];
	      XX[VIEW][ntrack[VIEW]][0]=zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],0,axis);
	      YY[VIEW][ntrack[VIEW]][0]=clcenter[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]];
	      trape[VIEW][ntrack[VIEW]]+=clpe[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]];

	      Yel[VIEW][ntrack[VIEW]][0]=YY[VIEW][ntrack[VIEW]][0] -Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],clchi[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]],axis,0);
	      Yeh[VIEW][ntrack[VIEW]][0]=-YY[VIEW][ntrack[VIEW]][0]+Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],clchf[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]],axis,1);
	      Xe[VIEW][ntrack[VIEW]][0]=scithick(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],YY[VIEW][ntrack[VIEW]][0],axis,0);
	      
	      plane[VIEW][ntrack[VIEW]][1]=track_pln[VIEW][ntrack[VIEW]][0];
	      clus[VIEW][ntrack[VIEW]][1]=cellu[VIEW][track_pln[VIEW][ntrack[VIEW]][0]][track_cell[VIEW][ntrack[VIEW]][0]][track_dist[VIEW][ntrack[VIEW]][0]];
	      XX[VIEW][ntrack[VIEW]][1]=zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],0,axis);
	      YY[VIEW][ntrack[VIEW]][1]=clcenter[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]];
	      trape[VIEW][ntrack[VIEW]]+=clpe[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]];
	      
	      Yel[VIEW][ntrack[VIEW]][1]=YY[VIEW][ntrack[VIEW]][1] -Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],clchi[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]],axis,0);
	      Yeh[VIEW][ntrack[VIEW]][1]=-YY[VIEW][ntrack[VIEW]][1]+Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],clchf[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]],axis,1);
	      Xe[VIEW][ntrack[VIEW]][1]=scithick(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],YY[VIEW][ntrack[VIEW]][1],axis,0);


	      ntracell[VIEW][ntrack[VIEW]]=1;
	      PLN2=PLN-1;
	      int chi2old=3;
	      if(mod==15)chi2old=neibor_th;
	      while(PLN2>=0){
		dwst=true;
		DIST3=track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]-1];
		//for(DIST2=0;DIST2<2;DIST2++){
		for(DIST2=0;DIST2<6;DIST2++){
		  //if(DIST3==1&&DIST2==1)continue;
		  //if(DIST2==1&&PLN2==0)continue;
		  if(PLN2-DIST2<0)continue;

		  for(TMP=0;TMP<2;TMP++){
		    for(NEI=0;NEI<nnei[VIEW][PLN2][DIST2][DIST3];NEI++){
		      if(neid[VIEW][PLN2][NEI][DIST2][DIST3]==track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]-1]){
			if(value[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]+1==value[VIEW][PLN2+1][neid[VIEW][PLN2][NEI][DIST2][DIST3]][DIST3]||TMP!=0){
			  if(!dwst){
			    
			    float Xetmp=Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];
			    float YYtmp=YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];
			    float Yeltmp=Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];
			    float Yehtmp=Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];


			    YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=clcenter[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]];
			    Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] -Yerr(mod,VIEW,PLN2-DIST2,clchi[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]],axis,0);
			    Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=-YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]+Yerr(mod,VIEW,PLN2-DIST2,clchf[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]],axis,1);
			    Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=scithick(mod,VIEW,PLN2-DIST2,YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],axis,0);


			    //TGraph *graph2=new TGraph((ntracell[VIEW][ntrack[VIEW]]+2),XX[VIEW][ntrack[VIEW]],YY[VIEW][ntrack[VIEW]]);
			    TGraphAsymmErrors *graph2=new TGraphAsymmErrors((ntracell[VIEW][ntrack[VIEW]]+2),XX[VIEW][ntrack[VIEW]],YY[VIEW][ntrack[VIEW]],Xe[VIEW][ntrack[VIEW]],Xe[VIEW][ntrack[VIEW]],Yel[VIEW][ntrack[VIEW]],Yeh[VIEW][ntrack[VIEW]]);
			    TF1 *f2=new TF1("f2","[0]+[1]*x");
			    f2->SetParameters(YY[VIEW][ntrack[VIEW]][0]-XX[VIEW][ntrack[VIEW]][0]*(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])/(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]),(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])/(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]));
 
			    graph2->Fit("f2","Q");
			    tmp=f2->GetChisquare();

			    graph2->Delete();
			    f2->Delete();
			    
			    float upstpe1=clpe[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]][DIST2]];
			    float upstpe2=clpe[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]];
			    
			    //if((tmp>chi2||upstpe2<4.5)&&upstpe1>=4.5){
			    if(   (mod<14  && (tmp>chi2||upstpe2<4.5)&&upstpe1>=4.5 ) ||
			          (mod==16 && (tmp>chi2||upstpe2<4.5)&&upstpe1>=4.5 ) ||
			          (mod==15 && (tmp>chi2) && upstpe1>=1.5 )  )
			    {
			      Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=Xetmp;
			      YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=YYtmp;
			      Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=Yeltmp;
			      Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=Yehtmp;
			      continue;
			    }
			    
			  }
			  

			  track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]=neiu[VIEW][PLN2][NEI][DIST2][DIST3];
			  track_pln[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]=PLN2-DIST2;
			  track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]=DIST2;
			  
			  plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=track_pln[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]];
			  clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=cellu[VIEW][track_pln[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]][track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]][track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]];
			  XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],0,axis);			  
			  YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=clcenter[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]];

			  Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] -Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],clchi[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]],axis,0);
			  Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=-YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]+Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],clchf[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]],axis,1);
			  Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=scithick(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],axis,0);
			  
			  //TGraph *graph1=new TGraph((ntracell[VIEW][ntrack[VIEW]]+2),XX[VIEW][ntrack[VIEW]],YY[VIEW][ntrack[VIEW]]);
			  TGraphAsymmErrors *graph1=new TGraphAsymmErrors((ntracell[VIEW][ntrack[VIEW]]+2),XX[VIEW][ntrack[VIEW]],YY[VIEW][ntrack[VIEW]],Xe[VIEW][ntrack[VIEW]],Xe[VIEW][ntrack[VIEW]],Yel[VIEW][ntrack[VIEW]],Yeh[VIEW][ntrack[VIEW]]);
			  TF1 *f1=new TF1("f1","[0]+[1]*x");
			  f1->SetParameters(YY[VIEW][ntrack[VIEW]][0]-XX[VIEW][ntrack[VIEW]][0]*(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])/(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]),(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])/(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]));
			  
			  graph1->Fit("f1","Q");
			  chi2=f1->GetChisquare();
			  graph1->Delete();
			  f1->Delete();


			  if(
			     (mod<14&&((chi2/ntracell[VIEW][ntrack[VIEW]]>3&&ntracell[VIEW][ntrack[VIEW]]==1)||
					(chi2/ntracell[VIEW][ntrack[VIEW]]>2&&ntracell[VIEW][ntrack[VIEW]]>1)||
					((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 -(ntracell[VIEW][ntrack[VIEW]]-1)*0.02)))
			     ||
			     (mod==16&&((chi2/ntracell[VIEW][ntrack[VIEW]]>3&&ntracell[VIEW][ntrack[VIEW]]==1)||
					(chi2/ntracell[VIEW][ntrack[VIEW]]>2&&ntracell[VIEW][ntrack[VIEW]]>1)||
					((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 -(ntracell[VIEW][ntrack[VIEW]]-1)*0.015)))
			     ||
			     (mod==15&&((chi2/ntracell[VIEW][ntrack[VIEW]]>neibor_th&&ntracell[VIEW][ntrack[VIEW]]==1)||
					//(chi2/ntracell[VIEW][ntrack[VIEW]]>10&&ntracell[VIEW][ntrack[VIEW]]>1)||
					((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 +(ntracell[VIEW][ntrack[VIEW]]-1)*0.1)))
					//((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>10 +(ntracell[VIEW][ntrack[VIEW]]-1)*10)))
			     ){
			    if(value[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]>2&&TMP==0){
			      neibor[VIEW][PLN2][NEI][DIST2][DIST3]=false;
			    }
			    continue;
			  }

			  dwst=false;
			}
		      }//if(!dwst)break;
		    }//NEI
		    if(!dwst)break;
		  }//TMP  
		  if(!dwst)break;
		}
		if(dwst)break;
		PLN2=PLN2-track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]-1;
		chi2old=chi2/ntracell[VIEW][ntrack[VIEW]];
		trape[VIEW][ntrack[VIEW]]+=clpe[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]];
		ntracell[VIEW][ntrack[VIEW]]++;
	      }
	      if(ntracell[VIEW][ntrack[VIEW]]==1)continue;
	      ntracl[VIEW][ntrack[VIEW]]=ntracell[VIEW][ntrack[VIEW]]+1;
	      ntrack[VIEW]++;
	    }
	  }
	}
      }
    }

    //*****Track rank*****
    memset(trank[VIEW],0,sizeof(trank[VIEW]));
    memset(rank[VIEW],0,sizeof(rank[VIEW]));
    memset(ncltra[VIEW],0,sizeof(ncltra[VIEW]));
    ntrack2[VIEW]=0;
    for(TRACL=MAX_CELLVALUE;TRACL>0;TRACL--){
      for(TRA=0;TRA<ntrack[VIEW];TRA++){
        if(ntracl[VIEW][TRA]==TRACL){
          trank[VIEW][TRACL][ncltra[VIEW][TRACL]]=TRA;
          ncltra[VIEW][TRACL]++;
        }
      }
      for(CL=0;CL<ncltra[VIEW][TRACL];CL++){
        for(CL2=CL+1;CL2<ncltra[VIEW][TRACL];CL2++){
          if(trape[VIEW][trank[VIEW][TRACL][CL]]<trape[VIEW][trank[VIEW][TRACL][CL2]]){
            dummy=trank[VIEW][TRACL][CL];
            trank[VIEW][TRACL][CL]=trank[VIEW][TRACL][CL2];
            trank[VIEW][TRACL][CL2]=dummy;
          }
        }
        rank[VIEW][ntrack2[VIEW]]=trank[VIEW][TRACL][CL];
        ntrack2[VIEW]++;
      }
    }

    //*****True track selection*****
    memset(ttrack[VIEW],true,sizeof(ttrack[VIEW]));
    //memset(ttrack[VIEW],false,sizeof(ttrack[VIEW]));
    memset(hitcl[VIEW],false,sizeof(hitcl[VIEW]));
    memset(ntracl2[VIEW],0,sizeof(ntracl2[VIEW]));
    memset(ovcl[VIEW],false,sizeof(ovcl[VIEW]));//add
    ntrack2[VIEW]=0;
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      TRA2=rank[VIEW][TRA];
      //std::cout << "debug " << axis << " " << VIEW << " " << TRA2 << std::endl;
      for(TRACL=0;TRACL<ntracl[VIEW][TRA2];TRACL++){
	if(!hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]
	   //&&clpe[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]>4.5){
	   //&&clpe[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]>4.5
	   &&clused[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]==0){
	   //){
	  ntracl2[VIEW][TRA2]++;
	}
	//added 2015/6/28
	else if(mod==15 && TRACL!=ntracl[VIEW][TRA2]-1){
	  //ntracl2[VIEW][TRA2] = ntracl2[VIEW][TRA2] -3;
	  ntracl2[VIEW][TRA2] = ntracl2[VIEW][TRA2] -1;
	}
      }

/*
      if(ntracl2[VIEW][TRA2]==2){
    for(TRACL2=0;TRACL2<ntracl[VIEW][TRA2];TRACL2++){
      if(!hitcl[VIEW][plane[VIEW][TRA2][TRACL2]][clus[VIEW][TRA2][TRACL2]]
         &&clpe[VIEW][plane[VIEW][TRA2][TRACL2]][clus[VIEW][TRA2][TRACL2]]>4.5)break;
    }
    TRACL2++;
    if(hitcl[VIEW][plane[VIEW][TRA2][TRACL2]][clus[VIEW][TRA2][TRACL2]]
       //&&clpe[VIEW][plane[VIEW][TRA2][TRACL2]][clus[VIEW][TRA2][TRACL2]]>6.0
       )ntracl2[VIEW][TRA2]=1;
      }
      */
/*
      if(ntracl2[VIEW][TRA2]==2){
    if(ntracl[VIEW][TRA2]<4){
      if(!hitcl[VIEW][plane[VIEW][TRA2][0]][clus[VIEW][TRA2][0]]
         &&clpe[VIEW][plane[VIEW][TRA2][0]][clus[VIEW][TRA2][0]]>4.5
         &&!hitcl[VIEW][plane[VIEW][TRA2][ntracl[VIEW][TRA2]-1]][clus[VIEW][TRA2][ntracl[VIEW][TRA2]-1]]
         &&clpe[VIEW][plane[VIEW][TRA2][ntracl[VIEW][TRA2]-1]][clus[VIEW][TRA2][ntracl[VIEW][TRA2]-1]]>4.5)
        ntracl2[VIEW][TRA2]=1;
    }
      }
*/
      //if(ntracl2[VIEW][TRA2]>1){
      if(ntracl2[VIEW][TRA2]>N_long){
	for(TRACL=0;TRACL<ntracl[VIEW][TRA2];TRACL++){
	  if(hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]){
		  ovcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]=true;//add
	  }
	  hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]=true;
	}
	ntrack2[VIEW]++;
        //ttrack[VIEW][TRA2]=true;
      }
      else ttrack[VIEW][TRA2]=false;
    }

    //for(TRA=0;TRA<ntrack[VIEW];TRA++){
    //  if(ttrack[VIEW][TRA]){
    //    std::cout << "debug2 " << axis << " " << VIEW << " " << TRA << std::endl;
    //  }
    //}

    //*****Fit tracks*****
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(!ttrack[VIEW][TRA])continue;
      TGraphAsymmErrors *g=new TGraphAsymmErrors(ntracell[VIEW][TRA]+1,XX[VIEW][TRA],YY[VIEW][TRA],Xe[VIEW][TRA],Xe[VIEW][TRA],Yel[VIEW][TRA],Yeh[VIEW][TRA]);
      TF1 *f=new TF1("f","[0]+[1]*x");
      f->SetParameters(YY[VIEW][TRA][0]-XX[VIEW][TRA][0]*(YY[VIEW][TRA][ntracell[VIEW][TRA]]-YY[VIEW][TRA][0])/(XX[VIEW][TRA][ntracell[VIEW][TRA]]-XX[VIEW][TRA][0]),(YY[VIEW][TRA][ntracell[VIEW][TRA]]-YY[VIEW][TRA][0])/(XX[VIEW][TRA][ntracell[VIEW][TRA]]-XX[VIEW][TRA][0]));
      g->Fit("f","Q");
      TF1 *func=g->GetFunction("f");
      par[VIEW][TRA][0]=func->GetParameter(0);
      par[VIEW][TRA][1]=func->GetParameter(1);
      /*
      chi2=func->GetChisquare();
      int ndf=func->GetNDF();
      if(chi2/ndf>2500)ttrack[VIEW][TRA]=false;
      */
      g->Delete();
      f->Delete();
    }

    //*****upstream veto & edge channel cut*****
    memset(vetowtracking[VIEW],false,sizeof(vetowtracking[VIEW]));
    memset(edgewtracking[VIEW],false,sizeof(edgewtracking[VIEW]));
    memset(track_stop[VIEW],true,sizeof(track_stop[VIEW]));
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(!ttrack[VIEW][TRA])continue;
      if((mod<15&&plane[VIEW][TRA][ntracl[VIEW][TRA]-1]==0) || (mod==15&&plane[VIEW][TRA][ntracl[VIEW][TRA]-1]<=2) || (mod==16&&plane[VIEW][TRA][ntracl[VIEW][TRA]-1]<=1) ){vetowtracking[VIEW][TRA]=true;}
      if(plane[VIEW][TRA][0]==plnmax(mod, VIEW, plane[VIEW][TRA][0] , axis)-1){track_stop[VIEW][TRA]=false;}

#ifdef PMEDGE
      for(TRACL=ntracl[VIEW][TRA]-2;TRACL<ntracl[VIEW][TRA];TRACL++){
	for(CLHIT=0;CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]];CLHIT++){
	  if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==0){edgewtracking[VIEW][TRA]=true;}
	  else if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==chmax(mod,VIEW,plane[VIEW][TRA][TRACL],axis)-1){edgewtracking[VIEW][TRA]=true;}
	}
      }

      for(TRACL=0;TRACL<2;TRACL++){
	for(CLHIT=0;CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]];CLHIT++){
	  if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==0){track_stop[VIEW][TRA]=false;}
	  else if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==chmax(mod,VIEW,plane[VIEW][TRA][TRACL],axis)-1){track_stop[VIEW][TRA]=false;}
	}
      }
#else
      //if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]]))+par[VIEW][TRA][0]<100){edgewtracking[VIEW][TRA]=true;}
      //if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]]))+par[VIEW][TRA][0]>1100){edgewtracking[VIEW][TRA]=true;}
      //if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0]<200){edgewtracking[VIEW][TRA]=true;}
      //if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0]>1000){edgewtracking[VIEW][TRA]=true;}
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0]<200){edgewtracking[VIEW][TRA]=true;}
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0]>1000){edgewtracking[VIEW][TRA]=true;}

      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0]<100){track_stop[VIEW][TRA]=false;}
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0]>1100){track_stop[VIEW][TRA]=false;}

#endif

    }


    //*****side veto*****
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      vetodist[VIEW][TRA]=1e5;
      if(!ttrack[VIEW][TRA])continue;
      for(HIT=0;HIT<hit[VIEW];HIT++){
	if(pln[VIEW][HIT]>=plnmax(mod, VIEW, pln[VIEW][HIT], axis)){
	  if(par[VIEW][TRA][1]>0){
	    //if(vposixy(mod,pln[VIEW][HIT])==1255)continue;
	    if(vposixy(mod,pln[VIEW][HIT])>0)continue;
	    //if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    //if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod<15)continue;
	  }
	  else if(par[VIEW][TRA][1]<0){
	    //if(vposixy(mod,pln[VIEW][HIT])==-55)continue;
	    if(vposixy(mod,pln[VIEW][HIT])<0)continue;
	    //if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    //if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod<15)continue;
	  }
	  else continue;
	  if(mod==15)continue;
	  dis=fabs(par[VIEW][TRA][1]*vposiz(mod,ch[VIEW][HIT])-vposixy(mod,pln[VIEW][HIT])+par[VIEW][TRA][0])/sqrt((par[VIEW][TRA][1])*(par[VIEW][TRA][1])+1);
	  if(dis<80)vetowtracking[VIEW][TRA]=true;
	  if(dis<vetodist[VIEW][TRA])vetodist[VIEW][TRA]=dis;
	}
      }
    }


    //*****mod fc with veto*****
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(!ttrack[VIEW][TRA])continue;
      for(HIT=0;HIT<hit[VIEW];HIT++){
	if(pln[VIEW][HIT]>=plnmax(mod, VIEW, pln[VIEW][HIT], axis)){
	  if(par[VIEW][TRA][1]>0){
	    //if(vposixy(mod,pln[VIEW][HIT])==1255)continue;
	    if(vposixy(mod,pln[VIEW][HIT])<0)continue;
	    if(XX[VIEW][TRA][0]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][0]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	  }
	  else if(par[VIEW][TRA][1]<0){
	    //if(vposixy(mod,pln[VIEW][HIT])==-55)continue;
	    if(vposixy(mod,pln[VIEW][HIT])>0)continue;
	    if(XX[VIEW][TRA][0]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][0]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	  }
	  else continue;
	  if(mod==15)continue;
	  dis=fabs(par[VIEW][TRA][1]*vposiz(mod,ch[VIEW][HIT])-vposixy(mod,pln[VIEW][HIT])+par[VIEW][TRA][0])/sqrt((par[VIEW][TRA][1])*(par[VIEW][TRA][1])+1);
	  if(dis<80)track_stop[VIEW][TRA]=false;
	}
      }
    }



//    //*****cut low upstream hit****
//    for(int TRA=0;TRA<ntrack[VIEW];TRA++){
//      if(mod==MOD_WAGASCI)continue;
//      if(vetowtracking[VIEW][TRA]||edgewtracking[VIEW][TRA])continue;//new
//      if(ttrack[VIEW][TRA]&&clpe[VIEW][plane[VIEW][TRA][ntracl[VIEW][TRA]-1]][clus[VIEW][TRA][ntracl[VIEW][TRA]-1]]<4.5){
//  if(ntracl[VIEW][TRA]==3){
//    if(ovcl[VIEW][plane[VIEW][TRA][0]][clus[VIEW][TRA][0]]||
//       ovcl[VIEW][plane[VIEW][TRA][1]][clus[VIEW][TRA][1]]||
//       ovcl[VIEW][plane[VIEW][TRA][2]][clus[VIEW][TRA][2]]
//       ){
//      ttrack[VIEW][TRA]=false;
//      ntrack2[VIEW]--;
//      continue;
//    }
//    else continue;
//  }
//  ntracl[VIEW][TRA]--;
//  ntracell[VIEW][TRA]--;
//      }
//    }



    //*****cluster P.E.*****
    /*
    memset(cluspe[VIEW],0,sizeof(cluspe[VIEW]));
    memset(clpdg[VIEW],0,sizeof(clpdg[VIEW]));
    for(PLN=0;PLN<plnmax(mod);PLN++){
      for(CL=0;CL<ncl[VIEW][PLN];CL++){
    for(CLHIT=0;CLHIT<numcl[VIEW][PLN][CL];CLHIT++){
      cluspe[VIEW][PLN][CL]+=pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]];
      if(pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==13||pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==-13){
        pid=0;
        clpdg[VIEW][PLN][CL][pid]++;
      }
      else if(pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==2212||pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==-2212){
        pid=1;
        clpdg[VIEW][PLN][CL][pid]++;
      }
      else if(pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==211||pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==-211){
        pid=2;
        clpdg[VIEW][PLN][CL][pid]++;
      }
      else if(pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==11||pdg[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==-11){
        pid=3;
        clpdg[VIEW][PLN][CL][pid]++;
      }

    }
      }
    }
    */


    //****track P.E.****
    /*
    memset(trackpe[VIEW],0,sizeof(trackpe[VIEW]));
    memset(trackpdg[VIEW],0,sizeof(trackpdg[VIEW]));
    memset(len[VIEW],0,sizeof(len[VIEW]));
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]){
    //for(TRACELL=0;TRACELL<ntracell[VIEW][TRA]+1;TRACELL++){
    for(TRACELL=0;TRACELL<ntracell[VIEW][TRA];TRACELL++){//changed to reject most upstream hit
      if(!ovcl[VIEW][plane[VIEW][TRA][TRACELL]][clus[VIEW][TRA][TRACELL]]){
        trackpe[VIEW][TRA]+=cluspe[VIEW][plane[VIEW][TRA][TRACELL]][clus[VIEW][TRA][TRACELL]];
        len[VIEW][TRA]++;
        for(TMP=0;TMP<4;TMP++)
          trackpdg[VIEW][TRA][TMP]+=clpdg[VIEW][plane[VIEW][TRA][TRACELL]][clus[VIEW][TRA][TRACELL]][TMP];

      }
    }
      }
    }
    */


    //****Track hits****
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]){
	TMP=0;
	for(TRACL=0;TRACL<ntracl[VIEW][TRA];TRACL++){
	  for(CLHIT=0;CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]];CLHIT++){
	    hit_id[VIEW][TRA][TMP]=id[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]];
	    if(ovcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]])
	      isohit[VIEW][TRA][TMP]=false;
	    else	    
	      isohit[VIEW][TRA][TMP]=true;
	    TMP++;
	  }
	}
	hit_n[VIEW][TRA]=TMP;
      }
    }



   //Arrange reconstructed hits order in  a track
   for(TRA=0;TRA<ntrack[VIEW];TRA++){
     if(ttrack[VIEW][TRA]){
	int id_temp[Chit];
	for(int TMP=0; TMP<hit_n[VIEW][TRA]; TMP++){
	  id_temp[TMP]=hit_id[VIEW][TRA][hit_n[VIEW][TRA]-1 - TMP];
	}
	for(int TMP=0; TMP<hit_n[VIEW][TRA]; TMP++){
	  hit_id[VIEW][TRA][TMP]=id_temp[TMP];
        }
     }
   }


    //**Joint hit around track**
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]){
	TMP=hit_n[VIEW][TRA];
        //for(PLN=0;PLN<plnmax(mod, VIEW, PLN, axis);PLN++){
        for(PLN=plane[VIEW][TRA][ntracell[VIEW][TRA]]; PLN<=plane[VIEW][TRA][0]; PLN++){
        for(CL=0;CL<ncl[VIEW][PLN];CL++){
	  if( clused[VIEW][PLN][CL]==0 && !hitcl[VIEW][PLN][CL] ){
          //std::cout << "judge joint " <<  fabs(par[VIEW][TRA][1]*zposi(mod,VIEW,PLN,0,axis) -clcenter[VIEW][PLN][CL] +par[VIEW][TRA][0])/sqrt(par[VIEW][TRA][1]+1.)  << std::endl;
	    if( fabs(par[VIEW][TRA][1]*zposi(mod,VIEW,PLN,0,axis) -clcenter[VIEW][PLN][CL] +par[VIEW][TRA][0])/sqrt(par[VIEW][TRA][1]+1.) < joint_th ){
	      for(CLHIT=0;CLHIT<numcl[VIEW][PLN][CL];CLHIT++){
	        hit_id[VIEW][TRA][TMP]=id[VIEW][clhit[VIEW][PLN][CL][CLHIT]];
	        hitcl[VIEW][PLN][CL]=true;
	        TMP++;
	      }
	    }
	  }
	}
	}
	hit_n[VIEW][TRA]=TMP;
      }
    }

    //*****reconstructed track info****
    Track track;
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]){
        track.view = VIEW;
        track.fpln = plane[VIEW][TRA][0];
        track.ipln = plane[VIEW][TRA][ntracell[VIEW][TRA]];
        track.fxy  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0];
        track.ixy  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0];
        track.fz   = zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis);
        track.iz   = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis);
        track.intcpt=par[VIEW][TRA][0];
        track.slope= par[VIEW][TRA][1];
        track.ang  = atan(par[VIEW][TRA][1])*180/3.141592;
        track.veto = vetowtracking[VIEW][TRA];
        track.edge = edgewtracking[VIEW][TRA];
        track.stop = track_stop[VIEW][TRA];
        track.vetodist = vetodist[VIEW][TRA];
        //track.totpe = trackpe[VIEW][TRA];
        //track.isohit = len[VIEW][TRA];
        //for(TMP=0;TMP<4;TMP++)
        //track.pdg[TMP] = trackpdg[VIEW][TRA][TMP];
        track.hitid.clear();
        track.isohit.clear();
        for(TMP=0;TMP<hit_n[VIEW][TRA];TMP++){
           track.hitid.push_back(hit_id[VIEW][TRA][TMP]);
           track.isohit.push_back(isohit[VIEW][TRA][TMP]);
        }

        //axis==1
        if(axis==1){

          if(track.slope>=0){
            track.fz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0];
            track.iz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0];
            track.fxy  = zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis);
            track.ixy  = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis);

            INGRID_Dimension *fdim = new INGRID_Dimension();
            int pln,grid,gridch;
            int fch = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][0]][clus[VIEW][TRA][0]][0]];
            int fpln = plane[VIEW][TRA][0];
                fdim -> get_plnch_fromrecon_loli( mod, VIEW, fpln, fch, axis, &pln, &gridch, &grid); //mod view pln ch axis
            if(VIEW==0){
              track.fpln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0)track.fpln = pln*3+1;
              if(grid==1)track.fpln = pln*3+0;
              if(grid==2)track.fpln = pln*3+2;
            }
            int ich = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][ntracell[VIEW][TRA]]][clus[VIEW][TRA][ntracell[VIEW][TRA]]][0]];
            int ipln = plane[VIEW][TRA][ntracell[VIEW][TRA]];
            fdim -> get_plnch_fromrecon_loli( mod, VIEW, ipln, ich, axis, &pln, &gridch, &grid); //mod view pln ch axis
            if(VIEW==0){
              track.ipln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0)track.ipln = pln*3+1;
              if(grid==1)track.ipln = pln*3+0;
              if(grid==2)track.ipln = pln*3+2;
            }
            delete fdim;
          }//if(track.slope>=0)
          else{
            track.iz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0];
            track.fz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0];
            track.ixy   = zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis);
            track.fxy   = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis);

            INGRID_Dimension *fdim = new INGRID_Dimension();
            int pln,grid,gridch;
            int ich = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][0]][clus[VIEW][TRA][0]][0]];
            int ipln = plane[VIEW][TRA][0];
              fdim -> get_plnch_fromrecon_loli( mod, VIEW, ipln, ich, axis, &pln, &gridch, &grid); //mod view pln ch axis
            if(VIEW==0){
              track.ipln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0)track.ipln = pln*3+1;
              if(grid==1)track.ipln = pln*3+0;
              if(grid==2)track.ipln = pln*3+2;
            }

            int fch = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][ntracell[VIEW][TRA]]][clus[VIEW][TRA][ntracell[VIEW][TRA]]][0]];
            int fpln = plane[VIEW][TRA][ntracell[VIEW][TRA]];
            fdim -> get_plnch_fromrecon_loli( mod, VIEW, fpln, fch, axis, &pln, &gridch, &grid); //mod view pln ch axis
            if(VIEW==0){
              track.fpln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0)track.fpln = pln*3+1;
              if(grid==1)track.fpln = pln*3+0;
              if(grid==2)track.fpln = pln*3+2;
            }
            delete fdim;
          }//else

          if(fabs(track.slope)<0.000001)track.slope=100000000.;
          else track.slope = 1./track.slope;       //y=ax+b  a'=-1/a

          track.intcpt= -track.intcpt*track.slope; //y=ax+b  b'=-b/a
          track.ang   = atan(track.slope)*180/3.141592;

        }//axis==1

        //angle cut
        if(mod==15){
          if(axis==0 && (track.ang>60 || track.ang<-60) )continue;
          //if(axis==1 && (track.ang<60 && track.ang>-60) )continue;
        }
        alltrack.push_back(track);
      }
    }//for(TRA)
  }//for(VIEW)

  if(alltrack.size()>0) return true;
  else                  return false;

};


bool fTracking_loli(int mod){

  	alltrack.clear();
	int N_long=5;
	vector<Hit> hit_temp;
	//recon long track along z axis
	hit_temp=hitcls;
        INGRID_Dimension *fdim_temp = new INGRID_Dimension();
	for(int i=0; i<hitcls.size(); i++){
		int reconpln,reconch;
        	fdim_temp -> get_reconplnch_loli( mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 0, &reconpln, &reconch ); //mod view pln axis
		hitcls[i].pln  = reconpln;
		hitcls[i].ch   = reconch;
		hitcls[i].used = 0;
	}
	fTracking(mod,0,N_long);

	//recon long track along xy axis
	for(int i=0;i<alltrack.size();i++){
		//for(int j=0;j<alltrack[i].hitid.size();j++){
		for(int j=1;j<alltrack[i].hitid.size();j++){
      			vector<Hit>::iterator it;
      			for(it = hit_temp.begin() ; it != hit_temp.end(); it++){
				if(it->id==alltrack[i].hitid[j]){
					//hit_temp.erase(it);
					//it--;
					it->used=1;
				}
			}
		}
	}
	hitcls = hit_temp;
	for(int i=0; i<hitcls.size(); i++){
		int reconpln,reconch;
        	fdim_temp -> get_reconplnch_loli( mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 1, &reconpln, &reconch ); //mod view pln axis
		hitcls[i].pln= reconpln;
		hitcls[i].ch = reconch;
	}
	fTracking(mod,1,N_long);





	//recon along z axis
	for(int i=0;i<alltrack.size();i++){
		//for(int j=0;j<alltrack[i].hitid.size();j++){
		for(int j=1;j<alltrack[i].hitid.size();j++){
      			vector<Hit>::iterator it;
      			for(it = hit_temp.begin() ; it != hit_temp.end(); it++){
				if(it->id==alltrack[i].hitid[j]){
					//hit_temp.erase(it);
					//it--;
					it->used=1;
				}
			}
		}
	}
	hitcls = hit_temp;
	for(int i=0; i<hitcls.size(); i++){
		int reconpln,reconch;
        	fdim_temp -> get_reconplnch_loli( mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 0, &reconpln, &reconch ); //mod view pln axis
		hitcls[i].pln  = reconpln;
		hitcls[i].ch   = reconch;
	}
	fTracking(mod,0,1);
	//recon track along xy axis
	for(int i=0;i<alltrack.size();i++){
		//for(int j=0;j<alltrack[i].hitid.size();j++){
		for(int j=1;j<alltrack[i].hitid.size();j++){
      			vector<Hit>::iterator it;
      			for(it = hit_temp.begin() ; it != hit_temp.end(); it++){
				if(it->id==alltrack[i].hitid[j]){
					//hit_temp.erase(it);
					//it--;
					it->used=1;
				}
			}
		}
	}
	hitcls = hit_temp;
	for(int i=0; i<hitcls.size(); i++){
		int reconpln,reconch;
        	fdim_temp -> get_reconplnch_loli( mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 1, &reconpln, &reconch ); //mod view pln axis
		hitcls[i].pln= reconpln;
		hitcls[i].ch = reconch;
	}
	fTracking(mod,1,1);





	delete fdim_temp;
	hit_temp.clear();

    return true;
};

bool fFindNearestTracks(int* a, int* b, float angle_th, float dist_th){
    /* Function to find two nearest tracks to be connected for water module.
     * Conditions:
     * (1) Angle between two tracks should be less than "angle_th".
     * (2) One track's position should be lower than the other track.
     * (3) Two tracks should be near close each other.
     *     (extending two tracks, the distance between two tracks at the middle point should be less than "dist_th".)
     * If more than one candidate is found, calculate joilik=sqrt((angle/angle_th)^2+(dist/dist_th)^2) for each candidate and select minimum one.
     * If no candidate is found, return false.
     */
    cout << "Find nearest tracks" << endl;
    const int Ntrack=(int)alltrack.size();
    float joilik=1e-6;
    bool candidate=false;
    for(int i=0; i<Ntrack; i++){
        for(int j=0; j<i; j++){
            if(alltrack[i].view!=alltrack[j].view) continue;
            float angle;
            if(fabs(alltrack[i].ang-alltrack[j].ang)<90.0) angle=fabs(alltrack[i].ang-alltrack[j].ang);
            else angle=180.0-fabs(alltrack[i].ang-alltrack[j].ang);
            if(angle>angle_th) continue; // condidtion (1)

            vector<vector<vector<float> > > pos(2, vector<vector<float> >(2, vector<float>(2, 0.0)));;
            pos[0][0][0]=alltrack[i].ixy;
            pos[0][0][1]=alltrack[i].iz;
            pos[0][1][0]=alltrack[i].fxy;
            pos[0][1][1]=alltrack[i].fz;
            pos[1][0][0]=alltrack[j].ixy;
            pos[1][0][1]=alltrack[j].iz;
            pos[1][1][0]=alltrack[j].fxy;
            pos[1][1][1]=alltrack[j].fz;

            /*  Let theta to be mean of two tracks,
             *  if -45 < theta < 45 then fitaxis=0 (fit along z-axis),
             *  else fitaxis=1 (fit along xy-axis).
             */
            int fitaxis=(fabs(alltrack[i].ang-alltrack[j].ang)>90.0 || fabs(alltrack[i].ang+alltrack[j].ang)>90.0);
            cout << "fitaxis " << fitaxis << " angle i " << alltrack[i].ang << " angle j " << alltrack[j].ang << endl;
            if(fitaxis){
                swap(pos[0][0][0], pos[0][0][1]);
                swap(pos[0][1][0], pos[0][1][1]);
                swap(pos[1][0][0], pos[1][0][1]);
                swap(pos[1][1][0], pos[1][1][1]);
            }
            //** swap positions in the order
            if(pos[0][0][1]>pos[0][1][1]) swap(pos[0][0], pos[0][1]);
            if(pos[1][0][1]>pos[1][1][1]) swap(pos[1][0], pos[1][1]);
            if(pos[0][0][1]>pos[1][0][1]) swap(pos[0]   , pos[1]   );
            if(pos[0][1][1]>pos[1][0][1]) continue; // condidion (2)
            //** dist is the distance between two tracks at the middle points.
            float dist=fabs(((pos[0][0][0]-pos[0][1][0])/(pos[0][0][1]-pos[0][1][1])+(pos[1][0][0]-pos[1][1][0])/(pos[1][0][1]-pos[1][1][1]))*(pos[0][1][1]-pos[1][0][1])/2.0-pos[0][1][0]+pos[1][0][0]);
            if(dist>dist_th) continue; // condition (3)

            if(!candidate || joilik>sqrt(angle*angle/angle_th/angle_th+dist*dist/dist_th/dist_th)){
                //** if more than one candidate is found, ...
                *a=i; *b=j;
                joilik=sqrt(angle*angle/angle_th/angle_th+dist*dist/dist_th/dist_th);
                candidate=true;
            }
        }
    }
    return candidate;
};

void fGetHitFromTrack(vector<Hit> &hit, const vector<Track> &track){
    cout << "Get hit from track" << endl;
    vector<Hit>::iterator it;
    for(int i=0; i<(int)track.size(); i++){
        for(int j=0; j<(int)track[i].hitid.size(); j++){
   	   for(it=hitcls_for_joint.begin(); it!=hitcls_for_joint.end(); it++){
                if(it->id==track[i].hitid[j]){
                    hit.push_back(*it);
                    break;
                }
            }
        }
    }
};

void fFitAllHit(const vector<Hit> &hit, vector<float> &par, int axis, bool disp=false, int mod=15){
    cout << "Fit all hit" << endl;
    vector<float> x, y, xe, ye;
    const int Nhits=(int)hit.size();
    INGRID_Dimension* fdim=new INGRID_Dimension();
    for(int i=0; i<Nhits; i++){
        int reconpln, reconch;
        fdim->get_reconplnch_loli(mod, hit[i].view, hit[i].pln, hit[i].ch, axis, &reconpln, &reconch);
        x.push_back (zposi   (mod, hit[i].view, reconpln, reconch, axis));
        y.push_back (xyposi  (mod, hit[i].view, reconpln, reconch, axis));
        xe.push_back(scithick(mod, hit[i].view, reconpln, y.at(i), axis)/sqrt(3));
        ye.push_back(sciwidth(mod, hit[i].view, reconpln, reconch, axis)/sqrt(3));
	std::cout << mod <<  " " << hit[i].view <<  " " << hit[i].pln <<  " " << hit[i].ch << " " << reconpln <<  " " << reconch <<  " " << axis << " " << x[i] << " " << y[i] << std::endl;
    }
    delete fdim;
    TGraphErrors* g=new TGraphErrors(Nhits, &x[0], &y[0], &xe[0], &ye[0]);
    TF1 *f=new TF1("f", "[0]+[1]*x");
    const float deltamin=0.1;
    const float deltax=(fabs(x[Nhits-1]-x[0])<deltamin ? deltamin : (x[Nhits-1]-x[0])); //avoid initial parameter=0
    const float deltay=(fabs(y[Nhits-1]-y[0])<deltamin ? deltamin : (y[Nhits-1]-y[0]));
    f->SetParameters(y[0]-x[0]*deltay/deltax, deltay/deltax);
    cout << "par before fit: " << f->GetParameter(0) << ", " << f->GetParameter(1) << endl;
    g->Fit("f", "WQ"); //ignore errors. Any better solutions?
    //cout << "par after  fit: " << f->GetParameter(0) << ", " << f->GetParameter(1) << " / Chisquare " << f->GetChisquare() << endl;
    cout << "par after  fit: " << f->GetParameter(0) << ", " << f->GetParameter(1) << " / Chisquare " << f->GetChisquare() << " " << 180./3.14*atan(f->GetParameter(1)) << " " << 180./3.14*atan(1./f->GetParameter(1)) << endl;
    if(disp){
        TCanvas* c_allhit=new TCanvas("c_allhit", "c_allhit");
        g->Draw("AP");
        c_allhit->Update();
        cout << "Enter to continue...";
        getchar();
        delete c_allhit;
    }
    par.push_back(f->GetParameter(0));
    par.push_back(f->GetParameter(1));
    f->Delete();
    g->Delete();
};

Track fGetConnectedTrack(const vector<Hit> &hit, const vector<Track> &track, const vector<float> &par, int axis, int mod=15){
    cout << "Get connected track" << endl;
    Track tr;
    tr.clear();
    tr.mod    = mod;
    tr.view   = hit[0].view;
    //tr.intcpt = par[0];
    //tr.slope  = par[1];
    //tr.ang    = atan(par[1])*180/3.141592;

    //tr.clstime //blank
    //tr.vetodist //not used for water module

    for(int i=0; i<(int)track.size(); i++){
        if(track[i].veto) tr.veto=true;
        if(track[i].edge) tr.edge=true;
        if(track[i].stop) tr.stop=true;
    }

    //** intcpt, slope, angle, pln, xy, z **//
    INGRID_Dimension* fdim=new INGRID_Dimension();
    if(axis==0){
        tr.intcpt = par[0];
        tr.slope  = par[1];
        tr.ang    = atan(tr.slope)*180.0/3.141592;

        //** Initialize **//
        int tmppln, tmpch;
        fdim->get_reconplnch_loli(mod, hit[0].view, hit[0].pln, hit[0].ch, 0, &tmppln, &tmpch);
        tr.ipln = hit[0].pln;
        tr.fpln = hit[0].pln;
        tr.iz   = zposi(mod, hit[0].view, tmppln, tmpch);
        tr.fz   = zposi(mod, hit[0].view, tmppln, tmpch);
        tr.ixy  = par[0]+par[1]*tr.iz;
        tr.fxy  = par[0]+par[1]*tr.fz;

        for(int i=0; i<(int)hit.size(); i++){
            int reconpln, reconch;
            fdim->get_reconplnch_loli(mod, hit[i].view, hit[i].pln, hit[i].ch, 0, &reconpln, &reconch);
            if(tr.iz>zposi(mod, hit[i].view, reconpln, reconch)){
                tr.ipln = hit[i].pln;
                tr.iz   = zposi(mod, hit[i].view, reconpln, reconch);
                tr.ixy  = par[0]+par[1]*tr.iz;
            }
            if(tr.fz<zposi(mod, hit[i].view, reconpln, reconch)){
                tr.fpln = hit[i].pln;
                tr.fz   = zposi(mod, hit[i].view, reconpln, reconch);
                tr.fxy  = par[0]+par[1]*tr.fz;
            }
        }
    } //axis==0
    else{ //axis==1
        if(par[1]==0){
            tr.intcpt = -par[0]*1e+8;
            tr.slope  = 1e+8;
            tr.ang    = 90.0;
        }
        else{
            tr.intcpt = -par[0]/par[1];
            tr.slope  = 1.0/par[1];
            tr.ang    = atan(tr.slope)*180.0/3.141592;
        }

        //** Initialize **//
        int tmppln, tmpch;
        fdim->get_reconplnch_loli(mod, hit[0].view, hit[0].pln, hit[0].ch, 0, &tmppln, &tmpch);
        tr.ipln = hit[0].pln;
        tr.fpln = hit[0].pln;
        tr.ixy  = xyposi(mod, hit[0].view, tmppln, tmpch);
        tr.fxy  = xyposi(mod, hit[0].view, tmppln, tmpch);
        tr.iz   = par[0]+par[1]*tr.ixy;
        tr.fz   = par[0]+par[1]*tr.fxy;

        for(int i=0; i<(int)hit.size(); i++){
            int reconpln, reconch;
            fdim->get_reconplnch_loli(mod, hit[i].view, hit[i].pln, hit[i].ch, 0, &reconpln, &reconch);
            if(par[1]>0 != tr.ixy<xyposi(mod, hit[i].view, reconpln, reconch)){
                tr.ipln = hit[i].pln;
                tr.ixy  = xyposi(mod, hit[i].view, reconpln, reconch);
                tr.iz   = par[0]+par[1]*tr.ixy;
            }
            if(par[1]>0 != tr.fxy>xyposi(mod, hit[i].view, reconpln, reconch)){
                tr.fpln = hit[i].pln;
                tr.fxy  = xyposi(mod, hit[i].view, reconpln, reconch);
                tr.fz   = par[0]+par[1]*tr.fxy;
            }
        }
    } //axis==1
    delete fdim;

    //** hitid, isohit **//
    for(int i=0; i<(int)hit.size(); i++){
        tr.hitid.push_back(hit[i].id);
        bool isohit=true;
        for(int j=0; j<(int)alltrack.size(); j++){ //Here alltrack does not include tracks to be connected.
            for(int k=0; k<(int)alltrack[j].hitid.size(); k++){
                if(hit[i].id==alltrack[j].hitid[k]){
                    isohit=false;
                    break;
                }
            }
            if(!isohit) break;
        }
        tr.isohit.push_back(isohit);
    }

    return tr;
};

bool fConnectTracks(float angle_th, float dist_th){ //hosomi 160622
    /* Connect break two tracks due to gap between layers in water module.
     * Procedure:
     * - Find nearest tracks
     * - Get hit info from two tracks to be connected
     * - Fit all hits to be connected
     * - Get connected trak info
     * If no candidate is found, return false
     */
    int a, b;
    if(!fFindNearestTracks(&a, &b, angle_th, dist_th)) return false;
    for(int i=0; i<(int)alltrack[a].hitid.size(); i++) cout << "track a hitid: " << alltrack[a].hitid[i] << " isohit: " << alltrack[a].isohit[i] << endl;
    for(int i=0; i<(int)alltrack[b].hitid.size(); i++) cout << "track b hitid: " << alltrack[b].hitid[i] << " isohit: " << alltrack[b].isohit[i] <<endl;

    vector<Hit> hit_connect;
    vector<Track> track_connect;
    track_connect.push_back(alltrack[a]);
    track_connect.push_back(alltrack[b]);
    fGetHitFromTrack(hit_connect, track_connect);
    if(hit_connect.size()<2) return false;

    /*  Let theta to be mean of two tracks,
     *  if -45 < theta < 45 then fitaxis=0 (fit along z-axis),
     *  else fitaxis=1 (fit along xy-axis).
     */
    int fitaxis=(fabs(alltrack[a].ang-alltrack[b].ang)>90.0 || fabs(alltrack[a].ang+alltrack[b].ang)>90.0);
    cout << "fitaxis " << fitaxis << " angle a " << alltrack[a].ang << " angle b " << alltrack[b].ang << endl;
    vector<float> par;
    fFitAllHit(hit_connect, par, fitaxis);

    alltrack.erase(alltrack.begin()+a);
    alltrack.erase(alltrack.begin()+b);
    alltrack.push_back(fGetConnectedTrack(hit_connect, track_connect, par, fitaxis));

    return true;
};

void fGetTrackInfo(char* buf, int ievt, int Nconnect){
    int ntracks=(int)alltrack.size();
    int nview=0;
    for(int i=0; i<ntracks; i++) nview += alltrack[i].view;
    //char buf[64];
    sprintf(buf, "%d\t%d\t%d\t%d", ievt, ntracks, nview, Nconnect);
}



#endif

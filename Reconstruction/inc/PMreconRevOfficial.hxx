#ifndef __PMRECON_HXX__
#define __PMRECON_HXX__

//#define PMEDGE

// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
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



//__________________________________________________________


class Hit{
public:
  Int_t       id;   //b/w Hit class and IngridHitSummary
  Int_t      mod;
  Int_t      pln;
  Int_t     view;
  Int_t       ch;
  //Int_t      adc;
  Float_t     pe;   //p.e.
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
vector<Hit> hitcls;


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


int chmax(int mod){
  int maxch;
  if(mod!=16)maxch=24;
  else maxch=32;
  return maxch;
}

int plnmax(int mod){
  int maxpln;
  if(mod!=16)maxpln=11;
  else maxpln=18;
  return maxpln;
}


int zposi(int mod,int view,int pln){
  int posiz;
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
  else{
    if(view==0){
      posiz=5+105*pln+1074.5+10;// ML 2017/06/27 -  INGRID modules start with vertical bars
    }
    else{
      posiz=5+105*pln+1074.5;// ML 2017/06/27
    }
  }  
  return posiz;
}

float xyposi(int mod,int pln,int ch){
  float posixy;

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
  else{
    posixy=ch*50+25;
  }
  return posixy;
}

float vposiz(int mod,int ch){
  float X;
  if(mod!=16)
    //X=9.5+50*ch+1020;
    X=5+9.5+50*ch+1075;
  else
    X=5+9.5+50*ch;
  return X;
}


int vposixy(int mod,int pln){
  int Y=0;
  if(mod==0){
    if(pln==11)Y=-105.75;  //right
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
    if(pln==11)Y=-105.75;  //right
    else if(pln==12)Y=1309;//left
    else if(pln==13)Y=-59; //botom
    else Y=1284;           //top
  }
  else if(mod<14){
    if(pln==11)Y=-105.75;  //right
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

float sciwidth(int mod,int pln,int ch){
  float sciw;

  if(mod!=16||pln==0||ch<8||ch>=24)sciw=25;
  else sciw=12.5;

  return sciw;
};

float scithick(int mod,int pln,float xy){
  float scith;

  if(mod!=16||pln==0||xy<400||xy>800)scith=5;
  else scith=6.5;

  return scith;
};


float Yerr(int mod,int pln,int ch,int edge){
  float yerr;
  if(edge==0)yerr=xyposi(mod,pln,ch)-sciwidth(mod,pln,ch);
  else       yerr=xyposi(mod,pln,ch)+sciwidth(mod,pln,ch);
  return yerr;
};






bool actpln[2][17];
int nactpln;
int fNactpln(int mod){
  nactpln=0;
  memset(actpln,false,sizeof(actpln));
  for(int i=0; i<hitcls.size(); i++){
    if(hitcls[i].pln<plnmax(mod))actpln[hitcls[i].view][hitcls[i].pln]=true;
  }
  for(int i=0; i<plnmax(mod); i++){
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
  for(int i=0; i<hitcls.size(); i++){
    if(hitcls[i].pln<plnmax(mod)){
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

  if(nhit <= maxhit)return false;
  for(Int_t i=0 ; i<nhit-maxhit; i++){
    //if( fabs( hit[i].time - hit[i+maxhit-1].time) < cTdcRsr ){
    if( fabs( hit[i].time - hit[i+maxhit].time) < cTdcRsr ){

      //long basetime = 0;
      //for( int j=0; j<maxhit; j++){
      //for( int j=0; j<maxhit+1; j++){
	//basetime += hit[i+j].time;
      //}
      //basetime = basetime/(maxhit);
      //basetime = basetime/(maxhit+1);
     
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
      hitcls.clear(); //Reset hit clster
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

int view,id[2][602],pln[2][602],ch[2][602],hit[2],hitnum[2][18][32];
float pe[2][602];
int PLN,PLN2,CH,CL,CL2,CLHIT,CLHIT2,CLHIT3,CELL,CELL2,NEI,TRA,TRA2,TRACELL,TRACL,TRACL2,HIT,DIST,DIST2,DIST3,dummy,TMP;
bool hitcl[2][18][16];
int clchi[2][18][16],clchf[2][18][16],ncl[2][18],clcenter[2][18][16],numcl[2][18][16],clhit[2][18][16][10],clpe[2][18][16];
int cellu[2][17][500][2],celld[2][17][500][2],ncell[2][17][2],value[2][17][500][2],nvalue[2][17][500][2];
bool neibor[2][16][500][2][2];
int neiu[2][16][500][2][2],neid[2][16][500][2][2],nnei[2][16][2][2];
int track_cell[2][500][17],track_pln[2][500][17],track_dist[2][500][17],ntracell[2][500],ntracl[2][500],ntracl2[2][500],ntrack[2],ntrack2[2],ntrack3[2],trank[2][19][50],ncltra[2][19],rank[2][500];
bool ttrack[2][500],ttrack2[2][500];
int plane[2][500][18],clus[2][500][18];
bool vetowtracking[2][500],edgewtracking[2][500],track_stop[2][500];
float vetodist[2][500];
float ang[2];
float hosei[2][500],dedx[2][500],cluspe[2][18][16],trackpe[2][500];//add
int pdg[2][1000],clpdg[2][18][16][4],trackpdg[2][500][4],trapdg[2][50],len[2][50];//add
bool ovcl[2][18][16];//add
float dis;
float XX[2][500][18],YY[2][500][18],Xe[2][500][18],Yel[2][500][18],Yeh[2][500][18];
float x[3],y[3],xe[3],yel[3],yeh[3];
float chi2,tmp;
float par[2][500][2];
float trape[2][500];
float ulimit;
bool upst,dwst;
//int select;
int nhit;
int ue,shita;
//bool actpln[8][2][18];
//float peactpln[8];
//int numactpln[8][2];
//bool hitcyc[8];
//int length[NMOD][2];
//int ncychit[8];
//int exptime[8];
//float timecl[8];
//float npecychit[8];
int star[2][50],fini[2][50];
float sep[2][50],kat[2][50];//,totpe[2][50];
bool shar[2][50];
int ntra[2],mtra[2];
//int hitpln[2][50][18];
//float hitpe[2][50][18];
int pid;
int hit_n[2][500],hit_id[2][500][200];
bool isohit[2][500][200];

bool fTracking(int mod, int VetoUpstreamCriteria, double VetoEdgeCriteria, double FVCriteria){

  memset(hit,0,sizeof(hit));
  memset(hitnum,-1,sizeof(hitnum));
  alltrack.clear();
  for(int i=0; i<hitcls.size(); i++){
    view=hitcls[i].view;
    id[view][hit[view]]=hitcls[i].id;
    pln[view][hit[view]]=hitcls[i].pln;
    ch[view][hit[view]]=hitcls[i].ch;
    pe[view][hit[view]]=hitcls[i].pe;
    //pdg[view][hit[view]]=hitcls[i].pdg;
    if(pln[view][hit[view]]<plnmax(mod))hitnum[view][pln[view][hit[view]]][ch[view][hit[view]]]=hit[view];
    hit[view]++;    
  }

  for(int VIEW=0;VIEW<2;VIEW++){

    //*****Define clusters*****
    memset(ncl[VIEW],0,sizeof(ncl[VIEW]));
    memset(numcl[VIEW],0,sizeof(numcl[VIEW]));
    memset(clhit[VIEW],0,sizeof(clhit[VIEW]));
    for(PLN=0;PLN<plnmax(mod);PLN++){
      CH=0;
      while(1){
	if(hitnum[VIEW][PLN][CH]>=0){
	  clchi[VIEW][PLN][ncl[VIEW][PLN]]=CH;
	  clhit[VIEW][PLN][ncl[VIEW][PLN]][numcl[VIEW][PLN][ncl[VIEW][PLN]]]=hitnum[VIEW][PLN][CH];
	  while(1){
	    numcl[VIEW][PLN][ncl[VIEW][PLN]]++;
	    if(CH==chmax(mod)-1)break;
	    if(hitnum[VIEW][PLN][CH+1]<0)break;
	    CH++;
	    clhit[VIEW][PLN][ncl[VIEW][PLN]][numcl[VIEW][PLN][ncl[VIEW][PLN]]]=hitnum[VIEW][PLN][CH];
	  }
	  clchf[VIEW][PLN][ncl[VIEW][PLN]]=CH;
	  ncl[VIEW][PLN]++;
	}//if
	if(CH==chmax(mod)-1)break;
	CH++;
      }//while
    }//for PLN
	
    //*****Fill pe and center of clusters*****
    memset(clpe[VIEW],0,sizeof(clpe[VIEW]));
    memset(clcenter[VIEW],0,sizeof(clcenter[VIEW]));
    for(PLN=0;PLN<plnmax(mod);PLN++){
      for(CL=0;CL<ncl[VIEW][PLN];CL++){
	for(CLHIT=0;CLHIT<numcl[VIEW][PLN][CL];CLHIT++){
	  clpe[VIEW][PLN][CL]+=pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]];
	  clcenter[VIEW][PLN][CL]+=pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]]*xyposi(mod,PLN,clchi[VIEW][PLN][CL]+CLHIT);
	}
	clcenter[VIEW][PLN][CL]=clcenter[VIEW][PLN][CL]/clpe[VIEW][PLN][CL];
      }
    }

    
    //*****Define cells*****
    memset(ncell[VIEW],0,sizeof(ncell[VIEW]));
    for(PLN=0;PLN<plnmax(mod)-1;PLN++){
      for(CL=0;CL<ncl[VIEW][PLN];CL++){
	for(DIST=0;DIST<2;DIST++){
	  if(PLN==plnmax(mod)-2&&DIST==1)continue;
	  for(CL2=0;CL2<ncl[VIEW][PLN+DIST+1];CL2++){
	    //if(DIST==1&&fabs(clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+2][CL2])>92)continue;
	    if(DIST==1&&fabs(clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+2][CL2])>150)continue;
	    cellu[VIEW][PLN][ncell[VIEW][PLN][DIST]][DIST]=CL;
	    celld[VIEW][PLN][ncell[VIEW][PLN][DIST]][DIST]=CL2;
	    ncell[VIEW][PLN][DIST]++;
	  }
	}
      }
    }
    

    //*****Define neiborhoods*****
    memset(nnei[VIEW],0,sizeof(nnei[VIEW]));
    for(PLN=0;PLN<plnmax(mod)-2;PLN++){
      for(DIST=0;DIST<2;DIST++){
	if(PLN==0&&DIST==1)continue;
	for(CELL=0;CELL<ncell[VIEW][PLN-DIST][DIST];CELL++){
	  for(DIST2=0;DIST2<2;DIST2++){
	    if(PLN==plnmax(mod)-3&&DIST2==1)continue;
	    if(DIST==1&&DIST2==1)continue;
	    for(CELL2=0;CELL2<ncell[VIEW][PLN+1][DIST2];CELL2++){
	      if(celld[VIEW][PLN-DIST][CELL][DIST]==cellu[VIEW][PLN+1][CELL2][DIST2]){
		
		x[0]=zposi(mod,VIEW,PLN-DIST);
		x[1]=zposi(mod,VIEW,PLN+1);
		x[2]=zposi(mod,VIEW,PLN+2+DIST2);
		
		y[0]=clcenter[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]];
		y[1]=clcenter[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]];
		y[2]=clcenter[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]];

		xe[0]=scithick(mod,PLN-DIST,y[0]);
		xe[1]=scithick(mod,PLN+1,y[1]);
		xe[2]=scithick(mod,PLN+2+DIST2,y[2]);

		yel[0]=y[0]-Yerr(mod,PLN-DIST,clchi[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]],0);
		yel[1]=y[1]-Yerr(mod,PLN+1,clchi[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]],0);
		yel[2]=y[2]-Yerr(mod,PLN+2+DIST2,clchi[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]],0);
		
		yeh[0]=-y[0]+Yerr(mod,PLN-DIST,clchf[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]],1);
		yeh[1]=-y[1]+Yerr(mod,PLN+1,clchf[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]],1);
		yeh[2]=-y[2]+Yerr(mod,PLN+2+DIST2,clchf[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]],1);
		
		TGraphAsymmErrors *graph=new TGraphAsymmErrors(3,x,y,xe,xe,yel,yeh);
		TF1 *f=new TF1("f","[0]+[1]*x");
		f->SetParameters(y[0]-x[0]*(y[2]-y[0])/(x[2]-x[0]),(y[2]-y[0])/(x[2]-x[0]));
		
		graph->Fit("f","Q");
		chi2=f->GetChisquare();
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
		
		if(DIST==0&&DIST2==0)ulimit=3;
		else                 ulimit=1.5;

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
    for(int jndex=0;jndex<18;jndex++){
      for(PLN=0;PLN<plnmax(mod)-2;PLN++){
	for(DIST=0;DIST<2;DIST++){
	  for(DIST2=0;DIST2<2;DIST2++){
	    if(PLN==0&&DIST==1)continue;
	    if(PLN==plnmax(mod)-3&&DIST2==1)continue;
	    if(DIST==1&&DIST2==1)continue;	    
	    for(NEI=0;NEI<nnei[VIEW][PLN][DIST][DIST2];NEI++){
	      if(value[VIEW][PLN-DIST][neiu[VIEW][PLN][NEI][DIST][DIST2]][DIST]==value[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]){
		nvalue[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]=value[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]+1;
	      }
	    }
	  } 
	}
      }
      for(PLN=0;PLN<plnmax(mod)-1;PLN++){	   
	for(DIST=0;DIST<2;DIST++){
	  if(PLN==plnmax(mod)-2&&DIST==1)continue;
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
    for(PLN=1;PLN<plnmax(mod)-1;PLN++){
      for(DIST=0;DIST<2;DIST++){
	if(PLN==plnmax(mod)-2&&DIST==1)continue;
	for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
	  if(value[VIEW][PLN][CELL][DIST]>0){
	    if(PLN+DIST==plnmax(mod)-2){upst=true;}
	    else{
	      upst=true;
	      for(DIST2=0;DIST2<2;DIST2++){
		if(DIST==1&&DIST2==1)continue;
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
	      XX[VIEW][ntrack[VIEW]][0]=zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][0]);
	      YY[VIEW][ntrack[VIEW]][0]=clcenter[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]];
	      trape[VIEW][ntrack[VIEW]]+=clpe[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]];

	      Yel[VIEW][ntrack[VIEW]][0]=YY[VIEW][ntrack[VIEW]][0]-Yerr(mod,plane[VIEW][ntrack[VIEW]][0],clchi[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]],0);
	      Yeh[VIEW][ntrack[VIEW]][0]=-YY[VIEW][ntrack[VIEW]][0]+Yerr(mod,plane[VIEW][ntrack[VIEW]][0],clchf[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]],1);
	      Xe[VIEW][ntrack[VIEW]][0]=scithick(mod,plane[VIEW][ntrack[VIEW]][0],YY[VIEW][ntrack[VIEW]][0]);
	      
	      plane[VIEW][ntrack[VIEW]][1]=track_pln[VIEW][ntrack[VIEW]][0];
	      clus[VIEW][ntrack[VIEW]][1]=cellu[VIEW][track_pln[VIEW][ntrack[VIEW]][0]][track_cell[VIEW][ntrack[VIEW]][0]][track_dist[VIEW][ntrack[VIEW]][0]];
	      XX[VIEW][ntrack[VIEW]][1]=zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][1]);
	      YY[VIEW][ntrack[VIEW]][1]=clcenter[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]];
	      trape[VIEW][ntrack[VIEW]]+=clpe[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]];
	      
	      Yel[VIEW][ntrack[VIEW]][1]=YY[VIEW][ntrack[VIEW]][1]-Yerr(mod,plane[VIEW][ntrack[VIEW]][1],clchi[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]],0);
	      Yeh[VIEW][ntrack[VIEW]][1]=-YY[VIEW][ntrack[VIEW]][1]+Yerr(mod,plane[VIEW][ntrack[VIEW]][1],clchf[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]],1);
	      Xe[VIEW][ntrack[VIEW]][1]=scithick(mod,plane[VIEW][ntrack[VIEW]][1],YY[VIEW][ntrack[VIEW]][1]);


	      ntracell[VIEW][ntrack[VIEW]]=1;
	      PLN2=PLN-1;
	      int chi2old=3;
	      while(PLN2>=0){
		dwst=true;
		DIST3=track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]-1];
		for(DIST2=0;DIST2<2;DIST2++){
		  if(DIST3==1&&DIST2==1)continue;
		  if(DIST2==1&&PLN2==0)continue;
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
			    Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-Yerr(mod,PLN2-DIST2,clchi[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]],0);
			    Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=-YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]+Yerr(mod,PLN2-DIST2,clchf[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]],1);
			    Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=scithick(mod,PLN2-DIST2,YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]);


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
			    
			    
			    if((tmp>chi2||upstpe2<4.5)&&upstpe1>=4.5){
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
			  XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]);			  
			  YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=clcenter[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]];

			  Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-Yerr(mod,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],clchi[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]],0);
			  Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=-YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]+Yerr(mod,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],clchf[VIEW][plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]][clus[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]],1);
			  Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=scithick(mod,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]);
			  
			  //TGraph *graph1=new TGraph((ntracell[VIEW][ntrack[VIEW]]+2),XX[VIEW][ntrack[VIEW]],YY[VIEW][ntrack[VIEW]]);
			  TGraphAsymmErrors *graph1=new TGraphAsymmErrors((ntracell[VIEW][ntrack[VIEW]]+2),XX[VIEW][ntrack[VIEW]],YY[VIEW][ntrack[VIEW]],Xe[VIEW][ntrack[VIEW]],Xe[VIEW][ntrack[VIEW]],Yel[VIEW][ntrack[VIEW]],Yeh[VIEW][ntrack[VIEW]]);
			  TF1 *f1=new TF1("f1","[0]+[1]*x");
			  f1->SetParameters(YY[VIEW][ntrack[VIEW]][0]-XX[VIEW][ntrack[VIEW]][0]*(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])/(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]),(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])/(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]));
			  
			  graph1->Fit("f1","Q");
			  chi2=f1->GetChisquare();
			  graph1->Delete();
			  f1->Delete();


			  if(
			     (mod!=16&&((chi2/ntracell[VIEW][ntrack[VIEW]]>3&&ntracell[VIEW][ntrack[VIEW]]==1)||
					(chi2/ntracell[VIEW][ntrack[VIEW]]>2&&ntracell[VIEW][ntrack[VIEW]]>1)||
					((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 -(ntracell[VIEW][ntrack[VIEW]]-1)*0.02)))
			     ||
			     (mod==16&&((chi2/ntracell[VIEW][ntrack[VIEW]]>3&&ntracell[VIEW][ntrack[VIEW]]==1)||
					(chi2/ntracell[VIEW][ntrack[VIEW]]>2&&ntracell[VIEW][ntrack[VIEW]]>1)||
					((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 -(ntracell[VIEW][ntrack[VIEW]]-1)*0.015)))
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
	      if(ntrack[VIEW]==500) break;//ML 2017/05/27
	    }
	  }
	}// for Cell
	if(ntrack[VIEW]==500) break;//ML 2017/05/27
      }// for Dist
      if(ntrack[VIEW]==500) break;//ML 2017/05/27
    }//for Pln


    
    //*****Track rank*****
    memset(trank[VIEW],0,sizeof(trank[VIEW]));
    memset(rank[VIEW],0,sizeof(rank[VIEW]));
    memset(ncltra[VIEW],0,sizeof(ncltra[VIEW]));
    ntrack2[VIEW]=0;
    for(TRACL=18;TRACL>0;TRACL--){
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
    
    //cout<<mod<<" "<<ntrack2[VIEW]<<endl;

    //*****True track selection*****
    memset(ttrack[VIEW],true,sizeof(ttrack[VIEW]));
    memset(hitcl[VIEW],false,sizeof(hitcl[VIEW]));
    memset(ntracl2[VIEW],0,sizeof(ntracl2[VIEW]));
    memset(ovcl[VIEW],false,sizeof(ovcl[VIEW]));//add
    ntrack2[VIEW]=0;
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      TRA2=rank[VIEW][TRA];
      for(TRACL=0;TRACL<ntracl[VIEW][TRA2];TRACL++){
	//cout<<allhit.size()<<" "<<VIEW<<" "<<TRA2<<" "<<TRACL<<" "<<plane[VIEW][TRA2][TRACL]<<" "<<clus[VIEW][TRA2][TRACL]<<endl;
	if(!hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]
	   &&clpe[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]>4.5){
	  ntracl2[VIEW][TRA2]++;
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

      if(ntracl2[VIEW][TRA2]==2){
	if(ntracl[VIEW][TRA2]<4){
	  if(!hitcl[VIEW][plane[VIEW][TRA2][0]][clus[VIEW][TRA2][0]]
	     &&clpe[VIEW][plane[VIEW][TRA2][0]][clus[VIEW][TRA2][0]]>4.5
	     &&!hitcl[VIEW][plane[VIEW][TRA2][ntracl[VIEW][TRA2]-1]][clus[VIEW][TRA2][ntracl[VIEW][TRA2]-1]]
	     &&clpe[VIEW][plane[VIEW][TRA2][ntracl[VIEW][TRA2]-1]][clus[VIEW][TRA2][ntracl[VIEW][TRA2]-1]]>4.5)
	    ntracl2[VIEW][TRA2]=1;
	}
      }

      if(ntracl2[VIEW][TRA2]>1){
	for(TRACL=0;TRACL<ntracl[VIEW][TRA2];TRACL++){
	  if(hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]==true)ovcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]=true;//add
	  hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]=true;
	}
	ntrack2[VIEW]++;
      }
      else ttrack[VIEW][TRA2]=false;
    }



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
    for(int TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]==false)continue;
      //if(plane[VIEW][TRA][ntracl[VIEW][TRA]-1]==0){vetowtracking[VIEW][TRA]=true;}
      if((mod!=16&&plane[VIEW][TRA][ntracl[VIEW][TRA]-1]<=VetoUpstreamCriteria)||(mod==16&&plane[VIEW][TRA][ntracl[VIEW][TRA]-1]<=VetoUpstreamCriteria)){vetowtracking[VIEW][TRA]=true;}
      if(plane[VIEW][TRA][0]==plnmax(mod)-1){track_stop[VIEW][TRA]=false;}

#ifdef PMEDGE
      for(int TRACL=ntracl[VIEW][TRA]-2;TRACL<ntracl[VIEW][TRA];TRACL++){
	for(CLHIT=0;CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]];CLHIT++){
	  if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==0){edgewtracking[VIEW][TRA]=true;}
	  else if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==chmax(mod)-1){edgewtracking[VIEW][TRA]=true;}
	}
      }

      for(int TRACL=0;TRACL<2;TRACL++){
	for(CLHIT=0;CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]];CLHIT++){
	  if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==0){track_stop[VIEW][TRA]=false;}
	  else if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==chmax(mod)-1){track_stop[VIEW][TRA]=false;}
	}
      }
#else
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]]))+par[VIEW][TRA][0]<FVCriteria){edgewtracking[VIEW][TRA]=true;}
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]]))+par[VIEW][TRA][0]>1200-FVCriteria){edgewtracking[VIEW][TRA]=true;}

      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0]))+par[VIEW][TRA][0]<100){track_stop[VIEW][TRA]=false;}
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0]))+par[VIEW][TRA][0]>1100){track_stop[VIEW][TRA]=false;}

#endif
	

    }


    //*****side veto*****
    // loop over all hits in the side veto planes to see if one could be the beginning of the track
    for(int TRA=0;TRA<ntrack[VIEW];TRA++){
      vetodist[VIEW][TRA]=1e5;
      if(ttrack[VIEW][TRA]==false)continue;
      for(HIT=0;HIT<hit[VIEW];HIT++){
	if(pln[VIEW][HIT]>=plnmax(mod)){
	  // cout<<"side veto is hit, mod="<<mod<<" plane="<<pln[VIEW][HIT]<<endl;
	  //if >=2 layers (46*2 in PM or 107*2 in I) between the track starting point 
	  // and what would be the entering point of the track on the side of the module,
	  // the track really starts inside the module -> no need to look for hits in side vetos
	  if(par[VIEW][TRA][1]>0){
	    if(vposixy(mod,pln[VIEW][HIT])>0)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	  }
	  else if(par[VIEW][TRA][1]<0){
	    if(vposixy(mod,pln[VIEW][HIT])<0)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	  }
	  else continue;
	  // dis is the distance of the hit to the track (orthogonal projection)
	  dis=fabs(par[VIEW][TRA][1]*vposiz(mod,ch[VIEW][HIT])-vposixy(mod,pln[VIEW][HIT])+par[VIEW][TRA][0])/sqrt((par[VIEW][TRA][1])*(par[VIEW][TRA][1])+1);
	  if(dis<VetoEdgeCriteria)vetowtracking[VIEW][TRA]=true;
	  if(dis<vetodist[VIEW][TRA])vetodist[VIEW][TRA]=dis;
	}
      }
    }


    //*****mod fc with veto*****
    for(int TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]==false)continue;
      for(HIT=0;HIT<hit[VIEW];HIT++){
	if(pln[VIEW][HIT]>=plnmax(mod)){
	  if(par[VIEW][TRA][1]>0){
	    if(vposixy(mod,pln[VIEW][HIT])<0)continue;
	    if(XX[VIEW][TRA][0]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][0]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	  }
	  else if(par[VIEW][TRA][1]<0){
	    if(vposixy(mod,pln[VIEW][HIT])>0)continue;
	    if(XX[VIEW][TRA][0]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2&&mod==16)continue;
	    if(XX[VIEW][TRA][0]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2&&mod!=16)continue;
	  }
	  else continue;
	  dis=fabs(par[VIEW][TRA][1]*vposiz(mod,ch[VIEW][HIT])-vposixy(mod,pln[VIEW][HIT])+par[VIEW][TRA][0])/sqrt((par[VIEW][TRA][1])*(par[VIEW][TRA][1])+1);
	  if(dis<VetoEdgeCriteria)track_stop[VIEW][TRA]=false;
	}
      }
    }

    //*****cut low upstream hit****
    for(int TRA=0;TRA<ntrack[VIEW];TRA++){
      if(vetowtracking[VIEW][TRA]||edgewtracking[VIEW][TRA])continue;//new
      if(ttrack[VIEW][TRA]&&clpe[VIEW][plane[VIEW][TRA][ntracl[VIEW][TRA]-1]][clus[VIEW][TRA][ntracl[VIEW][TRA]-1]]<4.5){
	if(ntracl[VIEW][TRA]==3){
	  if(ovcl[VIEW][plane[VIEW][TRA][0]][clus[VIEW][TRA][0]]||
	     ovcl[VIEW][plane[VIEW][TRA][1]][clus[VIEW][TRA][1]]||
	     ovcl[VIEW][plane[VIEW][TRA][2]][clus[VIEW][TRA][2]]
	     ){
	    ttrack[VIEW][TRA]=false;
	    ntrack2[VIEW]--;
	    continue;
	  }
	  else continue;
	}
	ntracl[VIEW][TRA]--;
	ntracell[VIEW][TRA]--;
      }
    }



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
	for(int TRACL=0;TRACL<ntracl[VIEW][TRA];TRACL++){

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




    //*****reconstructed track info****
    Track track;
    for(int TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA]){
	track.view = VIEW;
	track.fpln = plane[VIEW][TRA][0];
	track.ipln = plane[VIEW][TRA][ntracell[VIEW][TRA]];
	track.fxy  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0]))+par[VIEW][TRA][0];
	track.ixy  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]]))+par[VIEW][TRA][0];
	track.fz   = zposi(mod,VIEW,plane[VIEW][TRA][0]);
	track.iz   = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]]);
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
	alltrack.push_back(track);
      }
    }
  }

  if(alltrack.size()>0) return true;
  else                  return false;

};


#endif

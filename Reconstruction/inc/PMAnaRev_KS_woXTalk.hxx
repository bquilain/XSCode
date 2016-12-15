#ifndef __PMANA_HXX__
#define __PMANA_HXX__

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
#include "IngridHitSummary.h"
#include "PMrecon.hxx"

//__________________________________________________________

int pln_th  =2;
float ch_th =150;
int diff_th=4;
float ang_th=35;
float pos_th=85;

double TTCL;//for MuCL
int   nCL;  //for MuCL


class Hits{
public:

  Int_t    mod;
  Int_t    view;
  Int_t    pln;
  Int_t    ch;
  Int_t    pdg;
  Float_t  pe;
  Float_t  lope;
  Bool_t   isohit;
  Float_t  xy;
  Float_t  z;
  Float_t  time;

  void clear(){

    mod    = -1;
    view   = -1;
    pln    = -1;
    ch     = -1;
    pdg    = -1;
    pe     = -1.e-5;
    lope   = -1.e-5;
    isohit = false;
    xy = -1e-5;
    z = -1e-5;
    time = -1e-5;
  }
};



class TrackIng{
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
  vector<Hits> hit;  
  //Float_t totpe;
  //Int_t isohit;
  //vector<Int_t> hitid;
  //vector<Bool_t> isohit;
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
    hit.clear();
    //hitid.clear();
    //isohit.clear();
    //totpe   =  -1.e-5;
    //isohit   =  -1;
    //memset(pdg,-1,sizeof(pdg));
  }
};



class TrackPM{
public:
  
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
  Int_t    ing_imod;
  Int_t    ing_fmod;
  Int_t    ing_ipln;
  Int_t    ing_fpln;
  Bool_t   ing_trk;
  Bool_t   pm_stop;
  Bool_t   ing_stop;
  Int_t    iron_pene;
  vector<Hits> hit;  
  
  void clear(){
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
    ing_imod = -1;
    ing_fmod = -1;
    ing_ipln = -1;
    ing_fpln = -1;
    ing_trk = false;
    pm_stop = false;
    ing_stop = false;
    iron_pene = -1;
    hit.clear();
  }
};


class Trk{
public:

  Float_t  x;
  Float_t  y;
  Float_t  z;
  Float_t  zx;
  Float_t  zy;

  Int_t    startxpln;
  Int_t    startypln;
  Float_t    startxch;
  Float_t    startych;
  Int_t    endxpln;
  Int_t    endypln;
  Float_t    endxch;
  Float_t    endych;
  Float_t  thetax;
  Float_t  thetay;
  Float_t  angle;
  Int_t    ing_startmod;
  Int_t    ing_endmod;
  Int_t    ing_startpln;
  Int_t    ing_endpln;
  Bool_t   ing_trk;
  Bool_t   pm_stop;
  Bool_t   ing_stop;
  Float_t  sci_range;
  Float_t  iron_range;
  Int_t    iron_pene;
  Bool_t   vetowtracking; // Upstream VETO
  Bool_t   edgewtracking; // Fiducial CUT
  Int_t    pdg;
  Float_t  trkpe;
  Float_t  slopex;
  Float_t  slopey;
  Float_t  intcptx;
  Float_t  intcpty;
  Float_t    mucl;
  vector<Hits> hit;

  void clear(){

    x=  -1.e-5;
    y=  -1.e-5;
    z=  -1.e-5;
    zx=  -1.e-5;
    zy=  -1.e-5;
    startxpln= -1;
    startypln= -1;
    startxch=  -1.e-5;
    startych=  -1.e-5;
    endxpln= -1;
    endypln= -1;
    endxch=  -1.e-5;
    endych=  -1.e-5;
    thetax=  -1.e-5;
    thetay=  -1.e-5;
    angle=  -1.e-5;
    ing_startmod= -1;
    ing_endmod= -1;
    ing_startpln= -1;
    ing_endpln= -1;
    ing_trk= false;
    pm_stop= false;
    ing_stop= false;
    sci_range=  -1.e-5;
    iron_range=  -1.e-5;
    iron_pene= -1;
    vetowtracking=false; // Upstream VETO
    edgewtracking=false; // Fiducial CUT
    pdg = -1;
    trkpe =  -1.e-5;
    slopex =  -1.e-5;
    slopey =  -1.e-5;
    intcptx =  -1.e-5;
    intcpty =  -1.e-5;
    mucl    =  -1.e-5;
    hit.clear();
  }
};




class PMTrack{
public:

  Int_t Ntrack;
  Int_t Ningtrack;
  Float_t  clstime;
  Bool_t   vetowtracking; // Upstream VETO
  Bool_t   edgewtracking; // Fiducial CUT

  vector<Trk>    trk;

  void clear(){
    Ntrack = -1;
    Ningtrack = -1;
    clstime = -1.e-5;
    vetowtracking=false; // Upstream VETO
    edgewtracking=false; // Fiducial CUT
    trk.clear();
  }
};

vector<PMTrack> pmtrack;
vector<TrackPM> htrack;
vector<TrackPM> vtrack;
vector<TrackIng> hingtrack;
vector<TrackIng> vingtrack;


bool withend(const TrackPM& left, const TrackPM& right){
  if(left.fpln != right.fpln)
    return left.fpln > right.fpln;
  else if(left.ipln != right.ipln)
    return left.ipln < right.ipln;
  else
    return fabs(left.ang) < fabs(right.ang);
};

bool withcenter(const TrackIng& left, const TrackIng& right){
  if(abs(left.mod-3) != abs(right.mod-3))
    return abs(left.mod-3) < abs(right.mod-3);
  else if(left.ipln != right.ipln)
    return left.ipln < right.ipln;
  else if(left.fpln != right.fpln)
    return left.fpln > right.fpln;
  else
    return fabs(left.ang) < fabs(right.ang);
};

bool fSortTrack(vector<TrackPM> &a){
  std::stable_sort(a.begin(), a.end(), withend);
};

bool fIngSortTrack(vector<TrackIng> &a){
  std::stable_sort(a.begin(), a.end(), withcenter);
};

bool fIngPMJoint(vector<TrackIng> &itrk, vector<TrackPM> &ptrk, bool vertical){
  float diff_ang, diff_pos,joilik=-1e-5;
  int joitra=-1;
  bool jointed;
  bool hasingtrk=false;
  for(int j=0;j<itrk.size();j++){

    jointed=false;

    for(int i=0;i<ptrk.size();i++){
      
      if(itrk[j].mod==3){
	if(itrk[j].ipln>1)continue;
	if(ptrk[i].fpln<16)continue;
      }
      else if(vertical){
	if(itrk[j].ipln>1&&(!itrk[j].veto)&&(!itrk[j].edge))continue;
	if(ptrk[i].fpln<16&&ptrk[i].stop)continue;
      }

      diff_ang=ptrk[i].ang-itrk[j].ang;
      diff_pos=(ptrk[i].intcpt+ptrk[i].slope*946.75)-(itrk[j].intcpt+itrk[j].slope*946.75);//zpos = (1079.5+814)/2
      

      if(fabs(diff_ang)<ang_th&&fabs(diff_pos)<pos_th){

	if(jointed){
	  if(joilik>sqrt(fabs(diff_ang)*fabs(diff_ang)/35/35+fabs(diff_pos)*fabs(diff_pos)/85/85))	   
	    ptrk[joitra].ing_trk=false;

	  else
	    continue;
	}


	if(ptrk[i].ing_trk == false){
	  ptrk[i].ing_imod   = itrk[j].mod;
	  ptrk[i].ing_fmod   = itrk[j].mod;
	  ptrk[i].ing_ipln   = itrk[j].ipln;
	  ptrk[i].ing_fpln   = itrk[j].fpln;
	  ptrk[i].ing_trk    = true;
	  ptrk[i].ing_stop   = itrk[j].stop;
	  
	  if(hingtrack[j].fpln==10)
	    ptrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln-1;
	  else
	    ptrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln;
	}
	
	else{
	  if(abs(itrk[j].mod-3) <= abs(ptrk[i].ing_fmod-3))continue;
	  if(itrk[j].ipln < ptrk[i].ing_fpln)continue;
	  
	  ptrk[i].ing_fmod   = itrk[j].mod;
	  ptrk[i].ing_fpln   = itrk[j].fpln;
	  ptrk[i].ing_stop   = itrk[j].stop;
	  
	  if(hingtrack[j].fpln==10)
	    ptrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln-1;
	  else
	    ptrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln;
	}
	for(int ihit=0;ihit<itrk[j].hit.size();ihit++){
	  ptrk[i].hit.push_back(itrk[j].hit[ihit]);
	  //cout<<"all right added="<<itrk[j].hit[ihit].mod<<endl;
	}

	jointed=true;
	joitra=i;
	joilik=sqrt(fabs(diff_ang)*fabs(diff_ang)/35/35+fabs(diff_pos)*fabs(diff_pos)/85/85);
	hasingtrk=true;

      }//if
    }//ptrk
  }//itrk
  return hasingtrk;
};


float PE(float pe,float lope, int mod, int pln, int ch){
  float Pe;
  if(1.422*lope+pe<100)
    Pe=pe;
  else
    Pe=1.422*lope;
  if(mod==16&&pln>0&&ch>=8&&ch<24)
    Pe=Pe/2;
  return Pe;
  /*
  if(mod==16&&pln>0&&ch>=8&&ch<24)
    return pe/1.31;
  else
    return pe;
  */
};

int pdg2num(int pdg){
  int num;
  if(abs(pdg)==13)        num=0;
  else if(abs(pdg)==211)  num=1;
  else if(abs(pdg)==321)  num=2;
  else if(abs(pdg)==2212) num=3;
  else if(abs(pdg)==2112) num=4;
  else if(abs(pdg)==11)   num=5;
  else                    num=6;
  return num;
};

int num2pdg(int *trkpdg){
  int num=0;
  int pdg;
  //for(int i=1;i<5;i++){
  for(int i=1;i<6;i++){
    if(trkpdg[num]<=trkpdg[i])num=i;
  }
  //if(trkpdg[num]==0&&trkpdg[5]>0)num=5;
  //else if(trkpdg[num]==0)num=6;
  if(trkpdg[num]==0)num=6;
  if(num==0)      pdg = 13;
  else if(num==1) pdg = 211;
  else if(num==2) pdg = 321;
  else if(num==3) pdg = 2212;
  else if(num==4) pdg = 2112;
  else if(num==5) pdg = 11;
  else if(num==6) pdg = 0;
  return pdg;
};

TFile *fmucl = new TFile("/home/bquilain/CC0pi_XS/Reconstruction/app/mucl.root");
TH1F* ling = (TH1F*)fmucl -> Get("ling");
TH1F* lsci = (TH1F*)fmucl -> Get("lsci");

/*Benjamin*/
double CLi(float pe, bool scibar){
  double cl;
  int nbin=500,ipe;
  float nmin=0,nmax=500;
  ipe=nbin*(pe-nmin)/(nmax-nmin);
  if(ipe<0)ipe=0;
  if(ipe>nbin-1)ipe=nbin-1;
  if(scibar){
    cl=lsci->GetBinContent(ipe+1);
  }
  else{
    cl=ling->GetBinContent(ipe+1);
  }
  return cl;
};
/*
*/
//void initMuCL(float &TTCL, int &nCL]){
void initMuCL(){
  TTCL=1;
  nCL=0;
};

//void addMuCL(vector<Hits> &allhit, float ang,float &TTCL, int &nCL,int ipln, float intcpt, float slope){
/*Benjamin*/
void addMuCL(vector<Hits> &allhit, float ang,int ipln, float intcpt, float slope){
  bool hitpln[18];
  int stype[2][18];
  float plnpe[18];
  float corr=1;
  memset(hitpln,false,sizeof(hitpln));
  memset(plnpe,0,sizeof(plnpe));
  memset(stype,0,sizeof(stype));
  cout<<"new CL, number of hits="<<allhit.size()<<endl;
  cout<<setprecision(4)<<"(intcpt,grad)="<<intcpt<<", "<<slope<<endl;
  for(int k=0;k<allhit.size();k++){
    if(!allhit[k].isohit)continue;//test
    if(allhit[k].pln==ipln+1)continue;//test
    if(allhit[k].pln==ipln)continue;    
    corr=exp(-(1200-(intcpt+slope*zposi(allhit[k].mod,allhit[k].view,allhit[k].pln)))/2417);
    hitpln[allhit[k].pln]=true;
    plnpe[allhit[k].pln]+=PE(allhit[k].pe,allhit[k].lope,allhit[k].mod,allhit[k].pln,allhit[k].ch)*cos(ang*3.14159265/180)/corr;
    if(allhit[k].pln>0&&allhit[k].ch>=8&&allhit[k].ch<24) stype[0][allhit[k].pln]++;
    else stype[1][allhit[k].pln]++;
    cout<<allhit[k].mod<<", pln="<<allhit[k].pln<<", pe="<<allhit[k].pe<<", pecorr="<<allhit[k].pe*cos(ang*3.14159265/180)/corr<<", PETot="<<PE(allhit[k].pe,allhit[k].lope,allhit[k].mod,allhit[k].pln,allhit[k].ch)*cos(ang*3.14159265/180)/corr<<", dx="<<1/cos(ang*3.14159265/180)<<", att="<<corr<<", L="<<(1200-(intcpt+slope*zposi(allhit[k].mod,allhit[k].view,allhit[k].pln)))/10.<<", z="<<zposi(allhit[k].mod,allhit[k].view,allhit[k].pln)<<endl;
  }
  for(int i=0;i<18;i++){
    if(hitpln[i]){
      nCL++;
      cout<<"pln="<<i<<", pln tot pe="<<plnpe[i]<<endl;
      if(stype[1][i]<=stype[0][i]){
	TTCL*=CLi(plnpe[i],true);
      }
      else{
	TTCL*=CLi(plnpe[i],false);
      }
    }
  }
};
/**/

//float calcMuCL(float &TTCL, int &nCL){
double calcMuCL(){
  double lncli=-log(TTCL);
  double mucl=0;
  double kaijo=1;
  for(int m=0;m<nCL;m++){
    if(m!=0)kaijo*=m;
    mucl+=pow(lncli,m)/kaijo;
  };
  cout<<"TTCL="<<TTCL<<", mucl correction="<<mucl<<endl;
  mucl=TTCL*mucl;
  return mucl;
};

float veract(IngridHitSummary* inghitsum, int xpln, int ypln, float xch, float ych){
  int view = inghitsum->view;
  int ch   = inghitsum->ch;
  int pln  = inghitsum->pln;
  float pe = PE(inghitsum->pe,inghitsum->lope,16,pln,ch);
  if(view==0){
    if(abs(pln-xpln)>2||fabs(xyposi(16,pln,ch)-xch)>155)pe=0;
  }
  else{
    if(abs(pln-ypln)>2||fabs(xyposi(16,pln,ch)-ych)>155)pe=0;
  }
  return pe;
};


void fAddPE(vector<Hits> &allhit, float ang, float &totalpe, int &totalhit, int *trkpdg,int ipln, float intcpt, float slope){
  bool hitpln[18];
  float corr=1;
  memset(hitpln,false,sizeof(hitpln));
  for(int k=0;k<allhit.size();k++){
    if(!allhit[k].isohit)continue;//test
    if(allhit[k].pln==ipln+1)continue;//test
    if(allhit[k].pln==ipln)continue;    
    corr=exp(-(1200-(intcpt+slope*zposi(allhit[k].mod,allhit[k].view,allhit[k].pln)))/2417);
    totalpe+=PE(allhit[k].pe,allhit[k].lope,allhit[k].mod,allhit[k].pln,allhit[k].ch)*cos(ang*3.14159265/180)/corr;
    trkpdg[pdg2num(allhit[k].pdg)]++;
    hitpln[allhit[k].pln]=true;
  }
  for(int i=0;i<18;i++){
    if(hitpln[i])totalhit++;
  }
};

void fTrackMatch(Trk &trk, TrackPM &htrk, TrackPM &vtrk){

  trk.startxpln=htrk.ipln;
  trk.startypln=vtrk.ipln;
  trk.startxch=htrk.ixy;
  trk.startych=vtrk.ixy;
  trk.x=htrk.ixy;
  trk.y=vtrk.ixy;
  trk.endxpln=htrk.fpln;
  trk.endypln=vtrk.fpln;
  trk.endxch=htrk.fxy;
  trk.endych=vtrk.fxy;	
  
  trk.thetax=htrk.ang;
  trk.thetay=vtrk.ang;

  for(int ihit=0;ihit<htrk.hit.size();ihit++){
    trk.hit.push_back(htrk.hit[ihit]);
  }
  for(int ihit2=0;ihit2<vtrk.hit.size();ihit2++){
    trk.hit.push_back(vtrk.hit[ihit2]);
  }  

  trk.intcptx=htrk.intcpt;
  trk.intcpty=vtrk.intcpt;
  trk.slopex=htrk.slope;
  trk.slopey=vtrk.slope;
  
  float trkang=180/3.14159265*atan(sqrt(pow(tan(htrk.ang*3.14159265/180),2)+pow(tan(vtrk.ang*3.14159265/180),2)));
  
  trk.angle=trkang;
  trk.vetowtracking=htrk.veto||vtrk.veto;
  trk.edgewtracking=htrk.edge||vtrk.edge;
  
  if(abs(vtrk.ing_imod-3)<abs(htrk.ing_imod-3))
    trk.ing_startmod=vtrk.ing_imod;
  else
    trk.ing_startmod=htrk.ing_imod;
  
  if(abs(vtrk.ing_fmod-3)>abs(htrk.ing_fmod-3))
    trk.ing_endmod=vtrk.ing_fmod;
  else
    trk.ing_endmod=htrk.ing_fmod;
  
  if(vtrk.ing_ipln<htrk.ing_ipln)
    trk.ing_startpln=vtrk.ing_ipln;
  else
    trk.ing_startpln=htrk.ing_ipln;
  
  if(vtrk.ing_fpln>htrk.ing_fpln)
    trk.ing_endpln=vtrk.ing_fpln;
  else
    trk.ing_endpln=htrk.ing_fpln;
  
  trk.ing_trk=(htrk.ing_trk||vtrk.ing_trk);
  trk.pm_stop=(htrk.stop&&vtrk.stop);
  trk.ing_stop=(htrk.ing_stop&&vtrk.ing_stop);
  
  if(!trk.ing_trk){
    trk.iron_pene=0;
    trk.iron_range=0;
  }
  else if(vtrk.iron_pene>htrk.iron_pene){
    trk.iron_pene=vtrk.iron_pene;
    trk.iron_range=vtrk.iron_pene/cos(trkang*3.14159265/180);
  }
  else{
    trk.iron_pene=htrk.iron_pene;
    trk.iron_range=htrk.iron_pene/cos(trkang*3.14159265/180);
  }
  
  if((vtrk.fpln-vtrk.ipln)>(htrk.fpln-htrk.ipln))
    trk.sci_range=(vtrk.fpln-vtrk.ipln)/cos(trkang*3.14159265/180);
  else
    trk.sci_range=(htrk.fpln-htrk.ipln)/cos(trkang*3.14159265/180);
  
  float totalpe=0;
  int totalhit=0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,vtrk.intcpt,vtrk.slope);
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,htrk.intcpt,htrk.slope);
  if(totalhit>0){
    trk.trkpe=totalpe/totalhit;
  }
  else{
    trk.trkpe=0;
  }
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  //Benjamin
  addMuCL(htrk.hit,trkang,htrk.ipln,vtrk.intcpt,vtrk.slope);
  //Benjamin
  addMuCL(vtrk.hit,trkang,vtrk.ipln,htrk.intcpt,htrk.slope);
  //Benjamin
  trk.mucl=calcMuCL();
  cout<<"Total MuCL="<<trk.mucl<<endl;
};



void fTrackMatchX(Trk &trk,PMTrack &pmtrk, TrackPM &htrk){

  trk.startxpln=htrk.ipln;
  trk.startypln=pmtrk.trk[0].startypln;
  trk.startxch=htrk.ixy;
  trk.startych=pmtrk.trk[0].startych;
  trk.x=htrk.ixy;
  trk.y=pmtrk.trk[0].y;
  trk.endxpln=htrk.fpln;
  trk.endypln=htrk.fpln;
  trk.endxch=htrk.fxy;
  trk.endych=pmtrk.trk[0].endych;
  trk.thetax=htrk.ang;
  trk.thetay=pmtrk.trk[0].thetay;

  for(int ihit=0;ihit<htrk.hit.size();ihit++){
    trk.hit.push_back(htrk.hit[ihit]);
  }
  for(int ihit2=0;ihit2<pmtrk.trk[0].hit.size();ihit2++){
    if(pmtrk.trk[0].hit[ihit2].mod==16 && pmtrk.trk[0].hit[ihit2].pln<=trk.endypln && pmtrk.trk[0].hit[ihit2].view==1) trk.hit.push_back(pmtrk.trk[0].hit[ihit2]);
    else if(htrk.ing_trk && pmtrk.trk[0].hit[ihit2].mod==htrk.ing_imod && pmtrk.trk[0].hit[ihit2].pln<=htrk.ing_fpln && pmtrk.trk[0].hit[ihit2].view==1) trk.hit.push_back(pmtrk.trk[0].hit[ihit2]);
  }

  float trkang=180/3.14159265*atan(sqrt(pow(tan(htrk.ang*3.14159265/180),2)+pow(tan(pmtrk.trk[0].thetay*3.14159265/180),2)));
  
  trk.angle=trkang;
  trk.vetowtracking=htrk.veto||pmtrk.trk[0].vetowtracking;
  trk.edgewtracking=htrk.edge||pmtrk.trk[0].edgewtracking;

  trk.ing_trk=false;
  trk.pm_stop=htrk.stop;

  trk.iron_pene=0;
  trk.iron_range=0;

  trk.sci_range=(htrk.fpln-htrk.ipln)/cos(trkang*3.14159265/180);

  float totalpe=0;
  int totalhit=0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);
  if(totalhit>0){
    trk.trkpe=totalpe/totalhit;
  }
  else{
    trk.trkpe=0;
  }
  trk.pdg=num2pdg(trkpdg);

  //Benjamin
  initMuCL();
  //Benjamin
  addMuCL(htrk.hit,trkang,htrk.ipln,pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);
  //Benjamin
  trk.mucl=calcMuCL();
};



void fTrackMatchY(Trk &trk,PMTrack &pmtrk, TrackPM &vtrk){

  trk.startxpln=pmtrk.trk[0].startxpln;
  trk.startypln=vtrk.ipln;
  trk.startxch=pmtrk.trk[0].startxch;
  trk.startych=vtrk.ixy;
  trk.x=pmtrk.trk[0].x;
  trk.y=vtrk.ixy;
  trk.endxpln=vtrk.fpln;
  trk.endypln=vtrk.fpln;
  trk.endxch=pmtrk.trk[0].endxch;
  trk.endych=vtrk.fxy;
  trk.thetax=pmtrk.trk[0].thetax;
  trk.thetay=vtrk.ang;
  for(int ihit=0;ihit<vtrk.hit.size();ihit++){
    trk.hit.push_back(vtrk.hit[ihit]);
  }
  for(int ihit2=0;ihit2<pmtrk.trk[0].hit.size();ihit2++){
    if(pmtrk.trk[0].hit[ihit2].mod==16 && pmtrk.trk[0].hit[ihit2].pln<=trk.endxpln && pmtrk.trk[0].hit[ihit2].view==0) trk.hit.push_back(pmtrk.trk[0].hit[ihit2]);
    else if(vtrk.ing_trk && pmtrk.trk[0].hit[ihit2].mod==vtrk.ing_imod && pmtrk.trk[0].hit[ihit2].pln<=vtrk.ing_fpln && pmtrk.trk[0].hit[ihit2].view==0) trk.hit.push_back(pmtrk.trk[0].hit[ihit2]);

  }

  float trkang=180/3.14159265*atan(sqrt(pow(tan(vtrk.ang*3.14159265/180),2)+pow(tan(pmtrk.trk[0].thetax*3.14159265/180),2)));
  
  trk.angle=trkang;
  trk.vetowtracking=vtrk.veto||pmtrk.trk[0].vetowtracking;
  trk.edgewtracking=vtrk.edge||pmtrk.trk[0].edgewtracking;

  trk.ing_trk=false;
  trk.pm_stop=vtrk.stop;

  trk.iron_pene=0;
  trk.iron_range=0;

  trk.sci_range=(vtrk.fpln-vtrk.ipln)/cos(trkang*3.14159265/180);

  float totalpe=0;
  int totalhit=0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);
  if(totalhit>0){
    trk.trkpe=totalpe/totalhit;
  }
  else{
    trk.trkpe=0;
  }
  trk.pdg=num2pdg(trkpdg);

  //Benjamin
  initMuCL();
  //Benjamin
  addMuCL(vtrk.hit,trkang,vtrk.ipln,pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);
  //Benjamin
  trk.mucl=calcMuCL();
};



bool fPMAna(){

  PMTrack track;
  Trk     trk;
  pmtrack.clear();

  if(htrack.size()==0||vtrack.size()==0)return false;

  fIngSortTrack(hingtrack);
  fIngSortTrack(vingtrack);
  fSortTrack(htrack);
  fSortTrack(vtrack);

  if(!fIngPMJoint(hingtrack,htrack,false))return false;
  if(!fIngPMJoint(vingtrack,vtrack,true))return false;

  vector<int> id_h,id_v,track_h,track_v;
  vector<bool> tracked_h,tracked_v,used_h,used_v;


  tracked_h.clear();tracked_v.clear();
  for(int i=0;i<htrack.size();i++)tracked_h.push_back(false);
  for(int j=0;j<vtrack.size();j++)tracked_v.push_back(false);
  for(int dif=0;dif<diff_th;dif++){
    for(int pln=0;pln<plnmax(16)-1;pln++){//16 means PM
      id_h.clear();id_v.clear();
      used_h.clear();used_v.clear();
      for(int i=0;i<htrack.size();i++){
	if(!htrack[i].ing_trk)continue;
	if(tracked_h[i])continue;
	if((htrack[i].ipln-pln)>dif||(htrack[i].ipln-pln)<0)continue;
	id_h.push_back(i);
	used_h.push_back(false);
      }
      for(int j=0;j<vtrack.size();j++){
	if(!vtrack[j].ing_trk)continue;
	if(tracked_v[j])continue;
	if((vtrack[j].ipln-pln)>dif||(vtrack[j].ipln-pln)<0)continue;
	id_v.push_back(j);
	used_v.push_back(false);
      }

      track_h.clear();track_v.clear();
      for(int ddif=0;ddif<plnmax(3)-1;ddif++){
	for(int dpln=plnmax(3)-1;dpln>=0;dpln--){
	  for(int i=0;i<id_h.size();i++){
	    if((htrack[id_h[i]].ing_fpln-dpln)>ddif||(htrack[id_h[i]].ing_fpln-dpln)<0)continue;
	    for(int j=0;j<id_v.size();j++){
	      if(htrack[id_h[i]].clstime!=vtrack[id_v[j]].clstime)continue;
	      if(used_h[i])continue;
	      if(used_v[j])continue;
	      if((vtrack[id_v[j]].ing_fpln-dpln)>ddif||(vtrack[id_v[j]].ing_fpln-dpln)<0)continue;
	      track_h.push_back(id_h[i]);
	      track_v.push_back(id_v[j]);
	      used_h[i]=true;
	      used_v[j]=true;
	    }
	  }
	}
      }//ddif

      for(int k=0;k<track_h.size();k++){
	int h=track_h[k];
	int v=track_v[k];
	
	track.clear();
	
	tracked_h[h]=true;
	tracked_v[v]=true;

	fTrackMatch(trk,htrack[h],vtrack[v]);
	
	track.trk.push_back(trk);
	track.clstime=(htrack[h].clstime+htrack[v].clstime)/2;
	track.vetowtracking=htrack[h].veto||vtrack[v].veto;
	track.edgewtracking=htrack[h].edge||vtrack[v].edge;
	track.Ntrack=1;
	track.Ningtrack=1;
	
	pmtrack.push_back(track);
	
      }//for
      
    }//pln
  }//dif
    
    
  //for(int ddif=0;ddif<plnmax(3)-1;ddif++){
  for(int i=0;i<pmtrack.size();i++){
    for(int j=i+1;j<pmtrack.size();j++){

      if(pmtrack[i].Ntrack == 0||pmtrack[j].Ntrack == 0)continue;

      if(abs((pmtrack[i].trk[0].startxpln)-(pmtrack[j].trk[0].startxpln))+abs((pmtrack[i].trk[0].startypln)-(pmtrack[j].trk[0].startypln))>pln_th)continue;
      if(fabs((pmtrack[i].trk[0].y)-(pmtrack[j].trk[0].y))+fabs((pmtrack[i].trk[0].x)-(pmtrack[j].trk[0].x))>ch_th)continue;


      bool former = false;

      if((pmtrack[i].vetowtracking||pmtrack[i].edgewtracking)&&
	 (pmtrack[j].vetowtracking||pmtrack[j].edgewtracking)&&
	 (pmtrack[i].trk[0].endxpln+pmtrack[i].trk[0].endypln)<(pmtrack[j].trk[0].endxpln+pmtrack[j].trk[0].endypln)){
	former=false;
      }
      else if((pmtrack[i].vetowtracking||pmtrack[i].edgewtracking)&&
	      (pmtrack[j].vetowtracking||pmtrack[j].edgewtracking)){
	former=true;
      }     
      else if(pmtrack[j].vetowtracking||pmtrack[j].edgewtracking){
	former=false;
      }
      else if(pmtrack[i].vetowtracking||pmtrack[i].edgewtracking){
	former=true;
      }
      else if((pmtrack[i].trk[0].endxpln+pmtrack[i].trk[0].endypln)<(pmtrack[j].trk[0].endxpln+pmtrack[j].trk[0].endypln)){
	former=false;
      }
      else{
	former=true;
      }

      if(former){
	pmtrack[i].vetowtracking = pmtrack[i].vetowtracking || pmtrack[j].vetowtracking;
	pmtrack[i].edgewtracking = pmtrack[i].edgewtracking || pmtrack[j].edgewtracking;
	pmtrack[i].Ntrack += pmtrack[j].Ntrack;
	pmtrack[j].Ntrack =0;
	pmtrack[i].Ningtrack += pmtrack[j].Ningtrack;
	pmtrack[j].Ningtrack =0;
	for(int t=0;t<pmtrack[j].trk.size();t++)pmtrack[i].trk.push_back(pmtrack[j].trk[t]);
	pmtrack[j].trk.clear();
      }
      else{
	pmtrack[j].vetowtracking = pmtrack[i].vetowtracking || pmtrack[j].vetowtracking;
	pmtrack[j].edgewtracking = pmtrack[i].edgewtracking || pmtrack[j].edgewtracking;
	pmtrack[j].Ntrack += pmtrack[i].Ntrack;
	pmtrack[i].Ntrack =0;
	pmtrack[j].Ningtrack += pmtrack[i].Ningtrack;
	pmtrack[i].Ningtrack =0;
	for(int t=0;t<pmtrack[i].trk.size();t++)pmtrack[j].trk.push_back(pmtrack[i].trk[t]);
	pmtrack[i].trk.clear();
      }
    }
  }
  
    
  for(int pdif=0;pdif<4;pdif++){
    for(int k=0;k<pmtrack.size();k++){
      if(pmtrack[k].Ntrack==0)continue;
      for(int h=0;h<htrack.size();h++){
	if(tracked_h[h])continue;
	for(int v=0;v<vtrack.size();v++){
	  if(tracked_v[v])continue;	  
	  if((htrack[h].fpln-vtrack[v].fpln>pdif+1)||(vtrack[v].fpln-htrack[h].fpln>pdif))continue;
	  if(abs((pmtrack[k].trk[0].startxpln)-(htrack[h].ipln))+abs((pmtrack[k].trk[0].startypln)-(vtrack[v].ipln))>pln_th)continue;
	  if(fabs((pmtrack[k].trk[0].y)-(vtrack[v].ixy))+fabs((pmtrack[k].trk[0].x)-(htrack[h].ixy))>ch_th)continue;
	  
	  tracked_h[h]=true;
	  tracked_v[v]=true;
	  
	  trk.clear();
	  
	  fTrackMatch(trk,htrack[h],vtrack[v]);
	  
	  pmtrack[k].trk.push_back(trk);
	  pmtrack[k].vetowtracking = pmtrack[k].vetowtracking || trk.vetowtracking;
	  pmtrack[k].edgewtracking = pmtrack[k].edgewtracking || trk.edgewtracking;
	  if(htrack[h].ing_trk||vtrack[v].ing_trk)pmtrack[k].Ningtrack++;
	  pmtrack[k].Ntrack++;
	}
      }
    }
  }
  
  //check track match x and y
  for(int k=0;k<pmtrack.size();k++){
    if(pmtrack[k].Ntrack==0)continue;

    for(int h=0;h<htrack.size();h++){
      if(tracked_h[h])continue;
      if(abs((pmtrack[k].trk[0].startxpln)-(htrack[h].ipln))>pln_th)continue;
      if(fabs((pmtrack[k].trk[0].x)-(htrack[h].ixy))>ch_th)continue;
      tracked_h[h]=true;
      trk.clear();
      fTrackMatchX(trk,pmtrack[k],htrack[h]);
      pmtrack[k].trk.push_back(trk);
      pmtrack[k].vetowtracking = pmtrack[k].vetowtracking || trk.vetowtracking;
      pmtrack[k].edgewtracking = pmtrack[k].edgewtracking || trk.edgewtracking;
      if(htrack[h].ing_trk)pmtrack[k].Ningtrack++;
      pmtrack[k].Ntrack++;
    }
    for(int v=0;v<vtrack.size();v++){
      if(tracked_v[v])continue;
      if(abs((pmtrack[k].trk[0].startypln)-(vtrack[v].ipln))>pln_th)continue;
      if(fabs((pmtrack[k].trk[0].y)-(vtrack[v].ixy))>ch_th)continue;
      tracked_v[v]=true;
      trk.clear();
      fTrackMatchY(trk,pmtrack[k],vtrack[v]);
      pmtrack[k].trk.push_back(trk);
      pmtrack[k].vetowtracking = pmtrack[k].vetowtracking || trk.vetowtracking;
      pmtrack[k].edgewtracking = pmtrack[k].edgewtracking || trk.edgewtracking;
      if(vtrack[v].ing_trk)pmtrack[k].Ningtrack++;
      pmtrack[k].Ntrack++;
    }
  }

    
  return true;
};


#endif

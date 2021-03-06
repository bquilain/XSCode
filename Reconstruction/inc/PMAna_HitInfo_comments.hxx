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
#include "PMrecon.hxx"

//__________________________________________________________


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
    time   = -1e-5;
    xy     = -1e-5;
    z      = -1e-5;
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



bool withend(const TrackPM& left, const TrackPM& right){//l'opération pour classer les traces 2D. On les rangent par plans finaux du plus grand au plus petit. SI le plan final est le meme, on les trient du plan initial le plus petit au plus grand (donc de la trace la plus grande a la plus petite). Et si les 2 plans des traces sont les meme, on les trie par angle du plus petit au plus grand.
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

bool fIngPMJoint(vector<TrackIng> &itrk, vector<TrackPM> &ptrk, bool vertical){//continuer la lecture, mais devient vraie si on trouve une trace PM et INGRID qui matchent
  float diff_ang, diff_pos,joilik=-1e-5;//supposons que l'on ait 2 traces horiz dans le PM et une horiz dans Ingrid, et que ce soit 2 évenements indépendants. Déja, ne pourrait t'on matcher leur cluster en temps histoire de voir si rien ne cloche?
  //si les traces ne matchent pas on va avoir 2 continues => on retourne false => et dans PMAna, no reconstruction => On a sélectionné que les évents ou au moins une trace continuait dans INGRID => vire la plupart des NC non?
  //Si les traces muons sortent du PM avant le plan 16 => pk 16 et pas 17???. Surtout si elles sortent avant pour aller dans d'autres modules, on va jeter l'évenement non?
  //NON en fait on ne les a pas viré car le continues sont la si et seulement si on tape le module 3 ou les modules verticaux
  //si aucube trace n'a d'ami dans INGRID, par raison géométrique par exemple, alors on retourne false. ET dans ce cas, dans la recons, on passe alors à une recons suivante => ON ne reconstruit jamais une recons dans laquelle on n'a pas au moins une des traces ayant une suite dans INGRID


  int joitra=-1;
  bool jointed;
  bool hasingtrk=false;//vaut true dès qu'au moins une trace a des amis dans INGRID.
  for(int j=0;j<itrk.size();j++){//INGRID track loop

    jointed=false;

    for(int i=0;i<ptrk.size();i++){//PM Track loop
      
      if(itrk[j].mod==3){//si la trace n'est pas continue du PM vers le Mod3 alors on vire => c'est pas mal mais que faire si on a un plan inactif? (c'est quand meme rare de louper les 2 scintis du plan....)
	if(itrk[j].ipln>1)continue;
	if(ptrk[i].fpln<16)continue;
      }
      else if(vertical){//Inte dans les modules verticaux => avant le PM. SI la trace INGRID est plus loin que le premier plan et qu'elle est véto (donc véto coté) ou edge, on dégage. De meme si le plan final du PM est différent du dernier et que l'on se stoppe dans le PM. CORRIGER CE QUE VEUT DIRE VERTICAL
	if(itrk[j].ipln>1&&(!itrk[j].veto)&&(!itrk[j].edge))continue;
	if(ptrk[i].fpln<16&&ptrk[i].stop)continue;
      }

      diff_ang=ptrk[i].ang-itrk[j].ang;//angle entre traces reco
      diff_pos=(ptrk[i].intcpt+ptrk[i].slope*946.75)-(itrk[j].intcpt+itrk[j].slope*946.75);//zpos = (1079.5+814)/2. Distance entre le début de la trace INGRID et celle PM => a check 
      
      if(fabs(diff_ang)<35&&fabs(diff_pos)<85){//on peut jouer la dessus, ce sont les contraintes sur l'angle de matching PM et INGRID

	if(jointed){//cas ou on a déja la trace INGRID jointe avec une trace PM
	  if(joilik>sqrt(fabs(diff_ang)*fabs(diff_ang)/35/35+fabs(diff_pos)*fabs(diff_pos)/85/85))//on teste ?. Critère pour définir si les traces matchent les mieux	   
	    ptrk[joitra].ing_trk=false;//

	  else
	    continue;
	}


	if(ptrk[i].ing_trk == false){//cas initial: la pmtrack n'a pas encore d'ami ingrid => on lui donne l'ami que l'on vient de sélectionner
//cas ou on a pas réussi à matcher avec la trace INGRID car elle était déja occupé par meilleure qu'elle. On n'aura pas de hastrack
	  ptrk[i].ing_imod   = itrk[j].mod;
	  ptrk[i].ing_fmod   = itrk[j].mod;
	  ptrk[i].ing_ipln   = itrk[j].ipln;
	  ptrk[i].ing_fpln   = itrk[j].fpln;
	  ptrk[i].ing_trk    = true;
	  ptrk[i].ing_stop   = itrk[j].stop;//qui as t'il dans itrk[j].stop et pmtrk[i].stop???? Remplis en true sauf si on touche le dernier scinti du bord d'un module, ou alors, si on touche le dernier plan
	  
	  if(hingtrack[j].fpln==10)
	    ptrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln-1;
	  else
	    ptrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln;
	}
	
	else{//cas ou elle a déja un ami
	  if(abs(itrk[j].mod-3) <= abs(ptrk[i].ing_fmod-3))continue;//teste si notre différence entre protontrack et ingrid track stopped est bien inférieure à ?
	  if(itrk[j].ipln < ptrk[i].ing_fpln)continue;//
	  
	  ptrk[i].ing_fmod   = itrk[j].mod;
	  ptrk[i].ing_fpln   = itrk[j].fpln;
	  ptrk[i].ing_stop   = itrk[j].stop;
	  
	  if(hingtrack[j].fpln==10)
	    ptrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln-1;
	  else
	    ptrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln;
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


void fAddPE(vector<Hits> &allhit, float ang, float &totalpe, int &totalhit, int *trkpdg,int ipln){
  bool hitpln[18];
  memset(hitpln,false,sizeof(hitpln));
  for(int k=0;k<allhit.size();k++){
    if(!allhit[k].isohit)continue;
    if(allhit[k].pln==ipln)continue;
    totalpe+=PE(allhit[k].pe,allhit[k].lope,allhit[k].mod,allhit[k].pln,allhit[k].ch)*cos(ang*3.14159265/180);
    trkpdg[pdg2num(allhit[k].pdg)]++;
    hitpln[allhit[k].pln]=true;
  }
  for(int i=0;i<18;i++){
    if(hitpln[i])totalhit++;
  }
};


void fTrackMatch(Trk &trk, TrackPM &htrk, TrackPM &vtrk){//fait quasi rien a part rentrer le plan de départ et de fin des traces. ET construire les variables de Kikawa pour son analyse (ça on s'en fout)
  trk.clear();
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
  //cout<<"In Track, 2D 1st Track size="<<htrk.hit.size()<<" , second one="<<vtrk.hit.size()<<endl;
  for(int ihit=0;ihit<htrk.hit.size();ihit++){
    trk.hit.push_back(htrk.hit[ihit]);
  }
  for(int ihit2=0;ihit2<vtrk.hit.size();ihit2++){
    trk.hit.push_back(vtrk.hit[ihit2]);
  }
  //cout<<"In Track, #Hits="<<trk.hit.size()<<endl;

  float trkang=180/3.14159265*atan(sqrt(pow(tan(htrk.ang*3.14159265/180),2)+pow(tan(vtrk.ang*3.14159265/180),2)));
  
  trk.angle=trkang;
  trk.vetowtracking=htrk.veto||vtrk.veto;
  trk.edgewtracking=htrk.edge||vtrk.edge;
  
  if(abs(vtrk.ing_imod-3)<abs(htrk.ing_imod-3))//le module de début d'INGRID est storé (on prend celui le moins éloigné du module central, le 3)
    trk.ing_startmod=vtrk.ing_imod;
  else
    trk.ing_startmod=htrk.ing_imod;
  
  if(abs(vtrk.ing_fmod-3)>abs(htrk.ing_fmod-3))//le module final d'INGRID est stocké (cette fois, on prends celui le plus loin du 3. Pourquoi?)
    trk.ing_endmod=vtrk.ing_fmod;
  else
    trk.ing_endmod=htrk.ing_fmod;
  
  if(vtrk.ing_ipln<htrk.ing_ipln)//le plan INGRID initial est le plus petit des 2 => normal
    trk.ing_startpln=vtrk.ing_ipln;
  else
    trk.ing_startpln=htrk.ing_ipln;
  
  if(vtrk.ing_fpln>htrk.ing_fpln)//et le final, le plus grand évidemment
    trk.ing_endpln=vtrk.ing_fpln;
  else
    trk.ing_endpln=htrk.ing_fpln;
  
  trk.ing_trk=(htrk.ing_trk||vtrk.ing_trk);//on est ingrid track si l'une des 2 l'est (l'horiz ou la vert a une trace INGRID)
  trk.pm_stop=(htrk.stop&&vtrk.stop);//si les 2 sont stoppées dans le PM on est PM stop
  trk.ing_stop=(htrk.ing_stop&&vtrk.ing_stop);//si aucune dans le PM, on est INGRID stop
  
  if(!trk.ing_trk){//cas PM uniquement
    trk.iron_pene=0;
    trk.iron_range=0;
  }
  else if(vtrk.iron_pene>htrk.iron_pene){//ensuite c'est les variables de Kikawa
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
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln);
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln);
  if(totalhit>0){
    trk.trkpe=totalpe/totalhit;
  }
  else{
    trk.trkpe=0;
  }
  trk.pdg=num2pdg(trkpdg);
 
};


bool fPMAna(){//fonction principale appelée par PMAna.cc. Si le résultat renvoyé est faux, on n'aura pas de reconstruction 3D

  int pln_th  =2;
  float ch_th =150;

  PMTrack track;
  Trk     trk;
  pmtrack.clear();

  if(htrack.size()==0||vtrack.size()==0)return false;//vérifie qu'on ait forcément au moins une trace horiz ou vert

  fIngSortTrack(hingtrack);//trient les traces de la plus backward à la plus upward
  fIngSortTrack(vingtrack);
  fSortTrack(htrack);
  fSortTrack(vtrack);

  if(!fIngPMJoint(hingtrack,htrack,false))return false;//on vire toutes les traces n'ayant pas de Ingrid.... ça c'est un problème pour notre recons de plein de traces => TO CHANGE
  if(!fIngPMJoint(vingtrack,vtrack,true))return false;

  vector<int> id_h,id_v,track_h,track_v;
  vector<bool> tracked_h,tracked_v,used_h,used_v;
  tracked_h.clear();tracked_v.clear();
  for(int i=0;i<htrack.size();i++)tracked_h.push_back(false);
  for(int j=0;j<vtrack.size();j++)tracked_v.push_back(false);
  //Partie PM

  //toute cette mascarade dit juste qu'on va regarder toute les traces 2D vert et horiz, mais en partant de celles qui matchent dans le meme plan de départ (vertex) et qui ont un plan de départ le plus en amont

  for(int dif=0;dif<4;dif++){//qu'est ce que diff? différence entre le plan étudié (par ordre croissant) et le plan de début de la trace. Au premier coup, on cherche donc toutes les traces dont l'origine est située sur le meme plan qu celui étudié. que signifie le 4? C'est que l'on va tolérer jusqu'à 3 plans d'écart entre le plan étudié et le plan de début de la trace => dans quel but? voir la suite
    //on autorise jusqu'a 3 plans d'écart entre les débuts de 2 traces, une en XZ, l'autre en YZ
    for(int pln=0;pln<plnmax(16)-1;pln++){//16 means PM. C'est le plan de début de la trace!!!!
      id_h.clear();id_v.clear();
      used_h.clear();used_v.clear();
      for(int i=0;i<htrack.size();i++){
	if(!htrack[i].ing_trk)continue;
	if(tracked_h[i])continue;
	if((htrack[i].ipln-pln)>dif||(htrack[i].ipln-pln)<0)continue;//si la trace démarre bien sur le plan étudié (ou plus loin quand diff augmente) on ne s'arrete pas
	id_h.push_back(i);//alors on ajoute l'identifiant de la trace, i
	used_h.push_back(false);//on ajoute une notification false pour chaque trace ajoutée
      }
      for(int j=0;j<vtrack.size();j++){//de meme pour les verticales
	if(!vtrack[j].ing_trk)continue;
	if(tracked_v[j])continue;
	if((vtrack[j].ipln-pln)>dif||(vtrack[j].ipln-pln)<0)continue;
	id_v.push_back(j);
	used_v.push_back(false);
      }

      track_h.clear();track_v.clear();
      //toute cette mascarade dit juste qu'on va regarder toutes les traces 2D qu'on a listé juste avant. Les boucles sont là juste pour dire qu'on regarde d'abord celles qui s'arretent le plus loin dans INGRID
      for(int ddif=0;ddif<plnmax(3)-1;ddif++){//ddif va de 0 à 10 (les plans du module 3)
	for(int dpln=plnmax(3)-1;dpln>=0;dpln--){//dpln va de 10 à 0
	  for(int i=0;i<id_h.size();i++){//on parcourt les identifiant des traces PM qui sont passés plus haut=> dans le plan du PM étudié (ou à dist près)
	    if((htrack[id_h[i]].ing_fpln-dpln)>ddif||(htrack[id_h[i]].ing_fpln-dpln)<0)continue;//on checke si le plan final de notre trace PMM/INGRID est
	    for(int j=0;j<id_v.size();j++){//parcourt identifiant traces vert du PM
	      if(htrack[id_h[i]].clstime!=vtrack[id_v[j]].clstime)continue;//si tps entre les 2 traces ne matche pas => passer à 2 autres.
	      if(used_h[i])continue;//si déja matchées avec d'autres? ça c'est un PB pour le PM qui peut avoir 2 traces dans un plan et confondues dans l'autre
	      if(used_v[j])continue;
	      if((vtrack[id_v[j]].ing_fpln-dpln)>ddif||(vtrack[id_v[j]].ing_fpln-dpln)<0)continue;//on checke que le plan INGRID final de la trace PM est placé avant le plan étudié(de 10 à 0) => veut juste dire qu'on va commencer par étudier les traces s'arretant le plus loin dans INGRID
	      track_h.push_back(id_h[i]);//on rempli alors les track_h et track_v de ces traces, et on les marquent comme déja utilisées. Ce déja utilisé sera remis à 0 uniquement lorsqye l'on changera de plan de début de la trace. Autrement dit, si notre proton commence au meme plan que notre muon, on ne le reconstruit pas
	      track_v.push_back(id_v[j]);
	      used_h[i]=true;
	      used_v[j]=true;
	    }
	  }
	}
      }//ddif

      for(int k=0;k<track_h.size();k++){// on parcourt ces traces horiz. Au premier coup (pln et dif) on regardera par exemple les traces qui commencent dans le PM au premier plan et dont le vertex XZ et YZ est exactement dans le meme plan. 
	int h=track_h[k];
	int v=track_v[k];
	
	track.clear();
	
	tracked_h[h]=true;
	tracked_v[v]=true;
	trk.clear();
	//cout<<"Track Number=("<<h<<","<<v<<")"<<endl;
	fTrackMatch(trk,htrack[h],vtrack[v]);//c'est ici que ça se joue, avant on récolte juste les numéros des traces. Par contre sinon ça fait quasiment rien à part rentrer des infos sur le plan de début et fin de la trace 3D et merger les 2 traces qu'on a mise en une trace 3D (pas de cuts appliqués)
	
	track.trk.push_back(trk);
	track.clstime=(htrack[h].clstime+htrack[v].clstime)/2;//set du timing du cluster: moyenne des 2 clusters vert et horiz
	track.vetowtracking=htrack[h].veto||vtrack[v].veto;
	track.edgewtracking=htrack[h].edge||vtrack[v].edge;
	track.Ntrack=1;
	track.Ningtrack=1;
	
	pmtrack.push_back(track);//on rentre alors ces traces 3d dans pmtrack
	
      }//for
      
    }//pln
  }//dif
    
  //on se retrouve ici avec une liste de traces 3D qui, dans l'ordre croissant: les premières commencent le plus tot dans le PM, ont leur vertex XZ et YZ qui sont parfaitement dans le meme plan, et finissent le plus tard dans INGRID. Elles ont leur cluster en temps qui est le meme. Leur plan final est dans INGRID. Pour améliorer tout ça, il serait bon de virer le fait qu'il faille etre dans INGRID! virer aussi le fait que quand une trace 2d est utilisée avec une autre, elle ne peut plus l'etre avec une autre encore!

  //for(int ddif=0;ddif<plnmax(3)-1;ddif++){
  for(int i=0;i<pmtrack.size();i++){//on parcourt ces nouvelles traces 3d
    for(int j=i+1;j<pmtrack.size();j++){//et les uivantes

      if(pmtrack[i].Ntrack == 0||pmtrack[j].Ntrack == 0)continue;

      if(abs((pmtrack[i].trk[0].startxpln)-(pmtrack[j].trk[0].startxpln))+abs((pmtrack[i].trk[0].startypln)-(pmtrack[j].trk[0].startypln))>pln_th)continue;//si les 2 plans starts des traces reco sont trop proches => on dégage la trace (on ne la stocke pas) => pb pour le proton ça non?
      if(fabs((pmtrack[i].trk[0].y)-(pmtrack[j].trk[0].y)+fabs((pmtrack[i].trk[0].x)-(pmtrack[j].trk[0].x)))>ch_th)continue;//Si les (x,y) de début de la trace sont trop proches => on dégage la trace. de meme, pb pour le proton non?


      bool former = false;

      if((pmtrack[i].vetowtracking||pmtrack[i].edgewtracking)&&
	 (pmtrack[j].vetowtracking||pmtrack[j].edgewtracking)&&
	 (pmtrack[i].trk[0].endxpln+pmtrack[i].trk[0].endypln)<(pmtrack[j].trk[0].endxpln+pmtrack[j].trk[0].endypln)){//si les 2 traces sont soit véto, soit edge et que le plan FINAL de la ieme est en amont de celui de la jeme => former = false
	former=false;
      }
      else if((pmtrack[i].vetowtracking||pmtrack[i].edgewtracking)&&//si les 2 sont deges ou véto mais que c'est la iè qui descend la PLUS EN AVAL dans le PM, former=true
	      (pmtrack[j].vetowtracking||pmtrack[j].edgewtracking)){
	former=true;
      }     
      else if(pmtrack[j].vetowtracking||pmtrack[j].edgewtracking){//si juste la je est véto ou ege: former = false
	former=false;
      }
      else if(pmtrack[i].vetowtracking||pmtrack[i].edgewtracking){// si juste la ieme est véto ou edge: former = true
	former=true;
      }
      else if((pmtrack[i].trk[0].endxpln+pmtrack[i].trk[0].endypln)<(pmtrack[j].trk[0].endxpln+pmtrack[j].trk[0].endypln)){//si aucune n'est véto et que la jeme descend la plus bas: former = false
	former=false;
      }
      else{//si auxune véto et la ieme descend le plus bas => former=true
	former=true;
      }

      if(former){//si former<=>. On remplit alors la trace Ieme des détails additionnels de la je (véto ou edge) + on augmente son nombre de trace de +1! Les pmtracks sont des reco, pas des traces 3D attention donc!
	pmtrack[i].vetowtracking = pmtrack[i].vetowtracking || pmtrack[j].vetowtracking;
	pmtrack[i].edgewtracking = pmtrack[i].edgewtracking || pmtrack[j].edgewtracking;
	pmtrack[i].Ntrack += pmtrack[j].Ntrack;
	pmtrack[j].Ntrack =0;
	pmtrack[i].Ningtrack += pmtrack[j].Ningtrack;
	pmtrack[j].Ningtrack =0;
	for(int t=0;t<pmtrack[j].trk.size();t++)pmtrack[i].trk.push_back(pmtrack[j].trk[t]);
	pmtrack[j].trk.clear();
      }
      else{//si on est former false, c'est les détails de la ième qu'on va ajouter à ceux de la jème
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


  for(int k=0;k<pmtrack.size();k++){//parcourt les recos? => oui
    if(pmtrack[k].Ntrack==0)continue;//si aucune trace, passer à la reco suivante
    for(int h=0;h<htrack.size();h++){//parcourir les traces horiz. Lesquelles. C'est celles qui doivent avoir du INGRID non?
      if(tracked_h[h])continue;//?
      for(int v=0;v<vtrack.size();v++){
	if(tracked_v[v])continue;

	if(abs((pmtrack[k].trk[0].startxpln)-(htrack[h].ipln))+abs((pmtrack[k].trk[0].startypln)-(vtrack[v].ipln))>pln_th)continue;
	if(fabs((pmtrack[k].trk[0].y)-(vtrack[v].ixy)+fabs((pmtrack[k].trk[0].x)-(htrack[h].ixy)))>ch_th)continue;

	tracked_h[h]=true;
	tracked_v[v]=true;

	trk.clear();

	fTrackMatch(trk,htrack[h],vtrack[v]);

	pmtrack[k].trk.push_back(trk);
	pmtrack[k].vetowtracking = pmtrack[k].vetowtracking || trk.vetowtracking;
	pmtrack[k].edgewtracking = pmtrack[k].vetowtracking || trk.edgewtracking;
	pmtrack[k].Ntrack ++;
	
      }
    }
  }




  return true;
};


#endif

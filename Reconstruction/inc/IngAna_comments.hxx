#ifndef __INGANA_HXX__
#define __INGANA_HXX__

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
int diff_th=4;

class IngTrack{
public:
  Int_t    startxpln;
  Int_t    startypln;
  Int_t    startxch;
  Int_t    startych;
  Int_t    endxpln;
  Int_t    endypln;
  Int_t    endxch;
  Int_t    endych;
  Int_t    vertexz;
  Int_t    vertexxz;
  Int_t    vertexyz;
  Float_t  thetax;
  Float_t  thetay;
  Float_t  angle;
  Float_t  clstime;
  Float_t  x;
  Float_t  y;
  Float_t  z;
  Float_t  zx;
  Float_t  zy;
  Bool_t   vetowtracking; // Upstream VETO
  Bool_t   edgewtracking; // Fiducial CUT
  Bool_t   vtxtrk;
  Float_t    vetodist;
  Int_t Ntrack;

  void clear(){
    startxpln=  -1;
    startypln=  -1;
    startxch=  -1;
    startych=  -1;
    endxpln=  -1;
    endypln=  -1;
    endxch=  -1;
    endych=  -1;
    vertexz=  -1;
    vertexxz=  -1;
    vertexyz=  -1;
    thetax=  -1.e-5;
    thetay=  -1.e-5;
    angle=  -1.e-5;
    clstime=  -1.e-5;
    x=  -1.e-5;
    y=  -1.e-5;
    z=  -1.e-5;
    zx=  -1.e-5;
    zy=  -1.e-5;
    vetowtracking = false; // Upstream VETO
    edgewtracking = false; // Fiducial CUT
    vtxtrk = true;
    vetodist=  -1.e-5;
    Ntrack=-1;
  }
};

vector<IngTrack> ingtrack;
vector<Track> htrack;
vector<Track> vtrack;


void fIngAna(int mod){
  IngTrack track;
  ingtrack.clear();
  if(htrack.size()==1&&vtrack.size()==1){//cas ou, dans l'évenement, module, cycle, on n'a qu'une trace 2D vert et une 2D horiz=> On met une seule trace 3D dans track. Les conditions vérifiées: que le plan de départ ne soit pas éloigné de plus de 4plans entre horiz et vertcial+que le plan final de chacune soit plus loin que le début de trace de l'autre
    if(abs(htrack[0].ipln-vtrack[0].ipln)<4&&//dans track, on va mettre les infos comme si c'était une trace 3D avec pour check seulement que on ait 1,2, ou 3 plans de différences entre les 2 début des 2 traces
       htrack[0].fpln>vtrack[0].ipln&&
       vtrack[0].fpln>htrack[0].ipln){
      
      track.startxpln     = htrack[0].ipln;
      track.startypln     = vtrack[0].ipln;
      track.startxch      = (int)htrack[0].ixy/50;
      track.startych      = (int)vtrack[0].ixy/50;
      track.endxpln       = htrack[0].fpln;
      track.endypln       = vtrack[0].fpln;
      track.endxch        = (int)htrack[0].fxy/50;
      track.endych        = (int)vtrack[0].fxy/50;
      
      if(htrack[0].ipln<vtrack[0].ipln){
	track.vertexz       = htrack[0].ipln;
	track.z             = htrack[0].iz;
	}
      else{
	track.vertexz       = vtrack[0].ipln;
	track.z             = vtrack[0].iz;
      }
      
      track.vertexxz      = htrack[0].ipln;
      track.vertexyz      = vtrack[0].ipln;
      track.zx            = htrack[0].iz;
      track.zy            = vtrack[0].iz;
      track.y             = htrack[0].ixy;
      track.x             = vtrack[0].ixy;
      track.thetax        = htrack[0].ang;
      track.thetay        = vtrack[0].ang;
      track.angle         = 180/3.14159265*atan(sqrt(pow(tan(htrack[0].ang*3.14159265/180),2)+pow(tan(vtrack[0].ang*3.14159265/180),2)));
      track.clstime       = htrack[0].clstime;
      track.vetowtracking = (htrack[0].veto||vtrack[0].veto);
      track.edgewtracking = (htrack[0].edge||vtrack[0].edge);
      track.vtxtrk = true;
      if(htrack[0].vetodist<vtrack[0].vetodist)
	track.vetodist = htrack[0].vetodist;
      else
	track.vetodist = vtrack[0].vetodist;

      ingtrack.push_back(track);      
    }
  }

  else if(htrack.size()>=1&&vtrack.size()>=1){//cas ou plein de traces 2D
  vector<int> id_h,id_v,track_h,track_v;
  vector<bool> tracked_h,tracked_v,used_h,used_v;
  tracked_h.clear();tracked_v.clear();
  for(int i=0;i<htrack.size();i++)tracked_h.push_back(false);//replit la vzriable tracked de false
  for(int j=0;j<vtrack.size();j++)tracked_v.push_back(false);
  for(int dif=0;dif<diff_th;dif++){//au max, diff_th=4. c'est le nombre de plans limite d'écart avec le début de la trace
    for(int pln=plnmax(mod)-1;pln>=0;pln--){//du plan le plus up au plus down: mais pourquoi? le vertex?
      id_h.clear();id_v.clear();
      used_h.clear();used_v.clear();
      for(int i=0;i<htrack.size();i++){//parcourt les traces horizontales préstockées
	if(tracked_h[i])continue;//ce qui n'est pas le cas au début, mais ensuite, évite qu'une trace 2D soit utilisée deux fois => est ce forcément mal dans des cas très inélastiques????
	//if(htrack[i].ipln!=pln)continue;
	if((htrack[i].ipln-pln)>dif||(htrack[i].ipln-pln)<0)continue;//si le plan testé est trop loin du plan de dpéart d'une trace 2D=> on va a la suite
       	id_h.push_back(i);//numéro de la trace 2D
	used_h.push_back(false);
      }
      for(int j=0;j<vtrack.size();j++){//de meme que pour traces horizontales
	if(tracked_v[j])continue;
	//if(abs(vtrack[j].ipln-pln)>dif)continue;
	if((vtrack[j].ipln-pln)>dif||(vtrack[j].ipln-pln)<0)continue;
	id_v.push_back(j);
	used_v.push_back(false);
      }

      track_h.clear();track_v.clear();
      for(int ddif=0;ddif<plnmax(mod)-1;ddif++){
	for(int dpln=plnmax(mod)-1;dpln>=0;dpln--){//ça va etre seulement de la mise en commun de 2 traces 2D en une 3D. On peut avoir plusieurs traces 3D, mais elles seront mises dans un track différent
	  for(int i=0;i<id_h.size();i++){
	    if((htrack[id_h[i]].fpln-dpln)>ddif||(htrack[id_h[i]].fpln-dpln)<0)continue;//vérifie cette fois le plan de la fin
	    for(int j=0;j<id_v.size();j++){
	      if(htrack[id_h[i]].clstime!=vtrack[id_v[j]].clstime)continue;//vérifie compatibilité en timing
	      if(htrack[id_h[i]].fpln<=vtrack[id_v[j]].ipln||
		 vtrack[id_v[j]].fpln<=htrack[id_h[i]].ipln)continue;//compatibilité que la fin d'une trace doit etre au dessus du début d'une autre
	      if(used_h[i])continue;
	      if(used_v[j])continue;
	      if((vtrack[id_v[j]].fpln-dpln)>ddif||(vtrack[id_v[j]].fpln-dpln)<0)continue;
	      track_h.push_back(id_h[i]);//track_h et track_v contiennent donc les traces 2D appariées par numéro (track_h(5) et track_v(5) vont ensembles)
	      track_v.push_back(id_v[j]);
	      used_h[i]=true;
	      used_v[j]=true;
	    }
	  }
	}
      }


      for(int k=0;k<track_h.size();k++){//donc plusieurs traces 3D=différentes track
	int h=track_h[k];//contient le numéro des traces 2D qui sont appareillées
	int v=track_v[k];
	
	tracked_h[h]=true;
	tracked_v[v]=true;
	track.clear();

	track.startxpln     = htrack[h].ipln;
	track.startypln     = vtrack[v].ipln;
	track.startxch      = (int)htrack[h].ixy/50;
	track.startych      = (int)vtrack[v].ixy/50;
	track.endxpln       = htrack[h].fpln;
	track.endypln       = vtrack[v].fpln;
	track.endxch        = (int)htrack[h].fxy/50;
	track.endych        = (int)vtrack[v].fxy/50;

	if(track.startxch<0) track.startxch=0;
	if(track.startych<0) track.startych=0;
	if(track.endxch<0)   track.endxch=0;
	if(track.endych<0)   track.endych=0;
	if(track.startxch>23) track.startxch=23;
	if(track.startych>23) track.startych=23;
	if(track.endxch>23)   track.endxch=23;
	if(track.endych>23)   track.endych=23;
	
	if(htrack[h].ipln<vtrack[v].ipln){
	  track.vertexz       = htrack[h].ipln;
	  track.z             = htrack[h].iz;
	}
	else{
	  track.vertexz       = vtrack[v].ipln;
	  track.z             = vtrack[v].iz;
	}
	
	track.vertexxz      = htrack[h].ipln;
	track.vertexyz      = vtrack[v].ipln;
	track.zx            = htrack[h].iz;
	track.zy            = vtrack[v].iz;
	track.y             = htrack[h].ixy;
	track.x             = vtrack[v].ixy;
	track.thetax        = htrack[h].ang;
	track.thetay        = vtrack[v].ang;
	track.angle         = 180/3.14159265*atan(sqrt(pow(tan(htrack[h].ang*3.14159265/180),2)+pow(tan(vtrack[v].ang*3.14159265/180),2)));
	track.clstime       = htrack[h].clstime;
	track.vetowtracking = (htrack[h].veto||vtrack[v].veto);
	track.edgewtracking = (htrack[h].edge||vtrack[v].edge);
	if(htrack[h].vetodist<vtrack[v].vetodist)
	  track.vetodist =htrack[h].vetodist;
	else
	  track.vetodist =vtrack[v].vetodist;


	track.vtxtrk = true;//de base vtxtrk est true

	ingtrack.push_back(track);
	
      }



    }
  }
  }
};


#endif

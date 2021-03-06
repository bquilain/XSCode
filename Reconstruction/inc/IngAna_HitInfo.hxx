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
  vector < vector<int> > hitnum;

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
    for(int i=0;i<hitnum.size();i++){
      hitnum[i].clear();
    }
    hitnum.clear();
  }
};

vector<IngTrack> ingtrack;
vector<Track> htrack;
vector<Track> vtrack;


void fIngAna(int mod){
   IngTrack track;
  ingtrack.clear();
  track.clear();

  if(htrack.size()==1&&vtrack.size()==1){
    if(abs(htrack[0].ipln-vtrack[0].ipln)<4&&
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
      track.hitnum.push_back(htrack[0].hitid);
      for(int i=0;i<vtrack[0].hitid.size();i++){
	track.hitnum[0].push_back(vtrack[0].hitid[i]);
      }
      /*cout<<"In Ana : 1 Track 3D"<<endl;
      for(int i=0;i<track.hitnum.size();i++){//parcours des hits/trace 2D
	for(int j=0;j<track.hitnum[i].size();j++){
	  cout<<track.hitnum[i][j]<<"   ";
	}
      }
      cout<<endl;*/

      if(htrack[0].vetodist<vtrack[0].vetodist)
	track.vetodist = htrack[0].vetodist;
      else
	track.vetodist = vtrack[0].vetodist;

      ingtrack.push_back(track);      
    }
  }

  else if(htrack.size()>=1&&vtrack.size()>=1){
  vector<int> id_h,id_v,track_h,track_v;
  vector<bool> tracked_h,tracked_v,used_h,used_v;
  tracked_h.clear();tracked_v.clear();
  for(int i=0;i<htrack.size();i++)tracked_h.push_back(false);
  for(int j=0;j<vtrack.size();j++)tracked_v.push_back(false);
  for(int dif=0;dif<diff_th;dif++){
    for(int pln=plnmax(mod)-1;pln>=0;pln--){
      id_h.clear();id_v.clear();
      used_h.clear();used_v.clear();
      for(int i=0;i<htrack.size();i++){
	if(tracked_h[i])continue;
	//if(htrack[i].ipln!=pln)continue;
	if((htrack[i].ipln-pln)>dif||(htrack[i].ipln-pln)<0)continue;
       	id_h.push_back(i);
	used_h.push_back(false);
      }
      for(int j=0;j<vtrack.size();j++){
	if(tracked_v[j])continue;
	//if(abs(vtrack[j].ipln-pln)>dif)continue;
	if((vtrack[j].ipln-pln)>dif||(vtrack[j].ipln-pln)<0)continue;
	id_v.push_back(j);
	used_v.push_back(false);
      }

      track_h.clear();track_v.clear();
      for(int ddif=0;ddif<plnmax(mod)-1;ddif++){
	for(int dpln=plnmax(mod)-1;dpln>=0;dpln--){
	  for(int i=0;i<id_h.size();i++){
	    if((htrack[id_h[i]].fpln-dpln)>ddif||(htrack[id_h[i]].fpln-dpln)<0)continue;
	    for(int j=0;j<id_v.size();j++){
	      if(htrack[id_h[i]].clstime!=vtrack[id_v[j]].clstime)continue;
	      if(htrack[id_h[i]].fpln<=vtrack[id_v[j]].ipln||
		 vtrack[id_v[j]].fpln<=htrack[id_h[i]].ipln)continue;
	      if(used_h[i])continue;
	      if(used_v[j])continue;
	      if((vtrack[id_v[j]].fpln-dpln)>ddif||(vtrack[id_v[j]].fpln-dpln)<0)continue;
	      track_h.push_back(id_h[i]);
	      track_v.push_back(id_v[j]);
	      used_h[i]=true;
	      used_v[j]=true;
	    }
	  }
	}
      }


      for(int k=0;k<track_h.size();k++){
	//track.clear();
	//cout<<track.hitnum

	int h=track_h[k];
	int v=track_v[k];
	
	tracked_h[h]=true;
	tracked_v[v]=true;
	
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
	track.hitnum.push_back(htrack[h].hitid);
	for(int i=0;i<vtrack[v].hitid.size();i++){
	  track.hitnum.back().push_back(vtrack[v].hitid[i]);
	} 

	if(htrack[h].vetodist<vtrack[v].vetodist)
	  track.vetodist =htrack[h].vetodist;
	else
	  track.vetodist =vtrack[v].vetodist;

	//for(int i=0;i<track.hitnum.
	/*cout<<"In Ana: several tracks"<<endl;
	for(int i=0;i<track.hitnum.size();i++){//parcours des hits/trace 2D                                                                          
    	  for(int j=0;j<track.hitnum[i].size();j++){
	    cout<<track.hitnum[i][j]<<"   ";
	  }
	}
	cout<<endl;
	for(int i=0;i<htrack[h].hitid.size();i++){
          cout<<htrack[h].hitid[i]<<"   ";
        }
        for(int i=0;i<vtrack[v].hitid.size();i++){
          cout<<vtrack[v].hitid[i]<<"   ";
        }

	cout<<endl;*/
	track.vtxtrk = true;

	ingtrack.push_back(track);

	//hitid.clear();
	
      }



    }
  }
  }
};


#endif

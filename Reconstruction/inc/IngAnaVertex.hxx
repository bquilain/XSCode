#ifndef __INGANAVERTEX_HXX__
#define __INGANAVERTEX_HXX__

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

int pln_th  =2;
float ch_th =150;

void fIngVertex(int mod){

  for(int i=0;i<ingtrack.size();i++){
    ingtrack[i].Ntrack=1;
  }

  for(int i=0;i<ingtrack.size();i++){
    for(int j=i+1;j<ingtrack.size();j++){
      /*
      if(abs((ingtrack[i].startxpln)-(ingtrack[j].startxpln))>pln_th)continue;
      if(abs((ingtrack[i].startypln)-(ingtrack[j].startypln))>pln_th)continue;
      if(fabs((ingtrack[i].x)-(ingtrack[j].x))>ch_th)continue;
      if(fabs((ingtrack[i].y)-(ingtrack[j].y))>ch_th)continue;
      */
      
      if(abs((ingtrack[i].startxpln)-(ingtrack[j].startxpln))+abs((ingtrack[i].startypln)-(ingtrack[j].startypln))>pln_th)continue;
      if(fabs((ingtrack[i].y)-(ingtrack[j].y)+fabs((ingtrack[i].x)-(ingtrack[j].x)))>ch_th)continue;

      if(!(ingtrack[i].vtxtrk)||!(ingtrack[j].vtxtrk))continue;

      if((ingtrack[i].vetowtracking||ingtrack[i].edgewtracking)&&
	 (ingtrack[j].vetowtracking||ingtrack[j].edgewtracking)&&
	 (ingtrack[i].endxpln+ingtrack[i].endypln)<(ingtrack[j].endxpln+ingtrack[j].endypln)){
	ingtrack[i].vtxtrk = false;
	ingtrack[j].Ntrack += ingtrack[i].Ntrack;
      }
      else if((ingtrack[i].vetowtracking||ingtrack[i].edgewtracking)&&
	      (ingtrack[j].vetowtracking||ingtrack[j].edgewtracking)){
	ingtrack[j].vtxtrk = false;
	ingtrack[i].Ntrack += ingtrack[j].Ntrack;
      }     
      else if(ingtrack[j].vetowtracking||ingtrack[j].edgewtracking){
	ingtrack[i].vtxtrk = false;
	ingtrack[j].Ntrack += ingtrack[i].Ntrack;
      }
      else if(ingtrack[i].vetowtracking||ingtrack[i].edgewtracking){
	ingtrack[j].vtxtrk = false;
	ingtrack[i].Ntrack += ingtrack[j].Ntrack;
      }
      else if((ingtrack[i].endxpln+ingtrack[i].endypln)<(ingtrack[j].endxpln+ingtrack[j].endypln)){
	ingtrack[i].vtxtrk = false;
	ingtrack[j].Ntrack += ingtrack[i].Ntrack;
      }
      else{
	ingtrack[j].vtxtrk = false;
	ingtrack[i].Ntrack += ingtrack[j].Ntrack;
      }

    }
  }


};



/*


void fIngVertex(int mod){
  vector<IngTrack>::iterator it1;
  vector<IngTrack>::iterator it2;

  for(it1=ingtrack.begin();it1<ingtrack.end();it1++){
    for(it2=it1+1;it2<ingtrack.end();it2++){
      if(abs((it1->startxpln)-(it2->startxpln))>pln_th)continue;
      if(abs((it1->startypln)-(it2->startypln))>pln_th)continue;
      if(fabs((it1->x)-(it2->x))>ch_th)continue;
      if(fabs((it1->y)-(it2->y))>ch_th)continue;

      if((it1->vetowtracking||it1->edgewtracking)&&
	 (it2->vetowtracking||it2->edgewtracking)&&
	 (it1->endxpln+it1->endypln)<(it2->endxpln+it2->endypln)){
	it1 = ingtrack.erase(it1);
	it1--;
      }
      else if((it1->vetowtracking||it1->edgewtracking)&&
	      (it2->vetowtracking||it2->edgewtracking)){
	it2 = ingtrack.erase(it2);
	it2--;
      }     
      else if(it2->vetowtracking||it2->edgewtracking){
	it1 = ingtrack.erase(it1);
	it1--;
      }
      else if(it1->vetowtracking||it1->edgewtracking){
	it2 = ingtrack.erase(it2);
	it2--;
      }
      else if((it1->endxpln+it1->endypln)<(it2->endxpln+it2->endypln)){
	it1 = ingtrack.erase(it1);
	it1--;
      }
      else{
	it2 = ingtrack.erase(it2);
	it2--;
      }
    }
  }
};
*/

#endif

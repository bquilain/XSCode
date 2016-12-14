#ifndef __INGBASICRECON_HXX__
#define __INGBASICRECON_HXX__

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

//### INGRID-j data structure ####
//################################
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "IngridConstants.h"

const static int BREC_MAXLINE = 36;
const static float UVETOthreshold = 7.5;
class IngBasicRecon{
private:

  IngridHitSummary*   hit[INGRIDHIT_MAXHITS];
  //IngridHitSummary   hit;
  int    nhits;
  int    trgbit;

  float  layerpe;           //## sum of p.e. in active plane(tpl 1~10)/(active plane)
  float  layerpe_wo_tpl10;
  int    nactpln;
  int    nactpln_wo_tpl10;  //##
  bool   upstreamveto;
  bool        edgehit;
  int         vertexz;
  vector<int> vertexx;
  vector<int> vertexy;
  bool   isActive[nTPL+nVETO];

public:
  IngBasicRecon();


  ~IngBasicRecon(){};
  bool   SetHit(TRef inghit[]);
  int    fTrgbit();
  float  fLayerpe();
  int    fNactpln();
  bool   fUVETOHit();
  bool   fGetVertex();
  bool   fEdgeHit(){
    return edgehit;
  };
  int    fVertexz(){
    return vertexz;
  }
  vector<int> fVertexx(){
    return vertexx;
  }
  vector<int> fVertexy(){
    return vertexy;
  }

  //___________
  int fActInARow();
  //___________
  bool   Clear();
  bool   Print();
  bool   PrintAll();
  int    Nactpln(){return nactpln;}

  //### for rock muon study - 2
  float veto_xz; //### hit position when upstream VETO selection
  float veto_x;  //### hit position when upstream VETO selection
  float veto_yz; //### hit position when upstream VETO selection
  float veto_y;  //### hit position when upstream VETO selection

  //### for rock muon study
  int           most_upstream_tpl;
  int           most_downstream_tpl;
  vector<float> mupyz;
  vector<float> mupy;
  vector<float> mupxz;
  vector<float> mupx;
  vector<float> mdownyz;
  vector<float> mdowny;
  vector<float> mdownxz;
  vector<float> mdownx;
  TF1*          fLineX[BREC_MAXLINE];
  TF1*          fLineY[BREC_MAXLINE];
  TF1*          LineX(int i){
    if(i<nlx&&-1<i)return fLineX[i];
    else return 0;
  }
  TF1*          LineY(int i){
    if(i<nly&&-1<i)return fLineY[i];
    else return 0;
  }
  int           nlx;
  int           nly;
  bool GetMU_MD();
  bool GetLine();

  //### temporary for test
  int    inarowbit;
  vector<float> vetope;
  float vpe;
  int   vview ;
  int   vpln;
  int   vch;
};

#endif

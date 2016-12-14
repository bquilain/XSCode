#ifndef _INGSYSSEVERALFID_HXX__
#define _INGSYSSEVERALFID_HXX__

// ROOT includes
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TBox.h"
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
#include "INGRID_Dimension.hxx"
#include "IngridConstants.h"
FileStat_t fs;

const int nfid            = 6;
const int fidverxmin[nfid]  = {  2,  3,  4,  2,  2,  2};
const int fidverxmax[nfid]  = { 21, 20, 19, 21, 21, 21};
const int fidverymin[nfid]  = {  2,  3,  4,  2,  2,  2};
const int fidverymax[nfid]  = { 21, 20, 19, 21, 21, 21};
const int fidverzmin[nfid]  = {  1,  1,  1,  1,  4,  6};
const int fidverzmax[nfid]  = {  9,  9,  9,  3,  5,  9};

char EvtSelect[300];
char FidCut   [300];

IngridEventSummary* evt = new IngridEventSummary();
//TFile*              rfile;
//TTree*              tree;
//TBranch*            EvtBr;
const int nfile = 1;//7;
char* FileName[nfile] = 
  { 
    //"/home/ingrid/data/MC/10c/v1/numu/fe/ana/10c_nd3_numu_FE_1_trk.root",
    "/home/ingrid/data/MC/10c/v3/all/merged100.root",
    //"/home/daq/data/dst/merge.root",
    //"/home/daq/data/dst/merged_29_onlybrec.root",
    //"/home/daq/data/dst/merged_30_track.root",
    //"/home/daq/data/dst/merged_31_track.root",
    //"/home/daq/data/dst/merged_32_track.root",
    //"/home/daq/data/dst/merged_33_track.root",
    //"/home/daq/data/dst/merged_33_2_track.root",
    //"/home/daq/data/dst/merged_34_track.root"

  };

//const int nfile = 7;
//char* FileName[nfile] = 
//{ 


//  };

#endif

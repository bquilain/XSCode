#ifndef __INGANACOSMIC_HXX__
#define __INGANACOSMIC_HXX__

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


class IngAnaCosmic{
private:

  IngridHitSummary*   hit[INGRIDHIT_MAXHITS];
  int    nhits;


public:
  IngAnaCosmic();
  ~IngAnaCosmic(){};
  void        Clear();
  bool        SetHit(TRef inghit[]);
  vector<int> fTdcHit(int view, int pln);
  vector<int> fTrigger(int view, int pln, float thr);
  
};

#endif

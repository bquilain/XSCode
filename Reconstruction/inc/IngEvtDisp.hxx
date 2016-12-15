#ifndef _INGEVTDISP_HXX__
#define _INGEVTDISP_HXX__
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


class IngEvtDisp{
private:
  TStyle*         fStyle;
  TPad*           fPMod;
  TPad*           fPSideView;
  TPad*           fPTopView;
  TPad*           fPInfo;
  TH1F*           fHSideView;
  TH1F*           fHTopView;
  TBox*           fBXScinti;
  TLine*          TLscinti[Nscinti];
  TBox*           fBXIron[NIron];
  TBox*           fBYScinti;
  TBox*           fBYIron[NIron];
  TBox*           fBXGap[NIron + 1][2]; //Gap Between Iron and TPL                
  TBox*           fBYGap[NIron + 1][2]; //Gap Between Iron and TPL                
  TLatex*         fLatexModInfo;
  TLatex*         fLatexEvtInfo;
  TLatex*         fLatexBeamInfo;
  TMarker*        fTMdepX;
  TMarker*        fTMdepY;
  vector<TMarker> fTMVecX;
  vector<TMarker> fTMVecY;

  TF1*            fFtrack_x[10];
  TF1*            fFtrack_y[10];
  TF1*            mutrkx  [10];
  TF1*            mutrky  [10];
  int ntrack_x;
  int ntrack_y;


public:
  IngEvtDisp(){
  };
  ~IngEvtDisp(){};
  bool     Draw_Module( TCanvas& canvas, int mod);
  bool     Draw_Hit   ( Ingrid1stReducSummary& reduc,
			double msize);
  bool     Draw_BeamInfo   ( BeamInfoSummary& beam);
  

  //bool     Draw_Line  ( TF1& line, int view );
  bool     Draw_Line  ( TF1& line, int view, float iz, float fz ,int color, int size);
		        
  /*
  template <class T> bool Draw_Hit_A(T& sum, 
				     double msize,
				     int color,
				     int flag);
  */
  template <class T> bool Draw_Hit_A(T& sum, 
				     double msize,
				     int color,
				     int flag);
  template <class T> bool Print_Hit_A(T& sum);
 
  bool Print_Muon (IngridEventSummary& evt); 
  bool Draw_Muon  (IngridEventSummary& evt); 

				
				
				


  //### temprorary for rock muon study //####
  vector<TF1>     fTFlineX;
  vector<TF1>     fTFlineY;
  bool     Draw_Line_all  ();
  bool reset_track();
};



#endif

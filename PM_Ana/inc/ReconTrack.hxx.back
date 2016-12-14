#ifndef _RECONTRACK_HXX
#define _RECONTRACK_HXX

#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "setup.hxx"

//### INGRID-j data structure ####
//################################
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"

const static int  THETA_RESOLUTION =              200; // # of bins of theta
const static int  RHO_RESOLUTION   =             1300; // # of bins of rho
const static int  RHO_RESOLUTION_H =  RHO_RESOLUTION/2;//
const static int  isReconTrackN    =                3; // 
const static int  ReconTolerance   =               10;
const static int  NHitofTrack      =                4;
const static int  NHitLyr          =                4;

const static int  N_BREAKHOUGH     =                3; // because at least 3 hits
const static double            pai = 6*asin(0.5)/THETA_RESOLUTION;
const static int  max_ntry         =                1;
const static int  max_ntrack       =               30;

const static int  reset_region_theta = 20;
const static int  reset_region_rho   = 20;
//_______________________________________________________________
//_______________________________________________________________
class ReconTrack{
private:
  double          sn[THETA_RESOLUTION];   //### sin array for quick calc. 
  double          cs[THETA_RESOLUTION];   //### cos array for quick calc. 
  int             theta;                  //### variables for theta-rho space
  int             rho;                    //### variables for theta-rho space
  TF1*            fFtrack_x[max_ntrack];  //### reconstructed x-tracks
  TF1*            fFtrack_y[max_ntrack];  //### reconstructed y-tracks

  //_______________________________________________________________

  TH2F*   fH2_theta_rho_x;
  TH2F*   fH2_theta_rho_y;
  //### when the entry is larger than     //####
  //### isReconTrackN in theta-rho space, //####
  //### fill to these array               //####
  int             n_max[2];
  vector<int>     ent_max[2];
  vector<int>     theta_max[2];
  vector<int>     rho_max[2];
  //### After getting the array,   //#####
  //### reconstruct tracks         //#####
  int                    ntrack_x;          //### number of x-tracks
  int                    ntrack_y;          //### number of y-tracks
  int                    ntrack[2];
  IngridTrackSummary*    ftrk[2][max_ntrack];


  //_______________________________________________________________
public:
  ReconTrack();
  ~ReconTrack();
  IngridTrackSummary* trk_x(int i){
    if(i>=0 && i < ntrack[1])
       return ftrk[1][i];
     return 0;
  }
  IngridTrackSummary* trk_y(int i){
    if(i>=0 && i < ntrack[0])
       return ftrk[0][i];
     return 0;
  }
  TF1*     Track_x(int i){
    if(i>=0 && i < ntrack_x)
      return fFtrack_x[i];
    return 0;
  }
  TF1*     Track_y(int i){
    if(i>=0 && i < ntrack_y)
      return fFtrack_y[i];
    return 0;
  }
  int N_max_x(){return n_max[1];}
  int N_max_y(){return n_max[0];}
  int Ntrack_x(){return ntrack[1];}
  int Ntrack_y(){return ntrack[0];}


  //_______________________________________________________________

  bool HoughTrans(IngridBasicReconSummary& reduc);
  bool GetMaxInHoughSpace();
  bool GetTrack(IngridBasicReconSummary& reduc);
  bool isMu(IngridTrackSummary& trk){};

  //_______________________________________________________________
  int  NHitPln(IngridTrackSummary& trk);
  bool FitTrack(IngridTrackSummary& trk);

  //_______________________________________________________________

  void Draw_Hough(TCanvas& c1);  
  void debug();
  void Reset();    
  
};

#endif

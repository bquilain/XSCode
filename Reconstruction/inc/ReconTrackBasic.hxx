#ifndef _RECONTRACKBASIC_HXX
#define _RECONTRACKBASIC_HXX

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
//#include "setup.hxx"

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

const static int max_ntrack       =   5;
const static int TolOfPoint       =   1;
const static int TolOfTrkMatching =   1;
const static int TolOfFirstPoint  =   6;
const static int TolVetoFit       =   5; //cm
const static int isTrkNHit        =   3; 
const static int fidcosmicch[11][2]
={ {11,10}, //plane 0, no fiducial
   {11,10}, //plane 1  no fiducial
   { 4,19}, //plane 2
   { 2,21}, //plane 3
   { 1,22}, //plane 4
   {11,10}, //plane 5
   { 1,22}, //plane 6
   { 2,21}, //plane 7
   { 4,19}, //plane 8
   {11,10}, //plane 9
   {11,10}  //plane 10

};
const static float xyerror        = 0.5*ScintiWidth/sqrt(3);
//_______________________________________________________________

class ReconTrackBasic{
private:
  IngridTrackSummary*       ftrkxy[2][max_ntrack];
  IngridTrackSummary*         ftrk[max_ntrack];
  IngridHitSummary*           hit[INGRIDHIT_MAXHITS];

  int                             mod;
  int                           nhits;
  int                           nhitsxy[2];
  int                        md_pln[2]; // most downstream pln
  int                       md_pln2[2]; // next most downstream pln
  vector<pair<int, int> > trk_point[2]; // <plane, ch> **[view]
  int                         n_trk[2];
  int                          n_match;

  bool                    ftplistrk[2];
  bool                   uvetoistrk[2];
  bool                  uvetohashit[2];
  bool                   uedgeistrk[2];

  bool                       newfid[2]; // new fiducial volume 2010/5/12
  bool                 newfidcosmic[2]; // new fiducial volume for cosmic


  float                         vpe[2];

  int                      startpln[2];
  int                      starttpl[2];
  int                       startch[2];
  int                        endpln[2];
  int                        endtpl[2];
  int                         endch[2];

public:
  int                      NTrackhit(int v){return trk_point[v].size();};
  int                      NHit     (int v){return nhitsxy[v];};
  int                      StartPln (int v){return startpln[v];};
  int                      StartTPL (int v){return starttpl[v];};
  int                      StartCh  (int v){return startch [v];};
  int                      EndPln   (int v){return endpln  [v];};
  int                      EndTPL   (int v){return endtpl  [v];};
  int                      EndCh    (int v){return endch   [v];};

  bool                     startisveto[2];
  bool                       endisveto[2];

  ReconTrackBasic();
  ~ReconTrackBasic();

  //bool SetHit(TRef inghit[], int nhit);
  bool SetHit(TRef inghit[INGRIDHIT_MAXHITS], int nhit);

  bool Clear();
  bool ReconTrack      ( bool cosmic );
  bool ReconTrackXY    ( int v , bool cosmic);

  int  Last3TPL        ( int v ); //### if this function fails, return velue = -1;
                                  //### if sucsses, return value is 
                                  //### current most upstream trk layer 
  bool AddTrkPoint     ( int v, int cpln );
  bool AddVETOPoint    ( int v );
  bool AddDownStreamPoint  ( int v, int cpln );//### 2010/5/17, when study VETO efficiency with cosmic
  bool Point2Trk       ( int v );
  bool TrkMatching     (  );
  bool MargeTrkXY      ( int ntrk0, int ntrk1 );

  bool FirstTPLisTrack ( int v );
  bool UVETOisTrack    ( int v ); //### Not Ready
  bool UEdgeisTrack    ( int v );
  bool UVETOhasHit     ( int v );
  bool FitTrack        ( IngridTrackSummary& trk );
  int  GetInitialZ     ( int v );
  int  GetInitialXY    ( int v );
  int  GetLastZ        ( int v );
  int  GetLastXY       ( int v );

  vector<int>         fTdcHit  (int view, int pln);
  vector<float>       fPeHit   (int view, int pln);

  IngridTrackSummary* GetTrack (int v   , int i) const;
  bool                is1stTPL (int v)const {return ftplistrk[v];}
  bool                isUVETO  (int v)const {return uvetohashit[v];}
  bool                Newfidcosmic  (int v)const {return newfidcosmic[v];}
  //bool FillTrkSummary();
  //bool MakeTrack();
  float               Vpe(int v)const {return vpe[v];}

  void Print();
};

#endif

#ifndef PMANASUMMARY_H
#define PMANASUMMARY_H
#include <iostream>
#include "TObject.h"
#include "TClonesArray.h"

#include "TRefArray.h"
#include "TRef.h"
#include "vector"
using namespace std;

#include "IngridHitSummary.h"
#include "IngridTrackSummary.h"


#define INGRIDHIT_MAXHITS 1000  //## Temporary
#define INGRIDRECON_MAXTRACKS 100 //## Temporary

//......................................................................
class PMAnaSummary : public TObject{
 public:

  int       Ntrack;
  int       Ningtrack;
  float     clstime;      // time of cluster defined by most high p.e.
  float     clstimecorr;  // time after correction by measurement time of CT
  float     exptime;      // diff. from expected time
  int       hitcyc;       //
  bool      ontime;       //
  bool      vetowtracking; // Upstream VETO
  bool      edgewtracking; // Fiducial CUT

  vector<Float_t>  x;
  vector<Float_t>  y;
  vector<Float_t>  z;
  vector<Float_t>  zx;
  vector<Float_t>  zy;

  vector<Int_t>    startxpln;
  vector<Int_t>    startypln;
  vector<Float_t>    startxch;
  vector<Float_t>    startych;
  vector<Int_t>    endxpln;
  vector<Int_t>    endypln;
  vector<Float_t>    endxch;
  vector<Float_t>    endych;
  vector<Float_t>  thetax;
  vector<Float_t>  thetay;
  vector<Float_t>  angle;
  vector<Int_t>    ing_startmod;
  vector<Int_t>    ing_endmod;
  vector<Int_t>    ing_startpln;
  vector<Int_t>    ing_endpln;
  vector<Bool_t>   ing_trk;
  vector<Bool_t>   pm_stop;
  vector<Bool_t>   ing_stop;
  vector<Float_t>  sci_range;
  vector<Float_t>  iron_range;
  vector<Int_t>    iron_pene;
  vector<Bool_t>   veto; // Upstream VETO
  vector<Bool_t>   edge; // Fiducial CUT
  vector<Int_t>    pdg;
  vector<Float_t>  mucl;
  vector<Float_t>  trkpe;
 
  
  //###########################################
  //###########################################
  PMAnaSummary();
  PMAnaSummary(const PMAnaSummary& basicsum);
  virtual ~PMAnaSummary();
  void Clear   (Option_t* option="");
  //void ClearHits   (Option_t* option="");
  void Print();
  void AddIngridHit(IngridHitSummary* hit);
  IngridHitSummary* GetIngridHit(int i) const;
  void AddIngridHitTrk(IngridHitSummary* hit, int track);
  IngridHitSummary* GetIngridHitTrk(int i, int track) const;
  void AddTrack(IngridTrackSummary* trk);
  IngridTrackSummary* GetTrack(int i) const;
  int nhits;
  int nhitTs[INGRIDRECON_MAXTRACKS];
  
 private:
  
  TRef fIngridHit[INGRIDHIT_MAXHITS];
  TRef fIngridHitTrk[INGRIDHIT_MAXHITS][INGRIDRECON_MAXTRACKS];
  TRef fTrack[INGRIDRECON_MAXTRACKS];
 public:
  int ntracks;
  int Nhits(){return nhits;}
  int NhitTs(int trk){return nhitTs[trk];}
  int Ntracks(){return ntracks;}
  ClassDef(PMAnaSummary, 1) // Ingrid Hit Summary
    };

#endif // PMANASUMMARY_H
////////////////////////////////////////////////////////////////////////

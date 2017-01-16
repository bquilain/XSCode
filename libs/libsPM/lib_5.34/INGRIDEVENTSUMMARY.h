#ifndef INGRIDEVENTSUMMARY_H
#define INGRIDEVENTSUMMARY_H
#include <iostream>
#include "TObject.h"
#include "TClonesArray.h"
#include "TRef.h"
#include "TTimeStamp.h"

#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <algorithm>
using namespace std;
#include "IngridSimHitSummary.h"
#include "IngridHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "IngridTrackSummary.h"
#include "Ingrid1stReducSummary.h"
#include "NeutInfoSummary.h"
#include "PMReconSummary.h"
#include "PMAnaSummary.h"

//......................................................................

class IngridEventSummary : public TObject
{
public:
    IngridEventSummary();
    IngridEventSummary(const IngridEventSummary& evt);
    virtual ~IngridEventSummary();

    void NewTClonesArrays();
    void Clear     (Option_t* option="");

    void Print();

    IngridSimHitSummary* AddIngridSimHit (IngridSimHitSummary* hit);
    IngridSimHitSummary* GetIngridSimHit (int i) const;
    IngridHitSummary*    AddIngridHit    (IngridHitSummary* hit);
    IngridHitSummary*    GetIngridHit    (int i) const;
    IngridHitSummary*    AddIngridModHit (IngridHitSummary* hit, int nmod, int ncyc);
    IngridHitSummary*    GetIngridModHit (int i, int nmod, int ncyc) const;

    IngridSimVertexSummary*      AddSimVertex     (IngridSimVertexSummary* hit);
    IngridSimVertexSummary*      GetSimVertex     (int i) const;
    IngridSimParticleSummary*    AddSimParticle   (IngridSimParticleSummary* hit);
    IngridSimParticleSummary*    GetSimParticle    (int i) const;

    BeamInfoSummary*     AddBeamSummary     (BeamInfoSummary*  beamsummary);
    BeamInfoSummary*     GetBeamSummary     (int i) const;

    IngridBasicReconSummary*     AddBasicRecon     (IngridBasicReconSummary* basicrecon);
    IngridBasicReconSummary*     GetBasicRecon     (int i) const;
    IngridTrackSummary*          AddTrack     (IngridTrackSummary* trk);
    IngridTrackSummary*          GetTrack     (int i) const;


    Ingrid1stReducSummary*       Add1stReduc     (Ingrid1stReducSummary* reduc); 
    Ingrid1stReducSummary*       Get1stReduc     (int i) const;

    NeutInfoSummary*       AddNeut     (NeutInfoSummary* neut); 
    NeutInfoSummary*       GetNeut     (int i) const;


    PMReconSummary*     AddPMRecon     (PMReconSummary* pmrecon);
    PMReconSummary*     GetPMRecon     (int i) const;
    PMReconSummary*     AddPMModRecon  (PMReconSummary* hit, int nmod, int ncyc, int view);
    PMReconSummary*     GetPMModRecon  (int i, int nmod, int ncyc, int view) const;
    PMAnaSummary*       AddPMAna       (PMAnaSummary* pmana);
    PMAnaSummary*       GetPMAna       (int i) const;


public:
  
    unsigned int      run;        // run number
    unsigned int    event;        // event number
    int           runmode;        // run mode
    int             trgid;        // trgger id(1:Beam, 2:Calib. , 128:Cosmic)
    unsigned char version;        // data structure version
    int              date;
    int              time;
    int           trgtime;
    int       nd280nspill;        // spill # at nd280 
                                  // = ( spill # at beam line ) 0xffff  + 1
    bool     bunch_flag[23];

    int NIngridSimHits         ()    const {return ningridsimhits;}
    int NIngridHits            ()    const {return ningridhits;}  
    int NIngridSimVertexes     ()    const {return nsimvertexes;}
    int NIngridSimParticles    ()    const {return nsimparticles;}
    int NNeutInfo              ()    const {return nneutinfos;}    
    int NIngridBeamSummarys    ()    const {return nbeamsummarys;}
    int NIngridBasicRecons     ()    const {return nbasicrecons;}
    int NPMRecons              ()    const {return npmrecons;}  
    int NPMAnas                ()    const {return npmanas;}  
    int NIngridTracks          ()    const {return ningridtracks;}  
    int NIngrid1stReduc        ()    const {return n1streducs;}
    int NIngridModHits (int nmod, int ncyc ) const {return ningridmodhits[nmod][ncyc];}  
    int NPMModRecons   (int nmod, int ncyc , int view) const {return npmmodrecons[nmod][ncyc][view];}  
    //int NIngridTracks  () const {return ningridtracks;}

private:
    int ningridsimhits;   // number of IngridSimHitSummaries  in this event
    int ningridhits;      // number of IngridHitSummaries     in this event
    int nsimparticles;    // number of IngridSimParticles     in this event
    int nsimvertexes;     // number of IngridSimVertiexes     in this event
    int nbeamsummarys;    // number of BeamInfoSummarys       in this event  
    int nbasicrecons;     // number of IngridBasicReconSummarys  in this event    
    int npmrecons;        // number of PMReconSummarys  in this event    
    int npmanas;          // number of PMAnaSummarys  in this event    
    int ningridtracks;          // number of IngridTrackSummarys  in this event    
    int n1streducs;       // number of Ingrid1stReducs  in this event    
    int nneutinfos;       // number of Ingrid1stReducs  in this event    

    int ningridmodhits[17][23];
    vector<int> nidmodhits[17*23];

    int npmmodrecons[17][23][2];
    vector<int> nidpmmodrecons[17*23*2];


    TClonesArray* fIngridSimHit; 
    TClonesArray* fIngridHit;
    TClonesArray* fSimParticle;
    TClonesArray* fSimVertex;
    TClonesArray* fBeamSummary;
    TClonesArray* fBasicRecon;
    TClonesArray* fPMRecon;
    TClonesArray* fPMAna;
    TClonesArray* f1stReduc;
    TClonesArray* fIngridTrack;
    TClonesArray* fNeutInfo;
  
    ClassDef(IngridEventSummary, 5) // DST Summary of event
        
        };

#endif // INGRIDEVENTSUMMARY_H
////////////////////////////////////////////////////////////////////////


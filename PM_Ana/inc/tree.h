//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  1 14:43:20 2012 by ROOT version 5.23/04
// from TTree tree/tree
// found on file: /home/kikawa/scraid1/done.root
//////////////////////////////////////////////////////////

#ifndef tree_h
#define tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
   const Int_t kMaxfDefaultReco = 1;
   const Int_t kMaxfDefaultReco_fIngridSimHit = 1;
   const Int_t kMaxfDefaultReco_fIngridHit = 2069;
   const Int_t kMaxfDefaultReco_fSimParticle = 1;
   const Int_t kMaxfDefaultReco_fSimVertex = 1;
   const Int_t kMaxfDefaultReco_fBeamSummary = 1;
   const Int_t kMaxfDefaultReco_fBasicRecon = 7;
   const Int_t kMaxfDefaultReco_f1stReduc = 1;
   const Int_t kMaxfDefaultReco_fIngridTrack = 8;

class tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //IngridEventSummary *fDefaultReco.;
   UInt_t          fDefaultReco_TObject_fUniqueID;
   UInt_t          fDefaultReco_TObject_fBits;
   UInt_t          fDefaultReco_run;
   UInt_t          fDefaultReco_event;
   Int_t           fDefaultReco_runmode;
   Int_t           fDefaultReco_trgid;
   UChar_t         fDefaultReco_version;
   Int_t           fDefaultReco_date;
   Int_t           fDefaultReco_time;
   Int_t           fDefaultReco_trgtime;
   Int_t           fDefaultReco_nd280nspill;
   Bool_t          fDefaultReco_bunch_flag[23];
   Int_t           fDefaultReco_ningridsimhits;
   Int_t           fDefaultReco_ningridhits;
   Int_t           fDefaultReco_nsimparticles;
   Int_t           fDefaultReco_nsimvertexes;
   Int_t           fDefaultReco_nbeamsummarys;
   Int_t           fDefaultReco_nbasicrecons;
   Int_t           fDefaultReco_ningridtracks;
   Int_t           fDefaultReco_n1streducs;
   Int_t           fDefaultReco_ningridmodhits[17][23];
   vector<int>     fDefaultReco_nidmodhits[391];
   Int_t           fDefaultReco_fIngridSimHit_;
   UInt_t          fDefaultReco_fIngridSimHit_fUniqueID[kMaxfDefaultReco_fIngridSimHit];   //[fDefaultReco.fIngridSimHit_]
   UInt_t          fDefaultReco_fIngridSimHit_fBits[kMaxfDefaultReco_fIngridSimHit];   //[fDefaultReco.fIngridSimHit_]
   Float_t         fDefaultReco_fIngridSimHit_edeposit[kMaxfDefaultReco_fIngridSimHit];   //[fDefaultReco.fIngridSimHit_]
   Int_t           fDefaultReco_fIngridSimHit_trackid[kMaxfDefaultReco_fIngridSimHit];   //[fDefaultReco.fIngridSimHit_]
   Int_t           fDefaultReco_fIngridSimHit_pdg[kMaxfDefaultReco_fIngridSimHit];   //[fDefaultReco.fIngridSimHit_]
   Int_t           fDefaultReco_fIngridHit_;
   UInt_t          fDefaultReco_fIngridHit_fUniqueID[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   UInt_t          fDefaultReco_fIngridHit_fBits[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_mod[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_cyc[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_view[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_pln[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_ch[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_adc[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_loadc[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_pe[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_lope[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_pecorr[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_vise[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_visecorr[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Long_t          fDefaultReco_fIngridHit_tdc[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_time[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_tnearhit[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_timecorr[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_xy[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Float_t         fDefaultReco_fIngridHit_z[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_addbasicrecon[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_dummy[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Bool_t          fDefaultReco_fIngridHit_gocosmic[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Bool_t          fDefaultReco_fIngridHit_hitcosmic[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   Int_t           fDefaultReco_fIngridHit_nsimhits[kMaxfDefaultReco_fIngridHit];   //[fDefaultReco.fIngridHit_]
   TRef            fDefaultReco_fIngridHit_fIngridSimHit[1][kMaxfDefaultReco_fIngridHit];
   Int_t           fDefaultReco_fSimParticle_;
   UInt_t          fDefaultReco_fSimParticle_fUniqueID[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   UInt_t          fDefaultReco_fSimParticle_fBits[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Int_t           fDefaultReco_fSimParticle_trackid[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Int_t           fDefaultReco_fSimParticle_parentid[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Int_t           fDefaultReco_fSimParticle_pdg[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Float_t         fDefaultReco_fSimParticle_momentum[kMaxfDefaultReco_fSimParticle][4];   //[fDefaultReco.fSimParticle_]
   Float_t         fDefaultReco_fSimParticle_ipos[kMaxfDefaultReco_fSimParticle][4];   //[fDefaultReco.fSimParticle_]
   Float_t         fDefaultReco_fSimParticle_fpos[kMaxfDefaultReco_fSimParticle][4];   //[fDefaultReco.fSimParticle_]
   Int_t           fDefaultReco_fSimParticle_iposflag[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Int_t           fDefaultReco_fSimParticle_fposflag[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Float_t         fDefaultReco_fSimParticle_dir[kMaxfDefaultReco_fSimParticle][3];   //[fDefaultReco.fSimParticle_]
   Float_t         fDefaultReco_fSimParticle_length[kMaxfDefaultReco_fSimParticle];   //[fDefaultReco.fSimParticle_]
   Int_t           fDefaultReco_fSimVertex_;
   UInt_t          fDefaultReco_fSimVertex_fUniqueID[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   UInt_t          fDefaultReco_fSimVertex_fBits[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fSimVertex_mod[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fSimVertex_nutype[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_numomentum[kMaxfDefaultReco_fSimVertex][3];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_nuE[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fSimVertex_targeta[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fSimVertex_targetz[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fSimVertex_targettype[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_targetmomentum[kMaxfDefaultReco_fSimVertex][3];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_vnuclini[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_pfsurf[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fSimVertex_inttype[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_x[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_y[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_z[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_lx[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_ly[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_lz[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Float_t         fDefaultReco_fSimVertex_time[kMaxfDefaultReco_fSimVertex];   //[fDefaultReco.fSimVertex_]
   Int_t           fDefaultReco_fBeamSummary_;
   UInt_t          fDefaultReco_fBeamSummary_fUniqueID[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   UInt_t          fDefaultReco_fBeamSummary_fBits[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_nurun[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_spill_flag[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_good_spill_flag[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_trg_sec[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_spillnum[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_nd280spillnum[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBeamSummary_run_type[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Double_t        fDefaultReco_fBeamSummary_target_eff[kMaxfDefaultReco_fBeamSummary][3];   //[fDefaultReco.fBeamSummary_]
   Double_t        fDefaultReco_fBeamSummary_beam_time[kMaxfDefaultReco_fBeamSummary][5][9];   //[fDefaultReco.fBeamSummary_]
   Double_t        fDefaultReco_fBeamSummary_ct_np[kMaxfDefaultReco_fBeamSummary][5][9];   //[fDefaultReco.fBeamSummary_]
   Double_t        fDefaultReco_fBeamSummary_mumon[kMaxfDefaultReco_fBeamSummary][12];   //[fDefaultReco.fBeamSummary_]
   Double_t        fDefaultReco_fBeamSummary_hct[kMaxfDefaultReco_fBeamSummary][3][5];   //[fDefaultReco.fBeamSummary_]
   Double_t        fDefaultReco_fBeamSummary_otr[kMaxfDefaultReco_fBeamSummary][13];   //[fDefaultReco.fBeamSummary_]
   Bool_t          fDefaultReco_fBeamSummary_wohorn[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Bool_t          fDefaultReco_fBeamSummary_horn1[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Bool_t          fDefaultReco_fBeamSummary_whorn[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Bool_t          fDefaultReco_fBeamSummary_cutwotr[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Bool_t          fDefaultReco_fBeamSummary_horn250[kMaxfDefaultReco_fBeamSummary];   //[fDefaultReco.fBeamSummary_]
   Int_t           fDefaultReco_fBasicRecon_;
   UInt_t          fDefaultReco_fBasicRecon_fUniqueID[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   UInt_t          fDefaultReco_fBasicRecon_fBits[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_clstime[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_CTtimecorr[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_clstimecorr[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_exptime[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_nhitclster[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_nactpln[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_actinarow[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_layerpe[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_upstreamVETO[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_upstreamedge[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_newfid[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_newfidcosmic[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_vinternal[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_hitmod[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_hitcyc[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_spill_flag[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_bunch_flag[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_ontime[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_trgbit[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_inarowbit[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_vetowtracking[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_edgewtracking[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_hastrk[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_matchtrk[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Bool_t          fDefaultReco_fBasicRecon_modfc[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_penIron[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_muE[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_nuErec[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_nhitx[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_nhity[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_ntrackhitx[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_ntrackhity[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_retracktest[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_trg_sec[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_upstreamtpl[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_vertexz[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_vertexxz[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_vertexyz[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_angle[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_thetax[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_thetay[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   vector<int>     fDefaultReco_fBasicRecon_vertexx[kMaxfDefaultReco_fBasicRecon];
   vector<int>     fDefaultReco_fBasicRecon_vertexy[kMaxfDefaultReco_fBasicRecon];
   Bool_t          fDefaultReco_fBasicRecon_horn250[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_vpe[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_startxpln[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_startypln[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_startxch[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_startych[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_endxpln[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_endypln[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_endxch[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_endych[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_vertexx_true[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_vertexy_true[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_vertexz_true[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_lverx[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_lvery[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_lverxz[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Float_t         fDefaultReco_fBasicRecon_lveryz[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_fBasicRecon_nhits[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   TRef            fDefaultReco_fBasicRecon_fIngridHit[1000][kMaxfDefaultReco_fBasicRecon];
   TRef            fDefaultReco_fBasicRecon_fTrack[10][kMaxfDefaultReco_fBasicRecon];
   Int_t           fDefaultReco_fBasicRecon_ntracks[kMaxfDefaultReco_fBasicRecon];   //[fDefaultReco.fBasicRecon_]
   Int_t           fDefaultReco_f1stReduc_;
   UInt_t          fDefaultReco_f1stReduc_fUniqueID[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   UInt_t          fDefaultReco_f1stReduc_fBits[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Int_t           fDefaultReco_f1stReduc_hitmod[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Int_t           fDefaultReco_f1stReduc_hitcyc[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   UInt_t          fDefaultReco_f1stReduc_xhitbit[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   UInt_t          fDefaultReco_f1stReduc_yhitbit[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Int_t           fDefaultReco_f1stReduc_nhitxlyr[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Int_t           fDefaultReco_f1stReduc_nhitylyr[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Float_t         fDefaultReco_f1stReduc_xtotpe[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Float_t         fDefaultReco_f1stReduc_ytotpe[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Bool_t          fDefaultReco_f1stReduc_xtracklike[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   Bool_t          fDefaultReco_f1stReduc_ytracklike[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   vector<int>     fDefaultReco_f1stReduc_hity[kMaxfDefaultReco_f1stReduc];
   vector<int>     fDefaultReco_f1stReduc_hitx[kMaxfDefaultReco_f1stReduc];
   vector<int>     fDefaultReco_f1stReduc_hitxz[kMaxfDefaultReco_f1stReduc];
   vector<int>     fDefaultReco_f1stReduc_hityz[kMaxfDefaultReco_f1stReduc];
   Int_t           fDefaultReco_f1stReduc_nhits[kMaxfDefaultReco_f1stReduc];   //[fDefaultReco.f1stReduc_]
   TRef            fDefaultReco_f1stReduc_fIngridHit[1000][kMaxfDefaultReco_f1stReduc];
   Int_t           fDefaultReco_fIngridTrack_;
   UInt_t          fDefaultReco_fIngridTrack_fUniqueID[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   UInt_t          fDefaultReco_fIngridTrack_fBits[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_vtxi[kMaxfDefaultReco_fIngridTrack][3];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_vtxf[kMaxfDefaultReco_fIngridTrack][3];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_length[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_ekin[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_tx[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_ty[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_etx[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_ety[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_ex0[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_ey0[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_covx[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_covy[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_chi2x[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_chi2y[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_btheta[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_bphi[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Int_t           fDefaultReco_fIngridTrack_mrdhitid[kMaxfDefaultReco_fIngridTrack][2];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_mucl[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Float_t         fDefaultReco_fIngridTrack_vpe[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Int_t           fDefaultReco_fIngridTrack_view[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Int_t           fDefaultReco_fIngridTrack_nhits[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Int_t           fDefaultReco_fIngridTrack_nsimparticles[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   Int_t           fDefaultReco_fIngridTrack_nsimemshowers[kMaxfDefaultReco_fIngridTrack];   //[fDefaultReco.fIngridTrack_]
   TRef            fDefaultReco_fIngridTrack_fIngridHit[1000][kMaxfDefaultReco_fIngridTrack];
   TRef            fDefaultReco_fIngridTrack_fSimParticle[10][kMaxfDefaultReco_fIngridTrack];

   // List of branches
   TBranch        *b_fDefaultReco_TObject_fUniqueID;   //!
   TBranch        *b_fDefaultReco_TObject_fBits;   //!
   TBranch        *b_fDefaultReco_run;   //!
   TBranch        *b_fDefaultReco_event;   //!
   TBranch        *b_fDefaultReco_runmode;   //!
   TBranch        *b_fDefaultReco_trgid;   //!
   TBranch        *b_fDefaultReco_version;   //!
   TBranch        *b_fDefaultReco_date;   //!
   TBranch        *b_fDefaultReco_time;   //!
   TBranch        *b_fDefaultReco_trgtime;   //!
   TBranch        *b_fDefaultReco_nd280nspill;   //!
   TBranch        *b_fDefaultReco_bunch_flag;   //!
   TBranch        *b_fDefaultReco_ningridsimhits;   //!
   TBranch        *b_fDefaultReco_ningridhits;   //!
   TBranch        *b_fDefaultReco_nsimparticles;   //!
   TBranch        *b_fDefaultReco_nsimvertexes;   //!
   TBranch        *b_fDefaultReco_nbeamsummarys;   //!
   TBranch        *b_fDefaultReco_nbasicrecons;   //!
   TBranch        *b_fDefaultReco_ningridtracks;   //!
   TBranch        *b_fDefaultReco_n1streducs;   //!
   TBranch        *b_fDefaultReco_ningridmodhits;   //!
   TBranch        *b_fDefaultReco_nidmodhits;   //!
   TBranch        *b_fDefaultReco_fIngridSimHit_;   //!
   TBranch        *b_fDefaultReco_fIngridSimHit_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fIngridSimHit_fBits;   //!
   TBranch        *b_fDefaultReco_fIngridSimHit_edeposit;   //!
   TBranch        *b_fDefaultReco_fIngridSimHit_trackid;   //!
   TBranch        *b_fDefaultReco_fIngridSimHit_pdg;   //!
   TBranch        *b_fDefaultReco_fIngridHit_;   //!
   TBranch        *b_fDefaultReco_fIngridHit_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fIngridHit_fBits;   //!
   TBranch        *b_fDefaultReco_fIngridHit_mod;   //!
   TBranch        *b_fDefaultReco_fIngridHit_cyc;   //!
   TBranch        *b_fDefaultReco_fIngridHit_view;   //!
   TBranch        *b_fDefaultReco_fIngridHit_pln;   //!
   TBranch        *b_fDefaultReco_fIngridHit_ch;   //!
   TBranch        *b_fDefaultReco_fIngridHit_adc;   //!
   TBranch        *b_fDefaultReco_fIngridHit_loadc;   //!
   TBranch        *b_fDefaultReco_fIngridHit_pe;   //!
   TBranch        *b_fDefaultReco_fIngridHit_lope;   //!
   TBranch        *b_fDefaultReco_fIngridHit_pecorr;   //!
   TBranch        *b_fDefaultReco_fIngridHit_vise;   //!
   TBranch        *b_fDefaultReco_fIngridHit_visecorr;   //!
   TBranch        *b_fDefaultReco_fIngridHit_tdc;   //!
   TBranch        *b_fDefaultReco_fIngridHit_time;   //!
   TBranch        *b_fDefaultReco_fIngridHit_tnearhit;   //!
   TBranch        *b_fDefaultReco_fIngridHit_timecorr;   //!
   TBranch        *b_fDefaultReco_fIngridHit_xy;   //!
   TBranch        *b_fDefaultReco_fIngridHit_z;   //!
   TBranch        *b_fDefaultReco_fIngridHit_addbasicrecon;   //!
   TBranch        *b_fDefaultReco_fIngridHit_dummy;   //!
   TBranch        *b_fDefaultReco_fIngridHit_gocosmic;   //!
   TBranch        *b_fDefaultReco_fIngridHit_hitcosmic;   //!
   TBranch        *b_fDefaultReco_fIngridHit_nsimhits;   //!
   TBranch        *b_fDefaultReco_fIngridHit_fIngridSimHit;   //!
   TBranch        *b_fDefaultReco_fSimParticle_;   //!
   TBranch        *b_fDefaultReco_fSimParticle_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fSimParticle_fBits;   //!
   TBranch        *b_fDefaultReco_fSimParticle_trackid;   //!
   TBranch        *b_fDefaultReco_fSimParticle_parentid;   //!
   TBranch        *b_fDefaultReco_fSimParticle_pdg;   //!
   TBranch        *b_fDefaultReco_fSimParticle_momentum;   //!
   TBranch        *b_fDefaultReco_fSimParticle_ipos;   //!
   TBranch        *b_fDefaultReco_fSimParticle_fpos;   //!
   TBranch        *b_fDefaultReco_fSimParticle_iposflag;   //!
   TBranch        *b_fDefaultReco_fSimParticle_fposflag;   //!
   TBranch        *b_fDefaultReco_fSimParticle_dir;   //!
   TBranch        *b_fDefaultReco_fSimParticle_length;   //!
   TBranch        *b_fDefaultReco_fSimVertex_;   //!
   TBranch        *b_fDefaultReco_fSimVertex_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fSimVertex_fBits;   //!
   TBranch        *b_fDefaultReco_fSimVertex_mod;   //!
   TBranch        *b_fDefaultReco_fSimVertex_nutype;   //!
   TBranch        *b_fDefaultReco_fSimVertex_numomentum;   //!
   TBranch        *b_fDefaultReco_fSimVertex_nuE;   //!
   TBranch        *b_fDefaultReco_fSimVertex_targeta;   //!
   TBranch        *b_fDefaultReco_fSimVertex_targetz;   //!
   TBranch        *b_fDefaultReco_fSimVertex_targettype;   //!
   TBranch        *b_fDefaultReco_fSimVertex_targetmomentum;   //!
   TBranch        *b_fDefaultReco_fSimVertex_vnuclini;   //!
   TBranch        *b_fDefaultReco_fSimVertex_pfsurf;   //!
   TBranch        *b_fDefaultReco_fSimVertex_inttype;   //!
   TBranch        *b_fDefaultReco_fSimVertex_x;   //!
   TBranch        *b_fDefaultReco_fSimVertex_y;   //!
   TBranch        *b_fDefaultReco_fSimVertex_z;   //!
   TBranch        *b_fDefaultReco_fSimVertex_lx;   //!
   TBranch        *b_fDefaultReco_fSimVertex_ly;   //!
   TBranch        *b_fDefaultReco_fSimVertex_lz;   //!
   TBranch        *b_fDefaultReco_fSimVertex_time;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_fBits;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_nurun;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_spill_flag;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_good_spill_flag;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_trg_sec;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_spillnum;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_nd280spillnum;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_run_type;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_target_eff;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_beam_time;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_ct_np;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_mumon;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_hct;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_otr;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_wohorn;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_horn1;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_whorn;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_cutwotr;   //!
   TBranch        *b_fDefaultReco_fBeamSummary_horn250;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_fBits;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_clstime;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_CTtimecorr;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_clstimecorr;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_exptime;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_nhitclster;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_nactpln;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_actinarow;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_layerpe;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_upstreamVETO;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_upstreamedge;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_newfid;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_newfidcosmic;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vinternal;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_hitmod;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_hitcyc;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_spill_flag;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_bunch_flag;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_ontime;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_trgbit;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_inarowbit;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vetowtracking;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_edgewtracking;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_hastrk;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_matchtrk;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_modfc;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_penIron;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_muE;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_nuErec;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_nhitx;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_nhity;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_ntrackhitx;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_ntrackhity;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_retracktest;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_trg_sec;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_upstreamtpl;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexz;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexxz;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexyz;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_angle;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_thetax;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_thetay;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexx;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexy;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_horn250;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vpe;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_startxpln;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_startypln;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_startxch;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_startych;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_endxpln;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_endypln;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_endxch;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_endych;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexx_true;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexy_true;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_vertexz_true;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_lverx;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_lvery;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_lverxz;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_lveryz;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_nhits;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_fIngridHit;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_fTrack;   //!
   TBranch        *b_fDefaultReco_fBasicRecon_ntracks;   //!
   TBranch        *b_fDefaultReco_f1stReduc_;   //!
   TBranch        *b_fDefaultReco_f1stReduc_fUniqueID;   //!
   TBranch        *b_fDefaultReco_f1stReduc_fBits;   //!
   TBranch        *b_fDefaultReco_f1stReduc_hitmod;   //!
   TBranch        *b_fDefaultReco_f1stReduc_hitcyc;   //!
   TBranch        *b_fDefaultReco_f1stReduc_xhitbit;   //!
   TBranch        *b_fDefaultReco_f1stReduc_yhitbit;   //!
   TBranch        *b_fDefaultReco_f1stReduc_nhitxlyr;   //!
   TBranch        *b_fDefaultReco_f1stReduc_nhitylyr;   //!
   TBranch        *b_fDefaultReco_f1stReduc_xtotpe;   //!
   TBranch        *b_fDefaultReco_f1stReduc_ytotpe;   //!
   TBranch        *b_fDefaultReco_f1stReduc_xtracklike;   //!
   TBranch        *b_fDefaultReco_f1stReduc_ytracklike;   //!
   TBranch        *b_fDefaultReco_f1stReduc_hity;   //!
   TBranch        *b_fDefaultReco_f1stReduc_hitx;   //!
   TBranch        *b_fDefaultReco_f1stReduc_hitxz;   //!
   TBranch        *b_fDefaultReco_f1stReduc_hityz;   //!
   TBranch        *b_fDefaultReco_f1stReduc_nhits;   //!
   TBranch        *b_fDefaultReco_f1stReduc_fIngridHit;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_fUniqueID;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_fBits;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_vtxi;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_vtxf;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_length;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_ekin;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_tx;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_ty;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_etx;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_ety;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_ex0;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_ey0;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_covx;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_covy;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_chi2x;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_chi2y;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_btheta;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_bphi;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_mrdhitid;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_mucl;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_vpe;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_view;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_nhits;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_nsimparticles;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_nsimemshowers;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_fIngridHit;   //!
   TBranch        *b_fDefaultReco_fIngridTrack_fSimParticle;   //!

   tree(TTree *tree=0);
   virtual ~tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tree_cxx
tree::tree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/kikawa/scraid1/done.root");
      if (!f) {
         f = new TFile("/home/kikawa/scraid1/done.root");
      }
      tree = (TTree*)gDirectory->Get("tree");

   }
   Init(tree);
}

tree::~tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fDefaultReco.TObject.fUniqueID", &fDefaultReco_TObject_fUniqueID, &b_fDefaultReco_TObject_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.TObject.fBits", &fDefaultReco_TObject_fBits, &b_fDefaultReco_TObject_fBits);
   fChain->SetBranchAddress("fDefaultReco.run", &fDefaultReco_run, &b_fDefaultReco_run);
   fChain->SetBranchAddress("fDefaultReco.event", &fDefaultReco_event, &b_fDefaultReco_event);
   fChain->SetBranchAddress("fDefaultReco.runmode", &fDefaultReco_runmode, &b_fDefaultReco_runmode);
   fChain->SetBranchAddress("fDefaultReco.trgid", &fDefaultReco_trgid, &b_fDefaultReco_trgid);
   fChain->SetBranchAddress("fDefaultReco.version", &fDefaultReco_version, &b_fDefaultReco_version);
   fChain->SetBranchAddress("fDefaultReco.date", &fDefaultReco_date, &b_fDefaultReco_date);
   fChain->SetBranchAddress("fDefaultReco.time", &fDefaultReco_time, &b_fDefaultReco_time);
   fChain->SetBranchAddress("fDefaultReco.trgtime", &fDefaultReco_trgtime, &b_fDefaultReco_trgtime);
   fChain->SetBranchAddress("fDefaultReco.nd280nspill", &fDefaultReco_nd280nspill, &b_fDefaultReco_nd280nspill);
   fChain->SetBranchAddress("fDefaultReco.bunch_flag[23]", fDefaultReco_bunch_flag, &b_fDefaultReco_bunch_flag);
   fChain->SetBranchAddress("fDefaultReco.ningridsimhits", &fDefaultReco_ningridsimhits, &b_fDefaultReco_ningridsimhits);
   fChain->SetBranchAddress("fDefaultReco.ningridhits", &fDefaultReco_ningridhits, &b_fDefaultReco_ningridhits);
   fChain->SetBranchAddress("fDefaultReco.nsimparticles", &fDefaultReco_nsimparticles, &b_fDefaultReco_nsimparticles);
   fChain->SetBranchAddress("fDefaultReco.nsimvertexes", &fDefaultReco_nsimvertexes, &b_fDefaultReco_nsimvertexes);
   fChain->SetBranchAddress("fDefaultReco.nbeamsummarys", &fDefaultReco_nbeamsummarys, &b_fDefaultReco_nbeamsummarys);
   fChain->SetBranchAddress("fDefaultReco.nbasicrecons", &fDefaultReco_nbasicrecons, &b_fDefaultReco_nbasicrecons);
   fChain->SetBranchAddress("fDefaultReco.ningridtracks", &fDefaultReco_ningridtracks, &b_fDefaultReco_ningridtracks);
   fChain->SetBranchAddress("fDefaultReco.n1streducs", &fDefaultReco_n1streducs, &b_fDefaultReco_n1streducs);
   fChain->SetBranchAddress("fDefaultReco.ningridmodhits[17][23]", fDefaultReco_ningridmodhits, &b_fDefaultReco_ningridmodhits);
   fChain->SetBranchAddress("fDefaultReco.nidmodhits[391]", fDefaultReco_nidmodhits, &b_fDefaultReco_nidmodhits);
   fChain->SetBranchAddress("fDefaultReco.fIngridSimHit", &fDefaultReco_fIngridSimHit_, &b_fDefaultReco_fIngridSimHit_);
   fChain->SetBranchAddress("fDefaultReco.fIngridSimHit.fUniqueID", &fDefaultReco_fIngridSimHit_fUniqueID, &b_fDefaultReco_fIngridSimHit_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fIngridSimHit.fBits", &fDefaultReco_fIngridSimHit_fBits, &b_fDefaultReco_fIngridSimHit_fBits);
   fChain->SetBranchAddress("fDefaultReco.fIngridSimHit.edeposit", &fDefaultReco_fIngridSimHit_edeposit, &b_fDefaultReco_fIngridSimHit_edeposit);
   fChain->SetBranchAddress("fDefaultReco.fIngridSimHit.trackid", &fDefaultReco_fIngridSimHit_trackid, &b_fDefaultReco_fIngridSimHit_trackid);
   fChain->SetBranchAddress("fDefaultReco.fIngridSimHit.pdg", &fDefaultReco_fIngridSimHit_pdg, &b_fDefaultReco_fIngridSimHit_pdg);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit", &fDefaultReco_fIngridHit_, &b_fDefaultReco_fIngridHit_);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.fUniqueID", fDefaultReco_fIngridHit_fUniqueID, &b_fDefaultReco_fIngridHit_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.fBits", fDefaultReco_fIngridHit_fBits, &b_fDefaultReco_fIngridHit_fBits);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.mod", fDefaultReco_fIngridHit_mod, &b_fDefaultReco_fIngridHit_mod);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.cyc", fDefaultReco_fIngridHit_cyc, &b_fDefaultReco_fIngridHit_cyc);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.view", fDefaultReco_fIngridHit_view, &b_fDefaultReco_fIngridHit_view);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.pln", fDefaultReco_fIngridHit_pln, &b_fDefaultReco_fIngridHit_pln);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.ch", fDefaultReco_fIngridHit_ch, &b_fDefaultReco_fIngridHit_ch);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.adc", fDefaultReco_fIngridHit_adc, &b_fDefaultReco_fIngridHit_adc);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.loadc", fDefaultReco_fIngridHit_loadc, &b_fDefaultReco_fIngridHit_loadc);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.pe", fDefaultReco_fIngridHit_pe, &b_fDefaultReco_fIngridHit_pe);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.lope", fDefaultReco_fIngridHit_lope, &b_fDefaultReco_fIngridHit_lope);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.pecorr", fDefaultReco_fIngridHit_pecorr, &b_fDefaultReco_fIngridHit_pecorr);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.vise", fDefaultReco_fIngridHit_vise, &b_fDefaultReco_fIngridHit_vise);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.visecorr", fDefaultReco_fIngridHit_visecorr, &b_fDefaultReco_fIngridHit_visecorr);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.tdc", fDefaultReco_fIngridHit_tdc, &b_fDefaultReco_fIngridHit_tdc);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.time", fDefaultReco_fIngridHit_time, &b_fDefaultReco_fIngridHit_time);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.tnearhit", fDefaultReco_fIngridHit_tnearhit, &b_fDefaultReco_fIngridHit_tnearhit);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.timecorr", fDefaultReco_fIngridHit_timecorr, &b_fDefaultReco_fIngridHit_timecorr);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.xy", fDefaultReco_fIngridHit_xy, &b_fDefaultReco_fIngridHit_xy);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.z", fDefaultReco_fIngridHit_z, &b_fDefaultReco_fIngridHit_z);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.addbasicrecon", fDefaultReco_fIngridHit_addbasicrecon, &b_fDefaultReco_fIngridHit_addbasicrecon);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.dummy", fDefaultReco_fIngridHit_dummy, &b_fDefaultReco_fIngridHit_dummy);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.gocosmic", fDefaultReco_fIngridHit_gocosmic, &b_fDefaultReco_fIngridHit_gocosmic);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.hitcosmic", fDefaultReco_fIngridHit_hitcosmic, &b_fDefaultReco_fIngridHit_hitcosmic);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.nsimhits", fDefaultReco_fIngridHit_nsimhits, &b_fDefaultReco_fIngridHit_nsimhits);
   fChain->SetBranchAddress("fDefaultReco.fIngridHit.fIngridSimHit[1]", fDefaultReco_fIngridHit_fIngridSimHit, &b_fDefaultReco_fIngridHit_fIngridSimHit);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle", &fDefaultReco_fSimParticle_, &b_fDefaultReco_fSimParticle_);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.fUniqueID", &fDefaultReco_fSimParticle_fUniqueID, &b_fDefaultReco_fSimParticle_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.fBits", &fDefaultReco_fSimParticle_fBits, &b_fDefaultReco_fSimParticle_fBits);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.trackid", &fDefaultReco_fSimParticle_trackid, &b_fDefaultReco_fSimParticle_trackid);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.parentid", &fDefaultReco_fSimParticle_parentid, &b_fDefaultReco_fSimParticle_parentid);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.pdg", &fDefaultReco_fSimParticle_pdg, &b_fDefaultReco_fSimParticle_pdg);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.momentum[4]", &fDefaultReco_fSimParticle_momentum, &b_fDefaultReco_fSimParticle_momentum);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.ipos[4]", &fDefaultReco_fSimParticle_ipos, &b_fDefaultReco_fSimParticle_ipos);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.fpos[4]", &fDefaultReco_fSimParticle_fpos, &b_fDefaultReco_fSimParticle_fpos);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.iposflag", &fDefaultReco_fSimParticle_iposflag, &b_fDefaultReco_fSimParticle_iposflag);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.fposflag", &fDefaultReco_fSimParticle_fposflag, &b_fDefaultReco_fSimParticle_fposflag);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.dir[3]", &fDefaultReco_fSimParticle_dir, &b_fDefaultReco_fSimParticle_dir);
   fChain->SetBranchAddress("fDefaultReco.fSimParticle.length", &fDefaultReco_fSimParticle_length, &b_fDefaultReco_fSimParticle_length);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex", &fDefaultReco_fSimVertex_, &b_fDefaultReco_fSimVertex_);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.fUniqueID", &fDefaultReco_fSimVertex_fUniqueID, &b_fDefaultReco_fSimVertex_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.fBits", &fDefaultReco_fSimVertex_fBits, &b_fDefaultReco_fSimVertex_fBits);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.mod", &fDefaultReco_fSimVertex_mod, &b_fDefaultReco_fSimVertex_mod);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.nutype", &fDefaultReco_fSimVertex_nutype, &b_fDefaultReco_fSimVertex_nutype);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.numomentum[3]", &fDefaultReco_fSimVertex_numomentum, &b_fDefaultReco_fSimVertex_numomentum);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.nuE", &fDefaultReco_fSimVertex_nuE, &b_fDefaultReco_fSimVertex_nuE);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.targeta", &fDefaultReco_fSimVertex_targeta, &b_fDefaultReco_fSimVertex_targeta);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.targetz", &fDefaultReco_fSimVertex_targetz, &b_fDefaultReco_fSimVertex_targetz);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.targettype", &fDefaultReco_fSimVertex_targettype, &b_fDefaultReco_fSimVertex_targettype);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.targetmomentum[3]", &fDefaultReco_fSimVertex_targetmomentum, &b_fDefaultReco_fSimVertex_targetmomentum);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.vnuclini", &fDefaultReco_fSimVertex_vnuclini, &b_fDefaultReco_fSimVertex_vnuclini);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.pfsurf", &fDefaultReco_fSimVertex_pfsurf, &b_fDefaultReco_fSimVertex_pfsurf);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.inttype", &fDefaultReco_fSimVertex_inttype, &b_fDefaultReco_fSimVertex_inttype);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.x", &fDefaultReco_fSimVertex_x, &b_fDefaultReco_fSimVertex_x);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.y", &fDefaultReco_fSimVertex_y, &b_fDefaultReco_fSimVertex_y);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.z", &fDefaultReco_fSimVertex_z, &b_fDefaultReco_fSimVertex_z);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.lx", &fDefaultReco_fSimVertex_lx, &b_fDefaultReco_fSimVertex_lx);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.ly", &fDefaultReco_fSimVertex_ly, &b_fDefaultReco_fSimVertex_ly);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.lz", &fDefaultReco_fSimVertex_lz, &b_fDefaultReco_fSimVertex_lz);
   fChain->SetBranchAddress("fDefaultReco.fSimVertex.time", &fDefaultReco_fSimVertex_time, &b_fDefaultReco_fSimVertex_time);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary", &fDefaultReco_fBeamSummary_, &b_fDefaultReco_fBeamSummary_);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.fUniqueID", fDefaultReco_fBeamSummary_fUniqueID, &b_fDefaultReco_fBeamSummary_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.fBits", fDefaultReco_fBeamSummary_fBits, &b_fDefaultReco_fBeamSummary_fBits);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.nurun", fDefaultReco_fBeamSummary_nurun, &b_fDefaultReco_fBeamSummary_nurun);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.spill_flag", fDefaultReco_fBeamSummary_spill_flag, &b_fDefaultReco_fBeamSummary_spill_flag);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.good_spill_flag", fDefaultReco_fBeamSummary_good_spill_flag, &b_fDefaultReco_fBeamSummary_good_spill_flag);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.trg_sec", fDefaultReco_fBeamSummary_trg_sec, &b_fDefaultReco_fBeamSummary_trg_sec);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.spillnum", fDefaultReco_fBeamSummary_spillnum, &b_fDefaultReco_fBeamSummary_spillnum);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.nd280spillnum", fDefaultReco_fBeamSummary_nd280spillnum, &b_fDefaultReco_fBeamSummary_nd280spillnum);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.run_type", fDefaultReco_fBeamSummary_run_type, &b_fDefaultReco_fBeamSummary_run_type);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.target_eff[3]", fDefaultReco_fBeamSummary_target_eff, &b_fDefaultReco_fBeamSummary_target_eff);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.beam_time[5][9]", fDefaultReco_fBeamSummary_beam_time, &b_fDefaultReco_fBeamSummary_beam_time);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.ct_np[5][9]", fDefaultReco_fBeamSummary_ct_np, &b_fDefaultReco_fBeamSummary_ct_np);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.mumon[12]", fDefaultReco_fBeamSummary_mumon, &b_fDefaultReco_fBeamSummary_mumon);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.hct[3][5]", fDefaultReco_fBeamSummary_hct, &b_fDefaultReco_fBeamSummary_hct);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.otr[13]", fDefaultReco_fBeamSummary_otr, &b_fDefaultReco_fBeamSummary_otr);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.wohorn", fDefaultReco_fBeamSummary_wohorn, &b_fDefaultReco_fBeamSummary_wohorn);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.horn1", fDefaultReco_fBeamSummary_horn1, &b_fDefaultReco_fBeamSummary_horn1);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.whorn", fDefaultReco_fBeamSummary_whorn, &b_fDefaultReco_fBeamSummary_whorn);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.cutwotr", fDefaultReco_fBeamSummary_cutwotr, &b_fDefaultReco_fBeamSummary_cutwotr);
   fChain->SetBranchAddress("fDefaultReco.fBeamSummary.horn250", fDefaultReco_fBeamSummary_horn250, &b_fDefaultReco_fBeamSummary_horn250);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon", &fDefaultReco_fBasicRecon_, &b_fDefaultReco_fBasicRecon_);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.fUniqueID", fDefaultReco_fBasicRecon_fUniqueID, &b_fDefaultReco_fBasicRecon_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.fBits", fDefaultReco_fBasicRecon_fBits, &b_fDefaultReco_fBasicRecon_fBits);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.clstime", fDefaultReco_fBasicRecon_clstime, &b_fDefaultReco_fBasicRecon_clstime);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.CTtimecorr", fDefaultReco_fBasicRecon_CTtimecorr, &b_fDefaultReco_fBasicRecon_CTtimecorr);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.clstimecorr", fDefaultReco_fBasicRecon_clstimecorr, &b_fDefaultReco_fBasicRecon_clstimecorr);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.exptime", fDefaultReco_fBasicRecon_exptime, &b_fDefaultReco_fBasicRecon_exptime);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.nhitclster", fDefaultReco_fBasicRecon_nhitclster, &b_fDefaultReco_fBasicRecon_nhitclster);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.nactpln", fDefaultReco_fBasicRecon_nactpln, &b_fDefaultReco_fBasicRecon_nactpln);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.actinarow", fDefaultReco_fBasicRecon_actinarow, &b_fDefaultReco_fBasicRecon_actinarow);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.layerpe", fDefaultReco_fBasicRecon_layerpe, &b_fDefaultReco_fBasicRecon_layerpe);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.upstreamVETO", fDefaultReco_fBasicRecon_upstreamVETO, &b_fDefaultReco_fBasicRecon_upstreamVETO);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.upstreamedge", fDefaultReco_fBasicRecon_upstreamedge, &b_fDefaultReco_fBasicRecon_upstreamedge);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.newfid", fDefaultReco_fBasicRecon_newfid, &b_fDefaultReco_fBasicRecon_newfid);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.newfidcosmic", fDefaultReco_fBasicRecon_newfidcosmic, &b_fDefaultReco_fBasicRecon_newfidcosmic);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vinternal", fDefaultReco_fBasicRecon_vinternal, &b_fDefaultReco_fBasicRecon_vinternal);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.hitmod", fDefaultReco_fBasicRecon_hitmod, &b_fDefaultReco_fBasicRecon_hitmod);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.hitcyc", fDefaultReco_fBasicRecon_hitcyc, &b_fDefaultReco_fBasicRecon_hitcyc);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.spill_flag", fDefaultReco_fBasicRecon_spill_flag, &b_fDefaultReco_fBasicRecon_spill_flag);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.bunch_flag", fDefaultReco_fBasicRecon_bunch_flag, &b_fDefaultReco_fBasicRecon_bunch_flag);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.ontime", fDefaultReco_fBasicRecon_ontime, &b_fDefaultReco_fBasicRecon_ontime);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.trgbit", fDefaultReco_fBasicRecon_trgbit, &b_fDefaultReco_fBasicRecon_trgbit);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.inarowbit", fDefaultReco_fBasicRecon_inarowbit, &b_fDefaultReco_fBasicRecon_inarowbit);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vetowtracking", fDefaultReco_fBasicRecon_vetowtracking, &b_fDefaultReco_fBasicRecon_vetowtracking);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.edgewtracking", fDefaultReco_fBasicRecon_edgewtracking, &b_fDefaultReco_fBasicRecon_edgewtracking);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.hastrk", fDefaultReco_fBasicRecon_hastrk, &b_fDefaultReco_fBasicRecon_hastrk);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.matchtrk", fDefaultReco_fBasicRecon_matchtrk, &b_fDefaultReco_fBasicRecon_matchtrk);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.modfc", fDefaultReco_fBasicRecon_modfc, &b_fDefaultReco_fBasicRecon_modfc);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.penIron", fDefaultReco_fBasicRecon_penIron, &b_fDefaultReco_fBasicRecon_penIron);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.muE", fDefaultReco_fBasicRecon_muE, &b_fDefaultReco_fBasicRecon_muE);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.nuErec", fDefaultReco_fBasicRecon_nuErec, &b_fDefaultReco_fBasicRecon_nuErec);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.nhitx", fDefaultReco_fBasicRecon_nhitx, &b_fDefaultReco_fBasicRecon_nhitx);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.nhity", fDefaultReco_fBasicRecon_nhity, &b_fDefaultReco_fBasicRecon_nhity);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.ntrackhitx", fDefaultReco_fBasicRecon_ntrackhitx, &b_fDefaultReco_fBasicRecon_ntrackhitx);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.ntrackhity", fDefaultReco_fBasicRecon_ntrackhity, &b_fDefaultReco_fBasicRecon_ntrackhity);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.retracktest", fDefaultReco_fBasicRecon_retracktest, &b_fDefaultReco_fBasicRecon_retracktest);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.trg_sec", fDefaultReco_fBasicRecon_trg_sec, &b_fDefaultReco_fBasicRecon_trg_sec);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.upstreamtpl", fDefaultReco_fBasicRecon_upstreamtpl, &b_fDefaultReco_fBasicRecon_upstreamtpl);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexz", fDefaultReco_fBasicRecon_vertexz, &b_fDefaultReco_fBasicRecon_vertexz);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexxz", fDefaultReco_fBasicRecon_vertexxz, &b_fDefaultReco_fBasicRecon_vertexxz);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexyz", fDefaultReco_fBasicRecon_vertexyz, &b_fDefaultReco_fBasicRecon_vertexyz);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.angle", fDefaultReco_fBasicRecon_angle, &b_fDefaultReco_fBasicRecon_angle);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.thetax", fDefaultReco_fBasicRecon_thetax, &b_fDefaultReco_fBasicRecon_thetax);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.thetay", fDefaultReco_fBasicRecon_thetay, &b_fDefaultReco_fBasicRecon_thetay);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexx", fDefaultReco_fBasicRecon_vertexx, &b_fDefaultReco_fBasicRecon_vertexx);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexy", fDefaultReco_fBasicRecon_vertexy, &b_fDefaultReco_fBasicRecon_vertexy);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.horn250", fDefaultReco_fBasicRecon_horn250, &b_fDefaultReco_fBasicRecon_horn250);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vpe", fDefaultReco_fBasicRecon_vpe, &b_fDefaultReco_fBasicRecon_vpe);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.startxpln", fDefaultReco_fBasicRecon_startxpln, &b_fDefaultReco_fBasicRecon_startxpln);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.startypln", fDefaultReco_fBasicRecon_startypln, &b_fDefaultReco_fBasicRecon_startypln);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.startxch", fDefaultReco_fBasicRecon_startxch, &b_fDefaultReco_fBasicRecon_startxch);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.startych", fDefaultReco_fBasicRecon_startych, &b_fDefaultReco_fBasicRecon_startych);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.endxpln", fDefaultReco_fBasicRecon_endxpln, &b_fDefaultReco_fBasicRecon_endxpln);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.endypln", fDefaultReco_fBasicRecon_endypln, &b_fDefaultReco_fBasicRecon_endypln);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.endxch", fDefaultReco_fBasicRecon_endxch, &b_fDefaultReco_fBasicRecon_endxch);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.endych", fDefaultReco_fBasicRecon_endych, &b_fDefaultReco_fBasicRecon_endych);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexx_true", fDefaultReco_fBasicRecon_vertexx_true, &b_fDefaultReco_fBasicRecon_vertexx_true);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexy_true", fDefaultReco_fBasicRecon_vertexy_true, &b_fDefaultReco_fBasicRecon_vertexy_true);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.vertexz_true", fDefaultReco_fBasicRecon_vertexz_true, &b_fDefaultReco_fBasicRecon_vertexz_true);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.lverx", fDefaultReco_fBasicRecon_lverx, &b_fDefaultReco_fBasicRecon_lverx);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.lvery", fDefaultReco_fBasicRecon_lvery, &b_fDefaultReco_fBasicRecon_lvery);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.lverxz", fDefaultReco_fBasicRecon_lverxz, &b_fDefaultReco_fBasicRecon_lverxz);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.lveryz", fDefaultReco_fBasicRecon_lveryz, &b_fDefaultReco_fBasicRecon_lveryz);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.nhits", fDefaultReco_fBasicRecon_nhits, &b_fDefaultReco_fBasicRecon_nhits);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.fIngridHit[1000]", fDefaultReco_fBasicRecon_fIngridHit, &b_fDefaultReco_fBasicRecon_fIngridHit);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.fTrack[10]", fDefaultReco_fBasicRecon_fTrack, &b_fDefaultReco_fBasicRecon_fTrack);
   fChain->SetBranchAddress("fDefaultReco.fBasicRecon.ntracks", fDefaultReco_fBasicRecon_ntracks, &b_fDefaultReco_fBasicRecon_ntracks);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc", &fDefaultReco_f1stReduc_, &b_fDefaultReco_f1stReduc_);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.fUniqueID", &fDefaultReco_f1stReduc_fUniqueID, &b_fDefaultReco_f1stReduc_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.fBits", &fDefaultReco_f1stReduc_fBits, &b_fDefaultReco_f1stReduc_fBits);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.hitmod", &fDefaultReco_f1stReduc_hitmod, &b_fDefaultReco_f1stReduc_hitmod);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.hitcyc", &fDefaultReco_f1stReduc_hitcyc, &b_fDefaultReco_f1stReduc_hitcyc);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.xhitbit", &fDefaultReco_f1stReduc_xhitbit, &b_fDefaultReco_f1stReduc_xhitbit);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.yhitbit", &fDefaultReco_f1stReduc_yhitbit, &b_fDefaultReco_f1stReduc_yhitbit);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.nhitxlyr", &fDefaultReco_f1stReduc_nhitxlyr, &b_fDefaultReco_f1stReduc_nhitxlyr);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.nhitylyr", &fDefaultReco_f1stReduc_nhitylyr, &b_fDefaultReco_f1stReduc_nhitylyr);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.xtotpe", &fDefaultReco_f1stReduc_xtotpe, &b_fDefaultReco_f1stReduc_xtotpe);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.ytotpe", &fDefaultReco_f1stReduc_ytotpe, &b_fDefaultReco_f1stReduc_ytotpe);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.xtracklike", &fDefaultReco_f1stReduc_xtracklike, &b_fDefaultReco_f1stReduc_xtracklike);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.ytracklike", &fDefaultReco_f1stReduc_ytracklike, &b_fDefaultReco_f1stReduc_ytracklike);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.hity", &fDefaultReco_f1stReduc_hity, &b_fDefaultReco_f1stReduc_hity);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.hitx", &fDefaultReco_f1stReduc_hitx, &b_fDefaultReco_f1stReduc_hitx);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.hitxz", &fDefaultReco_f1stReduc_hitxz, &b_fDefaultReco_f1stReduc_hitxz);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.hityz", &fDefaultReco_f1stReduc_hityz, &b_fDefaultReco_f1stReduc_hityz);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.nhits", &fDefaultReco_f1stReduc_nhits, &b_fDefaultReco_f1stReduc_nhits);
   fChain->SetBranchAddress("fDefaultReco.f1stReduc.fIngridHit[1000]", &fDefaultReco_f1stReduc_fIngridHit, &b_fDefaultReco_f1stReduc_fIngridHit);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack", &fDefaultReco_fIngridTrack_, &b_fDefaultReco_fIngridTrack_);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.fUniqueID", fDefaultReco_fIngridTrack_fUniqueID, &b_fDefaultReco_fIngridTrack_fUniqueID);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.fBits", fDefaultReco_fIngridTrack_fBits, &b_fDefaultReco_fIngridTrack_fBits);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.vtxi[3]", fDefaultReco_fIngridTrack_vtxi, &b_fDefaultReco_fIngridTrack_vtxi);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.vtxf[3]", fDefaultReco_fIngridTrack_vtxf, &b_fDefaultReco_fIngridTrack_vtxf);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.length", fDefaultReco_fIngridTrack_length, &b_fDefaultReco_fIngridTrack_length);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.ekin", fDefaultReco_fIngridTrack_ekin, &b_fDefaultReco_fIngridTrack_ekin);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.tx", fDefaultReco_fIngridTrack_tx, &b_fDefaultReco_fIngridTrack_tx);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.ty", fDefaultReco_fIngridTrack_ty, &b_fDefaultReco_fIngridTrack_ty);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.etx", fDefaultReco_fIngridTrack_etx, &b_fDefaultReco_fIngridTrack_etx);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.ety", fDefaultReco_fIngridTrack_ety, &b_fDefaultReco_fIngridTrack_ety);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.ex0", fDefaultReco_fIngridTrack_ex0, &b_fDefaultReco_fIngridTrack_ex0);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.ey0", fDefaultReco_fIngridTrack_ey0, &b_fDefaultReco_fIngridTrack_ey0);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.covx", fDefaultReco_fIngridTrack_covx, &b_fDefaultReco_fIngridTrack_covx);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.covy", fDefaultReco_fIngridTrack_covy, &b_fDefaultReco_fIngridTrack_covy);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.chi2x", fDefaultReco_fIngridTrack_chi2x, &b_fDefaultReco_fIngridTrack_chi2x);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.chi2y", fDefaultReco_fIngridTrack_chi2y, &b_fDefaultReco_fIngridTrack_chi2y);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.btheta", fDefaultReco_fIngridTrack_btheta, &b_fDefaultReco_fIngridTrack_btheta);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.bphi", fDefaultReco_fIngridTrack_bphi, &b_fDefaultReco_fIngridTrack_bphi);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.mrdhitid[2]", fDefaultReco_fIngridTrack_mrdhitid, &b_fDefaultReco_fIngridTrack_mrdhitid);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.mucl", fDefaultReco_fIngridTrack_mucl, &b_fDefaultReco_fIngridTrack_mucl);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.vpe", fDefaultReco_fIngridTrack_vpe, &b_fDefaultReco_fIngridTrack_vpe);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.view", fDefaultReco_fIngridTrack_view, &b_fDefaultReco_fIngridTrack_view);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.nhits", fDefaultReco_fIngridTrack_nhits, &b_fDefaultReco_fIngridTrack_nhits);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.nsimparticles", fDefaultReco_fIngridTrack_nsimparticles, &b_fDefaultReco_fIngridTrack_nsimparticles);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.nsimemshowers", fDefaultReco_fIngridTrack_nsimemshowers, &b_fDefaultReco_fIngridTrack_nsimemshowers);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.fIngridHit[1000]", fDefaultReco_fIngridTrack_fIngridHit, &b_fDefaultReco_fIngridTrack_fIngridHit);
   fChain->SetBranchAddress("fDefaultReco.fIngridTrack.fSimParticle[10]", fDefaultReco_fIngridTrack_fSimParticle, &b_fDefaultReco_fIngridTrack_fSimParticle);
   Notify();
}

Bool_t tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_cxx

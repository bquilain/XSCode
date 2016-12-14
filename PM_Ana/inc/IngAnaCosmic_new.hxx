#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>

#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>

#include "TApplication.h"
//#include "setup.hxx"

//#include "root_setup.hxx"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "IngridTrackSummary.h"

#include "ReconTrackBasic.cxx"
#include "IngEvtDisp.cxx"
#include "IngridConstants.h"
#include "INGRID_Dimension.hxx"


#include <vector>
//#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
//#endif


const static int idtest[26][2] =
  { // {plane# ,  view#}
    { 0, 0}, { 0, 1}, 
    { 1, 0}, { 1, 1}, 
    { 2, 0}, { 2, 1}, 
    { 3, 0}, { 3, 1}, 
    { 4, 0}, { 4, 1}, 
    { 5, 0}, { 5, 1}, 
    { 6, 0}, { 6, 1}, 
    { 7, 0}, { 7, 1}, 
    { 8, 0}, { 8, 1}, 
    { 9, 0}, { 9, 1}, 
    {10, 0}, {10, 1}, 
    {11, 1},
    {12, 1},
    {13, 0},
    {14, 0}
  };
const static int  ntest = 26;

const static float CUT_PATH_LENZ = PlnIronThick * 5;


FileStat_t           fs;
TFile*            rfile;
TTree*            rtree;
TBranch*          EvtBr;
IngridEventSummary* evt;

int ReadFile(string str){


  //######## read root file #############################
  //#####################################################
  cout<<"reading "<<str<<"....."<<endl;
  if(gSystem->GetPathInfo(str.c_str(),fs)){
    cout<<"Cannot open file "<<str<<endl;
    return -1;
  }
  evt    = new IngridEventSummary();
  rfile  = new TFile(str.c_str(),"read");
  rtree  = (TTree*)rfile -> Get("tree");
  EvtBr  = rtree->GetBranch("fDefaultReco.");
  EvtBr  ->SetAddress(&evt);
  rtree  ->SetBranchAddress("fDefaultReco.", &evt);
  return (int)rtree -> GetEntries();


}

//### OUTPUT ROOT File ####
int             run;
int          numevt;
float         angle;
float          path;
int             mod;
int            view;
int             pln;

float  chi2x, chi2y;
float          chi2;
float            ax;
float            bx;
float            ay;
float            by;
int           expch;
float         expxy;
float          expz;
float         expzh;
float         expzl;
vector<int>      ch;
vector<float>    pe;
vector<float>    xy;
vector<float>  diff;

//### flag for selection
bool       trkmatching;
bool       trkendmatching;
bool       trklen;
bool       goodchi2;
bool       goVETO;
bool       notescape;
bool       trgsamepln;
bool       notstop;
 bool       anaflag;



class ch_info{
 public:
  int      ch;
  float    pe;
  float    xy;
  float  diff;
};
bool operator<(const ch_info& left, const ch_info& right){
  return left.diff < right.diff;
}
bool operator>(const ch_info& left, const ch_info& right){
  return left.diff > right.diff;
}
vector<ch_info> fch_info;



TFile*     fTFile;
TTree*       tree;
void Book(string str){
  gROOT->ProcessLine("#include<vector>");
  fTFile = new TFile(str.c_str(), "recreate");
  tree   = new TTree("tree", "tree");

  tree   ->Branch("run",          &run,          "run/I");
  tree   ->Branch("numevt",       &numevt,       "numevt/I");
  tree   ->Branch("angle",        &angle,        "angle/F");
  tree   ->Branch("path",         &path,         "path/F");
  tree   ->Branch("mod",          &mod,          "mod/I");
  tree   ->Branch("view",         &view,         "view/I");
  tree   ->Branch("pln",          &pln,          "pln/I");

  tree   ->Branch("chi2x",        &chi2x,        "chi2x/F");
  tree   ->Branch("chi2y",        &chi2y,        "chi2y/F");
  tree   ->Branch("ax",           &ax,           "ax/F");
  tree   ->Branch("bx",           &bx,           "bx/F");
  tree   ->Branch("ay",           &ay,           "ay/F");
  tree   ->Branch("by",           &by,           "by/F");

  tree   ->Branch("expch",        &expch,        "expch/I");
  tree   ->Branch("expxy",        &expxy,        "expxy/F");
  tree   ->Branch("expz",         &expz,         "expz/F");
  tree   ->Branch("expzh",        &expzh,        "expzh/F");
  tree   ->Branch("expzl",        &expzl,        "expzl/F");

  tree   ->Branch("trkmatching",   &trkmatching,   "trkmatching/O");
  tree   ->Branch("trkendmatching",&trkendmatching,"trkendmatching/O");
  tree   ->Branch("trklen",        &trklen,        "trklen/O");
  tree   ->Branch("goodchi2",      &goodchi2,      "goodchi2/O");
  tree   ->Branch("goVETO",        &goVETO,        "goVETO/O");
  tree   ->Branch("notescape",     &notescape,     "notescape/O");
  tree   ->Branch("trgsamepln",    &trgsamepln,    "trgsamepln/O");
  tree   ->Branch("notstop",       &notstop,       "notstop/O");
  tree   ->Branch("anaflag",       &anaflag,       "anaflag/O");


  tree   ->Branch("ch",           &ch);
  tree   ->Branch("pe",           &pe);
  tree   ->Branch("xy",           &xy);
  tree   ->Branch("diff",         &diff);

}




void Fill(){

  tree   ->Fill();
}

void Clear(){
  ch.      clear();
  pe.      clear();
  xy.      clear();
  diff.    clear();


}

void Write(){
  tree   -> Write();
  fTFile -> Write();
  fTFile -> Close();
}

//#### For Clear Cosmic Slection ####
class hit{
public:
  int  view;
  int   pln;
  int    ch;
  float  pe;
};
vector<hit> vhit;
bool fClearCosmic_1( IngridBasicReconSummary& basicrecon ){

  vhit.clear();
  for( int ihit=0; ihit < basicrecon. nhits; ihit++ ){
    IngridHitSummary* b = basicrecon. GetIngridHit(ihit);
    hit t;
    t.view = b-> view;
    t.pln  = b->  pln;
    t.ch   = b->   ch;
    t.pe   = b->   pe;
    vhit.push_back(t);
  }

  for(int i=0; i<vhit.size(); i++){
    for(int j=i+1; j<vhit.size(); j++){
      if( vhit[i].pe    > 10.5         &&
	  vhit[j].pe    > 10.5         &&
	  vhit[i].pln  == vhit[j].pln  &&
	  vhit[i].view == vhit[j].view &&
	  abs(vhit[i].ch - vhit[j].ch ) > 1 &&
	  vhit[i].pln < nTPL
	  )
	return false;
    }
  }
  return true;
}
ReconTrackBasic* frecontrack;
int nhits;

IngridBasicReconSummary* testbas    ;
IngridBasicReconSummary* trackingbas;
TRef  fIngridHit[INGRIDHIT_MAXHITS];
void AllClear(){
  nhits       = 0;
  testbas     -> Clear();
  trackingbas -> Clear();
  frecontrack -> Clear();
}

void CopyHit(IngridBasicReconSummary* basicrecon, 
	     IngridBasicReconSummary* testbas,
	     IngridBasicReconSummary* trackingbas,
	     int ipln, int iview){
  nhits=0;
  for(int ihit=0; ihit < basicrecon->nhits; ihit++){
    IngridHitSummary* b = basicrecon -> GetIngridHit(ihit);
    if(  ipln  == b->pln &&
	 ( iview == b->view || iview == -1 ) ){
      testbas -> AddIngridHit(b);
    
    }
    else{
      fIngridHit[nhits] = 0;
      fIngridHit[nhits] = TRef( (IngridHitSummary*)basicrecon -> GetIngridHit(ihit) );
      nhits++;
      trackingbas -> AddIngridHit( b );


    }
  }
  for(int ini=nhits; ini<INGRIDHIT_MAXHITS; ini++){
    fIngridHit[ini]=0;
  }

}

void SetHit2Basic(IngridEventSummary*      evt, 
		  IngridBasicReconSummary* basicrecon, 
		  int                      mod,
		  int                      cyc){
	    
  int nhit = evt -> NIngridModHits(mod, cyc);
  IngridHitSummary* b;
  for(int ihit=0; ihit < nhit; ihit++){
    b = (IngridHitSummary*) (evt->GetIngridModHit(ihit, mod, cyc));
    cout << "plane  : " << b->pln
	 << "\tview : " << b->view
	 << "\tch   : " << b->ch
	 << "\tpe   : " << b->pe
	 << endl;
    basicrecon -> AddIngridHit(b);
  }

}


//###############################################
//#### for VETO analysis ########################
//###############################################

bool SpecialForBVETO( float a, float b, int mod, int pln, int view, int cyc,IngridEventSummary* evt ){
  //cout << "-------------------------------------" << endl;
  //### calculate expected position at bottom module ###
  float expected_z_at_bMod = 1.0 * ( -35 - b )/a;


  //cout << "expected Z at bottom Module :" << expected_z_at_bMod << endl;
  if( StartModuleZ       < expected_z_at_bMod &&
      expected_z_at_bMod < EndModuleZ ){
    //cout << "bottom module should have hit" << endl;
  }
  else{
    //return true;
    return false;
  }

  int nhit = evt -> NIngridModHits(mod-1, cyc);

  bool hashit=false;
  for(int ihit=0; ihit < nhit; ihit++){
    //cout << nhit << "\t" << ihit << endl;
    IngridHitSummary* b = (IngridHitSummary*) (evt->GetIngridModHit(ihit, mod-1, cyc));
    if( b->pln != UVETO && 
	b->pe  > 6.5    &&
	fabs( b->z - expected_z_at_bMod ) < 15 &&
	b->view == 0
	){
      hashit = true;
    }
  }
  if(hashit){
    //cout << "bottom module has hit" << endl;
    return true;
  }
  return false;
  //cout << "-------------------------------------" << endl;
}


bool GetVetoPosZ(float a, float b, int mod, int ipln, float& exppos){
  if(ipln<11||ipln>14)return false;
  exppos = 1.0 * ( fVetoXY(mod, ipln) - b)/a;
  return true;
}
bool GetPosZ(float a, float b, int mod, int ipln, float& exppos){
  exppos = 1.0 * ( fVetoXY(mod, ipln) - b)/a;
  return true;
}
bool GetVetoPosXY(float a, float b, float posz, float& exppos){
  
  exppos = 1.0 * a * posz + b;
  return true;
}
bool GetPosXY(float a, float b, float posz, float& exppos){
  
  exppos = 1.0 * a * posz + b;
  return true;
}


IngridTrackSummary* retrk0;
IngridTrackSummary* retrk1;
float startxz   ;
float endxz     ;
float startyz   ;
float endyz     ;
int   startypln ;
int   startych  ;
int   endypln   ;
int   endych    ;
int   startxpln ;
int   startxch  ;
int   endxpln   ;
int   endxch    ;



int  diffverz;
bool  startisveto;
bool  endisveto;
float dax, dbx, day, dby;
void GetTrackInfo(){
  retrk0      = frecontrack -> GetTrack(0, 0);
  retrk1      = frecontrack -> GetTrack(1, 0);
  startxz     = retrk0      -> vtxi[2];
  endxz       = retrk0      -> vtxf[2];
  startyz     = retrk1      -> vtxi[2];
  endyz       = retrk1      -> vtxf[2];
  startypln   = frecontrack -> StartPln (FromX);
  startych    = frecontrack -> StartCh  (FromX);
  endypln     = frecontrack -> EndPln   (FromX);
  endych      = frecontrack -> EndCh    (FromX);
  startxpln   = frecontrack -> StartPln (FromY);
  startxch    = frecontrack -> StartCh  (FromY);
  endxpln     = frecontrack -> EndPln   (FromY);
  endxch      = frecontrack -> EndCh    (FromY);
  ax          = retrk1      -> tx; 
  bx          = retrk1      -> etx;
  ay          = retrk0      -> tx; 
  by          = retrk0      -> etx;
  diffverz    = frecontrack -> GetInitialZ( 0 ) - frecontrack -> GetInitialZ( 1 );
  startisveto = frecontrack -> startisveto[0] || frecontrack -> startisveto[1];
  endisveto   = frecontrack -> endisveto[0]   || frecontrack -> endisveto[1];
  day         = retrk0 -> covx;
  dby         = retrk0 -> covy;
  dax         = retrk1 -> covx;
  dbx         = retrk1 -> covy;
  chi2x       = retrk1 -> chi2x;
  chi2y       = retrk0 -> chi2x;
}

float expposz_thisview, expposxy_anotherview;
float expposz_thisview_low,expposz_thisview_high;
const static float TolOfOtherView = 2.5;
const static float TolOfThisView  = 2.5;
const static float TolOfChi2      = 5;
bool isVETO(int ipln){
  if( ipln >= nTPL )
    return true;
  else
    return false;
}

bool fTrgSamePln(int mod, int ipln, int iview){
  if( isVETO(ipln) )
    return true;
  float t;
  if(iview==FromX) GetPosXY( ax, bx, PlnZ[ipln], t);
  if(iview==FromY) GetPosXY( ay, by, PlnZ[ipln], t);

  if     ( ipln == 0 && iview == FromX && startxpln == 0 && startypln == 1)
    return true;
  else if( ipln == 0 && iview == FromY && startypln == 0 && startxpln == 1)
    return true;
  else if( ipln ==10 && iview == FromX && endxpln   ==10 && endypln   == 9)
    return true;
  else if( ipln ==10 && iview == FromY && endypln   ==10 && endxpln   == 9)
    return true;
  else
    return false;
}

bool fNotEscape(int mod, int ipln, int iview){
  if(  !isVETO(ipln) )
    return true;
  if     (  ipln == RVETO )
    GetVetoPosXY( ay, by, expposz_thisview ,  expposxy_anotherview );
  else if(  ipln == LVETO )
    GetVetoPosXY( ay, by, expposz_thisview , expposxy_anotherview );
  else if(  ipln == BVETO )
    GetVetoPosXY( ax, bx, expposz_thisview , expposxy_anotherview );
  else if(  ipln == UVETO )
    GetVetoPosXY( ax, bx, expposz_thisview , expposxy_anotherview );

  if( StartTPLXY + TolOfOtherView    < expposxy_anotherview         &&
      expposxy_anotherview           < EndTPLXY - TolOfOtherView   )
    return true;
  else 
    return false;
}

bool fNotStop(int mod, int ipln, int iview, int print=0){
  if(  !isVETO(ipln) )
    return true;

  float extxy;
  int   extpln;
  if     (  ipln == RVETO ){
    if     ( ax >  0 )
      extpln = startxpln - 1;
    else if( ax <  0 )
      extpln = endxpln   + 1;
    else if( ax == 0 )
      return false;
  }
  else if(  ipln == BVETO ){
    if     ( ay >  0 )
      extpln = startypln - 1;
    else if( ay <  0 )
      extpln = endypln   + 1;
    else if( ay == 0 )
      return false;
  }

  if( ipln==RVETO || ipln==LVETO )GetPosXY( ax, bx, PlnZ[extpln], extxy );
  if( ipln==BVETO || ipln==UVETO )GetPosXY( ay, by, PlnZ[extpln], extxy );
  if(print){
    cout << "extrapolate plane :" << extpln 
	 << " extrapolate xy   :" << extxy
	 << endl; 
  }
  if( StartTPLXY - 0.5 * ScintiWidth < extxy         &&
      extxy                          < EndTPLXY + 0.5 * ScintiWidth   )
    return false;
  else 
    return true;

  /*
  else if(  ipln == LVETO )
    GetVetoPosXY( ay, by, expposz_thisview , expposxy_anotherview );
  else if(  ipln == BVETO )
    GetVetoPosXY( ax, bx, expposz_thisview , expposxy_anotherview );
  else if(  ipln == UVETO )
    GetVetoPosXY( ax, bx, expposz_thisview , expposxy_anotherview );
  */

  if( StartTPLXY + TolOfOtherView    < expposxy_anotherview         &&
      expposxy_anotherview           < EndTPLXY - TolOfOtherView   )
    return true;
  else 
    return false;
}

bool fTrkLen(int print=0){
  if(print)
    cout << "endxz  :" << endxz   <<endl
	 << "startxz:" << startxz <<endl
	 << "diff   :" << endxz - startxz   <<endl
	 << "endyz  :" << endyz   <<endl
	 << "startyz:" << startyz <<endl
	 << "diff   :" << endyz - startyz   <<endl
	 << "thre.  :" << CUT_PATH_LENZ   <<endl
      ;

  if( ( endxz - startxz ) >= CUT_PATH_LENZ-0.1 &&
      ( endyz - startyz ) >= CUT_PATH_LENZ-0.1      //### Track length selection
      )
    return true;
  else
    return false;
}

bool fGoVETO(int mod, int ipln, int iview, int print=0){
  if( !isVETO(ipln) )
    return true;
  
  if     (  ipln == RVETO ) 
    GetVetoPosZ ( ax, bx, mod, ipln , expposz_thisview );
  else if(  ipln == LVETO ) 
    GetVetoPosZ ( ax, bx, mod, ipln , expposz_thisview );
  else if(  ipln == BVETO ) 
    GetVetoPosZ ( ay, by, mod, ipln , expposz_thisview );
  else if(  ipln == UVETO )
    GetVetoPosZ ( ay, by, mod, ipln , expposz_thisview );
    
  if( print )
    cout << "extrapolate z:" << expposz_thisview << endl;

  if(  StartVetoZ - 0.5 * ScintiWidth < expposz_thisview            &&
       expposz_thisview               < EndVetoZ + 0.5 * ScintiWidth  )
    return true;
  else
    return false;
}

bool fGoodChi2(){
  if( chi2x < TolOfChi2 && 
      chi2y < TolOfChi2 ) 
    return true;
  else
    return false;
}

bool fTrkEndMatching(int mod, int ipln, int iview){
  if(ipln!=10)
    return true;
  else if(ipln==10 &&
	  abs(endypln-endxpln)<2)
    return true;
  else
    return false;
}

bool FlagForLastTPL(int mod, int ipln, int iview){
  float t;
  if(iview==FromX) GetPosXY( ax, bx, PlnZ[ipln], t);
  if(iview==FromY) GetPosXY( ay, by, PlnZ[ipln], t);
  float t1;
  if(iview==FromX) GetPosXY( ay-day, by-dby, PlnZ[ipln], t1);
  if(iview==FromY) GetPosXY( ax-dax, bx-dbx, PlnZ[ipln], t1);
  float t2;
  if(iview==FromX) GetPosXY( ay+day, by+dby, PlnZ[ipln], t2);
  if(iview==FromY) GetPosXY( ax+dax, bx+dbx, PlnZ[ipln], t2);

  float tmax = max(t1, t2);
  float tmin = min(t1, t2);

  if     (  ipln      == 10   && iview == FromX && 
            endxpln   == 10   && 
            StartTPLXY <  t   && t < EndTPLXY   &&
            StartTPLXY < tmin && tmax < EndTPLXY )
    return true;
  else if(  ipln      == 10   && iview == FromY && 
	    endypln   == 10   &&
	    StartTPLXY <  t   && t < EndTPLXY   &&
	    StartTPLXY < tmin && tmax < EndTPLXY )
    return true;
  return false;
}
bool FlagForBVETO(int mod, int ipln, int iview){
  float expected;
  if(ay < 0){
    //GetPosXY( ay + day, by + dby, PlnZ[endypln] + PlnIronThick, expected );
    GetPosXY( ay , by , PlnZ[endypln] + PlnIronThick, expected );
  }
  if(ay > 0){
    //GetPosXY( ay - day, by + dby, PlnZ[startypln] - PlnIronThick, expected );
    GetPosXY( ay , by , PlnZ[startypln] - PlnIronThick, expected );
  }

  if( expected > StartTPLXY - 0.5 * ScintiWidth  )
    return false;

  return true;
}

bool FlagForRVETO(int mod, int ipln, int iview){
  float expected;
  if(ax < 0){
    GetPosXY( ax + dax, bx + dbx, PlnZ[endxpln] + PlnIronThick, expected );
  }
  if(ax > 0){
    GetPosXY( ax - dax, bx + dbx, PlnZ[startxpln] - PlnIronThick, expected );
  }
  if( expected > StartTPLXY - 0.5 * ScintiWidth  )
    return false;

  return true;
}

bool FlagForLVETO(int mod, int ipln, int iview){
  float expected;
  if(ax > 0){
    //GetPosXY( ax - dax, bx - dbx, PlnZ[endxpln] + PlnIronThick, expected );
    GetPosXY( ax , bx , PlnZ[endxpln] + PlnIronThick, expected );
  }
  if(ax < 0){
    //GetPosXY( ax + dax, bx - dbx, PlnZ[startxpln] - PlnIronThick, expected );
    GetPosXY( ax , bx , PlnZ[startxpln] - PlnIronThick, expected );
  }


  if( expected < EndTPLXY + 0.5 * ScintiWidth )
    return false;

  return true;
}


IngridBasicReconSummary*  basicrecon;

void Get(int print = 0){
  angle      = TMath::ATan( sqrt(ax*ax+ay*ay) ) / TMath::Pi() * 180;  	    

  path       = 1.0 / cos(angle/180*3.14);
  if( pln >= nTPL )
    path       = 1.0 / sin(angle/180*3.14);

  if( view == 0 )
    GetPosXY( ay, by, PlnZ[pln], expxy );
  else if( view == 1 )
    GetPosXY( ax, bx, PlnZ[pln], expxy );

  if( pln >= nTPL && view == 0 ) 
    GetPosZ ( ay, by, mod, pln, expz );
  else if( pln >= nTPL && view == 1 ) 
    GetPosZ ( ax, bx, mod, pln, expz );

  //####
  if( pln == RVETO ){
    if(ax<0){
      expzh = 1.0 * ( fVetoXY(mod, pln) - (bx+dbx))/(ax+dax);
      expzl = 1.0 * ( fVetoXY(mod, pln) - (bx-dbx))/(ax-dax);
    }
    if(ax>0){
      expzh = 1.0 * ( fVetoXY(mod, pln) - (bx-dbx))/(ax-dax);
      expzl = 1.0 * ( fVetoXY(mod, pln) - (bx+dbx))/(ax+dax);
    }
  }

  if( pln == BVETO ){
    if(ay<0){
      expzh = 1.0 * ( fVetoXY(mod, pln) - (by+dby))/(ay+day);
      expzl = 1.0 * ( fVetoXY(mod, pln) - (by-dby))/(ay-day);
    }
    if(ay>0){
      expzh = 1.0 * ( fVetoXY(mod, pln) - (by-dby))/(ay-day);
      expzl = 1.0 * ( fVetoXY(mod, pln) - (by+dby))/(ay+day);
    }
  }


  if( !isVETO(pln) ){
    expch   =  (int)( (expxy+0.5*ScintiWidth) / ScintiWidth) ;
    if( expch==-1 )expch=0;
    if( expch==24 )expch=23;
    if(print)cout << "expxy :" << expxy << "\texpch :" << expch << endl;
  }
  else {
    expch   =  (int)( (expz+1) / ScintiWidth) ;
    if(print)cout << "expz :" << expz << "\texpch :" << expch << endl;
  }
  //### channel Info ####
  fch_info.clear();


  for(int itestbrec = 0; itestbrec < testbas -> nhits; itestbrec++ ){
    IngridHitSummary* t = (IngridHitSummary*)testbas->GetIngridHit(itestbrec);
    if( t -> view == view ){


      ch_info tc;
      tc.ch   =  t->ch ;
      tc.pe   =  t->pe ;
      tc.xy   =  t->xy ;
      if      ( pln <  nTPL )
	tc.diff =  fabs( t->xy - expxy ) ;
      else if( pln >= nTPL ) 
	tc.diff =  fabs( t->z - expz ) ;
      fch_info.push_back(tc);
    }
  }
  sort( fch_info.begin(), fch_info.end() );
  if( fch_info.size()==0 ){
    ch_info tc;
    tc.ch   =  expch;
    tc.pe   =  -1e-5 ;
    tc.xy   =  130;
    tc.diff =  130;
    fch_info.push_back(tc);
  }

  for(int i=0; i<fch_info.size(); i++){
    ch.   push_back( fch_info[i].ch ); 
    pe.   push_back( fch_info[i].pe ); 
    xy.   push_back( fch_info[i].xy ); 
    diff. push_back( fch_info[i].diff ); 
  }

  ay         = retrk0 -> tx;  
  by         = retrk0 -> etx;
  ax         = retrk1 -> tx;
  bx         = retrk1 -> etx;



}


void Print(){
  cout << "---------------------- "  << endl;
  cout << "analysis module# :" <<  mod  << endl;
  if(view == 0)
    cout << "analysis          Side View" << endl; 
  if(view == 1)
    cout << "analysis          Top View" << endl; 
  cout << "analysis plane#  :" <<  pln              << endl
       << "difference       :" <<  fch_info[0].diff << endl
       << "trkmatching       " << trkmatching    << endl
       << "trkendmatching    " << trkendmatching << endl
       << "trklength         " << trklen         << endl
       << "good chi2         " << goodchi2       << endl
       << "goVETO            " << goVETO         << endl
       << "not escape        " << notescape      << endl
       << "trigger same pln  " << trgsamepln     << endl
       << "notstop           " << notstop        << endl
    ;

  cout << "ax                " << ax << "\t+-" << dax << endl
       << "dbx               " << bx << "\t+-" << dbx
       << endl;
  cout << "expected          " << expz  << endl;
  cout << "expected +        " << expzh << endl;
  cout << "expected -        " << expzl << endl;


  cout << "expch    :"         << expch << endl
    ;

}

void PrintHit(IngridBasicReconSummary* basicrecon){
  for(int ihit=0; ihit < basicrecon->nhits; ihit++){
    IngridHitSummary* b = basicrecon -> GetIngridHit(ihit);

    cout << "pln :" << b->pln  << "\t"
	 << "view:" << b->view << "\t"
	 << "ch  :" << b->ch   << "\t"
      	 << "pe  :" << b->pe   << "\t"
	 << "time:" << b->time << "\t"
	 << endl;
  }

}

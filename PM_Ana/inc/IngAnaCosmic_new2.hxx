//### C++ 
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>
//### ROOT
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
//### INGRID data class
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "IngridTrackSummary.h"
//### INGRID softwares
#include "ReconTrackBasic.cxx"
#include "IngEvtDisp.cxx"
#include "IngridConstants.h"
#include "INGRID_Dimension.hxx"
//### in order to use vector
#include <vector>
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;




FileStat_t           fs;
TFile*            rfile;
TTree*            rtree;
TBranch*          EvtBr;
IngridEventSummary* evt;

bool isVETO(int ipln){
  if( ipln >= nTPL )
    return true;
  else
    return false;
}

int ReadFile(string str){
  //######## read root file #############################
  //#####################################################
  cout<<"reading "<<str<<"....."<<endl;
  if(gSystem->GetPathInfo(str.c_str(),fs)){
    cout<<"Cannot open file "<<str<<endl;
    return -1;
  }
  evt    = new IngridEventSummary    ();
  rfile  = new TFile                 (str.c_str(),"read");
  rtree  = (TTree*)rfile -> Get      ("tree");
  EvtBr  = rtree         -> GetBranch("fDefaultReco.");
  EvtBr  ->SetAddress                (&evt);
  rtree  ->SetBranchAddress          ("fDefaultReco.", &evt);
  return (int)rtree -> GetEntries();
}

//### OUTPUT ROOT File ####
int              run;
int              ievt;
int              mod;
int              cyc;
int              view;
int              pln;

float            chi2x;
float            chi2y;
float            ax;
float            bx;
float            ay;
float            by;
float            angle;
float            path;
int              expch;

int              diffstartz;   //## difference start z[plane]
int              diffendz;     //## difference end   z[plane]
float            svpreVETOxy;  //## xy at plane extrapolated to VETO
float            svexpz;       //## expected z[cm]  at Same  View
float            ovexpxy;      //## expected xy[cm] at Other View
int              trklenx;      //## number of scintilator plane of track x
int              trkleny;      //## number of scintilator plane of track x

bool             fsvexpz;
bool             fovexpxy;
bool             fsvprevetoxy;
bool             ftrklenx;
bool             ftrkleny;
bool             ftrkstart;
bool             ftrkend;
bool             fana;

//### hit channel information 
vector<int>      ch;
vector<float>    pe;
vector<float>    xy;
vector<float>    diff;


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
  tree   ->Branch("ievt",         &ievt,         "ievt/I");
  tree   ->Branch("mod",          &mod,          "mod/I");
  tree   ->Branch("view",         &view,         "view/I");
  tree   ->Branch("pln",          &pln,          "pln/I");

  tree   ->Branch("chi2x",        &chi2x,        "chi2x/F");
  tree   ->Branch("chi2y",        &chi2y,        "chi2y/F");
  tree   ->Branch("ax",           &ax,           "ax/F");
  tree   ->Branch("bx",           &bx,           "bx/F");
  tree   ->Branch("ay",           &ay,           "ay/F");
  tree   ->Branch("by",           &by,           "by/F");
  tree   ->Branch("angle",        &angle,        "angle/F");
  tree   ->Branch("path",         &path,         "path/F");
  tree   ->Branch("expch",        &expch,        "expch/I");

  tree   ->Branch("svpreVETOxy",  &svpreVETOxy,  "svpreVETOxy/F");
  tree   ->Branch("diffstartz",   &diffstartz,   "diffstartz/I");
  tree   ->Branch("diffendz",     &diffendz,     "diffendz/I");
  tree   ->Branch("svexpz",       &svexpz,       "svexpz/F");
  tree   ->Branch("ovexpxy",      &ovexpxy,      "ovexpxy/F");
  tree   ->Branch("trklenx",      &trklenx,      "trklenx/I");
  tree   ->Branch("trkleny",      &trkleny,      "trkleny/I");

  tree   ->Branch("fsvexpz",      &fsvexpz,      "fsvexpz/O");
  tree   ->Branch("fovexpxy",     &fovexpxy,     "fovexpxy/O");
  tree   ->Branch("fsvprevetoxy", &fsvprevetoxy, "fsvprevetoxy/O");
  tree   ->Branch("ftrklenx",     &ftrklenx,     "ftrklenx/O");
  tree   ->Branch("ftrkleny",     &ftrkleny,     "ftrkleny/O");
  tree   ->Branch("ftrkstart",    &ftrkstart,    "ftrkstart/O");
  tree   ->Branch("ftrkend",      &ftrkend,      "ftrkend/O");
  tree   ->Branch("fana",         &fana,         "fana/O");


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


ReconTrackBasic* frecontrack;
int nhits;

IngridBasicReconSummary* basicrecon ;
IngridBasicReconSummary* testbas    ;
IngridBasicReconSummary* trackingbas;
TRef  fIngridHit[INGRIDHIT_MAXHITS];


void SetHit2Basic(IngridBasicReconSummary* bas, 
		  IngridBasicReconSummary* basicrecon, 
		  int                      imod,
		  int                      icyc,
		  int                     iview,
		  int                      ipln){
	    
  int nhit = bas -> Nhits();
  IngridHitSummary* b;
  for(int ihit=0; ihit < nhit; ihit++){
    b = (IngridHitSummary*) (bas->GetIngridHit(ihit));
    if( b->view == iview && b->pln == ipln )
      basicrecon -> AddIngridHit(b);
  }

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
int   diffverz    ;
float dax, dbx, day, dby;
void GetTrackInfo(){
  retrk0      = frecontrack -> GetTrack(0, 0);
  retrk1      = frecontrack -> GetTrack(1, 0);
  startxz     = retrk0      -> vtxi[2];
  endxz       = retrk0      -> vtxf[2];
  startyz     = retrk1      -> vtxi[2];
  endyz       = retrk1      -> vtxf[2];
  startypln   = frecontrack -> StartTPL (FromX);
  startych    = frecontrack -> StartCh  (FromX);
  endypln     = frecontrack -> EndTPL   (FromX);
  endych      = frecontrack -> EndCh    (FromX);
  startxpln   = frecontrack -> StartTPL (FromY);
  startxch    = frecontrack -> StartCh  (FromY);
  endxpln     = frecontrack -> EndTPL   (FromY);
  endxch      = frecontrack -> EndCh    (FromY);
  ax          = retrk1      -> tx; 
  bx          = retrk1      -> etx;
  ay          = retrk0      -> tx; 
  by          = retrk0      -> etx;
  diffverz    = frecontrack -> GetInitialZ( 0 ) - frecontrack -> GetInitialZ( 1 );
  day         = retrk0 -> covx;
  dby         = retrk0 -> covy;
  dax         = retrk1 -> covx;
  dbx         = retrk1 -> covy;
  chi2x       = retrk1 -> chi2x;
  chi2y       = retrk0 -> chi2x;

  diffstartz  = startxpln - startypln;
  diffendz    = endxpln   - endypln;
  trklenx     = endxpln   - startxpln;
  trkleny     = endypln   - startypln;
  angle       = TMath::ATan( sqrt(ax*ax+ay*ay) ) / TMath::Pi() * 180;  	    
  path        = 1.0 / cos(angle/180*3.14);
  if( isVETO(pln) )
    path       = 1.0 / sin(angle/180*3.14);



}

void GetHitChInfo(){
  expch   =  (int)( (svexpz+1) / ScintiWidth) ;
  fch_info.clear();
  for(int itestbrec = 0; itestbrec < testbas -> nhits; itestbrec++ ){
    IngridHitSummary* t = (IngridHitSummary*)testbas->GetIngridHit(itestbrec);
    ch_info tc;
    tc.ch   =  t->ch ;
    tc.pe   =  t->pe ;
    tc.xy   =  t->xy ;
    tc.diff =  fabs( t->z - svexpz ) ;
    fch_info.push_back(tc);
    
  }
  sort( fch_info.begin(), fch_info.end() );
  if( fch_info.size()==0 ){
    ch_info tc;
    tc.ch   =  -1;
    tc.pe   =  -1e-5 ;
    tc.xy   =  130;
    tc.diff =  130;
    fch_info.push_back(tc);
  }
  ch.clear();
  pe.clear();
  xy.clear();
  diff.clear();
  for(int i=0; i<fch_info.size(); i++){
    ch.   push_back( fch_info[i].ch ); 
    pe.   push_back( fch_info[i].pe ); 
    xy.   push_back( fch_info[i].xy ); 
    diff. push_back( fch_info[i].diff ); 
  }

}

float fSVPreVETOXY(){
  if( view == FromX && pln == BVETO ){ //### Bottom VETO
    if( ay > 0 )
      return ( ay * PlnZ[ startypln - 1 ] + by );
    if( ay < 0 )
      return ( ay * PlnZ[ endypln   + 1 ] + by );
  }

  if( view == FromX && pln == UVETO ){ //### Up VETO
    if( ay > 0 )
      return ( ay * PlnZ[ endypln     + 1 ] + by );
    if( ay < 0 )
      return ( ay * PlnZ[ startypln   - 1 ] + by );
  }

  if( view == FromY && pln == RVETO ){ //### Right VETO
    if( ax > 0 )
      return ( ax * PlnZ[ startxpln   - 1 ] + bx );
    if( ax < 0 )
      return ( ax * PlnZ[ endxpln     + 1 ] + bx );
  }

  if( view == FromY && pln == LVETO ){ //### Left VETO
    if( ax > 0 ){
      return ( ax * PlnZ[ endxpln     + 1 ] + bx );
    }
    if( ax < 0 ){
      return ( ax * PlnZ[ startxpln   - 1 ] + bx );
    }
  }
}

float fSVExpZ(){

  if( view == FromX ){
    return ( 1.0*( fVetoXY(mod, pln) - by )/ay );
  }
  if( view == FromY ){
    return ( 1.0*( fVetoXY(mod, pln) - bx )/ax );
  }

}

float fOVExpXY(){
  if( view == FromX ){
    return ( ax * svexpz + bx );
  }
  if( view == FromY ){
    return ( ay * svexpz + by );
  }
}

void Print(){
  cout << "evt#          :" << ievt         << endl;
  cout << "module        :" << mod          << endl;
  cout << "diffstartz    :" << diffstartz   << endl;
  cout << "diffendz      :" << diffendz     << endl;
  cout << "svpreVETOxy   :" << svpreVETOxy  << endl;
  cout << "svexpz        :" << svexpz       << endl;
  cout << "ovexpxy       :" << ovexpxy      << endl;
  cout << "trklenx       :" << trklenx      << endl;
  cout << "trkleny       :" << trkleny      << endl;
  cout << "diff          :" << diff[0]      << endl;
  cout << "startxpln     :" << startxpln    << endl;
  cout << "endxpln       :" << endxpln      << endl;
}


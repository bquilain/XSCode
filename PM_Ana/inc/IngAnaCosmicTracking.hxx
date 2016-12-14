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


float         angle;
int             pln;
int             mod;
int             cyc;
int           iview;
bool          sflag;
bool   clearcosmic1;

float         chi20;
float            a0;
float            b0;
float        expxy0;
int          expch0;
vector<int>     ch0;
vector<float>   pe0;
vector<float>   xy0;
vector<float> diff0;

float         chi21;
float            a1;
float            b1;
float        expxy1;
int          expch1;
vector<int>     ch1;
vector<float>   pe1;
vector<float>   xy1;
vector<float> diff1;

bool   newfidcosmic;

int        diffverz;


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
vector<ch_info> ch_info0;
vector<ch_info> ch_info1;


TFile*     fTFile;
TTree*       tree;
TH1F*      fHly[nModIng][nView][nTPL+nVETO][nCh];
int        nbinly=200;
int        minly =0;
int        maxly =100;
void MakeHist(){
  char name[300];
  for(int imod=0; imod<nModIng; imod++){
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nTPL+nVETO; ipln++){
	for(int ich=0; ich<nCh; ich++){
	  sprintf(name,"fH%02d%0d%02d%02d",imod,iview,ipln,ich);
	  fHly[imod][iview][ipln][ich] = new TH1F(name,name,
						  nbinly, minly, maxly);
	}
      }
    }
  }
}


ReconTrackBasic* frecontrack;
int nhits;
int nhits0;
int nhits1;
int nhits0high;
int nhits1high;
bool SpecialCosmicSelection_3(){
  if( nhits0 < 4 || nhits1 < 4 || nhits0 > 16 || nhits1 > 16 )
    return false;

}
bool SpecialCosmicSelection_4(){
  if( nhits0high + nhits1high > 12 )
    return true;
  else
    return false;

}
IngridBasicReconSummary* testbas    ;
IngridBasicReconSummary* trackingbas;
TRef  fIngridHit[INGRIDHIT_MAXHITS];
void AllClear(){
  nhits       = 0;
  nhits0      = 0;
  nhits1      = 0;
  nhits0high  = 0;
  nhits1high  = 0;
  testbas     -> Clear();
  trackingbas -> Clear();
  frecontrack -> Clear();
}

void CopyHit(IngridBasicReconSummary* basicrecon, 
	     IngridBasicReconSummary* testbas,
	     IngridBasicReconSummary* trackingbas,
	     int ipln, int iview){

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
      if(b->view==0){
	nhits0++;
	if( b -> pe > 8 )
	  nhits0high++;
      }
      if(b->view==1){
	nhits1++;
	if( b -> pe > 8 )
	  nhits1high++;
      }

    }
  }
  for(int ini=nhits; ini<INGRIDHIT_MAXHITS; ini++){
    fIngridHit[ini]=0;
  }
}



//###############################################
//#### for VETO analysis ########################
//###############################################

bool GetVetoPosZ(float a, float b, int mod, int ipln, float& exppos){
  if(ipln<11||ipln>14)return false;
  exppos = 1.0 * ( fVetoXY(mod, ipln) - b)/a;
  return true;
}
bool GetVetoPosXY(float a, float b, float posz, float& exppos){
  
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
float ax, bx, ay, by;

bool  startisveto;
bool  endisveto;
float da0, db0, da1, db1;
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
  da0        = retrk0 -> covx;
  db0        = retrk0 -> covy;
  da1        = retrk1 -> covx;
  db1        = retrk1 -> covy;

}
float expposz_thisview, expposxy_anotherview;
float expposz_thisview_low,expposz_thisview_high;
const static float TolOfOtherView = 22.5;
const static float TolOfThisView = 2.5;
bool FlagForVETOAnalysis(int mod, int ipln, int iview){


  if     (  ipln      ==  0 && iview == FromX && 
	    startxpln ==  0 && (startxch != 0 && startxch != 23 ) )
    return true;
  else if(  ipln      ==  0 && iview == FromY && 
	    startypln ==  0 && (startych != 0 && startych != 23 ) )
    return true;
  else if(  ipln      == 10 && iview == FromX && 
	    endxpln   == 10 && (endxch != 0 && endxch != 23) && endypln == 9)
    return true;
  else if(  ipln      == 10 && iview == FromY && 
	    endypln   == 10 && (endych != 0 && endych != 23) && endxpln == 9)
    return true;


  else if(  ipln >= nTPL ) {
    if(  ipln == RVETO ){ 
      GetVetoPosZ ( ax, bx, mod, ipln , expposz_thisview );
      GetVetoPosZ ( ax-da1, bx-db1, mod, ipln , expposz_thisview_low );
      GetVetoPosZ ( ax+da1, bx+db1, mod, ipln , expposz_thisview_high );
      GetVetoPosXY( ay, by, expposz_thisview , expposxy_anotherview );
     
    }
    if(  ipln == LVETO ){ 
      GetVetoPosZ ( ax, bx, mod, ipln , expposz_thisview );
      GetVetoPosZ ( ax-da1, bx-db1, mod, ipln , expposz_thisview_low );
      GetVetoPosZ ( ax+da1, bx+db1, mod, ipln , expposz_thisview_high );
      GetVetoPosXY( ay, by, expposz_thisview , expposxy_anotherview );
    }
    if(  ipln == BVETO ){ 
      GetVetoPosZ ( ay, by, mod, ipln , expposz_thisview );
      GetVetoPosZ ( ay-da0, by-db0, mod, ipln , expposz_thisview_low );
      GetVetoPosZ ( ay+da0, by+db0, mod, ipln , expposz_thisview_high );
      GetVetoPosXY( ax, bx, expposz_thisview , expposxy_anotherview );
    }
    if(  ipln == UVETO ){ 
      GetVetoPosZ ( ay, by, mod, ipln , expposz_thisview );
      GetVetoPosZ ( ay-da0, by-db0, mod, ipln , expposz_thisview_low );
      GetVetoPosZ ( ay+da0, by+db0, mod, ipln , expposz_thisview_high );
      GetVetoPosXY( ax, bx, expposz_thisview , expposxy_anotherview );
    }
    float maxz = max(expposz_thisview_low, expposz_thisview_high);
    float minz = min(expposz_thisview_low, expposz_thisview_high);
    /*
    if( ipln == RVETO || ipln == LVETO )
      cout << "expected z  on this    view: " << expposz_thisview     
	   << "( a: " << ax << "\t, b:" << bx << " )" << endl
	   << "expected xy on another view: " << expposxy_anotherview 
	   << "( a: " << ay << "\t, b:" << by << " )"  << endl;
    if( ipln == BVETO || ipln == UVETO )
      cout << "expected z  on this    view: " << expposz_thisview     
	   << "( a: " << ay << "\t, b:" << by << " )" << endl
	   << "expected xy on another view: " << expposxy_anotherview 
	   << "( a: " << ax << "\t, b:" << bx << " )"  << endl;
    */
    if( //StartVetoZ + TolOfThisView     < expposz_thisview &&
	//expposz_thisview               < EndVetoZ - TolOfThisView     &&
       StartVetoZ + TolOfThisView     < minz &&
       maxz                           < EndVetoZ - TolOfThisView     &&
	StartTPLXY + TolOfOtherView    < expposxy_anotherview &&
	expposxy_anotherview           < EndTPLXY - TolOfOtherView    
	)
      return true;
  }
  return false;
}
IngridBasicReconSummary*  basicrecon;

void Get(){
  chi20      = retrk0->chi2x;
  chi21      = retrk1->chi2x;

  ch_info0.clear();
  ch_info1.clear();
  mod        = basicrecon -> hitmod;

  a0         = retrk0 -> tx;  
  b0         = retrk0 -> etx;
  a1         = retrk1 -> tx;
  b1         = retrk1 -> etx;
  angle      = TMath::ATan( sqrt(a0*a0+a1*a1) ) / TMath::Pi() * 180;  	    


}

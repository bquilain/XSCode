//##### Standard C++ lib. ######
#include<iostream>
#include<sstream>
#include<fstream>
#include<math.h>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>
//##### Root Library ###########
#include <TROOT.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TMath.h>
#include <TLatex.h>
//##### INGRID Library #########
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridBasicReconSummary.h"
#include "IngridConstants.h"

float  nevt           [2][7];
float  nevt_verz_1_4  [2][7];
float  nevt_verz_5_8  [2][7];
float  nevt_verxy_3_20[2][7];
float  nevt_verxy_4_19[2][7];
float  nevt_error_verz_1_4  [2][7];
float  nevt_error_verz_5_8  [2][7];
float  nevt_error_verxy_3_20[2][7];
float  nevt_error_verxy_4_19[2][7];

float  nevt_fid_up0 [2][7]; //## with vertex z = 0~3, vertex xy = 3 ~ 20
float  nevt_fid_up1 [2][7]; //## with vertex z = 0~3, vertex xy = 4 ~ 19
float  nevt_fid_dn0 [2][7]; //## with vertex z = 4~8, vertex xy = 3 ~ 20
float  nevt_fid_dn1 [2][7]; //## with vertex z = 4~8, vertex xy = 4 ~ 19
float  nevt_ang     [2][7]; //## with vertex z = 4~8, vertex xy = 4 ~ 19
float  nevt_error         [2][7];
float  nevt_fid_up0_error [2][7]; //## with vertex z = 0~3, vertex xy = 3 ~ 20
float  nevt_fid_up1_error [2][7]; //## with vertex z = 0~3, vertex xy = 4 ~ 19
float  nevt_fid_dn0_error [2][7]; //## with vertex z = 4~8, vertex xy = 3 ~ 20
float  nevt_fid_dn1_error [2][7]; //## with vertex z = 4~8, vertex xy = 4 ~ 19
float  nevt_ang_error     [2][7]; //## with vertex z = 4~8, vertex xy = 4 ~ 19
float  xy[7];
float  xy_error[7];
void Reset(){
  for(int i=0; i<2; i++){
    for(int j=0; j<7; j++){
      xy       [j]         = ( j - 3 ) * DISTANCE_MODULE;
      xy_error [j]         = 0;
      nevt           [i][j]= 0;
      nevt_verz_1_4  [i][j]= 0;
      nevt_verz_5_8  [i][j]= 0;
      nevt_verxy_3_20[i][j]= 0;
      nevt_verxy_4_19[i][j]= 0;
      nevt           [i][j]= 0;
      nevt_fid_up0   [i][j]= 0;
      nevt_fid_up1   [i][j]= 0;
      nevt_fid_dn0   [i][j]= 0;
      nevt_fid_dn1   [i][j]= 0;
      nevt_ang       [i][j]= 0;
    }
  }
}

void Print(string str){
  ofstream wf(str.c_str());
  for(int i=0; i<2; i++){
    for(int j=0; j<7; j++){
      cout << i * 7 + j            << "\t" 
	   << nevt         [i][j]  << "\t" 
	   << nevt_verz_1_4 [i][j]  << "\t" 
	   << nevt_verz_5_8 [i][j]  << "\t" 
	   << nevt_verxy_3_20 [i][j]  << "\t" 
	   << nevt_verxy_4_19 [i][j]  << endl;

    }
  }
  wf.close();
}

void Count(IngridBasicReconSummary* basicrecon){
  int mod, id;
  if(basicrecon -> hitmod < 7){
    mod=basicrecon->hitmod; id=0;
  }
  else{
    mod = basicrecon->hitmod-7; id=1;
  }

  Bool_t  sigreg  = 
    (  basicrecon  -> layerpe > 6.5 ) &&
    (  basicrecon  -> nactpln >   2 ) &&
    (  basicrecon  -> hastrk        ) &&
    (  basicrecon  -> matchtrk      ) &&
    (  basicrecon  -> ontime        ) &&
    (!(basicrecon  -> vetowtracking)) &&
    (!(basicrecon  -> edgewtracking));

  if( sigreg ){


    nevt[id][mod]++;
    if( fabs( basicrecon -> angle ) < 20 ){
      nevt_ang[id][mod]++;
    }

    int verx  = basicrecon -> vertexx[0];
    int very  = basicrecon -> vertexy[0];
    int verxz = basicrecon -> vertexxz;
    int veryz = basicrecon -> vertexyz;

    if( 3 <= verx && verx <= 20 &&
	3 <= very && very <= 20  )
      nevt_verxy_3_20[id][mod]++;
    if( 4 <= verx && verx <= 19 &&
	4 <= very && very <= 19  )
      nevt_verxy_4_19[id][mod]++;

    if( verxz <=  4 && veryz <= 4)
      nevt_verz_1_4[id][mod]++;
    else if( verxz >=  5 && veryz >= 5 )
      nevt_verz_5_8[id][mod]++;
  }//sigreg
}


TFile*          fTFile;
TGraphErrors*   fGnevt           [2];
TGraphErrors*   fGnevt_verz_1_4  [2];
TGraphErrors*   fGnevt_verz_5_8  [2];
TGraphErrors*   fGnevt_verxy_3_20[2];
TGraphErrors*   fGnevt_verxy_4_19[2];

TH1F* fHhpro;
TH1F* fHvpro;
TH1F* fHhcenter;
TH1F* fHvcenter;


void Book(string str){

  for(int i=0; i<2; i++){
    for(int j=0; j<7; j++){
      cout << nevt[i][j] << endl;
      nevt_error            [i][j] = sqrt( nevt            [i][j] );
      nevt_error_verz_1_4   [i][j] = sqrt( nevt_verz_1_4   [i][j] );
      nevt_error_verz_5_8   [i][j] = sqrt( nevt_verz_5_8   [i][j] );
      nevt_error_verxy_3_20 [i][j] = sqrt( nevt_verxy_3_20 [i][j] );
      nevt_error_verxy_4_19 [i][j] = sqrt( nevt_verxy_4_19 [i][j] );
    }
  }

  fTFile = new TFile( str.c_str(), "recreate" );
  fHhpro = new TH1F("fHhpro",
		    "horizontal profile",
		    1000,-500,500);
  fHvpro = new TH1F("fHvpro",
		    "vertical profile",
		    1000,-500,500);
  fHhpro -> SetMaximum( nevt[1][3]*1.2 );
  fHvpro -> SetMaximum( nevt[1][3]*1.2 );
  fHhpro -> SetMinimum( 0 );
  fHvpro -> SetMinimum( 0 );
  fHhpro -> SetXTitle ("x[cm] from INGRID center");
  fHvpro -> SetXTitle ("y[cm] from INGRID center");
  fHhpro -> SetYTitle ("Number of events");
  fHvpro -> SetYTitle ("Number of events");
  fHhpro -> Write();
  fHvpro -> Write();

  for(int i=0; i<2; i++){
    fGnevt            [i] = new TGraphErrors(7, xy,    nevt[i],
					     xy_error, nevt_error        [i] );
    fGnevt_verz_1_4   [i] = new TGraphErrors(7, xy,    nevt_verz_1_4[i],
					     xy_error, nevt_error_verz_1_4[i] );
    fGnevt_verz_5_8   [i] = new TGraphErrors(7, xy,    nevt_verz_5_8[i],
					     xy_error, nevt_error_verz_5_8[i] );
    fGnevt_verxy_3_20 [i] = new TGraphErrors(7, xy,    nevt_verxy_3_20[i],
					     xy_error, nevt_error_verxy_3_20[i] );
    fGnevt_verxy_4_19 [i] = new TGraphErrors(7, xy,    nevt_verxy_4_19[i],
					     xy_error, nevt_error_verxy_4_19[i] );
    char Name[300];
    if(i==0)sprintf(Name, "h");
    if(i==1)sprintf(Name, "v");
    fGnevt            [i] -> SetName(Form("%spro", Name));
    fGnevt_verz_1_4   [i] -> SetName(Form("%sproup", Name));
    fGnevt_verz_5_8   [i] -> SetName(Form("%sprodown", Name));
    fGnevt_verxy_3_20 [i] -> SetName(Form("%sproin", Name));
    fGnevt_verxy_4_19 [i] -> SetName(Form("%sproinin", Name));
    fGnevt            [i] -> Write();
    fGnevt_verz_1_4   [i] -> Write();
    fGnevt_verz_5_8   [i] -> Write();
    fGnevt_verxy_3_20 [i] -> Write();
    fGnevt_verxy_4_19 [i] -> Write();
  
  }
}


TFile*              rfile;
TTree*              tree;
TBranch*            EvtBr;
IngridEventSummary* evt;
FileStat_t          fs;

int Read(string FileName){
  cout<<"reading "<<FileName<<"....."<<endl;
  if(gSystem->GetPathInfo(FileName.c_str(),fs)){
    cout<<"Cannot open file "<<FileName<<endl;
    return -1;
  }

  rfile = new TFile(FileName.c_str(),"read");
  tree  = (TTree*)rfile -> Get("tree");
  EvtBr = tree->GetBranch("fDefaultReco.");
  evt   = new IngridEventSummary();
  EvtBr         ->SetAddress(&evt);
  tree          ->SetBranchAddress("fDefaultReco.", &evt);
  return tree->GetEntries();
}

IngridBasicReconSummary* basicrecon;

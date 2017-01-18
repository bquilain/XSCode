//##### Standard C++ lib. ######
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
//##### Root Library ###########
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
#include <TSystem.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
//##### INGRID Library #########
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
//##### INGRID Software ########
#include "INGRID_Dimension.cxx"
#include "INGRID_BadCh_mapping.cxx"
INGRID_BadCh_mapping* IngBadChMap;

bool Is_Bad_Channel(IngridHitSummary* thit);


INGRID_Dimension* fINGRID_Dimension;

FileStat_t fs;

int ch_mapping[32]={
       4,12,28,20,
       5,13,21,29,
      30,22,14, 6,
      23,31,15, 7,
      24,16, 0, 8,
      25,17, 9, 1,
       2,10,18,26,
      11, 3,19,27
};

int pln1[16]={
	16,20,27,31, 2, 7, 8,13,
	 3, 6, 9,12,17,21,26,30
};

int pln2[16]={
	19,22,25,28, 1, 5,10,14,
	 0, 4,11,15,18,23,24,29
};

int easiroc_pln(int input, int x, int k){
  for(int i=0;i<16;i++){
	 if(pln1[i]==input){
	 	return 0+x*2;
	 } 
  }
  for(int i=0;i<16;i++){
	 if(pln2[i]==input){
	 	return 1+x*2;
	 } 
  }
  return -1;
}

int easiroc_ch(int input, int x, int k){
  if(easiroc_pln(input,x,k)%2==0){
    for(int i=0;i<8;i++){
      if(pln1[i]==input)return i    +k*8;
    }
    for(int i=8;i<12;i++){
      if(pln1[i]==input)return i-8  +k*4+40;
    }
    for(int i=12;i<16;i++){
      if(pln1[i]==input)return i-12 +k*4+60;
    }
  }
  else{
    for(int i=0;i<8;i++){
      if(pln2[i]==input)return i    +k*8;
    }
    for(int i=8;i<12;i++){
      if(pln2[i]==input)return i-8  +k*4+40;
    }
    for(int i=12;i<16;i++){
      if(pln2[i]==input)return i-12 +k*4+60;
    }
  }
  return -1;
}

int main(int argc,char *argv[]){
  Double_t NDUMMY= 4; //### Poisson mean value of # of noise hit at each module
  TROOT        root ("GUI","GUI");
  TApplication app  ("App",0,0);
  fINGRID_Dimension = new INGRID_Dimension();
  int  run_number;
  int  sub_run_number;
  int  c=-1;
  char FileName[300], Output[300];
  bool MC = false;
  while ((c = getopt(argc, argv, "f:o:")) != -1) {
    switch(c){
    case 'f':
      sprintf(FileName, "%s", optarg);
      break;
    case 'o':
      sprintf(Output, "%s", optarg);
      break;
    }
  }


  //#### Read file ######
  int bundle[4]={2,6,10,14};

  //TFile file(Form("./tree_plneff_2016_4_24_%d_%d_adc.root",bundle[0],bundle[1]));
  TFile file  (Form("../../../../data/cosmic_na_2016_4_26/tree_plneff_2016_4_24_%d_%d_adc.root",bundle[0],bundle[1]));
  TTree* tree=(TTree*)file.Get("treesum");
  Float_t pe1[40][40];
  tree->SetBranchAddress("pe",&pe1);
  //TFile file_2(Form("./tree_plneff_2016_4_24_%d_%d_adc.root",bundle[2],bundle[3]));
  TFile file_2(Form("../../../../data/cosmic_na_2016_4_26/tree_plneff_2016_4_24_%d_%d_adc.root",bundle[2],bundle[3]));
  TTree* tree_2=(TTree*)file_2.Get("treesum");
  Float_t pe2[40][40];
  tree_2->SetBranchAddress("pe",&pe2); 

  //#### Create file ######
  TFile*              wfile    = new TFile(Output,"recreate");
  TTree*              wtree    = new TTree("tree","tree");
  IngridEventSummary* wsummary = new IngridEventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary",  &wsummary,  64000,   99);

  int      startcyc = 4; //MC data is fille in cycle 4 
  int	   mod=15;
  IngridHitSummary* inghitsum;

  //for(int ievt = 0; ievt < 900000; ++ievt){
  for(int ievt = 0; ievt < 100000; ++ievt){
    if(ievt%1000==0)cout<<"\tadd noise event:"<<ievt<<endl;
    tree->GetEntry(ievt);
    tree_2->GetEntry(ievt);

    Float_t pe[40][40];
    for(int i=0;i<40;i++){
    	  pe[bundle[0]][i]=pe1[bundle[0]][i];
    	  pe[bundle[1]][i]=pe1[bundle[1]][i];
    	  pe[bundle[2]][i]=pe2[bundle[2]][i];
    	  pe[bundle[3]][i]=pe2[bundle[3]][i];
    }

    wsummary-> Clear();
    for(int k=0;k<4;k++){
    for(int i=0;i<32;i++){
	  if(pe[bundle[k]][i]<1.5){
		  continue;
	  }
	  inghitsum   = new IngridHitSummary();
	  int view, pln, ch;
	  view = 0;
	  pln  = easiroc_pln(i,2,k);
	  ch   = easiroc_ch (i,2,k);
	  //if(ch==15)std::cout << k << " " << i << " " << ch_mapping[i] << " " << pln << " " << ch << std::endl;	
	  inghitsum -> pe    = pe[bundle[k]][i];
	  //inghitsum -> lope  = tpe;
	  inghitsum -> pln   = pln;
	  inghitsum -> view  = view;
	  inghitsum -> ch    = ch;
	  inghitsum -> mod   = mod;
	  inghitsum -> cyc   = startcyc;
	  inghitsum -> time  = 0;

	  double xy, z;
	  fINGRID_Dimension -> get_pos_loli_xy( mod, view, pln, ch, &xy, &z);
	  inghitsum -> xy    = xy;
	  inghitsum ->  z    =  z;
	  inghitsum -> dummy = false;
          wsummary -> AddIngridModHit(inghitsum,mod, startcyc);
    }
    }

    wtree -> Fill();
  }


  wtree -> Write();
  wfile -> Write();
  wfile -> Close();
  //app.Run();
}
 


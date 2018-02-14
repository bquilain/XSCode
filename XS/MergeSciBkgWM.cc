#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>

#include "TChain.h"
#include "TTree.h"

//#include "SK__h1.h"
#include "Xsec.cc"

#include "INGRIDEVENTSUMMARY.h"
//#include "IngridSimVertexSummary.h"

// ML 2017/11/27
//   This macro is used for WM
//   First, run the WMMC with h2o input, then ch input
//   This macro merges the 2 output files selecting (ch in grid) + (h2o not in grid)



int main(int argc, char **argv)
{
  int c=-1;

  bool PM=false;

  char * h2o_FileName = new char[256];
  char * ch_FileName = new char[256];
  char * OutputFileName = new char[256];

  while ((c = getopt(argc, argv, "h:c:o:")) != -1) {
    switch(c){
    case 'h':
      h2o_FileName=optarg;
      break;
    case 'c':
      ch_FileName=optarg;
      break;
    case 'o':
      OutputFileName=optarg;
      break;
    }
  }

  Xsec* XS=new Xsec(PM);
  XS->Initialize();

  bool IsSand,IsAnti,IsNuE,IsBkgV,IsBkgH,IsSciBkg;

  IngridEventSummary* evt=new IngridEventSummary();  
  IngridSimVertexSummary* simver=new IngridSimVertexSummary();

  TChain *tree=new TChain("tree","tree");
  tree->Add(h2o_FileName);
  int nevt_h2o=tree->GetEntries();
  tree->Add(ch_FileName);
  tree->SetBranchAddress("fDefaultReco.",&evt);

  TTree *wtree=new TTree("tree","tree");
  IngridEventSummary* evt2 = new IngridEventSummary();
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary",
                                  &evt2,  64000,  99);
  
  int nevt=tree->GetEntries();
  cout<<nevt_h2o<<" H2O events, "<<nevt-nevt_h2o<<" CH events"<<endl;

  TChain * h1= new TChain("h1","h1");
  h1->Add(h2o_FileName);
  h1->Add(ch_FileName);

  tree->AddFriend(h1);
  TTree * wh1=(TTree*) h1->CloneTree(0);

  for (int ievt=0;ievt<nevt;ievt++){
    evt->Clear();
    tree->GetEntry(ievt);
    
    simver=evt->GetSimVertex(0);
    double position[3]={simver->xnu,simver->ynu,simver->znu};
    XS->DetermineNuType(IsSand,IsAnti,IsNuE,IsBkgH,IsBkgV,IsSciBkg,simver->nutype,simver->mod,position);

    if((ievt<nevt_h2o && IsSciBkg) || (ievt>=nevt_h2o && !IsSciBkg))  {
      //      cout<<"Event "<<ievt<<" rejected"<<endl;
      continue;
    } 
   
    evt2=evt;
    wtree->Fill();
    wh1->Fill();
  }

  cout<<wtree->GetEntries()<<" merged events"<<endl;
  cout<<wh1->GetEntries()<<" corresponding h1 entries"<<endl;

  cout<<"Writing in "<<OutputFileName<<endl;
  TFile * wfile=new TFile(OutputFileName,"recreate");
  wfile->cd();
  wtree->Write();
  wh1->Write();
  wfile->Close();

  return 0;


}

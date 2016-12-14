
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include <TError.h>
#include <TF1.h>
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
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector.h>
#include <TLegend.h>
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
#include <THStack.h>
#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"

int main(int argc, char **argv)
{
  bool Disp=false;
  int c=-1;
  bool MC=false;
  bool SystematicPE=false;
  int RandomIteration;
  string InputFileName;string OutputFileName;
  
  while ((c = getopt(argc, argv, "i:o:")) != -1) {
    switch(c){
    case 'i':
      InputFileName=optarg;
      break;
    case 'o':
      OutputFileName=optarg;
      break;
    }
  }

  cout<<"Welcome"<<endl;

  TFile * wfile = new TFile((OutputFileName).c_str(),"recreate");
  TTree*            wtree = new TTree("tree","tree");
  wtree -> SetMaxTreeSize(5000000000);
  IngridEventSummary* evt = new IngridEventSummary();
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary",
                                  &evt,  64000,  99);
  
  TFile * _file0 = new TFile((InputFileName).c_str(),"read");
  if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;
  else{cout<<"cannot open "<<_file0->GetName()<<endl; return 0;}
  TTree * tree=(TTree*) _file0->Get("tree");
  if(tree!=tree) return 0;
  int nevt=(int) tree->GetEntries();
  cout<<"Total Number Of Events="<<nevt<<endl;
  
  //   IngridEventSummary* evt = new IngridEventSummary();
  TBranch * Br=tree->GetBranch("fDefaultReco.");
  Br->SetAddress(&evt);
  tree->SetBranchAddress("fDefaultReco.",&evt);
  
  for(int ievt=0;ievt<nevt;ievt++){//loop over INGRID event (ingridsimvertex if MC, integration cycle of 580ns if data)
    if((ievt%100)==0) cout<<"Processing "<<ievt<<endl;
    evt->Clear();
    tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
    int NIngBasRec= evt->NPMAnas();
    if(NIngBasRec==0) continue;
    else cout<<NIngBasRec<<endl;
    wtree->Fill();
  }
  wfile->cd();
  wtree->Write();
  wfile->Close();
  
  return(0);
}
  

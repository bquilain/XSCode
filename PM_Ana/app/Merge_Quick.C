#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
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

//#include "root_setup.hxx"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"



//#include "IngBasicRecon.cxx"
//#include "IngEvtDisp.cxx"
//#include "IngridConstants.h"

int main(int argc, char** argv){

  char InputNeut[128];
  char InputSand[128];
  char Output[128];
  Int_t c=-1;
  int run_number;
  bool sand=false;
   while ((c = getopt(argc, argv, "h:n:s:o:")) != -1) {
    switch(c){
    case 'o':
      sprintf(Output,"%s",optarg);
      break;
    case 'n':
      sprintf(InputNeut,"%s",optarg);
     run_number=77;
      break;
    case 's':
      sprintf(InputSand,"%s",optarg);
      run_number=77;
      sand=true;
      break;
      /*    case 'h':
      optarg="1";
      cout<<"Print: ./ -n FirstFile.root -s SecondFile.root -o OutputFile.root"<<endl;
      return 0;*/
    }
   }

    /*  cout<<"Reading File "<<InputNeut<<endl;
   TFile * Neut = new TFile(InputNeut);
  TTree * TNeut = (TTree*) Neut->Get("tree");
  TBranch * BrNeut = TNeut->GetBranch("fDefaultReco.");//récupérée from le tree, pas le Tfile
  IngridEventSummary * evt = new IngridEventSummary;
  BrNeut->SetAddress(&evt);//pas directement setbrnch adress avec la branche (ne marche pas: setter l'adresse de la branche puis du tree)
  TNeut->SetBranchAddress("fDefaultReco.",&evt);
  int nevtNeut = TNeut->GetEntries();
  cout<<nevtNeut<<" Events in the Neutrino File"<<endl;
    */
  TFile * Out = new TFile(Output,"update");//ça l'ouvre également: on peut maintenant écrire dedans
  TTree * TOut = new TTree("tree","tree");
  IngridEventSummary * OutEvt = new IngridEventSummary;
  BrOut->SetAddress(&OutEvt);//pas directement setbrnch adress avec la branche (ne marche pas: setter l'adresse de la branche puis du tree)              
  TOut->SetBranchAddress("fDefaultReco.",&evt);
  int nevtOut = TOut->GetEntries();
  cout<<nevtOut<<" Events in the Neutrino File"<<endl;
  /*
  for(int ievt=0; ievt<nevtNeut; ievt++){
    if(ievt%100==0) cout<<"Number of Event "<<ievt<<endl;
    evt->Clear();
    OutEvt->Clear();
    TNeut->GetEntry(ievt);//mise automatiquement dans evt
    //puis mettre evt dans l'output Tree

    OutEvt = evt;
    TOut->Fill();//remplis la branche
    }*/

  if(sand){
    cout<<"Reading File "<<InputSand<<endl;
    TFile * Sand = new TFile(InputSand);
    TTree * TSand = (TTree*) Sand->Get("tree");
    TBranch * BrSand = TSand->GetBranch("fDefaultReco.");//récupérée from le tree, pas le Tfile
    BrSand->SetAddress(&evt);//pas directement setbrnch adress avec la branche (ne marche pas: setter l'adresse de la branche puis du tree)
    TSand->SetBranchAddress("fDefaultReco.",&evt);
    int nevtSand = TSand->GetEntries();
    cout<<nevtSand<<" Event in the Sand Muon File"<<endl;

    for(int ievt=0; ievt<nevtSand; ievt++){
      if(ievt%100==0) cout<<"Number of Event "<<ievt<<endl;
      evt->Clear();
      OutEvt->Clear();
      TSand->GetEntry(ievt);//mise automatiquement dans evt
      //puis mettre evt dans l'output Tree

      OutEvt = evt;
      TOut->Fill();//remplis la branche
    }
  }

  TOut -> Write();//écrire dans l'arbre
  Out -> Write();//écrire dans le TFile
  Out -> Close();//fermer le TFile
}

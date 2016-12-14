//##### Standard C++ lib. ######
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>
#include "math.h"
//##### Root Library ###########
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TBrowser.h>
//##### INGRID Library #########
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridBasicReconSummary.h"
#include "BeamInfoSummary.h"
#include "PMAnaSummary.h"
IngridBasicReconSummary* ingbasic;
BeamInfoSummary*         beaminfo;
IngridEventSummary*      evt;
//PMAnaSummary * ingbasic;

  //**Correction factors**//
  const float corr[14]=
    {8.71e-16,//module0
     1.06e-15,//module1
     1.13e-15,//module2
     1.22e-15,//module3
     1.32e-15,//module4
     1.22e-15,//module5
     9.23e-16,//module6
     8.12e-16,//module7
     1.00e-15,//module8
     1.05e-15,//module9
     1.22e-15,//module10
     1.19e-15,//module11
     1.09e-15,//module12
     7.14e-16}; //module13


float  getning(int imod){
  float ning=0;
  double ppb[8];
  for(int i=0;i<8;i++)
    ppb[i] = beaminfo->ct_np[4][i+1];
  /*
  for(int ibasic = 0; ibasic < evt->NPMAnas(); ibasic++){
    ingbasic     = (PMAnaSummary*)(evt->GetPMAna(ibasic));
    int cyc = ingbasic->hitcyc;

    //********  Event selection  ********
    if(ingbasic->hastrk &&           //Track reconstruction
       ingbasic->matchtrk &&         //Track matching
       ingbasic->ontime &&           //Timing cut
       !(ingbasic->vetowtracking) && //Veto cut
       !(ingbasic->edgewtracking) && //FV cut
       ingbasic->hitmod<14){         //Standard 14 module
      int mod = ingbasic->hitmod;
      if(imod!=mod) continue;
      else ning += (double)1/(1-(double)ppb[cyc-4]*corr[mod]);//Pileup correction
    }
  }
*/
  
  for(int ibasic = 0; ibasic < evt->NIngridBasicRecons(); ibasic++){
    ingbasic     = (IngridBasicReconSummary*)(evt->GetBasicRecon(ibasic));
    int cyc = ingbasic->hitcyc;

    //********  Event selection  ********
    if(ingbasic->hastrk &&           //Track reconstruction
       ingbasic->matchtrk &&         //Track matching
       ingbasic->ontime &&           //Timing cut
       !(ingbasic->vetowtracking) && //Veto cut
       !(ingbasic->edgewtracking) && //FV cut
       ingbasic->hitmod<14){         //Standard 14 module
      int mod = ingbasic->hitmod;
      if(imod!=mod) continue;
      else ning += (double)1/(1-(double)ppb[cyc-4]*corr[mod]);//Pileup correction
    }
  }
  
  return ning;
}


int main(int argc,char *argv[]){
  TROOT        root ("GUI","GUI");
  TApplication app  ("App",0,0);

  Int_t c=-1;  char FileName[300],Output[300];
  sprintf(Output,  "output.root");
  sprintf(FileName,"input.root");
  while ((c = getopt(argc, argv, "f:o:")) != -1) {
    switch(c){
    case 'f':
      sprintf(FileName,"%s",optarg);
      break;
    case 'o':
      sprintf(Output,"%s",optarg);
      break;
    }
  }

  FileStat_t fs;
  cout<<"reading and processsing "<<FileName<<"....."<<endl;
  if(gSystem->GetPathInfo(FileName,fs)){
    cout<<"Cannot open file "<<FileName<<endl;
    exit(1);
  }

  //#### Read ROOT File (Basic INGRID format)     ######
  //#####################################################
  TFile*            rfile = new TFile(FileName,"read");
  TTree*             tree = (TTree*)rfile -> Get("tree");
  TBranch*          Brutime = tree->GetBranch("utime");
  TBranch*          Brpot = tree->GetBranch("pot");
  TBranch*          BrspillID = tree->GetBranch("spillID");
  TBranch*          Brning = tree->GetBranch("ning[14]");

  float ningr[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int utimer;
  float potr;
  int spillIDr;

  Brutime -> SetAddress(&utimer);
  BrspillID -> SetAddress(&spillIDr);
  Brpot -> SetAddress(&potr);
  Brning -> SetAddress(&ningr);

  tree -> SetBranchAddress("utime",&utimer);
  tree -> SetBranchAddress("spillID",&spillIDr);
  tree -> SetBranchAddress("pot",&potr);
  tree -> SetBranchAddress("ning[14]",&ningr);
 int                nevt = (int)tree -> GetEntries();

  //#### Define ROOT File (LV format)     ######
  //############################################
  TFile*            wfile = new TFile(Output, "update");
  TTree*             wtree = (TTree*)wfile -> Get("tree");

  int utime;
  float pot;
  int spillFlag;
  float ning[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  TBranch*          BrutimeOut = wtree->GetBranch("utime");
  TBranch*          BrpotOut = wtree->GetBranch("pot");
  TBranch*          BrspillIDOut = wtree->GetBranch("spillID");
  TBranch*          BrningOut = wtree->GetBranch("ning[14]");

  BrutimeOut -> SetAddress(&utime);
  BrspillIDOut -> SetAddress(&spillFlag);
  BrpotOut -> SetAddress(&pot);
  BrningOut -> SetAddress(&ning);

  wtree -> SetBranchAddress("utime",&utime);
  wtree -> SetBranchAddress("spillID",&spillFlag);
  wtree -> SetBranchAddress("pot",&pot);
  wtree -> SetBranchAddress("ning[14]",&ning);
 

 /*
  wtree -> Branch("utime",&utime,"utime/I");
  wtree -> Branch("spillID",&spillFlag,"spillID/I");
  wtree -> Branch("pot",&pot,"pot/F");
  wtree -> Branch("ning[14]",&ning,"ning[14]/F");
  */
  cout<<"Entries="<<nevt<<endl;

  //######## Start Event Loop  #########
  for(int ievt=0; ievt<nevt; ievt++){
    if(ievt%100==0)cout << ievt << endl;
    tree      -> GetEntry(ievt);
    utime     =  utimer;
    pot       =  potr;
    spillFlag   =  spillIDr;
    //cout<<beaminfo -> spillnum<<endl;
    for(int imod=0;imod<14;imod++){
      //cout<<imod<<"   , #event="<<nIng[imod]<<endl;
      ning[imod] = ningr[imod];
    }
    wtree     -> Fill();
  }//Event Loop

  //######## Write and Close #############
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

}

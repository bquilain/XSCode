// ML 2017/11/22
// This macro is used to shrink the XSFormat files to keep only Detected Events (with at least one reconstruction)
//   The Flux and POT count histograms are kept unchanged to save the normalization of the file.

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include "INGRIDEVENTSUMMARY.h"
#include "IngridSimVertexSummary.h"
#include "PMAnaSummary.h"
#include "TTree.h"
#include "TFile.h"


#include "setup.h"

int ShrinkFile(char * fName, char * oName, int InteractionModule){
  TFile* ifile=new TFile(fName,"open");
  if(!ifile->IsOpen()) {cout<<" **no file** "; return 0;}
  TTree * itree=(TTree*) ifile->Get("tree");
  TFile * ofile=new TFile (oName,"recreate");
  TTree* otree=itree->CloneTree(0);
 
  int nevt=itree->GetEntries();
  IngridEventSummary* evt = new IngridEventSummary();
  itree->SetBranchAddress("fDefaultReco.",&evt);
  int nIngBasRec;
  IngridSimVertexSummary * SimVer;
    int IntMod;
 
  for(int i=0;i<nevt;i++){
    itree->GetEntry(i);
    nIngBasRec = evt->NPMAnas();
    SimVer = (IngridSimVertexSummary*)(evt->GetSimVertex(0));
    IntMod = SimVer->mod;
    if(i%100 == 0){
      cout<<i<<endl;
      cout<<"There is "<<nIngBasRec<<" events reconstructed" << ", Vertex mod=" << IntMod << endl;
    }
    if(IntMod == InteractionModule || nIngBasRec>0){//Remove Sand and events in INGRID only. Want to keep all events in PM or WM, even when not reconstructed for the normalization.
      otree->Fill();
    }
  }

  int oevt=otree->GetEntries();
  ofile->cd();
  otree->Write();
  ofile->Close();
  ifile->Close();
  /*
  otree->Delete();
  ofile->Delete();
  itree->Delete();
  ifile->Delete();
  */
  return oevt;
}



int main(int argc, char **argv)
{

  char * iBaseName = new char[256];
  char iName[256];
  char * oBaseName = new char[256];
  char oName[256];
  bool PM=true;
  int InteractionModule=16;//Default is PM
  
  int c=-1;
  while ((c = getopt(argc, argv, "i:o:w")) != -1) {
    switch(c){
    case 'i':
      iBaseName = optarg;
      break;
    case 'o':
      oBaseName = optarg;
      break;
    case 'w':
      PM=false;
      InteractionModule = 15;
      break;
    }
  }
      
  cout<<"shrinking "<<iBaseName<<" ...";
  int oevt=ShrinkFile(iBaseName,oBaseName,InteractionModule);
  cout<<" to "<<oevt<<" detected events: done!"<<endl;
  
  return 0;


}

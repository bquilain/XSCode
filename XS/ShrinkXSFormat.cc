// ML 2017/11/22
// This macro is used to shrink the XSFormat files to keep only Detected Events (with at least one reconstruction)
//   The Flux and POT count histograms are kept unchanged to save the normalization of the file.

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>

#include "TTree.h"
#include "TFile.h"


#include "setup.h"

int ShrinkFile(char * fName, char * oName){
  TFile* ifile=new TFile(fName,"open");
  if(!ifile->IsOpen()) {cout<<" **no file** "; return 0;}
  TTree * itree=(TTree*) ifile->Get("wtree");
  TFile * ofile=new TFile (oName,"recreate");
  TTree* otree=itree->CloneTree(0);

  TH1D* Flux=(TH1D*) ifile->Get("Flux");
  TH1D* POTCount=(TH1D*) ifile->Get("POTCount");

  int nevt=itree->GetEntries();

  bool IsDetected;
  itree->SetBranchAddress("IsDetected",&IsDetected);

  for(int i=0;i<nevt;i++){
    itree->GetEntry(i);
    if(IsDetected) otree->Fill();
  }

  int oevt=otree->GetEntries();
  ofile->cd();
  otree->Write();
  Flux->Write();
  POTCount->Write();
  ofile->Close();
  ifile->Close();

  /*otree->Delete();
  ofile->Delete();
  itree->Delete();
  ifile->Delete();
  */
  return oevt;
}



int main(int argc, char **argv)
{
  int c=-1;

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  bool PM=true;
  char DetName[2]; sprintf(DetName,(PM?"PM":"WM"));

  char iBaseName[256];
  char iName[256];
  char oBaseName[256];
  char oName[256];
  sprintf(iBaseName,"%s/XS/root_input/XSFormat_%s_redINGRID_Run1_%s_Plan.root",cINSTALLREPOSITORY,DetName,"%d");
  sprintf(oBaseName,"%s/XS/root_input/XSFormat_%s_redINGRID_Run1_%s_Plan_shrinked.root",cINSTALLREPOSITORY,DetName,"%d"); 

  for(int i=0;i<1000;i++){
    if(isBadIngridFile(i,PM)) continue;
    sprintf(iName,Form(iBaseName,i));
    sprintf(oName,Form(oBaseName,i));
    cout<<"shrinking "<<iName<<" ...";
    int oevt=ShrinkFile(iName,oName);
    cout<<" to "<<oevt<<" detected events: done!"<<endl;
  }

  return 0;


}

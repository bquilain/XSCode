
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
#include <TChain.h>
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

//#include "INGRIDEVENTSUMMARY.h"
//#include "IngridHitSummary.h"
//#include "IngridSimHitSummary.h"
//#include "IngridSimVertexSummary.h"
//#include "IngridSimParticleSummary.h"
//#include "BeamInfoSummary.h"
//#include "IngridBasicReconSummary.h"

int main(int argc, char **argv)
{
  int c=-1;
  char * NumuFileName = new char[256];
  char * NumubarFileName = new char[256];
  char * NueFileName = new char[256];
  char * WallFileName = new char[256];
  char * INGRIDHFileName = new char[256];
  char * INGRIDVFileName = new char[256];
  char * OutputFileName = new char[256];
  
  while ((c = getopt(argc, argv, "a:b:c:d:e:f:o:")) != -1) {
    switch(c){
    case 'a':
      NumuFileName=optarg;
      break;
    case 'b':
      NumubarFileName=optarg;
      break;
    case 'c':
      NueFileName=optarg;
      break;
    case 'd':
      WallFileName=optarg;
      break;
    case 'e':
      INGRIDHFileName=optarg;
      break;
    case 'f':
      INGRIDVFileName=optarg;
      break;
    case 'o':
      OutputFileName=optarg;
      break;
    }
  }

  TChain * h1 = new TChain("h1","h1");
  h1->Add(NumuFileName);
  h1->Add(NumubarFileName);
  h1->Add(NueFileName);
  h1->Add(WallFileName);
  h1->Add(INGRIDHFileName);
  h1->Add(INGRIDVFileName);

  TChain * tree = new TChain("tree","tree");
  tree->Add(NumuFileName);
  tree->Add(NumubarFileName);
  tree->Add(NueFileName);
  tree->Add(WallFileName);
  tree->Add(INGRIDHFileName);
  tree->Add(INGRIDVFileName);

  TFile * f = new TFile(OutputFileName,"recreate");
  h1->Write();
  tree->Write();
  f->Close();
  
  return(0);
}
  

#ifndef __INGRID_H__
#define __INGRID_H__
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include "ingrid.hxx"

#define NumMod 14
#define NumTPL 11
#define NumVPL 3
#define NumTFB 15
#define NumCh 48
#define NumCyc 23
const double cTdcBin=2.5;
const double cTrgTBin=10;

class INGRID_EVENT{
private:
  TFile *root_file;
  TTree *IngTree;

public:
  IngHdr *fIngHdr;
  IngEvt *fIngEvt;

  Bool_t hit[NumMod][NumTFB][NumCh][NumCyc];
  Int_t HighAdc[NumMod][NumTFB][NumCh][NumCyc];
  Int_t LowAdc[NumMod][NumTFB][NumCh][NumCyc];
  Long_t Tdc[NumMod][NumTFB][NumCh][NumCyc];

  Int_t Event;
  Int_t Spill;
  Int_t TrgId;
  Long_t TrgTime;

  INGRID_EVENT(char *file_name);
  ~INGRID_EVENT();
  
  void Book(Int_t numrun);
  void Fill();
  void Write();
  void Initialize();
};
#endif





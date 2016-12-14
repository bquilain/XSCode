#ifndef _ROOT_SETUP_H
#define _ROOT_SETUP_H

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>

#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>

class root_setup{
public:
  root_setup(){

  gStyle->SetOptFit(0000);
  gStyle->SetNdivisions(011,"X");
  gStyle->SetNdivisions(011,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.5,"Y");
  gStyle->SetLabelSize(0.06,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleSize(0.06,"Y");
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0);
  }
};


#endif

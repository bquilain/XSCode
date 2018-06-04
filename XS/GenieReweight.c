// This macro reads the neut and genie trees and computes the bin-by-bin Genie/NEUT ratio (true binning): 1st output
// It also gives the cross section curve for both generators: 2nd output

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std; 
#include <iomanip>
//#include <pair>
#include <sys/stat.h>
#include <cmath>
#include <TMinuit.h>
#include <TFitter.h>
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
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TFrame.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <THStack.h>
#include <TBox.h>
#include "setup.h"

int main(int argc, char ** argv){

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  bool water=false;  
  bool MC=true;
  int Selection = 1;//1=CC0pi, 2=CC1pi
  int c=-1;
  while ((c = getopt(argc, argv, "s:W")) != -1) {
    switch(c){
    case 's':
      Selection=atoi(optarg);
      break;      
    case 'W':
      water=true;
      break;      
    }
  }

  InitializeGlobal();
  char target[3];
  sprintf(target,(water?"H2O":"CH"));

  // opening the files
  char geniefile[256]; sprintf(geniefile,"%s/XS/genie/flat_genie12p2p4_numu_ingridFlux_%s_default_merge.root",cINSTALLREPOSITORY,target);
  char neutfile[256]; sprintf(neutfile,"%s/XS/genie/flat_neut5322_numu_ingridFlux_%s_SF_merge.root",cINSTALLREPOSITORY,target);

  TChain * tree=new TChain("FlatTree_VARS");
  tree->Add(geniefile);
  int nevt_genie=tree->GetEntries();
  tree->Add(neutfile);
  int nevt=tree->GetEntries();
  int nevt_neut=nevt-nevt_genie;
  cout<<"GENIE: "<<nevt_genie<<" events; \tNEUT: "<<nevt_neut<<" events"<<endl;  


  // branching the tree
  bool flagCC1pip,flagCC1pim;
  bool flagCC0pi;
  double fScaleFactor;
  TLorentzVector* mu_4mom=new TLorentzVector();
  tree->SetBranchAddress("flagCC0pi",&flagCC0pi);
  tree->SetBranchAddress("flagCC1pip",&flagCC1pip);
  tree->SetBranchAddress("flagCC1pim",&flagCC1pim);
  tree->SetBranchAddress("pmu_4mom",&mu_4mom);
  tree->SetBranchAddress("fScaleFactor",&fScaleFactor);

  //---------- output histograms ---------
  // to produce ratio in true binning
  TH2D * genieevents=new TH2D("genieevents","",NBinsTrueMom,0,NBinsTrueAngle,NBinsTrueAngle,0,NBinsTrueAngle);
  TH2D * neutevents=new TH2D("neutevents","",NBinsTrueMom,0,NBinsTrueAngle,NBinsTrueAngle,0,NBinsTrueAngle);
  // to produce ratio in much finer true binning
  TH2D * genieevents_fine=new TH2D("genieevents_fine","",50,0,5,90,0,90);
  TH2D * neutevents_fine=new TH2D("neutevents_fine","",50,0,5,90,0,90);

  //to produce the differential cross section curves
  TH1D * neut_dsigma_dpmu=new TH1D("neut_dsigma_dpmu","",100,0,5);
  TH1D * neut_dsigma_dthetamu=new TH1D("neut_dsigma_dthetamu","",180,0,90);
  TH1D * genie_dsigma_dpmu=new TH1D("genie_dsigma_dpmu","",100,0,5);
  TH1D * genie_dsigma_dthetamu=new TH1D("genie_dsigma_dthetamu","",180,0,90);
  

  // some global variables
  TVector3 beam(0,0,1); // z-axis is aligned with the beam in these files
  tree->GetEvent(nevt_genie);
  double AbsoluteScale=fScaleFactor;

  // loop over all events
  for(int ievt=0;ievt<nevt;ievt++){
    tree->GetEvent(ievt);
    if(ievt%100000==0) cout<<"event #"<<ievt<<endl;

    bool GENIE=(ievt<nevt_genie);

    bool IsCC1pi=flagCC1pip || flagCC1pim;
    bool IsCC0pi=flagCC0pi;
    double pmu=mu_4mom->P()/1000.;

    double theta=mu_4mom->Angle(beam)*180/TMath::Pi();    

    if((IsCC0pi && Selection == 1) || (IsCC1pi && Selection == 2)) {
      int BinTrueMom=0;
      int BinTrueAngle=0;
      
      for(int i=0;i<=NBinsTrueMom;i++){
	if(pmu<BinningTrueMom[i+1]){BinTrueMom=i;break;}
      }
      for(int i=0;i<=NBinsTrueAngle;i++){
	if(theta<BinningTrueAngle[i+1]){BinTrueAngle=i;break;}
      }
      //      cout<<pmu<<" "<<BinTrueMom<<" "<<theta<<" "<<BinTrueAngle<<endl;
      if(GENIE)  {
	genieevents->Fill(BinTrueMom,BinTrueAngle);
	genieevents_fine->Fill(pmu,theta);
      }
      else {
	neutevents->Fill(BinTrueMom,BinTrueAngle);
	neutevents_fine->Fill(pmu,theta);
      }

      if(GENIE){
	if(pmu>0.5 && pmu<5.) genie_dsigma_dthetamu->Fill(theta);
	if(theta<50) genie_dsigma_dpmu->Fill(pmu);
      }
      else {
	if(pmu>0.5 && pmu<5.) neut_dsigma_dthetamu->Fill(theta);
	if(theta<50) neut_dsigma_dpmu->Fill(pmu);
      }
    }
    
  }
  

  // postprocessing of histograms
  double geniescale=1.08*1e6/((double)2*nevt_genie);//genie overestimate the CC1pi cross section by ~8%
  double neutscale=1e6/((double)2*nevt_neut); //5e5 events/file -- normalize by nb of files
  genieevents->Scale(geniescale);
  neutevents->Scale(neutscale);
  genieevents_fine->Scale(geniescale);
  neutevents_fine->Scale(neutscale);
  neut_dsigma_dpmu->Scale(neutscale*AbsoluteScale/0.05);
  neut_dsigma_dthetamu->Scale(neutscale*AbsoluteScale/0.5);
  genie_dsigma_dpmu->Scale(geniescale*AbsoluteScale/0.05);
  genie_dsigma_dthetamu->Scale(geniescale*AbsoluteScale/0.5); 

  
  //---------- creating the output files ------------
  // file with ratio for reweighting
  TFile * ofile=new TFile(Form("true_events_fine_%s.root",target),"recreate");
  ofile->cd();
  genieevents->Write();
  neutevents->Write();
  TH2D * ratio=(TH2D*)genieevents->Clone("genie_neut_ratio");
  ratio->Divide(neutevents);
  ratio->Write();
  genieevents_fine->Write();
  neutevents_fine->Write();
  TH2D * ratio_fine=(TH2D*)genieevents_fine->Clone("genie_neut_ratio_fine");
  ratio_fine->Divide(neutevents_fine);
  ratio_fine->Write();
  ofile->Close();

  // file with differential cross section curves
  TFile * ofile2=new TFile(Form("true_xsection_%s.root",target),"recreate");
  ofile2->cd();
  neut_dsigma_dpmu->Write();
  neut_dsigma_dthetamu->Write();
  genie_dsigma_dpmu->Write();
  genie_dsigma_dthetamu->Write(); 
  ofile2->Close();

  return 0;
}

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>

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
#include <TProfile.h>
#include <TH2.h>
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
#include "setup.h"
#include "Xsec.cc"
Xsec * xs = new Xsec();

int main(int argc, char **argv){
  
  int c=-1;
  char * OutputFileName = new char[256];
  sprintf(OutputFileName,"/home/bquilain/CC0pi_XS/XS/files/PEXAngle.root");

  while ((c = getopt(argc, argv, "o:")) != -1) {
    switch(c){
    case 'o':
      OutputFileName=optarg;
      break;
    }
  }
  xs->Initialize();
  TApplication theApp("App",0,0);

  TH2D * MC_OV_ChargeDXCorrected_Ing = new TH2D("MC_OV_ChargeDXCorrected_Ing","MC_OV_ChargeDXCorrected_Ing",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * MC_OV_ChargeDXCorrected_PMIng = new TH2D("MC_OV_ChargeDXCorrected_PMIng","MC_OV_ChargeDXCorrected_PMIng",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * MC_OV_ChargeDXCorrected_PMSci = new TH2D("MC_OV_ChargeDXCorrected_PMSci","MC_OV_ChargeDXCorrected_PMSci",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_Ing = new TH2D("Data_OV_ChargeDXCorrected_Ing","Data_OV_ChargeDXCorrected_Ing",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_PMIng = new TH2D("Data_OV_ChargeDXCorrected_PMIng","Data_OV_ChargeDXCorrected_PMIng",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_PMSci = new TH2D("Data_OV_ChargeDXCorrected_PMSci","Data_OV_ChargeDXCorrected_PMSci",NBinsRecAngle,BinningRecAngle,200,0,200);

  MC_OV_ChargeDXCorrected_PMIng->Sumw2(); MC_OV_ChargeDXCorrected_PMSci->Sumw2();
  Data_OV_ChargeDXCorrected_PMIng->Sumw2(); Data_OV_ChargeDXCorrected_PMSci->Sumw2();
 
  int utime=0;
  double Charge=0;
  double ChargeCorrected=0;
  double ChargeDXCorrected=0;
  double dx=0;
  int Module=0;
  int Plane=0;
  int View=0;
  int Chan=0;
  int Channel=0;
  double Angle=0;
  double weight=0;
  bool SelectionFV=false;
  bool SelectionOV=false;
  int TrackSample=0;
  int Used=0;
  int MCSample=-1;

  TFile * fMC = new TFile("/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/MC_CalibrationPM.root");
  TTree * tMC=(TTree*) fMC->Get("wtree");
  int nevt_MC=(int) tMC->GetEntries();
  cout<<nevt_MC<<endl;
  
  TBranch* Br_MC_utime = tMC->GetBranch("utime");
  Br_MC_utime->SetAddress(&utime);
  tMC->SetBranchAddress("utime",&utime);

  TBranch* Br_MC_Angle = tMC->GetBranch("Angle");
  Br_MC_Angle->SetAddress(&Angle);
  tMC->SetBranchAddress("Angle",&Angle);

  TBranch* Br_MC_Charge = tMC->GetBranch("Charge");
  Br_MC_Charge->SetAddress(&Charge);
  tMC->SetBranchAddress("Charge",&Charge);

  TBranch* Br_MC_ChargeCorrected = tMC->GetBranch("ChargeCorrected");
  Br_MC_ChargeCorrected->SetAddress(&ChargeCorrected);
  tMC->SetBranchAddress("ChargeCorrected",&ChargeCorrected);

  TBranch* Br_MC_ChargeDXCorrected = tMC->GetBranch("ChargeDXCorrected");
  Br_MC_ChargeDXCorrected->SetAddress(&ChargeDXCorrected);
  tMC->SetBranchAddress("ChargeDXCorrected",&ChargeDXCorrected);

  TBranch* Br_MC_dx = tMC->GetBranch("dx");
  Br_MC_dx->SetAddress(&dx);
  tMC->SetBranchAddress("dx",&dx);

  TBranch* Br_MC_Channel = tMC->GetBranch("Channel");
  Br_MC_Channel->SetAddress(&Channel);
  tMC->SetBranchAddress("Channel",&Channel);

  TBranch* Br_MC_SelectionOV = tMC->GetBranch("SelectionOV");
  Br_MC_SelectionOV->SetAddress(&SelectionOV);
  tMC->SetBranchAddress("SelectionOV",&SelectionOV);

  TBranch* Br_MC_weight = tMC->GetBranch("weight");
  Br_MC_weight->SetAddress(&weight);
  tMC->SetBranchAddress("weight",&weight);

  for(int ievt=0;ievt<nevt_MC;ievt++){
    if(ievt%10000==0) cout<<ievt<<"/"<<nevt_MC<<endl;
    tMC->GetEntry(ievt);

    if(SelectionOV){
      if(Channel==0){
	MC_OV_ChargeDXCorrected_Ing->Fill(Angle,ChargeDXCorrected,weight);
      }
      else if(Channel==1){
	MC_OV_ChargeDXCorrected_PMIng->Fill(Angle,ChargeDXCorrected,weight);
      }
      else if(Channel==2){
	MC_OV_ChargeDXCorrected_PMSci->Fill(Angle,ChargeDXCorrected,weight);	
      }
    }
  }
  cout<<"end MC"<<endl;
  
  TFile * fData = new TFile("/home/bquilain/CC0pi_XS/XS/files_MCDataComparison/Data_CalibrationPM.root");
  TTree * tData=(TTree*) fData->Get("wtree");
  int nevt_Data=(int) tData->GetEntries();
  cout<<nevt_Data<<endl;
  
  TBranch* Br_Data_utime = tData->GetBranch("utime");
  Br_Data_utime->SetAddress(&utime);
  tData->SetBranchAddress("utime",&utime);

  TBranch* Br_Data_Angle = tData->GetBranch("Angle");
  Br_Data_Angle->SetAddress(&Angle);
  tData->SetBranchAddress("Angle",&Angle);

  TBranch* Br_Data_Charge = tData->GetBranch("Charge");
  Br_Data_Charge->SetAddress(&Charge);
  tData->SetBranchAddress("Charge",&Charge);

  TBranch* Br_Data_ChargeCorrected = tData->GetBranch("ChargeCorrected");
  Br_Data_ChargeCorrected->SetAddress(&ChargeCorrected);
  tData->SetBranchAddress("ChargeCorrected",&ChargeCorrected);

  TBranch* Br_Data_ChargeDXCorrected = tData->GetBranch("ChargeDXCorrected");
  Br_Data_ChargeDXCorrected->SetAddress(&ChargeDXCorrected);
  tData->SetBranchAddress("ChargeDXCorrected",&ChargeDXCorrected);

  TBranch* Br_Data_dx = tData->GetBranch("dx");
  Br_Data_dx->SetAddress(&dx);
  tData->SetBranchAddress("dx",&dx);

  TBranch* Br_Data_Channel = tData->GetBranch("Channel");
  Br_Data_Channel->SetAddress(&Channel);
  tData->SetBranchAddress("Channel",&Channel);

  TBranch* Br_Data_SelectionFV = tData->GetBranch("SelectionFV");
  Br_Data_SelectionFV->SetAddress(&SelectionFV);
  tData->SetBranchAddress("SelectionFV",&SelectionFV);

  TBranch* Br_Data_SelectionOV = tData->GetBranch("SelectionOV");
  Br_Data_SelectionOV->SetAddress(&SelectionOV);
  tData->SetBranchAddress("SelectionOV",&SelectionOV);

  TBranch* Br_Data_weight = tData->GetBranch("weight");
  Br_Data_weight->SetAddress(&weight);
  tData->SetBranchAddress("weight",&weight);

  for(int ievt=0;ievt<nevt_Data;ievt++){
    if(ievt%10000==0) cout<<ievt<<"/"<<nevt_Data<<endl;
    tData->GetEntry(ievt);
  
    if(SelectionOV){
      if(Channel==0){
	Data_OV_ChargeDXCorrected_Ing->Fill(Angle,ChargeDXCorrected,weight);	
      }
      else if(Channel==1){
	Data_OV_ChargeDXCorrected_PMIng->Fill(Angle,ChargeDXCorrected,weight);	
      }
      else if(Channel==2){
	Data_OV_ChargeDXCorrected_PMSci->Fill(Angle,ChargeDXCorrected,weight);	
      }
    }
  }
  
  
  TProfile * PEAngleData_PMIng = (TProfile*) Data_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleData_PMIng");
  TProfile * PEAngleData_PMSci = (TProfile*) Data_OV_ChargeDXCorrected_PMSci->ProfileX("PEAngleData_PMSci");
  
  TProfile * PEAngleMC_PMIng = (TProfile*) MC_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleMC_PMIng");
  TProfile * PEAngleMC_PMSci = (TProfile*) MC_OV_ChargeDXCorrected_PMSci->ProfileX("PEAngleMC_PMSci");
    
  TFile * wfile = new TFile(OutputFileName,"recreate");
  
  wfile->cd();
  
  PEAngleData_PMIng->Write();
  PEAngleData_PMSci->Write();
  PEAngleMC_PMIng->Write();
  PEAngleMC_PMSci->Write();
 

  wfile->Close();
  
  //Data_OV_ChargeDXCorrected_PMIng->Write();
  //Data_OV_ChargeDXCorrected_PMSci->Write();
  //theApp.Run();
  //TH1F * POTCount = new TH1F("POTCount","",2,0,2);
  return(0);
}
  

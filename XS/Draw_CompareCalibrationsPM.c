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

int main(int argc, char **argv){

  TApplication theApp("App",0,0);

  TH1D * MC_OV_Charge_Ing = new TH1D("MC_OV_Charge_Ing","MC_OV_Charge_Ing",200,0,200);
  TH1D * MC_OV_Charge_PMIng = new TH1D("MC_OV_Charge_PMIng","MC_OV_Charge_PMIng",200,0,200);
  TH1D * MC_OV_Charge_PMSci = new TH1D("MC_OV_Charge_PMSci","MC_OV_Charge_PMSci",200,0,200);
  TH1D * MC_FV_Charge_Ing = new TH1D("MC_FV_Charge_Ing","MC_FV_Charge_Ing",200,0,200);
  TH1D * MC_FV_Charge_PMIng = new TH1D("MC_FV_Charge_PMIng","MC_FV_Charge_PMIng",200,0,200);
  TH1D * MC_FV_Charge_PMSci = new TH1D("MC_FV_Charge_PMSci","MC_FV_Charge_PMSci",200,0,200);
  
  TH1D * Data_OV_Charge_Ing = new TH1D("Data_OV_Charge_Ing","Data_OV_Charge_Ing",200,0,200);
  TH1D * Data_OV_Charge_PMIng = new TH1D("Data_OV_Charge_PMIng","Data_OV_Charge_PMIng",200,0,200);
  TH1D * Data_OV_Charge_PMSci = new TH1D("Data_OV_Charge_PMSci","Data_OV_Charge_PMSci",200,0,200);
  TH1D * Data_FV_Charge_Ing = new TH1D("Data_FV_Charge_Ing","Data_FV_Charge_Ing",200,0,200);
  TH1D * Data_FV_Charge_PMIng = new TH1D("Data_FV_Charge_PMIng","Data_FV_Charge_PMIng",200,0,200);
  TH1D * Data_FV_Charge_PMSci = new TH1D("Data_FV_Charge_PMSci","Data_FV_Charge_PMSci",200,0,200);

  
  TH1D * MC_OV_ChargeCorrected_Ing = new TH1D("MC_OV_ChargeCorrected_Ing","MC_OV_ChargeCorrected_Ing",200,0,200);
  TH1D * MC_OV_ChargeCorrected_PMIng = new TH1D("MC_OV_ChargeCorrected_PMIng","MC_OV_ChargeCorrected_PMIng",200,0,200);
  TH1D * MC_OV_ChargeCorrected_PMSci = new TH1D("MC_OV_ChargeCorrected_PMSci","MC_OV_ChargeCorrected_PMSci",200,0,200);
  TH1D * MC_FV_ChargeCorrected_Ing = new TH1D("MC_FV_ChargeCorrected_Ing","MC_FV_ChargeCorrected_Ing",200,0,200);
  TH1D * MC_FV_ChargeCorrected_PMIng = new TH1D("MC_FV_ChargeCorrected_PMIng","MC_FV_ChargeCorrected_PMIng",200,0,200);
  TH1D * MC_FV_ChargeCorrected_PMSci = new TH1D("MC_FV_ChargeCorrected_PMSci","MC_FV_ChargeCorrected_PMSci",200,0,200);
  
  TH1D * Data_OV_ChargeCorrected_Ing = new TH1D("Data_OV_ChargeCorrected_Ing","Data_OV_ChargeCorrected_Ing",200,0,200);
  TH1D * Data_OV_ChargeCorrected_PMIng = new TH1D("Data_OV_ChargeCorrected_PMIng","Data_OV_ChargeCorrected_PMIng",200,0,200);
  TH1D * Data_OV_ChargeCorrected_PMSci = new TH1D("Data_OV_ChargeCorrected_PMSci","Data_OV_ChargeCorrected_PMSci",200,0,200);
  TH1D * Data_FV_ChargeCorrected_Ing = new TH1D("Data_FV_ChargeCorrected_Ing","Data_FV_ChargeCorrected_Ing",200,0,200);
  TH1D * Data_FV_ChargeCorrected_PMIng = new TH1D("Data_FV_ChargeCorrected_PMIng","Data_FV_ChargeCorrected_PMIng",200,0,200);
  TH1D * Data_FV_ChargeCorrected_PMSci = new TH1D("Data_FV_ChargeCorrected_PMSci","Data_FV_ChargeCorrected_PMSci",200,0,200);

 
  TH1D * MC_OV_ChargeCorrectedDX_Ing = new TH1D("MC_OV_ChargeCorrectedDX_Ing","MC_OV_ChargeCorrectedDX_Ing",200,0,200);
  TH1D * MC_OV_ChargeCorrectedDX_PMIng = new TH1D("MC_OV_ChargeCorrectedDX_PMIng","MC_OV_ChargeCorrectedDX_PMIng",200,0,200);
  TH1D * MC_OV_ChargeCorrectedDX_PMSci = new TH1D("MC_OV_ChargeCorrectedDX_PMSci","MC_OV_ChargeCorrectedDX_PMSci",200,0,200);
  TH1D * MC_FV_ChargeCorrectedDX_Ing = new TH1D("MC_FV_ChargeCorrectedDX_Ing","MC_FV_ChargeCorrectedDX_Ing",200,0,200);
  TH1D * MC_FV_ChargeCorrectedDX_PMIng = new TH1D("MC_FV_ChargeCorrectedDX_PMIng","MC_FV_ChargeCorrectedDX_PMIng",200,0,200);
  TH1D * MC_FV_ChargeCorrectedDX_PMSci = new TH1D("MC_FV_ChargeCorrectedDX_PMSci","MC_FV_ChargeCorrectedDX_PMSci",200,0,200);
  
  TH1D * Data_OV_ChargeCorrectedDX_Ing = new TH1D("Data_OV_ChargeCorrectedDX_Ing","Data_OV_ChargeCorrectedDX_Ing",200,0,200);
  TH1D * Data_OV_ChargeCorrectedDX_PMIng = new TH1D("Data_OV_ChargeCorrectedDX_PMIng","Data_OV_ChargeCorrectedDX_PMIng",200,0,200);
  TH1D * Data_OV_ChargeCorrectedDX_PMSci = new TH1D("Data_OV_ChargeCorrectedDX_PMSci","Data_OV_ChargeCorrectedDX_PMSci",200,0,200);
  TH1D * Data_FV_ChargeCorrectedDX_Ing = new TH1D("Data_FV_ChargeCorrectedDX_Ing","Data_FV_ChargeCorrectedDX_Ing",200,0,200);
  TH1D * Data_FV_ChargeCorrectedDX_PMIng = new TH1D("Data_FV_ChargeCorrectedDX_PMIng","Data_FV_ChargeCorrectedDX_PMIng",200,0,200);
  TH1D * Data_FV_ChargeCorrectedDX_PMSci = new TH1D("Data_FV_ChargeCorrectedDX_PMSci","Data_FV_ChargeCorrectedDX_PMSci",200,0,200);

  MC_OV_Charge_Ing->Sumw2(); MC_OV_Charge_PMIng->Sumw2(); MC_OV_Charge_PMSci->Sumw2();
  MC_FV_Charge_Ing->Sumw2(); MC_FV_Charge_PMIng->Sumw2(); MC_FV_Charge_PMSci->Sumw2();
  Data_OV_Charge_Ing->Sumw2(); Data_OV_Charge_PMIng->Sumw2(); Data_OV_Charge_PMSci->Sumw2();
  Data_FV_Charge_Ing->Sumw2(); Data_FV_Charge_PMIng->Sumw2(); Data_FV_Charge_PMSci->Sumw2();

  MC_OV_ChargeCorrected_Ing->Sumw2(); MC_OV_ChargeCorrected_PMIng->Sumw2(); MC_OV_ChargeCorrected_PMSci->Sumw2();
  MC_FV_ChargeCorrected_Ing->Sumw2(); MC_FV_ChargeCorrected_PMIng->Sumw2(); MC_FV_ChargeCorrected_PMSci->Sumw2();
  Data_OV_ChargeCorrected_Ing->Sumw2(); Data_OV_ChargeCorrected_PMIng->Sumw2(); Data_OV_ChargeCorrected_PMSci->Sumw2();
  Data_FV_ChargeCorrected_Ing->Sumw2(); Data_FV_ChargeCorrected_PMIng->Sumw2(); Data_FV_ChargeCorrected_PMSci->Sumw2();

  MC_OV_ChargeCorrectedDX_Ing->Sumw2(); MC_OV_ChargeCorrectedDX_PMIng->Sumw2(); MC_OV_ChargeCorrectedDX_PMSci->Sumw2();
  MC_FV_ChargeCorrectedDX_Ing->Sumw2(); MC_FV_ChargeCorrectedDX_PMIng->Sumw2(); MC_FV_ChargeCorrectedDX_PMSci->Sumw2();
  Data_OV_ChargeCorrectedDX_Ing->Sumw2(); Data_OV_ChargeCorrectedDX_PMIng->Sumw2(); Data_OV_ChargeCorrectedDX_PMSci->Sumw2();
  Data_FV_ChargeCorrectedDX_Ing->Sumw2(); Data_FV_ChargeCorrectedDX_PMIng->Sumw2(); Data_FV_ChargeCorrectedDX_PMSci->Sumw2();

  int utime=0;
  double Charge=0;
  double ChargeCorrected=0;
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

  TFile * fMC = new TFile("files_MCDataComparison/MC_CalibrationPM.root");
  TTree * tMC=(TTree*) fMC->Get("wtreee");
  int nevt_MC=(int) tMC->GetEntries();
  cout<<nevt_MC<<endl;
  
  TBranch* Br_MC_utime = tMC->GetBranch("utime");
  Br_MC_utime->SetAddress(&utime);
  tMC->SetBranchAddress("utime",&utime);

  TBranch* Br_MC_Charge = tMC->GetBranch("Charge");
  Br_MC_Charge->SetAddress(&Charge);
  tMC->SetBranchAddress("Charge",&Charge);

  TBranch* Br_MC_ChargeCorrected = tMC->GetBranch("ChargeCorrected");
  Br_MC_ChargeCorrected->SetAddress(&ChargeCorrected);
  tMC->SetBranchAddress("ChargeCorrected",&ChargeCorrected);

  TBranch* Br_MC_dx = tMC->GetBranch("dx");
  Br_MC_dx->SetAddress(&dx);
  tMC->SetBranchAddress("dx",&dx);

  TBranch* Br_MC_Channel = tMC->GetBranch("Channel");
  Br_MC_Channel->SetAddress(&Channel);
  tMC->SetBranchAddress("Channel",&Channel);

  TBranch* Br_MC_SelectionFV = tMC->GetBranch("SelectionFV");
  Br_MC_SelectionFV->SetAddress(&SelectionFV);
  tMC->SetBranchAddress("SelectionFV",&SelectionFV);

  TBranch* Br_MC_SelectionOV = tMC->GetBranch("SelectionOV");
  Br_MC_SelectionOV->SetAddress(&SelectionOV);
  tMC->SetBranchAddress("SelectionOV",&SelectionOV);

  TBranch* Br_MC_weight = tMC->GetBranch("weight");
  Br_MC_weight->SetAddress(&weight);
  tMC->SetBranchAddress("weight",&weight);

  for(int ievt=0;ievt<nevt_MC;ievt++){
    if(ievt%10000==0) cout<<ievt<<"/"<<nevt_MC<<endl;
    tMC->GetEntry(ievt);

    double ChargeCorrectedDX=ChargeCorrected;
    if(dx!=0) ChargeCorrected/=dx;
    
    if(SelectionOV){
      if(Channel==0){
	MC_OV_Charge_Ing->Fill(Charge,weight);
	MC_OV_ChargeCorrected_Ing->Fill(ChargeCorrected,weight);
	MC_OV_ChargeCorrectedDX_Ing->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==1){
	MC_OV_Charge_PMIng->Fill(Charge,weight);
	MC_OV_ChargeCorrected_PMIng->Fill(ChargeCorrected,weight);
	MC_OV_ChargeCorrectedDX_PMIng->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==2){
	MC_OV_Charge_PMSci->Fill(Charge,weight);
	MC_OV_ChargeCorrected_PMSci->Fill(ChargeCorrected,weight);
	MC_OV_ChargeCorrectedDX_PMSci->Fill(ChargeCorrectedDX,weight);	
      }
    }
    else if(SelectionFV){
      if(Channel==0){
	MC_FV_Charge_Ing->Fill(Charge,weight);
	MC_FV_ChargeCorrected_Ing->Fill(ChargeCorrected,weight);
	MC_FV_ChargeCorrectedDX_Ing->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==1){
	MC_FV_Charge_PMIng->Fill(Charge,weight);
	MC_FV_ChargeCorrected_PMIng->Fill(ChargeCorrected,weight);
	MC_FV_ChargeCorrectedDX_PMIng->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==2){
	MC_FV_Charge_PMSci->Fill(Charge,weight);
	MC_FV_ChargeCorrected_PMSci->Fill(ChargeCorrected,weight);
	MC_FV_ChargeCorrectedDX_PMSci->Fill(ChargeCorrectedDX,weight);	
      }
    }   
  }

  TFile * fData = new TFile("files_MCDataComparison/MC_CalibrationPM.root");
  TTree * tData=(TTree*) fData->Get("wtreee");
  int nevt_Data=(int) tData->GetEntries();
  cout<<nevt_Data<<endl;
  
  TBranch* Br_Data_utime = tData->GetBranch("utime");
  Br_Data_utime->SetAddress(&utime);
  tData->SetBranchAddress("utime",&utime);

  TBranch* Br_Data_Charge = tData->GetBranch("Charge");
  Br_Data_Charge->SetAddress(&Charge);
  tData->SetBranchAddress("Charge",&Charge);

  TBranch* Br_Data_ChargeCorrected = tData->GetBranch("ChargeCorrected");
  Br_Data_ChargeCorrected->SetAddress(&ChargeCorrected);
  tData->SetBranchAddress("ChargeCorrected",&ChargeCorrected);

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

    double ChargeCorrectedDX=ChargeCorrected;
    if(dx!=0) ChargeCorrected/=dx;
    
    if(SelectionOV){
      if(Channel==0){
	Data_OV_Charge_Ing->Fill(Charge,weight);
	Data_OV_ChargeCorrected_Ing->Fill(ChargeCorrected,weight);
	Data_OV_ChargeCorrectedDX_Ing->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==1){
	Data_OV_Charge_PMIng->Fill(Charge,weight);
	Data_OV_ChargeCorrected_PMIng->Fill(ChargeCorrected,weight);
	Data_OV_ChargeCorrectedDX_PMIng->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==2){
	Data_OV_Charge_PMSci->Fill(Charge,weight);
	Data_OV_ChargeCorrected_PMSci->Fill(ChargeCorrected,weight);
	Data_OV_ChargeCorrectedDX_PMSci->Fill(ChargeCorrectedDX,weight);	
      }
    }
    else if(SelectionFV){
      if(Channel==0){
	Data_FV_Charge_Ing->Fill(Charge,weight);
	Data_FV_ChargeCorrected_Ing->Fill(ChargeCorrected,weight);
	Data_FV_ChargeCorrectedDX_Ing->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==1){
	Data_FV_Charge_PMIng->Fill(Charge,weight);
	Data_FV_ChargeCorrected_PMIng->Fill(ChargeCorrected,weight);
	Data_FV_ChargeCorrectedDX_PMIng->Fill(ChargeCorrectedDX,weight);	
      }
      else if(Channel==2){
	Data_FV_Charge_PMSci->Fill(Charge,weight);
	Data_FV_ChargeCorrected_PMSci->Fill(ChargeCorrected,weight);
	Data_FV_ChargeCorrectedDX_PMSci->Fill(ChargeCorrectedDX,weight);	
      }
    }
  }

  TH1D * MC_POT = (TH1D*) fMC->Get("POTCount");
  double NPOT_MC=MC_POT->Integral();
  TH1D * Data_POT = (TH1D*) fData->Get("POTCount");
  double NPOT_Data=Data_POT->Integral();
  double Scaling=NPOT_Data; if(NPOT_MC!=0) Scaling/=NPOT_MC;


  TCanvas * cOV_Charge = new TCanvas("cOV_Charge");
  cOV_Charge->Divide(1,3);
  cOV_Charge->cd(1);
  MC_OV_Charge_Ing->Scale(Scaling);
  MC_OV_Charge_Ing->SetLineColor(kBlue);
  MC_OV_Charge_Ing->Draw();
  Data_OV_Charge_Ing->SetLineColor(kRed);
  Data_OV_Charge_Ing->Draw("same");

  
  theApp.Run();
  //TH1F * POTCount = new TH1F("POTCount","",2,0,2);
  return(0);
}
  

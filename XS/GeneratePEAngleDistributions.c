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
//#include "Xsec.cc"
//Xsec * xs = new Xsec();

int main(int argc, char **argv){
  
  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  int c=-1;
  char * OutputFileName = new char[256];
  char * BasisName = new char[64];sprintf(BasisName,"PEXAngle");
  bool PM=true;

  while ((c = getopt(argc, argv, "o:w")) != -1) {
    switch(c){
    case 'o':
      BasisName=optarg;
      break;
    case 'w':
      PM=false;
      break;
    }
  }
  
  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  sprintf(OutputFileName,"%s/XS/files/%s_%s.root",cINSTALLREPOSITORY,BasisName,DetName);
  //xs->Initialize();
  InitializeGlobal(PM);

  TApplication theApp("App",0,0);

  TH2D * MC_OV_ChargeDXCorrected_Ing = new TH2D("MC_OV_ChargeDXCorrected_Ing","MC_OV_ChargeDXCorrected_Ing",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * MC_OV_ChargeDXCorrected_PMIng = new TH2D("MC_OV_ChargeDXCorrected_PMIng","MC_OV_ChargeDXCorrected_PMIng",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * MC_OV_ChargeDXCorrected_PMSci = new TH2D("MC_OV_ChargeDXCorrected_PMSci","MC_OV_ChargeDXCorrected_PMSci",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_Ing = new TH2D("Data_OV_ChargeDXCorrected_Ing","Data_OV_ChargeDXCorrected_Ing",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_PMIng = new TH2D("Data_OV_ChargeDXCorrected_PMIng","Data_OV_ChargeDXCorrected_PMIng",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_PMSci = new TH2D("Data_OV_ChargeDXCorrected_PMSci","Data_OV_ChargeDXCorrected_PMSci",NBinsRecAngle,BinningRecAngle,200,0,200);

  TH2D * MC_OV_ChargeDXCorrected_WMPlan = new TH2D("MC_OV_ChargeDXCorrected_WMPlan","MC_OV_ChargeDXCorrected_WMPlan",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * MC_OV_ChargeDXCorrected_WMGridX = new TH2D("MC_OV_ChargeDXCorrected_WMGridX","MC_OV_ChargeDXCorrected_WMGridX",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * MC_OV_ChargeDXCorrected_WMGridY = new TH2D("MC_OV_ChargeDXCorrected_WMGridY","MC_OV_ChargeDXCorrected_WMGridY",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_WMPlan = new TH2D("Data_OV_ChargeDXCorrected_WMPlan","Data_OV_ChargeDXCorrected_WMPlan",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_WMGridX = new TH2D("Data_OV_ChargeDXCorrected_WMGridX","Data_OV_ChargeDXCorrected_WMGridX",NBinsRecAngle,BinningRecAngle,200,0,200);
  TH2D * Data_OV_ChargeDXCorrected_WMGridY = new TH2D("Data_OV_ChargeDXCorrected_WMGridY","Data_OV_ChargeDXCorrected_WMGridY",NBinsRecAngle,BinningRecAngle,200,0,200);

  MC_OV_ChargeDXCorrected_PMIng->Sumw2(); MC_OV_ChargeDXCorrected_PMSci->Sumw2();
  Data_OV_ChargeDXCorrected_PMIng->Sumw2(); Data_OV_ChargeDXCorrected_PMSci->Sumw2();

  MC_OV_ChargeDXCorrected_WMPlan->Sumw2(); MC_OV_ChargeDXCorrected_WMGridX->Sumw2(); MC_OV_ChargeDXCorrected_WMGridY->Sumw2();
  Data_OV_ChargeDXCorrected_WMPlan->Sumw2(); Data_OV_ChargeDXCorrected_WMGridX->Sumw2(); Data_OV_ChargeDXCorrected_WMGridY->Sumw2();

  MC_OV_ChargeDXCorrected_Ing->Sumw2();
  Data_OV_ChargeDXCorrected_Ing->Sumw2();

 
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

  TFile * fMC = new TFile(Form("%s/XS/files/MC_Calibration%s.root",cINSTALLREPOSITORY,DetName));
  TTree * tMC=(TTree*) fMC->Get("wtree");
  int nevt_MC=(int) tMC->GetEntries();
  cout<<nevt_MC<<endl;
  
  tMC->SetBranchAddress("utime",&utime);
  tMC->SetBranchAddress("Angle",&Angle);
  tMC->SetBranchAddress("Charge",&Charge);
  tMC->SetBranchAddress("ChargeCorrected",&ChargeCorrected);
  tMC->SetBranchAddress("ChargeDXCorrected",&ChargeDXCorrected);
  tMC->SetBranchAddress("dx",&dx);
  tMC->SetBranchAddress("Channel",&Channel);
  tMC->SetBranchAddress("SelectionOV",&SelectionOV);
  tMC->SetBranchAddress("weight",&weight);

  for(int ievt=0;ievt<nevt_MC;ievt++){
    if(ievt%10000==0) cout<<ievt<<"/"<<nevt_MC<<endl;
    tMC->GetEntry(ievt);

    if(SelectionOV){
      if(Channel==0){
	MC_OV_ChargeDXCorrected_Ing->Fill(Angle,ChargeDXCorrected,weight);
      }
      else if(PM){
	if(Channel==1) MC_OV_ChargeDXCorrected_PMIng->Fill(Angle,ChargeDXCorrected,weight);
	else if(Channel==2)MC_OV_ChargeDXCorrected_PMSci->Fill(Angle,ChargeDXCorrected,weight);	
      }
      else {
	if(Channel==1) MC_OV_ChargeDXCorrected_WMPlan->Fill(Angle,ChargeDXCorrected,weight);
	else if(Channel==2)MC_OV_ChargeDXCorrected_WMGridX->Fill(Angle,ChargeDXCorrected,weight);  
	else if(Channel==3)MC_OV_ChargeDXCorrected_WMGridY->Fill(Angle,ChargeDXCorrected,weight); 
      }
      
    }
  }
  cout<<"end MC"<<endl;
  
  TFile * fData = new TFile(Form("%s/XS/files/Data_Calibration%s.root",cINSTALLREPOSITORY,DetName));
  TTree * tData=(TTree*) fData->Get("wtree");
  int nevt_Data=(int) tData->GetEntries();
  cout<<nevt_Data<<endl;
  
  tData->SetBranchAddress("utime",&utime);
  tData->SetBranchAddress("Angle",&Angle);
  tData->SetBranchAddress("Charge",&Charge);
  tData->SetBranchAddress("ChargeCorrected",&ChargeCorrected);
  tData->SetBranchAddress("ChargeDXCorrected",&ChargeDXCorrected);
  tData->SetBranchAddress("dx",&dx);
  tData->SetBranchAddress("Channel",&Channel);
  tData->SetBranchAddress("SelectionFV",&SelectionFV);
  tData->SetBranchAddress("SelectionOV",&SelectionOV);
  tData->SetBranchAddress("weight",&weight);

  for(int ievt=0;ievt<nevt_Data;ievt++){
    if(ievt%10000==0) cout<<ievt<<"/"<<nevt_Data<<endl;
    tData->GetEntry(ievt);
  
    if(SelectionOV){
      if(Channel==0){
	Data_OV_ChargeDXCorrected_Ing->Fill(Angle,ChargeDXCorrected,weight);
      }
      else if(PM){
	if(Channel==1) Data_OV_ChargeDXCorrected_PMIng->Fill(Angle,ChargeDXCorrected,weight);
	else if(Channel==2)Data_OV_ChargeDXCorrected_PMSci->Fill(Angle,ChargeDXCorrected,weight);	
      }
      else {
	if(Channel==1) Data_OV_ChargeDXCorrected_WMPlan->Fill(Angle,ChargeDXCorrected,weight);
	else if(Channel==2)Data_OV_ChargeDXCorrected_WMGridX->Fill(Angle,ChargeDXCorrected,weight);  
	else if(Channel==3)Data_OV_ChargeDXCorrected_WMGridY->Fill(Angle,ChargeDXCorrected,weight); 
      }
      
    }
  }
  
  TProfile * PEAngleData_Ing = (TProfile*) Data_OV_ChargeDXCorrected_Ing->ProfileX("PEAngleData_Ing");
  TProfile * PEAngleMC_Ing = (TProfile*) MC_OV_ChargeDXCorrected_Ing->ProfileX("PEAngleMC_Ing");

  TProfile *PEAngleData_PMIng,*PEAngleData_PMSci,*PEAngleMC_PMIng,*PEAngleMC_PMSci;
  if(PM){
    PEAngleData_PMIng = (TProfile*) Data_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleData_PMIng");
    PEAngleData_PMSci = (TProfile*) Data_OV_ChargeDXCorrected_PMSci->ProfileX("PEAngleData_PMSci");
    PEAngleMC_PMIng = (TProfile*) MC_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleMC_PMIng");
    PEAngleMC_PMSci = (TProfile*) MC_OV_ChargeDXCorrected_PMSci->ProfileX("PEAngleMC_PMSci");
  }

  TProfile *PEAngleData_WMPlan,*PEAngleData_WMGridX,*PEAngleData_WMGridY,*PEAngleMC_WMPlan,*PEAngleMC_WMGridX,*PEAngleMC_WMGridY;
  if(!PM){
    PEAngleData_WMPlan = (TProfile*) Data_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleData_PMIng");
    PEAngleData_WMGridX = (TProfile*) Data_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleData_PMIng");
    PEAngleData_WMGridY = (TProfile*) Data_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleData_PMIng");
    PEAngleMC_WMPlan = (TProfile*) MC_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleMC_PMIng");
    PEAngleMC_WMGridX = (TProfile*) MC_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleMC_PMIng");
    PEAngleMC_WMGridY = (TProfile*) MC_OV_ChargeDXCorrected_PMIng->ProfileX("PEAngleMC_PMIng");
  }

  TFile * wfile = new TFile(OutputFileName,"recreate");
  
  wfile->cd();

  PEAngleData_Ing->Write();
  PEAngleMC_Ing->Write();

  if(PM){
    PEAngleData_PMIng->Write();
    PEAngleData_PMSci->Write();
    PEAngleMC_PMIng->Write();
    PEAngleMC_PMSci->Write();
  }
  else{
    PEAngleData_WMPlan->Write();
    PEAngleData_WMGridX->Write();
    PEAngleData_WMGridY->Write();
    PEAngleMC_WMPlan->Write();
    PEAngleMC_WMGridX->Write(); 
    PEAngleMC_WMGridY->Write();
  }

  wfile->Close();
  
  //Data_OV_ChargeDXCorrected_PMIng->Write();
  //Data_OV_ChargeDXCorrected_PMSci->Write();
  //theApp.Run();
  //TH1F * POTCount = new TH1F("POTCount","",2,0,2);
  return(0);
}
  

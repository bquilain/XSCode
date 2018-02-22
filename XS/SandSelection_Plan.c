#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
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
#include <TSpectrum.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <THStack.h>
//#include <TApplication.h>
#include <TBox.h>
#include <TVectorD.h>
#include "setup.h"
//#include "Xsec.cc"
//#include "Reconstruction.cc"
//#define DEBUG

//Xsec * XS = new Xsec();
char Type[32];
char Name[256];
char Name0[256];


void ProduceStack(TH1D * h[NFSIs], THStack * hStack){
  
  for(int fsi=0;fsi<NFSIs;fsi++){
    
    h[fsi]->GetYaxis()->SetTitleOffset(1.3);
    h[fsi]->Scale((5.86e20/1e21)*100/1000);
    
    if(fsi<3){//CC0Pi+/-/0
      if(fsi==0) h[0]->SetFillColor(kRed);
      else if(fsi==1) h[1]->SetFillColor(kOrange+10);
      else h[2]->SetFillColor(kOrange+9);
    }
    else{
      h[3]->SetFillColor(kBlue+2);//CC1Pi+/-
      h[4]->SetFillColor(kAzure+7);//CCNPi+/-
      h[5]->SetFillColor(kAzure+10);//CCpi0
      h[6]->SetFillColor(kGray);//NC
      if(fsi>6 && fsi<11){
	if(fsi==7) h[7]->SetFillColor(kYellow);//Sand
	if(fsi==8) h[8]->SetFillColor(kYellow+4);//Horizontal||Vertical Module Bkg
	if(fsi==9) h[9]->SetFillColor(kMagenta);//AntiNuMu
	if(fsi==10) h[10]->SetFillColor(kGreen);//NuE
      }
    }
    h[fsi]->SetLineColor(1);
    h[fsi]->SetLineWidth(2);
    hStack->Add(h[fsi]);
  }

}

void InitialiseTable(double ** DataSelected,double **** MCSelected,double ** BkgSelected, double ** Efficiency,double ** TotalCC0piEvent){
  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      Efficiency[c0][c1]=0;
      for(int e0=0;e0<NBinsRecMom;e0++){//loop over effect 0
	for(int e1=0;e1<NBinsRecAngle;e1++){//loop over effect 1
	  MCSelected[c0][c1][e0][e1]=0;
	  BkgSelected[e0][e1]=0;
	  DataSelected[e0][e1]=0;
	  TotalCC0piEvent[e0][e1]=0;
	}
      }
    }
  }
}

void CC0piDistributions(TChain * wtree,bool IsData,bool IsFlux,bool Plots,int Selection){
  //cout<<"hello, flux="<<Flux<<endl;
  int nevt=(int) wtree->GetEntries();
  cout<<"number of events="<<nevt<<endl;
  int err;
  int sig;

  int FSIInt;
  int Num_Int;
  int nTracks[10];
  double weight;
  bool IsFV;
  bool IsSand;
  bool IsAnti;
  bool IsBkgH;
  bool IsBkgV;
  bool IsNuE;
  bool IsDetected[10];
  bool SelectionFV[10];bool SelectionOV[10];
  double Enu;
  int nIngBasRec;
  double Errorweight;
  double TrueMomentumMuon, TrueAngleMuon;
  double POT;
  double TrackAngle[10][20];
  int TypeOfTrack[10][20];
  double CLMuon[10][20];// = new vector<double> [10][20];
  double CLMuon_KS[10][20];// = new vector<double> [10][20];
  double CLMuon_Plan[10][20];// = new vector<double> [10][20];
  double Momentum[10][20];
  double IronDistance[10][20];
  double PlasticDistance[10][20];
  int Sample[10][20];
  bool IsReconstructed[10][20];
  double ReWeight[175];
  double CriteriaAngleX[10][20];
  double CriteriaAngleY[10][20];
  double CriteriaHalfWayX[10][20];
  double CriteriaHalfWayY[10][20];
  double TrackWidth[10][20];
  int Spill;
  int GoodSpill;


  /////////////////////////////////PREPARE FOR READING/////////////////////////////////////////
  TBranch* Br_Spill = wtree->GetBranch("Spill");
  Br_Spill->SetAddress(&Spill);
  wtree->SetBranchAddress("Spill",&Spill);

  TBranch* Br_GoodSpill = wtree->GetBranch("GoodSpill");
  Br_GoodSpill->SetAddress(&GoodSpill);
  wtree->SetBranchAddress("GoodSpill",&GoodSpill);

  TBranch* Br_FSIInt = wtree->GetBranch("FSIInt");
  Br_FSIInt->SetAddress(&FSIInt);
  wtree->SetBranchAddress("FSIInt",&FSIInt);

  TBranch* Br_nIngBasRec = wtree->GetBranch("nIngBasRec");
  Br_nIngBasRec->SetAddress(&nIngBasRec);
  wtree->SetBranchAddress("nIngBasRec",&nIngBasRec);

  TBranch * Br_Num_Int = wtree->GetBranch("InteractionType");
  Br_Num_Int->SetAddress(&Num_Int);
  wtree->SetBranchAddress("InteractionType",&Num_Int);
  
  TBranch * Br_nTracks = wtree->GetBranch("nTracks[10]");
  Br_nTracks->SetAddress(&nTracks);
  wtree->SetBranchAddress("nTracks[10]",&nTracks);
 
  TBranch * Br_weight = wtree->GetBranch("weight");
  Br_weight->SetAddress(&weight);
  wtree->SetBranchAddress("weight",&weight);

  TBranch * Br_Enu = wtree->GetBranch("Enu");
  Br_Enu->SetAddress(&Enu);
  wtree->SetBranchAddress("Enu",&Enu);

  TBranch * Br_TrueMomentumMuon = wtree->GetBranch("TrueMomentumMuon");
  Br_TrueMomentumMuon->SetAddress(&TrueMomentumMuon);
  wtree->SetBranchAddress("TrueMomentumMuon",&TrueMomentumMuon);

  TBranch * Br_TrueAngleMuon = wtree->GetBranch("TrueAngleMuon");
  Br_TrueAngleMuon->SetAddress(&TrueAngleMuon);
  wtree->SetBranchAddress("TrueAngleMuon",&TrueAngleMuon);

  TBranch * Br_IsFV = wtree->GetBranch("IsFV");
  Br_IsFV->SetAddress(&IsFV);
  wtree->SetBranchAddress("IsFV",&IsFV);

  TBranch * Br_IsSand = wtree->GetBranch("IsSand");
  Br_IsSand->SetAddress(&IsSand);
  wtree->SetBranchAddress("IsSand",&IsSand);
  
  TBranch * Br_IsAnti = wtree->GetBranch("IsAnti");
  Br_IsAnti->SetAddress(&IsAnti);
  wtree->SetBranchAddress("IsAnti",&IsAnti);
  
  TBranch * Br_IsBkgH = wtree->GetBranch("IsBkgH");
  Br_IsBkgH->SetAddress(&IsBkgH);
  wtree->SetBranchAddress("IsBkgH",&IsBkgH);

  TBranch * Br_IsBkgV = wtree->GetBranch("IsBkgV");
  Br_IsBkgV->SetAddress(&IsBkgV);
  wtree->SetBranchAddress("IsBkgV",&IsBkgV);

  TBranch * Br_IsNuE = wtree->GetBranch("IsNuE");
  Br_IsNuE->SetAddress(&IsNuE);
  wtree->SetBranchAddress("IsNuE",&IsNuE);

  TBranch * Br_SelectionFV = wtree->GetBranch("SelectionFV[10]");
  Br_SelectionFV->SetAddress(&SelectionFV);
  wtree->SetBranchAddress("SelectionFV[10]",&SelectionFV);

  TBranch * Br_SelectionOV = wtree->GetBranch("SelectionOV[10]");
  Br_SelectionOV->SetAddress(&SelectionOV);
  wtree->SetBranchAddress("SelectionOV[10]",&SelectionOV);

  TBranch * Br_IsDetected = wtree->GetBranch("IsDetected[10]");
  Br_IsDetected->SetAddress(&IsDetected);
  wtree->SetBranchAddress("IsDetected[10]",&IsDetected);
  
  TBranch * Br_TrackAngle = wtree->GetBranch("TrackAngle[10][20]");
  Br_TrackAngle->SetAddress(TrackAngle);
  wtree->SetBranchAddress("TrackAngle[10][20]",TrackAngle);
 
 TBranch * Br_TypeOfTrack = wtree->GetBranch("TypeOfTrack[10][20]");
  Br_TypeOfTrack->SetAddress(TypeOfTrack);
  wtree->SetBranchAddress("TypeOfTrack[10][20]",TypeOfTrack);

  TBranch * Br_CLMuon = wtree->GetBranch("CLMuon[10][20]");
  Br_CLMuon->SetAddress(&CLMuon);
  wtree->SetBranchAddress("CLMuon[10][20]",&CLMuon);

  TBranch * Br_CLMuon_KS = wtree->GetBranch("CLMuon_KS[10][20]");
  Br_CLMuon_KS->SetAddress(&CLMuon_KS);
  wtree->SetBranchAddress("CLMuon_KS[10][20]",&CLMuon_KS);
  
  TBranch * Br_CLMuon_Plan = wtree->GetBranch("CLMuon_Plan[10][20]");
  Br_CLMuon_Plan->SetAddress(&CLMuon_Plan);
  wtree->SetBranchAddress("CLMuon_Plan[10][20]",&CLMuon_Plan);

  TBranch * Br_Momentum = wtree->GetBranch("Momentum[10][20]");
  Br_Momentum->SetAddress(Momentum);
  wtree->SetBranchAddress("Momentum[10][20]",Momentum);

  TBranch * Br_IronDistance = wtree->GetBranch("IronDistance[10][20]");
  Br_IronDistance->SetAddress(IronDistance);
  wtree->SetBranchAddress("IronDistance[10][20]",IronDistance);

  TBranch * Br_PlasticDistance = wtree->GetBranch("PlasticDistance[10][20]");
  Br_PlasticDistance->SetAddress(PlasticDistance);
  wtree->SetBranchAddress("PlasticDistance[10][20]",PlasticDistance);

  TBranch * Br_Sample = wtree->GetBranch("Sample[10][20]");
  Br_Sample->SetAddress(Sample);
  wtree->SetBranchAddress("Sample[10][20]",Sample);

  TBranch * Br_IsReconstructed = wtree->GetBranch("IsReconstructed[10][20]");
  Br_IsReconstructed->SetAddress(IsReconstructed);
  wtree->SetBranchAddress("IsReconstructed[10][20]",IsReconstructed);

  TBranch * Br_CriteriaAngleX = wtree->GetBranch("CriteriaAngleX[10][20]");
  Br_CriteriaAngleX->SetAddress(CriteriaAngleX);
  wtree->SetBranchAddress("CriteriaAngleX[10][20]",CriteriaAngleX);

  TBranch * Br_CriteriaAngleY = wtree->GetBranch("CriteriaAngleY[10][20]");
  Br_CriteriaAngleY->SetAddress(CriteriaAngleY);
  wtree->SetBranchAddress("CriteriaAngleY[10][20]",CriteriaAngleY);

  TBranch * Br_CriteriaHalfWayX = wtree->GetBranch("CriteriaHalfWayX[10][20]");
  Br_CriteriaHalfWayX->SetAddress(CriteriaHalfWayX);
  wtree->SetBranchAddress("CriteriaHalfWayX[10][20]",CriteriaHalfWayX);

  TBranch * Br_CriteriaHalfWayY = wtree->GetBranch("CriteriaHalfWayY[10][20]");
  Br_CriteriaHalfWayY->SetAddress(CriteriaHalfWayY);
  wtree->SetBranchAddress("CriteriaHalfWayY[10][20]",CriteriaHalfWayY);

  TBranch * Br_ReWeight = wtree->GetBranch("ReWeight[175]");
  Br_ReWeight->SetAddress(ReWeight);
  wtree->SetBranchAddress("ReWeight[175]",ReWeight);

  TBranch * Br_POT = wtree->GetBranch("POT");
  Br_POT->SetAddress(&POT);
  wtree->SetBranchAddress("POT",&POT);

  TBranch * Br_TrackWidth = wtree->GetBranch("TrackWidth[10][20]");
  Br_TrackWidth->SetAddress(TrackWidth);
  wtree->SetBranchAddress("TrackWidth[10][20]",TrackWidth);
  /*
  TH1D * FluxDistribution = new TH1D("FluxDistribution","",NBinsEnergyFlux,BinningEnergyFlux);
  for(int i=0;i<NBinsEnergyFlux;i++){
    if(IsFlux) FluxDistribution->SetBinContent(i+1,hNormNuMu_Binned->GetBinContent(i+1)+hNormNuMu_Binned->GetBinContent(i+1)*FluxVector[i]);
    else FluxDistribution->SetBinContent(i+1,0.);
  }
*/
  
  int Counter=0;
  int Counter2=0;
  int Counter3=0;
  double POTCount=0;

  double ** DataSelected = new double*[NBinsRecMom];
  double ** BkgSelected = new double*[NBinsRecMom];
  for(int i=0;i<NBinsRecMom;i++){
    DataSelected[i] = new double[NBinsRecAngle];
    BkgSelected[i] = new double[NBinsRecAngle];
  }

  double **** MCSelected = new double***[NBinsTrueMom];
  double ** Efficiency = new double*[NBinsTrueMom];
  double ** TotalCC0piEvent = new double*[NBinsTrueMom];
  for(int i=0;i<NBinsTrueMom;i++){
    MCSelected[i] = new double**[NBinsTrueAngle];
    Efficiency[i] = new double[NBinsTrueAngle];
    TotalCC0piEvent[i] = new double[NBinsTrueAngle];
    for(int j=0;j<NBinsTrueAngle;j++){
      MCSelected[i][j] = new double*[NBinsRecMom];
      for(int k=0;k<NBinsRecMom;k++){
	MCSelected[i][j][k] = new double[NBinsRecAngle];
      }
    }
  }


  InitialiseTable(DataSelected,MCSelected,BkgSelected,Efficiency,TotalCC0piEvent);

  TH2D * hMuCL_KSBNJ = new TH2D("hMuCL_KSBNJ","CL KS vs CL BNJ",20,0,1,20,0,1);
  TH2D * hMuCL_PlanBNJ = new TH2D("hMuCL_PlanBNJ","CL Plan vs CL BNJ",20,0,1,20,0,1);
  
  TH1D * hMuCL[NFSIs];
  TH1D * hMuCL_pvalue[NFSIs];
  TH1D * hMuCL_Lowest[NFSIs];
  TH1D * hNTracks[NFSIs];
  TH1D * hRecMom[NFSIs];
  TH1D * hRecAngle[NFSIs];
  
  if(Plots){
      char Name[256];char Title[256];

      for(int i=0;i<NFSIs;i++){
	sprintf(Name,"hMuCL%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi",i);
	hMuCL[i] = new TH1D(Name,Title,100,0,1);
	//hMuCL->Sumw2();
	hMuCL[i]->GetXaxis()->SetTitle("#mu_{CL}");    hMuCL[i]->GetYaxis()->SetTitle("Number of events");

	sprintf(Name,"hMuCL_pvalue%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi",i);
	hMuCL_pvalue[i] = new TH1D(Name,Title,100,0,1);
	//hMuCL_pvalue->Sumw2();
	hMuCL_pvalue[i]->GetXaxis()->SetTitle("#mu_{CL}");    hMuCL_pvalue[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hMuCL_Lowest%d",i);
	sprintf(Title,"Muon confidence level for the %dth-fsi, only for the lowest MuCL tracks if ntracks>1",i);
	hMuCL_Lowest[i] = new TH1D(Name,Title,20,0,1);
	//hMuCL_Lowest->Sumw2();
	hMuCL_Lowest[i]->GetXaxis()->SetTitle("#mu_{CL_Lowest}");    hMuCL_Lowest[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hNTracks%d",i);
	sprintf(Title,"Number of tracks for the %dth-fsi",i);
	hNTracks[i] = new TH1D(Name,Title,5,0,5);
	//hNTracks->Sumw2();
	hNTracks[i]->GetXaxis()->SetTitle("Number of reconstructed tracks");
	hNTracks[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecMom%d",i);
	sprintf(Title,"Distance in iron for the %dth-fsi",i);
	hRecMom[i] = new TH1D(Name,Title,20,0,100);
	//hRecMom->Sumw2();
	hRecMom[i]->GetXaxis()->SetTitle("Penetration in iron (cm)");
	hRecMom[i]->GetYaxis()->SetTitle("Number of events");
	
	sprintf(Name,"hRecAngle%d",i);
	sprintf(Title,"Reconstructed angle for the %dth-fsi",i);
	hRecAngle[i] = new TH1D(Name,Title,30,0,90);
	//hRecAngle->Sumw2();
	hRecAngle[i]->GetXaxis()->SetTitle("Angle (cm)");
	hRecAngle[i]->GetYaxis()->SetTitle("Number of events");
    }

  }    
  

  for(int ievt=0;ievt<nevt;ievt++){
    POTCount+=POT;
    //if(ievt%10000==0){cout<<ievt<<", "<<"POT processed="<<POTCount<<", xsec="<<IsXsec<<endl;
    if(ievt==(nevt-1)) cout<<"Final POT="<<POTCount<<endl;
    wtree->GetEntry(ievt);

    bool IsNuMu=((!IsSand) && (!IsAnti) && (!IsNuE) && (!IsBkgH) && (!IsBkgV));

      int BinTrueMom=0;
      int BinTrueAngle=0;

      for(int i=0;i<=NBinsTrueMom;i++){
	if(TrueMomentumMuon<BinningTrueMom[i+1]){BinTrueMom=i;break;}
      }
      for(int i=0;i<=NBinsTrueAngle;i++){
	if(TrueAngleMuon<BinningTrueAngle[i+1]){BinTrueAngle=i;break;}
      }

      if(FSIInt<3 && IsNuMu) TotalCC0piEvent[BinTrueMom][BinTrueAngle]+=weight;

      for(int irec=0;irec<nIngBasRec;irec++){

	if(!SelectionOV[irec]) continue;
	bool MuonFound=false;
	int MuonTrue;int MuonRec=0;
	int LowestMuCL=0;
	bool Trash=false;
	if(IsData) FSIInt=0;

	if(!IsDetected[irec]) continue;
	//cout<<FSIInt<<", int="<<Num_Int<<endl;
	  ////////////////////DETERMINE MUON TRACK/////////////////////////////////
	  for(int itrk=0;itrk<nTracks[irec];itrk++){
	    if(IsReconstructed[irec][itrk]){
	      if(CLMuon_Plan[irec][itrk]!=CLMuon_Plan[irec][itrk] || CLMuon_Plan[irec][itrk]==-1){ Trash=true; continue;}
	      cout<<TypeOfTrack[irec][itrk]<<endl;
	      if(TypeOfTrack[irec][itrk]==13) MuonTrue=itrk;
	      if(CLMuon_Plan[irec][itrk]>=CLMuon_Plan[irec][MuonRec]) MuonRec=itrk;
	      MuonFound=true;
	      if(CLMuon_Plan[irec][itrk]<CLMuon_Plan[irec][LowestMuCL]) LowestMuCL=itrk;
	    }
	  }
	  
	  if(!MuonFound){cout<<"whether NC, whether bug to solve"<<endl;continue;}
	  //if(!(Sample[irec][MuonRec]==5)) continue;
	  //if(!(Sample[irec][MuonRec]==3 || Sample[irec][MuonRec]==4 || Sample[irec][MuonRec]==5)) continue;

	    if(Plots){
	      hMuCL[FSIInt]->Fill(CLMuon_Plan[irec][MuonRec],weight);
	      hMuCL_KSBNJ->Fill(CLMuon_KS[irec][MuonRec],CLMuon[irec][MuonRec],weight);
	      hMuCL_PlanBNJ->Fill(CLMuon_Plan[irec][MuonRec],CLMuon[irec][MuonRec],weight);
	      
	      if(nTracks[irec]>=2) hMuCL_Lowest[FSIInt]->Fill(CLMuon_Plan[irec][LowestMuCL],weight);
	      hNTracks[FSIInt]->Fill(nTracks[irec],weight);
	      hRecMom[FSIInt]->Fill(IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio),weight);
	      hRecAngle[FSIInt]->Fill(TrackAngle[irec][MuonRec],weight);
	    }
	  ///////////////////THE ACTUAL SELECTION///////////////////////////////////////

	    int BinRecMom=0;
	    int BinRecAngle=0;
	    bool old=true;

	    for(int i=0;i<=NBinsRecMom;i++){
	      if((IronDistance[irec][MuonRec]+(PlasticDistance[irec][MuonRec]/IronCarbonRatio))<BinningRecMom[i+1]){BinRecMom=i;break;}
	    }
	    for(int i=0;i<=NBinsRecAngle;i++){
	      if(TrackAngle[irec][MuonRec]<BinningRecAngle[i+1]){BinRecAngle=i;break;}
	    }

	    DataSelected[BinRecMom][BinRecAngle]+=weight;
	    if(!IsData){
	      if(FSIInt<3 && IsNuMu){
		MCSelected[BinTrueMom][BinTrueAngle][BinRecMom][BinRecAngle]+=weight;
		Efficiency[BinTrueMom][BinTrueAngle]+=weight;
	      }
	      else{
		BkgSelected[BinRecMom][BinRecAngle]+=weight;
	      }
	    }
	}
      }

  for(int c0=0;c0<NBinsTrueMom;c0++){//loop over cause 0
    for(int c1=0;c1<NBinsTrueAngle;c1++){//loop over cause 1
      if(!IsData && TotalCC0piEvent[c0][c1]!=0) Efficiency[c0][c1]/=TotalCC0piEvent[c0][c1];
    }
  }


  if(Plots){
    if(!IsData){
      THStack * Stack_MuCL = new THStack("Stack_MuCL","");
      THStack * Stack_MuCL_pvalue = new THStack("Stack_MuCL_pvalue","");
      THStack * Stack_MuCL_Lowest = new THStack("Stack_MuCL_Lowest","");
      THStack * Stack_NTracks = new THStack("Stack_NTracks","");
      THStack * Stack_RecMom = new THStack("Stack_RecMom","");
      THStack * Stack_RecAngle = new THStack("Stack_RecAngle","");

      for(int fsi=0;fsi<NFSIs;fsi++){
	double Integral=hMuCL[fsi]->Integral();
	for(int ibinx=1;ibinx<=hMuCL_pvalue[fsi]->GetNbinsX();ibinx++){
	  double Value=hMuCL[fsi]->Integral(ibinx,hMuCL_pvalue[fsi]->GetNbinsX());
	  if(Integral!=0) Value/=Integral;
	  else Value=1;
	  hMuCL_pvalue[fsi]->SetBinContent(ibinx,1.-Value);
	  //double Error=;
	}
      }
        
      ProduceStack(hMuCL,Stack_MuCL);
      ProduceStack(hMuCL_pvalue,Stack_MuCL_pvalue);
      ProduceStack(hMuCL_Lowest,Stack_MuCL_Lowest);
      ProduceStack(hNTracks,Stack_NTracks);
      ProduceStack(hRecMom,Stack_RecMom);
      ProduceStack(hRecAngle,Stack_RecAngle);

      TFile * fMC = new TFile("plots/SandMCPlots_Plan.root","RECREATE");
      Stack_MuCL->Write();
      Stack_MuCL_pvalue->Write();
      Stack_MuCL_Lowest->Write();
      Stack_NTracks->Write();
      Stack_RecMom->Write();
      Stack_RecAngle->Write();
      hMuCL_KSBNJ->Write();
      hMuCL_PlanBNJ->Write(); 
      fMC->Close();
    }
    else{
      TFile * fData = new TFile("plots/SandDataPlots.root","RECREATE");
      hMuCL[0]->Write("Data_MuCL");
      hMuCL_Lowest[0]->Write("Data_MuCL_Lowest");
      hNTracks[0]->Write("Data_NTracks");
      hRecMom[0]->Write("Data_RecMom");
      hRecAngle[0]->Write("Data_RecAngle");
      fData->Close();
    }
    //ProduceStack(hNTracks,Stack_NTracks);
    //ProduceStack(hRecMom,Stack_RecMom);
    //ProduceStack(hRecAngle,Stack_RecAngle);

  }

  for(int i=0;i<NBinsRecMom;i++){
    delete DataSelected[i];
    delete BkgSelected[i];
  }
    delete DataSelected;
  delete BkgSelected;
  
  for(int i=0;i<NBinsTrueMom;i++){
    for(int j=0;j<NBinsTrueAngle;j++){
      for(int k=0;k<NBinsRecMom;k++){
	delete MCSelected[i][j][k];
      }
      delete MCSelected[i][j];
    }
    delete MCSelected[i];
    delete Efficiency[i];
    delete TotalCC0piEvent[i];
  }
  delete MCSelected;
  delete Efficiency;
  delete TotalCC0piEvent;

}

int main(int argc, char ** argv){

  cout<<"welcome, here is the CC0pi selection"<<endl;
  TApplication theApp("App",0,0);
  char fName[256];
  bool Data=false;
  int Sample=0;
  
  int c=-1;
  while ((c = getopt(argc, argv, "d")) != -1) {
    switch(c){
    case 'd':
      Data=true;
      break;
    case 's':
      Sample=atoi(optarg);//0 for sand muon, 1 for CC0pi
      break;
    }
  }
  
  cout<<"welcome"<<endl;
  //  XS->Xsec::Initialize();
  InitializeGlobal();

   TChain * chain = new TChain("wtree");
   if(Data){
     for(int r=14510;r<=14510;r++){
       for(int s=0;s=300;s++){
	 sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_%08d_%04d.root",r,s);
	 cout<<fName<<endl;
	 chain->Add(fName);
       }
     }
   }
   else{
    
     for(int i=1;i<=NMCfiles;i++){
       sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Sand_Run1_%d_Plan.root",i);
       //sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Sand_Run1_%d_Plan_RandomPE.root",i);
       cout<<fName<<endl;
       chain->Add(fName);
     }/*
     for(int i=1;i<=NMCfiles;i++){
       sprintf(fName,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Run1_%d.root",i);
       cout<<fName<<endl;
       chain->Add(fName);
       }  */

#ifdef DEBUG
     TH1D * hTest_PMIng;  TH1D * hTest_PMSci;  TH1D * hTest_Ing;
     TH1D * hTestImmediate_PMIng;  TH1D * hTestImmediate_PMSci;  TH1D * hTestImmediate_Ing;
     TF1 * f_PMIng_Plan;  TF1 * f_PMSci_Plan;  TF1 * f_Ing_Plan;
     TF1 * CL_PMIng_Plan;  TF1 * CL_PMSci_Plan;  TF1 * CL_Ing_Plan;
     TH2D * CLTest_PMIng;  TH2D * CLTest_PMSci;  TH2D * CLTest_Ing;

     char Name[256];
     TFile * f;
  for(int n=1;n<=NMCfiles;n++){
    sprintf(Name,"/home/bquilain/CC0pi_XS/XS/root_input/XSFormat_Sand_Run1_%d_Plan_RandomPE.root",n);
    f = new TFile(Name);
    if(f->IsOpen()) cout << f->GetName() <<" is open"<< endl ;
    else{cout<<"Not able to read the file"<<endl; continue;}
    
    if(n==1){
      hTest_PMIng = (TH1D*) ((f->Get("hTest_PMIng"))->Clone("hTest_PMIng"));
      f_PMIng_Plan = (TF1*) ((f->Get("f_PMIng_Plan"))->Clone("f_PMIng_Plan"));

      hTest_PMSci = (TH1D*) ((f->Get("hTest_PMSci"))->Clone("hTest_PMSci"));
      f_PMSci_Plan = (TF1*) ((f->Get("f_PMSci_Plan"))->Clone("f_PMSci_Plan"));
      
      hTest_Ing = (TH1D*) ((f->Get("hTest_Ing"))->Clone("hTest_Ing"));
      f_Ing_Plan = (TF1*) ((f->Get("f_Ing_Plan"))->Clone("f_Ing_Plan"));

      
      hTestImmediate_PMIng = (TH1D*) ((f->Get("hTestImmediate_PMIng"))->Clone("hTestImmediate_PMIng"));
      f_PMIng_Plan = (TF1*) ((f->Get("f_PMIng_Plan"))->Clone("f_PMIng_Plan"));

      hTestImmediate_PMSci = (TH1D*) ((f->Get("hTestImmediate_PMSci"))->Clone("hTestImmediate_PMSci"));
      f_PMSci_Plan = (TF1*) ((f->Get("f_PMSci_Plan"))->Clone("f_PMSci_Plan"));
      
      hTestImmediate_Ing = (TH1D*) ((f->Get("hTestImmediate_Ing"))->Clone("hTestImmediate_Ing"));
      f_Ing_Plan = (TF1*) ((f->Get("f_Ing_Plan"))->Clone("f_Ing_Plan"));
      
      CLTest_PMIng = (TH2D*) ((f->Get("CLTest_PMIng"))->Clone("CLTest_PMIng"));
      CL_PMIng_Plan = (TF1*) ((f->Get("CL_PMIng_Plan"))->Clone("CL_PMIng_Plan"));

      CLTest_PMSci = (TH2D*) ((f->Get("CLTest_PMSci"))->Clone("CLTest_PMSci"));
      CL_PMSci_Plan = (TF1*) ((f->Get("CL_PMSci_Plan"))->Clone("CL_PMSci_Plan"));
      
      CLTest_Ing = (TH2D*) ((f->Get("CLTest_Ing"))->Clone("CLTest_Ing"));
      CL_Ing_Plan = (TF1*) ((f->Get("CL_Ing_Plan"))->Clone("CL_Ing_Plan"));

    }
    else{
      hTest_PMIng->Add((TH1D*) (f->Get("hTest_PMIng")));
      hTest_PMSci->Add((TH1D*) (f->Get("hTest_PMSci")));
      hTest_Ing->Add((TH1D*) (f->Get("hTest_Ing")));

      hTestImmediate_PMIng->Add((TH1D*) (f->Get("hTestImmediate_PMIng")));
      hTestImmediate_PMSci->Add((TH1D*) (f->Get("hTestImmediate_PMSci")));
      hTestImmediate_Ing->Add((TH1D*) (f->Get("hTestImmediate_Ing")));

      CLTest_PMIng->Add((TH2D*) (f->Get("CLTest_PMIng")));
      CLTest_PMSci->Add((TH2D*) (f->Get("CLTest_PMSci")));
      CLTest_Ing->Add((TH2D*) (f->Get("CLTest_Ing")));

    }
    //f->Close();
    //f->Delete();
  }
  
  hTest_PMIng->Rebin(2);
  hTest_PMSci->Rebin(2);
  hTest_Ing->Rebin(2);

  hTest_PMIng->Scale(1./hTest_PMIng->GetMaximum());
  hTest_PMSci->Scale(1./hTest_PMSci->GetMaximum());
  hTest_Ing->Scale(1./hTest_Ing->GetMaximum());

  hTestImmediate_PMIng->Rebin(2);
  hTestImmediate_PMSci->Rebin(2);
  hTestImmediate_Ing->Rebin(2);

  hTestImmediate_PMIng->Scale(1./hTestImmediate_PMIng->GetMaximum());
  hTestImmediate_PMSci->Scale(1./hTestImmediate_PMSci->GetMaximum());
  hTestImmediate_Ing->Scale(1./hTestImmediate_Ing->GetMaximum());

  TCanvas * cPMIng = new TCanvas("cPMIng","PM INGRID type");
  hTest_PMIng->Draw("HIST");
  hTestImmediate_PMIng->SetLineColor(kBlue);
  hTestImmediate_PMIng->Draw("HISTsame");
  f_PMIng_Plan->SetLineColor(2);
  f_PMIng_Plan->Draw("same");

  TCanvas * cPMSci = new TCanvas("cPMSci","PM SciBar type");
  hTest_PMSci->Draw("HIST");
  hTestImmediate_PMSci->SetLineColor(kBlue);
  hTestImmediate_PMSci->Draw("HISTsame");
  f_PMSci_Plan->SetLineColor(2);
  f_PMSci_Plan->Draw("same");

  TCanvas * cIng = new TCanvas("cIng","INGRID type");
  hTest_Ing->Draw("HIST");
  hTestImmediate_Ing->SetLineColor(kBlue);
  hTestImmediate_Ing->Draw("HISTsame");
  f_Ing_Plan->SetLineColor(2);
  f_Ing_Plan->Draw("same");

  /*
  CLTest_PMIng->Rebin(2);
  CLTest_PMSci->Rebin(2);
  CLTest_Ing->Rebin(2);

  CLTest_PMIng->Scale(1./CLTest_PMIng->GetMaximum());
  CLTest_PMSci->Scale(1./CLTest_PMSci->GetMaximum());
  CLTest_Ing->Scale(1./CLTest_Ing->GetMaximum());
  */

  TProfile * pCLTest_PMIng = (TProfile*) CLTest_PMIng->ProfileX("pCLTest_PMIng",1,CLTest_PMIng->GetNbinsY());
  TCanvas * cCumulativePMIng = new TCanvas("cCumulativePMIng","PM INGRID type,CL");
  pCLTest_PMIng->Draw("HIST");
  CL_PMIng_Plan->SetLineColor(2);
  CL_PMIng_Plan->Draw("same");

    TProfile * pCLTest_PMSci = (TProfile*) CLTest_PMSci->ProfileX("pCLTest_PMSci",1,CLTest_PMSci->GetNbinsY());
  TCanvas * cCumulativePMSci = new TCanvas("cCumulativePMSci","PM INGRID type,CL");
  pCLTest_PMSci->Draw("HIST");
  CL_PMSci_Plan->SetLineColor(2);
  CL_PMSci_Plan->Draw("same");

    TProfile * pCLTest_Ing = (TProfile*) CLTest_Ing->ProfileX("pCLTest_Ing",1,CLTest_Ing->GetNbinsY());
  TCanvas * cCumulativeIng = new TCanvas("cCumulativeIng"," INGRID type,CL");
  pCLTest_Ing->Draw("HIST");
  CL_Ing_Plan->SetLineColor(2);
  CL_Ing_Plan->Draw("same");





  
    TH1D * projectedCLTest_PMIng = (TH1D*) CLTest_PMIng->ProjectionY("projectedCLTest_PMIng",1,CLTest_PMIng->GetNbinsY());
  TCanvas * cCLPMIng = new TCanvas("cCLPMIng","PM INGRID type,CL");
  projectedCLTest_PMIng->Draw("HIST");

    TH1D * projectedCLTest_PMSci = (TH1D*) CLTest_PMSci->ProjectionY("projectedCLTest_PMSci",1,CLTest_PMSci->GetNbinsY());
  TCanvas * cCLPMSci = new TCanvas("cCLPMSci","PM INGRID type,CL");
  projectedCLTest_PMSci->Draw("HIST");

    TH1D * projectedCLTest_Ing = (TH1D*) CLTest_Ing->ProjectionY("projectedCLTest_Ing",1,CLTest_Ing->GetNbinsY());
  TCanvas * cCLIng = new TCanvas("cCLIng"," INGRID type,CL");
  projectedCLTest_Ing->Draw("HIST");


  TH1D * projectedPETest_PMIng = (TH1D*) CLTest_PMIng->ProjectionX("projectedPETest_PMIng",1,CLTest_PMIng->GetNbinsX());
  TCanvas * cPEPMIng = new TCanvas("cPEPMIng","PM INGRID type,PE");
  projectedPETest_PMIng->Draw("HIST");

    TH1D * projectedPETest_PMSci = (TH1D*) CLTest_PMSci->ProjectionX("projectedPETest_PMSci",1,CLTest_PMSci->GetNbinsX());
  TCanvas * cPEPMSci = new TCanvas("cPEPMSci","PM INGRID type,PE");
  projectedPETest_PMSci->Draw("HIST");

    TH1D * projectedPETest_Ing = (TH1D*) CLTest_Ing->ProjectionX("projectedPETest_Ing",1,CLTest_Ing->GetNbinsX());
  TCanvas * cPEIng = new TCanvas("cPEIng"," INGRID type,PE");
  projectedPETest_Ing->Draw("HIST");

#endif
  
  }

   CC0piDistributions(chain,Data,false,true,Sample);
#ifdef DEBUG
  theApp.Run();
#endif

   return 0;
}

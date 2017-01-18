#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <utility>

#include <TFile.h>
#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TString.h>
#include <TRootEmbeddedCanvas.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TArrow.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TBox.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TColor.h>
#include <TGFileDialog.h>
#include "TSystem.h"
#include "TString.h"
#include "TGClient.h"
#include "TGWindow.h"
#include "TClass.h"
#include "THashList.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TEnv.h"
#include "TVirtualX.h"
#include "TImage.h"
#include "event_viewer.h"

// INGRID (PM) libraries
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h)
  : TGMainFrame(p, w, h)
{
  myfile=new TFile("PM_MC_Beam1_BirksCorrected_wNoise_ana.root");
  if(myfile->IsOpen()) cout << myfile->GetName() <<" is open"<< endl ;
  IsData=kFALSE;
  irec=0;
  tree=(TTree*) myfile->Get("tree");
  nevt=(int) tree->GetEntries();
  cout<<"Number of events = "<<nevt<<endl;   
  evt = new IngridEventSummary();
  Br=tree->GetBranch("fDefaultReco.");
  Br->SetAddress(&evt);
  tree->SetBranchAddress("fDefaultReco.",&evt);
  
  ievt=0;
  //a = new TH1D("a","a",100,-5,5);

  // bottom horizontal frame
  fHor1 = new TGHorizontalFrame(this, 60, 20, kFixedWidth);
  fExit = new TGTextButton(fHor1, "&Exit", "gApplication->Terminate(0)");
  fHor1->AddFrame(fExit, new TGLayoutHints(kLHintsCenterY | kLHintsRight 
					   , 4, 4, 4, 4));
  
  fOpenFile = new TGTextButton(fHor1, "&Open file");
  fOpenFile->Connect("Clicked()","MyMainFrame",this,"OpenFile()");
  fHor1->AddFrame(fOpenFile, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));

  fSaveAs = new TGTextButton(fHor1, "&Save as...");
  fSaveAs->Connect("Clicked()","MyMainFrame",this,"SaveAs()");
  fHor1->AddFrame(fSaveAs, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 4, 4, 4));
  //fHor1->AddFrame(fExit, new TGLayoutHints(kLHintsCenterY | kLHintsRight |
  //					   , 4, 4, 4, 4));
  AddFrame(fHor1, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 2, 2, 5, 2)); 

  /////////////////////////////////////////////
  
  fHor3 = new TGHorizontalFrame(this, 1300, 10, kFixedWidth);
  //============================================================ select evt
  fGframe = new TGGroupFrame(fHor3, "Selected event number");

  fNumber = new TGNumberEntry(fGframe, 0, 9, 999,
			      TGNumberFormat::kNESInteger,
			      TGNumberFormat::kNEANonNegative,
			      TGNumberFormat::kNELLimitMinMax,
			      0, 10000);
  fNumber->SetLimitValues(0,nevt-1);

  fNumber->Connect("ValueSet(Long_t)", "MyMainFrame", this, "GoNextEvent()");
  
  fGframe->AddFrame(fNumber, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));

  fHor3->AddFrame(fGframe, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  
  //==================================================== select reco

  fGframe4 = new TGGroupFrame(fHor3, "Selected Reconstruction");

  fHor5 = new TGHorizontalFrame(fGframe4);

  fNextReco = new TGTextButton(fHor5, "&Next Reco");
  fNextReco->Connect("Clicked()","MyMainFrame",this,"GoNextReco()");
  fNextReco->SetState(kButtonDisabled,kTRUE);

  fLabel = new TGLabel(fHor5, "0");

  fHor5->AddFrame(fNextReco,new TGLayoutHints(kLHintsCenterY, 4, 4, 4, 4)); 
  fHor5->AddFrame(fLabel, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 4, 4, 4, 4));

  fGframe4->AddFrame(fHor5, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2,2,2,2));

  fHor3->AddFrame(fGframe4, new TGLayoutHints(kLHintsCenterY, 80, 4, 4, 4));

  //================================================= MC Drawing options

  fGframeMCOptions = new TGGroupFrame(fHor3, "MC Drawing options");
  
  fHor4 = new TGHorizontalFrame(fGframeMCOptions);

  fVer1 = new TGVerticalFrame(fHor4);
  fHor4->AddFrame(fVer1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fButOption1 = new TGCheckButton(fVer1, "&Show vertex");
  fButOption1->Connect("Clicked()","MyMainFrame",this,"DrawOption1()");
  fButOption2 = new TGCheckButton(fVer1, "&Show nu beam");
  fButOption2->Connect("Clicked()","MyMainFrame",this,"DrawOption2()");
  fVer1->AddFrame(fButOption1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fVer1->AddFrame(fButOption2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));

  fVer2 = new TGVerticalFrame(fHor4);
  fHor4->AddFrame(fVer2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fButOption3 = new TGCheckButton(fVer2, "&Draw charged track(s)");
  fButOption3->Connect("Clicked()","MyMainFrame",this,"DrawOption3()");
  fButOption4 = new TGCheckButton(fVer2, "&Draw neutral track(s)");
  fButOption4->Connect("Clicked()","MyMainFrame",this,"DrawOption4()");
  fVer2->AddFrame(fButOption3, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fVer2->AddFrame(fButOption4, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));

  /*  fVer3 = new TGVerticalFrame(fHor4);
  fHor4->AddFrame(fVer3, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fButOption5 = new TGCheckButton(fVer3, "&Show extra hits");
  fButOption5->Connect("Clicked()","MyMainFrame",this,"DrawOption5()");
  fVer3->AddFrame(fButOption5, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));*/
		 
  fGframeMCOptions->AddFrame(fHor4, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,  15, 4, 4, 4)); 
  fHor3->AddFrame(fGframeMCOptions, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 80, 4, 4, 4)); 
  if(IsData){
    fButOption1->SetEnabled(kFALSE);
    fButOption2->SetEnabled(kFALSE);
    fButOption3->SetEnabled(kFALSE);
    fButOption4->SetEnabled(kFALSE);
  }
 

  //================================================ Reconstruction options

  fGframeROptions = new TGGroupFrame(fHor3, "Reconstruction options");

  fHor6 = new TGHorizontalFrame(fGframeROptions);

  fVer3 = new TGVerticalFrame(fHor6);
  fHor6->AddFrame(fVer3, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fButOption5 = new TGCheckButton(fVer3, "&Show extra hits");
  fButOption5->Connect("Clicked()","MyMainFrame",this,"DrawOption5()");
  fButOption6 = new TGCheckButton(fVer3, "&Show reconstructed vertex");
  fButOption6->Connect("Clicked()","MyMainFrame",this,"DrawOption6()");
  fVer3->AddFrame(fButOption5, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fVer3->AddFrame(fButOption6, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));

  /*fVer4 = new TGVerticalFrame(fHor6);
  fHor6->AddFrame(fVer4, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fButOption7 = new TGCheckButton(fVer4, "&Use veto cut");
  fButOption7->Connect("Clicked()","MyMainFrame",this,"DrawOption7()");
  fButOption8 = new TGCheckButton(fVer4, "&Use FV cut");
  fButOption8->Connect("Clicked()","MyMainFrame",this,"DrawOption8()");
  fVer4->AddFrame(fButOption7, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  fVer4->AddFrame(fButOption8, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 4, 4));
  */
  fGframeROptions->AddFrame(fHor6, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,  15, 4, 4, 4)); 
  fHor3->AddFrame(fGframeROptions, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 80, 4, 4, 4));  


  AddFrame(fHor3, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 4, 4, 4, 4)); 



  /////////////////////////////////////////////  histograms' frame



  fHor2 = new TGHorizontalFrame(this, 1300, 900, kFixedWidth);
  fGframe2 = new TGGroupFrame(fHor2, "Top view - Proton Module");
  fHor2->AddFrame(fGframe2, new TGLayoutHints(kLHintsTop | kLHintsLeft |
					      kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));
  
  fCanvasTopPM = new TRootEmbeddedCanvas("CanvasTopPM",fGframe2,300,300);
  fGframe2->AddFrame(fCanvasTopPM, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,4,4,12,12));
  fGframe3 = new TGGroupFrame(fHor2, "Top view - Indrid 3rd module");
  fHor2->AddFrame(fGframe3, new TGLayoutHints(kLHintsTop | kLHintsRight |
					      kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));
  
  fCanvasTopI = new TRootEmbeddedCanvas("CanvasTopI",fGframe3,300,300);
  fGframe3->AddFrame(fCanvasTopI, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,4,4,12,12));

  AddFrame(fHor2,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));

  //--------------------------------------------------------------------
  fHor7 = new TGHorizontalFrame(this, 1300, 900, kFixedWidth);
  fGframe5 = new TGGroupFrame(fHor7, "Side view - Proton Module");
  fHor7->AddFrame(fGframe5, new TGLayoutHints(kLHintsTop | kLHintsLeft |
					      kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));
  
  fCanvasSidePM = new TRootEmbeddedCanvas("CanvasSidePM",fGframe5,300,300);
  fGframe5->AddFrame(fCanvasSidePM, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,4,4,12,12));
  fGframe6 = new TGGroupFrame(fHor7, "Side view - Ingrid 3rd Module");
  fHor7->AddFrame(fGframe6, new TGLayoutHints(kLHintsTop | kLHintsRight |
					      kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));
  
  fCanvasSideI = new TRootEmbeddedCanvas("CanvasSideI",fGframe6,300,300);
  fGframe6->AddFrame(fCanvasSideI, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,4,4,12,12));

  AddFrame(fHor7,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));

  fIsCheckedOption1=false;
  fIsCheckedOption2=false;
  fIsCheckedOption3=false;
  fIsCheckedOption4=false;
  fIsCheckedOption5=false;
  fIsCheckedOption6=false;
  
  if(IsData) MyMainFrame::DoFillHistoData();
  else MyMainFrame::DoFillHisto();

  SetCleanup(kDeepCleanup);
  SetWindowName("PM-Ingrid event viewer");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}


///////////////////////////////////////////////////////////////////////////////////////////////////



void MyMainFrame::DoFillHisto()
{
  float pi=acos(-1);

  /*tree=(TTree*) myfile->Get("tree");
  nevt=(int) tree->GetEntries();
  cout<< nevt<<endl;   
  evt = new IngridEventSummary();
  Br=tree->GetBranch("fDefaultReco.");
  Br->SetAddress(&evt);
  tree->SetBranchAddress("fDefaultReco.",&evt);*/
  
  IngridHitSummary * Hit=new IngridHitSummary();

  //ievt = min((int) fNumber->GetNumberEntry()->GetIntNumber(),nevt-1);
  //irec=0;
  //if(ievt>=nevt) fNumber->SetIntNumber(ievt);

  evt->Clear();
  tree->GetEntry(ievt);
  //IngridSimVertexSummary * simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));

  nrec = evt->NPMAnas();
  while (nrec==0) {cout<<"Event "<<ievt<<" : no reconstruction"<<endl;
    evt->Clear();
    tree->GetEntry(++ievt);
    nrec=evt->NPMAnas();}
  //cout<<"Number of reconstructions="<<evt->NPMAnas()<<endl;
 
  if(nrec>1) fNextReco->SetEnabled();

  fNumber->SetIntNumber(ievt);

  /*std::cout << myNumber << std::endl;

  TH1D *a = new TH1D("a","a",100,-5,5);
  a->FillRandom("gaus",myNumber);
  a->Draw();
  TCanvas *CanvasLeft = fCanvasTopPM->GetCanvas();
  CanvasLeft->cd();
  CanvasLeft->Update();*/
  
  //======================================================PDG Code
  std::map< int, std::pair< std::string, bool > > myPDGCode;
  myPDGCode[ 2212] = std::make_pair("p",true);
  myPDGCode[-2212] = std::make_pair("#bar{p}",true);
  myPDGCode[   11] = std::make_pair("e^{-}",true);
  myPDGCode[  -11] = std::make_pair("e^{+}",true);
  myPDGCode[   12] = std::make_pair("#nu_{e}",false);
  myPDGCode[  -12] = std::make_pair("#bar{#nu}_{e}",false);
  myPDGCode[   22] = std::make_pair("#gamma",true);
  myPDGCode[ 2112] = std::make_pair("n",false);
  myPDGCode[-2112] = std::make_pair("#bar{n}",false);
  myPDGCode[  -13] = std::make_pair("#mu^{+}",true);
  myPDGCode[   13] = std::make_pair("#mu^{-}",true);
  myPDGCode[  130] = std::make_pair("K^{0}_{L}",false);
  myPDGCode[  211] = std::make_pair("#pi^{+}",true);
  myPDGCode[ -211] = std::make_pair("#pi^{-}",true);
  myPDGCode[  321] = std::make_pair("K^{+}",true);
  myPDGCode[ -321] = std::make_pair("K^{-}",true);
  myPDGCode[ 3122] = std::make_pair("#Lambda",false);
  myPDGCode[-3122] = std::make_pair("#bar{#Lambda}",false);
  myPDGCode[  310] = std::make_pair("K^{0}_{S}",false);
  myPDGCode[ 3112] = std::make_pair("#Sigma^{-}",true);
  myPDGCode[ 3222] = std::make_pair("#Sigma^{+}",true);
  myPDGCode[ 3212] = std::make_pair("#Sigma^{0}",false);
  myPDGCode[  111] = std::make_pair("#pi^{0}",false);
  myPDGCode[  311] = std::make_pair("K^{0}",false);
  myPDGCode[ -311] = std::make_pair("#bar{K}^{0}",false);
  myPDGCode[   14] = std::make_pair("#nu_{#mu}",false);
  myPDGCode[  -14] = std::make_pair("#bar#nu_{#mu}",false);
  myPDGCode[-3222] = std::make_pair("#bar{#Sigma}^{-}",true);
  myPDGCode[-3212] = std::make_pair("#bar{#Sigma}^{0}",false);
  myPDGCode[-3112] = std::make_pair("#bar{#Sigma}^{+}",true);
  myPDGCode[ 3322] = std::make_pair("#Xi^{0}",false);
  myPDGCode[-3322] = std::make_pair("#bar{#Xi}^{0}",false);
  myPDGCode[ 3312] = std::make_pair("#Xi^{-}",true);
  myPDGCode[-3312] = std::make_pair("#Xi^{+}",true);
  myPDGCode[ 3334] = std::make_pair("#Omega^{-}",true);
  myPDGCode[-3334] = std::make_pair("#bar{#Omega}",true);
  myPDGCode[  -15] = std::make_pair("#tau^{+}",true);
  myPDGCode[   15] = std::make_pair("#tau^{-}",true);
  myPDGCode[   16] = std::make_pair("#nu_{#tau}",false);
  myPDGCode[  -16] = std::make_pair("#bar{#nu}_{#tau}",false);
  myPDGCode[  411] = std::make_pair("D^{+}",true);
  myPDGCode[ -411] = std::make_pair("D^{-}",true);
  myPDGCode[  421] = std::make_pair("D^{0}",false);
  myPDGCode[ -421] = std::make_pair("#bar{D}^{0}",false);
  myPDGCode[  431] = std::make_pair("D_{s}^{+}",true);
  myPDGCode[ -431] = std::make_pair("D_{s}^{-}",true);
  myPDGCode[ 4122] = std::make_pair("#Lambda_{c}^{+}",true);
  myPDGCode[ 4232] = std::make_pair("#Xi_{c}^{+}",true);
  myPDGCode[ 4112] = std::make_pair("#Xi_{c}^{0}",false);
  myPDGCode[ 4322] = std::make_pair("#Xi'_{c}^{+}",true);
  myPDGCode[ 4312] = std::make_pair("#Xi'_{c}^{0}",false);
  myPDGCode[ 4332] = std::make_pair("#Omega_{c}^{0}",false);
  myPDGCode[-4122] = std::make_pair("#bar{#Lambda}_{c}^{-}",true);
  myPDGCode[-4232] = std::make_pair("#bar{#Xi}_{c}^{-}",true);
  myPDGCode[-4132] = std::make_pair("#bar{#Xi}_{c}^{0}",false);
  myPDGCode[-4322] = std::make_pair("#bar{#Xi'}_{c}^{-}",true);
  myPDGCode[-4312] = std::make_pair("#bar{#Xi'}_{c}^{0}",false);
  myPDGCode[-4332] = std::make_pair("#bar{#Omega}_{c}^{0}",false);
  myPDGCode[  221] = std::make_pair("#eta", false);
  
  // ============================================= definition of histograms' background
  
  //--------------------------------------------definition of bins
  double xbins[33];
  for(int i=0;i<8;i++) xbins[i] = 5.0*i;
  for(int i=8;i<24;i++) xbins[i] = 40.+2.5*(i-8);
  for(int i=24;i<=32;i++) xbins[i] = 80.+5.0*(i-24);

  double zbins[72];
  zbins[0] = 0.0;
  zbins[1] = 1.0;
  zbins[2] = 2.3;
  zbins[3] = 3.3;
  for(int i=0;i<34;i++) {
    zbins[2*i+4]=5.0+2.3*i;
    zbins[2*i+5]=6.3+2.3*i;
  }

  //for INGRID
  double zbinsI[33];
  for(int i=0;i<11;i++){
    zbinsI[3*i]=10.7*i;
    zbinsI[3*i+1]=10.7*i+1;
    zbinsI[3*i+2]=10.7*i+2;
  }

  double xbinsI[25];
  for(int i=0;i<25;i++) xbinsI[i]=5.*i;


  //------------------------------------------ definition of Veto plans
  myBoxVeto[0] = new TBox(0,0,1,120);
  myBoxVeto[1] = new TBox(2.3,0,3.3,120);
  myBoxVeto[0]->SetLineColor(kGray+1);
  myBoxVeto[0]->SetFillColor(kGray+1);
  myBoxVeto[0]->SetFillStyle(3004);
  myBoxVeto[1]->SetLineColor(kGray+1);
  myBoxVeto[1]->SetFillColor(kGray+1);
  myBoxVeto[1]->SetFillStyle(3004);
 
 //------------------------------------- definition of scintillator plans X and Y
  
  for(int i=0;i<17;i++) {
    myBoxScinX[i] = new TBox(5.0+4.6*i,0,6.3+4.6*i,120);
    myBoxScinX[i]->SetLineColor(kBlue-9);
    myBoxScinX[i]->SetFillStyle(0);

  }
  
  for(int i=0;i<17;i++) {
    myBoxScinY[i] = new TBox(7.3+4.6*i,0,8.6+4.6*i,120);
    myBoxScinY[i]->SetLineColor(kBlue-9);
    myBoxScinY[i]->SetFillStyle(0);
  }

  //for INGRID
  for(int i=0;i<11;i++) {
    myBoxScinIY[i] = new TBox(10.7*i,0,10.7*i+1,120);
    myBoxScinIY[i]->SetLineColor(kBlue-9);
    myBoxScinIX[i] = new TBox(10.7*i+1,0,10.7*i+2,120);
    myBoxScinIX[i]->SetLineColor(kBlue-9);
    myBoxScinIX[i]->SetFillStyle(0);
    myBoxScinIY[i]->SetFillStyle(0);
  }
  myBoxScinIX[0]->SetFillColor(kGray+1);
  myBoxScinIX[0]->SetFillStyle(3004);
  myBoxScinIY[0]->SetFillColor(kGray+1);
  myBoxScinIY[0]->SetFillStyle(3004);


  //======================================  initialization
  
  /*TFile * _file0=new TFile("PM_MC_Beam1_BirksCorrected_wNoise_ana.root");
    if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;*/


  /*TTree * tree=(TTree*) myfile->Get("tree");
  nevt=(int) tree->GetEntries();
  cout<< nevt<<endl;   
  IngridEventSummary* evt = new IngridEventSummary();
  TBranch * Br=tree->GetBranch("fDefaultReco.");
  Br->SetAddress(&evt);
  tree->SetBranchAddress("fDefaultReco.",&evt);
  IngridHitSummary * Hit=new IngridHitSummary();*/


  //======================================= real vertex, particles tracks, neutrino beam
  /*evt->Clear();
  tree->GetEntry(min(ievt,nevt-1));
  //IngridSimVertexSummary * simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));

  nrec = evt->NPMAnas();
  while (nrec==0) {cout<<ievt<<" Number of reconstructions="<<nrec<<endl;evt->Clear();tree->GetEntry(++ievt);nrec=evt->NPMAnas();}
  //cout<<"Number of reconstructions="<<evt->NPMAnas()<<endl;
 
  if(nrec>1) fNextReco->SetEnabled();


  fNumber->SetIntNumber(ievt);*/

  IngridSimVertexSummary * simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));
  //PMAnaSummary * recon = (PMAnaSummary*) evt->GetPMAna(irec);
  cout<<"\nEvent selected : "<< ievt<<endl;
  cout<<"Number of reconstructions = "<<nrec<< endl;
  cout<<"\nEnergy of neutrino = "<<simver->nuE<<endl; 
  
  
  //-------------------------------------- definition of neutrino beam track and label
  //TLatex* neutx;
  //TLatex* neuty;
  
  //........ what's the label ? ie  incoming neutrino type ? .....................
  if(simver->nutype ==1) {neutx = new TLatex((simver->znu+160)/2, simver->ynu+65, "#nu_{#mu}");neuty = new TLatex((simver->znu+160)/2, simver->xnu+65, "#nu_{#mu}");}
  else if (simver->nutype ==0) {neutx = new TLatex((simver->znu+160)/2, simver->ynu+65, "#bar{#nu}_{#mu}");neuty = new TLatex((simver->znu+160)/2, simver->xnu+65, "#bar{#nu}_{#mu}");}
  else if (simver->nutype ==3) {neutx = new TLatex((simver->znu+160)/2, simver->ynu+65, "#nu_{e}");neuty = new TLatex((simver->znu+160)/2, simver->xnu+65, "#nu_{e}");}
  else if (simver->nutype ==2) {neutx = new TLatex((simver->znu+160)/2, simver->ynu+65, "#bar{#nu}_{e}");neuty = new TLatex((simver->znu+160)/2, simver->xnu+65, "#bar{#nu}_{e}");}
    

  BeamX = new TArrow(0,simver->ynu+60, simver->znu+160-1, simver->ynu+60,15);
  BeamY = new TArrow(0,simver->xnu+60+(simver->znu+160)*tan(3.8/180*pi), simver->znu+160-1, simver->xnu+60,5);
    

  //--------------------------------------- definition of Vertex
  //  WARNING : x and y inverted for all simulated datas (vertex position, particle tracks etc)
  cout<<"Type of Interaction : "<<simver->inttype<<endl;
  
  //So as to set the origin of coordinnates on left botom corner of histograms, we have to make the following translations from the position given by IngridSimVertex (having swapped x and y):
  float Xvertex = simver->ynu + 60;
  float Yvertex = simver->xnu + 60;
  float Zvertex = simver->znu +160;

  nPart = evt->NIngridSimParticles();

  VertexX = new TMarker(Zvertex, Xvertex, 3);
  VertexY = new TMarker(Zvertex, Yvertex, 3);
  
  //-------------------------------------- definition of particles' tracks and labels      
 
  //string mycolorsname [7]={"red", "green","cyan", "orange", "magenta","blue", "brawn"};
  int mycolors [15] = {kRed, kGreen, kCyan, kOrange+1, kMagenta,kBlue,kOrange+4,kBlack,kGreen+3,kMagenta+2};

  float endTrack[3]; // end of simulated track of particle
  float nx, ny, nz;  // impulsion vector (normalized to 1)
  
  double length;
  IngridSimParticleSummary* Part [nPart]; //contains all outgoing particles of the reconstruction
  
  cout<<"\nNumber of particles =  "<<nPart<<endl;
  cout<<"Particles types : ";
    
  int J=0; //for short tracks (see end of loop over particles)

  for(int i=0;i<nPart;i++){
    Part[i] = evt->GetSimParticle(i);
    length = Part[i]->length;

    endTrack[0] = Xvertex + length*Part[i]->dir[1];  //Part->dir normalized to 1
    endTrack[1] = Yvertex + length*Part[i]->dir[0];
    endTrack[2] = Zvertex + length*Part[i]->dir[2];

    nx=Part[i]->dir[1];
    ny=Part[i]->dir[0];
    nz=Part[i]->dir[2];
 
    cout<<Part[i]->pdg<<" ";
    //cout<<"length="<<length<<endl;
    //cout<<nx<<endl;

    //............ particles tracks ....................................
    ParticlesX[i] = new TArrow(Zvertex,Xvertex,endTrack[2],endTrack[0]);
    ParticlesY[i] = new TArrow(Zvertex,Yvertex,endTrack[2],endTrack[1]);
    
    //............. neutral particles are plotted in dashed lines..........
    if(!myPDGCode[Part[i]->pdg].second) {
      ParticlesY[i]->SetLineStyle(7);
      ParticlesX[i]->SetLineStyle(7);
    }

    //-------------------------------------------- labels on XZ histogram
    // let's determine through what side of detector the particle is leaving first
    partTypeX[i]=new TLatex(endTrack[2],endTrack[0],myPDGCode[Part[i]->pdg].first.c_str());
    if(endTrack[2]>85) // end of track after zmax
      {
	if((-75-simver->znu)*nx/nz+Xvertex>125) {
	  partTypeX[i]->SetY(122); partTypeX[i]->SetX(Zvertex+(125-Xvertex)*nz/nx); 
	  //the track leaves the screen first through the top
	}
	else if ((85-Zvertex)*nx/nz+Xvertex<-5) {
	  partTypeX[i]->SetY(-9); partTypeX[i]->SetX(Zvertex+(-5-Xvertex)*nz/nx);
	  //the track leaves the screen first through the bottom
	}
	else {
	  partTypeX[i]->SetX(85); partTypeX[i]->SetY(2+Xvertex+(85-Zvertex)*nx/nz);
	  //the track leaves the screen first through the right side
	}
      }

    else if(endTrack[2]<-5) //end of track before zmin
      {
	if((-5-Zvertex)*nx/nz+Xvertex>125) {
	  partTypeX[i]->SetY(122); partTypeX[i]->SetX(Zvertex+(125-Xvertex)*nz/nx);
	  //the track leaves first through the top
	}
	else if ((-5-Zvertex)*nx/nz+Xvertex<-5) {
	  partTypeX[i]->SetY(-9); partTypeX[i]->SetX(Zvertex+(-5-Xvertex)*nz/nx);
	  //the track leaves first through the bottom
	}
	else {
	  partTypeX[i]->SetX(-8); partTypeX[i]->SetY(2+Xvertex+(-5-Zvertex)*nx/nz);
	  //the track leaves first through the front of detector (left side of screen)
	}
      }
    
    else if(endTrack[0]>125)  // the track leaves through the top with end of track in [zmin, zmax]
      {
	partTypeX[i]->SetY(122);partTypeX[i]->SetX(Zvertex+(125-Xvertex)*nz/nx);
      }
    else if(endTrack[0]<-5) //the track leaves through the bottom with end of track in [zmin, zmax]
      {
	partTypeX[i]->SetY(-9);partTypeX[i]->SetX(Zvertex+(-5-Xvertex)*nz/nx);
      }
    
    //------------------------------------------------- labels on YZ histogram
    partTypeY[i]=new TLatex(endTrack[2],endTrack[1],myPDGCode[Part[i]->pdg].first.c_str());
    if(endTrack[2]>85) // end of track after zmax
      {
	if((85-Zvertex)*ny/nz+Yvertex>125) {partTypeY[i]->SetY(122); partTypeY[i]->SetX(Zvertex+(125-Yvertex)*nz/ny);
	//the track leaves the screen first through the top
	}
	else if ((85-Zvertex)*ny/nz+Yvertex<-5) {
	  partTypeY[i]->SetY(-9); partTypeY[i]->SetX(Zvertex+(-5-Yvertex)*nz/ny);
	//the track leaves the screen first through the bottom
	}
	else {
	  partTypeY[i]->SetX(85); partTypeY[i]->SetY(2+Yvertex+(85-Zvertex)*ny/nz);
	//the track leaves the screen first through the right side
	}
      }
    
    else if(endTrack[2]<-5) //end of track before zmin
      {
	if((-5-Zvertex)*ny/nz+Yvertex>125) {partTypeY[i]->SetY(122); partTypeY[i]->SetX(Zvertex+(125-Yvertex)*nz/ny);
	//the track leaves first through the top
	}
	else if ((-5-Zvertex)*ny/nz+Yvertex<-5) {
	  partTypeY[i]->SetY(-9); partTypeY[i]->SetX(Zvertex+(-5-Yvertex)*nz/ny);
	//the track leaves first through the bottom
	}
	else {
	  partTypeY[i]->SetX(-8); partTypeY[i]->SetY(2+Yvertex+(-5-Zvertex)*ny/nz);
	//the track leaves first through the left side
	}
      }
    
    else if(endTrack[1]>125)  // the track leaves through the top with end of track in [zmin, zmax]
      {
	partTypeY[i]->SetY(122);partTypeY[i]->SetX(Zvertex+(125-Yvertex)*nz/ny);
      }
    else if(endTrack[1]<-5)//the track leaves through the bottom with end of track in [zmin, zmax]
      {
	partTypeY[i]->SetY(-9);partTypeY[i]->SetX(Zvertex+(-5-Yvertex)*nz/ny);
      }
    
    // for particles with short tracks (<5cm), label are drawn next to the vertex
    if(length<5) {
      partTypeX[i]->SetX(Zvertex-4*J-2); 
      partTypeX[i]->SetY(Xvertex-6);
      partTypeY[i]->SetX(Zvertex-4*J-2); 
      partTypeY[i]->SetY(Yvertex-6);

      if (Xvertex<7)  partTypeX[i]->SetY(Xvertex+2);
      if (Yvertex<7)  partTypeY[i]->SetY(Yvertex+2);

      J++;
    }
  }
  
  
  //===== histogram with all hits (even those that don't belong any reconstructed track)========
   
  HISTX=new TH2D("backgroundX","",71,zbins,32,xbins);
  HISTY=new TH2D("backgroundY","",71,zbins,32,xbins);
  HISTXI=new TH2D("backgroundXi","",32,zbinsI,24,xbinsI);
  HISTYI=new TH2D("backgroundYI","",32,zbinsI,24,xbinsI);


  int nHITS = evt->NIngridHits();
  for(int i=0;i<nHITS;i++)
    {
      Hit=evt->GetIngridHit(i);
      if(Hit->mod==16){
	if(Hit->view==0) {HISTX->SetBinContent(HISTX->GetXaxis()->FindBin(Hit->z+0.5),HISTX->GetYaxis()->FindBin(Hit->xy+1.25),(int) Hit->pe);}
	else {HISTY->SetBinContent(HISTY->GetXaxis()->FindBin(Hit->z+0.6),HISTY->GetYaxis()->FindBin(Hit->xy+1.25),(int) Hit->pe);}
      }
      if(Hit->mod==3){
	if(Hit->view==0) {HISTXI->SetBinContent(HISTXI->GetXaxis()->FindBin(Hit->z+0.5),HISTXI->GetYaxis()->FindBin(Hit->xy+1.25),(int) Hit->pe);}
	else {HISTYI->SetBinContent(HISTYI->GetXaxis()->FindBin(Hit->z+0.5),HISTYI->GetYaxis()->FindBin(Hit->xy+1.25),(int) Hit->pe);}
      }


    }
  HISTX->SetFillColor(15); //plotted in grey
  HISTY->SetFillColor(15);
  HISTXI->SetFillColor(15);
  HISTYI->SetFillColor(15);


  //================================================ canvas and  background of histograms
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetDrawBorder(1);
  gStyle->SetFillStyle(1001);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameFillStyle(1001);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(kSolid);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  
  
  
  //------------------------------------------------------ background of histogram XZ
    
  HISTX->GetXaxis()->SetLimits(0,84.8); 
  HISTX->GetXaxis()->SetTitle("z (cm)");
  HISTX->GetXaxis()->SetTitleOffset(1.5);
  HISTX->GetYaxis()->SetTitleOffset(1.5);
  HISTX->GetYaxis()->SetTitle("x (cm)");

  

  for(int i=0;i<nPart;i++){
    ParticlesX[i]->SetLineColor(4);
    ParticlesX[i]->SetLineWidth(2);
  }

  BeamX->SetLineStyle(10); BeamX->SetLineWidth(2); BeamX->SetLineColor(2);
   
  neutx->SetTextSize(0.05);neutx->SetTextColor(2);
  

  //-------------------------------------------------background of histogram YZ
   
  HISTY->GetXaxis()->SetLimits(0,84.8); 
  HISTY->GetXaxis()->SetTitle("z (cm)");
  HISTY->GetXaxis()->SetTitleOffset(1.5);
  HISTY->GetYaxis()->SetTitleOffset(1.5);
  HISTY->GetYaxis()->SetTitle("y (cm)");

  for(int i=0;i<nPart;i++){
    ParticlesY[i]->SetLineColor(4);
    ParticlesY[i]->SetLineWidth(2);
  }
 
  BeamY->SetLineStyle(10); BeamY->SetLineWidth(2);BeamY->SetLineColor(2);
  
  neuty->SetTextSize(0.05);neuty->SetTextColor(2);
  

  //======================================= fill histograms
 

  //if (evt->NPMAnas()!=0)                 // histos defined if there is at least one reconstruction
      //{
      PMAnaSummary * recon = (PMAnaSummary*) evt->GetPMAna(irec);
      
      nTracks=recon->Ntrack;
      cout<<"\nNumber of reconstructed tracks = "<<nTracks<<"\n"<<endl;

      //TH2D *histx [nTracks];
      //TH2D *histy [nTracks];
      double mx=90;double RVertX;
      double my=90;double RVertY;
      
      for(int itrk=0;itrk<nTracks;itrk++){
	
	int nHit = recon->NhitTs(itrk);
	
	histx[itrk] = new TH2D(Form("top view_%d",itrk),"Top view",71,zbins,32,xbins);
	histy[itrk] = new TH2D(Form("side view_%d",itrk),"Side view",71,zbins,32,xbins);
	histxI[itrk] = new TH2D(Form("top viewI_%d",itrk),"Top view",32,zbinsI,24,xbinsI);
	histyI[itrk] = new TH2D(Form("side viewI_%d",itrk),"Side view",32,zbinsI,24,xbinsI);

	
	for (int ihit=0;ihit<nHit;ihit++){
	  Hit=recon->GetIngridHitTrk(ihit,itrk);
	  if (Hit->mod==16){
	    if (Hit->view==0 ){
	      if (Hit->z>6 && Hit->z<mx){ mx=Hit->z; RVertX=Hit->xy; }
	      if(Hit->z>6) histx[itrk]->SetBinContent(histx[itrk]->GetXaxis()->FindBin(Hit->z+0.65),histx[itrk]->GetYaxis()->FindBin(Hit->xy+1.25), Hit->pe);
	    }
	    else {
	      if (Hit->z>6 && Hit->z<my) {my=Hit->z;RVertY=Hit->xy;}
	      if(Hit->z>6) histy[itrk]->SetBinContent(histy[itrk]->GetXaxis()->FindBin(Hit->z+0.65),histy[itrk]->GetYaxis()->FindBin(Hit->xy+1.25), Hit->pe);
	    }
	  }
	  if (Hit->mod==3){
	    if (Hit->view==0 ){
	       histxI[itrk]->SetBinContent(histxI[itrk]->GetXaxis()->FindBin(Hit->z+0.5),histxI[itrk]->GetYaxis()->FindBin(Hit->xy+1.25), Hit->pe);
	    }
	    else {
	       histyI[itrk]->SetBinContent(histyI[itrk]->GetXaxis()->FindBin(Hit->z+0.5),histyI[itrk]->GetYaxis()->FindBin(Hit->xy+1.25), Hit->pe);
	    }
	  }
	}
      }
    
      for(int i=0;i<nTracks;i++) {histx[i]->SetFillColor(mycolors[i]); //histx[i]->Draw("box same");
	histxI[i]->SetFillColor(mycolors[i]);
      }

      for(int i=0;i<nTracks;i++) {histy[i]->SetFillColor(mycolors[i]); //histy[i]->Draw("box same");
	histyI[i]->SetFillColor(mycolors[i]);
      }
      
      RVertexX = new TMarker(mx+0.65, RVertX+2, 30);
      RVertexY = new TMarker(my+0.65, RVertY+2, 30);
      RVertexX->SetMarkerSize(2);
      RVertexY->SetMarkerSize(2);
      //RVertexX->SetMarkerColor(29);
      //RVertexY->SetMarkerColor(29);


      histx[0]->GetXaxis()->SetLimits(0,84.8); 
      histx[0]->GetXaxis()->SetTitle("z (cm)");
      histx[0]->GetYaxis()->SetTitle("x (cm)");
      histx[0]->GetXaxis()->SetTitleOffset(1.5);
      histx[0]->GetYaxis()->SetTitleOffset(1.5);
      
      histy[0]->GetXaxis()->SetLimits(0,84.8); 
      histy[0]->GetXaxis()->SetTitle("z (cm)");
      histy[0]->GetYaxis()->SetTitle("y (cm)");
      histy[0]->GetXaxis()->SetTitleOffset(1.5);
      histy[0]->GetYaxis()->SetTitleOffset(1.5);

      histxI[0]->GetXaxis()->SetLimits(0,84.8); 
      histxI[0]->GetXaxis()->SetTitle("z (cm)");
      histxI[0]->GetYaxis()->SetTitle("x (cm)");
      histxI[0]->GetXaxis()->SetTitleOffset(1.5);
      histxI[0]->GetYaxis()->SetTitleOffset(1.5);
      
      histyI[0]->GetXaxis()->SetLimits(0,84.8); 
      histyI[0]->GetXaxis()->SetTitle("z (cm)");
      histyI[0]->GetYaxis()->SetTitle("y (cm)");
      histyI[0]->GetXaxis()->SetTitleOffset(1.5);
      histyI[0]->GetYaxis()->SetTitleOffset(1.5);






      
  //========================================== and finally labels of particles and vertex

  for(int i=0; i<nPart; i++) {
    partTypeX[i]->SetTextSize(0.05);
    VertexX->SetMarkerSize(3); 
  }
  
  
  for(int i=0; i<nPart; i++) {
    partTypeY[i]->SetTextSize(0.05);
    VertexY->SetMarkerSize(3); 
  }
  
   MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

///////////////////////////////////////////////////////////////////////////////////////////////////


void MyMainFrame::DoFillHistoData()
{
  float pi=acos(-1);
  int mycolors [15] = {kRed, kGreen, kCyan, kOrange+1, kMagenta,kBlue,kOrange+4,kBlack,kGreen+3,kMagenta+2};
  /*tree=(TTree*) myfile->Get("tree");
  nevt=(int) tree->GetEntries();
  cout<< nevt<<endl;   
  evt = new IngridEventSummary();
  Br=tree->GetBranch("fDefaultReco.");
  Br->SetAddress(&evt);
  tree->SetBranchAddress("fDefaultReco.",&evt);*/
  
  IngridHitSummary * Hit=new IngridHitSummary();

  //ievt = min((int) fNumber->GetNumberEntry()->GetIntNumber(),nevt-1);
  //irec=0;
  //if(ievt>=nevt) fNumber->SetIntNumber(ievt);

  evt->Clear();
  tree->GetEntry(ievt);
  //IngridSimVertexSummary * simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));

  nrec = evt->NPMAnas();
  while (nrec==0) {cout<<"Event "<<ievt<<" : no reconstruction"<<endl;
    evt->Clear();
    tree->GetEntry(++ievt);
    nrec=evt->NPMAnas();}
  //cout<<"Number of reconstructions="<<evt->NPMAnas()<<endl;
 
  if(nrec>1) fNextReco->SetEnabled();

  fNumber->SetIntNumber(ievt);

  /*std::cout << myNumber << std::endl;

  TH1D *a = new TH1D("a","a",100,-5,5);
  a->FillRandom("gaus",myNumber);
  a->Draw();
  TCanvas *CanvasLeft = fCanvasTopPM->GetCanvas();
  CanvasLeft->cd();
  CanvasLeft->Update();*/
  
 
  // ============================================= definition of histograms' background
  
  //--------------------------------------------definition of bins
  double xbins[33];
  for(int i=0;i<8;i++) xbins[i] = 5.0*i;
  for(int i=8;i<24;i++) xbins[i] = 40.+2.5*(i-8);
  for(int i=24;i<=32;i++) xbins[i] = 80.+5.0*(i-24);

  double zbins[72];
  zbins[0] = 0.0;
  zbins[1] = 1.0;
  zbins[2] = 2.3;
  zbins[3] = 3.3;
  for(int i=0;i<34;i++) {
    zbins[2*i+4]=5.0+2.3*i;
    zbins[2*i+5]=6.3+2.3*i;
  }

  //for INGRID
  double zbinsI[33];
  for(int i=0;i<11;i++){
    zbinsI[3*i]=10.7*i;
    zbinsI[3*i+1]=10.7*i+1;
    zbinsI[3*i+2]=10.7*i+2;
  }

  double xbinsI[25];
  for(int i=0;i<25;i++) xbinsI[i]=5.*i;


  //------------------------------------------ definition of Veto plans
  myBoxVeto[0] = new TBox(0,0,1,120);
  myBoxVeto[1] = new TBox(2.3,0,3.3,120);
  myBoxVeto[0]->SetLineColor(kGray+1);
  myBoxVeto[0]->SetFillColor(kGray+1);
  myBoxVeto[0]->SetFillStyle(3004);
  myBoxVeto[1]->SetLineColor(kGray+1);
  myBoxVeto[1]->SetFillColor(kGray+1);
  myBoxVeto[1]->SetFillStyle(3004);
 
 //------------------------------------- definition of scintillator plans X and Y
  
  for(int i=0;i<17;i++) {
    myBoxScinX[i] = new TBox(5.0+4.6*i,0,6.3+4.6*i,120);
    myBoxScinX[i]->SetLineColor(kBlue-9);
    myBoxScinX[i]->SetFillStyle(0);

  }
  
  for(int i=0;i<17;i++) {
    myBoxScinY[i] = new TBox(7.3+4.6*i,0,8.6+4.6*i,120);
    myBoxScinY[i]->SetLineColor(kBlue-9);
    myBoxScinY[i]->SetFillStyle(0);
  }

  //for INGRID
  for(int i=0;i<11;i++) {
    myBoxScinIY[i] = new TBox(10.7*i,0,10.7*i+1,120);
    myBoxScinIY[i]->SetLineColor(kBlue-9);
    myBoxScinIX[i] = new TBox(10.7*i+1,0,10.7*i+2,120);
    myBoxScinIX[i]->SetLineColor(kBlue-9);
    myBoxScinIX[i]->SetFillStyle(0);
    myBoxScinIY[i]->SetFillStyle(0);
  }
  myBoxScinIX[0]->SetFillColor(kGray+1);
  myBoxScinIX[0]->SetFillStyle(3004);
  myBoxScinIY[0]->SetFillColor(kGray+1);
  myBoxScinIY[0]->SetFillStyle(3004);


  //======================================  initialization
  
  /*TFile * _file0=new TFile("PM_MC_Beam1_BirksCorrected_wNoise_ana.root");
    if(_file0->IsOpen()) cout << _file0->GetName() <<" is open"<< endl ;*/


  /*TTree * tree=(TTree*) myfile->Get("tree");
  nevt=(int) tree->GetEntries();
  cout<< nevt<<endl;   
  IngridEventSummary* evt = new IngridEventSummary();
  TBranch * Br=tree->GetBranch("fDefaultReco.");
  Br->SetAddress(&evt);
  tree->SetBranchAddress("fDefaultReco.",&evt);
  IngridHitSummary * Hit=new IngridHitSummary();*/


  //======================================= real vertex, particles tracks, neutrino beam
  /*evt->Clear();
  tree->GetEntry(min(ievt,nevt-1));
  //IngridSimVertexSummary * simver = (IngridSimVertexSummary*)(evt->GetSimVertex(0));

  nrec = evt->NPMAnas();
  while (nrec==0) {cout<<ievt<<" Number of reconstructions="<<nrec<<endl;evt->Clear();tree->GetEntry(++ievt);nrec=evt->NPMAnas();}
  //cout<<"Number of reconstructions="<<evt->NPMAnas()<<endl;
 
  if(nrec>1) fNextReco->SetEnabled();


  fNumber->SetIntNumber(ievt);*/

  
  //PMAnaSummary * recon = (PMAnaSummary*) evt->GetPMAna(irec);
  cout<<"\nEvent selected : "<< ievt<<endl;
  cout<<"Number of reconstructions = "<<nrec<< endl;
  
  
  
  //-------------------------------------- definition of neutrino beam track and label
  //TLatex* neutx;
  //TLatex* neuty;
  
  
  /*
  BeamX = new TArrow(0,simver->ynu+60, simver->znu+160-1, simver->ynu+60,15);
  BeamY = new TArrow(0,simver->xnu+60+(simver->znu+160)*tan(3.8/180*pi), simver->znu+160-1, simver->xnu+60,5);
  */

 
  //===== histogram with all hits (even those that don't belong any reconstructed track)========
   
  HISTX=new TH2D("backgroundX","",71,zbins,32,xbins);
  HISTY=new TH2D("backgroundY","",71,zbins,32,xbins);
  HISTXI=new TH2D("backgroundXi","",32,zbinsI,24,xbinsI);
  HISTYI=new TH2D("backgroundYI","",32,zbinsI,24,xbinsI);


  int nHITS = evt->NIngridHits();
  for(int i=0;i<nHITS;i++)
    {
      Hit=evt->GetIngridHit(i);
      double pe;
      if(Hit->pe+Hit->lope>80) pe=Hit->lope;
      else pe=Hit->pe;
      if(Hit->mod==16){
	double xy;
	if(Hit->ch<8) xy=5*Hit->ch;
	else if(Hit->ch<24) xy=2.5*(Hit->ch-8)+40;
	else xy=5*(Hit->ch-24)+80;
	
	if(Hit->view==0) {HISTX->SetBinContent(HISTX->GetXaxis()->FindBin(Hit->z+0.5),HISTX->GetYaxis()->FindBin(xy+1.25),(int) pe);}
	else {HISTY->SetBinContent(HISTY->GetXaxis()->FindBin(Hit->z+0.6),HISTY->GetYaxis()->FindBin(xy+1.25),(int) pe);}
      }
      if(Hit->mod==3){
	if(Hit->view==0) {HISTXI->SetBinContent(HISTXI->GetXaxis()->FindBin(Hit->z+0.5),HISTXI->GetYaxis()->FindBin(Hit->ch*40+1.25),(int) pe);}
	else {HISTYI->SetBinContent(HISTYI->GetXaxis()->FindBin(Hit->z+0.5),HISTYI->GetYaxis()->FindBin(Hit->ch*40+1.25),(int) pe);}
      }


    }
  HISTX->SetFillColor(15); //plotted in grey
  HISTY->SetFillColor(15);
  HISTXI->SetFillColor(15);
  HISTYI->SetFillColor(15);


  //================================================ canvas and  background of histograms
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetDrawBorder(1);
  gStyle->SetFillStyle(1001);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameFillStyle(1001);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(kSolid);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  
  
  
  //------------------------------------------------------ background of histogram XZ
    
  HISTX->GetXaxis()->SetLimits(0,84.8); 
  HISTX->GetXaxis()->SetTitle("z (cm)");
  HISTX->GetXaxis()->SetTitleOffset(1.5);
  HISTX->GetYaxis()->SetTitleOffset(1.5);
  HISTX->GetYaxis()->SetTitle("x (cm)");

  

  BeamX->SetLineStyle(10); BeamX->SetLineWidth(2); BeamX->SetLineColor(2);
   
  

  //-------------------------------------------------background of histogram YZ
   
  HISTY->GetXaxis()->SetLimits(0,84.8); 
  HISTY->GetXaxis()->SetTitle("z (cm)");
  HISTY->GetXaxis()->SetTitleOffset(1.5);
  HISTY->GetYaxis()->SetTitleOffset(1.5);
  HISTY->GetYaxis()->SetTitle("y (cm)");

  for(int i=0;i<nPart;i++){
    ParticlesY[i]->SetLineColor(4);
    ParticlesY[i]->SetLineWidth(2);
  }
 
  BeamY->SetLineStyle(10); BeamY->SetLineWidth(2);BeamY->SetLineColor(2);
  
  neuty->SetTextSize(0.05);neuty->SetTextColor(2);
  

  //======================================= fill histograms
 

  //if (evt->NPMAnas()!=0)                 // histos defined if there is at least one reconstruction
      //{
      PMAnaSummary * recon = (PMAnaSummary*) evt->GetPMAna(irec);
      
      nTracks=recon->Ntrack;
      cout<<"\nNumber of reconstructed tracks = "<<nTracks<<"\n"<<endl;

      //TH2D *histx [nTracks];
      //TH2D *histy [nTracks];
      double mx=90;double RVertX;
      double my=90;double RVertY;
      
      for(int itrk=0;itrk<nTracks;itrk++){
	
	int nHit = recon->NhitTs(itrk);
	
	histx[itrk] = new TH2D(Form("top view_%d",itrk),"Top view",71,zbins,32,xbins);
	histy[itrk] = new TH2D(Form("side view_%d",itrk),"Side view",71,zbins,32,xbins);
	histxI[itrk] = new TH2D(Form("top viewI_%d",itrk),"Top view",32,zbinsI,24,xbinsI);
	histyI[itrk] = new TH2D(Form("side viewI_%d",itrk),"Side view",32,zbinsI,24,xbinsI);

	
	for (int ihit=0;ihit<nHit;ihit++){
	  Hit=recon->GetIngridHitTrk(ihit,itrk);
	  double pe;
	  if(Hit->pe+Hit->lope>80) pe=Hit->lope;

	  if (Hit->mod==16){
	    double xy;
	    if(Hit->ch<8) xy=5*Hit->ch;
	    else if(Hit->ch<24) xy=2.5*(Hit->ch-8)+40;
	    else xy=5*(Hit->ch-24)+80;
	    if (Hit->view==0 ){
	      if (Hit->z>6 && Hit->z<mx){ mx=Hit->z; RVertX=Hit->xy; }
	      if(Hit->z>6) histx[itrk]->SetBinContent(histx[itrk]->GetXaxis()->FindBin(Hit->z+0.65),histx[itrk]->GetYaxis()->FindBin(xy+1.25), pe);
	    }
	    else {
	      if (Hit->z>6 && Hit->z<my) {my=Hit->z;RVertY=Hit->xy;}
	      if(Hit->z>6) histy[itrk]->SetBinContent(histy[itrk]->GetXaxis()->FindBin(Hit->z+0.65),histy[itrk]->GetYaxis()->FindBin(Hit->xy+1.25), pe);
	    }
	  }
	  if (Hit->mod==3){
	    if (Hit->view==0 ){
	       histxI[itrk]->SetBinContent(histxI[itrk]->GetXaxis()->FindBin(Hit->z+0.5),histxI[itrk]->GetYaxis()->FindBin(Hit->ch*5+1.25), pe);
	    }
	    else {
	       histyI[itrk]->SetBinContent(histyI[itrk]->GetXaxis()->FindBin(Hit->z+0.5),histyI[itrk]->GetYaxis()->FindBin(Hit->ch*5+1.25), pe);
	    }
	  }
	}
      }
    
      for(int i=0;i<nTracks;i++) {histx[i]->SetFillColor(mycolors[i]); //histx[i]->Draw("box same");
	histxI[i]->SetFillColor(mycolors[i]);
      }

      for(int i=0;i<nTracks;i++) {histy[i]->SetFillColor(mycolors[i]); //histy[i]->Draw("box same");
	histyI[i]->SetFillColor(mycolors[i]);
      }
      
      RVertexX = new TMarker(mx+0.65, RVertX+2, 30);
      RVertexY = new TMarker(my+0.65, RVertY+2, 30);
      RVertexX->SetMarkerSize(2);
      RVertexY->SetMarkerSize(2);
      //RVertexX->SetMarkerColor(29);
      //RVertexY->SetMarkerColor(29);


      histx[0]->GetXaxis()->SetLimits(0,84.8); 
      histx[0]->GetXaxis()->SetTitle("z (cm)");
      histx[0]->GetYaxis()->SetTitle("x (cm)");
      histx[0]->GetXaxis()->SetTitleOffset(1.5);
      histx[0]->GetYaxis()->SetTitleOffset(1.5);
      
      histy[0]->GetXaxis()->SetLimits(0,84.8); 
      histy[0]->GetXaxis()->SetTitle("z (cm)");
      histy[0]->GetYaxis()->SetTitle("y (cm)");
      histy[0]->GetXaxis()->SetTitleOffset(1.5);
      histy[0]->GetYaxis()->SetTitleOffset(1.5);

      histxI[0]->GetXaxis()->SetLimits(0,84.8); 
      histxI[0]->GetXaxis()->SetTitle("z (cm)");
      histxI[0]->GetYaxis()->SetTitle("x (cm)");
      histxI[0]->GetXaxis()->SetTitleOffset(1.5);
      histxI[0]->GetYaxis()->SetTitleOffset(1.5);
      
      histyI[0]->GetXaxis()->SetLimits(0,84.8); 
      histyI[0]->GetXaxis()->SetTitle("z (cm)");
      histyI[0]->GetYaxis()->SetTitle("y (cm)");
      histyI[0]->GetXaxis()->SetTitleOffset(1.5);
      histyI[0]->GetYaxis()->SetTitleOffset(1.5);


      MyMainFrame::Draw(kFALSE,kFALSE,kFALSE, kFALSE,fIsCheckedOption5,fIsCheckedOption6);
}


//////////////////////////////////////////////////////////////////////////////////////////////////

void MyMainFrame::GoNextEvent()
{
  nrec=0; irec=0;fLabel->SetText(irec);
 
  bool GoUp = kTRUE;
  int iEntry=fNumber->GetNumberEntry()->GetIntNumber();
  if( ievt> iEntry) GoUp=kFALSE;

  //ievt =  fNumber->GetNumberEntry()->GetIntNumber();
  evt->Clear();
  tree->GetEntry(iEntry);
  if(GoUp) cout<<"Going to next event..."<<endl; else cout<<"Going to previous event..."<<endl;

  nrec = evt->NPMAnas();
  while (nrec==0) {cout<<"Event "<<iEntry<<" : no reconstruction"<<endl;evt->Clear();
    if (iEntry==0) GoUp=kTRUE; 
    else if (iEntry==nevt-1) GoUp=kFALSE;
    if (GoUp) {tree->GetEntry(++iEntry);}
    else {tree->GetEntry(--iEntry);}
    nrec=evt->NPMAnas(); }

  if(nrec>1) fNextReco->SetEnabled();
  else fNextReco->SetState(kButtonDisabled, kTRUE);
  
  fNumber->SetIntNumber(iEntry);ievt=iEntry;
  HISTX->Delete();
  HISTY->Delete();
  HISTXI->Delete();
  HISTYI->Delete();

  for(int i=0;i<nTracks;i++) {histx[i]->Delete(); histy[i]->Delete();
  histxI[i]->Delete(); histyI[i]->Delete();}
  if(IsData) MyMainFrame::DoFillHistoData();
  else MyMainFrame::DoFillHisto();  

}


void MyMainFrame::GoNextReco()
{
  std::cout << "Going to next reconstruction..." << std::endl; 
  if(irec==nrec-1) irec=0;
  else irec++;
  cout<<"reco number"<<irec<<endl;
  fLabel->SetText(irec);
  if(IsData) MyMainFrame::DoFillHistoData();
  else MyMainFrame::DoFillHisto();  

}

void MyMainFrame::DrawOption1()
{
  fIsCheckedOption1=!fIsCheckedOption1;
  if (fIsCheckedOption1) cout << "Plotting vertex of MC ..." << std::endl;
  else cout << "Plotting without vertex of MC ..." << std::endl;
  MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

void MyMainFrame::DrawOption2()
{
  fIsCheckedOption2=!fIsCheckedOption2;
  if (fIsCheckedOption2) cout << "Plotting nu beam ..." << std::endl;
  else cout << "Plotting without nu beam ..." << std::endl;
  MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

void MyMainFrame::DrawOption3()
{
  fIsCheckedOption3=!fIsCheckedOption3;
  if (fIsCheckedOption3) cout << "Plotting charged particles' trajectories ..." << std::endl;
  else cout << "Removing charged particles' trajectories ..." << std::endl;
  MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

void MyMainFrame::DrawOption4()
{
  fIsCheckedOption4=!fIsCheckedOption4;
  if (fIsCheckedOption4) cout << "Plotting neutral particles' trajectories ..." << std::endl;
  else cout << "Removing neutral particles' trajectories ..." << std::endl;
  MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

void MyMainFrame::DrawOption5()
{
  fIsCheckedOption5=!fIsCheckedOption5;
  if (fIsCheckedOption5) cout << "Plotting extra hits ..." << std::endl;
  else cout << "Removing extra hits ..." << std::endl;
  MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

void MyMainFrame::DrawOption6()
{
    fIsCheckedOption6=!fIsCheckedOption6;
  if (fIsCheckedOption5) std::cout << "Plotting reconstructed vertex ..." << std::endl; 
  else cout << "Removing reconstructed vertex ..." << std::endl;
  MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5,fIsCheckedOption6);
}

void MyMainFrame::DrawOption7()
{
  std::cout << "I am using veto cut" << std::endl; 
  fIsCheckedOption5=!fIsCheckedOption5;
  //MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5);
}

void MyMainFrame::DrawOption8()
{
  std::cout << "I am using FV cut" << std::endl; 
  fIsCheckedOption8=!fIsCheckedOption8;
  //MyMainFrame::Draw(fIsCheckedOption1,fIsCheckedOption2,fIsCheckedOption3,fIsCheckedOption4,fIsCheckedOption5);
}

void MyMainFrame::Draw(bool opt1,bool opt2, bool opt3, bool opt4, bool opt5, bool opt6) {

  fCanvasTopPM->GetCanvas()->cd();

  if(opt5) {HISTX->Draw("box");}
  else histx[0]->Draw("box");

  if(opt4) {for(int i=0;i<nPart;i++) {
      if(ParticlesX[i]->GetLineStyle()==7) ParticlesX[i]->Draw("same");}}
  if(opt3) {for(int i=0;i<nPart;i++) {
      if(ParticlesX[i]->GetLineStyle()!=7) ParticlesX[i]->Draw("same");}}
  
  for(int i=0;i<nTracks;i++) {histx[i]->Draw("box same");}
  myBoxVeto[0]->Draw("same");
  myBoxVeto[1]->Draw("same");
  for(int i=0;i<17;i++) {myBoxScinX[i]->Draw("same"); }

  if(opt4) {for(int i=0;i<nPart;i++) {
      if(ParticlesX[i]->GetLineStyle()==7) partTypeX[i]->Draw("same");}}
  if(opt3){for(int i=0;i<nPart;i++) {
      if(ParticlesX[i]->GetLineStyle()!=7) partTypeX[i]->Draw("same");}}
  if(opt6){RVertexX->Draw("same");}
  if(opt1) {VertexX->Draw("same");}
  if(opt2) {BeamX->Draw("same"); neutx->Draw("same");}

  //if(opt5) HISTX->Draw("text same");

  fCanvasTopPM->GetCanvas()->Update();

  //-------------------------------------------
  fCanvasSidePM->GetCanvas()->cd();

  if(opt5) {HISTY->Draw("box");}
  else histy[0]->Draw("box");

  if(opt3) {for(int i=0;i<nPart;i++) {
      if(ParticlesY[i]->GetLineStyle()!=7) ParticlesY[i]->Draw("same");}}
  if(opt4) {for(int i=0;i<nPart;i++) {
      if(ParticlesY[i]->GetLineStyle()==7) ParticlesY[i]->Draw("same");}}

  for(int i=0;i<nTracks;i++) {histy[i]->Draw("box same");}
  myBoxVeto[0]->Draw("same");
  myBoxVeto[1]->Draw("same");
  for(int i=0;i<17;i++) {myBoxScinY[i]->Draw("same"); }

  if(opt4) {for(int i=0;i<nPart;i++) {
      if(ParticlesY[i]->GetLineStyle()==7) partTypeY[i]->Draw("same");}}
  if(opt3){for(int i=0;i<nPart;i++) {
      if(ParticlesY[i]->GetLineStyle()!=7) partTypeY[i]->Draw("same");}}
  if(opt6) {RVertexY->Draw("same");}
  if(opt1) {VertexY->Draw("same");}
  if(opt2) {BeamY->Draw("same"); neuty->Draw("same");}
 
  //if(opt5) HISTY->Draw("text same");

 fCanvasSidePM->GetCanvas()->Update();

 //------------------------------------------------------
 fCanvasTopI->GetCanvas()->cd();
 HISTXI->Draw("box");
 for(int i=0;i<nTracks;i++) {histxI[i]->Draw("box same");}
 for(int i=0;i<11;i++) {myBoxScinIX[i]->Draw("same"); }
 //HISTXI->Draw("text same");
 fCanvasTopI->GetCanvas()->Update();


 //------------------------------------------------------
 fCanvasSideI->GetCanvas()->cd();
 HISTYI->Draw("box");
 for(int i=0;i<nTracks;i++) {histyI[i]->Draw("box same");}
 for(int i=0;i<11;i++) {myBoxScinIY[i]->Draw("same"); }
 //HISTYI->Draw("text same");
 fCanvasSideI->GetCanvas()->Update();



}

void MyMainFrame::OpenFile()
{
  TGFileInfo file_info_;
  const char *filetypes[] = {"ROOT files", "*.root", 0, 0};
  file_info_.fFileTypes = filetypes;
  file_info_.fIniDir = StrDup(".");
  TGFileDialog *dlg = new TGFileDialog(gClient->GetDefaultRoot(),
				       gClient->GetDefaultRoot(),kFDOpen, &file_info_);
  if( file_info_.fFilename ) {
    std::cout << "'" << file_info_.fFilename << "' selected." << std::endl;
  
    myfile = new TFile(file_info_.fFilename);
    tree=(TTree*) myfile->Get("tree");
    nevt=(int) tree->GetEntries();
    fNumber->SetLimitValues(0,nevt-1);
    cout<<"Number of events = "<<nevt<<endl;   
    evt = new IngridEventSummary();
    Br=tree->GetBranch("fDefaultReco.");
    Br->SetAddress(&evt);
    tree->SetBranchAddress("fDefaultReco.",&evt);
    ievt=0;irec=0;
    tree->GetEntry(0);
    if(evt->NIngridSimVertexes()==0) IsData=kTRUE; else IsData=kFALSE;
    if(IsData){
      MyMainFrame::DoFillHistoData();
      fButOption1->SetEnabled(kFALSE);
      fButOption2->SetEnabled(kFALSE);
      fButOption3->SetEnabled(kFALSE);
      fButOption4->SetEnabled(kFALSE);
      cout<<"Data file"<<endl;
    }
    else {
      fButOption1->SetEnabled();
      fButOption2->SetEnabled();
      fButOption3->SetEnabled();
      fButOption4->SetEnabled();
      MyMainFrame::DoFillHisto();
      cout<<"Monte-Carlo file"<<endl;
    }
  }
}

void MyMainFrame::SaveAs()
{
  TGFileInfo file_info_;
  const char *filetypes[] = {"PNG", "*.png", 0, 0};
  file_info_.fFileTypes = filetypes;
  file_info_.fIniDir = StrDup(".");
  TGFileDialog *dlg = new TGFileDialog(gClient->GetDefaultRoot(),
				       gClient->GetDefaultRoot(),kFDSave, &file_info_);
  if( file_info_.fFilename )
    std::cout << "Saving as...'" << file_info_.fFilename << "' selected." << std::endl;
  
  UInt_t nMainFrames = 0;
  TClass* clGMainFrame = TClass::GetClass("TGMainFrame");
  TGWindow* win = 0;
  TIter iWin(gClient->GetListOfWindows());
  while ((win = (TGWindow*)iWin())) {
    const TObject* winGetParent = win->GetParent();
    Bool_t winIsMapped = kFALSE;
    if (winGetParent == gClient->GetDefaultRoot())
      winIsMapped = win->IsMapped();
    if (winIsMapped && win->InheritsFrom(clGMainFrame)) {
      win->MapRaised();
      TString outfile = file_info_.fFilename;
      outfile += ".png";
      std::cout << outfile.Data() << std::endl;	  
      TImage *img = TImage::Create();
      win->RaiseWindow();
      img->FromWindow(win->GetId());
      img->WriteImage(outfile.Data());
      delete img;
    }
  }
}


int main(int argc, char **argv)
{
  


  TApplication *theApp = new TApplication("MyMainFrame", &argc, argv);
  MyMainFrame *myFrame = new MyMainFrame(gClient->GetRoot(), 50, 50);
  // terminate the application when closing the window
  myFrame->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  theApp->Run();
  return 0;

}



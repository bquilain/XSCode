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

#include "INGRIDEVENTSUMMARY.h"

class MyMainFrame : public TGMainFrame
{
 private:
  TGCompositeFrame    *fHor1;
  TGCompositeFrame    *fHor2;
  TGCompositeFrame    *fHor7;
  TGTextButton        *fExit;
  TGTextButton        *fExit2;
  TGTextButton        *fExit3;
  TGGroupFrame        *fGframe;
  TGGroupFrame        *fGframe2;
  TGGroupFrame        *fGframe3;
  TGGroupFrame        *fGframe4;
  TGGroupFrame        *fGframe5;
  TGGroupFrame        *fGframe6;


  TGNumberEntry       *fNumber;
  TGLabel             *fLabel;
  TGLabel             *fLabel2;
  TRootEmbeddedCanvas *fCanvasTopPM;
  TRootEmbeddedCanvas *fCanvasSidePM;
  TRootEmbeddedCanvas *fCanvasTopI;
  TRootEmbeddedCanvas *fCanvasSideI;
  TH1D                *myGauss;

  TH1D                *a;
  TGTextButton        *fNextReco;
  TGCompositeFrame    *fHor3;
  TGGroupFrame        *fGframeMCOptions;
  TGGroupFrame        *fGframeROptions;
  TGCompositeFrame    *fVer1;
  TGCompositeFrame    *fVer2;
  TGCompositeFrame    *fVer3;
  TGCompositeFrame    *fVer4;
  TGCheckButton       *fButOption1;
  TGCheckButton       *fButOption2;
  TGCheckButton       *fButOption3;
  TGCheckButton       *fButOption4;
  TGCheckButton       *fButOption5;
  TGCheckButton       *fButOption6;
  TGCheckButton       *fButOption7;
  TGCheckButton       *fButOption8;
  Bool_t              fIsCheckedOption1;
  Bool_t              fIsCheckedOption2;
  Bool_t              fIsCheckedOption3;
  Bool_t              fIsCheckedOption4;
  Bool_t              fIsCheckedOption5;
  Bool_t              fIsCheckedOption6;
  Bool_t              fIsCheckedOption7;
  Bool_t              fIsCheckedOption8;
  int                 ievt;
  int                 nevt;
  int                 irec;
  int                 nrec;
  TLatex              *neutx;
  TLatex              *neuty;
  TArrow              *BeamX;
  TArrow              *BeamY;
  TMarker             *VertexX;
  TMarker             *VertexY;
  TMarker             *RVertexX;
  TMarker             *RVertexY;
  TArrow              *ParticlesX [25];
  TArrow              *ParticlesY [25];
  TLatex              *partTypeX [25]; //labels of particles on XZ histogram
  TLatex              *partTypeY [25];
  TH2D                *HISTX;
  TH2D                *HISTY;
  TH2D                *HISTXI;
  TH2D                *HISTYI;

  TH2D                *histx [25];
  TH2D                *histy [25];
  TH2D                *histxI [25];
  TH2D                *histyI [25];

  int                 nPart;
  int                 nTracks;
  TBox                *myBoxVeto[2] ;
  TBox                *myBoxScinX[17] ;
  TBox                *myBoxScinY[17] ;
  TBox                *myBoxScinIX[11] ;
  TBox                *myBoxScinIY[11] ;

  TFile               *myfile;
  TGCompositeFrame    *fHor4;
  TGCompositeFrame    *fHor6;
  TGCompositeFrame    *fHor5;
  //  TGFileDialog        *inFile;
  TGTextButton        *fOpenFile;
  TGTextButton        *fSaveAs;
  TTree               *tree;
  TBranch             *Br;
  IngridEventSummary  *evt;
  Bool_t               IsData;


 public:
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~MyMainFrame() {}
  void DoFillHisto();
  void DoFillHistoData();
  void GoNextEvent();
  //void GoPreviousEvent();
  void GoNextReco();
  void DrawOption1();
  void DrawOption2();
  void DrawOption3();
  void DrawOption4();
  void DrawOption5();
  void DrawOption6();
  void DrawOption7();
  void DrawOption8();
  void Draw(bool opt1=kFALSE, bool opt2=kFALSE,bool opt3=kFALSE,bool opt4=kFALSE,bool opt5=kFALSE, bool opt6=kFALSE);
  void OpenFile();
  void SaveAs();
	   
};

void DoFillHisto();
void DoFillHistoData();
void GoNextEvent();
void GoNextReco();
void DrawOption1();
void DrawOption2();
void DrawOption3();
void DrawOption4();
void DrawOption5();
void DrawOption6();
void DrawOption7();
void DrawOption8();
void Draw(bool opt1=kFALSE, bool opt2=kFALSE, bool opt3=kFALSE, bool opt4=kFALSE, bool opt5=kFALSE, bool opt6=kFALSE);
void OpenFile();
void SaveAs();
//void GoPreviousEvent();

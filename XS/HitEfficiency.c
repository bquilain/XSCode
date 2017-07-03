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
#include <TF1.h>
#include <TH1.h>
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
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "INGRID_Dimension.cc"
#include "Reconstruction.cc"
#include "Xsec.cc"
#include "setup.h"

//#define DEBUG

int main(int argc, char **argv){

  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  int c=-1;
  bool MC=false;
  int NFiles=30;
  int IFiles=0;
  int IRuns=14510;
  int FRuns=14510;
  char Name[256];

  bool PM=true;

  while ((c = getopt(argc, argv, "wmi:f:r:t:")) != -1) {
    switch(c){
    case 'm':
      MC=true;
      break;
    case 'w':
      PM=false;
      break;
    case 'i':
      IFiles=atoi(optarg);
      break;
    case 'f':
      NFiles=atoi(optarg);
      break;
    case 'r':
      IRuns=atoi(optarg);
      break;
    case 't':
      FRuns=atoi(optarg);
      break;
    }
  }
  
  INGRID_Dimension * IngDim = new INGRID_Dimension();
  Reconstruction * Rec = new Reconstruction(PM);
  Xsec * xs = new Xsec(PM);
  xs->Xsec::Initialize();

  char DetName[2];sprintf(DetName,(PM?"PM":"WM"));
  char * cMCOUT = getenv((PM?"MCOUTPUTSTORAGE":"MCOUTPUTSTORAGE_WM"));
  char * cDATAOUT = getenv((PM?"DATAOUTPUTSTORAGE":"DATAOUTPUTSTORAGE_WM"));

  // for PM
  TH1D * hEffX_Ing = new TH1D("hEffX_Ing","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffY_Ing = new TH1D("hEffY_Ing","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffXNorm_Ing = new TH1D("hEffXNorm_Ing","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffYNorm_Ing = new TH1D("hEffYNorm_Ing","YZ distribution",NBinsRecAngle,BinningRecAngle);

  TH1D * hEffX_Sci = new TH1D("hEffX_Sci","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffY_Sci = new TH1D("hEffY_Sci","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffXNorm_Sci = new TH1D("hEffXNorm_Sci","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffYNorm_Sci = new TH1D("hEffYNorm_Sci","YZ distribution",NBinsRecAngle,BinningRecAngle);

  // for WM
  TH1D * hEffX_Plan = new TH1D("hEffX_Plan","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffY_Plan = new TH1D("hEffY_Plan","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffXNorm_Plan = new TH1D("hEffXNorm_Plan","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffYNorm_Plan = new TH1D("hEffYNorm_Plan","YZ distribution",NBinsRecAngle,BinningRecAngle);

  TH1D * hEffX_Grid = new TH1D("hEffX_Grid","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffY_Grid = new TH1D("hEffY_Grid","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffXNorm_Grid = new TH1D("hEffXNorm_Grid","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffYNorm_Grid = new TH1D("hEffYNorm_Grid","YZ distribution",NBinsRecAngle,BinningRecAngle);

  // for test
  TH1D* touched=new TH1D("touched","touched",60,-3,3);
  TH1D* nottouched=new TH1D("nottouched","not touched",60,-3,3);
  TH2D* gridhit=new TH2D("gridhit","expected vs true grid hit",20,0,20,20,0,20);

  TH1D * hHitsX = new TH1D("hHitsX","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hHitsY = new TH1D("hHitsY","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hHitsXNorm = new TH1D("hHitsXNorm","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hHitsYNorm = new TH1D("hHitsYNorm","YZ distribution",NBinsRecAngle,BinningRecAngle);

  double Icyc=4;
  double Ncyc=12;
  double Nmod=17;

  //TFile * _file0;
  char FileName[256];
  TTree * tree;
  IngridEventSummary* evt = new IngridEventSummary();
  IngridHitSummary * Hit = new IngridHitSummary(); 
  PMAnaSummary * recon = new PMAnaSummary();

  const int NPln=(PM?NPlnPM:NPlnWM);
  // following offsets to have the same coordinate system for both detectors
  //  (in particular, provides positive value for the hit position)  -- ML 2017/06/28
  float offsetZ=(PM?0.:40.95);
  float offsetXY=(PM?0.:60.);
  

  for(int R=IRuns;R<=FRuns;R++){
    for(int f=IFiles;f<=NFiles;f++){
      // ---- choice of the input files ----
      //if(MC && PM &&(f>=415 && f<=419)) continue; ML 2017/06/23
      if(MC) {
	sprintf(FileName,"%s/%sMC_Wall_Run1_%d_wNoise_ana.root",cMCOUT,DetName,f);
      }
      else {
	if(PM) sprintf(FileName,"/export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root",R,f);
      	else sprintf(FileName,"%s/ingrid_%08d_%04d_anadev.root",cDATAOUT,R,f);    
      }
      
      TFile * _file0=new TFile(FileName);
      
      if(_file0->IsOpen()) cout << _file0->GetName() << " is open"<<endl ;
      else continue;


      // ---- loading the content from the chosen file ----
      tree=(TTree*) _file0->Get("tree");
      double Nhits;
      int nevt=(int) tree->GetEntries();      
      tree->SetBranchAddress("fDefaultReco.",&evt);
      cout<<nevt<<endl;


      // ---- loop over the events ----
      for(int ievt=0;ievt<nevt;ievt++){
	if((ievt%1000)==0) cout<<ievt<<endl;
	evt->Clear();//vire l'evt précédent
	tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
	
	int NRec=evt->NPMAnas();
	for(int irec=0;irec<NRec;irec++){
	  recon=(PMAnaSummary*) evt->GetPMAna(irec);
	  int Mod=(PM?16:15);//recon->hitmod;
	  int nTracks=recon->Ntrack;

	  bool VSelectionOV=false;bool VSelectionFV=true;
	  Rec->Reconstruction::GetSelectionPM(&VSelectionFV,&VSelectionOV,recon,MC);

	  //bool sandMuon=VSelectionOV;// old definition
	  // new definition is following -- ML 2017/06/23
       	  //bool sandMuon=recon->Ntrack==1 && (min(recon->startxpln[0],recon->startypln[0])<(PM?1:2)) && (max(recon->endxpln[0],recon->endypln[0])>(PM?15:20));
	  bool sandMuon=recon->Ntrack==1 && VSelectionOV/* && (PM || recon->NhitTs(0)>5)*/;


	  if(sandMuon) {
	    for(int itrk=0;itrk<nTracks;itrk++){//loop on track  
	      double AngleX=TMath::Abs((recon->thetax)[itrk]);
	      double AngleY=TMath::Abs((recon->thetay)[itrk]);
	      double weight=1;

	      int MinplnX=NPln;
	      int MaxplnX=0;
	      int MinplnY=NPln;
	      int MaxplnY=0;                                                                                                     
	      int NHitsX=0;
	      int NHitsY=0;
	      int NHitPlnX[NPln];
	      int NHitPlnY[NPln];
	      bool IsHitGridX[NPln][loli_gridchnum];
	      bool IsHitGridY[NPln][loli_gridchnum];
	      for(int ipln=0;ipln<NPln;ipln++){
		NHitPlnX[ipln]=0;
		NHitPlnY[ipln]=0;
		for(int ich=0;ich<loli_gridchnum;ich++){
		  IsHitGridX[ipln][ich]=false;
		  IsHitGridY[ipln][ich]=false;
		}
	      }
#ifdef DEBUG
	      for(int ipln=0;ipln<NPln;ipln++){
		cout<<"ipln="<<ipln<<"  , Active Plane="<<NHitPlnX[ipln]<<endl;
	      }
#endif
	
	      // 1st step : get the fitted track
	      TF1 * fEffX = new TF1("fEffX","pol1",0,120);                                                                                                 
	      TF1 * fEffY = new TF1("fEffY","pol1",0,120);        
	   
	      TH1D * HistEffX = new TH1D("HistEffX","XZ distribution",120,0,120);
	      TH1D * HistEffY = new TH1D("HistEffY","YZ distribution",120,0,120);
	      TH1D * NormX = new TH1D("NormX","Norm",120,0,120);
	      TH1D * NormY = new TH1D("NormY","Norm",120,0,120);
	      HistEffX->Sumw2();
	      NormX->Sumw2();
	      HistEffY->Sumw2();
	      NormY->Sumw2();

	      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
		Hit=recon->GetIngridHitTrk(ihit,itrk);
		if(Hit->mod!=Mod) continue; // look at the target module

		// ML 2017/06/26 to get the recon pln -- unchanged for INGRID & PM
		int reconpln,reconch;
		IngDim->get_reconplnch_loli(Hit->mod,Hit->view,Hit->pln,Hit->ch,0,&reconpln,&reconch);
		
		if(Hit->view==0){
		  //cout<<"hit pos="<<Hit->xy<<", pln="<<Hit->pln<<endl;
		  HistEffX->Fill(Hit->z+offsetZ,(Hit->xy+offsetXY)*weight);
		  NormX->Fill(Hit->z+offsetZ,weight);

		  if(reconpln>MaxplnX) MaxplnX=reconpln;
		  if(reconpln<MinplnX) MinplnX=reconpln;
		  NHitsX++;
		  if(NHitPlnX[reconpln]==0) NHitPlnX[reconpln]=1;
		  if(reconpln%3 !=0) IsHitGridX[reconpln][reconch]=true;
#ifdef DEBUG
		  cout<<"NHits="<<NHitsX<<", Hit pos z="<<Hit->z<<", pos y="<<Hit->xy<<endl;
#endif
		}
		else{
		  HistEffY->Fill(Hit->z+offsetZ,(Hit->xy+offsetXY)*weight);
		  NormY->Fill(Hit->z+offsetZ,weight); 
		  if(reconpln>MaxplnY) MaxplnY=reconpln;
		  if(reconpln<MinplnY) MinplnY=reconpln;
		  NHitsY++;
		  if(NHitPlnY[reconpln]==0) NHitPlnY[reconpln]=1;
		  if(reconpln%3 !=1) IsHitGridY[reconpln][reconch]=true;
#ifdef DEBUG
		  cout<<"NHits="<<NHitsY<<", Hit pos z="<<Hit->z<<", pos x="<<Hit->xy<<endl;
#endif
		}  
	      }
	  	   
	      HistEffX->Divide(NormX);
	      HistEffY->Divide(NormY);
	      HistEffX->Fit("fEffX","RQ0");   
	      HistEffY->Fit("fEffY","RQ0");  
#ifdef DEBUG
	      for(int i=0; i<120;i++) if(HistEffX->GetBinContent(i)>0)cout<<"bin "<<i<<" (cm) : "<<HistEffX->GetBinContent(i)<<endl;
	      for(int i=0; i<120;i++) if(HistEffY->GetBinContent(i)>0)cout<<"bin "<<i<<" (cm) : "<<HistEffY->GetBinContent(i)<<endl;
	      cout<<"angleX="<<recon->thetax[itrk]<<", Param ="<<atan(fEffX->GetParameter(1))*180/TMath::Pi()<<endl;
	      cout<<"angleY="<<recon->thetay[itrk]<<", Param ="<<atan(fEffY->GetParameter(1))*180/TMath::Pi()<<endl;
#endif


	      // 2nd step : get the rate of planes actually hit 
	      //--for the PM
	      bool IsIngrid; 
	      double TotalIngridX=0;
	      double IngridX=0;
	      double TotalSciBarX=0;
	      double SciBarX=0;

	      double TotalIngridY=0;
	      double IngridY=0;
	      double TotalSciBarY=0;
	      double SciBarY=0;

	      //--for the WM
	      bool IsGrid; 
	      double TotalPlanX=0;
	      double PlanX=0;
	      double TotalGridX=0;
	      double GridX=0;

	      double TotalPlanY=0;
	      double PlanY=0;
	      double TotalGridY=0;
	      double GridY=0;


	      for(int reconpln=MinplnX;reconpln<=MaxplnX;reconpln++){
		double xy,z;
		IngDim->get_posi_lolirecon(Mod,0,reconpln,0,0,&xy,&z);
		double Pos=fEffX->Eval(z+offsetZ); // xy position where we look at a hit

		if(PM){
		  IsIngrid=(reconpln==0 || (Pos<40 || Pos>=80));
		  if(IsIngrid){
		    TotalIngridX++;
		    if(NHitPlnX[reconpln]==1) IngridX++; // the table is filled so that the possibe values are only 0 or 1
		  }
		  else{
		    TotalSciBarX++;
		    if(NHitPlnX[reconpln]==1)SciBarX++;
		  }
		}
		else {// WM
		  IsGrid=(reconpln%3 != 0);
		  if(IsGrid){
		    bool IsTouchable=false;
		    double DeltaPos=0.5*loli_scinti_width*fabs(fEffX->GetParameter(1));
		    //cout<<"slope="<<fEffX->GetParameter(1)<<" angle="<<recon->thetax[itrk]<<" Pos="<<Pos<<" Dpos="<<DeltaPos<<endl;
		    int reconch=0;
		    for(reconch=0;reconch<loli_gridchnum;reconch++){
		      IngDim->get_posi_lolirecon(15,0,reconpln,reconch,0,&xy,&z);
		      IsTouchable=( fabs(Pos-xy-offsetXY)<DeltaPos );
		      //cout<<"pln="<<reconpln<<" ch="<<reconch<<" y="<<xy<<" touchable? "<<IsTouchable<<endl;
		      if(IsTouchable) break;
		    }
		    if(IsTouchable){
		      if(IsHitGridX[reconpln][reconch]) gridhit->Fill(reconch,reconch);
		      else {
			for(int reconch2=0;reconch2<loli_gridchnum;reconch2++)
			  if(IsHitGridX[reconpln][reconch2]) {
			    gridhit->Fill(reconch,reconch2);
			    if(abs(reconch-reconch2)>2) {
			      cout<<ievt<<" view0 pln="<<reconpln<<" z="<<z+offsetZ<<" expch="<<reconch<<" truech="<<reconch2<<" angle="<<AngleY<<" slope="<<fEffY->GetParameter(1)<<" intcpt="<<fEffY->GetParameter(0)<<endl;
			      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
				Hit=recon->GetIngridHitTrk(ihit,itrk);
				if(Hit->mod!=15 || Hit->view!=0) continue;
				int _reconpln,_reconch;
				IngDim->get_reconplnch_loli(15,0,Hit->pln,Hit->ch,0,&_reconpln,&_reconch);
				if(reconpln==_reconpln && reconch2==_reconch) cout<<"ch="<<Hit->ch<<" xy="<<Hit->xy+offsetXY<<" z="<<Hit->z+offsetZ<<endl;
			      }
			    }
			    break;
			  }
		      }
		    }

		    if(IsTouchable){
		      TotalGridX++;
		      if(NHitPlnX[reconpln]==1){
			if(IsTouchable)GridX++;
			//cout<<"touched"<<endl;
			touched->Fill((Pos-xy-offsetXY)/DeltaPos);
		      }
		    }
		    else {
		      nottouched->Fill((Pos-xy-offsetXY)/DeltaPos);
		      // cout<<"not touched "<<NHitPlnX[reconpln]<<endl;
		    }
		  }
		  else{
		    TotalPlanX++;
		    if(NHitPlnX[reconpln]==1)PlanX++;
		  }
		}
	      }

	      for(int reconpln=MinplnY;reconpln<=MaxplnY;reconpln++){
		double xy,z;
		IngDim->get_posi_lolirecon(Mod,1,reconpln,0,0,&xy,&z);
		double Pos=fEffY->Eval(z+offsetZ); // xy position where we look at a hit
		  
		if(PM){
		  IsIngrid=(reconpln==0 || (Pos<40 || Pos>=80));
		  if(IsIngrid){
		    TotalIngridY++;
		    if(NHitPlnY[reconpln]==1) IngridY++;
		  }
		  else{
		    TotalSciBarY++;
		    if(NHitPlnY[reconpln]==1) SciBarY++;
		  }
		}
		else{// WM
		  IsGrid=(reconpln%3 != 1);
		  if(IsGrid){
		    bool IsTouchable=false;
		    double DeltaPos=0.5*loli_scinti_width*fabs(fEffY->GetParameter(1));
		    int reconch=0;
		    for(reconch=0;reconch<loli_gridchnum;reconch++){
		      IngDim->get_posi_lolirecon(15,1,reconpln,reconch,0,&xy,&z);
		      IsTouchable=( fabs(Pos-xy-offsetXY)<DeltaPos );
		      if(IsTouchable) break;
		    }
		    
		    if(IsTouchable){
		      if(IsHitGridY[reconpln][reconch])gridhit->Fill(reconch,reconch);
		      else {
			for(int reconch2=0;reconch2<loli_gridchnum;reconch2++)
			  if(IsHitGridY[reconpln][reconch2]) {
			    gridhit->Fill(reconch,reconch2);
			    if(abs(reconch-reconch2)>2) {
			      cout<<ievt<<" view1 pln="<<reconpln<<" z="<<z+offsetZ<<" expch="<<reconch<<" truech="<<reconch2<<" angle="<<AngleY<<" slope="<<fEffY->GetParameter(1)<<" intcpt="<<fEffY->GetParameter(0)<<endl;
			      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
				Hit=recon->GetIngridHitTrk(ihit,itrk);
				if(Hit->mod!=15 || Hit->view!=1) continue;
				int _reconpln,_reconch;
				IngDim->get_reconplnch_loli(15,1,Hit->pln,Hit->ch,0,&_reconpln,&_reconch);
				if(reconpln==_reconpln && reconch2==_reconch) cout<<"ch="<<Hit->ch<<" xy="<<Hit->xy+offsetXY<<" z="<<Hit->z+offsetZ<<endl;
			      }
			    }
			    break;
			  }
		      }
		    }
		    
		    if(IsTouchable){
		      TotalGridY++;
		      if(NHitPlnY[reconpln]==1){
			GridY++;
			touched->Fill((Pos-xy-offsetXY)/DeltaPos);
		      }
		      else nottouched->Fill((Pos-xy-offsetXY)/DeltaPos);
		    }
		  }
		  else{
		    TotalPlanY++;
		    if(NHitPlnY[reconpln]==1)PlanY++;
		  }
		}
	      }
	      //cout<<"Number of SciBar scintillator actually hit="<<SciBarX<<", And should be="<<TotalSciBarX<<endl;

	      if(TotalIngridX!=0) IngridX/=TotalIngridX;
	      if(TotalSciBarX!=0) SciBarX/=TotalSciBarX;
	      if(TotalIngridY!=0) IngridY/=TotalIngridY;
	      if(TotalSciBarY!=0) SciBarY/=TotalSciBarY;
	      if(TotalPlanX!=0) PlanX/=TotalPlanX;
	      if(TotalGridX!=0) GridX/=TotalGridX;
	      if(TotalPlanY!=0) PlanY/=TotalPlanY;
	      if(TotalGridY!=0) GridY/=TotalGridY;

	      //double InEffX=(double)1.- ((double)PlnX/(double)HitPlnX);
	      //double InEffY=(double)1.- ((double)PlnY/(double)HitPlnY);
	      /*
	      //int Minpln=MinplnX;
	      //int Maxpln=MaxplnX;
	      //if(MinplnX!=MinplnY || MaxplnX!=MaxplnY) cout<<"difference"<<endl;
	      int TotEffX=0;
	      int TotEffY=0;

	      double chfixed=0;//to change if we want also posxy
	      double * posz=new double();
	      double * posxy=new double();

	      for(int ipln=MinplnX;ipln<=MaxplnX;ipln++){
	      
	      HistEffX->Delete();
	      NormX->Delete();
	      HistEffX = new TH1D("HistEffX","XZ distribution",120,0,120);
	      NormX = new TH1D("NormX","Norm",120,0,120);
	      int EfficiencyX=0;
	     
	      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
	      Hit=recon->GetIngridHitTrk(ihit,itrk);
	      if(Hit->view!=0 || Hit->pln==ipln) continue; //remove hits of this plane or of different view                                                   
	      HistEffX->Fill(Hit->z,Hit->xy);                                                                                                                  
	      NormX->Fill(Hit->z);                                                                                                                                 
	      }
	     
	      HistEffX->Divide(NormX);
	      //now fit with a line
	      TF1 * EffX = new TF1("EffX","[0]+x*[1]",0,130);
	      HistEffX->Fit("EffX","RQI");
	      //cout<<"Angle X calculated="<<TMath::ATan(EffX->GetParameter(1))*180/TMath::Pi()<<endl;
	      //then check if there is a hit around
	      bool doneX=IngDim->INGRID_Dimension::get_posXY(Mod,0,ipln,chfixed,posxy,posz);
	      double ZscintiX=*posz;
	      double Xscinti=EffX->Eval(ZscintiX);
	      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
	      Hit=recon->GetIngridHitTrk(ihit,itrk);
	      if(Hit->view!=0 || Hit->pln!=ipln) continue;//now looking only hits in this plane
	      if(TMath::Abs(Hit->xy-Xscinti)<7.5) {
	      EfficiencyX++;
	      break;//not to double count if there are 2 hitsin near scinti
	      }
	      }
	      TotEffX+=EfficiencyX;
	      }
	   
	   
	   
	      for(int ipln=MinplnY;ipln<=MaxplnY;ipln++){
	     
	      HistEffY->Delete();
	      NormY->Delete();
	      HistEffY = new TH1D("HistEffY","YZ distribution",120,0,120);
	      NormY = new TH1D("NormY","Norm",120,0,120);
	      int EfficiencyY=0;
	     
	     
	      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
	      Hit=recon->GetIngridHitTrk(ihit,itrk);
	      if(Hit->view!=1 || Hit->pln==ipln) continue; //remove hits of this plane or of different view                                                   
	      HistEffY->Fill(Hit->z,Hit->xy);
	      NormY->Fill(Hit->z);
	      }
	     
	      HistEffY->Divide(NormY);	      
	      TF1 * EffY = new TF1("EffY","[0]+x*[1]",0,120);
	      HistEffY->Fit("EffY","RQI");
	      //cout<<"Angle Y="<<TMath::ATan(EffY->GetParameter(1))*180/TMath::Pi()<<endl;
	      //then check if there is a hit around
	      bool doneY=IngDim->INGRID_Dimension::get_posXY(Mod,1,ipln,chfixed,posxy,posz);
	      double ZscintiY=*posz;
	      double Yscinti=EffY->Eval(ZscintiY);
	      for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
	      Hit=recon->GetIngridHitTrk(ihit,itrk);
	      if(Hit->view!=1 || Hit->pln!=ipln) continue;//now looking only hits in this plane
	      if(TMath::Abs(Hit->xy-Yscinti)<7.5) {
	      EfficiencyY++;
	      break;//not to double count if there are 2 hits in near scinti
	      }
	      }
	      //	      cout<<"Plane="<<ipln<<"   , EfficiencyX="<<EfficiencyX<<"   , EfficiencyY="<<EfficiencyY<<endl;
	      TotEffY+=EfficiencyY;
	      }
	      */	   
	      //cout<<"Amount of hit planes, on X view="<<HitPlnX<<"  , on Y view="<<HitPlnY<<endl;	    
	      //cout<<"Amount of Planes Hits, in X="<<TotEffX<<" , in Y="<<TotEffY<<endl;   
	      //double TotEfficiencyX=(double) (1-(double)TotEffX/(double)HitPlnX);//inneffiency
	      //double TotEfficiencyY=(double) (1-(double)TotEffY/(double)HitPlnY);
	      //cout<<setprecision(4)<<"InEffiency Via Fit="<<TotEfficiencyX<<"   "<<TotEfficiencyY<<endl;
	      //cout<<"InEffiency="<<InEffX<<"   "<<InEffY<<endl<<endl;
	      //cout<<(double) (TotEffX/HitPlnX)<<"   "<<(double) (TotEffY/HitPlnY)<<endl;
	      /*if(TotEfficiencyX>0){
		cout<<"X loos efficiency"<<endl;
		for(int ihit=0;ihit<Vec.size();ihit++){
		if(Vec[ihit].view==0) cout<<Vec[ihit].x<<" ,"<<Vec[ihit].y<<" ,"<<Vec[ihit].z<<endl;
		}
		}
		if(TotEfficiencyY>0){
		cout<<"Y loos efficiency"<<endl;
		for(int ihit=0;ihit<Vec.size();ihit++){
		if(Vec[ihit].view==1) cout<<Vec[ihit].x<<" ,"<<Vec[ihit].y<<" ,"<<Vec[ihit].z<<endl;
		} 
		}*/


	      // 3rd step : get the histogram of this rate wrt 2D angle
	      // PM
	      if(TotalIngridX!=0){
		hEffX_Ing->Fill(AngleX,IngridX);
		hEffXNorm_Ing->Fill(AngleX,1.);
	      }
	      if(TotalSciBarX!=0){
		hEffX_Sci->Fill(AngleX,SciBarX);
		hEffXNorm_Sci->Fill(AngleX,1.);
	      }
	      if(TotalIngridY!=0){
		hEffY_Ing->Fill(AngleY,IngridY);
		hEffYNorm_Ing->Fill(AngleY,1.);
	      }
	      if(TotalSciBarY!=0){
		hEffY_Sci->Fill(AngleY,SciBarY);
		hEffYNorm_Sci->Fill(AngleY,1.);
	      }

	      // WM
	      if(TotalPlanX!=0){
		hEffX_Plan->Fill(AngleX,PlanX);
		hEffXNorm_Plan->Fill(AngleX,1.);
	      }
	      if(TotalGridX!=0){
		hEffX_Grid->Fill(AngleX,GridX);
		hEffXNorm_Grid->Fill(AngleX,1.);
	      }
	      if(TotalPlanY!=0){
		hEffY_Plan->Fill(AngleY,PlanY);
		hEffYNorm_Plan->Fill(AngleY,1.);
	      }
	      if(TotalGridY!=0){
		hEffY_Grid->Fill(AngleY,GridY);
		hEffYNorm_Grid->Fill(AngleY,1.);
	      }


	      hHitsX->Fill(AngleX,NHitsX);
	      hHitsXNorm->Fill(AngleX);
	      hHitsY->Fill(AngleY,NHitsY);
	      hHitsYNorm->Fill(AngleY);


	      delete fEffX;
	      delete fEffY;
	      delete HistEffX;
	      delete HistEffY;
	      delete NormX;
	      delete NormY;

	    }
	  }
	}
      }
      _file0->Close();
      delete _file0;
      //_file0->Delete();
    }
  }

  if(MC) sprintf(Name,"%s/XS/files/MCHitEfficiency_test_%s.root",cINSTALLREPOSITORY,DetName);
  else sprintf(Name,"%s/XS/files/DataHitEfficiency_test_%s.root",cINSTALLREPOSITORY,DetName);
  cout<<"writing in "<<Name<<endl;
  TFile * wfile = new TFile(Name,"recreate");

  if(PM){
    hEffX_Ing->Sumw2();
    hEffXNorm_Ing->Sumw2();
    hEffY_Ing->Sumw2();
    hEffYNorm_Ing->Sumw2();
    hEffX_Ing->Add(hEffY_Ing);
    hEffXNorm_Ing->Add(hEffYNorm_Ing);
    hEffX_Ing->Divide(hEffXNorm_Ing);
    hEffX_Ing->Write("Efficiency_PMIng");
  
    hEffX_Sci->Sumw2();
    hEffXNorm_Sci->Sumw2();
    hEffY_Sci->Sumw2();
    hEffYNorm_Sci->Sumw2();
    hEffX_Sci->Add(hEffY_Sci);
    hEffXNorm_Sci->Add(hEffYNorm_Sci);
    hEffX_Sci->Divide(hEffXNorm_Sci);
    hEffX_Sci->Write("Efficiency_PMSci");
  }
  else {
    hEffX_Plan->Sumw2();
    hEffXNorm_Plan->Sumw2();
    hEffY_Plan->Sumw2();
    hEffYNorm_Plan->Sumw2();
    hEffX_Plan->Add(hEffY_Plan);
    hEffXNorm_Plan->Add(hEffYNorm_Plan);
    hEffX_Plan->Divide(hEffXNorm_Plan);
    hEffX_Plan->Write("Efficiency_WMPlan");
  
    hEffX_Grid->Sumw2();
    hEffXNorm_Grid->Sumw2();
    hEffY_Grid->Sumw2();
    hEffYNorm_Grid->Sumw2();
    hEffX_Grid->Add(hEffY_Grid);
    hEffXNorm_Grid->Add(hEffYNorm_Grid);
    hEffX_Grid->Divide(hEffXNorm_Grid);
    hEffX_Grid->Write("Efficiency_WMGrid");

    // for test
    touched->Write();
    nottouched->SetLineColor(kRed);
    nottouched->Write();
    for(int i=1;i<=loli_gridchnum;i++){
      cout<<"expected channel "<<i<<": "<<gridhit->GetBinContent(i,i)*100./gridhit->Integral(i,i,1,20) <<" % of true channel "<<i<<endl;
    }
    gridhit->GetXaxis()->SetTitle("expected grid hit");
    gridhit->GetYaxis()->SetTitle("true grid hit(s)");
    gridhit->Write();
  }
  return 0;
}

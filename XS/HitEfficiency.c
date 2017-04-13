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
INGRID_Dimension * IngDim = new INGRID_Dimension();
Reconstruction * Rec = new Reconstruction();
Xsec * xs = new Xsec();


int main(int argc, char **argv){

  int c=-1;
  bool MC=false;
  int NFiles=30;
  int IFiles=0;
  int IRuns=14510;
  int FRuns=14510;
  char Name[256];

  while ((c = getopt(argc, argv, "mi:f:r:t:")) != -1) {
    switch(c){
    case 'm':
      MC=true;
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

  xs->Xsec::Initialize();


  TH1D * hEffX_Ing = new TH1D("hEffX_Ing","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffY_Ing = new TH1D("hEffY_Ing","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffXNorm_Ing = new TH1D("hEffXNorm_Ing","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffYNorm_Ing = new TH1D("hEffYNorm_Ing","YZ distribution",NBinsRecAngle,BinningRecAngle);

  TH1D * hEffX_Sci = new TH1D("hEffX_Sci","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffY_Sci = new TH1D("hEffY_Sci","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffXNorm_Sci = new TH1D("hEffXNorm_Sci","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hEffYNorm_Sci = new TH1D("hEffYNorm_Sci","YZ distribution",NBinsRecAngle,BinningRecAngle);

  TH1D * hHitsX = new TH1D("hHitsX","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hHitsY = new TH1D("hHitsY","YZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hHitsXNorm = new TH1D("hHitsXNorm","XZ distribution",NBinsRecAngle,BinningRecAngle);
  TH1D * hHitsYNorm = new TH1D("hHitsYNorm","YZ distribution",NBinsRecAngle,BinningRecAngle);

  double Icyc=4;
  double Ncyc=12;
  double Nmod=17;

  //TFile * _file0;
  char FileMC[256],FileData[256];
  TTree * tree;
  
  for(int R=IRuns;R<=FRuns;R++){
    for(int f=IFiles;f<=NFiles;f++){
      if(MC && (f>=415 && f<=419)) continue;
      if(MC) sprintf(FileMC,"/export/scraid2/data/bquilain/MCfiles/PMMC_Sand_Run1_%d_wNoise_anareduced.root",f);
      else sprintf(FileMC,"/export/scraid2/data/bquilain/DataNew/ingrid_%08d_%04d_pmmergedKSPManabsd_woXTalk.root",R,f);
      
       TFile * _file0=new TFile(FileMC);
      
      if(_file0->IsOpen()) cout << _file0->GetName() << "is open"<<endl ;
      else continue;
      
      tree=(TTree*) _file0->Get("tree");
      double Nhits;
      IngridHitSummary * Hit = new IngridHitSummary(); 
      int nevt=(int) tree->GetEntries();
      IngridEventSummary* evt = new IngridEventSummary();
      TBranch * Br=tree->GetBranch("fDefaultReco.");
      Br->SetAddress(&evt);
      tree->SetBranchAddress("fDefaultReco.",&evt);
      cout<<nevt<<endl;
      for(int ievt=0;ievt<nevt;ievt++){
	if((ievt%1000)==0) cout<<ievt<<endl;
	evt->Clear();//vire l'evt précédent
	tree->GetEntry(ievt);//charge l'evt grace au link avec la branche
	//for(int imod=0;imod<Nmod;imod++){
	//for(int icyc=Icyc;icyc<Ncyc;icyc++){
	//	Nhits=evt->NIngridModHits(imod,icyc);
	//	for(int ihit=0;ihit<Nhits;ihit++){

    int NRec=evt->NPMAnas();
    for(int irec=0;irec<NRec;irec++){
       PMAnaSummary * recon=(PMAnaSummary*) evt->GetPMAna(irec);
       int Mod=16;//recon->hitmod;
       int nTracks=recon->Ntrack;

       bool VSelectionOV=false;bool VSelectionFV=true;
       Rec->Reconstruction::GetSelectionPM(&VSelectionFV,&VSelectionOV,recon,MC);
       
       if(VSelectionOV) {
	 for(int itrk=0;itrk<nTracks;itrk++){//loop on track  
	   double AngleX=TMath::Abs((recon->thetax)[itrk]);
	   double AngleY=TMath::Abs((recon->thetay)[itrk]);
	   double weight=1;

	   int MinplnX=NPlnPM;
	   int MaxplnX=0;
	   int MinplnY=NPlnPM;
	   int MaxplnY=0;                                                                                                     
	   int NHitsX=0;
	   int NHitsY=0;
	   int NHitPlnX[NPlnPM];
	   int NHitPlnY[NPlnPM];
	   for(int ipln=0;ipln<NPlnPM;ipln++){
	     NHitPlnX[ipln]=0;
	     NHitPlnY[ipln]=0;
	   }

	   /* for(int ipln=0;ipln<NPlnPM;ipln++){
             cout<<"ipln="<<ipln<<"  , Active Plane="<<NHitPlnX[ipln]<<endl;
	     }*/
	   
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
             if(Hit->mod!=16) continue;

	     if(Hit->view==0){
	       //cout<<"hit pos="<<Hit->xy<<", pln="<<Hit->pln<<endl;
	       HistEffX->Fill(Hit->pln*4.6,Hit->xy*weight);
	       NormX->Fill(Hit->pln*4.6,weight);
	       if(Hit->pln>MaxplnX) MaxplnX=Hit->pln;
               if(Hit->pln<MinplnX) MinplnX=Hit->pln;
               NHitsX++;
               if(NHitPlnX[Hit->pln]==0) NHitPlnX[Hit->pln]=1;
#ifdef DEBUG
	       cout<<"NHits="<<NHitsX<<", Hit pos z="<<Hit->pln*4.6<<", pos y="<<Hit->xy<<endl;
#endif
	     }
	     else{
	       HistEffY->Fill(Hit->pln*4.6,Hit->xy*weight);
               NormY->Fill(Hit->pln*4.6,weight); 
	       if(Hit->pln>MaxplnY) MaxplnY=Hit->pln;
               if(Hit->pln<MinplnY) MinplnY=Hit->pln;
               NHitsY++;
               if(NHitPlnY[Hit->pln]==0) NHitPlnY[Hit->pln]=1;
             }  
	   }
	  	   
	   HistEffX->Divide(NormX);
	   HistEffY->Divide(NormY);
	   HistEffX->Fit("fEffX","RQ0");   
	   HistEffY->Fit("fEffY","RQ0");  
	   //cout<<"angleX="<<AngleX<<", Param ="<<atan(fEffX->GetParameter(1))*180/TMath::Pi()<<endl;
	   //cout<<"angleY="<<AngleY<<", Param ="<<atan(fEffY->GetParameter(1))*180/TMath::Pi()<<endl;
	   
	   bool IsIngrid; 
	   double TotalIngridX=0;
	   double IngridX=0;
	   double TotalSciBarX=0;
	   double SciBarX=0;

	   double TotalIngridY=0;
	   double IngridY=0;
	   double TotalSciBarY=0;
	   double SciBarY=0;

	   for(int ipln=MinplnX;ipln<=MaxplnX;ipln++){
	     double Pos=fEffX->Eval(ipln*4.6);
	     //cout<<"pln="<<ipln<<endl;
	     //cout<<"pos="<<Pos<<endl;
	     if(ipln==0 || (Pos<40 || Pos>=80)) IsIngrid=true;
	     else IsIngrid=false;
	     if(IsIngrid){
	       //cout<<"yes"<<endl;
	       TotalIngridX++;
	       if(NHitPlnX[ipln]==1) IngridX++;
	     }
	     else{
	       //cout<<"Pos="<<Pos<<endl;
	       TotalSciBarX++;
	       if(NHitPlnX[ipln]==1){SciBarX++;}
	     }
	   }

	   for(int ipln=MinplnY;ipln<=MaxplnY;ipln++){
	     double Pos=fEffY->Eval(ipln*4.6);
	     if(ipln==0 || (Pos<40 || Pos>=80)) IsIngrid=true;
	     else IsIngrid=false;
	     if(IsIngrid){
	       TotalIngridY++;
	       if(NHitPlnY[ipln]==1) IngridY++;
	     }
	     else{
	       TotalSciBarY++;
	       if(NHitPlnY[ipln]==1) SciBarY++;
	     }
	   }
	   //cout<<"Number of SciBar scintillator actually hit="<<SciBarX<<", And should be="<<TotalSciBarX<<endl;

	   if(TotalIngridX!=0) IngridX/=TotalIngridX;
	   if(TotalSciBarX!=0) SciBarX/=TotalSciBarX;
	   if(TotalIngridY!=0) IngridY/=TotalIngridY;
	   if(TotalSciBarY!=0) SciBarY/=TotalSciBarY;
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
	   /*
	   if(TotalIngridY!=0){
	     hEffY_Ing->Fill(AngleY,IngridY);
	     hEffYNorm_Ing->Fill(AngleY);
	   }
	   if(TotalSciBarY!=0){
	     hEffY_Sci->Fill(AngleY,SciBarY);
	     hEffYNorm_Sci->Fill(AngleY);
	   }
	   */
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

  if(MC) sprintf(Name,"/home/bquilain/CC0pi_XS/XS/files/MCHitEfficiency.root");
  else sprintf(Name,"/home/bquilain/CC0pi_XS/XS/files/DataHitEfficiency.root");
  cout<<"writing in "<<Name<<endl;
  TFile * wfile = new TFile(Name,"recreate");
  
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

  return 0;
}

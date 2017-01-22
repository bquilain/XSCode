

#ifndef Reconstruction_cc
#define Reconstruction_cc
/************CLibs********/
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
/********ROOT Libs*******/
#include <TF1.h>
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
#include <TMatrixD.h>
#include <TGraphErrors.h>
/**********INGRID Libs*************/
#include "TApplication.h"
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridSimVertexSummary.h"
#include "IngridSimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "IngridBasicReconSummary.h"
#include "INGRID_Dimension.cc"
/*********Hit Libs****************/
#include "Hit.h"
#include "Reconstruction.h"
#include "setup.h"
#include "PMrecon.hxx"
//#define DEBUG2
//
INGRID_Dimension * IngDimRec = new INGRID_Dimension();


vector <Hit3D> Reconstruction::ApplyPEError(vector <Hit3D> Vec, double angle){
    int BinAngle= (int) (angle/3);
  char FileName[32]; char HistName[32];
  sprintf(FileName,"Input/Landau/Landau[%d].root",BinAngle);
  sprintf(HistName,"Landau[%d]",BinAngle);
  TFile * _file = new TFile(FileName);
  TF1 * Landau = (TF1*) _file->Get(HistName);
  //Landau->SetNpx(1000);
  //cout<<"Maximum="<<Landau->GetMaximumX(5,200)<<", value at max="<<Landau->Eval(Landau->GetMaximumX(5,200))<<endl;
  //cout<<"function value at 20 p.e="<<Landau->Eval(20)<<endl;
  double Mean=Landau->GetMaximumX(5,200);
  //Landau->SaveAs("test.root");
  //cout<<"Mean="<<Mean<<endl;

  for(int ihit=0;ihit<Vec.size();ihit++){
    double Rand=Landau->GetRandom();
    //cout<<"PE="<<Vec[ihit].pe<<", Random="<<Rand<<", value of angle="<<angle<<", PE after="<<Vec[ihit].pe+(Rand-Mean)*(Vec[ihit].pe/Mean)<<endl;
    Vec[ihit].pe+=(Rand-Mean)*(Vec[ihit].pe/Mean);
  }

  return Vec;
}

/**Determine the track sample, depending on where it stops**/
/**0, 1, 2, 3, 4, 5 respectively corresponds to stopped in PM, out of PM without being geometrically able to reach INGRID, out of PM in the direction of INGRID w/o INGRID track, with and INGRID track but stopped in INGRID, with an INGRID track side escaping, with an INGRID track through going**/
int Reconstruction::SelectTrackSample(bool pm_stop, bool Geom, bool has_ingrid, bool ingrid_stop, int ing_last_pln){
  int TrackSample;

  if(has_ingrid && !ingrid_stop && ing_last_pln>=9) TrackSample=5;
  else if(has_ingrid && !ingrid_stop) TrackSample=4;
  else if(has_ingrid && ingrid_stop) TrackSample=3;
  else if(!has_ingrid && !pm_stop && Geom) TrackSample=2;
  else if(!has_ingrid && !pm_stop && !Geom) TrackSample=1;
  else if(!has_ingrid && pm_stop && !Geom) TrackSample=0;

  return TrackSample;
}
/***********************************************************/

/**Search for hits in INGRID that are aligned with the PM track, but are not enough to be reconstructed as an INGRID track**/
//To check
vector <Hit3D> Reconstruction::SearchIngridHit(vector <Hit3D> Vec, vector <Hit3D> VecAll, double thetaX, double thetaY, double TrackSample){

  if(TrackSample<=2){
    double ZX,ZY;
    double XMean=0;double YMean=0;
    double NClusterX=0;double NClusterY=0;
    IngridHitSummary * Hit = new IngridHitSummary();
    
    for(int ihit=0;ihit<Vec.size();ihit++){
      if(Vec[ihit].pln!=17) continue;
      else{
	if(Vec[ihit].view==0) {
	  XMean+=Vec[ihit].x;
	  ZX=Vec[ihit].z;
	  NClusterX++;
	}
	else{
	  YMean+=Vec[ihit].y;
	  ZY=Vec[ihit].z;
	  NClusterY++;
	}
      }
    }
    if(NClusterX!=0) XMean/=NClusterX;
    if(NClusterY!=0) YMean/=NClusterY;
    
    double Xpln0, Ypln0, Xpln1, Xpln2, Ypln1, Ypln2;
    for(int ihit=0;ihit<VecAll.size();ihit++){
      if(VecAll[ihit].mod==3 && VecAll[ihit].pln==0){
	Xpln0=TMath::Tan(thetaX)*(VecAll[ihit].z-ZX)+XMean;
	Ypln0=TMath::Tan(thetaY)*(VecAll[ihit].z-ZY)+YMean;
	if(Xpln0>120 || Xpln0<0) continue;
	if(Ypln0>120 || Ypln0<0) continue;
	bool Used=false;
	if(TMath::Abs(VecAll[ihit].x-Xpln0)<7.5 && TMath::Abs(VecAll[ihit].y-Ypln0)<7.5){
	  for(int ihit2=0;ihit2<Vec.size();ihit2++){
	    if(Vec[ihit2]==VecAll[ihit]) Used=true;
	  }
	  if(!Used){
	    Vec.push_back(VecAll[ihit]);
	    
	    for(int ihit3=0;ihit3<VecAll.size();ihit3++){
	      if(VecAll[ihit3].mod==VecAll[ihit].mod && VecAll[ihit3].pln==VecAll[ihit].pln && VecAll[ihit3].view==VecAll[ihit].view && TMath::Abs(VecAll[ihit].ch-VecAll[ihit3].ch)==1){
		bool Used2=false;
		for(int ihit2=0;ihit2<Vec.size();ihit2++){
		  if(Vec[ihit2]==VecAll[ihit3]) Used2=true;
		}
		if(!Used2) Vec.push_back(VecAll[ihit3]);
	      }
	    }
	  }
	}
      }
      
      
      else if(VecAll[ihit].mod==3 && VecAll[ihit].pln==1){
	Xpln1=TMath::Tan(thetaX)*(VecAll[ihit].z-ZX)+XMean;
	Ypln1=TMath::Tan(thetaY)*(VecAll[ihit].z-ZY)+YMean;
	if(Xpln1>120 || Xpln1<0) continue;
	if(Ypln1>120 || Ypln1<0) continue;
	bool Used=false;
	if(TMath::Abs(VecAll[ihit].x-Xpln1)<7.5 && TMath::Abs(VecAll[ihit].y-Ypln1)<7.5){
	  for(int ihit2=0;ihit2<Vec.size();ihit2++){
	    if(Vec[ihit2]==VecAll[ihit]) Used=true;
	  }
	  if(!Used){
	    Vec.push_back(VecAll[ihit]);
	    
	    for(int ihit3=0;ihit3<VecAll.size();ihit3++){
	      if(VecAll[ihit3].mod==VecAll[ihit].mod && VecAll[ihit3].pln==VecAll[ihit].pln && VecAll[ihit3].view==VecAll[ihit].view && TMath::Abs(VecAll[ihit].ch-VecAll[ihit3].ch)==1){
		bool Used2=false;
		for(int ihit2=0;ihit2<Vec.size();ihit2++){
		  if(Vec[ihit2]==VecAll[ihit3]) Used2=true;
		}
		if(!Used2) Vec.push_back(VecAll[ihit3]);
	      }
	    }
	  }
	}
      }
      //cout<<"Here in plane1. Angle are Thetax="<<thetaX*180/TMath::Pi()<<", Thetay="<<thetaY*180/TMath::Pi()<<", Last Hit in PM=("<<ZX<<","<<XMean<<","<<YMean<<") , Hit in plane=("<<VecAll[ihit].z<<","<<Xpln1<<","<<Ypln1<<") and real x is="<<VecAll[ihit].x<<","<<VecAll[ihit].y<<endl;
      
      else if(VecAll[ihit].mod==3 && VecAll[ihit].pln==2){
	Xpln2=TMath::Tan(thetaX)*(VecAll[ihit].z-ZX)+XMean;
	Ypln2=TMath::Tan(thetaY)*(VecAll[ihit].z-ZY)+YMean;
	if(Xpln2>120 || Xpln2<0) continue;
	if(Ypln2>120 || Ypln2<0) continue;
	bool Used=false;
	if(TMath::Abs(VecAll[ihit].x-Xpln2)<7.5 && TMath::Abs(VecAll[ihit].y-Ypln2)<7.5){
	  for(int ihit2=0;ihit2<Vec.size();ihit2++){
	    if(Vec[ihit2]==VecAll[ihit]) Used=true;
	  }
	  
	  if(!Used){
	  Vec.push_back(VecAll[ihit]);
	  
	  for(int ihit3=0;ihit3<VecAll.size();ihit3++){
	    if(VecAll[ihit3].mod==VecAll[ihit].mod && VecAll[ihit3].pln==VecAll[ihit].pln && VecAll[ihit3].view==VecAll[ihit].view && TMath::Abs(VecAll[ihit].ch-VecAll[ihit3].ch)==1){
	      bool Used2=false;
	      for(int ihit2=0;ihit2<Vec.size();ihit2++){
		if(Vec[ihit2]==VecAll[ihit3]) Used2=true;
	      }
	      if(!Used2) Vec.push_back(VecAll[ihit3]);
	    }
	    //cout<<"Here in plane2. Angle are Thetax="<<thetaX*180/TMath::Pi()<<", Thetay="<<thetaY*180/TMath::Pi()<<", Last Hit in PM=("<<ZX<<","<<XMean<<","<<YMean<<") , Hit in plane=("<<VecAll[ihit].z<<","<<Xpln2<<","<<Ypln2<<") and real x is="<<VecAll[ihit].x<<","<<VecAll[ihit].y<<endl;
	  }
	}
      }
    }
  }
  }



  return Vec;
}
/*************************************************************************************************************/

vector <Hit3D> Reconstruction::CountSharedHits(vector <Hit3D> Vec, vector< vector<Hit3D> > VecDouble, int Trk){
  vector <Hit3D> VecD(Vec);
  for(int itrk=0;itrk<VecDouble.size();itrk++){
    if(itrk==Trk) continue;
    else{
      for(int ihit=0;ihit<Vec.size();ihit++){
	//VecD[ihit]=Vec[ihit];
	for(int ihit2=0;ihit2<VecDouble[itrk].size();ihit2++){
	  if(Vec[ihit]==VecDouble[itrk][ihit2]){
#ifdef DEBUG2
	    cout<<"Original Track we compare="<<Trk<<" , Track with which we compare="<<itrk<<endl;
	    cout<<" Position=("<<Vec[ihit].x<<","<<Vec[ihit].y<<","<<Vec[ihit].z<<")   , NRJ="<<Vec[ihit].pe<<" , used="<<Vec[ihit].used<<endl;
	    cout<<" Posihittihiton=("<<VecDouble[itrk][ihit2].x<<","<<VecDouble[itrk][ihit2].y<<","<<VecDouble[itrk][ihit2].z<<")   , NRJ="<<VecDouble[itrk][ihit2].pe<<" , used="<<VecDouble[itrk][ihit2].used<<endl;
#endif
	  VecD[ihit].used++;
	  }
	  else continue;
	}
      }
    }
  }
  return VecD;
}


vector <HitTemp> Reconstruction::EraseDoubleHits(IngridBasicReconSummary * recon, int itrk, vector <HitTemp> HitV){
  //cout<<"hello"<<endl;
  HitV.clear();
  IngridHitSummary * Hit = new IngridHitSummary();
  HitTemp Coord;
  int ndouble=0;
  for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
    Hit=recon->GetIngridHitTrk(ihit,itrk);
    //Coord.x=Hit->view;
    //Coord.y=Hit->xy;
    //Coord.z=Hit->z;
    Coord.view=Hit->view;
    Coord.ch=Hit->ch;
    Coord.pln=Hit->pln; 
    Coord.trk=itrk;
    Coord.hit=ihit;
    Coord.mod=Hit->mod;
    HitV.push_back(Coord);

    for(int ihit2=0;ihit2<HitV.size()-1;ihit2++){
      if(HitV[ihit2]==Coord) {
	ndouble++;
	HitV.pop_back();
      }
    }
  }
  return HitV;
  
}

vector <HitTemp> Reconstruction::EraseDoubleHitsPM(PMAnaSummary * recon, int itrk, vector <HitTemp> HitV){
  HitTemp Coord;
  int ndouble;
  HitV.clear();
  IngridHitSummary * Hit = new IngridHitSummary();
#ifdef DEBUG2
  cout<<"New track, number of hits="<<recon->NhitTs(itrk)<<endl;
#endif
  for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
    //cout<<"hit="<<ihit<<endl;
    Hit=recon->GetIngridHitTrk(ihit,itrk);
    //cout<<"z position="<<Hit->z<<endl;
    //Coord.x=Hit->view;
    //Coord.y=Hit->xy;
    //Coord.z=Hit->z;
    Coord.view=Hit->view;
    Coord.ch=Hit->ch;
    Coord.pln=Hit->pln;
    Coord.trk=itrk;
    Coord.hit=ihit;
    Coord.mod=Hit->mod;
#ifdef DEBUG2
    cout<<"In erase, mod="<<Coord.mod<<", pln="<<Coord.pln<<", channel="<<Coord.ch<<", view="<<Coord.view<<endl;
#endif
    //Coord.time=Hit->time;
    HitV.push_back(Coord);
    //cout<<"HitV is filled"<<endl;

    for(int ihit2=0;ihit2<HitV.size()-1;ihit2++){
      if(HitV[ihit2]==Coord) {
        ndouble++;
        HitV.pop_back();
      }
    }
  }  
  //cout<<"return HitV"<<endl;
  return HitV;
}

vector <HitTemp> Reconstruction::EraseDoubleHitsAllTracks(IngridBasicReconSummary * recon, vector <HitTemp> HitV){
  //cout<<"hello"<<endl;
  HitV.clear();
  IngridHitSummary * Hit = new IngridHitSummary();
  HitTemp Coord;
  int ndouble;
  for(int itrk=0;itrk<recon->Ntrack;itrk++){
    for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
      Hit=recon->GetIngridHitTrk(ihit,itrk);
      //Coord.x=Hit->view;
      //Coord.y=Hit->xy;
      //Coord.z=Hit->z;
      Coord.view=Hit->view;
      Coord.ch=Hit->ch;
      Coord.pln=Hit->pln;
      Coord.trk=itrk;
      Coord.hit=ihit;
      HitV.push_back(Coord);
      
      for(int ihit2=0;ihit2<HitV.size()-1;ihit2++){
	if(HitV[ihit2]==Coord) {
	  ndouble++;
	  HitV.pop_back();
	}
      }
    }
  }
  
//    Hit=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit].hit,HitV[ihit].trk);
  //cout<<"HitV Size="<<HitV.size()<<endl;
  return HitV;
  // for(int ihit=0;ihit<HitV.size();ihit++){
}


vector <HitTemp> Reconstruction::EraseDoubleHitsAllTracksPM(PMAnaSummary * recon, vector <HitTemp> HitV){
  //cout<<"hello"<<endl;
  HitV.clear();
  IngridHitSummary * Hit = new IngridHitSummary();
  HitTemp Coord;
  int ndouble;
  for(int itrk=0;itrk<recon->Ntrack;itrk++){
    for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
      Hit=recon->GetIngridHitTrk(ihit,itrk);
      //Coord.x=Hit->view;
      //Coord.y=Hit->xy;
      //Coord.z=Hit->z;
      Coord.view=Hit->view;
      Coord.ch=Hit->ch;
      Coord.pln=Hit->pln;     
      Coord.trk=itrk;
      Coord.hit=ihit;
      HitV.push_back(Coord);

      for(int ihit2=0;ihit2<HitV.size()-1;ihit2++){
        if(HitV[ihit2]==Coord) {
          ndouble++;
          HitV.pop_back();
        }
      }
    }
  }

  return HitV;
}


vector <Hit3D> Reconstruction::Hit2DMatching( IngridEventSummary* evt, IngridBasicReconSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec, bool MC){
  Vec.clear();
  IngridHitSummary * hit=new IngridHitSummary();
  IngridHitSummary * hit2=new IngridHitSummary();
  Hit3D hit3d,hit3d2;
  IngridSimParticleSummary * SimPart2;
  //cout<<"hello"<<endl;
  double HitPln[NPlnPM+3][NView+1][NChPM+1];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ich=0;ich<NChPM;ich++){
        HitPln[ipln][iview][ich]=0;
      }
    }
  }
  for(int ihit=0;ihit<HitV.size();ihit++){
    hit=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit].hit,HitV[ihit].trk);

    if(hit->view==1) continue;
    HitTemp T;
    if(hit->cyc==-2) continue;
    for(int ihit2=0;ihit2<HitV.size();ihit2++){
      hit2=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit2].hit,HitV[ihit2].trk);
      if(hit->view==hit2->view) continue;
      if(hit->mod!=hit2->mod) continue;
      if(hit2->cyc==-2) continue;
      if(hit->pln==hit2->pln){

	hit3d.clear();
	hit3d2.clear();
	
        if(hit->mod==16) hit3d.z=(hit->pln-1)*4.6+5.0+2.3*hit->view+0.4;
        else hit3d.z=107.45+hit->pln*10.7+hit->view+1.0;
	//	if(hit->mod==16) hit3d.z=hit->z;
	//else hit3d.z=hit->z+107.45;

	if(MC)  hit3d.pe=hit->pe;
        else{
          if((hit->pe + hit->lope)/2.<39) hit3d.pe=hit->pe;
          else hit3d.pe=hit->lope;
        }
	hit3d.view=hit->view;
	hit3d.pln=hit->pln;
	hit3d.mod=hit->mod;

	hit3d.time=hit->time;	
	hit3d.used=1;
	hit3d.InTrk=true;

	//if(hit->NSimHits()>1)//cout<<"NOOOOOOOOOOOOOOOOOOOOOO"<<endl;
	//if(hit->NSimHits()>1)//cout<<"It shouldn't happen"<<endl;
	if(MC==true && hit->NSimHits()!=0) {
	  hit3d.truepe=((hit->GetIngridSimHit(0))->edeposit)*MeV2PEPM;
	  //cout<<(hit->GetIngridSimHit(0))->pdg<<endl;
	  hit3d.pdg=(hit->GetIngridSimHit(0))->pdg;
	  hit3d.trackid=(hit->GetIngridSimHit(0))->trackid;
	  /*
	  for(int ist=0;ist<evt->NIngridSimParticles();ist++){
	    SimPart = (IngridSimParticleSummary*) evt->GetSimParticle(ist);
	    if(SimPart->trackid==(hit->GetIngridSimHit(0))->trackid) {
	      hit3d.pdg=SimPart->pdg;
	      break;
	    }
	    }*/
	}
	else hit3d.truepe=-1;

	hit3d2.mod=hit2->mod;	
	//if(hit2->mod==16) hit3d2.z=hit2->z;
	//else hit3d2.z=hit2->z+107.45;
        if(hit2->mod==16) hit3d2.z=hit2->pln*4.6+2.3*hit2->view+0.4;
        else hit3d2.z=107.45+hit2->pln*10.7+hit2->view+1.0;

	if(MC) hit3d2.pe=hit2->pe;
	else{
	  if((hit2->pe + hit2->lope)/2.<39) hit3d2.pe=hit2->pe;
	  else hit3d2.pe=hit2->lope;
	}
	hit3d2.pln=hit2->pln;
	hit3d2.view=hit2->view;
	hit3d2.time=hit2->time;
	hit3d2.used=1;
	hit3d2.InTrk=true;
	if(MC && hit2->NSimHits()!=0) {
	  if(hit->NSimHits()>1)//cout<<"It can happen"<<endl;
	  hit3d2.truepe=((hit2->GetIngridSimHit(0))->edeposit)*MeV2PEPM;
	  //cout<<(hit2->GetIngridSimHit(0))->pdg<<endl;
	  hit3d2.pdg=(hit2->GetIngridSimHit(0))->pdg;
	  hit3d2.trackid=(hit2->GetIngridSimHit(0))->trackid;
	  /* for(int ist=0;ist<evt->NIngridSimParticles();ist++){
	    SimPart = (IngridSimParticleSummary*) evt->GetSimParticle(ist);
	    if(SimPart->trackid==(hit2->GetIngridSimHit(0))->trackid) {
	      hit3d2.pdg=SimPart->pdg;
	      break;
	    }
	    }*/
	}

	//if(hit->view==0) {
          hit3d.x=hit->xy;
          hit3d.y=hit2->xy;
          hit3d.ch=hit->ch;
          hit3d2.x=hit->xy;
          hit3d2.y=hit2->xy;
          hit3d2.ch=hit2->ch;
	  /*}
        else{
          hit3d.x=hit2->xy;
          hit3d.y=hit->xy;
          hit3d.ch=hit2->ch;
          hit3d2.x=hit->xy;
          hit3d2.y=hit2->xy;
          hit3d.ch=hit->ch;
	  }*/
	HitPln[hit3d.pln][hit3d.view][hit3d.ch]++;
	HitPln[hit3d2.pln][hit3d2.view][hit3d2.ch]++;
	Vec.push_back(hit3d);
	Vec.push_back(hit3d2);
	//cout<<" Vec, Position=("<<Vec[ihit].x<<","<<Vec[ihit].y<<","<<Vec[ihit].z<<")   , NRJ="<<Vec[ihit].pe<<endl;                                    
	//cout<<" hit3d #1, Position=("<<hit3d.x<<","<<hit3d.y<<","<<hit3d.z<<")   , NRJ="<<hit3d.pe<<endl;                                                    
	//cout<<" hit3d #2, Position=("<<hit3d2.x<<","<<hit3d2.y<<","<<hit3d2.z<<")   , NRJ="<<hit3d2.pe<<endl;                                                
	//cout<<"1st Hit params, view="<<hit->view<<" , xy="<<hit->xy<<" ,z="<<hit->z<<" , NRJ="<<hit->pe<<endl;                                                  
	//             //cout<<" In the hits, Position=("<<hit3d2.x<<","<<hit3d2.y<<","<<hit3d2.z<<")   , NRJ="<<hit3d2.pe<<endl;                                  
	//cout<<"2nd Hit params, view="<<hit2->view<<" , xy="<<hit2->xy<<" ,z="<<hit2->z<<" , NRJ="<<hit2->pe<<endl<<endl;                                        
	//double EnergyTest=FiberAttenuation(hit3d.view,hit3d.x,hit3d.y,hit3d.pe);                                                                                
	//cout<<hit3d.pe<<"   "<<EnergyTest<<endl;
      }
    }
  }

   for(int ihit=0;ihit<Vec.size();ihit++){
    //cout<<Vec[ihit].pln<<" "<<Vec[ihit].ch<<endl;
    //cout<<" Vec, Mod="<<Vec[ihit].mod<<", Position=("<<Vec[ihit].x<<","<<Vec[ihit].y<<","<<Vec[ihit].z<<")   , NRJ="<<Vec[ihit].pe<<" , Number of similar hits="<<HitPln[Vec[ihit].pln][Vec[ihit].view][Vec[ihit].ch]<<endl; 
    //Vec[ihit].pe/=HitPln[Vec[ihit].pln][Vec[ihit].view][Vec[ihit].ch];
    }
  sort(Vec.begin(),Vec.end());
  return Vec;
}


vector <Hit3D> Reconstruction::Hit2DMatchingPM( IngridEventSummary* evt, PMAnaSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec, bool MC){
  
  Vec.clear();
  IngridHitSummary * hit=new IngridHitSummary();
  IngridHitSummary * hit2=new IngridHitSummary();
  Hit3D hit3d,hit3d2;
  IngridSimParticleSummary * SimPart2;

  double HitPln[NPlnPM+3][NView+1][NChPM+1];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ich=0;ich<NChPM;ich++){
        HitPln[ipln][iview][ich]=0;
      }
    }
  }

  
  for(int ihit=0;ihit<HitV.size();ihit++){
    double thetax=(recon->thetax)[HitV[ihit].trk]; double thetay=(recon->thetay)[HitV[ihit].trk];
    double intcptx=(recon->zx)[HitV[ihit].trk]; double intcpty=(recon->zy)[HitV[ihit].trk];
    double gradx=TMath::ATan(thetax*TMath::Pi()/180.);
    double grady=TMath::ATan(thetay*TMath::Pi()/180.);
    
    
    hit=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit].hit,HitV[ihit].trk);
    HitTemp T;
    if(hit->cyc==-2) continue;//throw away the ?? hits
    hit3d.clear();

    hit3d.z=zposi(hit->mod,hit->view,hit->pln)/10.;

#ifdef DEBUG2
    cout<<"**************************************************"<<endl;
    cout<<"Test of Reconstruction::Hit2DMatching"<<endl;
    cout<<"Gradient x="<<gradx<<", intcpt="<<intcptx/10.<<", thetax="<<thetax<<endl;
    cout<<"Gradient y="<<grady<<", intcpt="<<intcpty/10.<<", thetay="<<thetay<<endl;
    cout<<"Hit pln="<<hit->pln<<", z position="<<zposi(hit->mod,hit->view,hit->pln)/10.<<", hit view="<<hit->view<<endl;
  cout<<"**************************************************"<<endl;
#endif

    
    if(MC) hit3d.pe=hit->pe;
    else{
      if((hit->pe + hit->lope)/2.<39) hit3d.pe=hit->pe;
      else hit3d.pe=hit->lope;
    }
    hit3d.view=hit->view;
    hit3d.pln=hit->pln;
    hit3d.mod=hit->mod;
    
    hit3d.time=hit->time;	
    hit3d.used=1;
    hit3d.InTrk=true;
    
    if(MC==true && hit->NSimHits()!=0) {
      if(hit->NSimHits()>1) cout<<"several simulated hits in this scintillator"<<endl;
      hit3d.truepe=((hit->GetIngridSimHit(0))->edeposit)*MeV2PEPM;
      hit3d.pdg=(hit->GetIngridSimHit(0))->pdg;
      hit3d.trackid=(hit->GetIngridSimHit(0))->trackid;
    }
    else hit3d.truepe=-1;
    
    if(hit->view==0){
      hit3d.x=hit->xy;
      hit3d.y=grady*hit3d.z+intcpty/10.;
    }
    else{
      hit3d.y=hit->xy;
      hit3d.x=gradx*hit3d.z+intcptx/10.;
    }

    hit3d.ch=hit->ch;
    HitPln[hit3d.pln][hit3d.view][hit3d.ch]++;
    Vec.push_back(hit3d);
  }

  
  sort(Vec.begin(),Vec.end());
  return Vec;
}

vector <Hit3D> Reconstruction::Hit2DMatchingAllTracksPM(PMAnaSummary * recon, bool MC){

  IngridHitSummary * hit =new IngridHitSummary();
  IngridHitSummary * hit2 =new IngridHitSummary();
  vector <Hit3D> VecAll;
  Hit3D hit3d,hit3d2;
  vector <HitTemp> HitV;

  double HitPln[NPlnPM+3][NView+1][NChPM+1];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ich=0;ich<NChPM;ich++){
        HitPln[ipln][iview][ich]=0;
      }
    }
  }

  IngridHitSummary * Hit = new IngridHitSummary();
  HitTemp Coord;
  int ndouble;

  for(int itrk=0;itrk<recon->Ntrack;itrk++){
    for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
      Hit=recon->GetIngridHitTrk(ihit,itrk);
      Coord.view=Hit->view;
      Coord.ch=Hit->ch;
      Coord.pln=Hit->pln;
      Coord.trk=itrk;
      Coord.hit=ihit;
      HitV.push_back(Coord);

      for(int ihit2=0;ihit2<HitV.size()-1;ihit2++){
        if(HitV[ihit2]==Coord) {
          ndouble++;
          HitV.pop_back();
        }
      }
    }
  }



  VecAll.clear();                                     

  for(int ihit=0;ihit<HitV.size();ihit++){
    //cout<<HitV[ihit].trk<<endl;
    double thetax=(recon->thetax)[HitV[ihit].trk]; double thetay=(recon->thetay)[HitV[ihit].trk];
    double intcptx=(recon->zx)[HitV[ihit].trk]; double intcpty=(recon->zy)[HitV[ihit].trk];
    double gradx=TMath::ATan(thetax*TMath::Pi()/180.);
    double grady=TMath::ATan(thetay*TMath::Pi()/180.);

    hit=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit].hit,HitV[ihit].trk);
    HitTemp T;
    if(hit->cyc==-2) continue;//throw away the ?? hits
    hit3d.clear();

    hit3d.z=zposi(hit->mod,hit->view,hit->pln)/10.;
    
    if(MC) hit3d.pe=hit->pe;
    else{
      if((hit->pe + hit->lope)/2.<39) hit3d.pe=hit->pe;
      else hit3d.pe=hit->lope;
    }
    hit3d.view=hit->view;
    hit3d.pln=hit->pln;
    hit3d.mod=hit->mod;
    
    hit3d.time=hit->time;	
    hit3d.used=1;
    hit3d.InTrk=true;
    
    if(MC==true && hit->NSimHits()!=0) {
      if(hit->NSimHits()>1) cout<<"several simulated hits in this scintillator"<<endl;
      hit3d.truepe=((hit->GetIngridSimHit(0))->edeposit)*MeV2PEPM;
      hit3d.pdg=(hit->GetIngridSimHit(0))->pdg;
      hit3d.trackid=(hit->GetIngridSimHit(0))->trackid;
    }
    else hit3d.truepe=-1;
    if(hit->view==0){
      hit3d.x=hit->xy;
      hit3d.y=grady*hit3d.z+intcpty/10;
    }
    else{
      hit3d.y=hit->xy;
      hit3d.x=gradx*hit3d.z+intcptx/10;
    }

    hit3d.ch=hit->ch;
    HitPln[hit3d.pln][hit3d.view][hit3d.ch]++;
    VecAll.push_back(hit3d);
  }
  
  sort(VecAll.begin(),VecAll.end());

  return VecAll;
}






vector <double> Reconstruction::GetTrackAngle(vector <Hit3D> Vec){
  double dx,AngleX,AngleY;

  int nx=0;
  int ny=0;
  for(int ihit=0;ihit<Vec.size();ihit++){
    if(Vec[ihit].view==0) nx++;
    else ny++;
  }
   double PosZX[nx];
   double PosX[nx];
   double ErrZX[nx];
   double ErrX[nx];

   double PosY[ny];
   double PosZY[ny];
  double ErrY[ny];
  double ErrZY[ny];

  int ix=0;
  int iy=0;

  for(int ihit=0;ihit<Vec.size();ihit++){
    if(Vec[ihit].view==0){
      PosZX[ix]=Vec[ihit].z;
      PosX[ix]=Vec[ihit].x;
      ErrX[ix]=2.5;
      ErrZX[ix]=1.;
      ix++;
    }
    else{
      PosZY[iy]=Vec[ihit].z;
      PosY[iy]=Vec[ihit].y;
      ErrY[iy]=2.5;
      ErrZY[iy]=1.;
      iy++;
      }
  }

  TGraphErrors * HistPosX=new TGraphErrors(nx,PosZX,PosX,ErrZX,ErrX);
  TGraphErrors * HistPosY=new TGraphErrors(ny,PosZY,PosY,ErrZY,ErrY);

  TF1 * fx = new TF1("fx","[0]+x*[1]",0,130);
  TF1 * fy = new TF1("fy","[0]+x*[1]",0,130);
  HistPosX->Fit("fx","Q");
  HistPosY->Fit("fy","Q");

  double SlopeX=fx->GetParameter(1);
  AngleX=TMath::ATan(SlopeX)*180/TMath::Pi();
  double SlopeY=fy->GetParameter(1);
  AngleY=TMath::ATan(SlopeY)*180/TMath::Pi();
  double bX=fx->GetParameter(0);
  double bY=fy->GetParameter(0);
  double ReducedChiSquareX=fx->GetChisquare()/fx->GetNDF();
  double ReducedChiSquareY=fy->GetChisquare()/fy->GetNDF();
  //cout<<"x chisquare="<<fx->GetChisquare()<<", NDF="<<fx->GetNDF()<<endl<<"nx="<<nx<<endl;
 
  //This formula for dx is wrong
  dx=TMath::Sqrt(SlopeY*SlopeY+SlopeX*SlopeX+1);
  vector <double> Out;
  Out.push_back(AngleX);
  Out.push_back(AngleY);
  Out.push_back(dx);
  Out.push_back(ReducedChiSquareX);
  Out.push_back(ReducedChiSquareY);
  Out.push_back(SlopeX);
  Out.push_back(bX);
  Out.push_back(SlopeY);
  Out.push_back(bY);
  HistPosX->Delete();
  HistPosY->Delete();
  fx->Delete();
  fy->Delete();
  return Out;
}


vector <double> Reconstruction::GetTrackAnglePM(vector <Hit3D> Vec, double AngleX, double AngleY, int TrackSample){
  vector <double> Out;
  Out.clear();

  if(TrackSample<5){
    double ZX=TMath::Cos(AngleX*TMath::Pi()/180);
    double X=TMath::Sin(AngleX*TMath::Pi()/180);
    double ZY=TMath::Cos(AngleY*TMath::Pi()/180);
    double Y=TMath::Sin(AngleY*TMath::Pi()/180);

    double X2=X/ZX;
    double Y2=Y/ZY;
    double Z2=1;
    double Norm0=TMath::Sqrt(X2*X2+Y2*Y2+Z2*Z2);
    X2/=Norm0;
    Y2/=Norm0;
    Z2/=Norm0;
    double Scalar2=X2*Beam[0]+Y2*Beam[1]+Z2*Beam[2];
    //double Norm=TMath::Sqrt((Fpos[0]-Ipos[0])*(Fpos[0]-Ipos[0])+(Fpos[1]-Ipos[1])*(Fpos[1]-Ipos[1])+(Fpos[2]-Ipos[2])*(Fpos[2]-Ipos[2]));
    double Angle2=TMath::ACos(Scalar2)*180/TMath::Pi();
    Out.push_back(Angle2);
  }

  else{
 int nx=0;int ny=0;
  int ix=0;int iy=0;
  for(int ihit=0;ihit<Vec.size();ihit++){
    if(Vec[ihit].view==0) nx++;
    else ny++;
  }

  double PosZX[nx];
  double PosX[nx];
  double ErrZX[nx];
  double ErrX[nx];
  
  double PosY[ny];
  double PosZY[ny];
  double ErrY[ny];
  double ErrZY[ny];

  for(int ihit=0;ihit<Vec.size();ihit++){
    //cout<<"Module="<<Vec[ihit].mod<<", x="<<Vec[ihit].x<<", y="<<Vec[ihit].y<<", z="<<Vec[ihit].z<<endl;
    if(Vec[ihit].mod==16){
      if(Vec[ihit].view==0){
	PosZX[ix]=Vec[ihit].z;
	PosX[ix]=Vec[ihit].x;
	if(Vec[ihit].ch<=7||Vec[ihit].ch>=24){
	  ErrX[ix]=2.5;
	  ErrZX[ix]=1.3;
	}
	else{
	  ErrX[ix]=5.;
	  ErrZX[ix]=1.;
	}
	ix++;
      }

      else{
	PosZY[iy]=Vec[ihit].z;
	PosY[iy]=Vec[ihit].y;
	if(Vec[ihit].ch<=7||Vec[ihit].ch>=24){
	  ErrY[iy]=2.5;
	  ErrZY[iy]=1.3;
	}
	else{
	  ErrY[iy]=5.;
	  ErrZY[iy]=1.;
	}
	iy++;
      }
    }
    else{
     if(Vec[ihit].view==0){
	PosZX[ix]=Vec[ihit].z;
	PosX[ix]=Vec[ihit].x;
	ErrX[ix]=5.;
	ix++;
      }
      else{
	PosZY[iy]=Vec[ihit].z;
	PosY[iy]=Vec[ihit].y;
	ErrY[iy]=5;
	iy++;
      }
    }
  }

  vector <double> MatchingPMINGRID;
  MatchingPMINGRID.clear();
   
  TGraphErrors * HistPosX=new TGraphErrors(nx,PosZX,PosX,ErrZX,ErrX);
  TGraphErrors * HistPosY=new TGraphErrors(ny,PosZY,PosY,ErrZY,ErrY);
  //TGraphErrors * HistPosX_Ing=new TGraphErrors(nx_Ing,PosZX_Ing,PosX_Ing,ErrZX_Ing,ErrX_Ing);
  //TGraphErrors * HistPosY_Ing=new TGraphErrors(ny_Ing,PosZY_Ing,PosY_Ing,ErrZY_Ing,ErrY_Ing);

  TF1 * fx = new TF1("fx","[0]+x*[1]",0,400);
  TF1 * fy = new TF1("fy","[0]+x*[1]",0,400);
  HistPosX->Fit("fx","Q");
  HistPosY->Fit("fy","Q");


  double Fpos[3];
  Fpos[0]=fx->Eval(400);
  Fpos[1]=fy->Eval(400);
  Fpos[2]=400;
  
  double Ipos[3];
  Ipos[0]=fx->Eval(0);
  Ipos[1]=fy->Eval(0);
  Ipos[2]=0;
  
  double Scalar=(Fpos[0]-Ipos[0])*Beam[0]+(Fpos[1]-Ipos[1])*Beam[1]+(Fpos[2]-Ipos[2])*Beam[2];
  double Norm=TMath::Sqrt((Fpos[0]-Ipos[0])*(Fpos[0]-Ipos[0])+(Fpos[1]-Ipos[1])*(Fpos[1]-Ipos[1])+(Fpos[2]-Ipos[2])*(Fpos[2]-Ipos[2]));
  double Angle=TMath::ACos(Scalar/Norm)*180/TMath::Pi();
  Out.push_back(Angle);

  } 
  return Out;
  /* TrueAngleMuon=TMath::ACos(Scalar/Norm)*180/TMath::Pi();

  double SlopeX=fx->GetParameter(1);
  double AngleX=TMath::ATan(SlopeX)*180/TMath::Pi();
  double SlopeY=fy->GetParameter(1);
  double AngleY=TMath::ATan(SlopeY)*180/TMath::Pi();
  double ReducedChiSquareX=fx->GetChisquare()/fx->GetNDF();
  double ReducedChiSquareY=fy->GetChisquare()/fy->GetNDF();

  vector <double> Out;
  Out.push_back(AngleX);
  Out.push_back(AngleY);
  Out.push_back(dx);

  Out.push_back(SlopeX);
  Out.push_back(bX);
  Out.push_back(Zi[0]);
  Out.push_back(Zf[0]);

  Out.push_back(SlopeY);
  Out.push_back(bY);
  Out.push_back(Zi[1]); 
  Out.push_back(Zf[1]);

  HistPosX->Delete();   
  HistPosY->Delete();
  Norm->Delete();
  fx->Delete();
  fy->Delete();
  return Out;
  */
}

/**Provides two track opening angle and coplanarity**/
vector <double> Reconstruction::GetKinematic(double ang1, double thetax1, double thetay1,double ang2, double thetax2, double thetay2){
  double * Track1 = new double[3];
  double * Track2 = new double[3];
  vector <double> Kinematic;
  Kinematic.clear();
  double AngleInTracks,Determinant;

  Track1[0]=TMath::Tan(thetax1)*TMath::Cos(ang1);
  Track1[1]=TMath::Tan(thetay1)*TMath::Cos(ang1);
  Track1[2]=TMath::Cos(ang1);
  double NormTrack1=TMath::Sqrt(Track1[0]*Track1[0]+Track1[1]*Track1[1]+Track1[2]*Track1[2]);

  Track2[0]=TMath::Tan(thetax2)*TMath::Cos(ang2);
  Track2[1]=TMath::Tan(thetay2)*TMath::Cos(ang2);
  Track2[2]=TMath::Cos(ang2);
  double NormTrack2=TMath::Sqrt(Track2[0]*Track2[0]+Track2[1]*Track2[1]+Track2[2]*Track2[2]);

  for(int i=0;i<3;i++){
    //cout<<i<<": (";
    Track1[i]/=NormTrack1;
    Track2[i]/=NormTrack2;
    //cout<< Track1[i]<<", "<< Track2[i]<<endl;
  }

  AngleInTracks=TMath::ACos(Track1[0]*Track2[0]+Track1[1]*Track2[1]+Track1[2]*Track2[2]);
  //cout<<"Angle between Tracks="<<AngleInTracks*180/TMath::Pi()<<endl;

  Kinematic.push_back(AngleInTracks);
  
  //Construct the coplanarity angle.
  //a. neutrino beam direction
  double * neutrinobeam = new double[3];
  neutrinobeam[0]=0;
  neutrinobeam[1]=-TMath::Sin(3.8*TMath::Pi()/180.);
  neutrinobeam[2]=TMath::Cos(3.8*TMath::Pi()/180.);
  double Nneutrinobeam=TMath::Sqrt(neutrinobeam[0]*neutrinobeam[0]+neutrinobeam[1]*neutrinobeam[1]+neutrinobeam[2]*neutrinobeam[2]);
  for(int i=0;i<3;i++) neutrinobeam[i]/=Nneutrinobeam;
  
  //b. Construct the projections of first & second tracks on the plan orthogonal to the beam direction
  double * firsttrackproj = new double[3];
  TMath::NormCross(neutrinobeam,Track1,firsttrackproj);
  double * secondtrackproj = new double[3];
  TMath::NormCross(neutrinobeam,Track2,secondtrackproj);

  //c. Determine the angle before these two projected vectors.
  double CoplanarityAngle=TMath::ACos(firsttrackproj[0]*secondtrackproj[0]+firsttrackproj[1]*secondtrackproj[1]+firsttrackproj[2]*secondtrackproj[2]);

#ifdef DEBUG2
  cout<<"**************************************************"<<endl;
  cout<<"Test of Reconstruction::GetKinematic"<<endl;
  cout<<"neutrino direction vector=("<<neutrinobeam[0]<<","<<neutrinobeam[1]<<","<<neutrinobeam[2]<<")"<<endl;
  cout<<"Test track 1, original vector=("<<Track1[0]<<","<<Track1[1]<<","<<Track1[2]<<") and projected=("<<firsttrackproj[0]<<","<<firsttrackproj[1]<<","<<firsttrackproj[2]<<")"<<endl;
  cout<<"Test track 2, original vector=("<<Track2[0]<<","<<Track2[1]<<","<<Track2[2]<<") and projected=("<<secondtrackproj[0]<<","<<secondtrackproj[1]<<","<<secondtrackproj[2]<<")"<<endl;
  cout<<"Coplanarity Angle="<<CoplanarityAngle*180/TMath::Pi()<<endl;
  cout<<"**************************************************"<<endl;
#endif
  /*
  //const double X[]={Beam[0],Track1[0],Track2[0],Beam[1],Track1[1],Track2[1],Beam[2],Track1[2],Track2[2]};
  //TMatrixD M(3,3);
  //M.SetMatrixArray(X);
  double ProdVec[3],ProdScal,Angle3D;

  ProdVec[0]=Track1[1]*Track2[2]-Track2[1]*Track1[2];
  ProdVec[1]=Track1[2]*Track2[0]-Track2[2]*Track1[0];
  ProdVec[2]=Track1[0]*Track2[1]-Track2[0]*Track1[1];
  double NormProdVec=TMath::Sqrt(ProdVec[0]*ProdVec[0]+ProdVec[1]*ProdVec[1]+ProdVec[2]*ProdVec[2]);
  double NormBeam=TMath::Sqrt(Beam[0]*Beam[0]+Beam[1]*Beam[1]+Beam[2]*Beam[2]);

  ProdScal=(Beam[0]*ProdVec[0]+Beam[1]*ProdVec[1]+Beam[2]*ProdVec[2])/NormProdVec;
  Angle3D=TMath::ACos(ProdScal);
  //Determinant=TMath::Abs(M.Determinant());
  //cout<<"Determinant"<<Determinant<<endl;
  //Kinematic.push_back(Determinant);
 
  //cout<<"Norm Track1="<<NormTrack1<<" , Norm Track2="<<NormTrack2<<endl;
  //cout<<"Norm Beam="<<NormBeam<<" , Norm ProdVec="<<NormProdVec<<endl;
  //cout<<"ProdScal="<<ProdScal<<endl;*/
  Kinematic.push_back(CoplanarityAngle);

  delete Track1;
  delete Track2;
  delete neutrinobeam;
  delete firsttrackproj;
  delete secondtrackproj;
  
  return Kinematic;
}
/********************************************************************************/

double Reconstruction::GetBeamAngle(double ang1, double thetax1, double thetay1){
  float Track1[3];

  double AngleBeamTrack;

  Track1[0]=TMath::Tan(thetax1)*TMath::Cos(ang1);
  Track1[1]=TMath::Tan(thetay1)*TMath::Cos(ang1);
  Track1[2]=TMath::Cos(ang1);
  double NormTrack1=TMath::Sqrt(Track1[0]*Track1[0]+Track1[1]*Track1[1]+Track1[2]*Track1[2]);

  for(int i=0;i<3;i++){
    Track1[i]/=NormTrack1;
  }

  AngleBeamTrack=TMath::ACos(Track1[0]*Beam[0]+Track1[1]*Beam[1]+Track1[2]*Beam[2]);

  return AngleBeamTrack;
}




vector <Hit3D> Reconstruction::SeveralHitsPlane(IngridBasicReconSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec){
  vector <Hit3D> Vec2;
  Hit3D hit3d, hit3d2, hitall;
  
  bool used[Vec.size()];
  for(int ihit=0;ihit<Vec.size();ihit++){
    used[ihit]=false;
  }
  
  for(int ihit=0;ihit<Vec.size();ihit++){
    if(used[ihit]==true) continue;
    hitall.x+=Vec[ihit].x;
    hitall.y+=Vec[ihit].y;
    hitall.z+=Vec[ihit].z;
    hitall.pe+=Vec[ihit].pe;
    hitall.pecorr+=Vec[ihit].pecorr;
    hitall.view=Vec[ihit].view;
    hitall.pln=Vec[ihit].pln;
    int IClster=1;
    for(int ihit2=ihit+1;ihit2<Vec.size();ihit2++){
      if(used[ihit2]==true) continue;
      if((Vec[ihit].pln==Vec[ihit2].pln)&&(Vec[ihit].view==Vec[ihit2].view)&&(TMath::Abs(Vec[ihit].x-Vec[ihit2].x)<10 && (TMath::Abs(Vec[ihit].y-Vec[ihit2].y)<10))){
	hitall.x+=Vec[ihit2].x;
	hitall.y+=Vec[ihit2].y;
	hitall.z+=Vec[ihit2].z;
	hitall.pe+=Vec[ihit2].pe;
	hitall.pecorr+=Vec[ihit2].pecorr;
	IClster++;
	used[ihit2]=true;
      }
  }
    hitall.x/=IClster;
    hitall.y/=IClster;
    hitall.z/=IClster;
    Vec2.push_back(hitall);
    hitall.clear();
  }
  
  hit3d.clear();
  hit3d2.clear();
  hitall.clear();
  return Vec2;
}


vector <Hit3D> Reconstruction::SeveralHitsPlanePM(PMAnaSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec){
  vector <Hit3D> Vec2;
  Hit3D hit3d, hit3d2, hitall;
  
  bool used[Vec.size()];
  for(int ihit=0;ihit<Vec.size();ihit++){
    used[ihit]=false;
  }
  
  for(int ihit=0;ihit<Vec.size();ihit++){
    if(used[ihit]==true) continue;
    hitall.x+=Vec[ihit].x;
    hitall.y+=Vec[ihit].y;
    hitall.z+=Vec[ihit].z;
    hitall.pe+=Vec[ihit].pe;
    hitall.pecorr+=Vec[ihit].pecorr;
    hitall.view=Vec[ihit].view;
    hitall.pln=Vec[ihit].pln;
    int IClster=1;
    for(int ihit2=ihit+1;ihit2<Vec.size();ihit2++){
      if(used[ihit2]==true) continue;
      if((Vec[ihit].pln==Vec[ihit2].pln)&&(Vec[ihit].view==Vec[ihit2].view)&&(TMath::Abs(Vec[ihit].x-Vec[ihit2].x)<10 && (TMath::Abs(Vec[ihit].y-Vec[ihit2].y)<10))){
	hitall.x+=Vec[ihit2].x;
	hitall.y+=Vec[ihit2].y;
	hitall.z+=Vec[ihit2].z;
	hitall.pe+=Vec[ihit2].pe;
	hitall.pecorr+=Vec[ihit2].pecorr;
	IClster++;
	used[ihit2]=true;
      }
  }
    hitall.x/=IClster;
    hitall.y/=IClster;
    hitall.z/=IClster;
    Vec2.push_back(hitall);
    hitall.clear();
  }
  
  hit3d.clear();
  hit3d2.clear();
  hitall.clear();
  return Vec2;
}




vector <Hit3D> Reconstruction::Hit2DMatchingCluster(IngridEventSummary* evt, IngridBasicReconSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec){
  //cout<<"Change Z and X/Y of hits in Hit2DMatchingClusterPM before to use"<<endl;
  int Mod=recon->hitmod;
  IngridHitSummary * hit =new IngridHitSummary();
  IngridHitSummary * hit2 =new IngridHitSummary();
  vector <Hit3D> VecAll;
  Hit3D hit3d,hit3d2;
  double TCluster= recon->clstime;

  double HitPln[NPlnPM+3][NView+1][NChPM+1];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ich=0;ich<NChPM;ich++){
        HitPln[ipln][iview][ich]=0;
      }
    }
  }

  VecAll.clear();                                                                                                   
  for(int ihit=0;ihit<evt->NIngridHits();ihit++){
    hit=(IngridHitSummary*) evt->GetIngridHit(ihit);//changer si data                                                                                           
    if((TMath::Abs(hit->time-TCluster)>50)) continue;
    if(hit->view==1) continue;
    HitTemp T;
    
    for(int ihit2=0;ihit2<evt->NIngridHits();ihit2++){
      hit2=(IngridHitSummary*) evt->GetIngridHit(ihit2);
      if(hit->view==hit2->view) continue;
      if((TMath::Abs(hit->time-TCluster)>50)) continue;
      if(hit->mod!=hit2->mod) continue;
      if(hit->pln!=hit2->pln) continue;
	   //cout<<"2nd Hit params, view="<<hit2->view<<" , xy="<<hit2->xy<<" ,z="<<hit2->z<<endl;                                                                 
	 hit3d.clear();
	 hit3d2.clear();
	 if(hit->mod==16) hit3d.z=hit->pln*4.6+2.3*hit->view+0.4;
	 else hit3d.z=107.45+hit->pln*10.7+hit->view+1.0;

	 if((hit->pe + hit->lope)/2.<39) hit3d.pe=hit->pe;
	 else hit3d.pe=hit->lope;     
	 hit3d.view=hit->view;
	 hit3d.pln=hit->pln;
	 hit3d.mod=hit->mod;	 
	 hit3d.time=hit->time;
	 hit3d.ch=hit->ch;

	 if(hit->mod==16) hit3d.z=hit->pln*4.6+2.3*hit->view+0.4;
	 else hit3d.z=107.45+hit->pln*10.7+hit->view+1.0;

	 if((hit2->pe + hit2->lope)/2.<39) hit3d2.pe=hit2->pe;
	 else hit3d2.pe=hit2->lope;
	 hit3d2.pln=hit2->pln;
	 hit3d2.view=hit2->view;
	 hit3d2.time=hit2->time;
	 hit3d2.mod=hit2->mod;
	 hit3d2.ch=hit2->ch;
	 //used[ihit2]==true;                                                                                                                                    
	 if(hit->mod==16){
	   if(hit->ch<=7) hit->xy=5*hit->ch;
	   else if(hit->ch<=23) hit->xy=5*8+2.5*(hit->ch-8);
	   else hit->xy=5*8+2.5*16+5*(hit->ch-24);
	 }

	 if(hit2->mod==16){
	   if(hit2->ch<=7) hit2->xy=5*hit2->ch;
	   else if(hit2->ch<=23) hit2->xy=5*8+2.5*(hit2->ch-8);
	   else hit2->xy=5*8+2.5*16+5*(hit2->ch-24);
	 }

	 if(hit->view==0) {
	   hit3d.x=hit->xy;
	   hit3d.y=hit2->xy;
	   hit3d2.x=hit->xy;
	   hit3d2.y=hit2->xy;
	 }
	 else{
	   hit3d.x=hit2->xy;
	   hit3d.y=hit->xy;
	   hit3d2.x=hit->xy;
	   hit3d2.y=hit2->xy;
	 }
	 HitPln[hit3d.pln][hit3d.view][hit3d.ch]++;
	 HitPln[hit3d2.pln][hit3d2.view][hit3d2.ch]++;
	 VecAll.push_back(hit3d);
	 VecAll.push_back(hit3d2);
    }
  }
  /*  for(int ihit=0;ihit<VecAll.size();ihit++){
    VecAll[ihit].pe/=HitPln[VecAll[ihit].pln][VecAll[ihit].view][VecAll[ihit].ch];
  }
  */
  sort(VecAll.begin(),VecAll.end());
  hit->Clear();
  hit2->Clear();
  return VecAll;
}



vector <Hit3D> Reconstruction::Hit2DMatchingClusterPM(IngridEventSummary* evt, PMAnaSummary * recon){
  //cout<<"Change Z and X/Y of hits in Hit2DMatchingClusterPM before to use"<<endl;
  IngridHitSummary * hit =new IngridHitSummary();
  IngridHitSummary * hit2 =new IngridHitSummary();
  vector <Hit3D> VecAll;
  Hit3D hit3d,hit3d2;
  double TCluster= recon->clstime;

  double HitPln[NPlnPM+3][NView+1][NChPM+1];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ich=0;ich<NChPM;ich++){
        HitPln[ipln][iview][ich]=0;
      }
    }
  }

  VecAll.clear();                                                                                                   
  for(int ihit=0;ihit<evt->NIngridHits();ihit++){
    hit=(IngridHitSummary*) evt->GetIngridHit(ihit);//changer si data                                                                                           
    if((TMath::Abs(hit->time-TCluster)>50)) continue;
    if(hit->view==1) continue;
    HitTemp T;
    
    for(int ihit2=0;ihit2<evt->NIngridHits();ihit2++){
      hit2=(IngridHitSummary*) evt->GetIngridHit(ihit2);
      if(hit->view==hit2->view) continue;
      if((TMath::Abs(hit->time-TCluster)>50)) continue;
      if(hit->mod!=hit2->mod) continue;
      if(hit->pln!=hit2->pln) continue;
	   //cout<<"2nd Hit params, view="<<hit2->view<<" , xy="<<hit2->xy<<" ,z="<<hit2->z<<endl;                                                                 
	 hit3d.clear();
	 hit3d2.clear();
	 //if(hit->mod==3 && hit2->pln==1) cout<<"View0, channel="<<hit->ch<<", view1 channel="<<hit2->ch<<endl;

	 if(hit2->mod==16) hit3d.z=hit->z;
	 else hit3d.z=hit->z+107.45;

	 if((hit->pe + hit->lope)/2.<39) hit3d.pe=hit->pe;
	 else hit3d.pe=hit->lope;     
	 hit3d.view=hit->view;
	 hit3d.pln=hit->pln;
	 hit3d.mod=hit->mod;	 
	 hit3d.time=hit->time;
	 hit3d.ch=hit->ch;

	 if(hit2->mod==16) hit3d2.z=hit2->z;
	 else hit3d2.z=hit2->z+107.45;
	 if((hit2->pe + hit2->lope)/2.<39) hit3d2.pe=hit2->pe;
	 else hit3d2.pe=hit2->lope;
	 hit3d2.pln=hit2->pln;
	 hit3d2.view=hit2->view;
	 hit3d2.time=hit2->time;
	 hit3d2.mod=hit2->mod;
	 hit3d2.ch=hit2->ch;
	 //used[ihit2]==true;                                                                                                                                    
	 if(hit->mod==16){
	   if(hit->ch<=7) hit->xy=5*hit->ch;
	   else if(hit->ch<=23) hit->xy=5*8+2.5*(hit->ch-8);
	   else hit->xy=5*8+2.5*16+5*(hit->ch-24);
	 }

	 if(hit2->mod==16){
	   if(hit2->ch<=7) hit2->xy=5*hit2->ch;
	   else if(hit2->ch<=23) hit2->xy=5*8+2.5*(hit2->ch-8);
	   else hit2->xy=5*8+2.5*16+5*(hit2->ch-24);
	 }

	 if(hit->view==0) {
	   hit3d.x=hit->xy;
	   hit3d.y=hit2->xy;
	   hit3d2.x=hit->xy;
	   hit3d2.y=hit2->xy;
	 }
	 else{
	   hit3d.x=hit2->xy;
	   hit3d.y=hit->xy;
	   hit3d2.x=hit->xy;
	   hit3d2.y=hit2->xy;
	 }
	 HitPln[hit3d.pln][hit3d.view][hit3d.ch]++;
	 HitPln[hit3d2.pln][hit3d2.view][hit3d2.ch]++;
	 VecAll.push_back(hit3d);
	 VecAll.push_back(hit3d2);
    }
  }
  /*  for(int ihit=0;ihit<VecAll.size();ihit++){
    VecAll[ihit].pe/=HitPln[VecAll[ihit].pln][VecAll[ihit].view][VecAll[ihit].ch];
  }
  */
  sort(VecAll.begin(),VecAll.end());
  hit->Clear();
  hit2->Clear();
  return VecAll;
}




vector <Hit3D> Reconstruction::ClusterPM(IngridEventSummary* evt, PMAnaSummary * recon, int nTracks){
  int Mod=16;
  IngridHitSummary * hit =new IngridHitSummary();
  vector <Hit3D> VecAll;
  vector <HitTemp> HitV;
  HitV.clear();

  Hit3D hit3d,hit3d2;
  double TCluster= recon->clstime;

  IngridHitSummary * hit2 = new IngridHitSummary();
  HitTemp Coord;
  int ndouble;

  VecAll.clear();                                                                                                   
  for(int ihit=0;ihit<evt->NIngridHits();ihit++){ 
    hit=(IngridHitSummary*) evt->GetIngridHit(ihit);//changer si data
    HitTemp T;
    //cout<<"Mod="<<hit->mod<<", pln="<<hit->pln<<", view="<<hit->view<<", ch="<<hit->ch<<endl;
    hit3d.clear();

    hit3d.mod=hit->mod;
    hit3d.z=hit->z;
    if((hit->pe + hit->lope)/2.<30) hit3d.pe=hit->pe;
    else hit3d.pe=hit->lope;
    hit3d.view=hit->view;
    hit3d.pln=hit->pln;
    hit3d.time=hit->time;
    hit3d.ch=hit->ch;

    for(int itrk=0;itrk<nTracks;itrk++){
      for(int ihit2=0;ihit2<recon->NhitTs(itrk);ihit2++){
	hit2=recon->GetIngridHitTrk(ihit2,itrk);
	if(hit3d.mod==hit2->mod && hit3d.view==hit2->view && hit3d.pln==hit2->pln && hit3d.ch==hit2->ch) {
	  hit3d.InTrk=true;
	  hit3d.RecTrk.push_back(itrk+1);
	}
      }
    }
    
    VecAll.push_back(hit3d);
  }

  sort(VecAll.begin(),VecAll.end());
  /*  for(int itrk=0;itrk<nTracks;itrk++){
    for(int ihit=0;ihit<VecAll.size();ihit++){
      if(VecAll[ihit].InTrk==false) continue;
      bool GoodTrk=false;
      for(int t=0;t<VecAll[ihit].RecTrk.size();t++){
	if(VecAll[ihit].RecTrk[t]==(itrk+1)) GoodTrk=true;
      }
      if(!GoodTrk) continue;
      cout<<"Hit is in Trk="<<itrk<<", in mod="<<VecAll[ihit].mod<<", in pln="<<VecAll[ihit].pln<<", in view="<<VecAll[ihit].view<<", in ch="<<VecAll[ihit].ch<<endl;
    }
    }*/
  hit->Clear();
  return VecAll;
}

vector <double> Reconstruction::GetAllTracksLength(vector <Hit3D> Vec){
  double XiIng=1000;double YiIng=1000; double ZiIng=1000;double ZiPM=1000; double YiPM=1000;double XiPM=1000;
  double XfIng=0;double YfIng=0; double ZfIng=0;double ZfPM=0; double YfPM=0;double XfPM=0;

  for(int ihit=0;ihit<Vec.size();ihit++){
    if(Vec[ihit].mod!=16){
      if(Vec[ihit].x<XiIng) XiIng=Vec[ihit].x;
      if(Vec[ihit].y<YiIng) YiIng=Vec[ihit].y;
      if(Vec[ihit].z<ZiIng) ZiIng=Vec[ihit].z;

      if(Vec[ihit].x>XfIng) XfIng=Vec[ihit].x;
      if(Vec[ihit].y>YfIng) YfIng=Vec[ihit].y;
      if(Vec[ihit].z>ZfIng) ZfIng=Vec[ihit].z;
    }
    else{
      if(Vec[ihit].x<XiPM) XiPM=Vec[ihit].x;
      if(Vec[ihit].y<YiPM) YiPM=Vec[ihit].y;
      if(Vec[ihit].z<ZiPM) ZiPM=Vec[ihit].z;

      if(Vec[ihit].x>XfPM) XfPM=Vec[ihit].x;
      if(Vec[ihit].y>YfPM) YfPM=Vec[ihit].y;
      if(Vec[ihit].z>ZfPM) ZfPM=Vec[ihit].z;
    }
  }
    vector <double> Length;
    Length.push_back(TMath::Sqrt(pow(XfPM-XiPM,2)+pow(YfPM-YiPM,2)+pow(ZfPM-ZiPM,2)));
    Length.push_back(TMath::Sqrt(pow(XfIng-XiIng,2)+pow(YfIng-YiIng,2)+pow(ZfIng-ZiIng,2)));
    return Length;
}

double Reconstruction::GetTrackLength(vector <Hit3D> Vec){
  double TrkLength=TMath::Sqrt(pow((Vec.front().x-Vec.back().x),2)+pow((Vec.front().y-Vec.back().y),2)+pow((Vec.front().z-Vec.back().z),2));
  return TrkLength;
}

vector <double> Reconstruction::GetVertex(vector <Hit3D> Vec){
  vector <double> Vertex;
  Vertex.push_back(Vec.front().x);
  Vertex.push_back(Vec.front().y);
  Vertex.push_back(Vec.front().z);
  return Vertex;
}



/**Estimates the total charge of the track, corrected by the track angle (dx)**/
double Reconstruction::GetTrackEnergyPerDistance(vector <Hit3D> Vec,double dx){
  double Energy=0;
  for(int ihit=0;ihit<Vec.size();ihit++){
    //if(Vec[ihit].used>1) continue;
    if(!IsINGRID(Vec[ihit].mod,Vec[ihit].pln,Vec[ihit].ch)) Energy+=Vec[ihit].pe/2;
       //   int pln,int ch)
    //if(Vec[ihit].ch<=7||Vec[ihit].ch>=24) Energy+=Vec[ihit].pe/2;
    else Energy+=Vec[ihit].pe;
  }
  if(dx!=0) Energy/=dx;
  return Energy;
}
/******************************************************************************/


vector <int> Reconstruction::TrackComposition(vector <Hit3D> Vec, double VertexX,double VertexY,double VertexZ){
  vector <int> Compo;
  int NMuons=0;
  int NPions=0;
  int NProtons=0;
  sort(Vec.begin(),Vec.end());
  //cout<<"New Track:"<<endl;
  for(int ihit=0/*Vec.size()-4*/;ihit<Vec.size();ihit++){
    //if(Vec[ihit].mod==16) continue;
    double Dist=TMath::Sqrt(pow(Vec[ihit].x-VertexX,2)+pow(Vec[ihit].y-VertexY,2)+pow(Vec[ihit].z-VertexZ,2));
    if(Vec[ihit].used>1) continue;//skip hits too close from vertex
    //cout<<"Particle Type="<<Vec[ihit].pdg<<" Distance from Vertex="<<Dist<<endl;
    if(Dist<10) continue;
    if(Vec[ihit].pdg==13) NMuons++;
    else if(Vec[ihit].pdg==211) NPions++;
    else if(Vec[ihit].pdg==2212) NProtons++;
  }

  if(NMuons==0 && NPions==0 && NProtons==0) {//case where all hits are shared
    //cout<<"All Hits are shared"<<endl;
    for(int ihit=0;ihit<Vec.size();ihit++){
      //if(Vec[ihit].mod==16) continue;
      double Dist=TMath::Sqrt(pow(Vec[ihit].x-VertexX,2)+pow(Vec[ihit].y-VertexY,2)+pow(Vec[ihit].z-VertexZ,2));                                                    
      //cout<<"Particle Type="<<Vec[ihit].pdg<<" Distance from Vertex="<<Dist<<endl;
      if(Vec[ihit].pdg==13) NMuons++;
      else if(Vec[ihit].pdg==211) NPions++;
      else if(Vec[ihit].pdg==2212) NProtons++;
    }
  }
  /*
  int Tot=NMuons+NPions+NProtons;
  if(Tot!=0) cout<<"Percentage of Muon="<<(float) 100*NMuons/Tot<<", Pions="<<(float) 100*NPions/Tot<<", Protons="<<(float) 100*NProtons/Tot<<endl;
  else cout<<"NO Hits of each type"<<endl;
  */
  Compo.push_back(NMuons);
  Compo.push_back(NPions);
  Compo.push_back(NProtons);
  
  return Compo;
}


vector <Hit3D> Reconstruction::Hit2DMatchingAllTracks(IngridBasicReconSummary * recon){
  int Mod=16;
  IngridHitSummary * hit =new IngridHitSummary();
  IngridHitSummary * hit2 =new IngridHitSummary();
  vector <Hit3D> VecAll;
  Hit3D hit3d,hit3d2;
  vector <HitTemp> HitV;

  double HitPln[NPlnPM+3][NView+1][NChPM+1];
  for(int ipln=0;ipln<NPlnPM;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ich=0;ich<NChPM;ich++){
        HitPln[ipln][iview][ich]=0;
      }
    }
  }

  IngridHitSummary * Hit = new IngridHitSummary();
  HitTemp Coord;
  int ndouble;
  // //cout<<"Ntracks="<<recon->Ntrack<<endl;
  for(int itrk=0;itrk<recon->Ntrack;itrk++){
    //if( recon->NhitTs(itrk)>100)//cout<<"NHits Track="<<recon->NhitTs(itrk)<<endl;
    for(int ihit=0;ihit<recon->NhitTs(itrk);ihit++){
      Hit=recon->GetIngridHitTrk(ihit,itrk);
      //Coord.x=Hit->view;
      //Coord.y=Hit->xy;
      //Coord.z=Hit->z;
      Coord.view=Hit->view;
      Coord.ch=Hit->ch;
      Coord.pln=Hit->pln;
      Coord.trk=itrk;
      Coord.hit=ihit;
      HitV.push_back(Coord);

      for(int ihit2=0;ihit2<HitV.size()-1;ihit2++){
        if(HitV[ihit2]==Coord) {
          ndouble++;
          HitV.pop_back();
        }
      }
    }
  }
  /*cout<<"Vector HitV:"<<endl;
   for(int ihit=0;ihit<HitV.size();ihit++){
   //cout<<HitV[ihit].x<<"  "<<HitV[ihit].y<<"  "<<HitV[ihit].z<<endl;
    }*/

  VecAll.clear();                                     
  for(int ihit=0;ihit<HitV.size();ihit++){
    hit=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit].hit,HitV[ihit].trk);
    if(hit->view==1) continue;
    HitTemp T;
    if(hit->cyc==-2) continue;
    for(int ihit2=0;ihit2<HitV.size();ihit2++){
      hit2=(IngridHitSummary*) recon->GetIngridHitTrk(HitV[ihit2].hit,HitV[ihit2].trk);
      if(hit->view==hit2->view) continue;
      if(hit->mod!=hit2->mod) continue;
      if(hit2->cyc==-2) continue;
      if(hit->pln==hit2->pln){

	hit3d.clear();
	hit3d2.clear();

	if(hit->mod==16) hit3d.z=hit->pln*4.6+2.3*hit->view+0.4;
	else hit3d.z=107.45+hit->pln*10.7+hit->view+1.0;
	/*	cout<<hit3d.z<<", plane="<<hit3d.pln<<", channel="<<hit3d.ch<<endl;
	if(hit2->mod==16) hit3d.z=hit->z;
	else hit3d.z=hit->z+107.45;
	cout<<"v2="<<hit3d.z<<", plane="<<hit3d.pln<<", channel="<<hit3d.ch<<endl;
	*/
	//if(hit->mod==16 && (hit->ch>7 && hit->ch<24)) cout<<"hit pe="<<hit->pe<<", hit lope="<<hit->lope<<endl;
	bool MC=false;
	if(MC) hit3d.pe=hit->pe;
        else{
          if((hit->pe + hit->lope)/2.<39) hit3d.pe=hit->pe;
          else hit3d.pe=hit->lope;
        }
	hit3d.view=hit->view;
	hit3d.pln=hit->pln;
	hit3d.mod=hit->mod;

	hit3d.time=hit->time;	

	//if(hit->NSimHits()>1)//cout<<"NOOOOOOOOOOOOOOOOOOOOOO"<<endl;
	//if(hit->NSimHits()>1)//cout<<"It shouldn't happen"<<endl;
	if(MC==true && hit->NSimHits()!=0) {
	  hit3d.truepe=((hit->GetIngridSimHit(0))->edeposit)*MeV2PEPM;
	  //cout<<(hit->GetIngridSimHit(0))->pdg<<endl;
	  hit3d.pdg=(hit->GetIngridSimHit(0))->pdg;
	  hit3d.trackid=(hit->GetIngridSimHit(0))->trackid;
	  /*
	  for(int ist=0;ist<evt->NIngridSimParticles();ist++){
	    SimPart = (IngridSimParticleSummary*) evt->GetSimParticle(ist);
	    if(SimPart->trackid==(hit->GetIngridSimHit(0))->trackid) {
	      hit3d.pdg=SimPart->pdg;
	      break;
	    }
	    }*/
	}
	else hit3d.truepe=-1;

	hit3d2.mod=hit2->mod;	
	if(hit2->mod==16) hit3d2.z=hit2->pln*4.6+2.3*hit2->view+0.4;
        else hit3d2.z=107.45+hit2->pln*10.7+hit2->view+1.0;
	//if(hit2->mod==16) hit3d2.z=hit2->z;
	//else hit3d2.z=hit2->z+107.45;

	if(MC) hit3d2.pe=hit2->pe;
	else{
	  if((hit2->pe + hit2->lope)/2.<39) hit3d2.pe=hit2->pe;
	  else hit3d2.pe=hit2->lope;
	}
	hit3d2.pln=hit2->pln;
	hit3d2.view=hit2->view;
	hit3d2.time=hit2->time;
	if(MC && hit2->NSimHits()!=0) {
	  if(hit->NSimHits()>1)//cout<<"It can happen"<<endl;
	  hit3d2.truepe=((hit2->GetIngridSimHit(0))->edeposit)*MeV2PEPM;
	  //cout<<(hit2->GetIngridSimHit(0))->pdg<<endl;
	  hit3d2.pdg=(hit2->GetIngridSimHit(0))->pdg;
	  hit3d2.trackid=(hit2->GetIngridSimHit(0))->trackid;
	  /* for(int ist=0;ist<evt->NIngridSimParticles();ist++){
	    SimPart = (IngridSimParticleSummary*) evt->GetSimParticle(ist);
	    if(SimPart->trackid==(hit2->GetIngridSimHit(0))->trackid) {
	      hit3d2.pdg=SimPart->pdg;
	      break;
	    }
	    }*/
	}

	//if(hit->view==0) {
	if(hit->mod==16){
	  if(hit->ch<=7) hit->xy=5*hit->ch;
	  else if(hit->ch<=23) hit->xy=5*8+2.5*(hit->ch-8);
	  else hit->xy=5*8+2.5*16+5*(hit->ch-24);
	}

	if(hit2->mod==16){
	  if(hit2->ch<=7) hit2->xy=5*hit2->ch;
	  else if(hit2->ch<=23) hit2->xy=5*8+2.5*(hit2->ch-8);
	  else hit2->xy=5*8+2.5*16+5*(hit2->ch-24);
	}

          hit3d.x=hit->xy;
          hit3d.y=hit2->xy;
          hit3d.ch=hit->ch;
          hit3d2.x=hit->xy;
          hit3d2.y=hit2->xy;
          hit3d2.ch=hit2->ch;
	  /*}
        else{
          hit3d.x=hit2->xy;
          hit3d.y=hit->xy;
          hit3d.ch=hit2->ch;
          hit3d2.x=hit->xy;
          hit3d2.y=hit2->xy;
          hit3d.ch=hit->ch;
	  }*/
                                                              
      HitPln[hit3d.pln][hit3d.view][hit3d.ch]++;
      HitPln[hit3d2.pln][hit3d2.view][hit3d2.ch]++;
      VecAll.push_back(hit3d);
      VecAll.push_back(hit3d2);
    }
  }
 }
  sort(VecAll.begin(),VecAll.end());
  /*
  for(int ihit=0;ihit<VecAll.size();ihit++){
    //cout<<"Hit/Pln="<<HitPln[VecAll[ihit].pln][VecAll[ihit].view]<<endl;
    VecAll[ihit].pe/=HitPln[VecAll[ihit].pln][VecAll[ihit].view][VecAll[ihit].ch];
  }
  */
  /* //cout<<"Vec All:"<<endl;
  for(int ihit=0;ihit<VecAll.size();ihit++){
   //cout<<VecAll[ihit].x<<"  "<<VecAll[ihit].y<<"  "<<VecAll[ihit].z<<"  "<<VecAll[ihit].pe<<endl;
    }*/
   return VecAll;
}

vector <Hit3D> Reconstruction::IsInTrk(vector <Hit3D> VecCluster, vector <Hit3D> VecAllTracks){
  for(int ihit=0;ihit<VecCluster.size();ihit++){
    for(int ihit2=0;ihit2<VecAllTracks.size();ihit2++){
      if( VecCluster[ihit].mod==VecAllTracks[ihit2].mod && VecCluster[ihit].pln==VecAllTracks[ihit2].pln && VecCluster[ihit].view==VecAllTracks[ihit2].view && VecCluster[ihit].ch==VecAllTracks[ihit2].ch && VecCluster[ihit].time==VecAllTracks[ihit2].time) VecCluster[ihit].InTrk=true;
    }
  }
  return VecCluster;
}

vector <double> Reconstruction::TrackRelativeLength(double posx, double posy, double posz, double tx, double ty, double TrkLg){//feed it by vertex position and angles in xz and yz planes (in radian), and then, by track length
  double pln=17;
  double ch=31;
  double posiz;
  int view=1;
  
  if(view==0){
    if(pln==0)posiz=5;
    else posiz=46*pln+9;
  }
  else{
    if(pln==0)posiz=28;
    else posiz=46*pln+32;
  }
  double posiz0=5;

  float posixy;

  if(pln==0){
    posixy=ch*50+25;
  }
  else{
    if(ch<8)posixy=ch*50+25;
    else if(ch<24)posixy=412.5+(ch-8)*25;
    else posixy=(ch-8)*50+25;
  }

  //  posixy-=posixy/2;
  double xend=((tx<0)?0:1)*posixy/10;//extremity of Ingrid
  double yend=((ty<0)?0:1)*posixy/10;
  double zendX=posz+(xend-posx)/TMath::Tan(tx);
  double zendY=posz+(yend-posy)/TMath::Tan(ty);
  double zend=(posiz-posiz0)/10;
  //cout<<"New Track"<<endl;
  //cout<<"Module ends= "<<xend<<"   "<<yend<<"   "<<zend<<endl;
  //cout<<"zends= "<<zendX<<"   "<<zendY<<"   "<<zend<<endl;
  //if(zendX<0||zendY<0)//cout<<"Problem!!!!!!!!!!!!!!!"<<endl;
  if(zend>1.3+zendX || zend>zendY){
    if(zendX<zendY) {
      zend=zendX;
      yend=(zendX-posz)*TMath::Tan(ty)+posy;//goes out by the side
    }
    else{
      zend=zendY;
      xend=(zendY-posz)*TMath::Tan(tx)+posx;//goes out by top/bottom
    }
  }
  else {
    xend=TMath::Tan(tx)*(zend-posz)+posx;
    yend=TMath::Tan(ty)*(zend-posz)+posy;
  }

  //cout<<"Vertex Position=("<<posx<<", "<<posy<<", "<<posz<<")  , ThetaX="<<tx*180/TMath::Pi()<<" , ThetaY="<<ty*180/TMath::Pi()<<endl;
  //cout<<"Ending Position=("<<xend<<", "<<yend<<", "<<zend<<")"<<endl;
  /*
  int StartPln=(int) posz/2.3;
  int EndPln=(int) zend/2.3;
  //L in scinti depends if INGRID or SciBar...
  double SlopeX=TMath::Tan(tx);
  double SlopeY=TMath::Tan(ty);
  double dx=TMath::Sqrt(SlopeY*SlopeY+SlopeX*SlopeX+1);
  */


  double GeomLength=TMath::Sqrt(pow(xend-posx,2)+pow(yend-posy,2)+pow(zend-posz,2));
  //cout<<"Ideal Track Length="<<GeomLength<<" ,True Track Length="<<TrkLg<<endl;

  //on a le début et la fin de la trace dans les coordonnées PM: il suffit maintenant de regarder de quels plan s'agit il => 


  vector <double> Length;
  Length.clear();
  Length.push_back(GeomLength);
  Length.push_back(TrkLg/GeomLength);
  return Length;
}


vector <double> Reconstruction::TheoreticalTrack(int startplnx, double startchx, double thetax, int startplny, double startchy, double thetay){//feed it by vertex position and angles in xz and yz planes (in radian), and then, by track length

  double Mod[17][3];
  
  Mod[3][0]=-0.000863*100;
  Mod[3][1]=-17.55557*100;
  Mod[3][2]=277.36844*100;
  
  Mod[10][0]=-0.038371*100;
  Mod[10][1]=-17.37957*100;//to check
  Mod[10][2]=273.35956*100;

  Mod[16][0]=0.001*100;
  Mod[16][1]=-17.563*100;
  Mod[16][2]=276.168*100;
  
  //Center[16][0]=60.15;
  //Center[16][1]=60.15;
  //Center[16][2]=39.1;

  /*conversion towards PM coordinates*/

  Mod[3][0]-=Mod[16][0];
  Mod[3][1]-=Mod[16][1];
  Mod[3][2]-=Mod[16][2];
  Mod[10][0]-=Mod[16][0];
  Mod[10][1]-=Mod[16][1];
  Mod[10][2]-=Mod[16][2];

  for(int imod=0;imod<NMod;imod++){
    if(imod<7){
      Mod[imod][0]=Mod[3][0]+(imod-3)*150;
      Mod[imod][1]=Mod[3][1];
      Mod[imod][2]=Mod[3][2];
    }
    else{
      Mod[imod][0]=Mod[10][0];
      Mod[imod][1]=Mod[10][1]+(imod-10)*150;
      Mod[imod][2]=Mod[10][2];
    }
  }
  
  double X[3];
  int StartPln=0;
  int EndPln=0;
  double StartChX=0; double StartChY=0; double EndChX=0; double EndChY=0;
  int StartPln_Ing=0;
  int EndPln_Ing=0;
  double StartChX_Ing,StartChY_Ing,EndChX_Ing,EndChY_Ing;
  double *Xx=new double();
  double *Zx=new double();
  double*Yy=new double();
  double *Zy=new double();
  double *x_ini=new double();
  double *z_iniX=new double();
  bool done1=IngDimRec->INGRID_Dimension::get_posXY(16,0,startplnx,0,x_ini,z_iniX);
  double *y_ini=new double();
  double *z_iniY=new double();
  bool done2=IngDimRec->INGRID_Dimension::get_posXY(16,1,startplnx,0,y_ini,z_iniY);
  *x_ini=startchx/10;
  *y_ini=startchy/10;
  StartPln=min(startplnx,startplny);
  StartChX=startchx;
  StartChY=startchy;

  for(int ipln=StartPln+1;ipln<NPlnPM;ipln++){
    //cout<<"Plane="<<ipln<<", startChX="<<startchx<<endl;
    bool done3=IngDimRec->INGRID_Dimension::get_posXY(16,0,ipln,0,Xx,Zx);
    X[0]=*x_ini+TMath::Tan(thetax)*(*Zx-*z_iniX);
    bool done4=IngDimRec->INGRID_Dimension::get_posXY(16,1,ipln,0,Yy,Zy);
    X[1]=*y_ini+TMath::Tan(thetay)*(*Zy-*z_iniY);

    if(X[0]<120.3 && X[0]>0 && X[1]<120.3 && X[1]>0) {
      //cout<<"Theoretical is in Module"<<endl;
      EndPln=ipln;
      EndChX=X[0]*10;
      EndChY=X[1]*10;
    }
    else break;
  }
  //cout<<"end of PMM"<<endl;

  bool start=false;
  bool HitMod=false;

  for(int imod=0;imod<7;imod++){
    if(!HitMod){
      if(imod==6) HitMod=true;
      for(int ipln=0;ipln<NPln;ipln++){
	//cout<<"ipln Ingrid="<<ipln<<endl;
	bool done3=IngDimRec->INGRID_Dimension::get_posXY(imod,0,ipln,1,Xx,Zx);
	X[0]=*x_ini+TMath::Tan(thetax)*(*Zx+107.45-*z_iniX);//X in the PM coordinates
	bool done4=IngDimRec->INGRID_Dimension::get_posXY(imod,1,ipln,1,Yy,Zy);
	X[1]=*x_ini+TMath::Tan(thetay)*(*Zy+107.45-*z_iniY);
	//cout<<"Z difference ="<<*Zy+107.45-*z_iniY<<endl;
	//cout<<"jist before"<<endl;	
	X[0]-=60.15+Mod[imod][0];
	X[1]-=60.15+Mod[imod][1];
	//cout<<"X and Y in INGRID="<<X[0]<<", "<<X[1]<<endl;
	if(TMath::Abs(X[0])<60.15 && TMath::Abs(X[1])<60.15) {
	  if(start==false){
	    StartPln_Ing=ipln;
	    StartChX_Ing=X[0]*10;
	    StartChY_Ing=X[1]*10;
	    start=true;
	  }
	  //cout<<"In INGRID"<<endl;
	  EndPln_Ing=ipln;
	  EndChX_Ing=X[0]*10;
	  EndChY_Ing=X[1]*10;
	  HitMod=true;
	  //cout<<"In module="<<imod<<endl;
	}
	//else break;
      }  
    }
  }
  if(StartPln_Ing>EndPln_Ing) cout<<"*********************************************************************"<<endl;
  vector <double> Track;
  //  double T[]={StartPln,StartChX,EndPln,EndChX, StartPln,StartChY,EndPln,EndChY, StartPln_Ing,EndPln_Ing};
  Track.push_back(startplnx);
  Track.push_back(StartChX);
  Track.push_back(thetax);
  Track.push_back(startplny);
  Track.push_back(StartChY);
  Track.push_back(thetay);
  Track.push_back(EndPln);
  Track.push_back(EndChX);
  Track.push_back(EndPln);
  Track.push_back(EndChY);
  Track.push_back(StartPln_Ing);
  Track.push_back(EndPln_Ing);
  ////cout<<StartPln_Ing<<", "<<EndPln_Ing<<endl;
  return Track;
}



/*Center of modules*/

//on va rentrer le intcptx,intpcty et thetax, thetay
//sortir si on a un module traversé avec une chance de créer une trace: renvoyer un booléen: oui et non

bool Reconstruction::HasGeomTrack(int mod, int startplnx, int startchx, double thetax,int startplny, int startchy,double thetay){
  double Mod[NMod][3];
  double Scinti[NCh][2];
  double *posXxy=new double();
  double *posXz=new double();
  double *posYxy=new double();
  double *posYz=new double();
  startchx/=10;
  startchy/=10;
  bool HasGeom=false;

  //cout<<endl;
  //cout<<"Mod="<<mod<<" , startchx="<<startchx<<" , startplnZ="<<startplnx<<endl;
  bool doneX=IngDimRec->INGRID_Dimension::get_posXY(mod,0,startplnx,0,posXxy,posXz);
  bool doneY=IngDimRec->INGRID_Dimension::get_posXY(mod,1,startplny,0,posYxy,posYz);

  //cout<<"Starting Point: PosX="<<startchx<<" , PosY="<<startchy<<" , PosZ="<<*posXz<<endl;
  //cout<<"Track Angle X="<<thetax*180/TMath::Pi()<<" , Angle Y="<<thetay*180/TMath::Pi()<<endl;
  Mod[3][0]=-0.000863*100;
  Mod[3][1]=-17.55557*100;
  Mod[3][2]=277.36844*100;
  
  Mod[10][0]=-0.038371*100;
  Mod[10][1]=-17.37957*100;//to check
  Mod[10][2]=273.35956*100;

  Mod[16][0]=0.001*100;
  Mod[16][1]=-17.563*100;
  Mod[16][2]=276.168*100;
  
  //Center[16][0]=60.15;
  //Center[16][1]=60.15;
  //Center[16][2]=39.1;

  /*conversion towards PM coordinates*/
  
  Mod[3][0]-=Mod[16][0];
  Mod[3][1]-=Mod[16][1];
  Mod[3][2]-=Mod[16][2];
  Mod[10][0]-=Mod[16][0];
  Mod[10][1]-=Mod[16][1];
  Mod[10][2]-=Mod[16][2];

  for(int imod=0;imod<NMod;imod++){
    if(imod<7){
      Mod[imod][0]=Mod[3][0]+(imod-3)*150;
      Mod[imod][1]=Mod[3][1];
      Mod[imod][2]=Mod[3][2];
    }
    else{
      Mod[imod][0]=Mod[10][0];
      Mod[imod][1]=Mod[10][1]+(imod-10)*150;
      Mod[imod][2]=Mod[10][2];
    }
  }
  
  //  pos[0]=posX - 60.15 + Mod[16][0]*100;
  //pos[1]=posY - 60.15 + Mod[16][1]*100;
  //pos[2]=posZ - 39.1 + Mod[16][2]*100;
  double PlnZx, PlnZy;
  double X[2];
  int SucPlane[NView][NMod];
  bool Successive[NView][NMod];
  //cout<<"Plane="<<ipln<<", Z pos="<<PlnZ<<endl;
  //cout<<endl;
  //cout<<"Position of Module 3 in PM coordinates="<<Mod[3][0]<<", "<<Mod[3][1]<<", "<<Mod[3][2]<<endl;
  PlnZx=PlnThick_PM*17;
  //cout<<"In the end of PM, here are positions="<<startchx+TMath::Tan(thetax)*(PlnZx-*posXz)<<", "<<startchy+TMath::Tan(thetay)*(PlnZx-*posYz)<<", "<<PlnZx<<endl; 
  double PlnZIng=(PlnThick+IronThick)*1+107.45; 
  //cout<<"At Ingrid Start, positions="<<startchx+TMath::Tan(thetax)*(PlnZIng-*posXz)<<", "<<startchy+TMath::Tan(thetay)*(PlnZIng-*posYz)<<", "<<PlnZIng<<endl;
  for(int imod=0;imod<7;imod++){
    //cout<<"Module"<<imod<<endl;
      for(int iview=0;iview<NView;iview++){//on va regarder pour le smodules horizontaux, les plans successifs. Puis on regarde si pour le z alors fixé, dans quel module est le Xs'il est dans un module. On regarde alors dans quel module est le y en ce meme point pour voir si on va déposer en XZ. Puis faire de meme pour les plans pour view=1: dans celui la, on ne peut savoir dans quel module on est avec le seul y=> si on a le x donné par la droite
	//cout<<"View="<<iview<<endl;
	Successive[iview][imod]=false;
	//if(iview==0)//cout<<"XZ Plane:"<<endl;
	//else//cout<<"YZ Plane:"<<endl;
	SucPlane[iview][imod]=0;
	//cout<<"difference in position="<<Mod[3][2]<<endl; 
	for(int ipln=0;ipln<NPln;ipln++){
	  //evaluer la droite XZ en Z=scinti vert, puis la droite YZ
	  //On obtient alors X,Z => Tester si on est dans un scinti d'INGRID
	  //Pour cela: on a scinti[j]: parcourir les j et voir si on est bien
	  //d'abord faudrait avoir le module après avoir testé les z =>
	  
	    //cout<<"XZ Plane:"<<endl;
	  PlnZx=(PlnThick+IronThick)*ipln+107.45+1;//(Mod[3][2]-Mod[16][2]);//différencier si view =0 ou 1. C'est donc dans le system de coordonnées du PM là (de son début d'ailleurs) 
	  X[0]=startchx+TMath::Tan(thetax)*(PlnZx-*posXz);
	  //cout<<"YZ Plane:"<<endl;
	  PlnZy=(PlnThick+IronThick)*ipln+107.45;//(Mod[3][2]-Mod[16][2])+1.;//différencier si view =0 ou 1
	  X[1]=startchy+TMath::Tan(thetay)*(PlnZy-*posYz);
	  //cout<<"X0="<<X[0]<<" , X1="<<X[1]<<endl;
	  
	  //	 //cout<<"Module Coordinates in PM referential=";
	  X[0]-=60.15+Mod[imod][0];
	  X[1]-=60.15+Mod[imod][1];
	  if(TMath::Abs(X[0])<60.15 && TMath::Abs(X[1])<60.15) {
	    //cout<<"In Module"<<imod<<endl;
	    //cout<<"Plane="<<ipln<<", X0="<<X[0]<<" , X1="<<X[1]<<endl;
	    for(int ich=0;ich<NCh;ich++){
	      Scinti[ich][1]=(ScintiWidth)*ich+0.743-60.15;//position to the center of module
	      Scinti[ich][0]=(ScintiWidth)*ich-0.186-60.15;
	      //cout<<"X="<<X[iview]<<" , Scvinti="<<Scinti[ich][iview]<<" ,Difference="<<TMath::Abs(X[iview]-Mod[imod][iview]-Scinti[ich][iview])<<endl;
	      if(TMath::Abs(X[iview]-Mod[imod][iview]-Scinti[ich][iview])<2.5){
		SucPlane[iview][imod]++;
		if(SucPlane[iview][imod]>=3) Successive[iview][imod]=true;
		//cout<<"it's in Plane="<<ipln<<" ,in scintillators="<<ich<<" and in Module="<<imod<<" , SucPlane="<<SucPlane[iview][imod]<<" , Bool Succ="<<Successive[iview][imod]<<endl;
		break;
	      }
	    }//ch
	  }
	}//mod
      }//view
      //droite XZ: a chaque scinti[i]
  }//planes
  
  for(int imod=0;imod<7;imod++){
    if(Successive[0][imod] && Successive[1][imod]) {
      HasGeom=true;
      //cout<<"Has Geom in Module="<<imod<<endl;
    }
  }

  return HasGeom;
}


vector <double> Reconstruction::IngridTrack(int mod, int startplnx, int startchx, double thetax,int startplny, int startchy,double thetay){
  double Mod[NMod][3];
  double Scinti[NCh][2];
  double *posXxy=new double();
  double *posXz=new double();
  double *posYxy=new double();
  double *posYz=new double();
  startchx/=10;
  startchy/=10;
  bool HasGeom=false;

  //cout<<endl;
  //cout<<"Mod="<<mod<<" , startchx="<<startchx<<" , startplnZ="<<startplnx<<endl;
  bool doneX=IngDimRec->INGRID_Dimension::get_posXY(mod,0,startplnx,0,posXxy,posXz);
  bool doneY=IngDimRec->INGRID_Dimension::get_posXY(mod,1,startplny,0,posYxy,posYz);

  //cout<<"Starting Point: PosX="<<startchx<<" , PosY="<<startchy<<" , PosZ="<<*posXz<<endl;
  //cout<<"Track Angle X="<<thetax*180/TMath::Pi()<<" , Angle Y="<<thetay*180/TMath::Pi()<<endl;
  Mod[3][0]=-0.000863*100;
  Mod[3][1]=-17.55557*100;
  Mod[3][2]=273.35956*100;
  
  Mod[10][0]=-0.038371*100;
  Mod[10][1]=-17.37957*100;//to check
  Mod[10][2]=273.35956*100;

  Mod[16][0]=0.001*100;
  Mod[16][1]=-17.563*100;
  Mod[16][2]=276.168*100;
  
  //Center[16][0]=60.15;
  //Center[16][1]=60.15;
  //Center[16][2]=39.1;

  /*conversion towards PM coordinates*/
  
  Mod[3][0]-=Mod[16][0];
  Mod[3][1]-=Mod[16][1];
  Mod[3][2]-=Mod[16][2];
  Mod[10][0]-=Mod[16][0];
  Mod[10][1]-=Mod[16][1];
  Mod[10][2]-=Mod[16][2];

  for(int imod=0;imod<NMod;imod++){
    if(imod<7){
      Mod[imod][0]=Mod[3][0]+(imod-3)*150;
      Mod[imod][1]=Mod[3][1];
      Mod[imod][2]=Mod[3][2];
    }
    else{
      Mod[imod][0]=Mod[10][0];
      Mod[imod][1]=Mod[10][1]+(imod-10)*150;
      Mod[imod][2]=Mod[10][2];
    }
  }
  
  //  pos[0]=posX - 60.15 + Mod[16][0]*100;
  //pos[1]=posY - 60.15 + Mod[16][1]*100;
  //pos[2]=posZ - 39.1 + Mod[16][2]*100;
  double PlnZ;
  double X[2];
  //cout<<"Plane="<<ipln<<", Z pos="<<PlnZ<<endl;
  //cout<<endl;
  int StartPln=100;
  int EndPln=0;
  for(int imod=0;imod<7;imod++){
      for(int iview=0;iview<NView;iview++){//on va regarder pour le smodules horizontaux, les plans successifs. Puis on regarde si pour le z alors fixé, dans quel module est le Xs'il est dans un module. On regarde alors dans quel module est le y en ce meme point pour voir si on va déposer en XZ. Puis faire de meme pour les plans pour view=1: dans celui la, on ne peut savoir dans quel module on est avec le seul y=> si on a le x donné par la droite
	for(int ipln=0;ipln<NPln;ipln++){
	  //evaluer la droite XZ en Z=scinti vert, puis la droite YZ
	  //On obtient alors X,Z => Tester si on est dans un scinti d'INGRID
	  //Pour cela: on a scinti[j]: parcourir les j et voir si on est bien
	  //d'abord faudrait avoir le module après avoir testé les z =>
	  
	  if(iview==0) {
	    //cout<<"XZ Plane:"<<endl;
	    PlnZ=(PlnThick+IronThick)*ipln+120.044;//différencier si view =0 ou 1
	    X[0]=startchx+TMath::Tan(thetax)*(PlnZ-*posXz);
	  }
	  else{
	    //cout<<"YZ Plane:"<<endl;
	    PlnZ=(PlnThick+IronThick)*ipln+120.044+1.;//différencier si view =0 ou 1
	    X[1]=startchy+TMath::Tan(thetay)*(PlnZ-*posYz);
	  }
	  
	  if(TMath::Abs(X[0]-Mod[imod][0])<60.15 && TMath::Abs(X[1]-Mod[imod][1])<60.15) {
	    for(int ich=0;ich<NCh;ich++){
	      Scinti[ich][1]=(ScintiWidth)*ich+0.743-60.15;
	      Scinti[ich][0]=(ScintiWidth)*ich-0.186-60.15;
	      if(TMath::Abs(X[1]-Mod[imod][1]-Scinti[ich][1])<2.5){
		if(ipln<StartPln) StartPln=ipln;
		if(ipln>EndPln) EndPln=ipln;
	      }
	    }//ch
	  }
	}//mod
      }//view
      //droite XZ: a chaque scinti[i]
  }//planes

 //cout<<"Has Ingrid Track="<<StartPln<<", Ending="<<EndPln<<endl;  
  vector <double> Pln;
  Pln.push_back(StartPln);
  Pln.push_back(EndPln);
  return Pln;
}

vector <double> Reconstruction::TrackPenetration(int Mod, int pln_iniX, double ch_iniX, double thetax,int pln_iniY, double ch_iniY, double thetay, int pln_finX, double ch_finX, int pln_finY, double ch_finY, double dx_Ing){
  //1plane in PM: careful if the channel is SciBar or Ingrid. 1plane = 2channels to cross. Look if the channel is a scibar or an Ingrid one?
  //for normal track: enter zini and zfinal? No plane ini and plane final:
  //dist=(plane final - plane ini + 1)*dxofchanneltype
  //case1: Ing scinti => how to know it? angle of tracks is known => identify channels it passes in. For this purpose, separate XZ and YZ:
  //enter pln_iniX,pln_iniY,pln_finX,pln_finY
  //dx_Ing, dx_Sci
  ch_iniX/=10;
  ch_finX/=10;
  ch_iniY/=10;
  ch_finY/=10;
  double dx_Iron=dx_Ing*6.5;  

  double *x_ini=new double();
  double *z_iniX=new double();
  bool done1=IngDimRec->INGRID_Dimension::get_posXY(Mod,0,pln_iniX,ch_iniX,x_ini,z_iniX);
  double *y_ini=new double();
  double *z_iniY=new double();
  bool done2=IngDimRec->INGRID_Dimension::get_posXY(Mod,1,pln_iniY,ch_iniY,y_ini,z_iniY);
  //cout<<"In penetration fonction"<<endl;
  double expx,expy,expzX,expzY;
  int xch,ych;
  double DistCarbon=0;
  double DistIron=0;


  vector <double> Dist;
  ////cout<<"Pln Ini="<<pln_iniX<<", Pln final="<<pln_finX<<endl;
  for(int ipln=pln_iniX;ipln<=pln_finX;ipln++){
    expzX=(ipln-pln_iniX)*PlnThick+*z_iniX;
    expx=ch_iniX+((thetax>0)?1:-1)*expzX*TMath::Tan(thetax);//not forgot before to convert thetax in deg->rad
    for(int numch=0;numch<48;numch++){
      double diffxy=expx-numch*ScintiWidth;
      if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	xch=numch;
      }
    }
    DistCarbon+=dx_Ing;//ok now the question: IngDim is in which coordinate system? the Module one!
    if(ipln<=9) DistIron+=dx_Iron;
    
  }
  //cout<<"X is over"<<endl;

  for(int ipln=pln_iniY;ipln<=pln_finY;ipln++){
    expzY=(ipln-pln_iniY)*PlnThick_PM+*z_iniY;
    expy=ch_iniY+((thetay>0)?1:-1)*expzY*TMath::Tan(thetay);//not forgot before to convert thetax in deg->rad
    for(int numch=0;numch<48;numch++){
      double diffxy=expx-numch*ScintiWidth;
      if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	ych=numch;
      }
    }
    DistCarbon+=dx_Ing;//ok now the question: IngDim is in which coordinate system? the Module one!
  }
  Dist.push_back(DistCarbon);
  Dist.push_back(DistIron);
  //cout<<"Distance in Carbon="<<DistCarbon<<endl;
  //cout<<"Distance in Iron="<<DistIron<<endl;
  double EqSciLength=DistIron*55.85/1.03+DistCarbon;//equivalent length of scintillators crossed by the particle
  Dist.push_back(EqSciLength);

  return Dist;
}


vector <double> Reconstruction::TrackPenetrationPM(int pln_iniX, double ch_iniX, double thetax,int pln_iniY, double ch_iniY, double thetay, int pln_finX, double ch_finX, int pln_finY, double ch_finY, int IngMod, int pln_ini_Ing, int pln_fin_Ing, double dx_Ing, bool PMStop){
  //1plane in PM: careful if the channel is SciBar or Ingrid. 1plane = 2channels to cross. Look if the channel is a scibar or an Ingrid one?
  //for normal track: enter zini and zfinal? No plane ini and plane final:
  //dist=(plane final - plane ini + 1)*dxofchanneltype
  //case1: Ing scinti => how to know it? angle of tracks is known => identify channels it passes in. For this purpose, separate XZ and YZ:
  //enter pln_iniX,pln_iniY,pln_finX,pln_finY
  //dx_Ing, dx_Sci
  ch_iniX/=10;
  ch_finX/=10;
  ch_iniY/=10;
  ch_finY/=10;
  double dx_Sci=dx_Ing*1.3;
  double dx_Iron=dx_Ing*6.5;  

  double *x_ini=new double();
  double *z_iniX=new double();
  bool done1=IngDimRec->INGRID_Dimension::get_posXY(16,0,pln_iniX,ch_iniX,x_ini,z_iniX);
  double *y_ini=new double();
  double *z_iniY=new double();
  bool done2=IngDimRec->INGRID_Dimension::get_posXY(16,1,pln_iniY,ch_iniY,y_ini,z_iniY);
  //cout<<"In penetration fonction"<<endl;
  double expx,expy,expzX,expzY;
  int xch,ych;
  double DistCarbon=0;
  vector <double> Dist;
  ////cout<<"Pln Ini="<<pln_iniX<<", Pln final="<<pln_finX<<endl;
  for(int ipln=pln_iniX;ipln<pln_finX;ipln++){
    expzX=(ipln-pln_iniX)*PlnThick_PM+*z_iniX;
    expx=ch_iniX+((thetax>0)?1:-1)*expzX*TMath::Tan(thetax);//not forgot before to convert thetax in deg->rad
    for(int numch=0;numch<48;numch++){
      double diffxy=expx-numch*ScintiWidth;
      if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	xch=numch;
      }
    }
    if(xch<8||xch>23) DistCarbon+=dx_Ing;//ok now the question: IngDim is in which coordinate system? the Module one!
    else DistCarbon+=dx_Sci;
  }
  //cout<<"X is over"<<endl;

  for(int ipln=pln_iniY;ipln<pln_finY;ipln++){
    expzY=(ipln-pln_iniY)*PlnThick_PM+*z_iniY;
    expy=ch_iniY+((thetay>0)?1:-1)*expzY*TMath::Tan(thetay);//not forgot before to convert thetax in deg->rad
    for(int numch=0;numch<48;numch++){
      double diffxy=expx-numch*ScintiWidth;
      if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	ych=numch;
      }
    }
    if(xch<8||xch>23) DistCarbon+=dx_Ing;//ok now the question: IngDim is in which coordinate system? the Module one!
    else DistCarbon+=dx_Sci;
  }
  //cout<<"Y is over"<<endl;
  //cout<<"DistCarbon before Ingrid="<<DistCarbon<<endl;
  double DistIron=0;
  if(!PMStop){
    DistCarbon+=dx_Ing*(pln_fin_Ing-pln_ini_Ing+1);
    if(pln_fin_Ing<=9) DistIron=dx_Iron*(pln_fin_Ing-pln_ini_Ing);
    else DistIron=dx_Iron*9;
  }
  //if(DistIron!=DistIron)//cout<<"************************************************************************/"<<endl<<"dx="<<dx_Ing<<" , Plane ="<<pln_fin_Ing<<endl;
  
  Dist.push_back(DistCarbon);
  Dist.push_back(DistIron);
  //cout<<"Distance in Carbon="<<DistCarbon<<endl;
  //cout<<"Distance in Iron="<<DistIron<<endl;
  double EqSciLength=DistIron*55.85/1.03+DistCarbon;//equivalent length of scintillators crossed by the particle
  Dist.push_back(EqSciLength);

  return Dist;
}

vector <double> Reconstruction::ConvertTruePM(float ipos[4], float fpos[4]){
    ipos[0]+=120.3/2;
    ipos[1]+=120.3/2;
    ipos[2]+=120.00+40.95+1.55;

    fpos[0]+=120.3/2;
    fpos[1]+=120.3/2;
    fpos[2]+=120.00+40.95+1.55;


   //cout<<"ipos="<<ipos[0]<<", "<<ipos[1]<<", "<<ipos[2]<<endl;
   //cout<<"fpos="<<fpos[0]<<", "<<fpos[1]<<", "<<fpos[2]<<endl;
    int startxpln=(int) ipos[2]/4.6;
    int startypln=(int) ipos[2]/4.6;
    int startxch, startych;
    if(ipos[0] < 8*5 +2.5) startxch=(int) ipos[0]/5;
    else if(ipos[0] < 8*5 + 16*2.5 + 2.5) startxch=(int) 8 + (ipos[0]-8*5)/2.5;
    else startxch=(int) 16 + 8 + (ipos[0]-8*5-16*2.5)/5;
    if(ipos[1] < 8*5 +2.5) startych=(int) ipos[1]/5;
    else if(ipos[1] < 8*5 + 16*2.5 + 2.5) startych=(int) 8 + (ipos[1]-8*5)/2.5;
    else startych=(int) 16 + 8 + (ipos[1]-8*5-16*2.5)/5;

   //cout<<"Startxpln="<<startxpln<<", Startypln="<<startypln<<", Startxch="<<startxch<<", Startych="<<startych<<endl;
    double thetaX=TMath::ATan((fpos[0]-ipos[0])/(fpos[2]-ipos[2]));
    double thetaY=TMath::ATan((fpos[1]-ipos[1])/(fpos[2]-ipos[2]));



    double posiz=4.6*NPlnPM+3.2;//careful: I assumed that XZ and YZ scintillators planes are at the same distance. I shouldn't!!!!! CORRECT THIS! 
    double xend=((thetaX<0)?0:1)*(16*2.5+16*5);//extremity of Ingrid                                                                                      
    double yend=((thetaY<0)?0:1)*(16*2.5+16*5);
    double zendX=ipos[2]+(xend-ipos[0])/TMath::Tan(thetaX);
    double zendY=ipos[2]+(yend-ipos[1])/TMath::Tan(thetaY);
    double zend=posiz;//careful different from in relativetrack function
    //cout<<"New Track"<<endl;                                                                                                                           
    //cout<<"Module ends= "<<xend<<"   "<<yend<<"   "<<zend<<endl;                                                                                       
    //cout<<"zends= "<<zendX<<"   "<<zendY<<"   "<<zend<<endl;                                                                                           
    //if(zendX<0||zendY<0)//cout<<"Problem!!!!!!!!!!!!!!!"<<endl;                                                                                         
    if(zend>1.3+zendX || zend>zendY){
      if(zendX<zendY) {
	zend=zendX;
	yend=(zendX-ipos[2])*TMath::Tan(thetaY)+ipos[1];//goes out by the side                                                                                     
      }
      else{
	zend=zendY;
	xend=(zendY-ipos[2])*TMath::Tan(thetaX)+ipos[0];//goes out by top/bottom                                                                                   
      }
    }
    else {
      xend=TMath::Tan(thetaX)*(zend-ipos[2])+ipos[0];
      yend=TMath::Tan(thetaY)*(zend-ipos[2])+ipos[1];
    }

   //cout<<"xend="<<xend<<", yend="<<yend<<", zend="<<zend<<endl;
    int endxpln;
    int endypln;
    int endxch, endych;
    bool StopPM=false;

    if(fpos[2]<zend){
      zend=fpos[2];
      xend=fpos[0];
      yend=fpos[1];
      StopPM=true;
    }      

      if(ipos[0] < 8*5 +2.5) endxch=(int) ipos[0]/5;
      else if(ipos[0] < 8*5 + 16*2.5 + 2.5) endxch=(int) 8 + (ipos[0]-8*5)/2.5;
      else endxch=(int) 16 + 8 + (ipos[0]-8*5-16*2.5)/5;
      if(ipos[1] < 8*5 +2.5) endych=(int) ipos[1]/5;
      else if(ipos[1] < 8*5 + 16*2.5 + 2.5) endych=(int) 8 + (ipos[1]-8*5)/2.5;
      else endych=(int) 16 + 8 + (ipos[1]-8*5-16*2.5)/5;
    
    endxpln=(int) zend/4.6;
    endypln=(int) zend/4.6;
   //cout<<"Endxpln="<<endxpln<<", Endypln="<<endypln<<", Endxch="<<endxch<<", Endych="<<endych<<endl;      
    //cout<<"StartPln="<<startxpln<<", StartCh="<<startxch<<", ending plane="<<endxpln<<endl;
    //idea: pass in each Ingrid Module
    //or apply the theoretical track construction on this point! => and then
    //Use HasGeomTrack and getout the number of the module, and possible first and last plane hit
    int ipln_ing,fpln_ing;
    double xIngIni, yIngIni, xIngFin, yIngFin, zIngIni, zIngFin;
    int IngMod;
    bool CrossINGRID=false;
    if((fpos[2]-(120.00+40.95+1.55)+54.5)>0){//Has INGRID Track
      fpos[2]-=(120.+40.95+1.55-54.5);
      ipos[2]-=(120.+40.95+1.55-54.5);
      if(fpos[2]>82) {
	CrossINGRID=true;
	//In this case, evaluate position at start and ending in z of Ingrid modules => if not in a module, iteration of planes to find where it enters and where it goes out
      }
     //cout<<"CrossIngrid="<<CrossINGRID<<endl;
      for(int ipln=0;ipln<NPln;ipln++){
	double ztest=.5+10.5*ipln;//vérifier .5
	xIngIni=(ztest-ipos[2])*TMath::Tan(thetaX)+ipos[0];
	yIngIni=(ztest-ipos[2])*TMath::Tan(thetaY)+ipos[1];
	for(int imod=0;imod<7;imod++){
	  if(TMath::Abs(xIngIni-(60.15+(imod-3)*150))<60.15 && TMath::Abs(yIngIni-60.15)<60.15){
	  IngMod=imod;
	  zIngIni=ztest;
	  ipln_ing=ipln;
	  //put a goto
	  goto EndLoop;
	  }
	}
      }
      //goto
    EndLoop:
     //cout<<"Ingrid Module="<<IngMod<<endl;
      //now test that in the same module, from the zIngIni to the end, which is the last plane
      //cases where it remains in the same module
      bool StayInModule=false;
      if(!CrossINGRID){
	if(TMath::Abs(fpos[0]-(60.15+(IngMod-3)*150))<60.15 && TMath::Abs(fpos[1]-60.15)<60.15){
	  fpln_ing=(int) (fpos[2]-.5)/10.5;
	  StayInModule=true;
	}
      }
      else{
	double zfin=.5+10.5*10;
	xIngFin=(zfin-ipos[2])*TMath::Tan(thetaX)+ipos[0];
        yIngFin=(zfin-ipos[2])*TMath::Tan(thetaY)+ipos[1];
	if(TMath::Abs(xIngFin-(60.15+(IngMod-3)*150))<60.15 && TMath::Abs(yIngFin-60.15)<60.15){
	  fpln_ing=(int) (fpos[2]-.5)/10.5;
	  StayInModule=true;
	}
      }
      
      if(!StayInModule){
	//determine where it leaves and estimate number of planes crossed
	for(int ipln=ipln_ing+1;ipln<NPln;ipln++){
	  double ztest=.5+10.5*ipln;                                                                                                                          
	    xIngFin=(ztest-ipos[2])*TMath::Tan(thetaX)+ipos[0];
	    yIngFin=(ztest-ipos[2])*TMath::Tan(thetaY)+ipos[1];
	    if(TMath::Abs(xIngIni-(60.15+(IngMod-3)*150))<60.15 && TMath::Abs(yIngIni-60.15)<60.15){
		IngMod=IngMod;
		fpln_ing=ipln;
		continue;
	    }
	}
      }
     //cout<<"startpln="<<ipln_ing<<", Ending plane="<<fpln_ing<<endl;
      //case 1: it ends before last plane, and in the same module: fpos<=>endpln, endch....
      //case 2: it goes out the back of the module or the edge without entering another module
      //case 3: it enters another module => take the first module hitted for now
      
      //evaluate at fpos[2]=0 if we are in a module => if yes, check if we are also finishing in a module
    }

    double dx=TMath::Sqrt(TMath::Tan(thetaX)*TMath::Tan(thetaX)+TMath::Tan(thetaY)*TMath::Tan(thetaY)+1);

    vector <double> Out;
    Out.push_back(startxpln);
    Out.push_back(startxch);
    Out.push_back(thetaX);
    Out.push_back(startypln);
    Out.push_back(startych);
    Out.push_back(thetaY);
    Out.push_back(endxpln);
    Out.push_back(endxch);
    Out.push_back(endypln);
    Out.push_back(endych);
    Out.push_back(IngMod);
    Out.push_back(ipln_ing);
    Out.push_back(fpln_ing);
    Out.push_back(dx);
    Out.push_back(StopPM);
    //needed informations: startxpln/y, startchxy, angleX, angleY, +ending planex/y, endchxy + Stop PM + Has Ingrid Track + INGRID informations(mod, start/endpln)
    return Out;
}

double Reconstruction::GetFSI(IngridEventSummary * evt){
  int FSIInt=0;
  int FSIMuons=0;int FSIPions=0;int FSIProtons=0;int FSIPions0=0;
  IngridSimParticleSummary * SimPart2;
  
  for(int is=0;is<evt->NIngridSimParticles();is++){
    SimPart2=(IngridSimParticleSummary*) evt->GetSimParticle(is);

    if(TMath::Abs(SimPart2->pdg)==211) FSIPions++;
    else if(TMath::Abs(SimPart2->pdg)==2212) FSIProtons++;
    else if(TMath::Abs(SimPart2->pdg)==13) FSIMuons++;
    else if(TMath::Abs(SimPart2->pdg)==111) FSIPions0++;
  }

  if(FSIMuons==1){
    if(FSIPions==0){
      if(FSIProtons==0 ) FSIInt=1;
      else if(FSIProtons==1) FSIInt=1;
      else if(FSIProtons>1) FSIInt=2;
    }
    else if(FSIPions==1) FSIInt=3;
    else if(FSIPions>1) FSIInt=4;
  }
  else FSIInt=5;
  
  return FSIInt;
}

bool Reconstruction::InPMFV(IngridSimVertexSummary * simver){

  int mod=simver->mod;
  double posX;
  double posY;
  double posZ;
  posX=simver->xnu;
  posY=simver->ynu;
  posZ=simver->znu;
  //cout<<"posX="<<posX<<", posY="<<posY<<", posZ="<<posZ<<endl;
  if((mod==16) && (fabs(posX)<=50) && (fabs(posY)<=50) && posZ>=-156.4 && posZ<-85 /*posZ>-152 && (posZ<-87.5)*/) return true;
  else return false;
}

vector <double> Reconstruction::GetTrueMuonInformation(IngridEventSummary * evt){
  vector <double> MuonProp;
  MuonProp.clear();
  IngridSimParticleSummary * SimPart2;
  double MuonMom;
  double TrueAngleMuon;
  double TrueMomentumMuon;
  double MuonAngle;

  for(int is=0;is<evt->NIngridSimParticles();is++){
    SimPart2=(IngridSimParticleSummary*) evt->GetSimParticle(is);   
    if(SimPart2->pdg==13){
      MuonMom=TMath::Sqrt(SimPart2->momentum[0]*SimPart2->momentum[0]+SimPart2->momentum[1]*SimPart2->momentum[1]+SimPart2->momentum[2]*SimPart2->momentum[2]);
      double thetaX=TMath::ATan((SimPart2->fpos[0]-SimPart2->ipos[0])/(SimPart2->fpos[2]-SimPart2->ipos[2]));
      double thetaY=TMath::ATan((SimPart2->fpos[1]-SimPart2->ipos[1])/(SimPart2->fpos[2]-SimPart2->ipos[2]));
      TrueAngleMuon=TMath::ACos(1/(pow(TMath::Tan(thetaX),2)+pow(TMath::Tan(thetaY),2)+1))*180/TMath::Pi();
      TrueMomentumMuon=MuonMom;
      MuonAngle=TMath::ACos(1/(pow(TMath::Tan(thetaX),2)+pow(TMath::Tan(thetaY),2)+1))*180/TMath::Pi();
      //double Scalar=(SimPart2->fpos[2]-SimPart2->ipos[2])*1;
      double Scalar=(SimPart2->fpos[0]-SimPart2->ipos[0])*Beam[0]+(SimPart2->fpos[1]-SimPart2->ipos[1])*Beam[1]+(SimPart2->fpos[2]-SimPart2->ipos[2])*Beam[2];
  // double Scalar=(Fpos[0]-Ipos[0])*Beam[0])+(Fpos[1]-Ipos[1])*Beam[1])+(Fpos[2]-Ipos[2])*Beam[2];
      double Norm=TMath::Sqrt((SimPart2->fpos[0]-SimPart2->ipos[0])*(SimPart2->fpos[0]-SimPart2->ipos[0])+(SimPart2->fpos[1]-SimPart2->ipos[1])*(SimPart2->fpos[1]-SimPart2->ipos[1])+(SimPart2->fpos[2]-SimPart2->ipos[2])*(SimPart2->fpos[2]-SimPart2->ipos[2]));
      TrueAngleMuon=TMath::ACos(Scalar/Norm)*180/TMath::Pi();
      //cout<<MuonAngle<<", Angle Test="<<AngleTest<<endl;
    }
  }
  MuonProp.push_back(TrueMomentumMuon);
  MuonProp.push_back(TrueAngleMuon);
  //delete SimPart2;
  
  return MuonProp;
}
 
bool Reconstruction::IsFV(int mod, double posx, double posy, double posz){
  bool IsFV=false;
  //if((mod==16) && (fabs(posx)<=50) && (fabs(posy)<=50) && posz>-152 && (posz<-87.5)) IsFV=true;
  if((mod==16) && (fabs(posx)<=50) && (fabs(posy)<=50) && (posz>=-156.4) && (posz<-85)) return true;
  return IsFV;
}

int Reconstruction::GetTrackParticle(IngridEventSummary *evt, PMAnaSummary * recon, int itrk, double TrkLength){
  int PartNum=0;
  int NSimPart=evt->NIngridSimParticles();
  double Min=180;
  int Particle=0;
  vector <int> SimList;
  SimList.clear();
  IngridSimParticleSummary * SimPart2;
  
  for(int is=0;is<evt->NIngridSimParticles();is++){
    SimPart2=(IngridSimParticleSummary*) evt->GetSimParticle(is);
    if(SimPart2->length<8.6) continue;
    if(SimPart2->pdg==2112) continue;
    double thetaX=(TMath::ATan((SimPart2->fpos[0]-SimPart2->ipos[0])/(SimPart2->fpos[2]-SimPart2->ipos[2])))*TMath::Pi()/180.;
    double thetaY=(TMath::ATan((SimPart2->fpos[1]-SimPart2->ipos[1])/(SimPart2->fpos[2]-SimPart2->ipos[2])))*TMath::Pi()/180.;
    
    bool Used=false;
    for(int i=0;i<SimList.size();i++){
      if(is==SimList[i]) Used=true;
    }
    if(Used) continue;
    
    double Tx=thetaY-(recon->thetax)[itrk];
    double Ty=thetaX-(recon->thetay)[itrk];
    //double dx_Temp=1+pow(TMath::Tan(180*(Tx)/TMath::Pi()),2)+pow(TMath::Tan(180*(Ty)/TMath::Pi()),2);
    double Scalar=(SimPart2->fpos[2]-SimPart2->ipos[2])*1;
    double Norm=TMath::Sqrt((SimPart2->fpos[0]-SimPart2->ipos[0])*(SimPart2->fpos[0]-SimPart2->ipos[0])+(SimPart2->fpos[1]-SimPart2->ipos[1])*(SimPart2->fpos[1]-SimPart2->ipos[1])+(SimPart2->fpos[2]-SimPart2->ipos[2])*(SimPart2->fpos[2]-SimPart2->ipos[2]));
    double Ang3D=TMath::ACos(Scalar/Norm)*180/TMath::Pi()-(recon->angle)[itrk];
    //double Ang3D=TMath::ACos(1/dx_Temp);
    if(TMath::Abs(Ang3D)<Min && TrkLength<2*SimPart2->length){
      

#ifdef DEBUG
      cout<<"Is selected: Particle="<<SimPart2->pdg<<"       :";
      cout<<"here is the rec thetax, thetay="<<(recon->thetay)[itrk]<<", "<<(recon->thetax)[itrk]<<", Angle 3D="<<(recon->angle)[itrk]<<", And here are the angles="<<thetaX<<", "<<thetaY<<", 3D="<<TMath::ACos(Scalar/Norm)*180/TMath::Pi()<<endl;
      cout<<"rec trk length="<<TrkLength<<", Simpart length="<<SimPart2->length<<endl;
#endif
      
      bool Change=false;
      for(int itrk2=0;itrk2<recon->Ntrack;itrk2++){
	if(itrk2==itrk) continue;
	if(SimPart2->length<8.6) continue;
	if(SimPart2->pdg==2112) continue;
	double Tx2=thetaY-(recon->thetax)[itrk];
	double Ty2=thetaX-(recon->thetay)[itrk];
	//double dx_Temp2=1+pow(TMath::Tan(180*(Tx2)/TMath::Pi()),2)+pow(TMath::Tan(180*(Ty2)/TMath::Pi()),2);
	Scalar=(SimPart2->fpos[2]-SimPart2->ipos[2])*1;
	Norm=TMath::Sqrt((SimPart2->fpos[0]-SimPart2->ipos[0])*(SimPart2->fpos[0]-SimPart2->ipos[0])+(SimPart2->fpos[1]-SimPart2->ipos[1])*(SimPart2->fpos[1]-SimPart2->ipos[1])+(SimPart2->fpos[2]-SimPart2->ipos[2])*(SimPart2->fpos[2]-SimPart2->ipos[2]));
	double Ang3D2=TMath::ACos(Scalar/Norm)*180/TMath::Pi()-(recon->angle)[itrk2];
	
    //double Ang3D2=TMath::ACos(1/dx_Temp2);
	if(TMath::Abs(Ang3D2)<TMath::Abs(Ang3D)) Change=true;
      }
      
      if(Change) continue;
      
      PartNum=is;
      Min=TMath::Abs(Ang3D);
      Particle=TMath::Abs(SimPart2->pdg);
      //TrueParticleNRJ=SimPart2->momentum[3];
      //TrueParticleNRJ=TMath::Sqrt(SimPart2->momentum[0]*SimPart2->momentum[0]+SimPart2->momentum[1]*SimPart2->momentum[1]+SimPart2->momentum[2]*SimPart2->momentum[2]);
      
    }
    else if(Min==180 && TMath::Abs(Ang3D)<Min && SimPart2->length>8.6){

      
      bool Change=false;
      for(int itrk2=0;itrk2<recon->Ntrack;itrk2++){
	if(itrk2==itrk) continue;
	if(SimPart2->length<8.6) continue;
	if(SimPart2->pdg==2112) continue;
	double Tx2=thetaY-(recon->thetax)[itrk];
	double Ty2=thetaX-(recon->thetay)[itrk];
	Scalar=(SimPart2->fpos[2]-SimPart2->ipos[2])*1;
	Norm=TMath::Sqrt((SimPart2->fpos[0]-SimPart2->ipos[0])*(SimPart2->fpos[0]-SimPart2->ipos[0])+(SimPart2->fpos[1]-SimPart2->ipos[1])*(SimPart2->fpos[1]-SimPart2->ipos[1])+(SimPart2->fpos[2]-SimPart2->ipos[2])*(SimPart2->fpos[2]-SimPart2->ipos[2]));
	double Ang3D2=TMath::ACos(Scalar/Norm)*180/TMath::Pi()-(recon->angle)[itrk2];
	
	//double dx_Temp2=1+pow(TMath::Tan(180*(Tx2)/TMath::Pi()),2)+pow(TMath::Tan(180*(Ty2)/TMath::Pi()),2);
	//double Ang3D2=TMath::ACos(1/dx_Temp2);
	if(TMath::Abs(Ang3D2)<TMath::Abs(Ang3D)) Change=true;
      }
      
      if(Change) continue;
      
#ifdef DEBUG
      cout<<"Is selected: Particle="<<SimPart2->pdg<<"       :";
      cout<<"here is the rec thetax, thetay="<<(recon->thetay)[itrk]<<", "<<(recon->thetax)[itrk]<<", Angle 3D="<<(recon->angle)[itrk]<<", And here are the angles="<<thetaX<<", "<<thetaY<<", 3D="<<TMath::ACos(Scalar/Norm)*180/TMath::Pi()<<endl;
      cout<<"rec trk length="<<TrkLength<<", Simpart length="<<SimPart2->length<<endl;
#endif DEBUG
      
      
      PartNum=is;
      Min=TMath::Abs(Ang3D);
      Particle=SimPart2->pdg;
      //TrueParticleNRJ=SimPart2->momentum[3];
      //TrueParticleNRJ=TMath::Sqrt(SimPart2->momentum[0]*SimPart2->momentum[0]+SimPart2->momentum[1]*SimPart2->momentum[1]+SimPart2->momentum[2]*SimPart2->momentum[2]);
      
    }
    else continue;
  }
  
  SimList.push_back(PartNum);
  if(Particle==0){
    
    for(int is=0;is<NSimPart;is++){
      SimPart2=(IngridSimParticleSummary*) evt->GetSimParticle(is);
      if(SimPart2->length<8.6) continue;
      double thetaX=(TMath::ATan((SimPart2->fpos[0]-SimPart2->ipos[0])/(SimPart2->fpos[2]-SimPart2->ipos[2])))*TMath::Pi()/180;
      double thetaY=(TMath::ATan((SimPart2->fpos[1]-SimPart2->ipos[1])/(SimPart2->fpos[2]-SimPart2->ipos[2])))*TMath::Pi()/180;
      
    
      double Tx=thetaY-(recon->thetax)[itrk];
      double Ty=thetaX-(recon->thetay)[itrk];
      //double dx_Temp=1+pow(TMath::Tan(180*(Tx)/TMath::Pi()),2)+pow(TMath::Tan(180*(Ty)/TMath::Pi()),2);
      // double Ang3D=TMath::ACos(1/dx_Temp);
      double Scalar=(SimPart2->fpos[2]-SimPart2->ipos[2])*1;
    double Norm=TMath::Sqrt((SimPart2->fpos[0]-SimPart2->ipos[0])*(SimPart2->fpos[0]-SimPart2->ipos[0])+(SimPart2->fpos[1]-SimPart2->ipos[1])*(SimPart2->fpos[1]-SimPart2->ipos[1])+(SimPart2->fpos[2]-SimPart2->ipos[2])*(SimPart2->fpos[2]-SimPart2->ipos[2]));
    double Ang3D=TMath::ACos(Scalar/Norm)*180/TMath::Pi()-(recon->angle)[itrk];

      
      if(TMath::Abs(Ang3D)<Min && TrkLength<2*SimPart2->length){
	PartNum=is;
	Min=TMath::Abs(Ang3D);
	Particle=TMath::Abs(SimPart2->pdg);
	//TrueParticleNRJ=SimPart2->momentum[3];
	//TrueParticleNRJ=TMath::Sqrt(SimPart2->momentum[0]*SimPart2->momentum[0]+SimPart2->momentum[1]*SimPart2->momentum[1]+SimPart2->momentum[2]*SimPart2->momentum[2]);
	
      }
      else if(Min==180 && TMath::Abs(Ang3D)<Min && SimPart2->length>8.6){
	
	PartNum=is;
	Min=TMath::Abs(Ang3D);
	Particle=SimPart2->pdg;
	//TrueParticleNRJ=SimPart2->momentum[3];
	//TrueParticleNRJ=TMath::Sqrt(SimPart2->momentum[0]*SimPart2->momentum[0]+SimPart2->momentum[1]*SimPart2->momentum[1]+SimPart2->momentum[2]*SimPart2->momentum[2]);
	
      }
      else continue;
    }
  }
  return PartNum;
}

bool Reconstruction::IsINGRID(int mod,int pln,int ch){
  bool Ing;
  if(mod==16){
    if(pln==0) Ing=true;
    else{
      if(ch<=7||ch>=24) Ing=true;
      else Ing=false;
    }
  }
  else Ing=true;
  return Ing;
}

vector <double> Reconstruction::GetMatchingPMINGRID(vector <Hit3D> Vec){
  int nx=0;int ny=0;
  int nx_Ing=0;int ny_Ing=0;
  int ix=0;int iy=0;
  int ix_Ing=0;int iy_Ing=0;
  for(int ihit=0;ihit<Vec.size();ihit++){
    if(Vec[ihit].mod==16){
      if(Vec[ihit].view==0) nx++;
      else ny++;
    }
    else{
      if(Vec[ihit].view==0) nx_Ing++;
      else ny_Ing++;
    }
  }

  double PosZX[nx];
  double PosX[nx];
  double ErrZX[nx];
  double ErrX[nx];
  
  double PosY[ny];
  double PosZY[ny];
  double ErrY[ny];
  double ErrZY[ny];
  
  double PosZX_Ing[nx_Ing];
  double PosX_Ing[nx_Ing];
  double ErrZX_Ing[nx_Ing];
  double ErrX_Ing[nx_Ing];
  
  double PosY_Ing[ny_Ing];
  double PosZY_Ing[ny_Ing];
  double ErrY_Ing[ny_Ing];
  double ErrZY_Ing[ny_Ing];
  //cout<<"new call"<<endl;
  for(int ihit=0;ihit<Vec.size();ihit++){
    //cout<<"Module="<<Vec[ihit].mod<<", x="<<Vec[ihit].x<<", y="<<Vec[ihit].y<<", z="<<Vec[ihit].z<<endl;
    if(Vec[ihit].mod==16){
      if(Vec[ihit].view==0){
	PosZX[ix]=Vec[ihit].z;
	PosX[ix]=Vec[ihit].x;
	if(Vec[ihit].ch<=7||Vec[ihit].ch>=24){
	  ErrX[ix]=2.5;
	  ErrZX[ix]=1.3;
	}
	else{
	  ErrX[ix]=5.;
	  ErrZX[ix]=1.;
	}
	ix++;
      }

      else{
	PosZY[iy]=Vec[ihit].z;
	PosY[iy]=Vec[ihit].y;
	if(Vec[ihit].ch<=7||Vec[ihit].ch>=24){
	  ErrY[iy]=2.5;
	  ErrZY[iy]=1.3;
	}
	else{
	  ErrY[iy]=5.;
	  ErrZY[iy]=1.;
	}
	iy++;
      }
    }
    else{
     if(Vec[ihit].view==0){
	PosZX_Ing[ix_Ing]=Vec[ihit].z;
	PosX_Ing[ix_Ing]=Vec[ihit].x;
	if(Vec[ihit].ch<=7||Vec[ihit].ch>=24){
	  ErrX_Ing[ix_Ing]=2.5;
	  ErrZX_Ing[ix_Ing]=1.3;
	}
	else{
	  ErrX_Ing[ix_Ing]=5.;
	  ErrZX_Ing[ix_Ing]=1.;
	}
	ix_Ing++;
      }

      else{
	PosZY_Ing[iy_Ing]=Vec[ihit].z;
	PosY_Ing[iy_Ing]=Vec[ihit].y;
	if(Vec[ihit].ch<=7||Vec[ihit].ch>=24){
	  ErrY_Ing[iy_Ing]=2.5;
	  ErrZY_Ing[iy_Ing]=1.3;
	}
	else{
	  ErrY_Ing[iy_Ing]=5.;
	  ErrZY_Ing[iy_Ing]=1.;
	}
	iy_Ing++;
      }
    }
  }

  vector <double> MatchingPMINGRID;
  MatchingPMINGRID.clear();
   
  TGraphErrors * HistPosX=new TGraphErrors(nx,PosZX,PosX,ErrZX,ErrX);
  TGraphErrors * HistPosY=new TGraphErrors(ny,PosZY,PosY,ErrZY,ErrY);
  TGraphErrors * HistPosX_Ing=new TGraphErrors(nx_Ing,PosZX_Ing,PosX_Ing,ErrZX_Ing,ErrX_Ing);
  TGraphErrors * HistPosY_Ing=new TGraphErrors(ny_Ing,PosZY_Ing,PosY_Ing,ErrZY_Ing,ErrY_Ing);

  TF1 * fx = new TF1("fx","[0]+x*[1]",0,260);
  TF1 * fy = new TF1("fy","[0]+x*[1]",0,260);
  HistPosX->Fit("fx","Q");
  HistPosY->Fit("fy","Q");


  double SlopeX=fx->GetParameter(1);
  double AngleX=TMath::ATan(SlopeX)*180/TMath::Pi();
  double SlopeY=fy->GetParameter(1);
  double AngleY=TMath::ATan(SlopeY)*180/TMath::Pi();
  //double bX=fx->GetParameter(0);
  //double bY=fy->GetParameter(0);
  double HalfX=fx->Eval(94);
  double HalfY=fy->Eval(94);
  double ReducedChiSquareX=fx->GetChisquare()/fx->GetNDF();
  double ReducedChiSquareY=fy->GetChisquare()/fy->GetNDF();

  TF1 * fx_Ing = new TF1("fx_Ing","[0]+x*[1]",0,260);
  TF1 * fy_Ing = new TF1("fy_Ing","[0]+x*[1]",0,260);

  double SlopeX_Ing=SlopeX;//=fx_Ing->GetParameter(1)=0;
  double AngleX_Ing=AngleX;//=TMath::ATan(SlopeX_Ing)*180/TMath::Pi()=0;
  double SlopeY_Ing=SlopeY;//=fy_Ing->GetParameter(1)=0;
  double AngleY_Ing=AngleY;//=TMath::ATan(SlopeY_Ing)*180/TMath::Pi()=0;
  double HalfX_Ing=HalfX;//=fx_Ing->Eval(94)=0;
  double HalfY_Ing=HalfY;//=fy_Ing->Eval(94)=0;
  double ReducedChiSquareX_Ing=ReducedChiSquareX;//=fx_Ing->GetChisquare()/fx_Ing->GetNDF()=0;
  double ReducedChiSquareY_Ing=ReducedChiSquareY;//=fy_Ing->GetChisquare()/fy_Ing->GetNDF()=0;
  if(nx_Ing!=0 && ny_Ing!=0){
    HistPosX_Ing->Fit("fx_Ing","Q");
    HistPosY_Ing->Fit("fy_Ing","Q");
    SlopeX_Ing=fx_Ing->GetParameter(1);
    AngleX_Ing=TMath::ATan(SlopeX_Ing)*180/TMath::Pi();
    SlopeY_Ing=fy_Ing->GetParameter(1);
    AngleY_Ing=TMath::ATan(SlopeY_Ing)*180/TMath::Pi();
    // bX_Ing=fx_Ing->GetParameter(0);
    // bY_Ing=fy_Ing->GetParameter(0);
    HalfX_Ing=fx_Ing->Eval(94);
    HalfY_Ing=fy_Ing->Eval(94);
    ReducedChiSquareX_Ing=fx_Ing->GetChisquare()/fx_Ing->GetNDF();
    ReducedChiSquareY_Ing=fy_Ing->GetChisquare()/fy_Ing->GetNDF();
  }
#ifdef DEBUG
  cout<<"Angle X, PM="<<AngleX<<", INGRID="<<AngleX_Ing<<endl;
  cout<<"Angle Y, PM="<<AngleY<<", INGRID="<<AngleY_Ing<<endl;
  cout<<"Half X, PM="<<HalfX<<", INGRID="<<HalfX_Ing<<endl;
  cout<<"Half Y, PM="<<HalfY<<", INGRID="<<HalfY_Ing<<endl;
#endif
  MatchingPMINGRID.push_back(AngleX);
  MatchingPMINGRID.push_back(AngleY);
  MatchingPMINGRID.push_back(HalfX);
  MatchingPMINGRID.push_back(HalfY);

  MatchingPMINGRID.push_back(AngleX_Ing);
  MatchingPMINGRID.push_back(AngleY_Ing);
  MatchingPMINGRID.push_back(HalfX_Ing);
  MatchingPMINGRID.push_back(HalfY_Ing);

  return MatchingPMINGRID;
}

double Reconstruction::GetINGRIDTrackWidth(vector <Hit3D> Vec){
  //cout<<"Interaction="<<Num_Int<<endl;                                                                                        
  double HitMax_All[2][12]={{-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999},{-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999}};
  double HitMin_All[2][12]={{999,999,999,999,999,999,999,999,999,999,999,999},{999,999,999,999,999,999,999,999,999,999,999,999}};
  bool HitPln_All[2][12]={{false,false,false,false,false,false,false,false,false,false,false,false},{false,false,false,false,false,false,false,false,false,false,false,false}};
  //cout<<"Initialized"<<endl;                                                                                                  
  int StartXPln=11;
  int StartYPln=11;

  for(int ihit=0;ihit<Vec.size();ihit++){
    if(Vec[ihit].mod!=16){
      HitPln_All[Vec[ihit].view][Vec[ihit].pln]=true;
      if(Vec[ihit].ch>HitMax_All[Vec[ihit].view][Vec[ihit].pln]) HitMax_All[Vec[ihit].view][Vec[ihit].pln]=Vec[ihit].ch;
      if(Vec[ihit].ch<HitMin_All[Vec[ihit].view][Vec[ihit].pln]) HitMin_All[Vec[ihit].view][Vec[ihit].pln]=Vec[ihit].ch;
      if(Vec[ihit].pln<StartXPln && Vec[ihit].view==0) StartXPln=Vec[ihit].pln;
      if(Vec[ihit].pln<StartYPln && Vec[ihit].view==1) StartYPln=Vec[ihit].pln;
      //cout<<"Plane="<<Vec[ihit].pln<<", View="<<Vec[ihit].view<<", ch="<<Vec[ihit].ch<<endl;
      //cout<<"Max="<<HitMax_All[Vec[ihit].view][Vec[ihit].pln]<<", Min="<<HitMin_All[Vec[ihit].view][Vec[ihit].pln]<<endl;
    }
  }
  //cout<<"looped"<<endl;                                                                                                       
  double MeanHit=0;
  int ActPln=0;
  for(int ipln=min(StartXPln, StartYPln);ipln<NPln;ipln++){
    for(int iview=0;iview<2;iview++){
      if((HitPln_All[iview][ipln]==false)) continue;
      else{
	//cout<<"Mean Hit="<<MeanHit<<", ActPln="<<ActPln<<endl;                                                              
	MeanHit+=(HitMax_All[iview][ipln]-HitMin_All[iview][ipln])+1;
	ActPln++;
      }
    }
  }
  
  //cout<<"end of loop"<<endl;                                                                                                  
  MeanHit/=ActPln;
  return MeanHit;
}

vector <double> Reconstruction::GetLastINGRIDChannel(vector <Hit3D> Vec, double TrackSample){
  vector <double> LastChan; LastChan.clear();
  if(TrackSample>=2){
    int MaxPlnX=-1;
    int MaxPlnY=-1;
    int _LastChannelINGRIDX=11.5;
    int _LastChannelINGRIDY=11.5;
    for(int ihit=0;ihit<Vec.size();ihit++){
      if(Vec[ihit].mod!=16 && Vec[ihit].view==0 && Vec[ihit].pln>=MaxPlnX){
	if(Vec[ihit].pln>MaxPlnX) _LastChannelINGRIDX=Vec[ihit].ch;
	else{
	  if(TMath::Abs(Vec[ihit].ch-11.5)>TMath::Abs(_LastChannelINGRIDX-11.5)) _LastChannelINGRIDX=Vec[ihit].ch;
	}
	MaxPlnX=Vec[ihit].pln;
	
      }
      else if(Vec[ihit].mod!=16 && Vec[ihit].view==1 && Vec[ihit].pln>=MaxPlnY){
	if(Vec[ihit].pln>MaxPlnY) _LastChannelINGRIDY=Vec[ihit].ch;
	else{
	  if(TMath::Abs(Vec[ihit].ch-11.5)>TMath::Abs(_LastChannelINGRIDY-11.5)) _LastChannelINGRIDY=Vec[ihit].ch;
	}
	MaxPlnY=Vec[ihit].pln;
      }
    }
    LastChan.push_back(_LastChannelINGRIDX);
    LastChan.push_back(_LastChannelINGRIDY);
  }
  else{
    LastChan.push_back(-1);
    LastChan.push_back(-1);
  }    
  return LastChan;
}

void Reconstruction::GetSelectionPM(bool * SelectionFV, bool * SelectionOV, PMAnaSummary * recon, bool MC){
  *SelectionFV=false;*SelectionOV=false;
  if(MC){
    *SelectionFV=!(recon->vetowtracking) && !(recon->edgewtracking);
    *SelectionOV=(recon->vetowtracking);
  }
  else{
    *SelectionFV=!(recon->vetowtracking) && !(recon->edgewtracking) && recon->ontime;
    *SelectionOV=(recon->vetowtracking) && recon->ontime;
  }
}

					 
#endif

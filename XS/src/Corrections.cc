#ifndef Corrections_cc
#define Corrections_cc

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include "Corrections.h"
//mettre un setup avec NPln, NView
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
#include "Hit.h"
#include "INGRID_Dimension.cc"
#include "setup.h"
#include "Reconstruction.cc"

Corrections::Corrections(bool PM){
  _isPM=PM;
  IngDimCor = new INGRID_Dimension();
}

Corrections::~Corrections(){
  delete IngDimCor;
}

void Corrections::SetDetector(bool PM){
  _isPM=PM;
}
bool Corrections::GetDetector(){
  return _isPM;
}


double Corrections::GetMCCorrections(double Nsel, int mod){//lui donner a manger C[mod]
  // this correction is non-zero only for INGRID modules so WM and PM are treated the same
  double Ncorr=Nsel/(1+C[mod]/100);
  return Ncorr;
}

vector <Hit3D> Corrections::GetFiberAttenuation(vector <Hit3D> Vec){
  const double LAtt=(_isPM? 241.7: 497.); // WM: based on measured value by ML
  double Lmax=120., offset=(_isPM?0.:10.);

  for(int ihit=0;ihit<Vec.size();ihit++){
    double Lfiber;
    if(Vec[ihit].view==0) Lfiber=(Vec[ihit].x-offset); // ML 2017/06/13 swap x<->y
    else Lfiber=(Lmax-Vec[ihit].y-offset); // WM:transverse coord are in [10,110]

    if(!_isPM) Lfiber+=20.; // in wmmc the distance scinti-bundle was added for the attenuation

    Vec[ihit].pecorr=Vec[ihit].pe/TMath::Exp(-Lfiber/LAtt);
#ifdef DEBUG2
    cout<<"**************************************************"<<endl;
    cout<<"Test of Correction::GetFiberAttenuation"<<endl;
    cout<<"Hit view="<<Vec[ihit].view<<", estimated position=("<<Vec[ihit].x<<","<<Vec[ihit].y<<","<<Vec[ihit].z<<")"<<endl;
    cout<<"Charge before ="<<Vec[ihit].pe<<" and after correction="<<Vec[ihit].pecorr<<endl;
#endif
    //if(Vec[ihit].mod==16 && (Vec[ihit].ch>7 && Vec[ihit].ch<24)) cout<<"pe="<<Vec[ihit].pe<<", corr="<<Vec[ihit].pecorr<<endl;
  }
  return Vec;
}

double Corrections::dzWM(double angle, double theta, bool grid){
  double L=(grid?25.:3.);
  double tanThetaLim=(grid ? 3./25. : 25./3.);
  return L/cos(angle)/(1+fabs(tan(theta))/tanThetaLim);
}

void Corrections::GetDXCorrectionWM(vector <Hit3D> & Vec, double angle3D, double thetax, double thetay){
  for(int i=0;i<Vec.size();i++){
    if(Vec[i].mod==15){//WM
      int view=Vec[i].view;
      bool grid=Vec[i].ch>=40;
      Vec[i].pecorr=Vec[i].pecorr/Corrections::dzWM(angle3D,(view==0?thetax:thetay),grid)*3; // output is pe/3mm
    }
  }
}

vector <Hit3D> Corrections::GetDXCorrection(vector <Hit3D> Vec,double dx){

  Reconstruction * Reco = new Reconstruction(_isPM);
  
  for(int i=0;i<Vec.size();i++){
    if(Vec[i].mod==16){//PM
      if(Reco->Reconstruction::IsINGRID(Vec[i].mod,Vec[i].pln,Vec[i].ch)){
	Vec[i].pecorr=Vec[i].pecorr/dx;
      }
      else {
	Vec[i].pecorr=(Vec[i].pecorr/(1.3*dx))/INGRIDSCIBAR;
      }
    }
    else if(Vec[i].mod<14){//INGRID
      Vec[i].pecorr=Vec[i].pecorr/dx;
    }
  }

  delete Reco;  
  return Vec;
}


void Corrections::GetHitPlaneCorrection(vector <Hit3D> Vec){

  TF1 * fx = new TF1("fx","[0]+x*[1]",0,130);
  TF1 * fy = new TF1("fy","[0]+x*[1]",0,130);
  
  int HitNum[NPln][NView];
  for(int ipln=0;ipln<NPln;ipln++){
    for(int iview=0;iview<NView;iview++){
      HitNum[ipln][iview]=0;
    }
  }
  for(int ihit=0;ihit<Vec.size();ihit++){
    HitNum[Vec[ihit].pln][Vec[ihit].view]++;
  }

  int Mod=Vec[0].mod;
  for(int iview=0;iview<NView;iview++){
    for(int ipln=0;ipln<NPln;ipln++){
      double *posxy=new double();
      double *posz=new double();
      int chfixed=0;
      bool done=IngDimCor->INGRID_Dimension::get_posXY(Mod,iview,ipln,chfixed,posxy,posz);
      double Zini=*posz-.5;
      double Zfin=*posz+.5;
      double Xini, Xfin;
      if(iview==0){
        Xini=fx->Eval(Zini);
        Xfin=fx->Eval(Zfin);
      }
      else{
        Xini=fy->Eval(Zini);
        Xfin=fy->Eval(Zfin);
      }
      int Sini=(int) Xini/5;
      int Sfin=(int) Xfin/5;
    }
  }  

  for(int ipln=0;ipln<NPln;ipln++){
    for(int iview=0;iview<NView;iview++){
      for(int ihit=0;ihit<Vec.size();ihit++){
        if((Vec[ihit].pln!=ipln)||(Vec[ihit].view==iview)) continue;
        Vec[ihit].pecorr /= HitNum[ipln][iview];
      }
    }  
  } 
}

#endif

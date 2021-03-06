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

INGRID_Dimension * IngDimCor = new INGRID_Dimension();

double Corrections::GetMCCorrections(double Nsel, int mod){//lui donner a manger C[mod]
  double Ncorr=Nsel/(1+C[mod]/100);
  return Ncorr;
}

vector <Hit3D> Corrections::GetFiberAttenuation(vector <Hit3D> Vec){
  const double LAtt=244;
  for(int ihit=0;ihit<Vec.size();ihit++){
    double L;
    if(Vec[ihit].view==0) L=(120-Vec[ihit].x);
    else L=(120-Vec[ihit].y);//changer pour le PM...
    Vec[ihit].pecorr=Vec[ihit].pe/TMath::Exp(-L/LAtt);
  }
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

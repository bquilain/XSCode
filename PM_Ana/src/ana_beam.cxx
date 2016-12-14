#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>

#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "setup.hxx"

#include "TApplication.h"


#include "ana_beam.hxx"
#include "setup.hxx"

ana_beam::ana_beam(){
};

Bool_t ana_beam::set_pedestal_and_gain(int file_number){
  char FileName[400];
  sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid.gap.txt",gain_and_pedestal_folder,file_number);
  ifstream f_gain(FileName);
  FileStat_t fs;
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    return 0;
  }

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){

	f_gain>>nmod>>ntfb>>nch>>pedestal[nummod][numtfb][numch]>>pedestal_sigma[nummod][numtfb][numch]>>gain[nummod][numtfb][numch];

      }//numch
    }//numtfb
  }//nummod
  return 1;
};//end of set_pedestal_and_gain

void ana_beam::adc_to_pe(Int_t *Adc,Int_t numcyc){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	
	if(*(Adc+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)<-100){//0 suppress data
	  pe[nummod][numtfb][numch][numcyc]=-10;
	}
	else{
	  pe[nummod][numtfb][numch][numcyc]=1.0*(*(Adc+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)-pedestal[nummod][numtfb][numch])/gain[nummod][numtfb][numch];
	}
      }//numch
    }//numtfb
  }//nummod
};//end of adc_to_pe

void ana_beam::set_Tdc(Long_t *tdc,Int_t numcyc){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	Tdc[nummod][numtfb][numch][numcyc]=(*(tdc+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc));
      }//numch
    }//numtfb
  }//nummod

};//end of set tdc

void ana_beam::set_pe(Double_t *npe,Int_t numcyc){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	pe[nummod][numtfb][numch][numcyc]=(*(npe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc));
      }//numch
    }//numtfb
  }//nummod

}//end of set tdc

Bool_t ana_beam::ana_plane_activity(int numpl,int numcyc,double PE_CUT,Long_t TDC_CUT){
  Bool_t flag_x=false,flag_y=false;
  Tdc_x[0][numpl][numcyc]=60000000;
  Pe_x[0][numpl][numcyc]=-10;
  Tdc_y[0][numpl][numcyc]=60000000;
  Pe_y[0][numpl][numcyc]=-10;
  for(Int_t xch=0;xch<LayerNumCh;xch++){
    if(pe[0][numpl][xch][numcyc]>PE_CUT&&Tdc[0][numpl][xch][numcyc]<10000000){
      flag_x=true;
      hit_chx=xch;
      Pe_x[0][numpl][numcyc]=pe[0][numpl][xch][numcyc];
      Tdc_x[0][numpl][numcyc]=Tdc[0][numpl][xch][numcyc];
      for(Int_t ych=LayerNumCh;ych<UseNumCh;ych++){
	if(pe[0][numpl][ych][numcyc]>PE_CUT&&Tdc[0][numpl][ych][numcyc]<10000000){
	  flag_y=true;
	  hit_chy=ych;
	  Tdc_y[0][numpl][numcyc]=Tdc[0][numpl][ych][numcyc];
	  Pe_y[0][numpl][numcyc]=pe[0][numpl][ych][numcyc];
	  if(fabs(Tdc[0][numpl][hit_chx][numcyc]-Tdc[0][numpl][hit_chy][numcyc])<TDC_CUT){
	    return true;
	  }
	}//ify
      }//ych
    }//ifx
  }//xch
  return false;
}

Bool_t ana_beam::top_veto_activity(int numcyc,double PE_CUT){
  for(Int_t numch=0;numch<22;numch++){
    if(pe[0][11][numch][numcyc]>PE_CUT){
      return true;
    }
  }
  return false;
}
Bool_t ana_beam::bottom_veto_activity(int numcyc,double PE_CUT){
  for(Int_t numch=0;numch<22;numch++){
    if(pe[0][12][numch][numcyc]>PE_CUT){
      return true;
    }
  }
  return false;
}
Bool_t ana_beam::left_veto_activity(int numcyc,double PE_CUT){
  for(Int_t numch=23;numch<46;numch++){
    if(pe[0][12][numch][numcyc]>PE_CUT){
      return true;
    }
  }
  return false;
}







#ifndef _ANA_COSMIC_C
#define _ANA_COSMIC_C

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
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "TApplication.h"
#include "TString.h"
#include "TSystem.h"
#include "setup.hxx"

#include "ana_cosmic.hxx"
ana_cosmic::ana_cosmic(){
}

void ana_cosmic::set_pedestal_and_gain(Int_t file_number){
  char FileName[300];
  sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid.gap.txt",gain_and_pedestal_folder,file_number);
  FileStat_t fs;
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    cin>>NCyc;
    return;
  }
  ifstream f(FileName);
  Double_t temp;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	f>>temp>>temp>>temp>>pedestal[nummod][numtfb][numch]>>temp>>gain[nummod][numtfb][numch];

	if(gain[nummod][numtfb][numch]<1){
	  gain[nummod][numtfb][numch]=100;
	  pedestal[nummod][numtfb][numch]=200;
	}
      }//numch
    }//numtfb
  }//nummod
}

void ana_cosmic::ana_trigger_tpl_ver1(ULong64_t twl){
  Bool_t bit[16];
  Int_t a,b;
  //cout<<"trigger "<<twl<<endl;
  a=twl;
  for(Int_t i=0;i<16;i++){
    b=a%2; if(b==1)bit[i]=true;else bit[i]=false;
    a=a/2;
  }
  for(Int_t i=0;i<11;i++){
    if(bit[15-i]){
      tpl_trigger_flag[i]=true;
    }
    if(!bit[15-i]){
      tpl_trigger_flag[i]=false;
    }
  }
}

void ana_cosmic::ana_trigger_type(ULong64_t twl){
  int bit[16];
  Long_t a,b;
  a=twl;
  for(Int_t i=0;i<16;i++){
    bit[i]=(a%16); 
    a=a/16;
  }
  if(bit[15-2]==8)trigger_type=128;//Cosmic
  if(bit[15-3]==1)trigger_type=1;//Beam
  if(bit[15-3]==2)trigger_type=2;//Spill
}
void ana_cosmic::ana_event_number(ULong64_t twl){
  event_number = twl & 0x00000000FFFFFFFF;
}


Bool_t ana_cosmic::pattern_match(Int_t alg){
  Int_t a=alg;
  for(Int_t i=0;i<11;i++){
    if(a%2==1){pattern[10-i]=true;}
    else{pattern[10-i]=false;}
    a=a/2;
  }
  for(Int_t i=0;i<11;i++){
    if(pattern[i]){
      if(!tpl_trigger_flag[i]){
	return false;
      }
    }
  }
  return true;
}


void ana_cosmic::print_pattern(){
  cout<<"pattern"<<endl;
  for(Int_t i=0;i<11;i++){
    if(pattern[i])cout<<1;
  if(!pattern[i])cout<<0;
  }
  cout<<endl;
  cout<<"trigger"<<endl;
  for(Int_t i=0;i<11;i++){
    if(tpl_trigger_flag[i])cout<<1;
  if(!tpl_trigger_flag[i])cout<<0;
  }
  cout<<endl;
}

void ana_cosmic::convert_adc_to_pe(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  pe[nummod][numtfb][numch][numcyc]=1.0*(Adc[nummod][numtfb][numch][numcyc]-pedestal[nummod][numtfb][numch])/gain[nummod][numtfb][numch];
	}//numcyc
      }//numch
    }//numtfb
  }//nummod
}//convert_adc_to_pe


void ana_cosmic::cut_with_pe(Double_t cut){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  Hit[nummod][numtfb][numch][numcyc]=false;
	  if(pe[nummod][numtfb][numch][numcyc]>cut)Hit[nummod][numtfb][numch][numcyc]=true;
	}//numcyc
      }//numchc
    }//numtfb
  }//nummod
}//convert_adc_to_pe

void ana_cosmic::init_Sync_Hit(int numcyc){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	
	  Sync_Hit[nummod][numtfb][numch][numcyc]=false;
	
      }
    }
  }
}

void ana_cosmic::initialize_tript_data(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  Adc[nummod][numtfb][numch][numcyc]=0;
	  Tdc[nummod][numtfb][numch][numcyc]=60000000;
	  pe[nummod][numtfb][numch][numcyc]=0;
	  Hit[nummod][numtfb][numch][numcyc]=false;
	  Tdc_Hit[nummod][numtfb][numch][numcyc]=false;
	  Sync_Hit[nummod][numtfb][numch][numcyc]=false;
	  Act_Hit[nummod][numtfb][numch][numcyc]=false;
	}//numcyc
      }//numchc
    }//numtfb
  }//nummod
}//convert_adc_to_pe


static const Long_t TDC_MIN=50000;
static const Long_t TDC_MAX=80000;
static const Int_t TDC_STEP=5;
static const Int_t N_STEP=(TDC_MAX-TDC_MIN)/TDC_STEP;
void ana_cosmic::cut_with_tdc(){

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  if(Tdc[nummod][numtfb][numch][numcyc]<10000000){
	    Tdc_Hit[nummod][numtfb][numch][numcyc]=true;
	  }
	}//numcyc
      }//numchc
    }//numtfb
  }//nummod
 
}

void ana_cosmic::tdc_synchronize(Long_t cut){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  if(Hit[nummod][numtfb][numch][numcyc]&&Tdc_Hit[nummod][numtfb][numch][numcyc]){
	    for(Int_t numtfb2=0;numtfb2<UseNumTFB;numtfb2++){
	      for(Int_t numch2=0;numch2<UseNumCh;numch2++){
		//if(numtfb2!=numtfb&&numch2!=numch){
		if(!(numtfb2==numtfb&&numch2==numch)){
		  if(Hit[nummod][numtfb2][numch2][numcyc]&&Tdc_Hit[nummod][numtfb2][numch2][numcyc]){
		    if(fabs(Tdc[nummod][numtfb][numch][numcyc]-Tdc[nummod][numtfb2][numch2][numcyc])<cut){
		      Sync_Hit[nummod][numtfb][numch][numcyc]=true;
		    }//abs
		  }//Hit2
		}//numtfb2!=numtfb
	      }//numch2
	    }//numtfb2
	  }//Hit
	}//numcyc
      }//numch
    }//numtfb
  }//nummod
}

void ana_cosmic::tdc_synchronize_old(Long_t cut){
  //make 24(arbit number) group. each group has tdc and number of menbers
  //a hit channel has Tdc value t1
  //if there is tdc near t1, number of members increase
  //if there is not tdc near t1, make new group
  Long_t T[24],N[24];
  Int_t N_group=0;
  Int_t G_Number[NumMod][NumTFB][NumCh][NumCyc];
  //G_Number=-2 initial
  //G_Number=-1 no hit
  for(int i=0;i<24;i++){
    T[i]=0;N[i]=0;
  }
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  G_Number[nummod][numtfb][numch][numcyc]=-2;//initial
	  if(Sync_Hit[nummod][numtfb][numch][numcyc]){
	    if(N_group==0){//first group
	      T[0]=Tdc[nummod][numtfb][numch][numcyc];
	      N[0]++;
	      G_Number[nummod][numtfb][numch][numcyc]=N_group;
	      N_group++;
	    }//end of first group
	    if(N_group>0){//more then second group
	      for(Int_t n=0;n<N_group;n++){
		if(abs(Tdc[nummod][numtfb][numch][numcyc]-T[n])<cut){
		  G_Number[nummod][numtfb][numch][numcyc]=n;
		  N[n]++;
		}//if
	      }
	      if(G_Number[nummod][numtfb][numch][numcyc]==-2){//make new group
	      T[N_group]=Tdc[nummod][numtfb][numch][numcyc];
	      N[N_group]++;
	      G_Number[nummod][numtfb][numch][numcyc]=N_group;
	      N_group++;
	      }//end of new group
	    }//and of second group
	  }//if Sync_Hit
	  else{
	    G_Number[numtfb][numtfb][numch][numcyc]=-1;
	  }
	}//numcyc
      }//numchc
    }//numtfb
  }//nummod

  Int_t selected_group=-1;
  for(Int_t n=0;n<N_group;n++){
    if(N[n]>2){

      break;
    }
  }

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	}
      }
    }
  }
}//convert_adc_to_pe


void ana_cosmic::only_one_activity(){
  Int_t count_X,count_Y;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	count_X=0;
	for(Int_t numch=0;numch<24;numch++){

	  if(Sync_Hit[nummod][numtfb][numch][numcyc])count_X++;

	}//numch
	if(count_X==1){
	  for(Int_t numch=0;numch<24;numch++){
	    for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	      Act_Hit[nummod][numtfb][numch][numcyc]=Sync_Hit[nummod][numtfb][numch][numcyc];
	   
	    }//numch
	  }
	}
	else{
	  for(Int_t numch=0;numch<24;numch++){
	    Act_Hit[nummod][numtfb][numch][numcyc]=false;
	    
	  }//numch
	}


	count_Y=0;
	for(Int_t numch=24;numch<UseNumCh;numch++){
	  
	  if(Sync_Hit[nummod][numtfb][numch][numcyc])count_Y++;
	  
	}//numch
	if(count_Y==1){
	  for(Int_t numch=24;numch<UseNumCh;numch++){
	    Act_Hit[nummod][numtfb][numch][numcyc]=Sync_Hit[nummod][numtfb][numch][numcyc];
	    
	  }//numch
	}
	else{
	  for(Int_t numch=24;numch<UseNumCh;numch++){
	    
	    Act_Hit[nummod][numtfb][numch][numcyc]=false;
	    
	  }//numch
	}
	  
      }//
    }//numtfb
  }//nummod
}



Int_t ana_cosmic::number_of_activity_plane(Int_t numcyc){
  Int_t number=0;
  Bool_t flag_x,flag_y;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      flag_x=false;
      flag_y=false;
      for(Int_t numch=0;numch<24;numch++){
	//if(Sync_Hit[nummod][numtfb][numch][numcyc]){
	if(Act_Hit[nummod][numtfb][numch][numcyc]){
	  flag_x=true;
	}//
      }//numch_X
      for(Int_t numch=24;numch<UseNumCh;numch++){
	//if(Sync_Hit[nummod][numtfb][numch][numcyc]){
	if(Act_Hit[nummod][numtfb][numch][numcyc]){
	  flag_y=true;
	}//
      }//numch_Y
      if(flag_x&&flag_y)number++;
    }//numtfb
  }//nummod
  return number;
}

void ana_cosmic::fit_track(Int_t numcyc){

  const double Scinti_width=5.0;
  const double Scinti_thickness=1.0;
  const double Iron_thickness=6.4;
  const double Plane_thickness=3.8;
  const double Error=pow(3.0,2);
  TGraphErrors *g_x,*g_y;

  int hitn=16;
  Double_t x[hitn],x_error[hitn],y[hitn],y_error[hitn],z[hitn],z_error[hitn];
  Int_t hit_number_x=0,hit_number_y=0;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t numch=0;numch<24;numch++){
	if(Act_Hit[nummod][numtfb][numch][numcyc]){
	  x[hit_number_x]=numch*Scinti_width;
	  x_error[hit_number_x]=Scinti_width/Error;
	  z[hit_number_x]=(Plane_thickness+Iron_thickness)*numtfb;
	  x_error[hit_number_x]=Scinti_thickness/Error;
	  hit_number_x++;
	}//
      }//numch_X
    }//numtfb
  }//nummod
  //fit
  double a=0;
  if(hit_number_x>0){
    //if(hit_number_x>0){
    TF1 *pol1 = new TF1("pol1","pol1",-10000,10000);
    g_x=new TGraphErrors(hit_number_x,z,x,z_error,x_error);
    g_x->Fit("pol1","q","",z[0],z[hit_number_x-1]);
    //cos_x[numcyc]=TMath::ATan(pol1->GetParameter(1));
    cos_x[numcyc]=pol1->GetParameter(1);
    a=pol1->GetParameter(0);
    delete pol1;
  }
  else{
    cos_x[numcyc]=-100;
  }

  //calculate expect hit channel
  for(Int_t numtfb=0;numtfb<11;numtfb++){
    for(Int_t numch=0;numch<UseNumCh;numch++){
      Expect_Hit_Channel[0][numtfb][numcyc][0]=-10;
      Expect_Hit_Channel[0][numtfb][numcyc][1]=-10;
    }
  }
  for(Int_t numtfb=0;numtfb<11;numtfb++){
    double pos = 1.0*a + 1.0*cos_x[numcyc]*(Plane_thickness+Iron_thickness)*numtfb;
    //cout<<"ch:"<<numtfb+1<<"\t"<<pos<<endl;
    for(Int_t numch=0;numch<24;numch++){
      if(Scinti_width*(numch+0.5)>pos&&pos>=Scinti_width*(numch-0.5)){
	Expect_Hit_Channel[0][numtfb][numcyc][0]=numch;
      }
    }
  }
  //

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t numch=24;numch<UseNumCh;numch++){
	if(Act_Hit[nummod][numtfb][numch][numcyc]){
	  y[hit_number_y]=numch*Scinti_width;
	  y_error[hit_number_y]=Error;
	  z[hit_number_y]=(Plane_thickness+Iron_thickness)*numtfb;
	  z_error[hit_number_y]=Scinti_thickness/Error;
	  hit_number_y++;
	}//
      }//numch_Y
    }//numtfb
  }//nummod


  if(hit_number_y>0){
    TF1 *pol2 = new TF1("pol1","pol1",-10000,10000);
    g_y=new TGraphErrors(hit_number_y,z,y,z_error,y_error);
    g_y->Fit("pol1","q","",z[0],z[hit_number_y-1]);
    //cos_y[numcyc]=TMath::ATan(pol1->GetParameter(1));
    cos_y[numcyc] = pol2->GetParameter(1);
    a=pol2->GetParameter(0);
    delete pol2;
  }
  else{
    cos_y[numcyc]=-100;
  }
  //calculate expect hit channel
  for(Int_t numtfb=0;numtfb<11;numtfb++){
    double pos = a + cos_y[numcyc]*(Plane_thickness+Iron_thickness)*numtfb;
    for(Int_t numch=24;numch<UseNumCh;numch++){
      if(Scinti_width*(numch+0.5)>pos&&pos>=Scinti_width*(numch-0.5)){
	Expect_Hit_Channel[0][numtfb][numcyc][1]=numch;
      }
    }
  }
  //

  //caluculate path length
  Double_t x_0 = 0;
  Double_t y_0 = 0;
  Double_t x_1 = 1.0*cos_x[numcyc];
  Double_t y_1 = 1.0*cos_y[numcyc];
  path_length[numcyc]=sqrt(x_1*x_1+y_1*y_1+1.0);
  if(cos_y[numcyc]==-100||cos_x[numcyc]==-100)path_length[numcyc]=1;
  //delete g_x;
  //delete g_y;
}



void ana_cosmic::debug(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  if(numcyc==18){
	    cout<<numtfb<<"\t"<<numch<<"\t"<<numcyc<<"\t"<<Adc[nummod][numtfb][numch][numcyc]<<"\t";
	    if(Hit[nummod][numtfb][numch][numcyc])cout<<1<<endl;
	    if(!Hit[nummod][numtfb][numch][numcyc])cout<<0<<endl;
	  }
	}//numcyc
      }//numchc
    }//numtfb
  }//nummod
}//debug





#endif

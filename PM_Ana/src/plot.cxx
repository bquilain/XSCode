///////////////////////////////////////////////////////////
////this is the class to monitor
//// gain
//// pedestal
//// sigma of pedestal
////at offline
////made by Masashi Otani
////last update 2009/04/17
/////////////////////////////////////////////////////////// 

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
#include <TPostScript.h>

#include "TApplication.h"
#include "plot.hxx"
#include "root_setup.hxx"
#define MAXHIST 200
#define MINHIST 100
const int USEUseNumCh=48;




plot::plot(int ffile_number_1,int ffile_number_2){
  file_number_1=ffile_number_1;
  file_number_2=ffile_number_2;
  flag=true;
  root_setup froot_setup;
  c1 = new TCanvas("c1","c1",10,10,700,700);
  struct stat st;
  //read gain and pedestal

  sprintf(buff1,"%s/ingrid_%08d_0000.daq.mid.gap.txt",gain_and_pedestal_folder,file_number_1);
  if(stat(buff1,&st)==-1){
    cout<<"error:"<<buff1<<endl;
    flag=false;
  }
  ifstream gain_file_1(buff1);
  sprintf(buff1,"%s/ingrid_%08d_0000.daq.mid.gap.txt",gain_and_pedestal_folder,file_number_2);
  if(stat(buff1,&st)==-1){
    cout<<"error:"<<buff1<<endl;
    flag=false;
  }
  ifstream gain_file_2(buff1);
  sprintf(buff1,"%s/ingrid_%08d_0000.daq.mid.noise.txt",noise_folder,file_number_1);
  if(stat(buff1,&st)==-1){
    cout<<"error:"<<buff1<<endl;
    flag=false;
  }
  ifstream noise_file_1(buff1);
  sprintf(buff1,"%s/ingrid_%08d_0000.daq.mid.noise.txt",noise_folder,file_number_2);
  if(stat(buff1,&st)==-1){
    cout<<"error:"<<buff1<<endl;
    flag=false;
  }
  ifstream noise_file_2(buff1);

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB+2;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	gain_file_1>>Mod[nummod]>>TFB[numtfb]>>Ch[numch]>>pedestal_1[nummod][numtfb][numch]>>pedestal_sigma_1[nummod][numtfb][numch]>>gain_1[nummod][numtfb][numch];
	gain_file_2>>Mod[nummod]>>TFB[numtfb]>>Ch[numch]>>pedestal_2[nummod][numtfb][numch]>>pedestal_sigma_2[nummod][numtfb][numch]>>gain_2[nummod][numtfb][numch];

	gain_ratio[nummod][numtfb][numch]=(gain_2[nummod][numtfb][numch]-gain_1[nummod][numtfb][numch])/gain_1[nummod][numtfb][numch];
	//error message
	flag_gain_error[nummod][numtfb][numch]=false;
	if(fabs(gain_ratio[nummod][numtfb][numch])>ERROR_RATIO){
	  flag_gain_error[nummod][numtfb][numch]=true;
	  cout<<"TPL "<<numtfb+1<<" numch"<<numch+1<<" gain ratio "<<gain_ratio[nummod][numtfb][numch]<<endl;
	  if(fabs(gain_ratio[nummod][numtfb][numch])>MAX_RATIO){
	    if(gain_ratio[nummod][numtfb][numch]>0)gain_ratio[nummod][numtfb][numch]=MAX_RATIO-0.01;
	    if(gain_ratio[nummod][numtfb][numch]<0)gain_ratio[nummod][numtfb][numch]=-MAX_RATIO+0.01;
	  }
	}

	pedestal_sigma_ratio[nummod][numtfb][numch]=(pedestal_sigma_2[nummod][numtfb][numch]-pedestal_sigma_1[nummod][numtfb][numch])/pedestal_sigma_1[nummod][numtfb][numch];
	noise_file_1>>Mod[nummod]>>TFB[numtfb]>>Ch[numch]>>noise_1[nummod][numtfb][numch];
	noise_file_2>>Mod[nummod]>>TFB[numtfb]>>Ch[numch]>>noise_2[nummod][numtfb][numch];

	noise_ratio[nummod][numtfb][numch]=(noise_2[nummod][numtfb][numch]-noise_1[nummod][numtfb][numch])/noise_2[nummod][numtfb][numch];

	if(fabs(noise_ratio[nummod][numtfb][numch])>MAX_RATIO)noise_ratio[nummod][numtfb][numch]=MAX_RATIO-0.01;
      }
      for(Int_t numch=0;numch<NumCh-UseNumCh;numch++){
	gain_file_1>>temp>>temp>>temp>>temp>>temp>>temp;
	gain_file_2>>temp>>temp>>temp>>temp>>temp>>temp;
	noise_file_1>>temp>>temp>>temp>>temp;
	noise_file_2>>temp>>temp>>temp>>temp;
      }
      this->veto_correction(nummod,numtfb);
    }//numtfb
  }//nummod
}



void plot::plot_gain(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB+2;numtfb++){
      graph_gain_1[nummod][numtfb]=new TGraph(UseNumCh,Ch,gain_1[nummod][numtfb]);
      graph_gain_2[nummod][numtfb]=new TGraph(UseNumCh,Ch,gain_2[nummod][numtfb]);
      c1->Clear();
      TH1F *frame = gPad->DrawFrame(0,0,50,13);
      frame->SetXTitle("Ch");
      frame->SetYTitle("gain[ADC counts]");
      graph_gain_1[nummod][numtfb]->SetMarkerStyle(8);
      graph_gain_1[nummod][numtfb]->SetMarkerSize(1.2);
      graph_gain_1[nummod][numtfb]->SetMarkerColor(1);
      
      sprintf(buff1,"run%d:black",file_number_1);
      TText *ref=new TText(1,0.5,buff1);
      ref->SetTextColor(1);
      ref->SetBit(kCanDelete);
      ref->SetTextSize(0.03);      

      graph_gain_2[nummod][numtfb]->SetMarkerStyle(8);
      graph_gain_2[nummod][numtfb]->SetMarkerSize(1.2);
      graph_gain_2[nummod][numtfb]->SetMarkerColor(2);
      sprintf(buff1,"run%d:red",file_number_2);
      TText *now=new TText(1,1.2,buff1);
      now->SetTextColor(2);
      now->SetBit(kCanDelete);
      now->SetTextSize(0.03);

      sprintf(buff1,"TFB%02d",numtfb+1);
      frame->SetTitle(buff1);
      frame->Draw();
      ref->Draw();
      now->Draw();
      graph_gain_1[nummod][numtfb]->Draw("P");
      graph_gain_2[nummod][numtfb]->Draw("P");
      
      sprintf(buff1,"TPL%02d gain",numtfb+1);
      c1->SetName(buff1);
      c1->Update();
      c1->Write();

      delete ref,now;
    }
  }
}



void plot::plot_pedestal(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB+2;numtfb++){
      graph_pedestal_1[nummod][numtfb]=new TGraph(UseNumCh,Ch,pedestal_1[nummod][numtfb]);
      graph_pedestal_2[nummod][numtfb]=new TGraph(UseNumCh,Ch,pedestal_2[nummod][numtfb]);
      c1->Clear();
      TH1F *frame = gPad->DrawFrame(0,0,50,200);
      frame->SetXTitle("Ch");
      frame->SetYTitle("pedestal[ADC counts]");
      graph_pedestal_1[nummod][numtfb]->SetMarkerStyle(8);
      graph_pedestal_1[nummod][numtfb]->SetMarkerSize(1.2);
      graph_pedestal_1[nummod][numtfb]->SetMarkerColor(1);
      
      sprintf(buff1,"run%d:black",file_number_1);
      TText *ref=new TText(1,1,buff1);
      ref->SetTextColor(1);
      ref->SetBit(kCanDelete);
      ref->SetTextSize(0.03);      

      graph_pedestal_2[nummod][numtfb]->SetMarkerStyle(8);
      graph_pedestal_2[nummod][numtfb]->SetMarkerSize(1.2);
      graph_pedestal_2[nummod][numtfb]->SetMarkerColor(2);
      sprintf(buff1,"run%d:red",file_number_2);
      TText *now=new TText(1,10,buff1);
      now->SetTextColor(2);
      now->SetBit(kCanDelete);
      now->SetTextSize(0.03);

      sprintf(buff1,"TFB%02d",numtfb+1);
      frame->SetTitle(buff1);
      frame->Draw();
      ref->Draw();
      now->Draw();
      graph_pedestal_1[nummod][numtfb]->Draw("P");
      graph_pedestal_2[nummod][numtfb]->Draw("P");
      
      sprintf(buff1,"TPL%02d pedestal",numtfb+1);
      c1->SetName(buff1);
      c1->Update();
      c1->Write();
      delete ref,now;
    }
  }
}


void plot::plot_pedestal_sigma(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB+2;numtfb++){
      graph_pedestal_sigma_1[nummod][numtfb]=new TGraph(UseNumCh,Ch,pedestal_sigma_1[nummod][numtfb]);
      graph_pedestal_sigma_2[nummod][numtfb]=new TGraph(UseNumCh,Ch,pedestal_sigma_2[nummod][numtfb]);
      c1->Clear();
      TH1F *frame = gPad->DrawFrame(0,0,50,3);
      frame->SetXTitle("Ch");
      frame->SetYTitle("pedestal sigma");
      graph_pedestal_sigma_1[nummod][numtfb]->SetMarkerStyle(8);
      graph_pedestal_sigma_1[nummod][numtfb]->SetMarkerSize(1.2);
      graph_pedestal_sigma_1[nummod][numtfb]->SetMarkerColor(1);
      

      graph_pedestal_sigma_2[nummod][numtfb]->SetMarkerStyle(8);
      graph_pedestal_sigma_2[nummod][numtfb]->SetMarkerSize(1.2);
      graph_pedestal_sigma_2[nummod][numtfb]->SetMarkerColor(2);


      sprintf(buff1,"TFB%02d",numtfb+1);
      frame->SetTitle(buff1);
      frame->Draw();

      sprintf(buff1,"run%d:black",file_number_1);
      TText *ref=new TText(1,0.2,buff1);
      ref->SetTextColor(1);
      ref->SetBit(kCanDelete);
      ref->SetTextSize(0.03);      

      sprintf(buff1,"run%d:red",file_number_2);
      TText *now=new TText(1,1.2,buff1);
      now->SetTextColor(2);
      now->SetBit(kCanDelete);
      now->SetTextSize(0.03);

      ref->Draw();
      now->Draw();
      graph_pedestal_sigma_1[nummod][numtfb]->Draw("P");
      graph_pedestal_sigma_2[nummod][numtfb]->Draw("P");
      
      sprintf(buff1,"TPL%02d pedestal",numtfb+1);
      c1->SetName(buff1);
      c1->Update();
      c1->Write();
      delete ref,now;
    }
  }

}



void plot::hist_gain(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      sprintf(buff1,"Mod.%02d TFB.%02d gain 1",nummod+1,numtfb+1);	
      h_gain_1[nummod][numtfb]=new TH1F(buff1,buff1,100,0,60);
      sprintf(buff1,"Mod.%02d TFB.%02d gain 2",nummod+1,numtfb+1);	
      h_gain_2[nummod][numtfb]=new TH1F(buff1,buff1,100,0,60);
    }//numtfb
  }//nummod


  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){

	if(numch<USEUseNumCh)h_gain_1[nummod][numtfb]->Fill(gain_1[nummod][numtfb][numch]);

	if(numch<USEUseNumCh)h_gain_2[nummod][numtfb]->Fill(gain_2[nummod][numtfb][numch]);
      }
    }
  }
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      c1->Clear();
      h_gain_1[nummod][numtfb]->SetXTitle("gain[ADC counts]");
      h_gain_1[nummod][numtfb]->SetLineColor(1);
      cout<<file_number_1<<":black";
      h_gain_1[nummod][numtfb]->Draw();
      c1->Update();
      c1->Clear();

      h_gain_2[nummod][numtfb]->SetXTitle("gain[ADC counts]");
      h_gain_2[nummod][numtfb]->SetLineColor(2);
      cout<<file_number_2<<":red";
      h_gain_2[nummod][numtfb]->Draw();
      c1->Update();
      //sprintf(buff_1,"%s/gain_hist_%08d_%02d_%02d.pdf",pdf_folder,gain_file_number,nummod,numtfb);
      //c1->Print(buff_1);

    }
  }
}

void plot::hist_gain_all(){

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
 
	if(numtfb<11){
	  if(numch<USEUseNumCh){
	    h_gain_all->Fill(gain_1[nummod][numtfb][numch]);
	  }
	}

	if(numtfb==11){
	  if(numch<22){
	    h_gain_all->Fill(gain_1[nummod][numtfb][numch]);
	  }
	}
	if(numtfb==12){
	  if(numch<22||(numch>23&&numch<46)){
	    h_gain_all->Fill(gain_1[nummod][numtfb][numch]);
	  }
	}

      }//numch
    }//numtfb
  }//numtfb

  c1->Clear();
  h_gain_all->SetXTitle("gain[ADC counts]");
  h_gain_all->Draw();
  c1->Update();

}

void plot::hist_pedestal_sigma(){
  double temp;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      sprintf(buff1,"Mod.%02d TFB.%02d pedestal sigma",nummod+1,numtfb+1);	
      h_pedestal_sigma[nummod][numtfb]=new TH1F(buff1,buff1,10,0,10);
    }//nummod
  }//numtfb


  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){

	if(numch<USEUseNumCh)h_pedestal_sigma[nummod][numtfb]->Fill(pedestal_sigma_1[nummod][numtfb][numch]);
      }
    }
  }
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      c1->Clear();
      h_pedestal_sigma[nummod][numtfb]->SetXTitle("pedestal sigma");
      h_pedestal_sigma[nummod][numtfb]->Draw();
      c1->Update();
 
    }
  }
}

void plot::plot_adc(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	sprintf(buff1,"Mod.%02d TFB.%02d Ch.%02d",nummod+1,numtfb+1,numch+1);	
	hist_Adc[nummod][numtfb][numch]=new TH1F(buff1,buff1,MAXHIST-MINHIST,MINHIST,MAXHIST);
      }//numch
    }
  }  
  cout<<"ok"<<endl;
  sprintf(buff1,"%s/ingrid_%08d_0000.daq.mid.new.root",root_file_folder,file_number_2);

  TFile *f = new TFile(buff1);
  TTree *IngEvt = (TTree*)f->Get("IngEvt");
  Int_t Nevent = IngEvt->GetEntries();
  Int_t Adc[NumMod][NumTFB][UseNumCh][NumCyc];
  sprintf(buff1,"Adc[%d][%d][%d][%d]",NumMod,NumTFB,UseNumCh,NumCyc);
  IngEvt->SetBranchAddress(buff1,Adc);
  cout<<"fill hist..."<<endl;
  for(Int_t nevent=0;nevent<Nevent;nevent++){
    IngEvt->GetEntry(nevent);
    if(nevent+1%100==0)cout<<nevent+1<<"\tend"<<endl;
    for(Int_t nummod=0;nummod<NumMod;nummod++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
	for(Int_t numch=0;numch<UseNumCh;numch++){
	  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	    hist_Adc[nummod][numtfb][numch]->Fill(Adc[nummod][numtfb][numch][numcyc]);
	  }//numcyc
	}//numch
      }//numbtfb
    }//nummod
  }//nevent
  cout<<"fill hist end"<<endl;
  double temp;
  c1->Divide(2,2);

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	Int_t n_c=(numch+1)%4;
	if(n_c==0)n_c=4;
	c1->cd(n_c);
	sprintf(buff1,"Mod.%02d TFB.%02d Ch.%02d",nummod+1,numtfb+1,numch+1);
	hist_Adc[nummod][numtfb][numch]->SetName(buff1);
	hist_Adc[nummod][numtfb][numch]->SetXTitle("adc");
	hist_Adc[nummod][numtfb][numch]->SetYTitle("number");
	hist_Adc[nummod][numtfb][numch]->Draw();
	if((numch+1)%4==0){
	  c1->Update();
	  c1->Update();
	  c1->Clear();
	  c1->Divide(2,2);
	}
      }
    }
  }
}




void plot::plot_noise(){

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      graph_noise_1[nummod][numtfb]=new TGraph(UseNumCh,Ch,noise_1[nummod][numtfb]);
      graph_noise_2[nummod][numtfb]=new TGraph(UseNumCh,Ch,noise_2[nummod][numtfb]);
      c1->Clear();
      TH1F *frame = gPad->DrawFrame(0,0,50,1.0);
      frame->SetXTitle("Ch");
      frame->SetYTitle("noise[MHz]");
      graph_noise_1[nummod][numtfb]->SetMarkerStyle(8);
      graph_noise_1[nummod][numtfb]->SetMarkerSize(1.2);
      graph_noise_1[nummod][numtfb]->SetMarkerColor(1);
      graph_noise_2[nummod][numtfb]->SetMarkerStyle(8);
      graph_noise_2[nummod][numtfb]->SetMarkerSize(1.2);
      graph_noise_2[nummod][numtfb]->SetMarkerColor(2);
      sprintf(buff1,"TFB%02d",numtfb+1);
      frame->SetTitle(buff1);
      frame->Draw();

      sprintf(buff1,"run%d:black",file_number_1);
      TText *ref=new TText(1,0.2,buff1);
      ref->SetTextColor(1);
      ref->SetBit(kCanDelete);
      ref->SetTextSize(0.03);      

      sprintf(buff1,"run%d:red",file_number_2);
      TText *now=new TText(1,1.2,buff1);
      now->SetTextColor(2);
      now->SetBit(kCanDelete);
      now->SetTextSize(0.03);

      ref->Draw();
      now->Draw();

      graph_noise_1[nummod][numtfb]->Draw("P");
      graph_noise_2[nummod][numtfb]->Draw("P");
      sprintf(buff1,"TPL%02d noise",numtfb+1);
      c1->SetName(buff1);
      c1->Update();
      c1->Write();
      delete ref,now;

    }
  }
}

void plot::plot_gain_ratio(){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      graph_gain_ratio[nummod][numtfb]=new TGraph(UseNumCh,Ch,gain_ratio[nummod][numtfb]);
      c1->Clear();
      TH1F *frame = gPad->DrawFrame(0,-MAX_RATIO,50,MAX_RATIO);
      frame->SetXTitle("Ch");
      frame->SetYTitle("gain ratio");
      graph_gain_ratio[nummod][numtfb]->SetMarkerStyle(8);
      graph_gain_ratio[nummod][numtfb]->SetMarkerSize(1.2);
      graph_gain_ratio[nummod][numtfb]->SetMarkerColor(1);
      sprintf(buff1,"TFB%02d",numtfb+1);
      frame->SetTitle(buff1);
      frame->Draw();
      graph_gain_ratio[nummod][numtfb]->Draw("P");
      sprintf(buff1,"TPL%02d gain ratio",numtfb+1);
      c1->SetName(buff1);
      c1->Update();
      c1->Write();
  

     
    }
  }
}


void plot::plot_noise_ratio(){

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      graph_noise_ratio[nummod][numtfb]=new TGraph(UseNumCh,Ch,noise_ratio[nummod][numtfb]);
 
      c1->Clear();
      TH1F *frame = gPad->DrawFrame(0,-MAX_RATIO,50,MAX_RATIO);
      frame->SetXTitle("Ch");
      frame->SetYTitle("noise ratio");
      graph_noise_ratio[nummod][numtfb]->SetMarkerStyle(8);
      graph_noise_ratio[nummod][numtfb]->SetMarkerSize(1.2);
      graph_noise_ratio[nummod][numtfb]->SetMarkerColor(1);
      sprintf(buff1,"TFB%02d",numtfb+1);
      frame->SetTitle(buff1);
      frame->Draw();
      graph_noise_ratio[nummod][numtfb]->Draw("P");
      sprintf(buff1,"TPL%02d noise ratio",numtfb+1);
      c1->SetName(buff1);
      c1->Update();
      c1->Write();
     
    }
  }
}


void plot::evdisp(Double_t *npe){
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	  pe[nummod][numtfb][numch][numcyc]=(*(npe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc));
	}
      }//numch
    }//numtfb
  }//nummod         
  Int_t n;
  Int_t nummod =0;
  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
    cout<<"cycle"<<numcyc+1<<endl;
    cout<<"X"<<endl;
    for(Int_t numtfb=0;numtfb<11;numtfb++){
      for(Int_t numch=0;numch<24;numch++){
	n=pe[nummod][numtfb][numch][numcyc]=(*(npe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc));
	if(n>9)n=9;
	if(n<2)n=0;
	cout<<n;
      }//numch
      cout<<endl;
    }//numtfb
    cout<<"Y"<<endl;
    for(Int_t numtfb=0;numtfb<11;numtfb++){
      for(Int_t numch=24;numch<48;numch++){
	n=pe[nummod][numtfb][numch][numcyc];
	if(n>9)n=9;
	if(n<2)n=0;
	cout<<n;
      }//numch
      cout<<endl;
    }//numtfb
    cin>>n;
  }//nummod         
  

};


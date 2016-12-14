
//////////////////////////////////////////////////
//analysis code for cosmic event
//Made by Masashi Otani
//////////////////////////////////////////////////
/*
|||||||||||||||||||||||||||||||
consists 4 function
1. main
2. ProcessFile
3. Read
4. Analysis
main{
ProcessFile
->Read->Analysis(each event)
}
|||||||||||||||||||||||||||||||
*/
// ND280 software includes
#include "TMidasBank.hxx"
#include "TMidasFile.hxx"
#include "TMidasBankProxy.hxx"
#include "TMidasBankProxyRegistry.hxx"
#include "TND280RawEvent.hxx"
#include "TRawDataHeader.hxx"
// oaRawEvent includes
#include "TTripTDigitBank.hxx"
#include "TMidasTripTDigitItr.hxx"
#include "TMidasTripTDigit.hxx"
#include "TMCMBank.hxx"
#include "TTriggerBank.hxx"
// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
//C++ libraly includes
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>
#include <sys/stat.h>
#include <unistd.h> // using getopt      
using namespace std;
// INGRID includes
#include "setup.hxx"
#include "Beam_EVENT.cxx"
#include "ana_cosmic.cxx"
#include "ana_MPPC.cxx"
#include "ana_each_MidasBank.cxx"
#include "COSMIC_EVENT.cxx"
#include "hit.cxx"

//parameter setting___________________________________________________________________
int flag_txt=0;
int cosmic_flag=0;
ana_cosmic fana_cosmic;
FileStat_t fs;
//ana_MPPC *fana_mppc[NumMod][NumTFB][NumCh];
ana_MPPC *fana_mppc[NumMod*NumTFB*NumCh];
ana_each_MidasBank *fana_each_MidasBank;
Int_t NBEvent;
Int_t temperature;
Int_t HVtrim;
Long_t NumNum;//debug
COSMIC_EVENT *fCosEvt;
hit fHit;
//Int_t number_of_hit[NumCyc];

void Checkout(){
  //Tdc synchronized check
  Int_t hit_channel_x[NumTFB],hit_channel_y[NumTFB];
  Double_t cut=4.5;
  //Long_t diff_Tdc[NumTFB];
  Long_t Tdc_temp[NumTFB];
  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
    Int_t flag_x = fHit.one_hit_every_one_layer_x(&fana_cosmic.pe[0][0][0][0],numcyc,cut,&hit_channel_x[0]);
    Int_t flag_y = fHit.one_hit_every_one_layer_y(&fana_cosmic.pe[0][0][0][0],numcyc,cut,&hit_channel_x[0]);
    if(flag_x==1&&flag_y==1){
      for(Int_t numtfb=0;numtfb<11;numtfb++){
	for(Int_t numch=0;numch<24;numch++){
	  if(fana_cosmic.pe[0][numtfb][numch][numcyc]>cut){
	    Tdc_temp[numtfb]=fana_cosmic.Tdc[0][numtfb][numch][numcyc];
	  }
	}
      }
      for(Int_t numtfb=0;numtfb<11;numtfb++){
	fCosEvt->diff_Tdc[numtfb]=Tdc_temp[0]-Tdc_temp[numtfb];
	//cout<<diff_Tdc[numtfb]<<"\t";
      }
      fCosEvt->BookFill();
    }
  }

  //hit cycle check
  Int_t temp_number_of_hit[NumCyc];
  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
    temp_number_of_hit[numcyc]=0;;
  }
  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<UseNumCh;numch++){
	if(fana_cosmic.pe[0][numtfb][numch][numcyc]>10){
	  fCosEvt->Num_Hit->Fill(numcyc);
	  fCosEvt->number_of_hit[numcyc]++;
	  temp_number_of_hit[numcyc]++;
	}
      }//numch
    }//numtfb
  }//numcyc
  if(temp_number_of_hit[6]>10){
    cout<<"********"<<NumNum<<endl;
  }

}//Checkout

void Debug(){
  //fana_cosmic.debug();
  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      //if(numcyc==18||numcyc==19)cout<<fana_cosmic.Expect_Hit_Channel[0][numtfb][numcyc][0]<<endl;
      for(Int_t numch=0;numch<NumCh;numch++){
	//if(numcyc==18||numcyc==19){//for 1500nsec gate
	if(numcyc==18||numcyc==19){//for 800nsec gate



	  if(fana_cosmic.Sync_Hit[0][numtfb][numch][numcyc]){

	    if(numch<24){
	      if(fabs(fana_cosmic.Expect_Hit_Channel[0][numtfb][numcyc][0]-numch)<2){
		//cout<<numtfb<<"\t"<<numch<<"\t"<<numcyc<<"\t"<<fana_cosmic.pe[0][numtfb][numch][numcyc]<<"\t"<<fana_cosmic.cos_x[numcyc]<<"\t"<<fana_cosmic.cos_y[numcyc]<<endl;
		if(fana_cosmic.pe[0][numtfb][numch][numcyc]>10){
		  fCosEvt->cosmic[0][numtfb][numch]->Fill(fana_cosmic.pe[0][numtfb][numch][numcyc]);
		}
	      }
	    }
	    if(48>numch&&numch>=24){
	      if(fabs(fana_cosmic.Expect_Hit_Channel[0][numtfb][numcyc][1]-numch)<2){
		if(fana_cosmic.pe[0][numtfb][numch][numcyc]>10){

		  fCosEvt->cosmic[0][numtfb][numch]->Fill(fana_cosmic.pe[0][numtfb][numch][numcyc]);

		}
		//cout<<numtfb<<"\t"<<numch<<"\t"<<numcyc<<"\t"<<fana_cosmic.pe[0][numtfb][numch][numcyc]<<"\t"<<fana_cosmic.cos_x[numcyc]<<"\t"<<fana_cosmic.cos_y[numcyc]<<endl;
	      }
	    }
	    if(12>=numtfb&&numtfb>=11&&numch<48){
	      if(fana_cosmic.pe[0][numtfb][numch][numcyc]>10){
		fCosEvt->cosmic[0][numtfb][numch]->Fill(fana_cosmic.pe[0][numtfb][numch][numcyc]);

	      }
	    }
	  }//SyncHit


	  //cout<<numtfb<<"\t"<<numch<<"\t"<<numcyc<<"\t"<<fana_cosmic.pe[0][numtfb][numch][numcyc]<<"\t"<<sqrt(1+fana_cosmic.cos_x[numcyc]*fana_cosmic.cos_x[numcyc]+fana_cosmic.cos_y[numcyc]*fana_cosmic.cos_y[numcyc])<<endl;
	

	  
	}//numcyc
      }//numch
    }//numtfb
  //int temp;
  //cin>>temp;

  }
}

//Analysis___________________________________________________________________
void Analysis(){
  fana_cosmic.convert_adc_to_pe();
  fana_cosmic.cut_with_pe(2.5);
  fana_cosmic.cut_with_tdc();
  fana_cosmic.tdc_synchronize(12);
  fana_cosmic.only_one_activity();
 
  if(fana_cosmic.number_of_activity_plane(18)>3){
    fana_cosmic.fit_track(18);
  }
  else{
    fana_cosmic.init_Sync_Hit(18);
  }


  if(fana_cosmic.number_of_activity_plane(19)>3){
    fana_cosmic.fit_track(19);
  }
  else{
    fana_cosmic.init_Sync_Hit(19);
  }


}
//Aanalysis end___________________________________________________________________

//Read___________________________________________________________________
void Read(ND::TND280RawEvent* re) {
  //Block(re) consists of Header(header) and some Bank(triggerbank,triptbank and so on)
  //read each bank and put data into fana_cosmic
  fana_cosmic.initialize_tript_data();//initilize before read
  fana_each_MidasBank->ana_header(fana_cosmic,re);//read header 
  fana_each_MidasBank->ana_triggerbank(fana_cosmic,re);//read triggerbank
  fana_each_MidasBank->ana_triptBank(fana_cosmic,re);//read adc,tdc data
}
//Read end___________________________________________________________________



//ProcessFile ___________________________________________________________________
void ProcessFile(const char *FileName) {
  //cout<<"processing file..."<<endl;
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    return;
  }
  ND::TMidasFile mf;
  mf.Open(FileName);

  // Loop over events in this file
  // Midas data consists of Blocks(=re).
  // Promote(re->PromoteMidasBanks(false)),
  // read(Read(re))
  // and analysis(fana_cosmic.****) each Block 
  cout<<"loop all event..."<<endl;


  while ( ND::TND280RawEvent* re = mf.ReadRawEvent()) {
   
    if(NumNum==0){
      re->PromoteMidasBanks(false);//first event is header
    }
    re->PromoteMidasBanks(false);
    fana_each_MidasBank->ana_mcmBank(fana_cosmic,re);//read trigger information
    if(fana_cosmic.cosmic_or_not()){//if cosmic trigger read block
      Read(re);
    }
    else{
      //cout<<"event:"<<fana_cosmic.get_event_number()<<"is not cosmic triggerd."<<endl;
    }
    //if(fana_cosmic.event_number%100==0)cout<<"end to read\t"<<fana_cosmic.event_number<<"event"<<endl;
    delete re;

    //analysis
    Analysis();
    Debug();
    Checkout();
    NumNum++;

    
    if(NumNum%1000==0){
      cout<<NumNum<<endl;
      for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
	cout<<fCosEvt->number_of_hit[numcyc]<<"\t";
      }
      cout<<endl;
      
    }


   
  }// End loop over events
  //cout<<"end of loop all event..."<<endl;
}

//_____________________________________________________________________________

int main(int argc,char *argv[]){
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);

  int c=-1;
  char FileName[300],txt_file_name[300];
  int run_number;
  Int_t pedestal_and_gain_for_cosmic;
  Int_t pedestal_number=1976;
  while ((c = getopt(argc, argv, "hr:p:t:o:")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid",data_file_folder,run_number);
      cout<<"file name is :"<<FileName<<endl;

      break;
    case 't':
      temperature=atoi(optarg);
      cout<<"temperature is :"<<temperature<<endl;
      break;
    case 'p':
      pedestal_number=atoi(optarg);
      break;
    case 'o':
      HVtrim=atoi(optarg);
      cout<<"HVtrim is :"<<HVtrim<<endl;
      break;
    case 'h':
      cout<<"-r [run number]\t midas file to root file"<<endl;
 
      break;
    case '?':
      cout<<"Unknown option"<<endl;
      cout<<"-r [run number]\t midas file to root file"<<endl;
      cout<<"-b \t : get all beam event"<<endl;
      exit(1);
      break;
    }
  }//option end

  //write root file
  char root_file_name[300];
  sprintf(root_file_name,"%s/ingrid_%08d_cosmic.root",root_file_folder,run_number);
  fCosEvt = new COSMIC_EVENT(root_file_name);
  for(Int_t numcyc=0;numcyc<NumCyc;numcyc++){
    fCosEvt->number_of_hit[numcyc]=0;
  }
  fCosEvt->HIST();
  fCosEvt->Book();
  //MPPC setup
  fana_cosmic.set_pedestal_and_gain(pedestal_number);
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	fana_mppc[nummod*NumTFB*NumCh+numtfb*NumCh+numch] = new ana_MPPC();
	fana_mppc[nummod*NumTFB*NumCh+numtfb*NumCh+numch]->set_gain_and_pedestal(14,709,nummod,numtfb,numch);//(temperature,HVtrim,mod,tfb,ch)
      }
    }
  }

  NumNum=0;


  ProcessFile(FileName);
  fCosEvt->BookWrite();
  fCosEvt->HISTWrite();
  return 0;

}

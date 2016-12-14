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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>
#include <sys/stat.h>
using namespace std;
// INGRID includes
#include "setup.hxx"
#include "INGRID_EVENT.cxx"
#include "Beam_EVENT.cxx"
#include "ana_cosmic.cxx"
#include "ana_MPPC.cxx"




#include <unistd.h> // using getopt      
ofstream writing_file;
ofstream cosmic_file;
INGRID_EVENT *fIngEvt;
Beam_EVENT *fBeamEvt;
int flag_txt=0;
int cosmic_flag=0;
ana_cosmic fana_cosmic;
Bool_t beam_mode=false;
Bool_t normal_mode=false;
Bool_t beambackground_mode=false;
FileStat_t fs;
//ana_MPPC *fana_mppc[NumMod][NumTFB][NumCh];
ana_MPPC *fana_mppc[NumMod*NumTFB*NumCh];
Int_t NBEvent;
int run_number;


void Event(ND::TND280RawEvent* re) {
  fIngEvt->Initialize();

  //Time Stamp
  ND::TRawDataHeader header =re->GetHeader();
  fIngEvt->UTime=header.GetTimeStamp();


  //trigger type information 
  ND::THandle<ND::TMCMBank> mcmBank;
  while ( mcmBank = re->GetMidasBank<ND::TMCMBank>("IMCM",mcmBank) ) {
    ND::TMCMBank& mcm = re->UseMidasBank<ND::TMCMBank>("IMCM");
    ULong64_t twl=mcm.GetTriggerWord();
    fana_cosmic.ana_trigger_type(twl);
    fana_cosmic.ana_event_number(twl);
  }




  //trigger information for cosmic
  ND::THandle<ND::TTriggerBank> triggerBank;
  while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
    ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
    ULong64_t t =trigger.GetTriggerPatternHigh();
    ULong64_t twl = trigger.GetTriggerPatternLow();
    ULong64_t twm = trigger.GetTriggerPatternMid();
    ULong64_t twh = trigger.GetTriggerPatternHigh();
    fana_cosmic.ana_trigger_tpl_ver1(twl);
    fIngEvt->TType=fana_cosmic.get_trigger_type();
    fIngEvt->SNum=(trigger.GetTriggerWord()>>32)&0xffff;
  }

  // Loop over all banks of type TTripTHitBank
  ND::THandle<ND::TTripTDigitBank> triptBank;
  while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {
    // Create an iterator over digits
    ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
    while ( ! itr.EOD() ) {

      ND::TMidasTripTDigit digit(itr.Get());
      Int_t lowadc  =  digit.GetLowGainADC();
      Int_t highadc =  digit.GetHighGainADC();
      Int_t chan64  =  digit.GetChannelNum() + 16*digit.GetTripTNum();
      Int_t tfb     =  digit.GetTFBNum();
      Int_t time    =  digit.GetTimeOffset();
      Int_t cycle   =  digit.GetIntegrationNum();
      Long_t offset =  digit.GetTimeOffsetT0();



      const ND::TRawDataHeader& h = re->GetHeader();
      int evno = h.GetSeqNo();    // Event Sequence Number
      // Plot the event
      int iint = digit.GetIntegrationNum();    // = Capacitor number
      int icoff = triptBank->GetTFBStartCycle();
      int co = iint - icoff;
      if (co<0) co += 23;

      ////////////////////////////////////////////////////////
      //if(time>10000000)cout<<highadc<<"\t"<<time<<"\t"<<digit.GetTime()<<"\t"<<digit.GetTimeOffsetT0()<<endl;
      ////////////////////////////////////////////////////////

      fIngEvt->Adc[0][tfb][chan64][co]=highadc;
      fIngEvt->LoAdc[0][tfb][chan64][co]=lowadc;
      fIngEvt->Tdc[0][tfb][chan64][co]=time;
      fIngEvt->Time[0][tfb][chan64][co]=digit.GetTime();
      fIngEvt->Offset[0][tfb][chan64][co]=offset;
      //fIngEvt->T0[0][tfb][chan64][co]=digit.GetTimeOffsetT0();
      //if(digit.GetTimeOffsetTI()!=0)cout<<digit.GetTimeOffsetTI()<<endl;

      //writing cosmic event

     //pettern match
      Bool_t flag=false;
      //cout<<"event\t"<<fana_cosmic.get_event_number()<<endl;
      if(cosmic_flag==1&&time<16777201&&(co==cosmic_event_cycle_1||co==cosmic_event_cycle_2)){
	if(tfb==0){if(fana_cosmic.pattern_match(2816)){flag=true;}}
	if(tfb==1){if(fana_cosmic.pattern_match(2432)){flag=true;}}
	if(tfb==2){if(fana_cosmic.pattern_match(2240)||fana_cosmic.pattern_match(3584)){flag=true;}}
	if(tfb==3){if(fana_cosmic.pattern_match(2816)||fana_cosmic.pattern_match(2144)){flag=true;}}
	if(tfb==4){if(fana_cosmic.pattern_match(2432)||fana_cosmic.pattern_match(2096)){flag=true;}}
	if(tfb==5){if(fana_cosmic.pattern_match(2240)||fana_cosmic.pattern_match(3584)){flag=true;}}
	if(tfb==6){if(fana_cosmic.pattern_match(2144)||fana_cosmic.pattern_match(2060)){flag=true;}}
	if(tfb==7){if(fana_cosmic.pattern_match(2096)||fana_cosmic.pattern_match(2054)){flag=true;}}
	if(tfb==8){if(fana_cosmic.pattern_match(2072)||fana_cosmic.pattern_match(2051)){flag=true;}}
	if(tfb==9){if(fana_cosmic.pattern_match(2060)){flag=true;}}
	if(tfb==10){if(fana_cosmic.pattern_match(2054)){flag=true;}}
      }

      if(flag){
	//cout<<tfb<<"\t"<<chan64<<"\t"<<co<<"ok"<<endl;
	cosmic_file<<fana_cosmic.get_event_number()<<"\t";
	cosmic_file<<tfb<<"\t";
	cosmic_file<<chan64<<"\t";
	cosmic_file<<1.0*(highadc-fana_cosmic.pedestal[0][tfb][chan64])/fana_cosmic.gain[0][tfb][chan64]<<"\t";
	//cosmic_file<<time<<"\t";
	cosmic_file<<endl;
      }


      if(flag_txt==1){
	//digit.Print();

	writing_file<<highadc<<"\t";
	writing_file<<lowadc<<"\t";
	writing_file<<time<<"\t";
	writing_file<<digit.GetRMMNum()<<"\t";
	writing_file<<tfb<<"\t";
	writing_file<<cycle<<"\t";
	writing_file<<digit.GetTripTNum()<<"\t";
	writing_file<<chan64<<"\t";
	writing_file<<digit.GetFormat()<<endl;
	
      }
// A little more work is needed to develop the time variable.  The hardware records
// a 44-bit offset from the nearest hardware 1-sec reset. To compress this, the
// online subtracts a fixed offset and a further offset depending on the cycle
// number.  It fits it into 24 bits.  We have plans to later further compress the
// data and store this in a 12 bit number.  For now, the digit.GetTimeOffset();
// only returns the compressed time.  This should normally be OK.  We will add
// another two methods: GetFullTime() which will add back the offsets and give 
// the 44-bit number and GetTimeStatus() which will return codes for e.g.
// 'there was no TDC hit on this channel'.  At the moment 'no TDC hit on this channel'
// is flagged by digit.GetTimeOffset() returning the value RAWDPT_DVX_FMT0_TIMENOHIT.
// (which is 0x00FFFFF1).

// This is why there is a big spike of entries in the right hand bin on the 
// gTime plot we make below).

    }   // End of loop over digits in this bank

  }     // End of loop over banks of digits in this event
  fIngEvt->EvtNum=fana_cosmic.get_event_number();
  fIngEvt->BookFill();
  if(fana_cosmic.get_trigger_type()==1){
    
    cout<<"event:"<<""<<fana_cosmic.get_event_number()<<"\tTime:"<<fIngEvt->UTime<<"\t beam event"<<endl;
    //debug


  }

}


void BeamEvent(ND::TND280RawEvent* re,ana_MPPC *fana_mppc[]) {
  fBeamEvt->Initialize();
  Int_t SNum;
  Long_t UTime;
  ULong64_t twl;

  //trigger type information 
  ND::THandle<ND::TMCMBank> mcmBank;
  while ( mcmBank = re->GetMidasBank<ND::TMCMBank>("IMCM",mcmBank) ) {
    ND::TMCMBank& mcm = re->UseMidasBank<ND::TMCMBank>("IMCM");
    twl=mcm.GetTriggerWord();
    fana_cosmic.ana_trigger_type(twl);
    fana_cosmic.ana_event_number(twl);
  }

  Bool_t real_beam_event=false;//there are dummy beam triggers
  if(1==fana_cosmic.get_trigger_type()){//beam event
    //Time Stamp
    ND::TRawDataHeader header =re->GetHeader();
    UTime=header.GetTimeStamp();
    //trigger information 
    ND::THandle<ND::TTriggerBank> triggerBank;
    while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
      ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
      SNum=(trigger.GetTriggerWord()>>32)&0xffff;
    }
    char buff1[300];
    sprintf(buff1,"%s/ingrid_%08d_beam_shot_unixtime.txt",beam_folder,run_number);
    if ( gSystem->GetPathInfo(buff1,fs) ) {
      std::cerr << "Cannot find file: " << buff1 << std::endl;
      exit(1);
    }
    
    ifstream f_spill_and_unix(buff1);
    Long_t daq_UTime=-10;
    Long_t daq_SNum=-10;
    while(!f_spill_and_unix.eof()){
      f_spill_and_unix>>daq_SNum>>daq_UTime;
      Long_t daq_SNum_ = (daq_SNum+1) & 0xffff;
      if(abs(UTime-daq_UTime)<8&&SNum==daq_SNum_){
	NBEvent++;
	cout<<"***total "<<NBEvent<<"*******"<<endl;
	cout<<"Beam DAQ Spill\t"<<daq_SNum_<<"("<<daq_SNum<<")"<<endl;
	cout<<"INGRID Spill\t"<<SNum<<endl;
	cout<<"DAQ time\t"<<daq_UTime<<endl;
	cout<<"INGRID time\t"<<UTime<<endl;
	cout<<"event:"<<fana_cosmic.get_event_number()<<endl;
	cout<<"****************************"<<endl;
	real_beam_event=true;
      }//if(UTime-daq_UTime)
    }//f_spill_and_unix.eof
    f_spill_and_unix.close();
  }//!f_spill_and_unix.eof()


  if(real_beam_event){
  ND::THandle<ND::TTripTDigitBank> triptBank;
  while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {
    // Create an iterator over digits
    ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
    while ( ! itr.EOD() ) {
      ND::TMidasTripTDigit digit(itr.Get());
      Int_t lowadc  =  digit.GetLowGainADC();
      Int_t highadc =  digit.GetHighGainADC();
      Int_t chan64  =  digit.GetChannelNum() + 16*digit.GetTripTNum();
      Int_t tfb     =  digit.GetTFBNum();
      Int_t time    =  digit.GetTimeOffset();
      Int_t cycle   =  digit.GetIntegrationNum();
      Long_t offset =  digit.GetTimeOffsetT0();
      const ND::TRawDataHeader& h = re->GetHeader();
      int evno = h.GetSeqNo();    // Event Sequence Number
      // Plot the event
      int iint = digit.GetIntegrationNum();    // = Capacitor number
      int icoff = triptBank->GetTFBStartCycle();
      int co = iint - icoff;
      if (co<0) co += 23;
      fBeamEvt->pe[0][tfb][chan64][co]=(highadc-fana_mppc[tfb*NumCh+chan64]->get_pedestal())/fana_mppc[tfb*NumCh+chan64]->get_gain();
      fBeamEvt->Adc[0][tfb][chan64][co]=highadc;
      fBeamEvt->LoAdc[0][tfb][chan64][co]=lowadc;
      fBeamEvt->Tdc[0][tfb][chan64][co]=time;
      fBeamEvt->SNum=SNum;
      fBeamEvt->UTime=UTime;
    }   // End of loop over digits in this bank
  }     // End of loop over banks of digits in this event
  //fIngEvt->BookFill();
  fBeamEvt->BookFill();  
  }//BeamEvent

}

//____________________________________________________________________________


//____________________________________________________________________________

void ProcessFile(const char *FileName) {
  
  cout<<"processing file..."<<endl;
  ND::TMidasBankProxyRegistry::Instance().Print();
  
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    return;
  }

  ND::TMidasFile mf;
  mf.Open(FileName);

  // Loop over events in this file
  int Nevent=0;

  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
      for(Int_t numch=0;numch<NumCh;numch++){
	fana_mppc[nummod*NumTFB*NumCh+numtfb*NumCh+numch] = new ana_MPPC();
	fana_mppc[nummod*NumTFB*NumCh+numtfb*NumCh+numch]->set_gain_and_pedestal(run_number,nummod,numtfb,numch);
//// now tempolary use temp=14,MPPC_HV=709
//// 090429 Masashi Otani
//// you should consider stability of pedestal and 
//// the difference of temperature at each run
	
      }
    }
  }

  while ( ND::TND280RawEvent* re = mf.ReadRawEvent() ) {
    if(Nevent==0)re->PromoteMidasBanks(false);//first event is header
    if(beam_mode){if(Nevent%10000==0)cout<<"start to read\t"<<Nevent<<"event"<<endl;}
    if(normal_mode){if(Nevent%100==0)cout<<"start to read\t"<<Nevent<<"event"<<endl;}

    Nevent++;
    re->PromoteMidasBanks(false);
    if(normal_mode)Event(re);
    if(beam_mode)BeamEvent(re,fana_mppc);
    if(beambackground_mode)BeamEvent(re,fana_mppc);
    delete re;
  }  // End loop over events
}

//_____________________________________________________________________________

int main(int argc,char *argv[]){
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);

  int c=-1;
  char root_file_name[300],FileName[300],txt_file_name[300];

  Int_t pedestal_and_gain_for_cosmic;
  beam_mode=false;
  normal_mode=false;
  beambackground_mode=false;
  while ((c = getopt(argc, argv, "hb:r:c:t:g")) != -1) {
    switch(c){
    case 'r':
      normal_mode=true;
      run_number=atoi(optarg);
      sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid",data_file_folder,run_number);
      cout<<"file name is :"<<FileName<<endl;
      sprintf(root_file_name,"%s/ingrid_%08d_0000.daq.mid.new.root",output_file_folder,run_number);
      cout<<"out put file name is :"<<root_file_name<<endl;
      break;
    case 'c':
      cosmic_flag=1;
      pedestal_and_gain_for_cosmic=atoi(optarg);
      cout<<"comic analysis mode..."<<endl;
      break;
    case 't':
      cout<<"output txt file"<<endl;
      flag_txt=1;
      break;
    case 'b':
      beam_mode=true;//get beam and
      run_number=atoi(optarg);
      sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid",data_file_folder,run_number);
      cout<<"file name is :"<<FileName<<endl;
      sprintf(root_file_name,"%s/ingrid_%08d_0000.beam.root",beam_folder,run_number);
      cout<<"out put file name is :"<<root_file_name<<endl;
      NBEvent=0;
      break;
    case 'g':
      beambackground_mode=true;//get beam background
      run_number=atoi(optarg);
      sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid",data_file_folder,run_number);
      cout<<"file name is :"<<FileName<<endl;
      sprintf(root_file_name,"%s/ingrid_%08d_0000.beam.root",beam_folder,run_number);
      cout<<"out put file name is :"<<root_file_name<<endl;
      NBEvent=0;
      break;
    case 'h':
      cout<<"-r [run number]\t midas file to root file"<<endl;
      cout<<"-b \t : get all beam event"<<endl;
      cout<<"-g \t : back ground "<<endl;
      cout<<"-c \t : get all cosmic analysis mode"<<endl;
      cout<<"-t \t : output txt file"<<endl;
      break;
    case '?':
      cout<<"Unknown option"<<endl;
      cout<<"-r [run number]\t midas file to root file"<<endl;
      cout<<"-b \t : get all beam event"<<endl;
      cout<<"-c \t : get all cosmic analysis mode"<<endl;
      cout<<"-t \t : output txt file"<<endl;
      exit(1);
      break;
    }
  }//option end

  if(flag_txt==1){//output txt file
    sprintf(txt_file_name,"%s/ingrid_%08d_0000.daq.mid.txt",output_file_folder,run_number);
    writing_file.open(txt_file_name);
  }

  if(normal_mode){//midas to root for each run(normal mode)
    cout<<"normal mode..."<<endl;
    fIngEvt = new INGRID_EVENT(root_file_name);
    fIngEvt->Book();
    if(cosmic_flag==1){
      char cosmic_file_name[300];
      sprintf(cosmic_file_name,"%s/ingrid_%08d_cosmic.txt",cosmic_txt_folder,run_number);
      cosmic_file.open(cosmic_file_name);
      fana_cosmic.set_pedestal_and_gain(pedestal_and_gain_for_cosmic);
    }
    cout<<"start to midas to root..."<<endl;
    ProcessFile(FileName);
    fIngEvt->BookWrite();
    if(flag_txt)writing_file.close(); 
    if(cosmic_flag==1)cosmic_file.close();
  }//normal mode
  char buff1[300];
  FileStat_t fs;
  if(beam_mode||beambackground_mode){//beam mode
    if ( gSystem->GetPathInfo(FileName,fs) ) {
      std::cerr << "Cannot find file: " << FileName << std::endl;
      return false;
    }
    fBeamEvt = new Beam_EVENT(root_file_name);
    fBeamEvt->Book();
    cout<<FileName<<" read start"<<endl;
    ProcessFile(FileName);
    fBeamEvt->BookWrite();   
  }



  return 0;

}

#ifndef _ANA_EACH_MIDASBANK_C
#define _ANA_EACH_MIDASBANK_C

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


#include "ana_each_MidasBank.hxx"


ana_each_MidasBank::ana_each_MidasBank(){};
ana_each_MidasBank::~ana_each_MidasBank(){};

void ana_each_MidasBank::ana_triggerbank(ana_cosmic &fana_cosmic,ND::TND280RawEvent *re){
  ND::THandle<ND::TTriggerBank> triggerBank;
  while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
    ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
    ULong64_t t =trigger.GetTriggerPatternHigh();
    ULong64_t twl = trigger.GetTriggerPatternLow();
    ULong64_t twm = trigger.GetTriggerPatternMid();
    ULong64_t twh = trigger.GetTriggerPatternHigh();
    fana_cosmic.ana_trigger_tpl_ver1(twl);
  }
}

void ana_each_MidasBank::ana_header(ana_cosmic &fana_cosmic,ND::TND280RawEvent *re){
  ND::TRawDataHeader header =re->GetHeader();
  fana_cosmic.NSeq = header.GetSeqNo();    // Event Sequence Number
  fana_cosmic.UTime=header.GetTimeStamp();
}

void ana_each_MidasBank::ana_triptBank(ana_cosmic &fana_cosmic,ND::TND280RawEvent *re){
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

      int iint = digit.GetIntegrationNum();    // = Capacitor number
      int icoff = triptBank->GetTFBStartCycle();
      int co = iint - icoff;
      if (co<0) co += 23;
      fana_cosmic.NCyc=co;

      fana_cosmic.Adc[0][tfb][chan64][co]=highadc;
      fana_cosmic.Tdc[0][tfb][chan64][co]=time;
    }
  }
}

void ana_each_MidasBank::ana_mcmBank(ana_cosmic &fana_cosmic,ND::TND280RawEvent *re){
  ND::THandle<ND::TMCMBank> mcmBank;
  while ( mcmBank = re->GetMidasBank<ND::TMCMBank>("IMCM",mcmBank) ) {
    ND::TMCMBank& mcm = re->UseMidasBank<ND::TMCMBank>("IMCM");
    ULong64_t twl=mcm.GetTriggerWord();
    fana_cosmic.ana_trigger_type(twl);
    fana_cosmic.ana_event_number(twl);
  }
}
#endif

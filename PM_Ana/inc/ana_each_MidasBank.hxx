#ifndef _ANA_EACH_MIDASBANK_H
#define _ANA_EACH_MIDASBANK_H

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

#include "ana_cosmic.hxx"
//#include "setup.hxx"
//#include "ana_MPPC.cxx"

class ana_each_MidasBank{
public:
  ana_each_MidasBank();
  ~ana_each_MidasBank();
  void ana_triggerbank(ana_cosmic &fana_cosmic,ND::TND280RawEvent *re);
  void ana_header(ana_cosmic &fana_cosmic,ND::TND280RawEvent *re);
  void ana_triptBank(ana_cosmic &ana_cosmic,ND::TND280RawEvent *re);
  void ana_mcmBank(ana_cosmic &ana_cosmic,ND::TND280RawEvent *re);


};

#endif

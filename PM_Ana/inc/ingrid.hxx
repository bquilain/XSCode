////////////////////////////////////
///////// INGRID data class
///////// ver1.9
////////////////////////////////////


#ifndef __INGRID_HXX__
#define __INGRID_HXX__
#include "TObject.h"
#include "TClonesArray.h"

#define VetoRight 11
#define VetoLeft 12
#define VetoBottom 13
#define VetoUp 14
#define XLayer 0
#define YLayer 1

class IngHit : public TObject {
public:
  Int_t          HighAdc;             //High gain ADC
  Int_t          LowAdc;              //Low  gain ADC
  Int_t          HighPe;              //Photoelectron from High gain ADC
  Int_t          LowPe;               //Photoelectron from Low  gain ADC
  Long_t         Time;                //Time from Trigger[nsec]
  //channel ID
  Int_t          View;                //Side View(X layer):0 Top View(Y laer):1
  Int_t          Pln;                 //0~10:TPL, 11~14:Veto
  Int_t          Ch;                  //0~23:TPL, 0~21:Veto
  IngHit();
  ~IngHit();
  ClassDef(IngHit,1)
}; 

class ModEvt : public TObject {
public:
  //Int_t          Mod;
  Int_t          NHitX;
  TClonesArray*  fXHits;
  Int_t          NHitXVeto;
  TClonesArray*  fXVetoHits;
  Int_t          NHitY;
  TClonesArray*  fYHits;
  Int_t          NHitYVeto;
  TClonesArray*  fYVetoHits;
  void AddHit(Int_t  highadc,
              Int_t  lowadc,
              Long_t tdc, 
              Int_t  view, 
              Int_t  pln, 
              Int_t  ch    );

  ModEvt();
  void InitHit();
  ~ModEvt();
  ClassDef(ModEvt,1)
}; 

class IngEvt : public TObject {
public:

  Int_t          Event;                //Event number
  Int_t          Spill;                //Beam Spill number
  Int_t          Cycle;                //Cycle number(0~22)
  Int_t          TrgId;                //Trigger mode
                                       //1:Beam
                                       //2:Periodic
                                       //128:Cosmic
  Long_t         TrgTime;              //Trigger Time[10nsec]
  Long_t         UTime;                //Header Time[10nsec]
  Int_t          NHit;                 //Number of Hit,
  ModEvt         fModEvt[16];
  IngEvt();
  ~IngEvt();
  void Make();
  void AddHit(Int_t  highadc,
              Int_t  lowadc,
              Long_t tdc, 
              Int_t  mod, 
              Int_t  view, 
              Int_t  pln, 
              Int_t  ch    );
  void InitHits();
  void Print();
  ClassDef(IngEvt,1)
}; 

class IngHdr : public TObject {
public:
  Int_t           Run;
  IngHdr();
  ~IngHdr(); 
  ClassDef(IngHdr,1)
}; 

#endif

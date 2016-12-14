/////////////////////////////////////
////////  INGRID data class
////////  ver1.9
/////////////////////////////////////

#ifndef __INGRID_CXX__
#define __INGRID_CXX__
#include "TROOT.h"
#include "iostream"
#include "sstream"
#include <new>
using namespace std;
#include "ingrid.hxx"
ClassImp(IngHit)
ClassImp(ModEvt)
ClassImp(IngEvt)
ClassImp(IngHdr)

  
IngHit::IngHit(){};
IngHit::~IngHit(){}

ModEvt::ModEvt(){
  fXHits          =  new TClonesArray("IngHit",264);
  fXVetoHits      =  new TClonesArray("IngHit",44);
  fYHits          =  new TClonesArray("IngHit",264);
  fYVetoHits      =  new TClonesArray("IngHit",44);
};


void ModEvt::InitHit(){
  fXHits          -> Delete();
  fXVetoHits      -> Delete();
  fYHits          -> Delete();
  fYVetoHits      -> Delete();
  NHitX           =  0;
  NHitXVeto       =  0;
  NHitY           =  0;
  NHitYVeto       =  0;
};
void ModEvt::AddHit(Int_t highadc, Int_t lowadc, Long_t time, 
                     Int_t view, Int_t pln, Int_t ch){
  if(view==0){
    if(pln==VetoBottom||
       pln==VetoUp){
      TClonesArray& IngHits = *(fXVetoHits);
      IngHit*       hit     = new(IngHits[NHitXVeto++]) IngHit();
      hit  -> HighAdc       =  highadc;
      hit  -> LowAdc        =  lowadc;
      hit  -> Time          =  time;
      hit  -> View          =  view;
      hit  -> Pln           =  pln;
      hit  -> Ch            =  ch;
    }
    else{
      TClonesArray& IngHits = *(fXHits);
      IngHit*       hit     = new(IngHits[NHitX++]) IngHit();
      hit  -> HighAdc       =  highadc;
      hit  -> LowAdc        =  lowadc;
      hit  -> Time          =  time;
      hit  -> View          =  view;
      hit  -> Pln           =  pln;
      hit  -> Ch            =  ch;
    }
  }
  if(view==1){
    if(pln==VetoRight||
       pln==VetoLeft){
      TClonesArray& IngHits = *(fYVetoHits);
      IngHit*       hit     = new(IngHits[NHitYVeto++]) IngHit();
      hit  -> HighAdc       =  highadc;
      hit  -> LowAdc        =  lowadc;
      hit  -> Time          =  time;
      hit  -> View          =  view;
      hit  -> Pln           =  pln;
      hit  -> Ch            =  ch;
    }
    else{
      TClonesArray& IngHits = *(fYHits);
      IngHit*       hit     = new(IngHits[NHitY++]) IngHit();
      hit  -> HighAdc       =  highadc;
      hit  -> LowAdc        =  lowadc;
      hit  -> Time          =  time;
      hit  -> View          =  view;
      hit  -> Pln           =  pln;
      hit  -> Ch            =  ch;
    }
  }
};

ModEvt::~ModEvt(){
};


IngEvt::IngEvt(){
  for(Int_t nummod=0;nummod<16;nummod++){
    //fModEvt[nummod].Mod=nummod;
    //fModEvt[nummod]=ModEvt();
    fModEvt[nummod].InitHit();
  } 
};
void IngEvt::Make(){
};

void IngEvt::AddHit(Int_t highadc, Int_t lowadc, Long_t time, 
                    Int_t mod, Int_t view, Int_t pln, Int_t ch){
  fModEvt[mod].AddHit(highadc, lowadc, time, view, pln, ch);
};
void IngEvt::InitHits(){
  for(Int_t nummod=0;nummod<16;nummod++)fModEvt[nummod].InitHit(); 
};
void IngEvt::Print(){
};

IngEvt::~IngEvt(){};

IngHdr::IngHdr(){};
IngHdr::~IngHdr(){};

#endif

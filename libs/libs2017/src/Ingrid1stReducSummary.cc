#include "Ingrid1stReducSummary.h"


//......................................................................

Ingrid1stReducSummary::Ingrid1stReducSummary():
  nhits(0)
{ 
  hitmod   =    -1;
  hitcyc   =    -1;
  xhitbit  = 0x000;
  xhitbit  = 0x000;
  nhitxlyr =    -1;
  nhitylyr =    -1;
  xtotpe   = -1e-5;
  ytotpe   = -1e-5;
  xtracklike = false;
  ytracklike = false;

  hitx.clear();
  hity.clear();
  hitxz.clear();
  hityz.clear();
}
//......................................................................
Ingrid1stReducSummary::Ingrid1stReducSummary(const Ingrid1stReducSummary& reduc) :
  nhits(0)
{ 
  hitmod    =  reduc.hitmod;
  hitcyc    =  reduc.hitcyc;
  xhitbit   =  reduc.xhitbit;
  yhitbit   =  reduc.yhitbit;
  nhitxlyr  =  reduc.nhitxlyr;
  nhitylyr  =  reduc.nhitylyr;
  xtotpe    =  reduc.xtotpe;
  ytotpe    =  reduc.ytotpe;
  xtracklike=  reduc.xtracklike;
  ytracklike=  reduc.ytracklike;
  nhits     = 0;
  for (int i=0; i < INGRIDHIT_MAX_1STREDUC; ++i) {
    fIngridHit[i] = TRef(NULL);
  }
  for (int i=0; i < reduc.nhits; ++i) 
    AddIngridHit(reduc.GetIngridHit(i));

  hitx  = reduc.hitx;
  hity  = reduc.hity;
  hitxz = reduc.hitxz;
  hityz = reduc.hityz;
  //for (int i=0; i < reduc.hitx.size(); ++i)
  //hitx.push_back( reduc.hitx[i] );



}

//......................................................................



Ingrid1stReducSummary::~Ingrid1stReducSummary() 
{

}

//......................................................................


void Ingrid1stReducSummary::Clear(Option_t* option)
{
  for (int i=0; i<INGRIDHIT_MAX_1STREDUC; ++i)
    fIngridHit[i] = TRef(NULL);
  nhits = 0;

}

//......................................................................

void Ingrid1stReducSummary::Print()
{

}
//......................................................................
void Ingrid1stReducSummary::AddIngridHit(IngridHitSummary* hit) 
{

    if (nhits < INGRIDHIT_MAX_1STREDUC) {
        fIngridHit[nhits] = TRef((IngridHitSummary*) hit);
        ++nhits;
    }

}
//......................................................................


IngridHitSummary* Ingrid1stReducSummary::GetIngridHit(int i) const
{ 

  return (IngridHitSummary*)fIngridHit[i].GetObject();
}


//......................................................................


ClassImp(Ingrid1stReducSummary)

////////////////////////////////////////////////////////////////////////

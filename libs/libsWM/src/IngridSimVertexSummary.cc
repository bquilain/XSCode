#include "IngridSimVertexSummary.h"

//......................................................................

IngridSimVertexSummary::IngridSimVertexSummary() 
{ 
  nutype= -1;
  nuE   =  0;
  norm = -1.e5;
  for(int i=0;i<3;i++) numomentum[i] = -1.e5;

  mod   = -1;
  targeta = -1;
  targetz = -1;
  targettype = -1;
  vnuclini = -1.e5;
  pfsurf = -1.e5;
  inttype = -1;
  totcrsne = -1.e5;
  xnu = -1.e5;
  ynu = -1.e5;
  znu = -1.e5;
}

//......................................................................

IngridSimVertexSummary::~IngridSimVertexSummary() 
{ 
}

//......................................................................

void IngridSimVertexSummary::Clear(Option_t* option)
{    
  nutype= -1;
  nuE   =  0;
  norm = -1.e5;
  for(int i=0;i<3;i++) numomentum[i] = -1.e5;

  mod   = -1;
  targeta = -1;
  targetz = -1;
  targettype = -1;
  vnuclini = -1.e5;
  pfsurf = -1.e5;
  inttype = -1;
  totcrsne = -1.e5;
  xnu = -1.e5;
  ynu = -1.e5;
  znu = -1.e5;

  ng = -1;
  gpid.clear();
  gmec.clear();
  gcosbm.clear();
  gposx.clear();
  gposy.clear();
  gposz.clear();
  gmomx.clear();
  gmomy.clear();
  gmomz.clear();
}

//......................................................................

void IngridSimVertexSummary::Print(Option_t* option) const
{
}

//......................................................................

ClassImp(IngridSimVertexSummary)

////////////////////////////////////////////////////////////////////////



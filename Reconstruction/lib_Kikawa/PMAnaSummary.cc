#include "PMAnaSummary.h"

//......................................................................

PMAnaSummary::PMAnaSummary():
  nhits(0), ntracks(0)
{ 
  Ntrack      = -1; 
  Ningtrack   = -1; 
  clstime     = -1.e5;  
  clstimecorr = -1.e5;
  exptime     = -1.e5;
  hitcyc      = -1;
  ontime      =  false;
  vetowtracking=  false;
  edgewtracking=  false;

  x.clear();
  y.clear();
  z.clear();
  zx.clear();
  zy.clear();
  startxpln.clear();
  startypln.clear();
  startxch.clear();
  startych.clear();
  endxpln.clear();
  endypln.clear();
  endxch.clear();
  endych.clear();
  thetax.clear();
  thetay.clear();
  angle.clear();
  ing_startmod.clear();
  ing_endmod.clear();
  ing_startpln.clear();
  ing_endpln.clear();
  ing_trk.clear();
  pm_stop.clear();
  ing_stop.clear();
  sci_range.clear();
  iron_range.clear();
  iron_pene.clear();
  veto.clear();
  edge.clear();
  pdg.clear();
  mucl.clear();
  trkpe.clear();


}

//......................................................................
PMAnaSummary::PMAnaSummary(const PMAnaSummary& basicsum) :
  nhits(0), ntracks(0)
{
  Ntrack      = basicsum. Ntrack ;
  Ningtrack   = basicsum. Ningtrack ;
  clstime     = basicsum. clstime ;
  clstimecorr = basicsum. clstimecorr;
  exptime     = basicsum. exptime;
  hitcyc      = basicsum. hitcyc;
  ontime      = basicsum. ontime;
  vetowtracking  = basicsum. vetowtracking ;
  edgewtracking  = basicsum. edgewtracking ;

  x = basicsum. x;
  y = basicsum. y;
  z = basicsum. z;
  zx= basicsum. zx;
  zy= basicsum. zy;
  startxpln = basicsum. startxpln;
  startypln = basicsum. startypln;
  startxch = basicsum. startxch;
  startych = basicsum. startych;
  endxpln = basicsum. endxpln;
  endypln = basicsum. endypln;
  endxch = basicsum. endxch;
  endych = basicsum. endych;
  thetax = basicsum. thetax;
  thetay = basicsum. thetay;
  angle = basicsum. angle;
  ing_startmod = basicsum. ing_startmod;
  ing_endmod = basicsum. ing_endmod;
  ing_startpln = basicsum. ing_startpln;
  ing_endpln = basicsum. ing_endpln;
  ing_trk = basicsum. ing_trk;
  pm_stop = basicsum. pm_stop;
  ing_stop = basicsum. ing_stop;
  sci_range = basicsum. sci_range;
  iron_range = basicsum. iron_range;
  iron_pene = basicsum. iron_pene;
  veto = basicsum. veto;
  edge = basicsum. edge;
  pdg = basicsum. pdg;
  mucl = basicsum. mucl;
  trkpe = basicsum. trkpe;


  ntracks = 0;
  nhits = 0;
  

  for (int i=0; i<INGRIDHIT_MAXHITS; ++i) {
    fIngridHit[i] = TRef(NULL);
  }
  for (int i=0; i < basicsum.nhits; ++i) 
    AddIngridHit(basicsum.GetIngridHit(i));

  for (int i=0; i<INGRIDRECON_MAXTRACKS; ++i) {
    fTrack[i] = TRef(NULL);
  }
  for (int i=0; i < basicsum.ntracks; ++i) 
    AddTrack(basicsum.GetTrack(i));

}

//......................................................................



PMAnaSummary::~PMAnaSummary() 
{

}

//......................................................................


void PMAnaSummary::Clear(Option_t* option)
{
  for (int i=0; i<INGRIDHIT_MAXHITS; ++i)
    fIngridHit[i] = TRef(NULL);
  nhits = 0;
  
  for (int i=0; i<INGRIDRECON_MAXTRACKS; ++i)
    fTrack[i] = TRef(NULL);
  ntracks = 0;
}

//......................................................................

void PMAnaSummary::Print()
{

}
//......................................................................
void PMAnaSummary::AddIngridHit(IngridHitSummary* hit) 
{
  if (nhits < INGRIDHIT_MAXHITS) {
    fIngridHit[nhits] = TRef((IngridHitSummary*) hit);
    ++nhits;
  }
}
//......................................................................


IngridHitSummary* PMAnaSummary::GetIngridHit(int i) const
{ 
  return (IngridHitSummary*)fIngridHit[i].GetObject();
}

//......................................................................

void PMAnaSummary::AddTrack(IngridTrackSummary* trk) 
{
  if (ntracks < INGRIDRECON_MAXTRACKS) {
    fTrack[ntracks] = TRef((IngridTrackSummary*) trk);
    ++ntracks;
  }
}
//......................................................................


IngridTrackSummary* PMAnaSummary::GetTrack(int i) const
{ 

  return (IngridTrackSummary*)fTrack[i].GetObject();
}

//......................................................................


ClassImp(PMAnaSummary)

////////////////////////////////////////////////////////////////////////

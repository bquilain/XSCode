#include "./INGRIDEVENTSUMMARY.h"

#include <vector>
#include <algorithm>
using namespace std;

//......................................................................

IngridEventSummary::IngridEventSummary():
  fIngridSimHit(0),fIngridHit(0), 
  fSimParticle(0), fSimVertex(0),
  fBeamSummary(0), fBasicRecon(0), fPMRecon(0), fPMAna(0),
  fIngridTrack(0),f1stReduc(0),
  fNeutInfo(0)//,
{ 
    this->NewTClonesArrays();
    this->Clear("C");



}

//......................................................................

IngridEventSummary::IngridEventSummary(const IngridEventSummary& evt) 
{ 

    this->NewTClonesArrays();
    this->Clear("C");

    run      = evt.run;
    event    = evt.event;
    runmode  = evt.runmode;
    trgid    = evt.trgid;
    version  = evt.version;
    date     = evt.date;
    time     = evt.time;
    trgtime  = evt.trgtime;
  
    for(int i=0; i<23;i++)bunch_flag[i] = evt.bunch_flag[i];

    ningridsimhits = 0;
    ningridhits    = 0;
    nbasicrecons   = 0;
    npmrecons   = 0;
    npmanas   = 0;
    nbeamsummarys  = 0;
    nsimvertexes   = 0;
    nsimparticles  = 0;
    n1streducs     = 0;
    nneutinfos      = 0;

    for(int mod=0; mod<17; mod++){
      for(int cyc=0; cyc<23; cyc++){
	nidmodhits    [mod * 23 + cyc].clear();
	ningridmodhits[mod][cyc] = 0;
      }
    }

    for(int mod=0; mod<17; mod++){
      for(int cyc=0; cyc<23; cyc++){
	for(int view=0; view<2; view++){
	  nidpmmodrecons    [(mod * 23 + cyc )*2 + view].clear();
	  npmmodrecons[mod][cyc][view] = 0;
	}
      }
    }

    //ningridtracks  = 0;

    // Copy all arrays of detector summaries
    int i;
    for (i=0; i < evt.ningridsimhits;    ++i) {
      AddIngridSimHit(evt.GetIngridSimHit(i));
    }
    for (i=0; i < evt.ningridhits;       ++i) {
      IngridHitSummary*  tingsum = evt.GetIngridHit(i);
      int                    mod = tingsum -> mod;
      int                    cyc = tingsum -> cyc;
      AddIngridModHit   ( evt.GetIngridHit(i), mod, cyc);
    }
    for (i=0; i < evt.nbeamsummarys;    ++i) {
      AddBeamSummary(evt.GetBeamSummary(i));
    }
    for (i=0; i < evt.nbasicrecons;    ++i) {
      AddBasicRecon(evt.GetBasicRecon(i));
    }
    /*
    for (i=0; i < evt.npmrecons;    ++i) {
      AddPMRecon(evt.GetPMRecon(i));
    }
    */
    for (i=0; i < evt.npmrecons;    ++i) {
      PMReconSummary*  tpmsum = evt.GetPMRecon(i);
      int                    mod = tpmsum -> hitmod;
      int                    cyc = tpmsum -> hitcyc;
      int                   view = tpmsum -> view;
      AddPMModRecon   ( evt.GetPMRecon(i), mod, cyc, view);
    }
    for (i=0; i < evt.npmanas;    ++i) {
      AddPMAna(evt.GetPMAna(i));
    }
    for (i=0; i < evt.ningridtracks;    ++i) {
      AddTrack(evt.GetTrack(i));
    }
    for (i=0; i < evt.n1streducs;    ++i) {
      Add1stReduc(evt.Get1stReduc(i));
    }
    for (i=0; i < evt.nneutinfos; ++i) {
      AddNeut(evt.GetNeut(i));
    }

}

//......................................................................

IngridEventSummary::~IngridEventSummary() 
{ 
  this->Clear("C");

  if (fIngridSimHit) {
    delete fIngridSimHit;
    fIngridSimHit = 0;
  }

  if (fIngridHit) {
    delete fIngridHit;
    fIngridHit = 0;
  }

  if (fSimParticle) {
    delete fSimParticle;
      fSimParticle = 0;
  }

  if (fSimVertex) {
    delete fSimVertex;
      fSimVertex = 0;
  }

  if (fBeamSummary) {
    delete fBeamSummary;
      fBeamSummary = 0;
  }
  if (fBasicRecon) {
    delete fBasicRecon;
    fBasicRecon = 0;
  }
  if (fPMRecon) {
    delete fPMRecon;
    fPMRecon = 0;
  }
  if (fPMAna) {
    delete fPMAna;
    fPMAna = 0;
  }
  if (fIngridTrack) {
    delete fIngridTrack;
    fIngridTrack = 0;
  }
  if (f1stReduc) {
    delete f1stReduc;
    f1stReduc = 0;
  }
  if (fNeutInfo) {
    delete fNeutInfo;
    nneutinfos = 0;
  }

}

//......................................................................

void IngridEventSummary::NewTClonesArrays()  
{ 

    fIngridSimHit = new TClonesArray("IngridSimHitSummary", 1000);
    fIngridHit    = new TClonesArray("IngridHitSummary"   , 1000);
    fSimVertex    = new TClonesArray("IngridSimVertexSummary", 1000);
    fSimParticle  = new TClonesArray("IngridSimParticleSummary"   , 1000);
    fBeamSummary  = new TClonesArray("BeamInfoSummary"   , 1000);
    fBasicRecon   = new TClonesArray("IngridBasicReconSummary"   , 1000);
    fPMRecon      = new TClonesArray("PMReconSummary"   , 1000);
    fPMAna        = new TClonesArray("PMAnaSummary"   , 1000);
    fIngridTrack  = new TClonesArray("IngridTrackSummary"   , 1000);
    f1stReduc     = new TClonesArray("Ingrid1stReducSummary"   , 1000);
    fNeutInfo     = new TClonesArray("NeutInfoSummary"   , 1000);

}

//......................................................................

void IngridEventSummary::Clear(Option_t* option) 
{ 

    if (fIngridSimHit) fIngridSimHit -> Clear(option);
    if (fIngridHit)    fIngridHit    -> Clear(option);
    if (fSimParticle)  fSimParticle  -> Clear(option);
    if (fSimVertex)    fSimVertex    -> Clear(option);
    if (fBeamSummary)  fBeamSummary  -> Clear(option);
    if (fBasicRecon)   fBasicRecon   -> Clear(option);
    if (fPMRecon)      fPMRecon      -> Clear(option);
    if (fPMAna)        fPMAna        -> Clear(option);    
    if (fIngridTrack)  fIngridTrack  -> Clear(option);
    if (f1stReduc)     f1stReduc     -> Clear(option);
    if (fNeutInfo)     fNeutInfo     -> Clear(option);

    ningridsimhits  = 0;
    ningridhits     = 0;
    nsimparticles   = 0;
    nsimvertexes    = 0;
    nbeamsummarys   = 0;
    nbasicrecons    = 0;
    npmrecons       = 0;
    npmanas         = 0;
    ningridtracks   = 0;
    n1streducs      = 0;
    nneutinfos      = 0;
    for(int mod=0; mod<17; mod++){
      for(int cyc=0; cyc<23; cyc++){
	nidmodhits[mod*23+cyc].clear();
	ningridmodhits[mod][cyc] = 0;
      }
    }
    for(int mod=0; mod<17; mod++){
      for(int cyc=0; cyc<23; cyc++){
	for(int view=0; view<2; view++){
	  nidpmmodrecons    [(mod * 23 + cyc )*2 + view].clear();
	  npmmodrecons[mod][cyc][view] = 0;
	}
      }
    }


    // initialize rest of data members
    run     = 0;
    event   = 0;
    runmode = -10;
    trgid   = -10;
    version = ' ';
    date = -1;
    time = -1;
    trgtime = -1;
    for(int i=0;i<23;i++)bunch_flag[i]=false;
}

//......................................................................

void IngridEventSummary::Print()
{
    std::cout << "Event summary: " << std::endl;

}

//......................................................................


IngridSimHitSummary* IngridEventSummary::AddIngridSimHit(IngridSimHitSummary* ingridhitsum) 
{

    TClonesArray &ingridsimhit_s = *fIngridSimHit;
    
//    use memberwise copy
    new(ingridsimhit_s[ningridsimhits++]) IngridSimHitSummary(*ingridhitsum);

    return (IngridSimHitSummary*)(fIngridSimHit->At(ningridsimhits-1));

}


//......................................................................

IngridSimHitSummary* IngridEventSummary::GetIngridSimHit(int i) const
{ 

    if (i < ningridsimhits && i>=0 ) return (IngridSimHitSummary*)(fIngridSimHit->At(i));
    return 0;

}


//......................................................................

IngridHitSummary* IngridEventSummary::AddIngridHit(IngridHitSummary* ingridhitsum) 
{

    TClonesArray &ingridhit_s = *fIngridHit;
    
//    use memberwise copy
    new(ingridhit_s[ningridhits++]) IngridHitSummary(*ingridhitsum);

    return (IngridHitSummary*)(fIngridHit->At(ningridhits-1));

}



//......................................................................

IngridHitSummary* IngridEventSummary::GetIngridHit(int i) const
{ 

    if (i < ningridhits && i>=0 ) return (IngridHitSummary*)(fIngridHit->At(i));
    return 0;

}



IngridHitSummary* IngridEventSummary::AddIngridModHit(IngridHitSummary* ingridhitsum, int nmod, int ncyc) 
{

    TClonesArray &ingridhit_s = *fIngridHit;
    
    //##### push back ningridhits to hit # ID vector ######
    nidmodhits    [ nmod * 23 + ncyc ].push_back( ningridhits );
    ningridmodhits[ nmod ][ ncyc ]++;


    new(ingridhit_s[ningridhits++]) IngridHitSummary(*ingridhitsum);
    return (IngridHitSummary*)(fIngridHit->At(ningridhits-1));

}


IngridHitSummary* IngridEventSummary::GetIngridModHit(int i, int nmod, int ncyc ) const
{ 


    if (i < ningridmodhits[nmod][ncyc] && i>=0 ) return (IngridHitSummary*)(fIngridHit->At( nidmodhits[nmod * 23 + ncyc][i] ));
    return 0;


}

IngridSimVertexSummary* IngridEventSummary::AddSimVertex(IngridSimVertexSummary* simvertex) 
{

    TClonesArray &simvertex_s = *fSimVertex;
    
//    use memberwise copy
    new(simvertex_s[nsimvertexes++]) IngridSimVertexSummary(*simvertex);

    return (IngridSimVertexSummary*)(fSimVertex->At(nsimvertexes-1));
}


//......................................................................

IngridSimVertexSummary* IngridEventSummary::GetSimVertex(int i) const
{ 

    if (i < nsimvertexes && i>=0 ) return (IngridSimVertexSummary*)(fSimVertex->At(i));
    return 0;

}


IngridSimParticleSummary* IngridEventSummary::AddSimParticle(IngridSimParticleSummary* simparticle) 
{

    TClonesArray &simparticle_s = *fSimParticle;
    
//    use memberwise copy
    new(simparticle_s[nsimparticles++]) IngridSimParticleSummary(*simparticle);

    return (IngridSimParticleSummary*)(fSimParticle->At(nsimparticles-1));
}


//......................................................................

IngridSimParticleSummary* IngridEventSummary::GetSimParticle(int i) const
{ 

    if (i < nsimparticles && i>=0 ) return (IngridSimParticleSummary*)(fSimParticle->At(i));
    return 0;

}
//......................................................................
BeamInfoSummary* IngridEventSummary::AddBeamSummary(BeamInfoSummary* beamsummary) 
{

  TClonesArray &beamsummary_s = *fBeamSummary;
    
//    use memberwise copy
    new(beamsummary_s[nbeamsummarys++]) BeamInfoSummary(*beamsummary);
    return (BeamInfoSummary*)(fBeamSummary->At(nbeamsummarys-1));
}


//......................................................................

BeamInfoSummary* IngridEventSummary::GetBeamSummary(int i) const
{ 

    if (i < nbeamsummarys && i>=0 ) return (BeamInfoSummary*)(fBeamSummary->At(i));
    return 0;

}

//......................................................................
IngridBasicReconSummary* IngridEventSummary::AddBasicRecon(IngridBasicReconSummary* basicrecon) 
{

  TClonesArray &basicrecon_s = *fBasicRecon;
    
//    use memberwise copy
    new(basicrecon_s[nbasicrecons++]) IngridBasicReconSummary(*basicrecon);
    return (IngridBasicReconSummary*)(fBasicRecon->At(nbasicrecons-1));
}


//......................................................................

IngridBasicReconSummary* IngridEventSummary::GetBasicRecon(int i) const
{ 

    if (i < nbasicrecons && i>=0 ) return (IngridBasicReconSummary*)(fBasicRecon->At(i));
    return 0;

}




//......................................................................
PMReconSummary* IngridEventSummary::AddPMRecon(PMReconSummary* pmrecon) 
{

  TClonesArray &pmrecon_s = *fPMRecon;
    
//    use memberwise copy
    new(pmrecon_s[npmrecons++]) PMReconSummary(*pmrecon);
    return (PMReconSummary*)(fPMRecon->At(npmrecons-1));
}


//......................................................................

PMReconSummary* IngridEventSummary::GetPMRecon(int i) const
{ 

    if (i < npmrecons && i>=0 ) return (PMReconSummary*)(fPMRecon->At(i));
    return 0;

}




PMReconSummary* IngridEventSummary::AddPMModRecon(PMReconSummary* pmreconsum, int nmod, int ncyc, int nview) 
{

    TClonesArray &pmrecon_s = *fPMRecon;
    
    //##### push back ningridhits to hit # ID vector ######
    nidpmmodrecons    [ (nmod * 23 + ncyc )*2 + nview ].push_back( npmrecons );
    npmmodrecons[ nmod ][ ncyc ][ nview ]++;


    new(pmrecon_s[npmrecons++]) PMReconSummary(*pmreconsum);
    return (PMReconSummary*)(fPMRecon->At(npmrecons-1));

}


PMReconSummary* IngridEventSummary::GetPMModRecon(int i, int nmod, int ncyc, int nview ) const
{ 


  if (i < npmmodrecons[nmod][ncyc][nview] && i>=0 ) return (PMReconSummary*)(fPMRecon->At( nidpmmodrecons[(nmod * 23 + ncyc)*2 + nview][i] ));
  return 0;

}




PMAnaSummary* IngridEventSummary::AddPMAna(PMAnaSummary* pmana) 
{
  TClonesArray &pmana_s = *fPMAna;
  //use memberwise copy
  new(pmana_s[npmanas++]) PMAnaSummary(*pmana);
  return (PMAnaSummary*)(fPMAna->At(npmanas-1));
}


//......................................................................

PMAnaSummary* IngridEventSummary::GetPMAna(int i) const
{ 
  if (i < npmanas && i>=0 ) return (PMAnaSummary*)(fPMAna->At(i));
  return 0;
}



//......................................................................
Ingrid1stReducSummary* IngridEventSummary::Add1stReduc(Ingrid1stReducSummary* reduc) 
{
  TClonesArray &reduc_s = *f1stReduc;
    
  //    use memberwise copy
  new(reduc_s[n1streducs++]) Ingrid1stReducSummary(*reduc);
  return (Ingrid1stReducSummary*)(f1stReduc->At(n1streducs-1));
}


//......................................................................

Ingrid1stReducSummary* IngridEventSummary::Get1stReduc(int i) const
{ 

    if (i < n1streducs && i>=0 ) return (Ingrid1stReducSummary*)(f1stReduc->At(i));
    return 0;

}


//......................................................................
NeutInfoSummary* IngridEventSummary::AddNeut(NeutInfoSummary* neut) 
{
  TClonesArray &neutinfo_s = *fNeutInfo;

//    use memberwise copy
  new(neutinfo_s[nneutinfos++]) NeutInfoSummary(*neut);
  return (NeutInfoSummary*)(fNeutInfo->At(nneutinfos-1));
}


//......................................................................

NeutInfoSummary* IngridEventSummary::GetNeut(int i) const
{ 
    if (i < nneutinfos && i>=0 ) return (NeutInfoSummary*)(fNeutInfo->At(i));
    return 0;
}



IngridTrackSummary* IngridEventSummary::AddTrack(IngridTrackSummary* trk)
{
    TClonesArray &track_s = *fIngridTrack;
    
//    use memberwise copy
    new(track_s[ningridtracks++]) IngridTrackSummary(*trk);
    return (IngridTrackSummary*)(fIngridTrack->At(ningridtracks-1));

}


//......................................................................

IngridTrackSummary* IngridEventSummary::GetTrack(int i) const
{ 
    if (i < ningridtracks && i>=0 ) return (IngridTrackSummary*)(fIngridTrack->At(i));
    return 0;

 }



ClassImp(IngridEventSummary)


////////////////////////////////////////////////////////////////////////

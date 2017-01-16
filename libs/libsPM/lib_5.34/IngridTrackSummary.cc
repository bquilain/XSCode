#include "IngridTrackSummary.h"

//......................................................................

IngridTrackSummary::IngridTrackSummary() :
    nhits(0),nsimparticles(0),nsimemshowers(0)
{
    for (unsigned int i=0; i<3; ++i) {
        vtxi[i] = -1.e5;
        vtxf[i] = -1.e5;
    }
    length = -1.e5;
    ekin = -1.e5;
    tx = -1.e5;
    ty = -1.e5;
    etx = -1.e5;
    ety = -1.e5;
    ex0 = -1.e5;
    ey0 = -1.e5;
    covx = -1.e5;
    covy = -1.e5;
    chi2x = -1.e5;
    chi2y = -1.e5;
    btheta = -1.e5;
    bphi = -1.e5;

    mrdhitid[0] = -1;
    mrdhitid[1] = -1;
    mucl = -1.e5;
    vpe  = -1.e5;
    view = -1;
}

//......................................................................

IngridTrackSummary::IngridTrackSummary(const IngridTrackSummary& trk) :
    nhits(0),nsimparticles(0),nsimemshowers(0)
{ 
    for (unsigned int i=0; i<3; ++i) {
        vtxi[i] = trk.vtxi[i];
        vtxf[i] = trk.vtxf[i];
    }
    length = trk.length;
    ekin = trk.ekin;
    tx = trk.tx;
    ty = trk.ty;
    etx = trk.etx;
    ety = trk.ety;
    ex0 = trk.ex0;
    ey0 = trk.ey0;
    covx = trk.covx;
    covy = trk.covy;
    chi2x = trk.chi2x;
    chi2y = trk.chi2y;
    btheta = trk.btheta;
    bphi = trk.bphi;
    vpe  = trk.vpe;
    view = trk.view;

    mrdhitid[0] = trk.mrdhitid[0];
    mrdhitid[1] = trk.mrdhitid[1];


    for (int i=0; i<TRACK_MAXHITS; ++i)
        fIngridHit[i] = TRef(NULL);
    for (int i=0; i < trk.nhits; ++i) 
        AddIngridHit(trk.GetIngridHit(i));

    for (int i=0; i<TRACK_MAXSIMPARTICLES; ++i)
        fSimParticle[i] = TRef(NULL);
    for (int i=0; i < trk.nsimparticles; ++i) 
        AddSimParticle(trk.GetSimParticle(i));



}

//......................................................................

IngridTrackSummary::~IngridTrackSummary() 
{ 
}


//......................................................................


void IngridTrackSummary::Clear(Option_t* option)
{
    for (int i=0; i<TRACK_MAXHITS; ++i)
        fIngridHit[i] = TRef(NULL);
    nhits = 0;

    for (int i=0; i<TRACK_MAXSIMPARTICLES; ++i)
        fSimParticle[i] = TRef(NULL);
    nsimparticles = 0;


}


//......................................................................

void IngridTrackSummary::Print()
{
    std::cout << "number of hits = " << nhits << ", length = " << 
        length << std::endl;
}

//......................................................................

void IngridTrackSummary::AddIngridHit(IngridHitSummary* sbhitsum) 
{
    if (nhits < TRACK_MAXHITS) {
        fIngridHit[nhits] = TRef((IngridHitSummary*) sbhitsum);
        ++nhits;
    }
}

//......................................................................


IngridHitSummary* IngridTrackSummary::GetIngridHit(int i) const
{ 
    return (IngridHitSummary*)fIngridHit[i].GetObject();
}

//......................................................................


void IngridTrackSummary::AddSimParticle(IngridSimParticleSummary* sbsimpart) 
{
    if (nsimparticles < TRACK_MAXSIMPARTICLES) {
        fSimParticle[nsimparticles] = TRef((IngridSimParticleSummary*) sbsimpart);
        ++nsimparticles;
    }
}

//......................................................................


IngridSimParticleSummary* IngridTrackSummary::GetSimParticle(int i) const
{ 
    return (IngridSimParticleSummary*)fSimParticle[i].GetObject();
}

//......................................................................




ClassImp(IngridTrackSummary)

////////////////////////////////////////////////////////////////////////


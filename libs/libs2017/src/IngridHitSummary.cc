#include "IngridHitSummary.h"

//......................................................................

IngridHitSummary::IngridHitSummary():
  nsimhits(0)
{ 
    view     = -1;
    mod      = -1;
    pln      = -1;
    ch       = -1;
    adc      = -1;
    loadc    = -1;
    cyc      = -1;
    pe       = -1.e5;
    lope     = -1.e5;
    pecorr   = -1.e5;
    vise     = -1.e5;
    visecorr = -1.e5;
    tdc      = -1.e5;
    time     = -1.e5;
    tnearhit = -1.e5;
    timecorr = -1.e5;
    xy       = -1.e5;
    z        = -1.e5;
    addbasicrecon  = 0;
    dummy          = 0;
    gocosmic = false;
    hitcosmic= false;
    isohit   = false;
    gridcell_id_x1=0; //for cross talk study 2016/1/13
    gridcell_id_x2=0; //for cross talk study 2016/1/13
    gridcell_id_y1=0; //for cross talk study 2016/1/13
    gridcell_id_y2=0; //for cross talk study 2016/1/13
    pe_cross=0; //for cross talk study 2016/1/13
}

//......................................................................

IngridHitSummary::IngridHitSummary(const IngridHitSummary& hit) :
  nsimhits(0)
{ 

    mod      = hit.mod;
    view     = hit.view;
    pln      = hit.pln;
    ch       = hit.ch;
    cyc      = hit.cyc;
    adc      = hit.adc;
    loadc    = hit.loadc;
    pe       = hit.pe;
    lope     = hit.lope;
    pecorr   = hit.pecorr;
    vise     = hit.vise;
    visecorr = hit.visecorr;
    tdc      = hit.tdc;
    time     = hit.time;
    tnearhit = hit.tnearhit;
    timecorr = hit.timecorr;
    xy       = hit.xy;
    z        = hit.z;
    addbasicrecon  = hit.addbasicrecon;
    dummy          = hit.dummy;
    gocosmic       = hit.gocosmic;
    hitcosmic      = hit.hitcosmic;
    isohit         = hit.isohit;
    gridcell_id_x1 = hit.gridcell_id_x1; //for cross talk study 2016/1/13
    gridcell_id_x2 = hit.gridcell_id_x2; //for cross talk study 2016/1/13
    gridcell_id_y1 = hit.gridcell_id_y1; //for cross talk study 2016/1/13
    gridcell_id_y2 = hit.gridcell_id_y2; //for cross talk study 2016/1/13
    pe_cross=hit.pe_cross; //for cross talk study 2016/1/13

    for (int i=0; i<INGRIDHIT_MAXSIMHITS; ++i) {
        fIngridSimHit[i] = TRef(NULL);
    }

    for (int i=0; i < hit.nsimhits; ++i) 
        AddIngridSimHit(hit.GetIngridSimHit(i));
}

//......................................................................

IngridHitSummary::~IngridHitSummary() 
{
}

//......................................................................


void IngridHitSummary::Clear(Option_t* option)
{
    for (int i=0; i<INGRIDHIT_MAXSIMHITS; ++i)
        fIngridSimHit[i] = TRef(NULL);
    nsimhits = 0;
}

//......................................................................

void IngridHitSummary::Print()
{
  std::cout << "-------------------" <<std::endl
	    << "Module# = " << mod   <<std::endl
	    << "View    = " << view  <<std::endl
	    << "Plane # = " << pln   <<std::endl
	    << "Ch. #   = " << ch    <<std::endl
	    << "ADC     = " << adc   <<std::endl
	    << "TDC     = " << tdc   <<std::endl
	    << "xy      = " << xy    <<std::endl
	    << "z       = " << z     <<std::endl
	    << std::endl;
}


//......................................................................

void IngridHitSummary::AddIngridSimHit(IngridSimHitSummary* sbhitsum) 
{
    if (nsimhits < INGRIDHIT_MAXSIMHITS) {
        fIngridSimHit[nsimhits] = TRef((IngridSimHitSummary*) sbhitsum);
        ++nsimhits;
    }
}
//......................................................................


IngridSimHitSummary* IngridHitSummary::GetIngridSimHit(int i) const
{ 
    return (IngridSimHitSummary*)fIngridSimHit[i].GetObject();
}

//......................................................................


ClassImp(IngridHitSummary)

////////////////////////////////////////////////////////////////////////

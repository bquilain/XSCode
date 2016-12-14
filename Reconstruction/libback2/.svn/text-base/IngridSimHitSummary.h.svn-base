#ifndef INGRIDSIMHITSUMMARY_H
#define INGRIDSIMHITSUMMARY_H
#include <iostream>
#include "TObject.h"


//......................................................................

class IngridSimHitSummary : public TObject
{
public:
    IngridSimHitSummary();
    virtual ~IngridSimHitSummary();
    
    void Clear   (Option_t* option="");
    void Print();

public:
    
    float edeposit;             // Energy deposit (MeV)
    int   trackid;              // ID of GEANT4 track responsible for hit
    int   pdg;                  // PDG particle ID of GEANT4 track responsible for hit

private:

    ClassDef(IngridSimHitSummary, 2) // Ingrid Hit Sim Summary
};

#endif // INGRIDSIMHITSUMMARY_H
////////////////////////////////////////////////////////////////////////

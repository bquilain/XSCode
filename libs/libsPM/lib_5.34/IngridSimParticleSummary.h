#ifndef INGRIDSIMPARTICLESUMMARY_H
#define INGRIDSIMPARTICLESUMMARY_H
#include <iostream>
#include "TObject.h"

//......................................................................

class IngridSimParticleSummary : public TObject
{
public:
    IngridSimParticleSummary();
    virtual ~IngridSimParticleSummary();
    
    void Clear   (Option_t* option="");
    void Print();
    
public:

    int trackid;             // GEANT4 track ID number
    int parentid;            // GEANT4 parent track ID number
    int pdg;                 // PDG particle numbering scheme
    float momentum[4];       // particle's 4-momentum (GeV)
    float ipos[4];           // particle's initial position/time (cm, ns)
    float fpos[4];           // particle's final position/time (cm, ns)
    int iposflag;            // particle's initial position flag
    int fposflag;            // particle's final position flag
    float dir[3];            // particle's direction
    float length;
private:

    ClassDef(IngridSimParticleSummary, 3) //  Simulation (detector mc) particle Summary
        };

#endif // IngridSIMPARTICLESUMMARY_H
////////////////////////////////////////////////////////////////////////

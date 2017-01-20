#ifndef INGRID_RESPONSE
#define INGRID_RESPONSE

#include "G4ThreeVector.hh"
#include "G4EmCalculator.hh"
#include "G4EmSaturation.hh"
#include "G4Track.hh"

class IngridResponse {
public:
  IngridResponse();
  IngridResponse(G4double);
  ~IngridResponse();

  double ApplyConversion(G4double* edep);
  void ApplyScintiResponse(G4double* edep, G4Track* aTrack);
  void ApplyLightCollection(G4double* edep, G4int mod, G4int view, G4ThreeVector pos);//added
  void ApplyFiberResponse(G4double* edep, G4double* time, G4int view, G4ThreeVector pos);
  void ApplyFiberResponseV(G4double* edep, G4double* time, G4int pln, G4ThreeVector pos);
  void ApplyMPPCResponse(G4double edep, G4int mod, G4double* pe);
  void ApplyADCResponse(G4double* pe, G4double* lope, G4int* adc, G4int* loadc);
  void ApplyTDCResponse(G4double time, G4int* tdc);

  void ApplyScintiResponse2(G4double* edep, G4double* steplength, G4Track* aTrack);
private:
  G4EmCalculator emcal;
  G4double CBIRKS;

  void BirksSaturation(G4double* edeposit, G4Track* aTrack);
  void BirksSaturation2(G4double* edeposit, G4double * steplength, G4Track* aTrack);


};

#endif

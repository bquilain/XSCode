#ifndef IngridPrimaryGeneratorAction_h
#define IngridPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"

#include "IngridRunAction.hh"
#include "IngridEventAction.hh"

#include "Neut.hh"

class G4ParticleGun;
class G4Event;

class IngridPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  //IngridPrimaryGeneratorAction(Neut*, IngridRunAction*, IngridEventAction*, int, int);
  IngridPrimaryGeneratorAction(Neut*, IngridRunAction*, IngridEventAction*, int, int, int);
  ~IngridPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void SetParticleGun(bool, int);

private:
  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;

  int module_mode;
  int neutrino_flavor;

  IngridRunAction* runaction;
  IngridEventAction* eventaction;

  // for option -p : propagate only a single type of particle
  bool isParticleGun;
  int particleGun_pdg;
  
  Neut *neut_fe;
};

#endif

#ifndef IngridHLayerHit_h
#define IngridHLayerHit_h 1
 
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class IngridHLayerHit : public G4VHit 
{
public:
  //IngridHLayerHit(G4int, G4int, G4int, G4double, G4double, const G4ThreeVector&, G4double);
  IngridHLayerHit(G4int, G4int, G4int, G4double, G4double, G4ThreeVector, G4double);
  IngridHLayerHit(G4int, G4double, G4int, const G4ThreeVector&);
  IngridHLayerHit(G4int, G4double, G4int);
  IngridHLayerHit(G4int, G4double);
  
  ~IngridHLayerHit();
  IngridHLayerHit(const IngridHLayerHit &right);
  const IngridHLayerHit& operator=(const IngridHLayerHit &right);
  G4int operator==(const IngridHLayerHit &right) const;
  G4int CompareID(const IngridHLayerHit right) ;
  G4int CompareP(const IngridHLayerHit right);
  G4int isFaster(const IngridHLayerHit right);
  G4int LargerEdep(const IngridHLayerHit right);

  // new/delete operators
  inline void *operator new(size_t);
  inline void operator delete(void  *aHit);

  void Draw();
  void Print();
  void Print_Hit();

private:
    G4int detID;
    G4int trackID;
    G4int mod;
    G4int view;
    G4int pln;
    G4int ch;
    G4double edep;
    G4double edep_q;
    G4double pe;
    G4double lope;
    G4int Particle;
    G4ThreeVector position;
    G4ThreeVector posinmod;
    G4int eventID;
    G4double time;
    G4double delay_time;
    G4int gridcell_id_x1;//added by koga 2016/1/7
    G4int gridcell_id_x2;//added by koga 2016/1/7
    G4int gridcell_id_y1;//added by koga 2016/1/7
    G4int gridcell_id_y2;//added by koga 2016/1/7
  
  //set/get functions 
public:
		//
    inline void AddEdep(G4double e, G4double eq) { edep += e; edep_q += eq; }
    inline void SetEventID(G4int event_ID) { eventID = event_ID;}
    inline void SetPE(G4double p) { pe = p; }
    inline void SetLOPE(G4double p) { lope = p; }
    inline void SetTime(G4double t) { time = t; }
    inline void SetParticle(G4int p) { Particle = p; }
    inline void SetDelayTime(G4double dt) { delay_time = dt; }
    inline void SetCellidX1(G4int id) { gridcell_id_x1 = id; }//added by koga 2016/1/7
    inline void SetCellidX2(G4int id) { gridcell_id_x2 = id; }//added by koga 2016/1/7
    inline void SetCellidY1(G4int id) { gridcell_id_y1 = id; }//added by koga 2016/1/7
    inline void SetCellidY2(G4int id) { gridcell_id_y2 = id; }//added by koga 2016/1/7

		//
    inline G4int GetDetID() { return detID; }// ID of scintillator
    inline G4int GetTrackID() { return trackID; }
    inline G4int GetParticle() { return Particle; }
    inline G4double GetEdep() { return edep; }
    inline G4double GetEdepQ() { return edep_q; }
    inline G4ThreeVector GetPosition() { return position; }
    inline G4ThreeVector GetPosInMod() { return posinmod; }
    inline G4double GetTime() { return time; }
    inline G4double GetDelayTime() { return delay_time; }
    inline G4int GetEventID() { return eventID; }
    inline G4int GetMod() { return mod; }
    inline G4int GetView() { return view; }
    inline G4int GetPln() { return pln; }
    inline G4int GetCh() { return ch; }
    inline G4double GetPE() { return pe; }
    inline G4double GetLOPE() { return lope; }
    inline G4int GetCellidX1() { return gridcell_id_x1; }//added by koga 2016/1/7
    inline G4int GetCellidX2() { return gridcell_id_x2; }//added by koga 2016/1/7
    inline G4int GetCellidY1() { return gridcell_id_y1; }//added by koga 2016/1/7
    inline G4int GetCellidY2() { return gridcell_id_y2; }//added by koga 2016/1/7

};


// ====================================================================
// inline functions
// ====================================================================

// externally instanciated.
typedef G4THitsCollection<IngridHLayerHit> IngridHLayerHitsCollection;

extern G4Allocator<IngridHLayerHit> IngridHLayerHitAllocator; 

inline void* IngridHLayerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) IngridHLayerHitAllocator.MallocSingle();
  return aHit;
}

inline void IngridHLayerHit::operator delete(void *aHit)
{
  IngridHLayerHitAllocator.FreeSingle((IngridHLayerHit*) aHit);
}

#endif

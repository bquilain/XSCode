//#include "cfortran/cfortran.h"
//#include "cfortran/hbook.h"
//comment out by akira090918

#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4EmSaturation.hh"
#include "G4EmCalculator.hh"
#include "IngridHLayerSD.hh"
#include "IngridResponse.hh"

//#define DEBUG

IngridHLayerSD::IngridHLayerSD(const G4String& name, double cbirks)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("hlayerHitCollection");

  ingresp = new IngridResponse(cbirks);
}


IngridHLayerSD::~IngridHLayerSD()
{
    if(ingresp!=NULL) delete ingresp;
}


void IngridHLayerSD::Initialize(G4HCofThisEvent* HCTE)
{
  hlayerHitCollection = 
    new IngridHLayerHitsCollection(SensitiveDetectorName,collectionName[0]);
  TotalhlayerDep = 0.;
  
  static int HCID = -1;
  if(HCID<0) HCID = GetCollectionID(0); 
  HCTE->AddHitsCollection(HCID,hlayerHitCollection);
}

//
G4bool IngridHLayerSD::ProcessHits(G4Step* aStep, 
				G4TouchableHistory* ROhist)
{
  G4EmSaturation emsat;
  G4EmCalculator emcal2;

  // only when there is energy deposit in a sensitive detector,
  // create a new hit.   
  G4Track* track = aStep->GetTrack();
  G4double TrackEnergy = track->GetKineticEnergy()/MeV;
  G4int PDG = track->GetDefinition()->GetPDGEncoding();
  G4double hittime = aStep->GetPreStepPoint()->GetGlobalTime();
  G4double edep = aStep->GetTotalEnergyDeposit()/MeV;
  G4double length = aStep->GetStepLength()/cm;
  //G4double edep2 = emsat.G4EmSaturation::VisibleEnergyDeposition(aStep);

  G4double              kineticE = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4Material*           material = track->GetMaterial();
  G4double edep3=0;
  //if(particle->GetPDGCharge()!=0) edep3 = emcal2.GetDEDX(kineticE, particle, material)/(MeV/cm);
  /*
  if((PDG==2212||PDG==2112) && TrackEnergy<5){
    track->SetTrackStatus(fStopAndKill);
    return true;
  }
  
  if(PDG==2112 && TrackEnergy<5) {
    track->SetTrackStatus(fStopAndKill);
    //G4cout<<"Track is Killed="<<track->GetTrackID()<<G4endl;
    // track->SetTrackStatus(fKillTrackAndSecondaries);
    return true;
  }
  */

 if(edep==0.) return false;
 //if(TrackEnergy==0 && PDG==11) return false;
 //if(TrackEnergy<5) edep=0;
 
  TotalhlayerDep += edep;

  // volume information must be extracted from Touchable of "PreStepPoint"
  const G4VTouchable* Touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int detID = Touchable->GetVolume(0)->GetCopyNo();
  G4int mod = detID/264;
  if(mod>=16) mod=16;
  G4int trackID = track->GetTrackID();
  G4ThreeVector hitPos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector hitPosIni;
  if(aStep->IsFirstStepInVolume()==true) hitPosIni=aStep->GetPreStepPoint()->GetPosition();
  G4double LengthMaterial = std::sqrt(pow(hitPos.getX()-hitPosIni.getX(),2) + pow(hitPos.getY()-hitPosIni.getY(),2) + pow(hitPos.getZ()-hitPosIni.getZ(),2));


  /*#ifdef DEBUG
  G4cout<<"Track Energy="<<TrackEnergy/MeV<<"MeV, Particle="<<PDG<<G4endl;
  G4cout<<"Track ID="<<trackID<<", ParentID="<<track->GetParentID()<<", Hit timing="<<hittime<<", energy left="<<edep<<", step length="<<length<<G4endl;
  #endif
*/

#ifdef DEBUG
  G4cout<<"Track ID="<<track->GetTrackID()<<",Particle type="<<PDG<<", ParentID="<<track->GetParentID()<<G4endl;
  G4cout<<"Kinetic Energy ="<<TrackEnergy<<", Track Energy="<<TrackEnergy/MeV<<"MeV, Length in Material = "<<LengthMaterial/mm<<" mm"<<G4endl;
  G4cout<<", energy deposit="<<edep<<", step length="<<length<<G4endl;
  G4cout<<"Hit Position="<<hitPos.getX()<<", "<<hitPos.getY()<<", "<<hitPos.getZ()<<G4endl<<G4endl;
  G4cout<<"Module="<<detID/264<<G4endl;
  //G4cout<<"edep with attenuation="<<edep2<<", test of kinetic="<<kineticE<<", last dedx="<<edep3<<G4endl;
#endif

      //if(track->GetTrackLength()/cm<.1 && TrackEnergy<2) return true;
  //if(track->GetTrackLength()/cm<.05 && (PDG==2112||PDG>10000)) return false;

	//apply quenching effect
	G4double edep_q = edep;
	//ingresp->ApplyScintiResponse(&edep_q,track);
	ingresp->ApplyScintiResponse2(&edep_q,&length,track);
	//edep_q=edep2;
	 //
	IngridHLayerHit* aHit 
		= new IngridHLayerHit(detID,PDG,trackID,edep,edep_q,hitPos,hittime);
    
	IngridHLayerHit* bHit;

	for(int k=0;k<hlayerHitCollection->entries();k++){
		bHit = (*hlayerHitCollection)[k];

		if(bHit->CompareID(*aHit)){
		  bHit->AddEdep(edep,edep_q);
#ifdef DEBUG
		  //		  G4cout <<"E Ini="<<bHit->GetEdep()<<", Timing="<<bHit->GetTime()<<", Particle Type="<<bHit->GetParticle()<< G4endl;
		  //		  G4cout <<"E Added="<<aHit->GetEdep()<<", Timing="<<aHit->GetTime()<<", Particle Type="<<aHit->GetParticle()<< G4endl<<G4endl;
#endif
			if(bHit->isFaster(*aHit)) { 
			  bHit->SetTime(aHit->GetTime()); 
			}
			if(bHit->LargerEdep(*aHit)) { 
			  bHit->SetParticle(aHit->GetParticle()); 
			}
			return true;
		}
	}

	hlayerHitCollection->insert( aHit );
#ifdef DEBUG
	G4cout<<G4endl<<G4endl;
#endif

	return true;

}

void IngridHLayerSD::EndOfEvent(G4HCofThisEvent* HCTE)
{
#ifdef DEBUG
  //G4cout << "hlayer hits : " << fHIT << "\n"; // number of Hits
#endif
  IngridHLayerHit *cHit;

  G4double edep_tmp;
  G4double time_tmp;
  G4ThreeVector posinmod;
  G4int mod;
  G4int view;
  G4int adc;
  G4int loadc;
  G4double pe;
  G4double lope;

	//
  // apply ingrid response
  for(G4int k=0;k<hlayerHitCollection->entries();k++) {
		cHit = (*hlayerHitCollection)[k];

		//
    edep_tmp = cHit->GetEdepQ();
    time_tmp = cHit->GetTime();
    posinmod = cHit->GetPosInMod();
    mod = cHit->GetMod();
    view = cHit->GetView();


    //apply light collection
    ingresp->ApplyLightCollection(&edep_tmp,mod,view,posinmod);

    //apply fiber attenuation
    ingresp->ApplyFiberResponse(&edep_tmp,&time_tmp,view,posinmod);

    //convert edep -> p.e. &  p.e. statistic & cross-talk & after-pulse
    ingresp->ApplyMPPCResponse(edep_tmp,mod,&pe);

    //apply ADC responce
    ingresp->ApplyADCResponse(&pe,&lope,&adc,&loadc);

    //fill variable to hitcollection
    cHit->SetPE(pe);
    cHit->SetLOPE(lope);
    cHit->SetDelayTime(time_tmp);
  }

}

void IngridHLayerSD::DrawAll()
{
  for(G4int k=0; k < hlayerHitCollection->entries(); k++)
   (*hlayerHitCollection)[k]->Draw(); 
}

void IngridHLayerSD::PrintAll()
{
   for(G4int k=0; k < hlayerHitCollection->entries(); k++)
     (*hlayerHitCollection)[k]->Print(); 
   //hlayerHitCollection-> PrintAllHits();
}


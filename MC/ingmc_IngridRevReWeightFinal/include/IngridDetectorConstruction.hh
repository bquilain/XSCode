#ifndef IngridDetectorConstruction_H
#define IngridDetectorConstruction_H 1

#include "IngridHLayerSD.hh"
#include "IngridVetoSD.hh"
#include "IngridVLayerSD.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

#include "G4VUserDetectorConstruction.hh"

#define HALL_Z 1.0
#define HALL_XY 1.5
#define nSciVertex 12
#define nSciVertex2 8

class IngridDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  void DefineMaterial();
  void DefineSpace();
  void DefineStructures();

  G4Material *Air, *Fe, *Scinti, *Water, *Concrete, *ABS, *Paint;
  G4double WorldX, WorldY, WorldZ;
  G4double HallDirtRadiusMin, HallDirtRadiusMax, HallDirtHeight,
           HallDirtSPhi, HallDirtDPhi; //## Dirt Dimension
  G4double HallDirtX, HallDirtY, HallDirtZ; //## Dirt position
  G4double HallX, HallZ; //## center position of Hall

  G4double HorizonX,HorizonY, HorizonZ, PHorizonZ;
  G4double ModuleX, ModuleY, ModuleZ;
  G4double ModuleSpace, ModuleStart;
  G4double iron_z,iron_xy, Gap, iron_start;
  G4double Plastic_z,Plastic_xy;
  G4double scibar_x, scibar_y,scibar_z, scibar_start, scibar_xy_start,Pscibar_start;
  G4double Lveto_x, Lveto_y, Lveto_z, veto_start, Pveto_start; 
  G4double Sveto_x, Sveto_y, Sveto_z;
  G4double Niron_start, Niron_z, Nscibar_start, NGap, NModuleSpace;
  G4double SciVertex_x[nSciVertex], SciVertex_y[nSciVertex];
  G4double SciPaint_x[nSciVertex], SciPaint_y[nSciVertex];
  G4double SciVertex2_x[nSciVertex], SciVertex2_y[nSciVertex];
  G4double SciVertex3_x[nSciVertex], SciVertex3_y[nSciVertex];



  //for Proton Module added by kikawa
  G4double PModuleZ;
  G4double scibar2_y,scibar2_z;//scibar_start, scibar_xy_start;
  G4double PLveto_x, PLveto_y, PLveto_z;
  G4double PSveto_x, PSveto_y, PSveto_z;
  G4double dist_first,dist;


  G4double CBIRKS;
  IngridVetoSD* avetoSD;
  IngridHLayerSD* ahlayerSD;  
  IngridVLayerSD* avlayerSD; 
 
  int flag;
  int mode;

  IngridDetectorConstruction(int, G4double = 0.0208);
  //IngridDetectorConstruction();
  ~IngridDetectorConstruction();
  
  
  G4VPhysicalVolume* Construct();
 
};

#endif

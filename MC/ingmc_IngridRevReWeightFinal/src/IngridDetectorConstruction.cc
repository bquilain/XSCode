#include <stdio.h>
#include <iostream>
#include <vector>

#include "IngridDetectorConstruction.hh"
#include "IngridHLayerSD.hh"
#include "IngridVLayerSD.hh"
#include "IngridVetoSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "globals.hh"

#define PROTON 2
#define HORIZONTAL 3
#define VERTICAL 4

//#define MOD3

IngridDetectorConstruction::IngridDetectorConstruction(int MODE, G4double cbirks)
{
    mode = MODE;
    CBIRKS=cbirks;
}

IngridDetectorConstruction::~IngridDetectorConstruction()
{
}


G4VPhysicalVolume* IngridDetectorConstruction::Construct()
{

  //Initialise materials
  DefineMaterial();

  //Initialise volumes
  DefineSpace();

  //Initialise detector elements parameters
  DefineStructures();

  //
  int start_mod,stop_mod;
  if(mode==HORIZONTAL) { start_mod=0; stop_mod=7; }
  else if(mode==VERTICAL) { start_mod=7; stop_mod=14; }
  else if(mode==PROTON) { start_mod=0; stop_mod=17; }
  else { start_mod =0; stop_mod=0; }

  //
	G4RotationMatrix *xrot = new G4RotationMatrix(G4ThreeVector(1,0,0),90.*degree);
	G4RotationMatrix *yrot = new G4RotationMatrix(G4ThreeVector(0,1,0),90.*degree);


  // 7 modules
  // Each module:
  // 9 "Iron" blocks (120x120x6.5 cm3)
  // 1 air block (120x120x6.5 cm3)
  // 11 scintillator planes = 11x2x24 scintis (120x5x1 cm3)
  // 
  // 3 veto planes per module, 4 for the last one, each plane: 22 scintis                
  //
  // Distance between two modules  D = 30cm
  
  //World volume
  G4Box* experimentalHall_box = new G4Box("Hall",WorldX,WorldY,WorldZ);
  G4LogicalVolume* worldLV = new G4LogicalVolume(experimentalHall_box, Air, "hall_log",0,0,0);
  G4VPhysicalVolume* worldPV    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),worldLV,"hall",0,false,0);

  //Hall Dirt volume
  //### added by M.Otani 2010/08/02
  G4Tubs* HallDirtSLD = 
    new G4Tubs("halldirt", HallDirtRadiusMin, HallDirtRadiusMax, HallDirtHeight, HallDirtSPhi, HallDirtDPhi);
  G4LogicalVolume* HallDirtLV = new G4LogicalVolume(HallDirtSLD, Concrete, "HallDirt",0,0,0);
  G4VPhysicalVolume* HallDirtPV  = 
    new G4PVPlacement(xrot,G4ThreeVector(HallDirtX,HallDirtY,HallDirtZ),HallDirtLV,"HallDirt",worldLV,false,0);

  G4VisAttributes* HallDirtVisAtt = new G4VisAttributes(G4Color(1.,1.,1.));
  HallDirtVisAtt->SetForceWireframe(true);
  HallDirtLV->SetVisAttributes(HallDirtVisAtt);

  //Horizontal volume (for the 7 horizontal modules)
  G4Box* horizontalHall_box = new G4Box("HorizontalHall",HorizonX, HorizonY, HorizonZ);
  G4LogicalVolume* horizontalLV =new G4LogicalVolume(horizontalHall_box, Air, "horizon_log",0,0,0);
  //G4VPhysicalVolume* horizontalPV  = new G4PVPlacement(0,G4ThreeVector(-0.001863*m, 0.00743*m, 0.00044*m),horizontalLV,"horizontal",worldLV,false,0);
  //G4VPhysicalVolume* horizontalPV  = new G4PVPlacement(0,G4ThreeVector(0*m, 0.00523*m, 0.00044*m),horizontalLV,"horizontal",worldLV,false,0);
  G4VPhysicalVolume* horizontalPV  = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),horizontalLV,"horizontal",worldLV,false,0);

  //Vertical volume (for the 7 vertical modules)
  G4Box* verticalHall_box =  new G4Box("VerticalHall",HorizonY, HorizonX, HorizonZ);
  G4LogicalVolume* verticalLV =new G4LogicalVolume(verticalHall_box, Air, "vertical_log",0,0,0);
  // untill 10a
  //G4VPhysicalVolume* verticalPV = new G4PVPlacement(0,G4ThreeVector(0.,18.0*cm,-4.0*m),verticalLV,"vertical",worldLV,false,0);
  //from jnubeam 10a
  G4VPhysicalVolume* verticalPV = new G4PVPlacement(0,G4ThreeVector(0.,0.,-4.0*m),verticalLV,"vertical",worldLV,false,0);


  //Proton Moduel volume added by kikawa
  G4Box* protonHall_box= new G4Box("ProtonHall",HorizonY, HorizonY, PHorizonZ);
  G4LogicalVolume* ProtonLV =new G4LogicalVolume(protonHall_box, Air, "proton_log",0,0,0);
  G4VPhysicalVolume* ProtonPV = new G4PVPlacement(0,G4ThreeVector(0.,0.,-1.2*m),ProtonLV,"proton",worldLV,false,0);




//////////////////////////////////////////////////////////////
//// We first create modules, than put materials inside /////
//////////////////////////////////////////////////////////////

  // MODULE ====================================
  //Modules objects:
  G4Box* module_box = new G4Box("Module",ModuleX, ModuleY, ModuleZ);
  //Proton Module object added by kikawa
  G4Box* Pmodule_box= new G4Box("PModule",ModuleX, ModuleY, PModuleZ);
  G4LogicalVolume* moduleLV[17];//changed 14->17 by kikawa
  
  for (int k=start_mod;k<stop_mod;k++){

#ifdef MOD3
    if(k!=3) continue;  // for only one module test
#endif
  //Positions of modules
  //G4double PosX=-6.3; G4double PosY=2.7;//Run 1                                                                                           
  //G4double PosX=-0.2; G4double PosY=1.4;//Run 2                                                                                           
  //G4double PosX=8.35; G4double PosY=3.6;//Run 3                                                                                          
  //G4double PosX=1.9; G4double PosY=2.7;                                                                                                   
  G4double PosX=0; G4double PosY=0;

  //Set 7 horizontal modules:
  if( k>=0 && k<7 ) {
    moduleLV[k] = new G4LogicalVolume(module_box, Air, "module_log_h");
    char moduleName[14];
    sprintf(moduleName,"module%d",k);

    G4double x = ModuleStart + ModuleSpace*k-PosX/100;
    new G4PVPlacement(0,G4ThreeVector(x,0.,0.),moduleLV[k],moduleName,horizontalLV,false,k);
  }


  //Set 7 vertical modules
  else if( k>=7 && k<14 ) {
    moduleLV[k] = new G4LogicalVolume(module_box, Air, "module_log_v");
    char moduleName[14];
    sprintf(moduleName,"module%d",k);
    G4double y = ModuleStart + ModuleSpace*(k-7)-PosY/100;    
    new G4PVPlacement(0,G4ThreeVector(0.,y,0.),moduleLV[k],moduleName,verticalLV,false,k);     
  }



  //Set proton module added by kikawa  
  else if(k==16){
  moduleLV[16] = new G4LogicalVolume(Pmodule_box, Air, "module_log_p");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),moduleLV[k],"module16",ProtonLV,false,k);
  }
  }



  // IRON BLOCK =========================================  
   
  G4Box* iron_block = new G4Box("Iron",iron_xy,iron_xy,iron_z);
  G4LogicalVolume* ironLV = new G4LogicalVolume(iron_block,Fe,"ironLV");
  
  G4VisAttributes* ironVisAtt = new G4VisAttributes(G4Color(0.7,0.,0.7)); // magenta
  //ironVisAtt->SetForceWireframe(true);
  ironVisAtt->SetForceSolid(true);
  ironLV->SetVisAttributes(ironVisAtt);

  //PLASTIC PLANES ADDED BY B.QUILAIN (2013/10/10)
  G4Box* Plastic_plane = new G4Box("ABS",Plastic_xy,Plastic_xy,Plastic_z);
  G4LogicalVolume* PlasticLV = new G4LogicalVolume(Plastic_plane,ABS,"PlasticLV");

  G4VisAttributes* PlasticVisAtt = new G4VisAttributes(G4Color(0.0,1.0,1.0)); // cyan 
                                                                        
  PlasticVisAtt->SetForceSolid(true);
  PlasticLV->SetVisAttributes(PlasticVisAtt);


  //SCINTILLATORS FOR TRACKING PLANES =======================================

  // vertical-layer scintillator (internal)
  std::vector<G4TwoVector> vdim;
  // horizontal-layer scintillator (internal)
  std::vector<G4TwoVector> hdim;

  // vertical-layer scintillator (internal)                                                   
  std::vector<G4TwoVector> vdim_Paint;
  // horizontal-layer scintillator (internal)                                                 
  std::vector<G4TwoVector> hdim_Paint;

  //SciBar type for Proton module added by kikawa
  // vertical-layer scintillator (internal)
  std::vector<G4TwoVector> svdim;
  // horizontal-layer scintillator (internal)
  std::vector<G4TwoVector> shdim;

  //INGRID type for Proton module added by kikawa
  // vertical-layer scintillator (internal)
  std::vector<G4TwoVector> ivdim;
  // horizontal-layer scintillator (internal)
  std::vector<G4TwoVector> ihdim;

  for(int iver=0; iver<nSciVertex; iver++){
    vdim.push_back( G4TwoVector( SciVertex_x[iver], SciVertex_y[iver] ) );
    hdim.push_back( G4TwoVector( SciVertex_y[nSciVertex-1-iver], 
    SciVertex_x[nSciVertex-1-iver] ) );
    vdim_Paint.push_back( G4TwoVector( SciPaint_x[iver], SciPaint_y[iver] ) );
    hdim_Paint.push_back( G4TwoVector( SciPaint_y[nSciVertex-1-iver], SciPaint_x[nSciVertex-1-iver] ) );
  }

 
  //SciBar Scintillator dimension
  SciVertex2_x[0] = -11.7*mm; 
  SciVertex2_x[1] = -12.25*mm; 
  SciVertex2_x[2] = -12.25*mm; 
  SciVertex2_x[3] = -11.7*mm; 
  SciVertex2_x[4] =  11.7*mm; 
  SciVertex2_x[5] =  12.25*mm; 
  SciVertex2_x[6] =  12.25*mm; 
  SciVertex2_x[7] =  11.7*mm;

  SciVertex2_y[0] = - 6.5 *mm; 
  SciVertex2_y[1] = - 3.5*mm; 
  SciVertex2_y[2] =   3.5*mm; 
  SciVertex2_y[3] =   6.5 *mm; 
  SciVertex2_y[4] =   6.5 *mm; 
  SciVertex2_y[5] =   3.5*mm; 
  SciVertex2_y[6] =  -3.5*mm; 
  SciVertex2_y[7] =  -6.5 *mm; 

  //for SciBar type
  for(int iver=0; iver<nSciVertex2; iver++){
    svdim.push_back( G4TwoVector( SciVertex2_x[iver], SciVertex2_y[iver] ) );
    shdim.push_back( G4TwoVector( SciVertex2_y[nSciVertex2-1-iver], 
    SciVertex2_x[nSciVertex2-1-iver] ) );
  }
 

  //INGRID Scintillator dimension for PM
  SciVertex3_x[0] = -23.65*mm; 
  SciVertex3_x[1] = -24.15*mm; 
  SciVertex3_x[2] = -24.15*mm; 
  SciVertex3_x[3] = -23.65*mm; 
  SciVertex3_x[4] =  23.65*mm; 
  SciVertex3_x[5] =  24.15*mm; 
  SciVertex3_x[6] =  24.15*mm; 
  SciVertex3_x[7] =  23.65*mm;

  SciVertex3_y[0] = - 5. *mm; 
  SciVertex3_y[1] = - 0.5*mm; 
  SciVertex3_y[2] =   0.5*mm; 
  SciVertex3_y[3] =   5. *mm; 
  SciVertex3_y[4] =   5. *mm; 
  SciVertex3_y[5] =   0.5*mm; 
  SciVertex3_y[6] =  -0.5*mm; 
  SciVertex3_y[7] =  -5. *mm; 


  //for INGRID type of PM
  for(int iver=0; iver<nSciVertex2; iver++){
    ivdim.push_back( G4TwoVector( SciVertex3_x[iver], SciVertex3_y[iver] ) );
    ihdim.push_back( G4TwoVector( SciVertex3_y[nSciVertex2-1-iver], 
    SciVertex3_x[nSciVertex2-1-iver] ) );
  }

  std::vector<G4ExtrudedSolid::ZSection> zsec;
  zsec.push_back( G4ExtrudedSolid::ZSection(-600*mm, G4TwoVector(0*mm,0*mm), 1) );
  zsec.push_back( G4ExtrudedSolid::ZSection(600*mm, G4TwoVector(0*mm,0*mm), 1) );
  
  //INGRID type for INGRID
  //INGRID vertical scintillator
  G4ExtrudedSolid* vscint_tmp = new G4ExtrudedSolid("vscint_tmp", vdim, zsec);
  G4EllipticalTube* vsci_hole = new G4EllipticalTube("vsci_hole",1.3*mm,1.95*mm,600*mm);
  G4SubtractionSolid* vscint_int = new G4SubtractionSolid("vscint_int",  vscint_tmp, vsci_hole);
  G4LogicalVolume* vscint_intLV = new G4LogicalVolume(vscint_int,Scinti,"vscint_intLV");
  G4VisAttributes* vscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));
  G4double maxStep = 0.00001;
  G4UserLimits * fStepLimit = new G4UserLimits(maxStep);
  vscint_intLV->SetUserLimits(fStepLimit);

  //vscint_intVisAtt->SetForceWireframe(true);
  vscint_intVisAtt->SetForceSolid(true);
  vscint_intLV->SetVisAttributes(vscint_intVisAtt);
  G4ExtrudedSolid * vscintPaint_tmp = new G4ExtrudedSolid("vscintPaint_tmp",vdim_Paint,zsec);
  G4SubtractionSolid* vscintPaint_int = new G4SubtractionSolid("vscintPaint_int",  vscintPaint_tmp, vscint_tmp);
  G4LogicalVolume* vscintPaint_intLV = new G4LogicalVolume(vscintPaint_int,Paint,"vscintPaint_intLV");
  G4VisAttributes* vscintPaint_intVisAtt = new G4VisAttributes(G4Color(0.,0.,0.));
  //vscintPaint_intVisAtt->SetForceWireframe(true);
  vscintPaint_intVisAtt->SetForceWireframe(true);
  vscintPaint_intLV->SetVisAttributes(vscintPaint_intVisAtt);
  
  //INGRID horizontal scintillator
  G4ExtrudedSolid* hscint_tmp = new G4ExtrudedSolid("hscint_tmp", hdim, zsec);
  G4EllipticalTube* hsci_hole = new G4EllipticalTube("hsci_hole",1.95*mm,1.3*mm,600*mm);
  G4SubtractionSolid* hscint_int = new G4SubtractionSolid("hscint_int",  hscint_tmp, hsci_hole);
  G4LogicalVolume* hscint_intLV = new G4LogicalVolume(hscint_int,Scinti,"hscint_intLV");
  G4VisAttributes* hscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  hscint_intLV->SetUserLimits(fStepLimit);

  //hscint_intVisAtt->SetForceWireframe(true);
  hscint_intVisAtt->SetForceSolid(true);
  hscint_intLV->SetVisAttributes(hscint_intVisAtt);
  G4ExtrudedSolid * hscintPaint_tmp = new G4ExtrudedSolid("hscintPaint_tmp",hdim_Paint,zsec);
  G4SubtractionSolid* hscintPaint_int = new G4SubtractionSolid("hscintPaint_int",  hscintPaint_tmp, hscint_tmp);
  G4LogicalVolume* hscintPaint_intLV = new G4LogicalVolume(hscintPaint_int,Paint,"hscintPaint_intLV");
  G4VisAttributes* hscintPaint_intVisAtt = new G4VisAttributes(G4Color(0.,0.,0.));
  //hscintPaint_intVisAtt->SetForceWireframe(true);
  hscintPaint_intVisAtt->SetForceWireframe(true);
  hscintPaint_intLV->SetVisAttributes(hscintPaint_intVisAtt);
  

  //SciBar type for proton module added by kikawa
  //SciBar vertical scintillator
  G4ExtrudedSolid* vscint2_tmp = new G4ExtrudedSolid("vscint2_tmp", svdim, zsec);
  G4EllipticalTube* vsci_hole2 = new G4EllipticalTube("vsci_hole2",0.9*mm,0.9*mm,600*mm);
  G4SubtractionSolid* vscint2_int = new G4SubtractionSolid("vscint2_int",  vscint2_tmp, vsci_hole2);
  G4LogicalVolume* vscint2_intLV = new G4LogicalVolume(vscint2_int,Scinti,"vscint2_intLV");
  G4VisAttributes* vscint2_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));

  vscint2_intLV->SetUserLimits(fStepLimit);

  //vscint2_intVisAtt->SetForceWireframe(true);
  vscint2_intVisAtt->SetForceSolid(true);
  vscint2_intLV->SetVisAttributes(vscint2_intVisAtt);

  //SciBar horizontal scintillator
  G4ExtrudedSolid* hscint2_tmp = new G4ExtrudedSolid("hscint2_tmp", shdim, zsec);
  G4EllipticalTube* hsci_hole2 = new G4EllipticalTube("hsci_hole2",0.9*mm,0.9*mm,600*mm);
  G4SubtractionSolid* hscint2_int = new G4SubtractionSolid("hscint2_int",  hscint2_tmp, hsci_hole2);
  G4LogicalVolume* hscint2_intLV = new G4LogicalVolume(hscint2_int,Scinti,"hscint2_intLV");
  G4VisAttributes* hscint2_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  hscint2_intLV->SetUserLimits(fStepLimit);

  //hscint2_intVisAtt->SetForceWireframe(true);
  hscint2_intVisAtt->SetForceSolid(true);
  hscint2_intLV->SetVisAttributes(hscint2_intVisAtt);




  //INGRID type for proton module added by kikawa
  //INGRID vertical scintillator for PM
  G4ExtrudedSolid* vscint3_tmp = new G4ExtrudedSolid("vscint3_tmp", ivdim, zsec);
  G4EllipticalTube* vsci_hole3 = new G4EllipticalTube("vsci_hole3",1.3*mm,1.95*mm,600*mm);
  G4SubtractionSolid* vscint3_int = new G4SubtractionSolid("vscint3_int",  vscint3_tmp, vsci_hole3);
  G4LogicalVolume* vscint3_intLV = new G4LogicalVolume(vscint3_int,Scinti,"vscint3_intLV");
  G4VisAttributes* vscint3_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));
  //vscint3_intVisAtt->SetForceWireframe(true);
  vscint3_intVisAtt->SetForceSolid(true);
  vscint3_intLV->SetVisAttributes(vscint3_intVisAtt);

  //INGRID horizontal scintillator for PM
  G4ExtrudedSolid* hscint3_tmp = new G4ExtrudedSolid("hscint3_tmp", ihdim, zsec);
  G4EllipticalTube* hsci_hole3 = new G4EllipticalTube("hsci_hole3",1.95*mm,1.3*mm,600*mm);
  G4SubtractionSolid* hscint3_int = new G4SubtractionSolid("hscint3_int",  hscint3_tmp, hsci_hole3);
  G4LogicalVolume* hscint3_intLV = new G4LogicalVolume(hscint3_int,Scinti,"hscint3_intLV");
  G4VisAttributes* hscint3_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hscint3_intVisAtt->SetForceWireframe(true);
  hscint3_intVisAtt->SetForceSolid(true);
  hscint3_intLV->SetVisAttributes(hscint3_intVisAtt);



  // horizontal-layer scintillator
  /*
  G4Box* hlayer = new G4Box("hlayer",scibar_x,scibar_y,scibar_z);
  G4LogicalVolume* hlayerLV = new G4LogicalVolume(hlayer,Scinti,"hlayerLV");
  G4VisAttributes* hlayerVisAtt = new G4VisAttributes(G4Color(0.,0.7,0.));
  hlayerVisAtt->SetForceWireframe(true);
  hlayerLV->SetVisAttributes(hlayerVisAtt);

    //vertical-layer scintillator
  G4Box* vlayer = new G4Box("vlayer",scibar_y,scibar_x,scibar_z);
  G4LogicalVolume* vlayerLV = new G4LogicalVolume(vlayer,Scinti,"vlayerLV");
  G4VisAttributes* vlayerVisAtt = new G4VisAttributes(G4Color(0.,0.7,0.));
  vlayerVisAtt->SetForceWireframe(true);
  vlayerLV->SetVisAttributes(vlayerVisAtt);
  */



  //SciBar type for proton module added by kikawa
  //SciBar vertical scintillator
  /*
  G4Box* hlayer2 = new G4Box("hlayer2",scibar_x,scibar2_y,scibar2_z);
  G4LogicalVolume* hlayer2LV = new G4LogicalVolume(hlayer2,Scinti,"hlayer2LV");
  G4VisAttributes* hlayer2VisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
  hlayer2VisAtt->SetForceWireframe(true);
  hlayer2LV->SetVisAttributes(hlayer2VisAtt);

  //SciBar horizontal scintillator
  G4Box* vlayer2 = new G4Box("vlayer2",scibar2_y,scibar_x,scibar2_z);
  G4LogicalVolume* vlayer2LV = new G4LogicalVolume(vlayer2,Scinti,"vlayer2LV");
  G4VisAttributes* vlayer2VisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
  vlayer2VisAtt->SetForceWireframe(true);
  vlayer2LV->SetVisAttributes(vlayer2VisAtt);
  */



  
  //SCINTILLATORS FOR VETO PLANES ==========================================

  G4VisAttributes* vetoVisAtt = new G4VisAttributes(G4Color(0.,0.6,0.5));
  vetoVisAtt->SetForceWireframe(true);

  //Long Veto plane
  G4Box* Lveto_box = new G4Box("Lveto_box",Lveto_x,Lveto_y,Lveto_z);
  G4LogicalVolume *LvetoLV = new G4LogicalVolume(Lveto_box,Scinti,"LvetoLV");

  LvetoLV->SetVisAttributes(G4VisAttributes::Invisible);
  //LvetoLV->SetVisAttributes(vetoVisAtt);
  
  //Short Veto plane
  G4Box* Sveto_box = new G4Box("Sveto_box",Sveto_x,Sveto_y,Sveto_z);
  G4LogicalVolume *SvetoLV = new G4LogicalVolume(Sveto_box,Scinti,"SvetoLV");

  SvetoLV->SetVisAttributes(G4VisAttributes::Invisible);
  //SvetoLV->SetVisAttributes(vetoVisAtt);



  //Proton module Long Veto plane added by kikawa
  G4Box* PLveto_box = new G4Box("PLveto_box",PLveto_x,PLveto_y,PLveto_z);
  G4LogicalVolume *PLvetoLV = new G4LogicalVolume(PLveto_box,Scinti,"PLvetoLV");

  PLvetoLV->SetVisAttributes(G4VisAttributes::Invisible);
  //PLvetoLV->SetVisAttributes(vetoVisAtt);
  
  //Proton module Short Veto plane added by kikawa
  G4Box* PSveto_box = new G4Box("PSveto_box",PSveto_x,PSveto_y,PSveto_z);
  G4LogicalVolume *PSvetoLV = new G4LogicalVolume(PSveto_box,Scinti,"PSvetoLV");

  PSvetoLV->SetVisAttributes(G4VisAttributes::Invisible);
  //PSvetoLV->SetVisAttributes(vetoVisAtt);




  //Rotation matrix for veto planes
  G4RotationMatrix *TB2LR = new G4RotationMatrix(G4ThreeVector(0,0,1.),90.*degree);


  //POSITIONNING OF ALL THE ELEMENTS ================================================

  // 7 horizontal Modules
  for(int k=start_mod;k<stop_mod;k++){
    

    //if(k!=16){
    if(k<14){
    //Module HO3 has gap 0.2cm shorter than other modules
    G4double Gap2;
    G4double iron_start2,scibar_start2;

    if(k==3){
      Gap2=Gap-0.2*cm; 
      iron_start2=iron_start+0.9*cm;
      scibar_start2=scibar_start+1*cm;
    }
    else{
      Gap2=Gap;
      iron_start2=iron_start;
      scibar_start2=scibar_start;
    }
    
#ifdef MOD3
    if(k!=3) continue; // for only one module test
#endif

    // 9 iron-blocs per module -----------------------------------------------------
    for(int i=0;i<9;i++){
      char ironname[30];
      sprintf(ironname,"iron[%d][%d]",k,i);
      G4double z = iron_start2 + Gap2*i;
      new G4PVPlacement(0,G4ThreeVector(0.,0.,z),ironLV,ironname,moduleLV[k],false,i+k*9 ); // the world in Module! Hereiron platyes are added     
    }
    for(int i=0;i<11;i++){
      char Plasticname_U[30];
      char Plasticname_D[30];
      sprintf(Plasticname_U,"Plastic_U[%d][%d]",k,i);
      sprintf(Plasticname_D,"Plastic_D[%d][%d]",k,i);
      G4double z1 = scibar_start2 + Gap2*i - Plastic_z;
      G4double z2 = scibar_start2 + Gap2*i + scibar_z + 2*scibar_z + scibar_z + Plastic_z;
      //new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),PlasticLV,Plasticname_U,moduleLV[k],false,i*2+k*22 );
      //new G4PVPlacement(0,G4ThreeVector(0.,0.,z2),PlasticLV,Plasticname_D,moduleLV[k],false,1+i*2+k*22 );
    }

    // 11 planes of scintillator per module-----------------------------------------
    for(int i=0;i<11;i++){
      G4double z = scibar_start2 + scibar_z + Gap2*i;
      for(int j=0;j<24;j++){ 
	char name[30]; 
	char name_Paint[30];

	sprintf(name,"vlayer[%d][%d][%d]",k,i,j);
	sprintf(name_Paint,"vlayer_Paint[%d][%d][%d]",k,i,j);

	G4double x_scinti =  scibar_xy_start + 2*j*scibar_y; 
	//new G4PVPlacement(0,G4ThreeVector(x_scinti,0.,z),vlayerLV,name,moduleLV[k],false,j+i*24+k*264); // in Module        
	new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z),vscint_intLV,name,moduleLV[k],false,j+i*24+k*264); // in Module        
	//new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z),vscintPaint_intLV,name_Paint,moduleLV[k],false,j+i*24+k*264);
      
	sprintf(name,"hlayer[%d][%d][%d]",k,i,j);
	sprintf(name_Paint,"vlayer_Paint[%d][%d][%d]",k,i,j);

	G4double y_scinti =  scibar_xy_start + 2*j*scibar_y;     
	//new G4PVPlacement(0,G4ThreeVector(0.,y_scinti,z+2*scibar_z),hlayerLV,name,moduleLV[k],false,j+i*24+k*264); // in Module
	new G4PVPlacement(yrot,G4ThreeVector(0.,y_scinti,z+2*scibar_z),hscint_intLV,name,moduleLV[k],false,j+i*24+k*264); // in Module
	//new G4PVPlacement(yrot,G4ThreeVector(0.,y_scinti,z+2*scibar_z),hscintPaint_intLV,name_Paint,moduleLV[k],false,j+i*24+k*264); // in Module

      }
    }

    //4 veto planes in the first horizontal modules-----------------------------------------------
    if(k==0){
        // 22 veto-bars per veto-plane (4 veto-planes for this module)
      for(int i=0;i<22;i++){
		G4double z = veto_start + 2*i*Sveto_z;

		char vetoname[4][22];

		sprintf(vetoname[0],"veto[%d][0][%d]",k,i);
		sprintf(vetoname[1],"veto[%d][1][%d]",k,i);
		sprintf(vetoname[2],"veto[%d][2][%d]",k,i);
		sprintf(vetoname[3],"veto[%d][3][%d]",k,i);

		// Precise positions in module from Oscar
		new G4PVPlacement(0,G4ThreeVector(5.9*cm,68.4*cm,z),LvetoLV,vetoname[3],moduleLV[k],false, i+66);  // UP
		new G4PVPlacement(0,G4ThreeVector(-0.9*cm,-65.9*cm,z),SvetoLV,vetoname[2],moduleLV[k],false, i+44); // DOWN
		new G4PVPlacement(TB2LR,G4ThreeVector(70.9*cm,0.3*cm,z),LvetoLV,vetoname[1],moduleLV[k],false, i+22); // LEFT
		new G4PVPlacement(TB2LR,G4ThreeVector(-70.575*cm,-3.7*cm,z),LvetoLV,vetoname[0],moduleLV[k],false, i); // RIGTH
	
		//new G4PVPlacement(TB2LR,G4ThreeVector(75.0*cm,0.3*cm,z),LvetoLV,vetoname[1],moduleLV[k],false, i+22); // RIGTH
		//new G4PVPlacement(TB2LR,G4ThreeVector(80.0*cm,0.3*cm,z),LvetoLV,vetoname[1],moduleLV[k],false, i+22); // RIGTH


      }  
    }
    
    //3 veto planes on the 6 horizontal modules-------------------------------------------
    if (k>0 && k<7){
      // 22 veto-bars per veto-plane (3 veto-plane / module)
      for(int i=0;i<22;i++){
		G4double z = veto_start + 2*i*Sveto_z;
		
		char vetoname[3][22];

		sprintf(vetoname[0],"veto[%d][1][%d]",k,i);
		sprintf(vetoname[1],"veto[%d][2][%d]",k,i);
		sprintf(vetoname[2],"veto[%d][3][%d]",k,i);

		// Precise positions in module from Oscar
		new G4PVPlacement(0,G4ThreeVector(5.9*cm,68.4*cm,z),LvetoLV,vetoname[2],moduleLV[k],false,88*k + i+66); // UP
		new G4PVPlacement(0,G4ThreeVector(-0.9*cm,-65.9*cm,z),SvetoLV,vetoname[1],moduleLV[k],false,88*k + i+44); // DOWN
		new G4PVPlacement(TB2LR,G4ThreeVector(70.9*cm,0.3*cm,z),LvetoLV,vetoname[0],moduleLV[k],false,88*k + i+22); // LEFT

      } 
    }

    //4 veto planes on the first vertical module ----------------------------------------------------
    if (k==7){
      // 22 veto-bars per veto-plane (4 veto-planes for this module)
      for(int i=0;i<22;i++){
		G4double z = veto_start + 2*i*Sveto_z;
	
		char vetoname[4][22];

		sprintf(vetoname[0],"veto[%d][0][%d]",k,i);
		sprintf(vetoname[1],"veto[%d][1][%d]",k,i);
		sprintf(vetoname[2],"veto[%d][2][%d]",k,i);
		sprintf(vetoname[3],"veto[%d][3][%d]",k,i);

		// Precise positions in module from Oscar
		new G4PVPlacement(0,G4ThreeVector(5.9*cm,68.4*cm,z),LvetoLV,vetoname[3],moduleLV[k],false, 88*k + i+66);  // UP
		new G4PVPlacement(0,G4ThreeVector(-0.9*cm,-65.9*cm,z),SvetoLV,vetoname[2],moduleLV[k],false, 88*k + i+44); // DOWN
		new G4PVPlacement(TB2LR,G4ThreeVector(70.9*cm,0.3*cm,z),LvetoLV,vetoname[1],moduleLV[k],false, 88*k + i+22); // LEFT
		new G4PVPlacement(TB2LR,G4ThreeVector(-70.575*cm,-3.7*cm,z),LvetoLV,vetoname[0],moduleLV[k],false, 88*k + i); // RIGTH

      }  

    }

    //3 veto planes on the remaining vertical modules
    if(k>7){
      // 22 veto-bars per veto-plane (4 veto-planes for this module)
      for(int i=0;i<22;i++){
		G4double z = veto_start + 2*i*Sveto_z;

		char vetoname[3][22];
	
		sprintf(vetoname[0],"veto[%d][1][%d]",k,i);
		sprintf(vetoname[1],"veto[%d][2][%d]",k,i);
		sprintf(vetoname[2],"veto[%d][3][%d]",k,i);

	
		// Precise positions in module from Oscar
		new G4PVPlacement(0,G4ThreeVector(5.9*cm,68.4*cm,z),LvetoLV,vetoname[2],moduleLV[k],false, 88*k + i+66);  // UP
		new G4PVPlacement(TB2LR,G4ThreeVector(70.9*cm,0.3*cm,z),LvetoLV,vetoname[1],moduleLV[k],false, 88*k + i+22); // LEFT
		new G4PVPlacement(TB2LR,G4ThreeVector(-70.575*cm,-3.7*cm,z),LvetoLV,vetoname[0],moduleLV[k],false,88*k + i); // RIGTH

      }
    }
    }

  //Proton Module added by kikawa
  
    else if(k==16){
  for(int i=0;i<18;i++){
    G4double z;
    if (i==0){
      z = Pscibar_start + scibar_z;
      G4cout << i << " "<<z<<G4endl;
    }
    else{
      z = Pscibar_start + scibar_z + dist_first + dist + 2*dist*(i-1);
      G4cout << i << " "<<z<<G4endl;
    }

    //First plane similar to a INGRID module TKP, will be a veto plane
    if (i==0){
      for(int j=0;j<24;j++){ 
        char name[10]; 
	
        sprintf(name,"x%d",j+k*264);
        G4double x_scinti =  scibar_xy_start + 2*j*scibar_y;
	new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z+dist),vscint3_intLV,name,moduleLV[16],false,j+k*264); // in Module
	
        sprintf(name,"y%d",j+k*264);
        G4double y_scinti =   scibar_xy_start + 2*j*scibar_y;
	new G4PVPlacement(yrot,G4ThreeVector(0.,y_scinti,z),hscint3_intLV,name,moduleLV[16],false,j+k*264); // in Module     
	
	}
    }
    else{
      //---------- 17 TKP

      // 8 INGRID type scintillators
      for(int j=0;j<8;j++){ 
        char name[10];
	
	sprintf(name,"vlayer[%d][%d][%d]",k,i,j);
        G4double x_scinti = scibar_xy_start + 2*j*scibar_y; 
	new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z+dist),vscint3_intLV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // vertical bars
	
	sprintf(name,"hlayer[%d][%d][%d]",k,i,j);
        G4double y_scinti = scibar_xy_start + 2*j*scibar_y;       
	new G4PVPlacement(yrot,G4ThreeVector(0,y_scinti,z),hscint3_intLV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // horizontal bars      

	
      }

      //16 SciBar type scintillators

      for(int j=8;j<24;j++){ 
        char name[10]; 
	
	sprintf(name,"vlayer[%d][%d][%d]",k,i,j);
        G4double x_scinti = scibar_xy_start + 15*scibar_y +(2*(j-8)+1)*scibar2_y; 
        //new G4PVPlacement(0,G4ThreeVector(x_scinti,0.,z+dist),vlayer2LV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // vertical bars
	new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z+dist),vscint2_intLV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // vertical bars
	
	sprintf(name,"hlayer[%d][%d][%d]",k,i,j);
        G4double y_scinti = scibar_xy_start + 15*scibar_y +(2*(j-8)+1)*scibar2_y; 
	//new G4PVPlacement(0,G4ThreeVector(0,y_scinti,z),hlayer2LV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // horizontal bars
	new G4PVPlacement(yrot,G4ThreeVector(0,y_scinti,z),hscint2_intLV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // horizontal bars
	
      }

      // 8 INGRID type scintillators 
      for(int j=24;j<32;j++){ 
        char name[10]; 
	
	sprintf(name,"vlayer[%d][%d][%d]",k,i,j);
        G4double x_scinti = scibar_xy_start + 15*scibar_y +2*16*scibar2_y+(2*(j-24)+1)*scibar_y; 
	new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z+dist),vscint3_intLV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // vertical bars

	sprintf(name,"hlayer[%d][%d][%d]",k,i,j);
        G4double y_scinti = scibar_xy_start + 15*scibar_y +2*16*scibar2_y+(2*(j-24)+1)*scibar_y; 
        new G4PVPlacement(yrot,G4ThreeVector(0,y_scinti,z),hscint3_intLV,name,moduleLV[16],false,j+32*(i-1)+k*264+24); // horizontal bars


      }
    }

  }

    //Veto Planes  
      for(int i=0;i<17;i++){
	G4double z = Pveto_start + Sveto_z + 2*i*Sveto_z;
	char vetoname[4][22];

	sprintf(vetoname[0],"veto[%d][0][%d]",k,i);
	sprintf(vetoname[1],"veto[%d][1][%d]",k,i);
	sprintf(vetoname[2],"veto[%d][2][%d]",k,i);
	sprintf(vetoname[3],"veto[%d][3][%d]",k,i);

	// Precise positions in module from Oscar

	new G4PVPlacement(0,G4ThreeVector(-0.5*cm,65.5*cm,z),PLvetoLV,vetoname[3],moduleLV[k],false, 88*k+i+17*2);  // UP
	new G4PVPlacement(0,G4ThreeVector(1.5*cm,-65.5*cm,z),PSvetoLV,vetoname[2],moduleLV[k],false, 88*k+i); // DOWN
	new G4PVPlacement(TB2LR,G4ThreeVector(65.5*cm,-2.5*cm,z),PLvetoLV,vetoname[1],moduleLV[k],false, 88*k+i+17); // LEFT
	new G4PVPlacement(TB2LR,G4ThreeVector(-65.5*cm,0.5*cm,z),PLvetoLV,vetoname[0],moduleLV[k],false, 88*k+i+17*3); // RIGTH
	      }
   
  }
  
  }







  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  G4String hlayerSDname = "Ingrid/hlayerSD";
  ahlayerSD = new IngridHLayerSD( hlayerSDname , CBIRKS);
  SDman->AddNewDetector( ahlayerSD );
  //hlayerLV->SetSensitiveDetector( ahlayerSD );
  //hlayer2LV->SetSensitiveDetector( ahlayerSD );
  hscint_intLV->SetSensitiveDetector( ahlayerSD );
  hscint2_intLV->SetSensitiveDetector( ahlayerSD );//SciBar type added by kikawa
  hscint3_intLV->SetSensitiveDetector( ahlayerSD );//INGRID type added by kikawa
 
  G4String vlayerSDname = "Ingrid/vlayerSD";
  avlayerSD = new IngridVLayerSD( vlayerSDname, CBIRKS );
  SDman->AddNewDetector( avlayerSD );
  //vlayerLV->SetSensitiveDetector( avlayerSD );
  //vlayer2LV->SetSensitiveDetector( avlayerSD );
  vscint_intLV->SetSensitiveDetector( avlayerSD );
  vscint2_intLV->SetSensitiveDetector( avlayerSD );//SciBar type added by kikawa
  vscint3_intLV->SetSensitiveDetector( avlayerSD );//INGRID type added by kikawa
 
  G4String vetoSDname = "Ingrid/vetoSD";
  avetoSD = new IngridVetoSD( vetoSDname , CBIRKS);
  SDman->AddNewDetector( avetoSD );

  SvetoLV->SetSensitiveDetector( avetoSD );
  LvetoLV->SetSensitiveDetector( avetoSD );
  PSvetoLV->SetSensitiveDetector( avetoSD );//Veto for Proton Module added by kikawa
  PLvetoLV->SetSensitiveDetector( avetoSD );//Veto for Proton Module added by kikawa

  return worldPV;

}

//___________________________________________________________________________________________________________

void IngridDetectorConstruction::DefineSpace()
{
  // The length at World volume, volume of 7 modules, and module volume means the half length.
  //WorldX=5.50*m;
  //WorldY=5.50*m;
  //WorldZ=5.50*m;
  WorldX=20.0*m;
  WorldY=10.0*m;
  WorldZ=20.0*m;

  //The dimension of Hall dirt volume
  HallDirtRadiusMin =  8.5*m;
  HallDirtRadiusMax = (10. + 3.2) *m;
  HallDirtHeight    = 10. *m;
  HallDirtSPhi      = 0;
  HallDirtDPhi      = 3.1415926535897932384626*2;
  HallDirtX         = -2.167 *m;
  HallDirtY         = 0. *m;
  //HallDirtZ         = -11.55*m;
  HallDirtZ         = 1.7*m;
  HallX             = HallDirtX;
  HallZ             = HallDirtZ;

  //Volume of 7 modules
  HorizonX=5.25*m;
  HorizonY=1.0*m;
  //HorizonZ=90.0*cm;
  HorizonZ=68.0*cm;


  //Volume of Proton module added by kikawa
  PHorizonZ=45.0*cm;

  //Module volume
  //ModuleX=0.82*m;  // overlap of logical volumes is not good !!!
  //ModuleY=0.82*m;
  ModuleX=0.75*m;
  ModuleY=0.75*m;
  //ModuleZ=70.0*cm;
  ModuleZ=65.0*cm;

  //Proton Module volume
  PModuleZ=44.0*cm;
  

  ModuleSpace= 1.50*m;     // 120 + 30  distance between the centers of two consecutive modules
  ModuleStart=-4.50*m;     // -150*3 center of the first module
}

//___________________________________________________________________________________________________________

void IngridDetectorConstruction::DefineStructures()
{
  //Iron blocks
  iron_z=3.25*cm;
  iron_xy= 62.0*cm;
  Gap= 10.7*cm;           // 6.5 + 2 +2*1.1 means " scinti-scinti-gap-iron-gap-scinti-scinti-... "
  iron_start= -48.15*cm;  // Position of the center of the first iron block  
	    		// = -54.5cm(1st plane surface) + 2cm(width of 1 plane) + 1.1cm(gap of scinti-iron) +3.25cm(half of iron)
  

  //Plastic Planes
  Plastic_z=.1*cm;
  Plastic_xy=60.25*cm;

  //Same as before, but doubles instead of G4doubles for use by other functions
  Niron_start=-48.15;
  Niron_z=3.25;
  Nscibar_start=-54.5;
  NGap=10.7;
  NModuleSpace=150;

  //Scintillators
  scibar_x= 0.6*m;        // 120cm long
  scibar_y=2.5*cm;        // 5cm wide
  scibar_z= 5.0*mm;       // 1cm Thick

  //Scibar type for proton module added by kikawa
  scibar2_y=1.25*cm;        // 2.5cm wide
  scibar2_z= 6.5*mm;       // 1.3cm thick


  //Proton Module
  Pscibar_start=-40.95*cm;  
  Pveto_start=-40.95*cm-1.55*cm;  


  scibar_start=-54.5*cm;  // in Module , which is surface of 1st tracking x-layer
  scibar_xy_start= -57.5*cm; // in Module   57.5 = 12mai*5cm - 5cm/2

  //Long veto planes
  Lveto_x= 0.65*m;        // 130cm long
  Lveto_y= 0.5*cm;        // 1cm thick
  Lveto_z= 2.5*cm;        // 5cm wide 
  veto_start=-52.5*cm;

  //Short veto planes
  Sveto_x=0.56*m;         // 112cm long
  Sveto_y=0.5*cm;         // 1cm thick
  Sveto_z=2.5*cm;         // 5cm wide 

  //Protn Module Long veto planes added by kikawa
  PLveto_x= 0.625*m;        // 125cm long
  PLveto_y= 0.5*cm;        // 1cm thick
  PLveto_z= 2.5*cm;        // 5cm wide 

  //Proton Module Short veto planes added by kikawa
  PSveto_x=0.6*m;         // 120cm long
  PSveto_y=0.5*cm;         // 1cm thick
  PSveto_z=2.5*cm;         // 5cm wide 


  //Distance between planes of Proton Module added by kikawa
  dist_first=2.7*cm;
  dist=2.3*cm;


  //Scintillator dimension (changed)
  /*
  SciVertex_x[0] = -22.7*mm; 
  SciVertex_x[1] = -24.2*mm; 
  SciVertex_x[2] = -24.2*mm; 
  SciVertex_x[3] = -22.7*mm; 
  SciVertex_x[4] =  22.7*mm; 
  SciVertex_x[5] =  24.2*mm; 
  SciVertex_x[6] =  24.2*mm; 
  SciVertex_x[7] =  22.7*mm;

  SciVertex_y[0] = - 5. *mm; 
  SciVertex_y[1] = - 0.5*mm; 
  SciVertex_y[2] =   0.5*mm; 
  SciVertex_y[3] =   5. *mm; 
  SciVertex_y[4] =   5. *mm; 
  SciVertex_y[5] =   0.5*mm; 
  SciVertex_y[6] =  -0.5*mm; 
  SciVertex_y[7] =  -5. *mm; 
  */

  //The good simulation of the scintillator
  
  SciVertex_x[0] = -21.5*mm;                                                                                        
  SciVertex_x[1] = -23.4*mm;                                                                                                     
  SciVertex_x[2] = -24.5*mm;                                                                                                 
  SciVertex_x[3] = -23.4*mm;
  SciVertex_x[4] = -23.4*mm;
  SciVertex_x[5] = -21.5*mm;                                                                                                          
  SciVertex_x[6] =  21.0*mm;                                                                                                       
  SciVertex_x[7] =  23.9*mm;                                                                                                     
  SciVertex_x[8] =  24.0*mm;                                                                                                  
  SciVertex_x[9] =  24.0*mm;
  SciVertex_x[10] =  23.9*mm;
  SciVertex_x[11] =  21.0*mm;                                                                                                          

  SciVertex_y[0] = -5.0 *mm;                                                                                                           
  SciVertex_y[1] = - 4.0*mm;                                                                                                        
  SciVertex_y[2] = -1.0*mm;                                                                                                            
  SciVertex_y[3] = 1.0*mm;                                                                                                       
  SciVertex_y[4] = 4.0 *mm;                                                                                                     
  SciVertex_y[5] =   5.0*mm;
  SciVertex_y[6] = 5.0 *mm;                                                                                                           
  SciVertex_y[7] = 4.0*mm;                                                                                                           
  SciVertex_y[8] = 1.0*mm;                                                                                                             
  SciVertex_y[9] = -1.0*mm;                                                                                                         
  SciVertex_y[10] = -4.0 *mm; 
  SciVertex_y[11] =   -5.0*mm;
  
  /*
  SciVertex_x[0] = -24.7*mm;                                                                                        
  SciVertex_x[1] = -24.8*mm;                                                                                                     
  SciVertex_x[2] = -24.9*mm;                                                                                                 
  SciVertex_x[3] = -24.9*mm;
  SciVertex_x[4] = -24.8*mm;
  SciVertex_x[5] = -24.7*mm;

  SciVertex_x[6] = 24.7*mm;                                                                                        
  SciVertex_x[7] = 24.8*mm;                                                                                                     
  SciVertex_x[8] = 24.9*mm;                                                                                                 
  SciVertex_x[9] = 24.9*mm;
  SciVertex_x[10] = 24.8*mm;
  SciVertex_x[11] = 24.7*mm;
                                                                                                          

  SciVertex_y[0] = -4.5 *mm;                                                                                                           
  SciVertex_y[1] = - 4.0*mm;                                                                                                        
  SciVertex_y[2] = -1.0*mm;                                                                                                            
  SciVertex_y[3] = 1.0*mm;                                                                                                       
  SciVertex_y[4] = 4.0 *mm;                                                                                                     
  SciVertex_y[5] =   4.5*mm;
  SciVertex_y[6] = 4.5 *mm;                                                                                                           
  SciVertex_y[7] = 4.0*mm;                                                                                                           
  SciVertex_y[8] = 1.0*mm;                                                                                                             
  SciVertex_y[9] = -1.0*mm;                                                                                                         
  SciVertex_y[10] = -4.0 *mm; 
  SciVertex_y[11] =   -4.5*mm;
*/
  /*
  SciVertex_x[0] = -24.7*mm;                                                                                        
  SciVertex_x[1] = -24.8*mm;                                                                                                     
  SciVertex_x[2] = -24.9*mm;                                                                                                 
  SciVertex_x[3] = -24.9*mm;
  SciVertex_x[4] = -24.8*mm;
  SciVertex_x[5] = -24.7*mm;

  SciVertex_x[6] = 24.7*mm;                                                                                        
  SciVertex_x[7] = 24.8*mm;                                                                                                     
  SciVertex_x[8] = 24.9*mm;                                                                                                 
  SciVertex_x[9] = 24.9*mm;
  SciVertex_x[10] = 24.8*mm;
  SciVertex_x[11] = 24.7*mm;
                                                                                                          

  SciVertex_y[0] = -5.0 *mm;                                                                                                           
  SciVertex_y[1] = - 4.0*mm;                                                                                                        
  SciVertex_y[2] = -1.0*mm;                                                                                                            
  SciVertex_y[3] = 1.0*mm;                                                                                                       
  SciVertex_y[4] = 4.0 *mm;                                                                                                     
  SciVertex_y[5] =   5.0*mm;
  SciVertex_y[6] = 5.0 *mm;                                                                                                           
  SciVertex_y[7] = 4.0*mm;                                                                                                           
  SciVertex_y[8] = 1.0*mm;                                                                                                             
  SciVertex_y[9] = -1.0*mm;                                                                                                         
  SciVertex_y[10] = -4.0 *mm; 
  SciVertex_y[11] =   -5.0*mm;
*/
  
  for(int i=0;i<nSciVertex;i++){
    if(SciVertex_x[i]<0) SciPaint_x[i]=SciVertex_x[i]-0*mm;
    else SciPaint_x[i]=SciVertex_x[i]+0*mm;
    if(SciVertex_y[i]<0) SciPaint_y[i]=SciVertex_y[i]-0*mm;
    else SciPaint_y[i]=SciVertex_y[i]+0*mm;
    //G4cout<<"i="<<i<<", x position="<<SciPaint_x[i]<<", y position="<<SciPaint_y[i]<<G4endl;
  }

  //**original
  /*  SciVertex_x[0] = -23.8*mm; 
  SciVertex_x[1] = -24.3*mm; 
  SciVertex_x[2] = -24.3*mm; 
  SciVertex_x[3] = -23.8*mm; 
  SciVertex_x[4] =  23.8*mm; 
  SciVertex_x[5] =  24.3*mm; 
  SciVertex_x[6] =  24.3*mm; 
  SciVertex_x[7] =  23.8*mm;

  SciVertex_y[0] = - 5. *mm; 
  SciVertex_y[1] = - 0.5*mm; 
  SciVertex_y[2] =   0.5*mm; 
  SciVertex_y[3] =   5. *mm; 
  SciVertex_y[4] =   5. *mm; 
  SciVertex_y[5] =   0.5*mm; 
  SciVertex_y[6] =  -0.5*mm; 
  SciVertex_y[7] =  -5. *mm; 
  */
 

}

//___________________________________________________________________________________________________________

void IngridDetectorConstruction::DefineMaterial()
{
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4String name, symbol;
  G4int nel;
  G4double iz;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
  a = 28.1*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);
  a = 40.1*g/mole;
  G4Element* elCa = new G4Element(name="Calusium", symbol="Ca", z=20., a);
  a = 23.0*g/mole;
  G4Element* elNa = new G4Element(name="Sodium", symbol="Na", z=11., a);
  a = 55.8*g/mole;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z=26., a);
  a = 27.0*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium", symbol="Al", z=13., a);
  a = 47.9*g/mole;
  G4Element* elTi = new G4Element(name="Titanium", symbol="Ti", z=22., a);

  //Air
  density = 1.29*mg/cm3;
  Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  //Iron
  a = 55.845*g/mole;
  density = 7.86*g/cm3;
  Fe = new G4Material(name="Iron", z=26., a, density);

  //ABS == Plastic (acrylonitrile butdadiene styrene)
  density= 1.05*g/cm3;
  ABS = new G4Material(name="ABS",density,nel=3);
  ABS->AddElement(elC,15);
  ABS->AddElement(elH,17);
  ABS->AddElement(elN,1);

  //Painting (TiO2) Added by B.Quilain
  density=4.250*g/cm3;
  Paint = new G4Material(name="Paint",density,nel=2);
  Paint->AddElement(elO,2);
  Paint->AddElement(elTi,1);

  //Water
  density = 1.000*g/cm3;
  Water = new G4Material(name="Water",density,nel=2);
  Water->AddElement(elH,2);
  Water->AddElement(elO,1);

  //Scintillator
  density = 1.032*g/cm3;
  Scinti = new G4Material(name="Scintillator", density, nel=2);
  Scinti->AddElement(elC, 9);
  Scinti->AddElement(elH, 10);
  Scinti->GetIonisation()->SetBirksConstant(0.208*mm/MeV);

  //Concrete
  density = 2.2*g/cm3;
  Concrete = new G4Material(name="Concrete", density, nel=6);

  Concrete->AddElement(elO, .53);
  Concrete->AddElement(elSi, .335);
  Concrete->AddElement(elCa, 0.06);
  Concrete->AddElement(elNa, 0.015);
  Concrete->AddElement(elFe, 0.02);
  Concrete->AddElement(elAl, 0.04);

}

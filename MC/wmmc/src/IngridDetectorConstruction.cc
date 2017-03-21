#include <stdio.h>
#include <iostream>
#include <vector>

#include "IngridDetectorConstruction.hh"
#include "IngridHLayerSD.hh"
#include "IngridVLayerSD.hh"
#include "IngridVetoSD.hh"
#include "INGRID_Dimension.hh"

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

#include "globals.hh"

#include "Const.hh"

#define PROTON 2
#define HORIZONTAL 3
#define VERTICAL 4
#define WAGASCI 5
#define WAGASCIBG 6

//#define MOD3

IngridDetectorConstruction::IngridDetectorConstruction(int MODE)
{
    mode = MODE;
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
  //if(mode==HORIZONTAL) { start_mod=0; stop_mod=7; }
  if(mode==HORIZONTAL) { start_mod=0; stop_mod=17; }//activate all
  //else if(mode==VERTICAL) { start_mod=7; stop_mod=14; }
  else if(mode==VERTICAL) { start_mod=0; stop_mod=17; }//activate all
  else if(mode==PROTON) { start_mod=0; stop_mod=17; }
  else if(mode==WAGASCI) { start_mod=0; stop_mod=21; }
  else if(mode==WAGASCIBG) { start_mod=0; stop_mod=21; }
  else { start_mod =0; stop_mod=0; }

  //rotation reverse junnjo
	G4RotationMatrix *xrot  = new G4RotationMatrix(G4ThreeVector(1,0,0),90.*degree);
	G4RotationMatrix *xrot2 = new G4RotationMatrix(G4ThreeVector(1,0,0),-90.*degree);
	G4RotationMatrix *yrot  = new G4RotationMatrix(G4ThreeVector(0,1,0),90.*degree);
	G4RotationMatrix *yrot2 = new G4RotationMatrix(G4ThreeVector(0,1,0),-90.*degree);

	//for Y grid layer
	G4RotationMatrix *rotgridY_v = new G4RotationMatrix(G4ThreeVector(0,0,1),180.*degree);
	rotgridY_v->rotateY(-90*degree);
	rotgridY_v->rotateX(-90*degree);
	G4RotationMatrix *rotgridY_h = new G4RotationMatrix(G4ThreeVector(1,0,0),-90.*degree);
	rotgridY_h->rotateY(-90*degree);

	//for Xgrid layer
	G4RotationMatrix *rotgridX_v = new G4RotationMatrix(G4ThreeVector(0,1,0),-90.*degree);
	rotgridX_v->rotateX(-90*degree);
	G4RotationMatrix *rotgridX_h = new G4RotationMatrix(G4ThreeVector(0,0,1),180.*degree);
	rotgridX_h->rotateX(-90*degree);
	rotgridX_h->rotateY(-90*degree);



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
  G4VPhysicalVolume* horizontalPV  = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),horizontalLV,"horizontal",worldLV,false,0);

  //Vertical volume (for the 7 vertical modules)
  G4Box* verticalHall_box =  new G4Box("VerticalHall",HorizonY, HorizonX, HorizonZ);
  G4LogicalVolume* verticalLV =new G4LogicalVolume(verticalHall_box, Air, "vertical_log",0,0,0);
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

  //Set 7 horizontal modules:
  if( k>=0 && k<7 ) {
    moduleLV[k] = new G4LogicalVolume(module_box, Air, "module_log_h");
    char moduleName[14];
    sprintf(moduleName,"module%d",k);

    G4double x = ModuleStart + ModuleSpace*k;
    new G4PVPlacement(0,G4ThreeVector(x,0.,0.),moduleLV[k],moduleName,horizontalLV,false,k);
  }


  //Set 7 vertical modules
  else if( k>=7 && k<14 ) {
    moduleLV[k] = new G4LogicalVolume(module_box, Air, "module_log_v");
    char moduleName[14];
    sprintf(moduleName,"module%d",k);
    G4double y = ModuleStart + ModuleSpace*(k-7);    
    new G4PVPlacement(0,G4ThreeVector(0.,y,0.),moduleLV[k],moduleName,verticalLV,false,k);     
  }



  //Set proton module added by kikawa  common for WAGASCI
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


  // Water tank (SUS304) BLOCK =========================================  

  G4Box* SUS_box = new G4Box("SUS_box",(62.8+1.2)*cm,(62.8+1.6)*cm,(25.+0.4)*cm);	//2015/6/15
  G4LogicalVolume *SUSLV = new G4LogicalVolume(SUS_box,SUS304,"SUSLV");
  G4VisAttributes* SUSVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
  SUSLV->SetVisAttributes(SUSVisAtt);


  // Water BLOCK =========================================  

  //G4Box* Water_box = new G4Box("Water_box",62.8*cm,62.8*cm,23.3*cm);	//2015/6/15
  G4Box* Water_box = new G4Box("Water_box",62.8*cm,62.8*cm,25.*cm);	//2017/03/15
  G4LogicalVolume *WaterLV = new G4LogicalVolume(Water_box,Water,"WaterLV");
  G4VisAttributes* WaterVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
  WaterLV->SetVisAttributes(WaterVisAtt);


  // Al BLOCK ============================================
  G4Box* Al_box = new G4Box("Al_box",72.8*cm,72.8*cm,30.3*cm);	//2015/6/15
  G4LogicalVolume *AlLV = new G4LogicalVolume(Al_box,Al,"AlLV");
  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
  AlLV->SetVisAttributes(AlVisAtt);


  //SCINTILLATORS FOR TRACKING PLANES =======================================

  // vertical-layer scintillator (internal)
  std::vector<G4TwoVector> vdim;
  // horizontal-layer scintillator (internal)
  std::vector<G4TwoVector> hdim;

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
  }

 
  //SciBar Scintillator dimension
/*//original
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
*/
//kikawa rev
  SciVertex2_x[0] = -11.672*mm; 
  SciVertex2_x[1] = -12.21*mm; 
  SciVertex2_x[2] = -12.21*mm; 
  SciVertex2_x[3] = -11.672*mm; 
  SciVertex2_x[4] =  11.672*mm; 
  SciVertex2_x[5] =  12.21*mm; 
  SciVertex2_x[6] =  12.21*mm; 
  SciVertex2_x[7] =  11.672*mm;

  SciVertex2_y[0] = - 6.17*mm; 
  SciVertex2_y[1] = - 3.5*mm; 
  SciVertex2_y[2] =   3.5*mm; 
  SciVertex2_y[3] =   6.17*mm; 
  SciVertex2_y[4] =   6.17*mm; 
  SciVertex2_y[5] =   3.5*mm; 
  SciVertex2_y[6] =  -3.5*mm; 
  SciVertex2_y[7] =  -6.17*mm; 




  //for SciBar type
  for(int iver=0; iver<nSciVertex; iver++){
    svdim.push_back( G4TwoVector( SciVertex2_x[iver], SciVertex2_y[iver] ) );
    shdim.push_back( G4TwoVector( SciVertex2_y[nSciVertex-1-iver], 
    SciVertex2_x[nSciVertex-1-iver] ) );
  }
 

  //INGRID Scintillator dimension for PM
  SciVertex3_x[0] = -23.616*mm; 
  SciVertex3_x[1] = -24.389*mm; 
  SciVertex3_x[2] = -24.389*mm; 
  SciVertex3_x[3] = -23.616*mm; 
  SciVertex3_x[4] =  23.616*mm; 
  SciVertex3_x[5] =  24.389*mm; 
  SciVertex3_x[6] =  24.389*mm; 
  SciVertex3_x[7] =  23.616*mm;

  SciVertex3_y[0] = - 4.71*mm; 
  SciVertex3_y[1] = - 0.5 *mm; 
  SciVertex3_y[2] =   0.5 *mm; 
  SciVertex3_y[3] =   4.71*mm; 
  SciVertex3_y[4] =   4.71*mm; 
  SciVertex3_y[5] =   0.5 *mm; 
  SciVertex3_y[6] =  -0.5 *mm; 
  SciVertex3_y[7] =  -4.71*mm; 


  //for INGRID type of PM
  for(int iver=0; iver<nSciVertex; iver++){
    ivdim.push_back( G4TwoVector( SciVertex3_x[iver], SciVertex3_y[iver] ) );
    ihdim.push_back( G4TwoVector( SciVertex3_y[nSciVertex-1-iver], 
    SciVertex3_x[nSciVertex-1-iver] ) );
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
  //vscint_intVisAtt->SetForceWireframe(true);
  vscint_intVisAtt->SetForceSolid(true);
  vscint_intLV->SetVisAttributes(vscint_intVisAtt);
  
  //INGRID horizontal scintillator
  G4ExtrudedSolid* hscint_tmp = new G4ExtrudedSolid("hscint_tmp", hdim, zsec);
  G4EllipticalTube* hsci_hole = new G4EllipticalTube("hsci_hole",1.95*mm,1.3*mm,600*mm);
  G4SubtractionSolid* hscint_int = new G4SubtractionSolid("hscint_int",  hscint_tmp, hsci_hole);
  G4LogicalVolume* hscint_intLV = new G4LogicalVolume(hscint_int,Scinti,"hscint_intLV");
  G4VisAttributes* hscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hscint_intVisAtt->SetForceWireframe(true);
  hscint_intVisAtt->SetForceSolid(true);
  hscint_intLV->SetVisAttributes(hscint_intVisAtt);


  //SciBar type for proton module added by kikawa
  //SciBar vertical scintillator
  G4ExtrudedSolid* vscint2_tmp = new G4ExtrudedSolid("vscint2_tmp", svdim, zsec);
  G4EllipticalTube* vsci_hole2 = new G4EllipticalTube("vsci_hole2",0.9*mm,0.9*mm,600*mm);
  G4SubtractionSolid* vscint2_int = new G4SubtractionSolid("vscint2_int",  vscint2_tmp, vsci_hole2);
  G4LogicalVolume* vscint2_intLV = new G4LogicalVolume(vscint2_int,Scinti,"vscint2_intLV");
  G4VisAttributes* vscint2_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));
  //vscint2_intVisAtt->SetForceWireframe(true);
  vscint2_intVisAtt->SetForceSolid(true);
  vscint2_intLV->SetVisAttributes(vscint2_intVisAtt);

  //SciBar horizontal scintillator
  G4ExtrudedSolid* hscint2_tmp = new G4ExtrudedSolid("hscint2_tmp", shdim, zsec);
  G4EllipticalTube* hsci_hole2 = new G4EllipticalTube("hsci_hole2",0.9*mm,0.9*mm,600*mm);
  G4SubtractionSolid* hscint2_int = new G4SubtractionSolid("hscint2_int",  hscint2_tmp, hsci_hole2);
  G4LogicalVolume* hscint2_intLV = new G4LogicalVolume(hscint2_int,Scinti,"hscint2_intLV");
  G4VisAttributes* hscint2_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
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



  //SCINTILLATORS FOR WAGASCI ===============================================


  //for WAGASCI simple box
  G4Box* vlayer6 = new G4Box("vlayer6",1.25*cm,0.50*m,0.15*cm);
  //G4Box* vlayer6 = new G4Box("vlayer6",Vlayer6X,Vlayer6Y,Vlayer6Z);//2014/6/11
  G4LogicalVolume* vlayerLV6 = new G4LogicalVolume(vlayer6,Scinti,"vlayerLV6");
  G4VisAttributes* vlayerVisAtt6 = new G4VisAttributes(G4Color(0.,0.7,0.));
  vlayerVisAtt6->SetForceWireframe(true);
  vlayerLV6->SetVisAttributes(vlayerVisAtt6);
  
  //for WAGASCI simple box
  G4Box* hlayer6 = new G4Box("hlayer6",0.50*m,1.25*cm,0.15*cm);
  //G4Box* hlayer6 = new G4Box("hlayer6",Hlayer6X,Hlayer6Y,Hlayer6Z);//2014/6/11
  G4LogicalVolume* hlayerLV6 = new G4LogicalVolume(hlayer6,Scinti,"hlayerLV6");
  G4VisAttributes* hlayerVisAtt6 = new G4VisAttributes(G4Color(0.,0.7,0.));
  hlayerVisAtt6->SetForceWireframe(true);
  hlayerLV6->SetVisAttributes(hlayerVisAtt6);

  //Scintillator dimension for WAGASCI
  std::vector<G4TwoVector> wvdim;
  std::vector<G4TwoVector> whdim;

//  SciVertex4_x[0] = -10.8*mm; 
//  SciVertex4_x[1] = -12.0*mm; 
//  SciVertex4_x[2] = -12.0*mm; 
//  SciVertex4_x[3] = -10.8*mm; 
//  SciVertex4_x[4] =  10.8*mm; 
//  SciVertex4_x[5] =  12.0*mm; 
//  SciVertex4_x[6] =  12.0*mm; 
//  SciVertex4_x[7] =  10.8*mm;
//
//  SciVertex4_y[0] = - 1.5*mm;
//  SciVertex4_y[1] = - 1.0*mm; 
//  SciVertex4_y[2] =   1.0*mm; 
//  SciVertex4_y[3] =   1.5*mm; 
//  SciVertex4_y[4] =   1.5*mm; 
//  SciVertex4_y[5] =   1.0*mm; 
//  SciVertex4_y[6] =  -1.0*mm; 
//  SciVertex4_y[7] =  -1.5*mm; 

  SciVertex4_x[0] = -10.4*mm; 
  SciVertex4_x[1] = -11.8*mm; 
  SciVertex4_x[2] = -11.8*mm; 
  SciVertex4_x[3] = -10.4*mm; 
  SciVertex4_x[4] =  10.9*mm; 
  SciVertex4_x[5] =  11.8*mm; 
  SciVertex4_x[6] =  11.8*mm; 
  SciVertex4_x[7] =  10.9*mm;

  SciVertex4_y[0] = - 1.5*mm;
  SciVertex4_y[1] = - 0.75*mm; 
  SciVertex4_y[2] =   0.75*mm; 
  SciVertex4_y[3] =   1.5*mm; 
  SciVertex4_y[4] =   1.5*mm; 
  SciVertex4_y[5] =   1.0*mm; 
  SciVertex4_y[6] =  -1.0*mm; 
  SciVertex4_y[7] =  -1.5*mm; 





  //for WAGASCI
  for(int iver=0; iver<nSciVertex; iver++){
    wvdim.push_back( G4TwoVector( SciVertex4_x[iver], SciVertex4_y[iver] ) );
    whdim.push_back( G4TwoVector( SciVertex4_y[nSciVertex-1-iver], SciVertex4_x[nSciVertex-1-iver] ) );
  }

  std::vector<G4ExtrudedSolid::ZSection> wzsec;
  wzsec.push_back( G4ExtrudedSolid::ZSection(-500*mm, G4TwoVector(0*mm,0*mm), 1) );
  wzsec.push_back( G4ExtrudedSolid::ZSection( 500*mm, G4TwoVector(0*mm,0*mm), 1) );
  
  //vertical scintillator
  G4ExtrudedSolid* vwscint_tmp = new G4ExtrudedSolid("vwscint_tmp", wvdim, wzsec);
  G4Box* vglue = new G4Box("vglue",0.61*mm,0.61*mm,500.1*mm); //add 0.1mm to avoid bug
  G4SubtractionSolid* vwscint_int = new G4SubtractionSolid("vwscint_int", vwscint_tmp, vglue, 0, G4ThreeVector(-3.9*mm,-0.9*mm,0*mm));
  G4LogicalVolume* vwscint_intLV = new G4LogicalVolume(vwscint_int,Scinti,"vwscint_intLV");
  G4VisAttributes* vwscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  vwscint_intVisAtt->SetForceSolid(true);
  vwscint_intLV->SetVisAttributes(vwscint_intVisAtt);
  
  //horizontal scintillator
  G4ExtrudedSolid* hwscint_tmp = new G4ExtrudedSolid("hwscint_tmp", whdim, wzsec);
  G4Box* hglue = new G4Box("hglue",0.61*mm,0.61*mm,500.1*mm); //add 0.1mm to avoid bug
  G4SubtractionSolid* hwscint_int = new G4SubtractionSolid("hwscint_int", hwscint_tmp, hglue, 0, G4ThreeVector(0.9*mm,-3.9*mm,0*mm));
  G4LogicalVolume* hwscint_intLV = new G4LogicalVolume(hwscint_int,Scinti,"hwscint_intLV");
  G4VisAttributes* hwscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  hwscint_intVisAtt->SetForceSolid(true);
  hwscint_intLV->SetVisAttributes(hwscint_intVisAtt);

 
  //vertical grid scintillator
  //G4VSolid* vwgridscint_int = new G4Box("vwgridscint_int",12*mm,1.5*mm,500*mm);
  //G4VSolid* vwgridscint_int = new G4ExtrudedSolid("vwgridscint_int", wvdim, wzsec);
  G4VSolid* vwgridscint_tmp = new G4ExtrudedSolid("vwgridscint_tmp", wvdim, wzsec);
  G4VSolid* vwgridscint_int = new G4SubtractionSolid("vwgridscint_int", vwgridscint_tmp, vglue, 0, G4ThreeVector(-3.9*mm,-0.9*mm,0*mm));  
  G4VSolid* vcut = new G4Box("vcut",13.1/2.*mm,3.1/2.*mm,3.5/2.*mm); //add 0.1mm to avoid bug
  for(int i=0;i<20;i++){
    //  	G4VSolid* vwgridscint_tmp2 = new G4SubtractionSolid("", vwgridscint_int, vcut, 0, G4ThreeVector(5.5*mm,0*mm,(loli_offsetxy_grid + loli_cutwidth/2. + i*loli_cutgap)*cm));
    G4VSolid* vwgridscint_tmp2 = new G4SubtractionSolid("", vwgridscint_int, vcut, 0, G4ThreeVector(5.5*mm,0*mm,(loli_offsetxy_grid + loli_scinti_thick/2. + i*loli_cutgap)*cm));
    vwgridscint_int = vwgridscint_tmp2;
  }
  G4LogicalVolume* vwgridscint_intLV     = new G4LogicalVolume(vwgridscint_int,Scinti,"vwgridscint_intLV");
  G4VisAttributes* vwgridscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  vwgridscint_intVisAtt->SetForceSolid(true);
  vwgridscint_intLV->SetVisAttributes(vwgridscint_intVisAtt);
 
  //horizontal grid scintillator
  //G4VSolid* hwgridscint_int = new G4Box("hwgridscint_int",1.5*mm,12*mm,500*mm);
  //G4VSolid* hwgridscint_int = new G4ExtrudedSolid("hwgridscint_int", whdim, wzsec);
  G4VSolid* hwgridscint_tmp = new G4ExtrudedSolid("hwgridscint_tmp", whdim, wzsec);
  G4VSolid* hwgridscint_int = new G4SubtractionSolid("hwgridscint_int", hwgridscint_tmp, hglue, 0, G4ThreeVector(0.9*mm,-3.9*mm,0*mm));
  G4VSolid* hcut = new G4Box("hcut",3.1/2.*mm,13./2.*mm,3.5/2.*mm); //add 0.1mm to avoid bug
  for(int i=0;i<20;i++){
        G4VSolid* hwgridscint_tmp2 = new G4SubtractionSolid("", hwgridscint_int, hcut, 0, G4ThreeVector(0*mm,5.5*mm,(loli_offsetxy_grid + loli_scinti_thick/2. + i*loli_cutgap)*cm));
  	hwgridscint_int = hwgridscint_tmp2;
  }
  G4LogicalVolume* hwgridscint_intLV = new G4LogicalVolume(hwgridscint_int,Scinti,"hwgridscint_intLV");
  G4VisAttributes* hwgridscint_intVisAtt = new G4VisAttributes(G4Color(0.,0.,1.));
  hwgridscint_intVisAtt->SetForceSolid(true);
  hwgridscint_intLV->SetVisAttributes(hwgridscint_intVisAtt);

 
  
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
      new G4PVPlacement(0,G4ThreeVector(0.,0.,z),ironLV,ironname,moduleLV[k],false,i+k*9 ); // the world in Module     
    }

    // 11 planes of scintillator per module-----------------------------------------
    for(int i=0;i<11;i++){
      G4double z = scibar_start2 + scibar_z + Gap2*i;
      for(int j=0;j<24;j++){ 
	char name[30]; 

      sprintf(name,"vlayer[%d][%d][%d]",k,i,j);
      G4double x_scinti =  scibar_xy_start + 2*j*scibar_y; 
      //new G4PVPlacement(0,G4ThreeVector(x_scinti,0.,z),vlayerLV,name,moduleLV[k],false,j+i*24+k*264); // in Module        
      new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z),vscint_intLV,name,moduleLV[k],false,j+i*24+k*264); // in Module        
      //new G4PVPlacement(xrot,G4ThreeVector(x_scinti,0.,z),vscint3_intLV,name,moduleLV[k],false,j+i*24+k*264); // in Module        

      sprintf(name,"hlayer[%d][%d][%d]",k,i,j);
      G4double y_scinti =  scibar_xy_start + 2*j*scibar_y;     
      //new G4PVPlacement(0,G4ThreeVector(0.,y_scinti,z+2*scibar_z),hlayerLV,name,moduleLV[k],false,j+i*24+k*264); // in Module
      new G4PVPlacement(yrot,G4ThreeVector(0.,y_scinti,z+2*scibar_z+0.1*mm),hscint_intLV,name,moduleLV[k],false,j+i*24+k*264); // in Module
      //new G4PVPlacement(yrot,G4ThreeVector(0.,y_scinti,z+2*scibar_z),hscint3_intLV,name,moduleLV[k],false,j+i*24+k*264); // in Module

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
  
    //else if(k==16){
    else if(k==16 && mode==PROTON){
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
    }//else(i!=0)

  }//for i

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
   
  }//else if(proton)
 

    //WAGASCI prototype added by koga 2015/1/23
  
    else if(k==20 && mode==WAGASCI){

	//2017/03/15 
	new G4PVPlacement(0,G4ThreeVector(0,0,0),SUSLV,  "sustarget"  ,moduleLV[16],false,0,1); //water tank
	new G4PVPlacement(0,G4ThreeVector(0,0,0),WaterLV,"watertarget",SUSLV       ,false,0,1); //water target


        INGRID_Dimension *fdim = new INGRID_Dimension();
	for(int i=0;i<8;i++){
              for(int j=0;j<40;j++){ 
      		char name[100];
		double x,y,z;

		//vertical layer
        	fdim -> get_pos_loli( 15, 1, i, j, 0, &x, &y, &z); //mod view pln ch grid
                sprintf(name,"vlayer[%d][%d][%d]",k,i,j);
                //new G4PVPlacement(0,G4ThreeVector(x*cm,y*cm,z*cm), vlayerLV6,name,WaterLV,false,j + (i)*CHMAX + k*264,1); //test by simple box (can be visualized) 
                new G4PVPlacement(xrot2,G4ThreeVector(x*cm,y*cm,z*cm),vwscint_intLV,name,WaterLV,false,j + (i)*CHMAX + k*264);        

		//horizontal layer
        	fdim -> get_pos_loli( 15, 0, i, j, 0, &x, &y, &z); //mod view pln ch grid
                sprintf(name,"hlayer[%d][%d][%d]",k,i,j); 
                //new G4PVPlacement(0,G4ThreeVector(x*cm,y*cm,z*cm), hlayerLV6,name,WaterLV,false,j + (i)*CHMAX + k*264,1); //test by simple box (can be visualized)
                new G4PVPlacement(yrot2,G4ThreeVector(x*cm,y*cm,z*cm),hwscint_intLV,name,WaterLV,false,j + (i)*CHMAX + k*264);

		if(j<20){
			//vertical grid layer
        		fdim -> get_pos_loli( 15, 1, i, j, 1, &x, &y, &z); //mod view pln ch grid
        	        sprintf(name,"vlayer1[%d][%d][%d]",k,i,j);
        	        new G4PVPlacement(rotgridX_v,G4ThreeVector(x*cm,y*cm,z*cm),vwgridscint_intLV,name,WaterLV,false, 40 + j + (i)*CHMAX + k*264);
        		fdim -> get_pos_loli( 15, 1, i, j, 2, &x, &y, &z); //mod view pln ch grid
        	        sprintf(name,"vlayer2[%d][%d][%d]",k,i,j);
        	        new G4PVPlacement(rotgridY_v,G4ThreeVector(x*cm,y*cm,z*cm),vwgridscint_intLV,name,WaterLV,false, 60 + j + (i)*CHMAX + k*264);

			//horizontal grid layer
        		fdim -> get_pos_loli( 15, 0, i, j, 1, &x, &y, &z); //mod view pln ch grid
        	        sprintf(name,"hlayer1[%d][%d][%d]",k,i,j); 
        	        new G4PVPlacement(rotgridX_h,G4ThreeVector(x*cm,y*cm,z*cm),hwgridscint_intLV,name,WaterLV,false, 40 + j + (i)*CHMAX + k*264);
        		fdim -> get_pos_loli( 15, 0, i, j, 2, &x, &y, &z); //mod view pln ch grid
        	        sprintf(name,"hlayer2[%d][%d][%d]",k,i,j); 
        	        new G4PVPlacement(rotgridY_h,G4ThreeVector(x*cm,y*cm,z*cm),hwgridscint_intLV,name,WaterLV,false, 60 + j + (i)*CHMAX + k*264);
		}
              }
        }

	    
     }//else if(wagasci)
 
    //WAGASCI prototype BG study added by koga 2015/1/23
    else if(k==20 && mode==WAGASCIBG){
	new G4PVPlacement(0,G4ThreeVector(0,0,0),AlLV,"alflame",moduleLV[16],false,0,1); //al frame
	new G4PVPlacement(0,G4ThreeVector(0,0,0),WaterLV,"watertarget",AlLV,false,0,1); //water target
    }

  }//for k







  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  G4String hlayerSDname = "Ingrid/hlayerSD";
  ahlayerSD = new IngridHLayerSD( hlayerSDname );
  SDman->AddNewDetector( ahlayerSD );
  hscint_intLV->SetSensitiveDetector ( ahlayerSD );
  hscint2_intLV->SetSensitiveDetector( ahlayerSD );//SciBar type added by kikawa
  hscint3_intLV->SetSensitiveDetector( ahlayerSD );//INGRID type added by kikawa
  hlayerLV6->SetSensitiveDetector    ( ahlayerSD );    //wagasci simple box
  hwscint_intLV->SetSensitiveDetector( ahlayerSD ); //wagasci scinti h
  hwgridscint_intLV->SetSensitiveDetector( ahlayerSD ); //wagasci scinti h
 
  G4String vlayerSDname = "Ingrid/vlayerSD";
  avlayerSD = new IngridVLayerSD( vlayerSDname );
  SDman->AddNewDetector( avlayerSD );
  vscint_intLV->SetSensitiveDetector ( avlayerSD );
  vscint2_intLV->SetSensitiveDetector( avlayerSD );//SciBar type added by kikawa
  vscint3_intLV->SetSensitiveDetector( avlayerSD );//INGRID type added by kikawa
  vlayerLV6->SetSensitiveDetector    ( avlayerSD );    //wagasci simple box
  vwscint_intLV->SetSensitiveDetector( avlayerSD ); //wagasci scinti v
  vwgridscint_intLV->SetSensitiveDetector( avlayerSD ); //wagasci scinti v
 
  G4String vetoSDname = "Ingrid/vetoSD";
  avetoSD = new IngridVetoSD( vetoSDname );
  SDman->AddNewDetector( avetoSD );

  SvetoLV->SetSensitiveDetector ( avetoSD );
  LvetoLV->SetSensitiveDetector ( avetoSD );
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
  WorldY=15.0*m;
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
  //PHorizonZ=45.0*cm;
  PHorizonZ=48.0*cm;

  //Module volume
  //ModuleX=0.82*m;  // overlap of logical volumes is not good !!!
  //ModuleY=0.82*m;
  ModuleX=0.75*m;
  ModuleY=0.75*m;
  //ModuleZ=70.0*cm;
  ModuleZ=65.0*cm;

  //Proton Module volume
  //PModuleZ=44.0*cm;
  PModuleZ=47.0*cm;
  
  //Positions of modules
  ModuleSpace= 1.50*m;     // 120 + 30  distance between the centers of two consecutive modules
  ModuleStart=-4.50*m;     // -150*3 center of the first module
}

//___________________________________________________________________________________________________________

void IngridDetectorConstruction::DefineStructures()
{
  //Iron blocks
  iron_z= 3.25*cm;
  iron_xy= 62.0*cm;
  Gap= 10.7*cm;           // 6.5 + 2 +2*1.1 means " scinti-scinti-gap-iron-gap-scinti-scinti-... "
  iron_start= -48.15*cm;  // Position of the center of the first iron block  
	    		// = -54.5cm(1st plane surface) + 2cm(width of 1 plane) + 1.1cm(gap of scinti-iron) +3.25cm(half of iron)
  
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

  //B2
  offset_mod5_v = -500.0*mm;
  offset_mod5_h = -500.0*mm;

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

//kikawa rev
  SciVertex_x[0] = -23.616*mm; 
  SciVertex_x[1] = -24.389*mm; 
  SciVertex_x[2] = -24.389*mm; 
  SciVertex_x[3] = -23.616*mm; 
  SciVertex_x[4] =  23.616*mm; 
  SciVertex_x[5] =  24.389*mm; 
  SciVertex_x[6] =  24.389*mm; 
  SciVertex_x[7] =  23.616*mm;

  SciVertex_y[0] = - 4.71*mm; 
  SciVertex_y[1] = - 0.5*mm; 
  SciVertex_y[2] =   0.5*mm; 
  SciVertex_y[3] =   4.71*mm; 
  SciVertex_y[4] =   4.71*mm; 
  SciVertex_y[5] =   0.5*mm; 
  SciVertex_y[6] =  -0.5*mm; 
  SciVertex_y[7] =  -4.71*mm; 



  //original
/*
  SciVertex_x[0] = -23.8*mm; 
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

  
/*  
  SciVertex_x[0] = -24.8*mm; 
  SciVertex_x[1] = -24.9*mm; 
  SciVertex_x[2] = -24.9*mm; 
  SciVertex_x[3] = -24.8*mm; 
  SciVertex_x[4] =  24.8*mm; 
  SciVertex_x[5] =  24.9*mm; 
  SciVertex_x[6] =  24.9*mm; 
  SciVertex_x[7] =  24.8*mm;

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
  a = 58.69*g/mole;
  G4Element* elNi = new G4Element(name="Nickel", symbol="Ni", z=28., a);
  a = 51.99*g/mole;
  G4Element* elCr = new G4Element(name="Chromium", symbol="Cr", z=24., a);

  //Air
  density = 1.29*mg/cm3;
  Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  //Iron
  a = 55.845*g/mole;
  density = 7.86*g/cm3;
  Fe = new G4Material(name="Iron", z=26., a, density);

  //Water
  density = 1.000*g/cm3;
  Water = new G4Material(name="Water",density,nel=2);
  Water->AddElement(elO,1);
  Water->AddElement(elH,2);

  //Al
  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  Al = new G4Material(name="Aluminum", z=13., a, density);

  //Scintillator
  density = 1.032*g/cm3;
  Scinti = new G4Material(name="Scintillator", density, nel=2);
  Scinti->AddElement(elC, 9);
  Scinti->AddElement(elH, 10);

  //Concrete
  density = 2.2*g/cm3;
  Concrete = new G4Material(name="Concrete", density, nel=6);
  Concrete->AddElement(elO, .53);
  Concrete->AddElement(elSi, .335);
  Concrete->AddElement(elCa, 0.06);
  Concrete->AddElement(elNa, 0.015);
  Concrete->AddElement(elFe, 0.02);
  Concrete->AddElement(elAl, 0.04);

  //SUS304
  density = 7.93*g/cm3;
  SUS304 = new G4Material(name="SUS304", density, nel=3);
  SUS304->AddElement(elFe, 0.72);
  SUS304->AddElement(elCr, 0.19);
  SUS304->AddElement(elNi, 0.09);

}

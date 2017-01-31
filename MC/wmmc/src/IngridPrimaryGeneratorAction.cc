#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandGauss.h"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

#include "IngridPrimaryGeneratorAction.hh"
#include "IngridDetectorConstruction.hh"
#include "IngridRunAction.hh"
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <math.h>

#define PROTON 2
#define HORIZONTAL 3
#define VERTICAL 4
#define WAGASCI 5
#define WAGASCIBG 6

#define MOD_CEN_GAP 150.	// the gap of each module center = 60(half of module)*2cm + 30(module gap)cm
#define OFFSET 400.		// Gap between horizotal modules and vertical modules (cm)
#define OFFSETPM 120. 	// Gap between horizotal modules and proton modules (cm)
//#define OFFSETLOLI 120.95 	// Gap between horizotal modules and proton modules (cm)
//#define OFFSETLOLI 144.95 	//120.95 + 24.0 Gap between horizotal modules and proton modules (cm)
#define OFFSETLOLI 143.3 	//120. + 23.3 Gap between horizotal modules and proton modules (cm)
#define OFFSETLOLIBG 150.3 	//120. + 30.3 Gap between horizotal modules and proton modules (cm)

// IRON 120x120x6.5 -------------------------

//

/* new value ?
const double total_mass_fe = 6.6165e3; //kg
const double total_mass_sci = 3.1738e2; //kg
*/

const double total_mass_fe = 99.54; //ton
const double total_mass_sci = 3.74; //ton
const double scinti_start = -54.5;  // the surface of the first scintillator by akira
const double iron_start = scinti_start + 2.0 + 1.1;  // the surface fo the first scinti + 2 scnti + air gap  by akira
const double width_fe = 6.5; // cm
const double width_sci = 2.0; // 1.0cm * 2
const double GAP = 10.7; // 6.5 + 2 + 2.2

const double HallX = -216.7; //cm 
const double HallZ = 170.;    //cm
//const double HallRadiusMax = 1320.;//cm
const double HallRadiusMax = 1150.;//cm
const double HallRadiusMin = 850.;//cm

//for Proton Module added by kikawa
//const double total_mass_sci_pm = 0.556 ;//ton (total mass)
const double total_mass_sci_pm = 0.56904 ;//ton (total mass)
const double total_mass_front_pm = 0.028848 ;//ton (total mass of front veto)
const double ingrid_width = 1 ;//(cm) (total mass of ingrid type)
const double scibar_width = 1.3; //(cm) (total mass of scibar type)
const double width_ingrid =1.0; //INGRID type 
const double width_scibar =1.3; //SciBar type 
const double scibar_region=44.2;//(cm)
const double sciing_region=39.1;//(cm)
const double ingrid_region=34;//(cm)
const double Pscinti_start=-40.95;//(cm)
const double distance_first=2.7;//(cm)
const double distance_pln=2.3;//(cm)
const double diff=-0.15;//(cm) difference of scibar and ingrid type start

//for Mod03
const double GAP03 = GAP - 0.2;
const double scinti_start03 = scinti_start + 1.0;
const double iron_start03 = iron_start + 0.9;

//#define MOD3
//#define DEBUG 
//#define DEBUG2 

//
IngridPrimaryGeneratorAction::
IngridPrimaryGeneratorAction(Neut *neut0,IngridRunAction* rac, IngridEventAction* evt,int nd,int flavor0,int seed=-1)
  :runaction(rac)
{
  eventaction = evt;
  neut_fe = neut0;
  module_mode = nd;
  neutrino_flavor = flavor0;
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  isParticleGun=false;
  particleGun_pdg=0;

  if(seed!=-1)CLHEP::HepRandom::setTheSeed(seed);
  runaction->NotEntry = 0;
}

//
IngridPrimaryGeneratorAction::~IngridPrimaryGeneratorAction()
{
  if( particleGun!=NULL ) { delete particleGun; particleGun=NULL; }
}

// option -p : propagate only a single type of particle
void IngridPrimaryGeneratorAction::SetParticleGun(bool isGun, int pdg){
  isParticleGun=isGun;
  particleGun_pdg=pdg;
}


//
void IngridPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  SecondaryVector Secondary;
  Neut *neut = neut_fe;
  int fdid = 0;
  int mode = 0;
  float pos[3];
  int ID=-1;
  float direction[3];
  int vertex_flag=0;
  int flayer;
  double prob;//for Proton Module
  int scitype;//scintillator type


  //cosmic mode
  if(neutrino_flavor==-1){
      
      G4ParticleDefinition* particle;
      particle = particleTable->FindParticle(13);
      particleGun->SetParticleDefinition(particle);

      //set momentum
      G4float mass  = particle->GetPDGMass();
      //G4float energy = 3*GeV; // refer to PDG2015
      //particleGun->SetParticleEnergy(energy);

      //load momentum distribution file
      //std::ifstream ifs("/export/scbn25/data1/taichiro/T2K/work/basesoft/git/watermodule/wmmc/cosmic_mom.txt");
      //double x[1000],y[1000];
      //int i=0;
      //while(ifs >> x[i] >> y[i]){
      //	i++;
      //}
      //for(int j=0;j<i;j++){
      //	y[j]=y[j]/pow(x[j],2.7);
      //}
      //TGraph *gr = new TGraph(i,x,y);
      //TF1 *f = new TF1("f","[0]+[1]*pow(x,-1)+[2]*pow(x,-2)+[3]*pow(x,-3)+[4]*pow(x,-4)");
      //gr->Fit("f","","",0.5,100);
      //double fitparam[10];
      //for(int j=0;j<5;j++){
      //  fitparam[j]=f->GetParameter(j);
      //}

      double fitparam[5]={0.01371,-2.407,149.7,-247.0,131.9};
      G4float mumom;
      float mumom_calc;
      float rand_mom_calc;
      float func_calc;
      while(1){
        mumom_calc = G4UniformRand()*99+1; //1-100GeV
        func_calc  = fitparam[0] + fitparam[1]*pow(mumom_calc,-1) + fitparam[2]*pow(mumom_calc,-2) + fitparam[3]*pow(mumom_calc,-3) + fitparam[4]*pow(mumom_calc,-4);
        rand_mom_calc = G4UniformRand()*40;
	if(func_calc>rand_mom_calc){
		mumom=mumom_calc*GeV;
		break;
	}
      }
      G4float energy = sqrt(mumom*mumom+mass*mass);
      particleGun->SetParticleEnergy(energy);




      //set theta and position
      float pi_calc = 4.0*atan(1.0);
      float thetamax_calc =70./180.*pi; // reffer to PDG2015
      float rand_theta_calc;
      float theta_calc,cos2theta_calc;
      float phi_calc;
      float xpos_calc;
      float ypos_calc;
      float zpos_calc;
      float px;
      float py;
      float pz;

      while(1){
	//set theta
        while(1){
          theta_calc = G4UniformRand()*pi_calc/2.;
          phi_calc = 2*pi_calc*G4UniformRand();
          cos2theta_calc = cos(theta_calc)*cos(theta_calc);
          rand_theta_calc = G4UniformRand();
          //if(cos2theta_calc>rand_theta_calc && theta_calc<thetamax_calc){
          if(cos2theta_calc>rand_theta_calc){
          	break;
          }
        }
        px = sin(phi_calc)*sin(theta_calc);
        py = -cos(theta_calc);
        pz = cos(phi_calc)*sin(theta_calc);

        //set position
        xpos_calc=(G4UniformRand()-0.5)*2000;
        ypos_calc=450;
        zpos_calc=(G4UniformRand()-0.5)*2000-125;

	//reject far event
	float point_0[3]={0,0,-125};
	float t = px*(point_0[0]-xpos_calc) + py*(point_0[1]-ypos_calc) + pz*(point_0[2]-zpos_calc);
	float point_1[3]={px*t+xpos_calc, py*t+ypos_calc, pz*t+zpos_calc};
	float dist = sqrt( pow(point_0[0]-point_1[0],2) + pow(point_0[1]-point_1[1],2) + pow(point_0[2]-point_1[2],2) );
        //if(1){
        if(dist<sqrt(50*50+50*50+30*30)){
          	break;
        }
      }

      particleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
      G4float xpos = xpos_calc*cm;
      G4float ypos = ypos_calc*cm;
      G4float zpos = zpos_calc*cm;
      particleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
      particleGun->SetParticleTime(0.0*ns);

      IngridSimVertexSummary* simvertex = new IngridSimVertexSummary();
      simvertex -> Clear();
      simvertex -> nutype   = neutrino_flavor;
      simvertex -> inttype  = -100;
      //simvertex -> nuE      = energy/GeV;
      simvertex -> nuE      = mumom/GeV;
      simvertex -> xnu      = xpos/cm;
      simvertex -> ynu      = ypos/cm;
      simvertex -> znu      = zpos/cm;
      simvertex -> gmomx.push_back(px);
      simvertex -> gmomy.push_back(py);
      simvertex -> gmomz.push_back(pz);
      //simvertex -> mod      = ID;
      //simvertex -> norm			= (neut->Vector).Neutrino.Norm;
      runaction  -> GetEvtSum() -> AddSimVertex( simvertex );

      particleGun->GeneratePrimaryVertex(anEvent);
      return;
  }


  // start loop of neut file
  while(1){
    //G4cout << G4UniformRand() << G4endl;
    //break;
    mode = neut->NtupleReadEvent(Secondary,direction);
    if( mode==-1111 ){
      G4cout <<"Aboart Run (mode =" << mode << G4endl;
      G4RunManager* runManager = G4RunManager::GetRunManager();
      eventaction->SetWriteFlag(-1); 
      runManager->AbortRun(1);
      return;
    }

    // Check flavor ====================
    //neutrino mode
    int  neutrino_flavor_tmp =  (int)(((neut->Vector).Neutrino.ProductionMode)/10);
    if( neutrino_flavor_tmp != neutrino_flavor ) {
#ifdef DEBUG
      G4cout << " === This neutrino id : " << neutrino_flavor_tmp
             << " ,not selected now === " << G4endl;
#endif
      continue;
    }

  // Define neutrino interaction vertex
    fdid = (neut->Vector).Neutrino.FDID;

    // X-Y vertex
    pos[0] = (neut->Vector).Neutrino.x;
    pos[1] = (neut->Vector).Neutrino.y;
    //pos[1] = (neut->Vector).Neutrino.y +15.;//for 15cm lower beam center

    

    // Z-Vertex for INGRID Signal MC
    if( fdid==3 || fdid==4 ) {

      // select vertex in Fe or scinti ?
      double ratio = total_mass_sci / (total_mass_fe + total_mass_sci);
      if( ratio > G4UniformRand() ) vertex_flag = 1; // flag = 1 -> vertex in scinti
      else vertex_flag = 0;
      
      //the origin is the center of the horizontal modules.
      // vertcical moduleds center <---4m---> horizontal modules center 
      
      // iron#/module is 0 ~ 8, scinti#/module is 0~10
      // unit of pos is [cm];

      if( fdid==3 && fabs(pos[0]) <= 60 && fabs(pos[1]) <= 60){
	if( vertex_flag==0 ) {
	  flayer = (int)(9*G4UniformRand());      // 0 < 9*rand < 9
	  pos[2] = width_fe*G4UniformRand();
	  pos[2] = pos[2] + iron_start03 + GAP03*flayer;    
	}
	else if( vertex_flag==1 ) {
	  flayer = (int)(11*G4UniformRand());     // 0 < 11*rand < 11
	  pos[2] = width_sci*G4UniformRand();
	  pos[2] = pos[2] + scinti_start03 + GAP03*flayer;
	}
      }
      else{
	if( vertex_flag==0 ) {
	  flayer = (int)(9*G4UniformRand());      // 0 < 9*rand < 9
	  pos[2] = width_fe*G4UniformRand();
	  pos[2] = pos[2] + iron_start + GAP*flayer;    
	}
	else if( vertex_flag==1 ) {
	  flayer = (int)(11*G4UniformRand());     // 0 < 11*rand < 11
	  pos[2] = width_sci*G4UniformRand();
	  pos[2] = pos[2] + scinti_start + GAP*flayer;
	}
      }


    }

    //Proton Module
    if( fdid==2 && module_mode==PROTON) {

      double front =  total_mass_front_pm / (total_mass_front_pm + total_mass_sci_pm);
    
      if( front > (G4UniformRand()) ){
	vertex_flag=0;
	prob=sciing_region/scibar_region;
      }
      else if(fabs(pos[0])<=20&&fabs(pos[1])<=20){
	vertex_flag=1;
	prob=1;
      }
      else if(fabs(pos[0])<=20){
	vertex_flag=2;
	prob=sciing_region/scibar_region;
      }      
      else if(fabs(pos[1])<=20){
	vertex_flag=3;
	prob=sciing_region/scibar_region;
      }      
      else{
	vertex_flag=4;
	prob=ingrid_region/scibar_region;
      }


      if(prob<(G4UniformRand()))continue;

      
      // unit of pos is [cm];
      if( vertex_flag==0 ) {
        flayer = (int)(2*(G4UniformRand()));     // 0 < 2*rand < 2
	pos[2] = width_ingrid*(G4UniformRand());
        pos[2] = pos[2] + Pscinti_start + flayer*distance_pln;    
      }
      else if( vertex_flag==1 ) {
        flayer = (int)(34*(G4UniformRand()));     // 0 < 34*rand < 34
        pos[2] = width_scibar*(G4UniformRand());
        pos[2] = pos[2] + Pscinti_start + distance_pln + diff + distance_first + flayer*distance_pln;
      }
      else if( vertex_flag==2){
	scitype= (int)((G4UniformRand())/scibar_width*(ingrid_width+scibar_width));
	flayer = (int)(17*(G4UniformRand()));     // 0 < 17*rand < 17
	if(scitype==0){
	  pos[2] = width_scibar*(G4UniformRand());
	  pos[2] = pos[2] + Pscinti_start + 2*distance_pln + diff + distance_first + flayer*distance_pln*2;
	}
	else{
	  pos[2] = width_ingrid*(G4UniformRand());
	  pos[2] = pos[2] + Pscinti_start + distance_pln + distance_first + flayer*distance_pln*2;
	}
      }
      else if( vertex_flag==3){
	scitype= (int)((G4UniformRand())/scibar_width*(ingrid_width+scibar_width));
	flayer = (int)(17*(G4UniformRand()));     // 0 < 17*rand < 17
	if(scitype==0){
	  pos[2] = width_scibar*(G4UniformRand());
	  pos[2] = pos[2] + Pscinti_start + distance_pln + diff + distance_first + flayer*distance_pln*2;
	}
	else{
	  pos[2] = width_ingrid*(G4UniformRand());
	  pos[2] = pos[2] + Pscinti_start + 2*distance_pln + distance_first + flayer*distance_pln*2;
	}
      }
      else if( vertex_flag==4){
	flayer = (int)(34*(G4UniformRand()));     // 0 < 34*rand < 34
        pos[2] = width_ingrid*(G4UniformRand());
        pos[2] = pos[2] + Pscinti_start + distance_pln + distance_first + flayer*distance_pln;
      }


    }


    if(fdid==2 && module_mode==WAGASCI){
        pos[2]=46.6*(G4UniformRand());     // 0 < 46.6*rand < 46.6
    }

    if(fdid==3 && module_mode==WAGASCIBG){
        pos[2]=60.6*(G4UniformRand());     // 0 < 60.6*rand < 60.6
    }


    // jnubeam ndid = 7 (INGRID-upstream surface for BG study)
    if(fdid==7){ // jnubeam ndid = 7(INGRID-upstream surface for BG study
      G4double lx = pos[0]- HallX;
      if( fabs(lx) < HallRadiusMin ){
 	//G4cout << HallRadiusMin * HallRadiusMin - lx * lx <<G4endl;
	pos[2] = -1.0 * sqrt( HallRadiusMin * HallRadiusMin - lx * lx ) - G4UniformRand()*(HallRadiusMax-HallRadiusMin)+HallZ;
	//G4cout << pos[2] << G4endl;
      }
      else{
	pos[2] = -1.0 *  G4UniformRand()*(HallRadiusMax-HallRadiusMin)+HallZ;
	//G4cout << pos[2] << G4endl;
	//pos[2] = -1.*G4UniformRand()*(HallRadiusMax-HallRadiusMin)/370. * (1261.7-fabs(lx))+HallZ;
      }
    }






    // define vertex module
    ID = -1;
#ifdef MOD3 
     // for one module test
	if(fabs(pos[0])<60  &&  fabs(pos[1])<60){
		ID = 3;
		goto NEXTSTEP;
	}
#endif

#ifndef MOD3 
    if( fdid == 7 ) goto NEXTSTEP;

    //


    if( fdid==2 && module_mode == PROTON ) {
      for( int m=0;m<1;m++ ) {
        if( fabs(pos[0]) <= 60 &&
            fabs(pos[1]) <= 60 ) {
          ID = m+16;
          goto NEXTSTEP;
        }
      }
    }
    //if( fdid==3 && module_mode == HORIZONTAL ) {
    if( fdid==3 && (module_mode == HORIZONTAL || module_mode == PROTON || module_mode == WAGASCI ) ) { //changed 2016/6/30 for bg study
      for( int m=0;m<7;m++ ) {
        if( fabs(pos[0]+MOD_CEN_GAP*(3-m)) <= 60 &&
            fabs(pos[1]) <= 60 ) {
          ID = m;
          goto NEXTSTEP;
        }
      }
    }
    if( fdid==4 && ( module_mode == VERTICAL || module_mode == PROTON || module_mode == WAGASCI ) ) {//changed ML 2017/01/25 for BkgV study
      for( int m=0;m<7;m++ ) {
        if( fabs(pos[0]) <= 60 &&
            fabs(pos[1]+MOD_CEN_GAP*(3-m)) <= 60 ) {
          ID = m+7;
          goto NEXTSTEP;
        }
      }
    }

    if( fdid==2 && module_mode == WAGASCI ) {
      for( int m=0;m<1;m++ ) {
        if( fabs(pos[0]) <= 125.6/2. &&
            fabs(pos[1]) <= 125.6/2. ) {
          ID = m+17;
          goto NEXTSTEP;
        }
      }
    }

    if( fdid==3 && module_mode == WAGASCIBG ) {
      for( int m=0;m<1;m++ ) {
        if( fabs(pos[0]) <= 145.6/2. &&
            fabs(pos[1]) <= 145.6/2. ) {
          ID = m+18;
          goto NEXTSTEP;
        }
      }
    }

#endif

    // count events which have vertex out of modules
#ifdef DEBUG
    G4cout << "##### Interaction vertex is out of modules #####" << G4endl;
    G4cout << "##### Skip this event                      #####" << G4endl;
#endif
    runaction->NotEntry++; 

  } // end while loop

//
NEXTSTEP:

// Gap between horizotal modules and vertical modules (cm) 
    if(ID > 6 && ID < 14)  pos[2] = pos[2] - OFFSET;
// Gap between horizotal modules and proton modules (cm) 
    else if(ID == 16)  pos[2] = pos[2] - OFFSETPM;
    else if(ID == 17)  pos[2] = pos[2] - OFFSETLOLI;
    else if(ID == 18)  pos[2] = pos[2] - OFFSETLOLIBG;

    // Input Neut file info to output ROOT class
    neut->ID = ID;
    for(int i=0;i<3;i++) (runaction->vertex)[i] = pos[i];

    IngridSimVertexSummary* simvertex = new IngridSimVertexSummary();
    simvertex -> Clear();
    simvertex -> nutype   = neutrino_flavor;
    simvertex -> inttype  = (neut->Vector).Primary.Mode;
    simvertex -> nuE      = (neut->Vector).Neutrino.Energy;
    simvertex -> xnu      = pos[0];
    simvertex -> ynu      = pos[1];
    simvertex -> znu      = pos[2];
    simvertex -> mod      = ID;
    simvertex -> norm			= (neut->Vector).Neutrino.Norm;
    simvertex -> totcrsne	= (neut->Vector).neutcrs.Totcrsne;

    //for Al density added by koga
    if(ID==18){
	simvertex -> norm			= (neut->Vector).Neutrino.Norm * 2.7;
    }

    //removed for 11b (before using t2kreweight)                                                
    /*
    simvertex -> ng = (neut->Vector).Neutrino.ancestor.ng;
    for(int i=0;i<(neut->Vector).Neutrino.ancestor.ng;i++) {
      simvertex -> gpid.push_back   ( (neut->Vector).Neutrino.ancestor.gpid[i] );
      simvertex -> gmec.push_back   ( (neut->Vector).Neutrino.ancestor.gmec[i] );
      simvertex -> gcosbm.push_back ( (neut->Vector).Neutrino.ancestor.gcosbm[i] );
      simvertex -> gposx.push_back  ( (neut->Vector).Neutrino.ancestor.gvx[i] );
      simvertex -> gposy.push_back  ( (neut->Vector).Neutrino.ancestor.gvy[i] );
      simvertex -> gposz.push_back  ( (neut->Vector).Neutrino.ancestor.gvz[i] );
      simvertex -> gmomx.push_back  ( (neut->Vector).Neutrino.ancestor.gpx[i] );
      simvertex -> gmomy.push_back  ( (neut->Vector).Neutrino.ancestor.gpy[i] );
      simvertex -> gmomz.push_back  ( (neut->Vector).Neutrino.ancestor.gpz[i] );
    }
    */

    runaction  -> GetEvtSum() -> AddSimVertex( simvertex );

    G4cout.precision( 3 );

#ifdef DEBUG
    G4cout << "\n=== Neutrino Information from Jnubeam ===" << G4endl;
		G4cout << "Norm: " <<  (neut->Vector).Neutrino.Norm << G4endl;
		G4cout << "Totcrsne: " <<  (neut->Vector).neutcrs.Totcrsne << G4endl;
    G4cout << "ParentID: " << (neut->Vector).Neutrino.ParentID;
    G4cout << "  Neut Production Mode: " << (neut->Vector).Neutrino.ProductionMode;
    G4cout << "  Neutrino.FDID: " << (neut->Vector).Neutrino.FDID << G4endl;
    G4cout << "Neut interaction Mode: " << (neut->Vector).Primary.Mode << G4endl;
    G4cout << "Energy[GeV]: " << (neut->Vector).Neutrino.Energy;
    G4cout << "  Direction: {" << direction[0] << "," << direction[1] << "," << direction[2] << "}" << G4endl;
    G4cout << "Vertex(cm): {" << pos[0] << ", "<< pos[1] << ", "<< pos[2] << "}";
    G4cout << "  Module: " << ID << "\n\n";
#endif

    particleGun->SetParticlePosition(G4ThreeVector(pos[0]*cm,pos[1]*cm,pos[2]*cm));

    // Input Neut info for T2KReWeight to SK__h1 class
    runaction -> numnu = (neut->Vector).Primary.NumParticle;
    runaction -> mode  = (neut->Vector).Primary.Mode;
    for ( int i = 0; i<50; i++ ) {
      runaction -> ipnu[i] = (neut->Vector).Primary.ParticleID[i];
      runaction -> pnu[i] = (neut->Vector).Primary.AbsMomentum[i];
      for ( int j = 0 ; j < 3 ; j++ ){
        runaction -> dirnu[i][j] = (neut->Vector).Primary.Momentum[i][j] / (neut->Vector).Primary.AbsMomentum[i];
      }
    }

    runaction -> Crsx   = (neut->Vector).Crs.Crsx;
    runaction -> Crsy   = (neut->Vector).Crs.Crsy;
    runaction -> Crsz   = (neut->Vector).Crs.Crsz;
    runaction -> Crsphi = (neut->Vector).Crs.Crsphi;

    runaction -> Nvert = (neut->Vector).Fsi.Nvert;
    for (int ivert=0; ivert<150; ivert++) {
      runaction -> Iflgvert[ivert] = (neut->Vector).Fsi.Iflgvert[ivert];
      for (int j=0; j<3; j++)
        runaction -> Posvert[ivert][j] = (neut->Vector).Fsi.Posvert[ivert][j];
    }

    runaction -> Nvcvert = (neut->Vector).Fsi.Nvcvert;
    for (int ip=0; ip<900; ip++) {

      runaction -> Abspvert[ip]  = (neut->Vector).Fsi.Abspvert[ip];
      runaction -> Abstpvert[ip] = (neut->Vector).Fsi.Abstpvert[ip];
      runaction -> Ipvert[ip]    = (neut->Vector).Fsi.Ipvert[ip];
      runaction -> Iverti[ip]    = (neut->Vector).Fsi.Iverti[ip];
      runaction -> Ivertf[ip]    = (neut->Vector).Fsi.Ivertf[ip];

      for (int j=0; j<3; j++)
        runaction -> Dirvert[ip][j] = (neut->Vector).Fsi.Dirvert[ip][j];
    }

    runaction -> Fsiprob = (neut->Vector).Fsi.Fsiprob;
    runaction -> Numbndn = (neut->Vector).target_info.Numbndn;
    runaction -> Numbndp = (neut->Vector).target_info.Numbndp;
    runaction -> Numfrep = (neut->Vector).target_info.Numfrep;
    runaction -> Numatom = (neut->Vector).target_info.Numatom;
    runaction -> Ibound  = (neut->Vector).Fsi.Ibound;
    runaction -> Npvc    = (neut->Vector).Secondary.NumParticle;
    for (int i=0; i<100; i++) {
      runaction -> Ipvc[i]    = (neut->Vector).Secondary.ParticleID[i];
      runaction -> Ichvc[i]   = (neut->Vector).Secondary.TrackingFlag[i];
      runaction -> Iorgvc[i]  = (neut->Vector).Secondary.ParentID[i];
      runaction -> Iflvc[i]   = (neut->Vector).Secondary.InteractionCode[i];
      runaction -> Abspvc[i]  = (neut->Vector).Secondary.AbsMomentum[i];
      for (int j=0; j<3; j++)
        runaction -> Pvc[i][j]     = (neut->Vector).Secondary.Momentum[i][j];
    }


  // #############################################################################
  // ### Fill primary state info of partcle generated at neutrino interaction
  // #############################################################################
    /*
      NeutInfoSummary* neutinfo = new NeutInfoSummary();
  neutinfo -> Clear();
  neutinfo -> Mode = (neut->Vector).Primary.Mode;
  neutinfo -> Numnu = (neut->Vector).Primary.NumParticle;
  for(int i=0;i<(neutinfo->Numnu);i++) {
    neutinfo -> Ipnu[i] = (neut->Vector).Primary.ParticleID[i];
    neutinfo -> Abspnu[i] = (neut->Vector).Primary.AbsMomentum[i];
    for(int j=0;j<3;j++) neutinfo -> Pnu[i][j] = (neut->Vector).Primary.Momentum[i][j];
  }
  runaction -> GetEvtSum() -> AddNeut( neutinfo );
    */
	for(int ipart=0; ipart<Secondary.NumParticle; ipart++) {
  // #############################################################################
	// ### consider only TrackingFlag for use non interacted particle in neucleus ###
  // #############################################################################
		if( Secondary.TrackingFlag[ipart]==1 ) {

#ifdef DEBUG2
	    G4cout << "Particle:" << (neut->Vector).Secondary.ParticleID[ipart];
	    G4cout << "  Index: " << ipart;
	    G4cout << "  Parent Index: " << (neut->Vector).Secondary.ParentID[ipart] -1 << "\n";
	    G4cout << "Tracking Flag: " << (neut->Vector).Secondary.TrackingFlag[ipart];
	    G4cout << "  Interaction code: " << (neut->Vector).Secondary.InteractionCode[ipart] << "\n";
	    G4cout << " Momentum[MeV/c]:";
	    for (int k=0;k<3;k++)   G4cout << (neut->Vector).Secondary.Momentum[ipart][k]*MeV << " ";
	    G4cout << "\n";
#endif

			G4ParticleDefinition* particle;
			particle = particleTable->FindParticle(Secondary.ParticleID[ipart]);

			double nvec[3];
			for(int ixyz=0; ixyz<3; ixyz++)
				nvec[ixyz] = Secondary.Momentum[ipart][ixyz]/ Secondary.AbsMomentum[ipart];
			G4ThreeVector dir(nvec[0], nvec[1], nvec[2]);

			G4double mass = particle->GetPDGMass();
			G4double mom = Secondary.AbsMomentum[ipart]*MeV;
			G4double energy = sqrt(mass*mass+mom*mom) - mass;

			//
			if(Secondary.ParticleID[ipart]==13)
			  (runaction->p_theta_muons)->Fill(mom/GeV,dir.z());
			else if(Secondary.ParticleID[ipart]==211)
			  (runaction->p_theta_piPos)->Fill(mom/GeV,dir.z());
			else if (Secondary.ParticleID[ipart]==-211)
			  (runaction->p_theta_piNeg)->Fill(mom/GeV,dir.z());
			else if(Secondary.ParticleID[ipart]==2212)
			  (runaction->p_theta_protons)->Fill(mom/GeV,dir.z());

			if(isParticleGun && Secondary.ParticleID[ipart]!=particleGun_pdg){}
			else{
			  particleGun->SetParticleDefinition(particle);
			  particleGun->SetParticleMomentumDirection(dir);
			  particleGun->SetParticleEnergy(energy);
			  particleGun->SetParticleTime(0.0*ns);			
			  particleGun->GeneratePrimaryVertex(anEvent);
			}


#ifdef DEBUG2
	    G4cout << "ipart: " << ipart << "\n";
	    G4cout << "PID:" << (neut->Vector).Secondary.ParticleID[ipart] << "\n";
	    G4cout << "Tracking Flag: " << (neut->Vector).Secondary.TrackingFlag[ipart] << "\n";
	    G4cout << "  Kinetic Energy[MeV]: " << energy << G4endl;;
	    G4cout << "  Momentum:";
	    for (int k=0;k<3;k++)   G4cout << (neut->Vector).Secondary.Momentum[ipart][k]*MeV << " ";
	    G4cout << " [MeV]\n";
#endif

		} // end of if
	} // end of for loop
}

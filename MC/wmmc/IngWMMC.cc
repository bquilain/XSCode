#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
//#include "QGSP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC.hh"
//#include "QGSP_BERT_HP.hh"
//#include "G4NeutronHPData.hh"

#include "INGRID_Dimension.hh"
#include "IngridDetectorConstruction.hh"
//#include "ND280mPhysicsList.hh"
#include "IngridPrimaryGeneratorAction.hh"
#include "IngridRunAction.hh"
#include "IngridEventAction.hh"
#include "IngridTrackingAction.hh"
#include "Neut.hh"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#define PROTON 2
#define HORIZONTAL 3
#define VERTICAL 4 

// For watermodule (measured position of the bars)-------------------
// filling is done in the end of this code
/*double position_xy[VIEWMAX][PLNMAX][CHMAX];
double position_z [VIEWMAX][PLNMAX][CHMAX];
void Initialize_INGRID_Dimension();
*/

// ====================================================================
//     main
// ====================================================================
int main(int argc, char** argv) 
{
	char neutfile[300];
	char output[300];
	char cmd[300];

	int nd = 0;
	int batch = 0;
	int flav = 0;
	// 1:numu, 2:numubar, 3:nue, 4:nuebar
	int seed=-1;//added by koga 2016/7/2
	
	bool produceKinDistributions=false;
 	bool isParticleGun=false;
	int particleGun_pdg=0;

	double CBIRKS_values[3]={0.0185,0.0208,0.0231}; // nominal 0.0208; Birks_Minus = 0.0185; Birks_Plus = 0.0231 
	int BirksIndex=1;

	int NumberOfEvent = 0;

	int c = -1;
	while ((c = getopt(argc, argv, "ho:i:m:b:f:r:dB:p:")) != -1) {
    switch(c){
			case 'o':
				sprintf(output,"%s",optarg);
				break;
			case 'i':
				sprintf(neutfile,"%s",optarg);
				break;
			case 'm':
				nd = atoi(optarg);
				break;
			case 'f':
				flav = atoi(optarg);
				break;
			case 'b':
				batch = 1;
				sprintf(cmd,"%s",optarg);
				break;
			case 'r':
				seed = atoi(optarg);
				break;
                        case 'B':
			        BirksIndex=atoi(optarg);
				if(abs(BirksIndex-1)>1) {
				  cout<<"Please select 0,1 or 2 for Birks constant"<<endl;
				  exit(1);
				}
			        break;
                        case 'd':
			        produceKinDistributions=true;
				break;
                        case 'p':
			       // option -p : propagate only a single type of particle
			       isParticleGun=true;
			       particleGun_pdg=atoi(optarg);
			        break;
			case 'h':
				std::cerr << "o:output root file name" << std::endl;
				std::cerr << "i:input neut file" << std::endl;
				std::cerr << "B:Birks mode: 0(Birks_Minus) 1(nominal) 2(Birks_Plus)" << std::endl;
				std::cerr << "m:2(nd2) or 3(nd3) or 4(nd4)" << std::endl;
				std::cerr << "b:batch command" << std::endl;
				exit(1);
		}
	}

  if( nd==0 ) {
    G4cout << "Select horizontal or vertical or proton module" << G4endl;
    exit(1);
  }
  if( flav==0 ) {
    G4cout << "Select neutrino flavor" << G4endl;
    exit(1);
  }
  
  cout<<"Birks constant is set to "<<CBIRKS_values[BirksIndex]<<endl;
  Initialize_INGRID_Dimension();

  // run manager
  G4RunManager* runManager= new G4RunManager;  //G4cout << G4endl;

  // Neut initialization
  Neut *neut = new Neut;
  if( flav==-1 ){
	NumberOfEvent = 5000;
  	G4cout << "NumberOfEvent :" << NumberOfEvent << G4endl;
  }
  else{
	NumberOfEvent = neut->NtupleReadInit(neutfile);
	G4cout << "NumberOfEvent :" << NumberOfEvent << G4endl;
  }
  // set mandatory user initialization classes...

  // detector setup
  runManager-> SetUserInitialization(new IngridDetectorConstruction(nd));
  //runManager-> SetUserInitialization(new IngridDetectorConstruction());
  G4cout << "Detector Init OK" << G4endl;

  // particles and physics processes
  //runManager-> SetUserInitialization(new ND280mPhysicsList);
  //runManager-> SetUserInitialization(new QGSP);
  runManager-> SetUserInitialization(new QGSP_BERT); // ML 2017/05/02
  //runManager-> SetUserInitialization(new QGSP_BIC);
  G4cout << "PhysicsList Init OK" << G4endl;

	IngridRunAction * rac = new IngridRunAction(output);
	rac->produceDistributions= produceKinDistributions;
	G4cout << "RunAction init OK" << (produceKinDistributions? " ... producing kinematic distributions":"")<< G4endl;

	IngridEventAction * evt = new IngridEventAction(rac);
	G4cout << "EventAction init OK" << G4endl;

	IngridTrackingAction * tra = new IngridTrackingAction(rac, evt);
	runManager->SetUserAction(tra);
	G4cout << "TrackingAction init OK" << G4endl;

	IngridPrimaryGeneratorAction *gen=new IngridPrimaryGeneratorAction(neut, rac, evt, nd, flav,seed);
	gen->SetParticleGun(isParticleGun,particleGun_pdg);
	runManager-> SetUserAction(gen);
	G4cout << "PrimaryGenerator init OK" ;
	if(isParticleGun) G4cout<<" ... propagating only particles of type " << particleGun_pdg << G4endl;
	else G4cout<<G4endl;
    
  // user action classes... (optional)
	runManager-> SetUserAction(rac);
	runManager-> SetUserAction(evt);

#ifdef G4VIS_USE
  // initialize visualization package
  G4VisManager* visManager= new G4VisExecutive;
  visManager-> Initialize();
  G4cout <<"visualization init OK" << G4endl;
#endif

  // Initialize G4 kernel
  runManager-> Initialize();
    
	// get the pointer to the UI manager and set verbosities
	G4UImanager* UI= G4UImanager::GetUIpointer();

	if(batch==1)
	// Define (G)UI terminal for interactive mode  
	{ 
			// G4UIterminal is a (dumb) terminal.
			G4UIsession * session = 0;
			G4String command = "/control/execute ";
			G4String macro = cmd;
			session = new G4UIterminal(new G4UItcsh);
			UI->ApplyCommand(command+macro);
			session->SessionStart();
			delete session;
	}

	else { // batch mode
	  			runManager->BeamOn(NumberOfEvent);
	}

  // terminating...
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;  //G4cout << G4endl;

  return 0;

}

//------------------------------------------------
/*
void Initialize_INGRID_Dimension(){
  TFile* f = TFile::Open("/export/scbn25/data1/taichiro/T2K/work/basesoft/git/watermodule/wmmc/assembly_checklist.root","READ");
  TTree* t = (TTree*) f->Get("tree");
  double xy1,xy3,xy0;
  int view,pln,ch;
  t -> SetBranchAddress("xy1",&xy1);
  t -> SetBranchAddress("xy3",&xy3);
  t -> SetBranchAddress("view",&view);
  t -> SetBranchAddress("pln",&pln);
  t -> SetBranchAddress("ch",&ch);
  for(int i=0;i<t->GetEntries();i++){
    t -> GetEntry(i);
    if(ch==0){
        xy0 = (xy1+xy3)/2./10. + 1.;//offset
    }
    if(ch<40){
      position_xy[view][pln][ch] = (xy1+xy3)/2./10.; //mm->cm
      position_xy[view][pln][ch] = position_xy[view][pln][ch] - xy0; //mm->cm
    }
    else{
      position_xy[view][pln][ch] = 0;
    }
  }
  t->Delete();
  f->Close();
  f->Delete();

  f = TFile::Open("/export/scbn25/data1/taichiro/T2K/work/basesoft/git/watermodule/wmmc/position_module_z.root","READ");
  t = (TTree*) f->Get("tree");
  double z;
  double z0;
  t -> SetBranchAddress("z",&z);
  t -> SetBranchAddress("view",&view);
  t -> SetBranchAddress("pln",&pln);
  for(int i=0;i<t->GetEntries();i++){
    t -> GetEntry(i);
    if(i==0){
        z0 = z/10.;//offset
    }
    for(int ch_temp=0;ch_temp<80;ch_temp++){
    	position_z [view][pln][ch_temp] = z/10.; //mm->cm
    	position_z [view][pln][ch_temp] = position_z [view][pln][ch_temp] - z0;
    }
  }
  t->Delete();
  f->Close();
  f->Delete();
  //  std::cout << "call DIMENSION" << std::endl;
}
*/

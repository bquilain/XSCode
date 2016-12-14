#include "IngridResponse.hh"
#include <CLHEP/Random/Randomize.h>
#include <algorithm>
#include <TRandom.h>

using namespace std;

const G4double scilen = 120.;  //cm
const G4double longlen = 130.;  //cm
const G4double longlen_pm = 120.;  //cm
const G4double shortlen = 112.;  //cm
const G4double attleng = 241.7; //cm
const G4double sciattleng = 10.46; //cm added
//const G4double SciBarFactor = 1.96;  //P.E. factor for SciBar scintillator
//const G4double SciBarFactor = 1.78;  //P.E. factor for SciBar scintillator
//const G4double SciBarFactor = 1.72;  //P.E. factor for SciBar scintillator
//Default: const G4double SciBarFactor = 1.77;  //P.E. factor for SciBar scintillator
const G4double SciBarFactor = 1.88;//1.95
const G4double CBIRKS = 0.0231;//0.0208; // used in SciBooNE MC
const G4double CBIRKS_LE = 10*CBIRKS;
//const G4double CBIRKS_Hadrons = 0.01;
//const G4double CBIRKS = 0.1;//Benjamin Test (other solution: dans les stepping action: if(track->getdef==Gpartneutron) trackstatus(stopandkill)
//const G4double CBIRKS_LE = 1;
const G4double TransTimeInFiber = 1./28.;  //  1cm/2.8e10[cm/s] * 1e9 [ns]
const G4double Pedestal = 0;//145;  // pedeltal of ADC counts
//const G4double Gain = 10;  // Gain ADC counts of high gain channel
//const G4double LowGain  = 1;  // Gain ADC counts of low gain channel
const G4double ElecNoise = 1.7; // sigma of high gain electronics noise
const G4double LowElecNoise = 1.2;  // sigma of low gain electronics noise

/****B.Quilain********/
const G4double Gain = 10.35;
const G4double LowGain = 1.09;
//const G4double ElecNoise = 0.34;
//const G4double LowElecNoise = 0.09;


const G4double PixelGainVari = 0.031;  // gain variation among pixels

//
//const G4double MeV2PE = 15./2.6;  // pe/MeV
//const G4double MeV2PE = 18.45/2.6;  // pe/MeV
//const G4double MeV2PE = 36.6;  // v1r1_3
//const G4double MeV2PE = 35.1;  // v2
//const G4double MeV2PE = 37.2;  // v2_1
//const G4double MeV2PE = 39.1;  // v2r2
//const G4double MeV2PE = 60.;  // v3
//const G4double MeV2PE = 60.6;  // v3r1
//const G4double MeV2PE = 49.47;  // v3r2
//const G4double MeV2PE = 40.130;  // v3r3
//const G4double MeV2PE = 45.1;  // v3r4
//const G4double MeV2PE = 48.8;  // v3r4
//const G4double MeV2PE = 46.1;  // v3r4
//Kikawa-san one's: const G4double MeV2PE = 44.8;  // v3r4
/***B.Quilain 01/05/2013******/
//const G4double MeV2PE = 36.9;
//const G4double MeV2PE = 35.2;
//const G4double MeV2PE = 46.0;
const G4double MeV2PE_PM = 40.0;
const G4double MeV2PE = 46.0;
//const G4double MeV2PE = 41.0;

const G4double MPPCPixel = 667.; //v3
//const G4double Eff_PDE = -0.292;  // @deltaV = 1.14 V
//const G4double Eff_PDE = -0.38;  // @deltaV = 1.4 V
const G4double Eff_PDE = -0.275;//*1.134;  // @deltaV = 1.1 V

//
const G4double ETC_CONST = 5.0;  // v3
const G4double rPDE = 1.7;  // @deltaV = 1.1 V
const G4double cross_after_rate = 0.09;//*1.051;  // @deltaV = 1.1 V



//
const G4int topview = 1;
const G4int sideview = 0;

//
IngridResponse::IngridResponse()
{
}

//
IngridResponse::~IngridResponse()
{
}

double IngridResponse::ApplyConversion(G4double* edep)
{
  return ((*edep)*MeV2PE);
}

//
void IngridResponse::ApplyScintiResponse(G4double* edep, G4Track* aTrack)
{
  //quenching
  BirksSaturation(edep,aTrack);
  
  return;
}

void IngridResponse::ApplyScintiResponse2(G4double* edep, G4double* length,G4Track* aTrack)
{
  //quenching                                                                                                                             
  BirksSaturation2(edep,length,aTrack);

  return;
}


//
void IngridResponse::BirksSaturation(G4double* edep, G4Track* aTrack)
{
  //const G4double CBIRKS = sbcard->Birks;
  
  G4double              kineticE = aTrack->GetKineticEnergy();
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4Material*           material = aTrack->GetMaterial();

/*
  G4cout << "=======" << G4endl;
  G4cout << "edep : " << (*edep) << G4endl;
  G4cout << "kineticE : " << kineticE << G4endl;
  G4cout << "particle : " << particle->GetParticleName() << G4endl;
  G4cout << "material : " << material->GetName() << G4endl;
*/

  if(particle->GetPDGCharge()==0) return;

  G4double dedx = emcal.GetDEDX(kineticE, particle, material)/(MeV/cm);

  /*
  G4cout << "dedx : " << dedx << G4endl;
  G4cout << G4endl;
  */
  //G4cout<<"Lepton Number ==================================================================================================================================================== "<<particle->GetLeptonNumber()<<" PDG code="<<particle->GetPDGEncoding()<<G4endl;
  //if(particle->GetLeptonNumber()!=0) (*edep) = (*edep)/(1. + CBIRKS*dedx);
  //else (*edep) = (*edep)/(1. + CBIRKS_Hadrons*dedx);
  /*  if(kineticE==0) (*edep) = (*edep)/(1. + CBIRKS_LE*dedx);
  else (*edep) = (*edep)/(1. + CBIRKS*dedx);
*/
  //G4cout<<"**********************************************************************************************************************************************************************************"<<G4endl<<"edep="<<(*edep)<<", dedx="<<dedx<<", track NRJ="<<kineticE/MeV<<G4endl;
  if(kineticE/MeV < 5) (*edep) = (*edep)/(1. + 1000*CBIRKS*dedx);
  else (*edep) = (*edep)/(1. + CBIRKS*dedx);
  //G4cout<<"edep after attenuation="<<(*edep)<<G4endl;
  /*
  G4double mmu = 105.65837;
  for (int i=0;i<100000;i++) {
  G4double p = mmu * (i+1)*0.1;
  G4double kineticE = (sqrt(mmu*mmu+p*p)-mmu)*MeV;
  G4double dedx = emcal.GetDEDX(kineticE, particle, material)/(MeV/cm)/1.032;
  
  G4cout << (i+1)*0.1 << " " << dedx << G4endl;
  }
  */
  return;
}



void IngridResponse::BirksSaturation2(G4double* edep, G4double* length, G4Track* aTrack)
{
  //const G4double CBIRKS = sbcard->Birks;
  
  G4double              kineticE = aTrack->GetKineticEnergy();
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4Material*           material = aTrack->GetMaterial();

  if(particle->GetPDGCharge()==0) return;


  G4double dedx = (*edep)/(*length);//emcal.GetDEDX(kineticE, particle, material)/(MeV/cm);

  /*
  G4cout << "dedx : " << dedx << G4endl;
  G4cout << G4endl;
  */
  //G4cout<<"Lepton Number ==================================================================================================================================================== "<<particle->GetLeptonNumber()<<" PDG code="<<particle->GetPDGEncoding()<<G4endl;
  //if(particle->GetLeptonNumber()!=0) (*edep) = (*edep)/(1. + CBIRKS*dedx);
  //else (*edep) = (*edep)/(1. + CBIRKS_Hadrons*dedx);
  /*  if(kineticE==0) (*edep) = (*edep)/(1. + CBIRKS_LE*dedx);
  else (*edep) = (*edep)/(1. + CBIRKS*dedx);
*/
  //if((*edep)*(MeV2PE)<20){
  //if(kineticE==0){
  /*  if(kineticE==0){
    G4cout << "particle : " << particle->GetParticleName() << G4endl;                                                                   
    G4cout << "material : " << material->GetName() << G4endl;
    G4cout<<"*******************************************************************************************************************************************************************************"<<G4endl<<"edep="<<(*edep)<<", dedx="<<dedx<<", de/dx from emcal="<<emcal.GetDEDX(kineticE, particle, material)/(MeV/cm)<<", track NRJ="<<kineticE/MeV<<", step length="<<*length<<G4endl;
    }*/
  //if(kineticE<1){G4cout<<"**********************************************************************************************************************************************************************************"<<G4endl<<"edep="<<(*edep)<<", dedx="<<dedx<<", track NRJ="<<kineticE/MeV<<", step length="<<*length<<G4endl;}

  (*edep) = (*edep)/(1. + CBIRKS*dedx);
  /*
  if((kineticE/MeV) !=0) (*edep) = (*edep)/(1. + CBIRKS*dedx);
  else (*edep) = (*edep)/(1. + 2*CBIRKS*dedx);
  
  if(dedx!=0) (*edep) = (*edep)/(1. + CBIRKS*dedx);
  else (*edep)=0;
  */
  /*  if(particle->GetPDGEncoding()>10000) (*edep)=0;
  else *(edep) = (*edep)/(1. + CBIRKS*dedx);
  */
  //if((*edep)*(MeV2PE)<20)
  //if(kineticE==0) G4cout<<"edep after attenuation="<<(*edep)<<G4endl;
  //if(kineticE==0)
  return;
}


void IngridResponse::ApplyLightCollection(G4double* edep, G4int mod, G4int view, G4ThreeVector pos){
  
  G4double x = 0.;
  G4int i = 0.;

  if( view==topview ) x = fabs(scilen/2. + pos[0]/cm);
  else if( view==sideview ) x = fabs(scilen/2. + pos[1]/cm);
  
  if(x<40||x>80||mod!=16){
    i=x/5;
    x=fabs(x-i*5-2.5);
    *edep *= exp(-1.*x/sciattleng);
  }
  else{
    i=x/2.5;
    x=fabs(x-i*2.5-1.25);
    *edep *= exp(-1.*x/sciattleng)*SciBarFactor;
  }


  return;
}



void IngridResponse::ApplyFiberResponse(G4double* edep, G4double* time, G4int view, G4ThreeVector pos)
{
  G4double x = 0.;
  
  if( view==topview ) x = fabs(scilen/2. - pos[1]/cm);
  else if( view==sideview ) x = fabs(scilen/2. + pos[0]/cm);
  
  /*Added by B.Quilain*/
  //G4double x2=0;
  //if( view==topview ) x2=fabs(scilen/2+pos[1]/cm);
  //else if( view==sideview ) x2=fabs(scilen/2-pos[0]/cm);
  //*edep = (*edep*exp(-1.*x/attleng)+*edep*exp(-1.*(2*x2+x)/attleng))*.5;
  // attenuation
  *edep *= exp(-1.*x/attleng);
  
  // delay in fiber
  *time += TransTimeInFiber*x;

  return;
}

//
void IngridResponse::ApplyFiberResponseV(G4double* edep, G4double* time, G4int pln, G4ThreeVector pos)
{
  G4double x = 0.;

  /*
    if( pln==11 ) x = fabs(longlen/2. - (pos[1]/cm + 3.7) );
    if( pln==12 ) x = fabs(longlen/2. + (pos[1]/cm - 0.3) );
    if( pln==13 ) x = fabs(shortlen/2. - (pos[0]/cm + 0.9) );
    if( pln==14 ) x = fabs(longlen/2. + (pos[0]/cm - 5.9) );
  */
  

  
  // attenuation
  *edep *= exp(-1.*x/attleng);
  
  // delay in fiber
  *time += TransTimeInFiber*x;
  
  return;
}

//
void IngridResponse::ApplyMPPCResponse(G4double edep, G4int mod, G4double* pe)
{
  G4double npe;
  
  // energy to p.e.
  //pe = (edep/MeV)*MeV2PE;
  if(mod==16) npe = edep*(MeV2PE_PM);
  else npe = edep*(MeV2PE);  
  // PDE 
  /*
    npe = edep * rPDE * ETC_CONST;
  */

  // MPPC linearity
  npe = MPPCPixel * (1. - exp( Eff_PDE * npe / MPPCPixel ));
  //*pe = npe;
  /*
    npe = MPPCPixel * (1. - exp( -1.* npe / MPPCPixel ));
  */
  
  // fake signal p.e. from cross-talk & after pulse
  npe = npe / (1. - cross_after_rate);


  // Poisson statistics & 1 pe resolution
  npe = CLHEP::RandPoisson::shoot(npe);

  npe = gRandom->Gaus(npe,npe*PixelGainVari);
  
  //
  *pe = npe;

  return;
}


//
void IngridResponse::ApplyADCResponse(G4double *pe, G4double *lope, G4int* adc, G4int* loadc)
{

  G4double adc_tmp, loadc_tmp,Q,loQ;
  
  //PE to ADC
  adc_tmp = Pedestal + (*pe)*Gain;
  loadc_tmp = Pedestal + (*pe)*LowGain;
  
  //Electronics noise
  adc_tmp = gRandom->Gaus(adc_tmp,ElecNoise);
  loadc_tmp = gRandom->Gaus(loadc_tmp,LowElecNoise);

  /*added by b.quilain*/
  //  adc_tmp2=gRandom->Gaus(adc_tmp,ChannelDiff);
  //loadc_tmp2 = gRandom->Gaus(loadc_tmp,ChannelDiff);
  //ADC to Charge
  //Q=(adc_tmp+53)/217;
  //loQ=(loadc_tmp+82)/26;
  Q=(adc_tmp)/135.5;
  loQ=(loadc_tmp)/14.29;

  //Non linearlity of high gain ADC
    if(Q<0.65)*adc=135.5*Q;
  else if(Q<3.2)*adc=217*Q-53;
  else if(Q<4.2)*adc=158.6*Q+133.9;
  else if(Q<14)*adc=5.1*Q+778.6;
  else *adc=850;

  //Non linearlity of low gain ADC
  if(loQ<7)*loadc=14.29*loQ;
  else if(loQ<27)*loadc=26*loQ-82;
  else if(loQ<35.5)*loadc=21.12*loQ+48.24;
  else if(loQ<178.4)*loadc=0.7*loQ+775.1;
  else *loadc=900;
  

  //ADC to PE
  //*pe = (float)((*adc) - Pedestal)/Gain;
  //*lope = (float)((*loadc) - Pedestal)/LowGain;
  //changed by B.Quilain to remove non linearity added by kikawa-san in MC:
  *pe = (float)((adc_tmp) - Pedestal)/Gain;
  *lope = (float)((loadc_tmp) - Pedestal)/LowGain;


/*
  const G4double VA_NOISE = 1.6;
  const G4double VA_OVERFLOW = 2048.0;

  G4double adc_tmp;

  IngridMap* sbmap= sbdb->GetSciBarMap();

  for (G4int view=0;view<2;view++) {
    for (G4int layer=0;layer<64;layer++) {
      for (G4int strip=2;strip<114;strip++) {

	G4int isci = sbmap->sbisci(view,layer,strip);

	// pe -> ADC conversion
	SciBarCalibration* sbcal= sbdb->GetSciBarCalibrationInfo();
	adc_tmp = pe[isci]/sbcal->pe[isci];

	// VA noise
	adc_tmp += CLHEP::RandGaussQ::shoot(0.0,VA_NOISE);

	// overflow
	if (adc_tmp>VA_OVERFLOW) adc_tmp = VA_OVERFLOW;

	adc[isci] = (int)adc_tmp;
      }
    }
  }
*/
}

//
void IngridResponse::ApplyTDCResponse(G4double time, G4int* tdc)
{
  *tdc = 0;
}

#ifndef Neut_h
#define Neut_h 1

#define NVEC 100
#define NG 20

#include <stdio.h>
#include <iostream>

#define JNUBEAM10C

//
typedef struct {
  int NumParticle, ParticleID[NVEC];
  int ParentID[NVEC], TrackingFlag[NVEC], InteractionCode[NVEC];
  float AbsMomentum[NVEC], Momentum[NVEC][3], Vertex[3];
} SecondaryVector;

//
typedef struct {
    int gipart;
    float gpos0[3],gvec0[3],gamom0;
} PrimaryParticle;

//
typedef struct {
    int spid;
    int pgen;
    float psi0,xsi0[3],nsi0[3],cossi0bm;
    float xsi[3];
    int smech;
    int intgt;
    float prvtx[3];
    int smed;
    int gppid;
    float xgpi0[3],xgpi[3],pgpi0;
    int gpmech, gpmed;
    int prmech, prmed;
    int prdght, sdght, gpdght;
    int chain;
} OptionalVariables;

//
typedef struct {
	int ng;
	float gpx[NG],gpy[NG],gpz[NG],gcosbm[NG];
	float gvx[NG],gvy[NG],gvz[NG];
	int gpid[NG],gmec[NG];
} Ancestor;

//
typedef struct {
	float Enusk;
	float normsk;
	float anorm;
} SKInfo;

//
typedef struct {
	int Numbndn;
	int Numbndp;
	int Numfrep;
	int Numatom;
} TargetInfo;

//
typedef struct {
	float Crsenergy;
	float Totcrsne;					// total cross-section
	float Difcrsne[8];
} NeutCrs;

//
class Neut
{
  int NumberOfEvent;

    typedef class ParentVector{
      float AbsMomentum, Vertex[3], Direction[3], CosBeam;
    } ;

public:
  struct{

    int EventNumber;

    struct{
      int Mode;
      int NumParticle, ParticleID[NVEC];
      float AbsMomentum[NVEC], Momentum[NVEC][3];
    } Primary;
    
    SecondaryVector Secondary;
        
    struct{
      float Energy;
      int ParentID, ProductionMode;
      ParentVector ParentDecay;
      float Norm;
      int nvtx0;
      ParentVector ParentProduction;
      float r, x, y, Direction[3];
      int FDID;
      PrimaryParticle primary_particle;
#ifdef JNBEAM10A
      OptionalVariables optional_variables;
#endif
#ifdef JNUBEAM10C
			Ancestor ancestor;
      SKInfo sk_info;
#endif
    } Neutrino; 

#ifdef JNUBEAM10C
    TargetInfo target_info;
		NeutCrs neutcrs;
#endif

  } Vector;

	int ID;

public:

  int NtupleReadInit(char*);
  int NtupleReadEvent(SecondaryVector&,float*);

};

#endif

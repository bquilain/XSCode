#ifndef __INGTREE_HH__
#define __INGTREE_HH__ 1

#include "TObject.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TROOT.h"

#include <vector>
#include <iostream>


#define MAX_HITS 1000
#define NMOD 16


//////////////////////////////////
class IngTrack : public TObject {

private:
    Int_t PID;			//particle code according to PDG
    Int_t ID;			//Track ID

    Int_t NHitX;		//# of Hits in xlayer-scinti
    Int_t NHitY;		//# of Hits in ylayer-scinti
    Int_t NHitXVeto;		//# of Hits in xlayer-veto-scinti
    Int_t NHitYVeto;		//# of Hits in ylayer-veto-scinti
    
    std::vector<Int_t> HitXID;		//array of Hit index in x-scinti
    std::vector<Int_t> HitYID;		//array of Hit index in y-scinti
    std::vector<Int_t> HitXVetoID;	//array of Hit index in v-scinti
    std::vector<Int_t> HitYVetoID;	//array of Hit index in v-scinti
    
    Float_t StartPos[3];	//Starting point of this track [x, y, z]
    Float_t StopPos[3];		//Stopping point of this track
    Float_t Pi[3];		//initial momentum
    Float_t Pf[3];		//final momentum

public:
    IngTrack() : TObject(),NHitX(0),NHitY(0),NHitXVeto(0),NHitYVeto(0) {};
    IngTrack(IngTrack &trk) : TObject(trk),NHitX(0),NHitY(0),NHitXVeto(0),NHitYVeto(0)
    {
	PID = trk.GetPID();
	ID = trk.GetID();
	for(int i=0;i<3;i++) {
	    StartPos[i] = trk.GetStartPos(i);
	    StopPos[i] = trk.GetStopPos(i);
	    Pi[i] = trk.GetFirstMom(i);
	    Pf[i] = trk.GetFinalMom(i);
	}
    };

   virtual ~IngTrack() { Clear();};
   void Clear(Option_t *option="");

    Int_t GetPID() const { return PID; }
    Int_t GetID() const {return ID; }
    Int_t GetNHitX() const  { return NHitX; }
    Int_t GetNHitY() const  { return NHitY; }
    Int_t GetNHitXVeto() const  { return NHitXVeto; }
    Int_t GetNHitYVeto() const  { return NHitYVeto; }

    Int_t GetHitID(int index,int view,int pln);

    Float_t GetStartPos(Int_t i=0) const { return (i<3)?StartPos[i]:0; }
    Float_t GetStopPos(Int_t i=0)  const { return (i<3)?StopPos[i]:0; }
    Float_t GetFirstMom(Int_t i=0) const { return (i<3)?Pi[i]:0; }
    Float_t GetFinalMom(Int_t i=0) const { return (i<3)?Pf[i]:0; }
    
    void SetPID(Int_t pid) { PID=pid; }
    void SetID(Int_t id) { ID=id; }
    void AddHitID(int id,int view,int pln);

    void SetStartPos(float *pos) { for(int i=0;i<3;i++) StartPos[i]=pos[i]; }
    void SetStopPos(float *pos) { for(int i=0;i<3;i++) StopPos[i]=pos[i]; }
    void SetFirstMom(float *mom) { for(int i=0;i<3;i++) Pi[i]=mom[i]; }
    void SetFinalMom(float *mom) { for(int i=0;i<3;i++) Pf[i]=mom[i]; }

   ClassDef(IngTrack,2) 
};

///////////////////////////////////////////
class TrueTrack : public TObject {

private:
   Int_t          NTrack;            //Number of tracks
   TClonesArray  *fTracks;            //->array with all tracks

    Long_t fSeed[10];
    Int_t fN;

   static TClonesArray *fgTracks;

public:
   TrueTrack();
   virtual ~TrueTrack(){ Clear();};
   void Clear(Option_t *option ="");
   void Reset(Option_t *option ="");
    void SetSeed(long seed){ fSeed[fN]=seed; fN++; }

   IngTrack* AddTrack(Int_t id, Int_t pid, Float_t *pos1, Float_t *pos2, Float_t *mom1, Float_t *mom2);
   IngTrack* AddTrack(IngTrack &track);

   Int_t GetNTrack() const { return NTrack; }
   TClonesArray *GetTracks() const { return fTracks; }
    IngTrack *GetIngTrack(Int_t id);

    Int_t CheckTrack(Int_t id);
    Int_t SetHitID(int id, int hitid, int view, int pln); 

    Long_t GetSeed(int i) const { return fSeed[i]; }

   ClassDef(TrueTrack,1)

};


///////////////////////////////////
class IngHit : public TObject {

private:

    Int_t Mod;
    Int_t View;    	//Top View(X layer):1 Side View(Y layer):0
    Int_t Pln;		//0~10:TPL, 11~14:Veto
    Int_t Ch;    	//0~23:TPL, 0~21:Veto

    Int_t HighAdc;
    Int_t LowAdc;
    Long_t Time;
    Int_t HighPe;
    Int_t LowPe;

    /// Special information of MC 
    Int_t ID;			// Hit ID
    Float_t Edep;		//energy deposite
    Float_t Position[3];		//Position
    Int_t PID;

public:
    IngHit() {};

    /// constructor for Raw Data
    IngHit(int mod,int view,int pln,int ch,int hadc,int ladc,long time) 
    : Mod(mod),View(view),Pln(pln),Ch(ch),HighAdc(hadc),LowAdc(ladc),Time(time),
    HighPe(0),LowPe(0),ID(0),Edep(0.)
    { 
	for(int i=0;i<3;i++) Position[i]=0.; 
    };
	
    /// constructor for MC
    IngHit(int mod,int view,int pln,int ch,int id,float *pos,float edep,int pid) 
    : Mod(mod),View(view),Pln(pln),Ch(ch),
    HighAdc(0),LowAdc(0),Time(0),HighPe(0),LowPe(0),ID(id),Edep(edep),PID(pid)
    { 
	for(int i=0;i<3;i++) Position[i]=pos[i]; 
    }; 

    virtual ~IngHit() {};

    void Clear(Option_t *option="");

    Int_t	GetModID() const { return Mod; }
    Int_t 	GetView() const { return View; }
    Int_t 	GetPln() const { return Pln; }
    Int_t	GetCh() const { return Ch; }

    Float_t	GetPos(Int_t i=0)  const { return (i<3)?Position[i]:0; }
    Float_t 	GetEdep() const { return Edep; }
    Int_t 	GetID() const { return ID; }
    Int_t 	GetPID() const { return PID; }
    
    void Init(int mod,int view,int pln,int ch,int hadc,int ladc,long time)
    {
	Mod=mod; View=view; Pln=pln; Ch=ch;
        HighAdc=hadc; LowAdc=ladc; Time=time;
    };

  void SetHighPe(int hpe)       { HighPe=hpe; }
  void SetLowPe(int lpe)        { LowPe=lpe; }

  void SetHighAdc(int hadc)     { HighAdc=hadc; }
  void SetLowAdc(int ladc)      { LowAdc=ladc; }
  void SetTime(long time)       { Time=time; }

  void SetMod(int mod)          { Mod=mod; }
  void SetView(int view)        { View=view; }
  void SetPln(int pln)          { Pln=pln; }
  void SetCh(int ch)            { Ch=ch; }

    void SetEdep(double edep)   { Edep=edep; }
    void SetPos(float *pos) { for(int i=0;i<3;i++) Position[i]=pos[i]; }
    void SetID(int id) { ID=id; }

   ClassDef(IngHit,3)  //A hit segment

};


///////////////////////////////////////////
class ModHit : public TObject {

private:
   Int_t	Mod; 
   Int_t        NHitX;            //Number of xlayer-scinti hits
   Int_t        NHitY;            //Number of ylayer-scinti hits
   Int_t        NHitXVeto;            //Number of veto-xlayer-scinti hits
   Int_t        NHitYVeto;            //Number of veto-ylayer-scinti hits

   TClonesArray  *fXHits;            //->array with all xlayer-scinti hits
   TClonesArray  *fYHits;            //->array with all ylayer-scinti hits
   TClonesArray  *fXVetoHits;            //->array with all veto-xlayer-scinti hits
   TClonesArray  *fYVetoHits;            //->array with all veto-ylayer-scinti hits

   static TClonesArray *fgXHits;
   static TClonesArray *fgYHits;
   static TClonesArray *fgXVetoHits;
   static TClonesArray *fgYVetoHits;

public:
   ModHit();
   ~ModHit();

    void Init();
   void   Clear(Option_t *option ="");
   void   Reset(Option_t *option ="");
    
   void SetMod(int mod) { Mod=mod; } 

    IngHit* AddHit(int mod,int view,int pln,int ch,int hadc,int ladc,long time);
    IngHit* AddHit(int mod,int view,int pln,int ch,int id,float *pos,float edep,int pid);

   Int_t GetMod() const { return Mod; }
   Int_t GetNHitX() const { return NHitX; }
   Int_t GetNHitY() const { return NHitY; }
   Int_t GetNHitXVeto() const { return NHitXVeto; }
   Int_t GetNHitYVeto() const { return NHitYVeto; }

    Int_t AddNHitX() { return NHitX++; }
    Int_t AddNHitY() { return NHitY++; }
    Int_t AddNHitXVeto() { return NHitXVeto++; }
    Int_t AddNHitYVeto() { return NHitYVeto++; }

    TClonesArray *GetXHits() const { return fXHits; }
    TClonesArray *GetYHits() const { return fYHits; }
    TClonesArray *GetXVetoHits() const { return fXVetoHits; }
    TClonesArray *GetYVetoHits() const { return fYVetoHits; }
    
    IngHit *GetXHit(Int_t i);
    IngHit *GetYHit(Int_t i);
    IngHit *GetXVetoHit(Int_t i);
    IngHit *GetYVetoHit(Int_t i);

    ClassDef(ModHit,2)  //

};



class EvtSum : public TObject {

private:
  Int_t          Event;                //Event number
  Int_t          Spill;                //Beam Spill number
  Int_t          Cycle;                //Cycle number(0~22)
  Int_t          TrgId;                //Trigger mode
                                       //1:Beam
                                       //2:Periodic
                                       //128:Cosmic
  Long_t         TrgTime;              //Trigger Time[nsec]
  Long_t         UTime;                //Header Time[sec]

  ModHit	fModHit[16];

public:
    
    EvtSum();

    virtual ~EvtSum() 
    { 
    };

    void Clear(Option_t *option="");
    void Reset(Option_t *option="");

    void SetEvent(int event) { Event=event; }
    void SetSpill(int spill) { Spill=spill; }
    void SetCycle(int cycle) { Cycle=cycle; }
    void SetTrgId(int trgid) { TrgId=trgid; }
    void SetTrgTime(long trgtime) { TrgTime=trgtime; }
    void SetUTime(long utime) { UTime=utime; }

    void AddHit(int mod,int view,int pln,int ch,int hadc,int ladc,long time);
    void AddHit(int mod,int view,int pln,int ch,int id,float *pos,float edep,int pid);

    Int_t GetEvent() const { return Event; }
    Int_t GetSpill() const { return Spill; }
    Int_t GetCycle() const { return Cycle; }
    Int_t GetTrgId() const { return TrgId; }
    Long_t GetTrgTime() const { return TrgTime; }
    Long_t GetUTime() const { return UTime; }

    ModHit* GetModHit(int index) { return &fModHit[index]; }

    ClassDef(EvtSum,1)

};

#if 0

/////////////////////////////////////////////////
class TNeut 
{
private:
  int ID;  	// Module ID
  float Enu;
  int NeutMode;
  float NeutDirection[3];
  float AbsMom;
  float Mom[3];
  //
    float Proton_AbsMom;
    float Proton_Mom[3];

  float Vertex[3];
  
public:
  TNeut();
  ~TNeut();
  void Clear();
  void SetID(int id){ID = id;}  
  void SetEnu(float enu){Enu = enu;}
  void SetNeutMode(int neutmode){NeutMode = neutmode;}
  void SetNeutDirection(float* neutdirection);
  void SetAbsMom(float absmom){AbsMom = absmom;}
  void SetMom(float* mom);
  // QEL
    void SetProton_AbsMom(float absmon) { Proton_AbsMom=absmon; }
    void SetProton_Mom(float *mom);
  void SetVertex(float* vertex);
  //
    //
  ////
  ////
  int GetID(){return ID;}
  float GetEnu(){return Enu;}
  int GetNeutMode(){return NeutMode;}
  float* GetNeutDirection(){return NeutDirection;}
  float GetAbsMom(){return AbsMom;}
  float* GetMom(){return Mom;}
  //
    float GetProton_AbsMom() { return Proton_AbsMom; }
    float* GetProton_Mom() { return Proton_Mom; }
  float* GetVertex(){return Vertex;}

    ClassDef(TNeut,1)
};


#endif

#define NVECT 100

/*
class ParentInfo
{
public:
    ParentInfo() {};
    ~ParentInfo() {};

    float AbsMomentum, Vertex[3], Direction[3], CosBeam;

    ClassDef(ParentInfo,1)
};
*/

class  PrimaryState
{
public:
    PrimaryState() {};
    ~PrimaryState() {};

    Int_t Mode;
    Int_t NumParticle, ParticleID[NVECT];
    Float_t AbsMomentum[NVECT], Momentum[NVECT][3];

    ClassDef(PrimaryState,1)
};

class FinalState
{
public:
    FinalState() {};
    ~FinalState() {};

    Int_t NumParticle, ParticleID[NVECT];
    Int_t ParentID[NVECT], TrackingFlag[NVECT], InteractionCode[NVECT];
    Float_t AbsMomentum[NVECT], Momentum[NVECT][3], Vertex[3];

    ClassDef(FinalState,1)
};

class NeutInfo
{
public:
    NeutInfo() {};
    ~NeutInfo() {};

    Float_t Energy;
    Int_t ParentID, ProductionMode;
//    ParentInfo ParentDecay;
    Float_t WeightingFactor;
    Int_t nvtx0;
//    ParentInfo ParentProduction;
    Float_t r, x, y, Direction[3];
    Int_t FDID;
    Int_t ModID;

    ClassDef(NeutInfo,1)
};

class TNeut
{
public:
    TNeut() {};
    ~TNeut() {};

    Int_t NumberOfEvent;

    PrimaryState primary;

    FinalState final;

    NeutInfo neutinfo;

    Float_t vertex[3];

    ClassDef(TNeut,1)
};

#endif

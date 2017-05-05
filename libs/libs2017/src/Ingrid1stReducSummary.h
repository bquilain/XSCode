#ifndef INGRID1STREDUCSUMMARY_H
#define INGRID1STREDUCSUMMARY_H
#include <iostream>
#include "TObject.h"
#include "TClonesArray.h"

#include "TRefArray.h"
#include "TRef.h"
#include "vector"
using namespace std;

#include "IngridHitSummary.h"
#define INGRIDHIT_MAX_1STREDUC 1000 //## Maximum # of IngridHitSummaries 

//......................................................................
class Ingrid1stReducSummary : public TObject{
public:
    int            hitmod;  //### hit module
    int            hitcyc;  //### hit cycle
    unsigned int  xhitbit;  //### hit X layer 11(TPL) bits
    unsigned int  yhitbit;  //### hit Y layer 11(TPL) bits
    int          nhitxlyr;  //### number of hit X layers
    int          nhitylyr;  //### number of hit Y layers
    float          xtotpe;  //### total p.e. in X layers
    float          ytotpe;  //### total p.e. in Y layers
    bool       xtracklike;  //### definition of track-like will be changed or removed
    bool       ytracklike;  //### now(2010/2/12), 3 layers hit in a row 

    vector<int>      hity;  //### temporary for stdudying external(?) redioactivity
    vector<int>      hitx;
    vector<int>      hitxz;
    vector<int>      hityz;


    Ingrid1stReducSummary();
    Ingrid1stReducSummary(const Ingrid1stReducSummary& basicsum);
    virtual ~Ingrid1stReducSummary();
    void Clear   (Option_t* option="");
    void Print();
    void AddIngridHit(IngridHitSummary* hit);
    IngridHitSummary* GetIngridHit(int i) const;
    int  nhits;   
private:
    TRef fIngridHit[INGRIDHIT_MAX_1STREDUC];
    ClassDef(Ingrid1stReducSummary, 4)
};

#endif // INGRIDHITSUMMARY_H
////////////////////////////////////////////////////////////////////////

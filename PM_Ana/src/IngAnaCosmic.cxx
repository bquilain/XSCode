#ifndef  __INGANACOSMIC_CXX__
#define  __INGANACOSMIC_CXX__
#include"IngAnaCosmic.hxx"
#include"TRef.h"

IngAnaCosmic::IngAnaCosmic(){

  Clear();
}
void IngAnaCosmic::Clear(){
  for(int i=0; i<INGRIDHIT_MAXHITS; i++){
    hit[i] = new IngridHitSummary();
  }
}

bool IngAnaCosmic::SetHit( TRef inghit[] ){
  nhits = 0;
  while(1){
    if( inghit[nhits] != 0 )
      nhits++;
    else
      break;
  }
  for( int i=0; i<nhits; i++ ){
    hit[i] = (IngridHitSummary*)inghit[i].GetObject(); 
  }
}

vector<int> IngAnaCosmic::fTrigger(int view, int pln, float thr){
  vector<int> hitch;
  hitch.clear();
  for(int i=0; i<nhits; i++){
    if( hit[i] -> pln == pln  &&
	hit[i] -> view== view &&
	hit[i] ->  pe >  thr )
      hitch.push_back( hit[i]->ch );
  }
  return hitch;
}
vector<int> IngAnaCosmic::fTdcHit(int view, int pln){
  vector<int> hitch;
  hitch.clear();
  for(int i=0; i<nhits; i++){
    if( hit[i] -> pln == pln &&
	hit[i] -> view == view)
      hitch.push_back( hit[i]->ch );
  }
  return hitch;
}

#endif

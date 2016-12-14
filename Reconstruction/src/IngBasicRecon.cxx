#ifndef  __INGBASICRECON_CXX__
#define  __INGBASICRECON_CXX__
#include"IngBasicRecon.hxx"
#include"TRef.h"

IngBasicRecon::IngBasicRecon(){
  for(int i=0; i<INGRIDHIT_MAXHITS; i++){
    hit[i] = new IngridHitSummary();
  }
  //### for rock muon study
  char buff[300];
  for(int i=0; i<BREC_MAXLINE; i++){
    sprintf(buff,"Xl%02d",i);
    fLineX[i] = new TF1(buff,"pol1",0,130);
    sprintf(buff,"Yl%02d",i);
    fLineY[i] = new TF1(buff,"pol1",0,130);
  }
  Clear();
}
bool IngBasicRecon::Clear(){
  nhits            = -1;
  trgbit           = -1;
  layerpe          = -1;  
  layerpe_wo_tpl10 = -1;
  nactpln          = -1;
  nactpln_wo_tpl10 = -1;
  for(int i=0; i<INGRIDHIT_MAXHITS; i++){
     hit[i]->Clear("C");
  }
  vertexz          = -1;
  vertexy.clear();
  vertexx.clear();

  veto_xz = -1;
  veto_x  = -1;
  veto_yz = -1;
  veto_y  = -1;
  mupyz.clear();
  mupy.clear();
  mupxz.clear();
  mupx.clear();
  mdownyz.clear();
  mdowny.clear();
  mdownxz.clear();
  mdownx.clear();
  vetope.clear();
}

bool IngBasicRecon::SetHit( TRef inghit[] ){
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

int IngBasicRecon::fTrgbit(){
  if( nhits  == -1 ) return -1;
  if( trgbit != -1 ) return trgbit;
  trgbit = 0;
  for( Int_t i = 0 ; i < nhits ; i++){
    for( Int_t j=i+1 ; j < nhits ; j++){
      if( hit[i]->pln == hit[j]->pln ){
	if( hit[i]->view != hit[j]->view){
	  if( hit[i]->mod != 16)
	    trgbit  =  trgbit | ((0x400)>>(10-hit[i]->pln));
	  else
	    trgbit  =  trgbit | ((0x20000)>>(17-hit[i]->pln));
	  isActive[hit[i]->pln] = true;
	  //if( MUpAct > hitclster[i].pln && hitclster[i].pln>=1 )
	}//view difference
      }//pln coince
    }//j
  }//i
  return trgbit;  
}

int IngBasicRecon::fNactpln(){
  if( nhits   == -1     ) return -1.e-5;
  if( nactpln != -1     ) return nactpln;
  if( trgbit  == -1     ) return -1;
  nactpln = 0;
  nactpln_wo_tpl10 = 0;
  //for(Int_t i=1;i<=1024;i=i<<1){
  for(Int_t i=1;i<=131072;i=i<<1){
    if(trgbit & i){
      if(i!=1)nactpln++;
      nactpln_wo_tpl10++;
    }
  }
  return nactpln;
}

float IngBasicRecon::fLayerpe(){
  if( nhits   == -1     ) return -1.e-5;
  if( layerpe != -1     ) return layerpe;
  if( trgbit  == -1     ) return -1.e-5;
  if( nactpln == -1     ) return -1.e-5;
  layerpe          = 0;
  layerpe_wo_tpl10 = 0;
  for(Int_t i = 0 ; i < nhits ; i++){

    if( ((hit[i]->mod != 16) &&(trgbit & (  (0x400)>>(10 - hit[i]->pln) )))||
	((hit[i]->mod == 16) &&(trgbit & (  (0x20000)>>(17 - hit[i]->pln) )))
	){
      if( hit[i]->pln != 0 )
	layerpe         += hit[i]->pe;
      layerpe_wo_tpl10  += hit[i]->pe;
    }
  }
  
  if( layerpe > 0 )layerpe  = layerpe / ( nactpln * 2 );   //#### Updated at 091229
  if( layerpe_wo_tpl10 > 0 )layerpe_wo_tpl10  = layerpe_wo_tpl10 / ( nactpln_wo_tpl10 * 2 );   //#### Updated at 091229

  return layerpe;
}


//____________________________________________________________________
//____________________________________________________________________

bool IngBasicRecon::fUVETOHit(){
  upstreamveto = false;

  vpe = 0;
  float tvpe;
  //### check most upstream z position of hit ch. of VETO ##
  //########################################################
  bool   flag_hitveto = false;
  float  veto_posz = 200;
  for( int i=0; i < nhits; i++ ){
    int   m  = hit[i]->mod;
    int   p  = hit[i]->pln;
    float z  = hit[i]->z;
    float pe = hit[i]->pe;
    if( ((m != 16 && ( p == 0 || p == 11 || p == 12 || p == 13 || p == 14)) ||
	 (m == 16 && ( p == 0 || p == 18 || p == 19 || p == 20 || p == 21)) )
	&& pe > UVETOthreshold){
      if( veto_posz > z ){
	veto_posz     = z;
	flag_hitveto  = true;
	tvpe = hit[i]->pe;
	if(tvpe>vpe){
	  vpe   = tvpe;
	  vview = hit[i]->view;
	  vch   = hit[i]->ch;
	  vpln  = p;
	}
	int v = hit[i]->view;
	if(v==0){//Y
	  veto_y  = hit[i]->xy;
	  veto_yz = hit[i]->z;
	}
	if(v==1){//Y
	  veto_x  = hit[i]->xy;
	  veto_xz = hit[i]->z;
	}
      }
    }
  }
  //### check most upstream z position of active TPLs   ####
  //########################################################
  bool   flag_hittpl = false;
  float  tpl_posz    =   200;
  for( int i=0; i < nhits; i++ ){
    int   p = hit[i]->pln;
    float z = hit[i]->z;

    if((((hit[i]->mod != 16)&&(trgbit & (  (0x400)>>(10 - hit[i]->pln) )))||
	((hit[i]->mod == 16)&&(trgbit & (  (0x20000)>>(17 - hit[i]->pln) ))))&&
       p != 0 ){
      if( tpl_posz > z );    
      tpl_posz    =    z;
      flag_hittpl = true;
    }
  }
  if( flag_hitveto && flag_hittpl &&
      veto_posz    <  tpl_posz
      )
    upstreamveto = true;
	

  return upstreamveto;
}

//____________________________________________________________________
//____________________________________________________________________

bool IngBasicRecon::fGetVertex(){
  //### vertex Z
  vertexz = 17;
  for( int i=0; i < nhits; i++ ){
    int   p = hit[i]->pln;
    if((((hit[i]->mod != 16)&&(trgbit & (  (0x400)>>(10 - hit[i]->pln) )))||
	((hit[i]->mod == 16)&&(trgbit & (  (0x20000)>>(17 - hit[i]->pln) ))))&&
       p != 0 ){
      if( vertexz > p )
	vertexz = p;
    }
  }
  //### vertex X and Y 
  vertexx.clear();
  vertexy.clear();
  edgehit = false;
  for( int i=0; i < nhits; i++ ){
    int   p = hit[i]->pln;
    if((((hit[i]->mod != 16)&&(trgbit & (  (0x400)>>(10 - hit[i]->pln) )))||
	((hit[i]->mod == 16)&&(trgbit & (  (0x20000)>>(17 - hit[i]->pln) ))))&&
       p != 0 ){
      if( p == vertexz  ){
	int c = hit[i]->ch;
	int v = hit[i]->view;
	if( v == 1 )
	  vertexx.push_back(c);
	if( v == 0 )
	  vertexy.push_back(c);
	if( ((hit[i]->mod != 16) && (c == 0 || c == 23))||
	    ((hit[i]->mod == 16) && (c == 0 || c == 31)))
	  edgehit = true;
      }
    }
  }
  return true;
}

//____________________________________________________________________
//____________________________________________________________________


bool IngBasicRecon::Print(){
  if( nhits == -1 ){
    cout << "Please SetHit(TRef(IngridHitSummary*))" << endl;
    return true;
  }
  else if( nhits >= 0 ){
    cout << "----------- # of hits = " << nhits 
	 << " ---------- "             
	 << endl;
    for( int i=0; i<nhits; i++ ){
      cout << hit[i] -> pe << endl;
    }
  }
}

bool IngBasicRecon::PrintAll(){
  for(int i=0; i<nhits; i++){
    cout << "mod:"     << hit[i]->mod 
	 << "\tview:"  << hit[i]->view
	 << "\tpln:"   << hit[i]->pln
	 << "\tch:"    << hit[i]->ch
	 << "\tpe:"    << hit[i]->pe
	 <<endl;
  }
}


//### for rock muon study
bool IngBasicRecon::GetMU_MD(){
  most_upstream_tpl   = 17;
  most_downstream_tpl =  0;
  for(Int_t i = 0 ; i < nhits ; i++){
    int p = hit[i]->pln;
    if( p == 0)
      continue;
    if(((hit[i]->mod != 16)&&(trgbit & (  (0x400)>>(10 - p) )))||
       ((hit[i]->mod == 16)&&(trgbit & (  (0x20000)>>(17 - p) )))){
      if( p < most_upstream_tpl )
	most_upstream_tpl = p;
      if( p > most_downstream_tpl )
	most_downstream_tpl = p;
    }
  }
  for(Int_t i = 0 ; i < nhits ; i++){
    int   p  = hit[i] -> pln;
    int   v  = hit[i] -> view;
    float z  = hit[i] -> z;
    float xy = hit[i] -> xy;
    if( p == most_upstream_tpl ){
      if( v == 0){
	mupx. push_back(xy);
	mupxz.push_back(z);
      }
      if( v == 1){
	mupy. push_back(xy);
	mupyz.push_back(z);
      }
    }
    if( p == most_downstream_tpl ){
      if( v == 0){
	mdownx. push_back(xy);
	mdownxz.push_back(z);
      }
      if( v == 1){
	mdowny. push_back(xy);
	mdownyz.push_back(z);
      }

  
    }
  }
}


bool IngBasicRecon::GetLine(){
  nlx = 0;
  nly = 0;
  for(int i=0; i<mupx.size(); i++){
    for(int j=0; j<mdownx.size(); j++){
      float z0    = mupxz[i];
      float z1    = mdownxz[j];
      float x0    = mupx[i];
      float x1    = mdownx[j];
      float slope = 1.0 * ( x1 - x0 )/( z1 - z0 );
      float cons  = x0 - slope * z0;
      fLineX[nlx] -> SetParameter(0, cons);
      fLineX[nlx] -> SetParameter(1, slope);
      nlx++;
    }
  }
  for(int i=0; i<mupy.size(); i++){
    for(int j=0; j<mdowny.size(); j++){
      float z0    = mupyz[i];
      float z1    = mdownyz[j];
      float x0    = mupy[i];
      float x1    = mdowny[j];
      float slope = 1.0 * ( x1 - x0 )/( z1 - z0 );
      float cons  = x0 - slope * z0;
      fLineY[nly] -> SetParameter(0, cons);
      fLineY[nly] -> SetParameter(1, slope);
      nly++;
    }
  }

}



int IngBasicRecon::fActInARow(){
  inarowbit=0;
  int inarowbit_temp=0;
  int nActInaRow_temp = 0;
  int nActInaRow      = 0;
  //for(Int_t i=1;i<=1024;i=i<<1){
  for(Int_t i=1;i<=131072;i=i<<1){
    if( ( trgbit & i ) &&
	i != 1 ){ // eliminate TPL0
      nActInaRow_temp++;
      inarowbit_temp = inarowbit_temp | trgbit;
    }
    else{
      if( nActInaRow < nActInaRow_temp){
	nActInaRow =  nActInaRow_temp;
	inarowbit  =  inarowbit_temp;
      }
      nActInaRow_temp = 0;
    }
  }
  if( nActInaRow < nActInaRow_temp){
    nActInaRow =  nActInaRow_temp;
    inarowbit  =  inarowbit_temp;
  }
  return nActInaRow;
}

#endif


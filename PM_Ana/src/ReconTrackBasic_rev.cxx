#ifndef _RECONTRACKBASIC_CXX
#define _RECONTRACKBASIC_CXX
//#define DEBUG111206

#include"ReconTrackBasic_rev.hxx"

ReconTrackBasic::ReconTrackBasic(){
  for(int i=0; i<INGRIDHIT_MAXHITS; i++){
    hit[i]     = new IngridHitSummary();
  }
  for(int i=0; i<max_ntrack; i++){
    for(int v=0; v<2; v++){
      ftrkxy[v][i] = new IngridTrackSummary();
    }
    ftrk[i] = new IngridTrackSummary();
  }
  Clear();
}

ReconTrackBasic::~ReconTrackBasic(){
}

bool ReconTrackBasic::Clear(){
  for(int v=0; v<2; v++){
    vpe[v] = 0;
    n_trk[v] = 0;
    ftplistrk[v]   = false;
    uvetohashit[v] = false;
    uvetoistrk[v]  = false;
    uedgeistrk[v]  = false;
  }

  for(int i=0; i<INGRIDHIT_MAXHITS; i++){
    hit[i]      = 0;
  }
  for(int i=0; i<max_ntrack; i++){
    for(int v=0; v<2; v++){
      ftrkxy[v][i]->Clear();
    }
    ftrk[i]->Clear();
  }
  nhits     = 0;
  n_match   = 0;
}
//bool ReconTrackBasic::SetHit( TRef inghit[] ){
//bool ReconTrackBasic::SetHit( TRef inghit[] , int nhit){
bool ReconTrackBasic::SetHit( TRef inghit[INGRIDHIT_MAXHITS] , int nhit){

  nhits = nhit;
  /*
  while(1){
    if( inghit[nhits] != 0 )
      nhits++;
    else
      break;
  }
  */
  //### Get IngridHitSummarys #######
  //#################################
  for(int iview=0; iview<2; iview++)nhitsxy[iview]=0;
  for( int i=0; i<nhits; i++ ){
    hit[i] = (IngridHitSummary*)inghit[i].GetObject(); 
    nhitsxy[hit[i]->view]++;
  }
  mod = hit[0]->mod;
  //### Get md_pln(most down stream active TPL) ###
  //###############################################
  for(int p = nTPL - 1; p >= 0; p-- ){
    vector<int> h = fTdcHit(0, p);
    vector<int> v = fTdcHit(1, p);
    if(h.size()>0 && v.size()>0){
      md_pln[0] = p;
      md_pln[1] = p;
      break;
    }
  }
  //### Get md_pln(next most down stream active TPL) ###
  //###############################################
  for(int p = md_pln[0]; p >= 0; p-- ){
    vector<int> h = fTdcHit(0, p);
    vector<int> v = fTdcHit(1, p);
    if(h.size()>0 && v.size()>0){
      md_pln2[0] = p;
      md_pln2[1] = p;
      break;
    }
  }

}

vector<int> ReconTrackBasic::fTdcHit(int view, int pln){
  vector<int> hitch;
  hitch.clear();
  for(int i=0; i<nhits; i++){
    if( hit[i] -> pln == pln &&
	hit[i] -> view == view)
      hitch.push_back( hit[i]->ch );
  }
  return hitch;
}

vector<float> ReconTrackBasic::fPeHit(int view, int pln){
  vector<float> hitpe;
  hitpe.clear();
  for(int i=0; i<nhits; i++){
    if( hit[i] -> pln == pln &&
	hit[i] -> view == view)
      hitpe.push_back( hit[i]->pe );
  }
  return hitpe;
}




int ReconTrackBasic::Last3TPL(int v){//Return last(most upstream) TPL
  trk_point[v].clear();
  //#### Get Last active 3 TPLs ####
  //################################
  int            p = md_pln[v];
  vector<int>   ch[3];
  vector<float> pe[3];
  int         mpln[3]; //### last 
  int temp =0;

  //#### Get Last 3 hit layers from    ###
  //#### most down stream active plane ###
  bool gap = false; //### gap flat is added at 2010/4/29
  while( 1 ){
    
    /*
    vector<int> tempx, tempy;
    tempx.clear(); tempy.clear();
    tempx = fTdcHit( FromX, p );
    tempy = fTdcHit( FromY, p );
    if( tempx.size()>0 && tempy.size()>0 ){
    */
    //##### modified 2010/4/28
    vector<int> tempch;
    tempch.clear();
    tempch = fTdcHit( v, p );
    if( tempch.size()>0 ){

      ch  [temp] = fTdcHit( v, p );
      pe  [temp] = fPeHit ( v, p );
      mpln[temp] = p;
      //cout << temp << "'s down stream pln# " << p << "\tch"; 
      //for(int ich=0; ich<ch[temp].size(); ich++)
      //cout << ch[temp][ich]<< " ";
      //cout << endl; 
      temp++;
    }
    else if(!gap){
      gap = true;
    }
    else if(gap){
      return -1;
    }
    p--;
    if( p < 0 ||
	temp > 2 )
      break;
  }
  int i,j,k;




  for(i = 0; i < ch[0].size(); i++){
    for(j = 0; j < ch[1].size(); j++){
      //int diff01 = ch[0][i] - ch[1][j];  //### modified 2010/4/28
      float diff01 = 1.0 * (ch[0][i] - ch[1][j])/(mpln[0]-mpln[1]);
      //cout << "diff01 " << diff01 << "\t" ; 
      for( k = 0; k < ch[2].size(); k++){
	//int diff12 = ch[1][j] - ch[2][k];  //### modified 2010/4/28
	float diff12 = 1.0 * (ch[1][j] - ch[2][k])/(mpln[1]-mpln[2]);
	//cout << "diff12 " << diff12 << endl;
	if( diff01 == diff12 ){
	  //### 2009/5/31 added in order to eliminate bad fit
	  int lopech = 0;
	  float meanpe = (pe[0][i]+pe[1][j]+pe[2][k])/3;
	  if( pe[0][i] < 6.5 )lopech++;
	  if( pe[1][j] < 6.5 )lopech++;
	  if( pe[2][k] < 6.5 )lopech++;
	  //cout << pe[0][i] << "\t"<< pe[1][j] << "\t"<< pe[2][k] <<endl;
	  //if(lopech<=1)

	  if(lopech<=1 || ((lopech>=2)&&fabs(diff01)<4) )
	    //#### 2009/5/31 added end
	    //goto Find;
	    //#### 2011/
	    {
	      if(!(this_3TPL_is_already_found(v,ch[0][i],ch[1][j],ch[2][k])))
		goto Find;
	    }
	}
      }
    }
  }

  for(i=0; i<ch[0].size(); i++){
    for(j=0; j<ch[1].size(); j++){
      //int diff01 = ch[0][i] - ch[1][j];
      float diff01 = 1.0 * (ch[0][i] - ch[1][j])/(mpln[0]-mpln[1]);
      //cout << "diff01 " << diff01 << "\t" ; 
      for(k=0; k<ch[2].size(); k++){
	//int diff12 = ch[1][j] - ch[2][k];
	float diff12 = 1.0 * (ch[1][j] - ch[2][k])/(mpln[1]-mpln[2]);
	//cout << "diff12 " << diff12 << endl;
	if( fabs( diff01 - diff12 ) <= TolOfPoint ){
	  //### 2009/5/31 added in order to eliminate bad fit
	  int   lopech = 0;
	  if( pe[0][i] < 6.5 )lopech++;
	  if( pe[1][j] < 6.5 )lopech++;
	  if( pe[2][k] < 6.5 )lopech++;

	  if(lopech<=1 || ((lopech>=2)&&fabs(diff01)<4) ){
	    //#### 2009/5/31 added end
	    //goto Find;
	    if(!(this_3TPL_is_already_found(v,ch[0][i],ch[1][j],ch[2][k])))
	      goto Find;

	  }

	}
      }
    }
  }
  
  return -1;
 Find:
  trk_point[v].push_back( pair<int, int>( mpln[0], ch[0][i] ) );
  trk_point[v].push_back( pair<int, int>( mpln[1], ch[1][j] ) );
  trk_point[v].push_back( pair<int, int>( mpln[2], ch[2][k] ) );
  return mpln[2];
}

bool ReconTrackBasic::AddTrkPoint(int v, int cpln){
  int       np = trk_point[v].size();
  //### calculate channel difference(slope) of current line ###
  //###########################################################

  float chdiff = 
    1.0 * 
    ( trk_point[v][0].second - trk_point[v][np-1].second )
    /  ( trk_point[v][0].first - trk_point[v][np-1].first );

  vector<int>   ch;
  ch = fTdcHit(v, cpln);
 
  int i;
  //### calculate channel difference b/w this cpln and last track point ###
  //#######################################################################
  for( i=0; i < ch.size(); i++ ){
    float diff = 
      1.0 * 
      ( trk_point[v][np-1].second - ch[i] )
      /  ( trk_point[v][np-1].first - cpln );
    if( fabs( diff - chdiff ) <= 0.5  )
      goto Find;
  }
  for( i=0; i < ch.size(); i++ ){
    float diff = 
      1.0 * 
      ( trk_point[v][np-1].second - ch[i] )
      /  ( trk_point[v][np-1].first - cpln );
    if( fabs( diff - chdiff ) <= TolOfPoint + 0.5  )
      goto Find;
  }
  return false;
 Find:
  trk_point[v].push_back( pair<int, int>( cpln, ch[i] ) );
  return true;
}


bool ReconTrackBasic::AddDownStreamPoint(int v, int cpln){
 

  int       np = trk_point[v].size();
  //### calculate channel difference(slope) of current line ###
  //###########################################################

  float chdiff = 
    1.0 * 
    ( trk_point[v][0].second - trk_point[v][np-1].second )
    /  ( trk_point[v][0].first - trk_point[v][np-1].first );
  vector<int>   ch;
 
  ch = fTdcHit(v, cpln);
  int i;
  //### calculate channel difference b/w this cpln and last track point ###
  //#######################################################################
  for( i=0; i < ch.size(); i++ ){
    float diff = 
      1.0 * 
      ( trk_point[v][0].second - ch[i] )
      /  ( trk_point[v][0].first - cpln );
 
    if( fabs( diff - chdiff ) <= 0.5  )
      goto Find;
  }
  for( i=0; i < ch.size(); i++ ){
    float diff = 
      1.0 * 
      ( trk_point[v][0].second - ch[i] )
      /  ( trk_point[v][0].first - cpln );
    if( fabs( diff - chdiff ) <= TolOfPoint + 0.5  )
      goto Find;
  }
  return false;
 Find:
 
  trk_point[v].push_back( pair<int, int>( cpln, ch[i] ) );
  return true;
}

bool ReconTrackBasic::Point2Trk(int v){

  float vtxf[3] = {0,0,-1000};
  float vtxi[3] = {0,0, 1000};

  for(int i=0; i<trk_point[v].size(); i++){
    for(int j=0; j<nhits; j++){
      if( hit[j] -> view != v )
	continue;
      if( trk_point[v][i].first  == hit[j] -> pln &&
	  trk_point[v][i].second == hit[j] -> ch ){
	ftrkxy[v][n_trk[v]] -> AddIngridHit( hit[j] );

	if( hit[j]->z > vtxf[2] ){                     //### last point
	  ftrkxy[v][n_trk[v]] -> vtxf[0] = hit[j] -> xy;
	  ftrkxy[v][n_trk[v]] -> vtxf[2] = hit[j] -> z;
	  vtxf[2]  = hit[j]->z;
	  vtxf[0]  = hit[j]->xy;
	  endpln[v]= hit[j]->pln;
	  if( hit[j]->pln < nTPL )
	    endtpl[v] = hit[j]->pln;

	  endch[v] = hit[j]->ch;


	  if( hit[j]->pln >= nTPL )  //### tpl or veto
	    endisveto[v] = true;
	  else
	    endisveto[v] = false;
	}
	if( hit[j]->z < vtxi[2] ){ //### initial point
	  ftrkxy[v][n_trk[v]] -> vtxi[0] = hit[j] -> xy;
	  ftrkxy[v][n_trk[v]] -> vtxi[2] = hit[j] -> z;
	  vtxi[2]     = hit[j]->z;
	  vtxi[0]     = hit[j]->xy;
	  vpe[v]      = hit[j]->pe;
	  startpln[v] = hit[j]->pln;
	  if( hit[j]->pln < nTPL )
	    starttpl[v] = hit[j]->pln;
	  startch[v]  = hit[j]->ch;
	  if( hit[j]->pln >= nTPL )  //### tpl or veto
	    startisveto[v] = true;
	  else
	    startisveto[v] = false;

	}

	break;
      }// it is track point
    }
  }


  (ftrkxy[v][n_trk[v]])->length = 
    1.0 * sqrt( ( vtxf[0] - vtxi[0] )*( vtxf[0] -vtxi[0] ) + ( vtxf[2] - vtxi[2] )*( vtxf[2] -vtxi[2] ) ) ; 
    //pow( (vtxf[0] - vtxi[0]),2 ) + pow( (vtxf[2] - vtxi[2]),2 );

  FitTrack(*ftrkxy[v][n_trk[v]]);
  n_trk[v]++;

  //#### New fiducial definition
  if( !startisveto[v] && 
      !endisveto[v]   &&
      fidcosmicch[startpln[v]][0] <= startch[v] && 
      startch[v] <= fidcosmicch[startpln[v]][1] &&
      fidcosmicch[endpln[v]][0]   <= endch[v]   && 
      endch[v] <= fidcosmicch[endpln[v]][1] 
      )
    newfidcosmic[v] = true;
  else
    newfidcosmic[v] = false;



}

bool ReconTrackBasic::FitTrack(IngridTrackSummary& trk){
  float Sxx=0, Sx=0, Sxy=0, Sy=0, S=0;
  int nhits = trk.NHits();
 
  for( Int_t i = 0 ; i < nhits ; i++){
    IngridHitSummary* hita = (IngridHitSummary*)trk.GetIngridHit(i);
    float              xy  = hita -> xy;
    float               z  = hita ->  z;
 
    S++;
    Sx  += z;
    Sxx += pow(z,2);
    Sy  += xy;
    Sxy += xy * z;
  }//i
  float a = ( S * Sxy - Sx * Sy )/( S * Sxx - pow(Sx, 2) );
  float b = ( Sxx * Sy - Sxy * Sx )/( S * Sxx - pow(Sx, 2) );
  float chi2=0;
  //cout << "*** slope ;" << a << endl;
  for( Int_t i = 0 ; i < nhits ; i++){
    IngridHitSummary* hita = (IngridHitSummary*)trk.GetIngridHit(i);
    float              xy  = hita -> xy;
    float               z  = hita ->  z;
    chi2 = chi2 + pow( xy - ( a * z + b ), 2 ); 
  }//i
  trk. tx    = a;
  trk. etx   = b;
  trk. chi2x = chi2/nhits;

  float da = xyerror * sqrt( S   / ( S * Sxx - Sx * Sx ) );
  float db = xyerror * sqrt( Sxx / ( S * Sxx - Sx * Sx ) );
  trk. covx  = da;
  trk. covy  = db;
}

bool ReconTrackBasic::FirstTPLisTrack(int v){
  int np = trk_point[v].size();
  for( int i=0; i<np; i++ ){
    if( trk_point[v][i].first == 0 )
      return true;
 
  }
  return false;
}
bool ReconTrackBasic::UVETOisTrack(int v){ 
  return uvetoistrk[v];

}
bool ReconTrackBasic::UVETOhasHit(int v){
  int z0     = GetInitialZ (v);
  int xy0    = GetInitialXY(v);
  float zpos = 1.0 * z0 * ( IronThick + PlnThick ); 
 
  //if( xy0 != 0 && xy0 != 23 )
  if( xy0 <= 1 && xy0 >= 22 )
    return false;
 
  for( int i = 0; i<nhits; i++ ){
    int   p = hit[i] -> pln;
    float z = hit[i] -> z;
    if( z < zpos ){
      if( xy0 == 0  ){ //#### Bottom or Right VETO 
	if( ( p == BVETO && v == FromX ) ||
	    ( p == RVETO && v == FromY ) )
	  return true;
	  }
      if( xy0 == 23 ){ //#### Upper or Left VETO 
	if( ( p == UVETO && v == FromX ) ||
	    ( p == LVETO && v == FromY ) )
	  return true;
      }
    }
  }
  return false;
}

bool ReconTrackBasic::UEdgeisTrack(int v){
  int z0     = GetInitialZ (v);
  int xy0    = GetInitialXY(v);
  if( xy0 == 0  ||
      xy0 == 1  ||
      xy0 == 22 ||
      xy0 == 23 ){
    uedgeistrk[v]=true;
    return true;
  }
  return false;
}

int ReconTrackBasic::GetInitialZ( int v ){
  int tpl = 11;
  for(int i=0; i<trk_point[v].size(); i++){
    int t = trk_point[v][i].first;
    if( tpl > t )
      tpl = t;
  }
  return tpl;
}
int ReconTrackBasic::GetLastZ( int v ){
  int tpl = -1;
  for(int i=0; i<trk_point[v].size(); i++){
    int t = trk_point[v][i].first;
    if( tpl < t )
      tpl = t;
  }
  return tpl;
}
int ReconTrackBasic::GetInitialXY( int v ){
  int tpl = 11;
  int  ch;
  for(int i=0; i<trk_point[v].size(); i++){
    int t = trk_point[v][i].first;
    if( tpl > t ){
      tpl = t;
      ch  = trk_point[v][i].second;
    }
  }
  return ch;
}
int ReconTrackBasic::GetLastXY( int v ){
  int tpl = -1;
  int  ch;
  for(int i=0; i<trk_point[v].size(); i++){
    int t = trk_point[v][i].first;
    if( tpl < t ){
      tpl = t;
      ch  = trk_point[v][i].second;
    }
  }
  return ch;
}

bool ReconTrackBasic::AddVETOPoint(int v){

  int       np = trk_point[v].size();
  //### calculate position of trk point ####
  //########################################
  float     z0 = 1.0 * GetInitialZ  (v) * ( IronThick + PlnThick ); 
  float     z1 = 1.0 * GetLastZ     (v) * ( IronThick + PlnThick ); 
  float    xy0 = 1.0 * GetInitialXY (v) * ScintiWidth;
  float    xy1 = 1.0 * GetLastXY    (v) * ScintiWidth;
  //### calculate line equation for expected VETO position ###
  //########################################################## 
  float  slope = ( xy1 - xy0 ) / ( z1 - z0 );
  float   cons = 0.5 * ( xy1 + xy0 ) - slope * 0.5 * ( z0 + z1 );

  //### calculate position of trk point ####
  //########################################
  for(int i=0; i<nhits; i++){
    int p = hit[i] -> pln;
    if( p < nTPL )continue; //### TPL
    if( ( v == FromX && ( p == BVETO || p == UVETO ) )|| 
	( v == FromY && ( p == RVETO || p == LVETO ) )   
	){

      float  z             = hit[i] -> z;
      float  xy            = hit[i] -> xy;
      //### Now(2010/4/1) common use VETO has not valid xy   ############
      //### so correction is done by hand here               ############
      //#################################################################
      if( ( 1 <= mod && mod <= 6  )&&( p == RVETO ) )   //### Horizontal common use VETO
	xy = -19.1 - 2.5;
      //if( ( 8 <= mod && mod <= 13 )&&( p == UVETO ) ) //### Vertical common use VETO
      if( ( 8 <= mod && mod <= 13 )&&( p == BVETO ) )   //### Vertical common use VETO
	xy = -21.6;  //### Bug fix 2010/4/28
	


      //### calculate the positions at each edge of VETO scintillator ###
      //#################################################################
      float  expected_xy_u = slope * ( z - 0.5 * ScintiWidth ) + cons;
      float  expected_xy_d = slope * ( z + 0.5 * ScintiWidth ) + cons;
      if( fabs( xy - expected_xy_u ) < TolVetoFit  ||
	  fabs( xy - expected_xy_d ) < TolVetoFit ){
	int ch = hit[i] -> ch;
	trk_point[v].push_back( pair<int, int>( p, ch ) );

	if( z < z0 ){
	  //cout <<"----- upstream VETO -----" << endl;
	  //cout <<"TPL:" <<  z0 <<"\t" <<"VETO:" << z << endl;
	  uvetoistrk[v] = true;
	}
      }
      /*
      float  expected_xy = slope * z + cons;
      cout<<"VEOT z =" << z << "\t xy="<<xy << "\texpxy="<<expected_xy << endl;
      if( fabs( xy - expected_xy) < TolVetoFit ){
	int ch = hit[i] -> ch;
	trk_point[v].push_back( pair<int, int>( p, ch ) );
      }
      */
    }

  }//nhits loop
}

bool ReconTrackBasic::TrkMatching( ){//### 100404
  int diffverz = abs( GetInitialZ(0) - GetInitialZ(1) );
  if( diffverz <= TolOfTrkMatching )
    return true;
  else
    return false;
}

bool ReconTrackBasic::MargeTrkXY(int ntrk0, int ntrk1){//### Not Ready
  float vtxi[3]={0,0, 1000};
  float vtxf[3]={0,0,-1000};
  for( int i=0; i<ftrkxy[FromX][ntrk1]->NHits();i++){
    IngridHitSummary* th = 
      (IngridHitSummary*)ftrkxy[FromX][n_match]->GetIngridHit(i);
    ftrk[n_match] -> AddIngridHit( th );
    float xy = th -> xy;
    float  z = th ->  z;
    if( vtxi[Z] > z ){
   
    }
    if( vtxf[Z] < z ){
    }
  }
  for( int i=0; i<ftrkxy[FromY][ntrk0]->NHits();i++){
    IngridHitSummary* th = 
      (IngridHitSummary*)ftrkxy[FromY][n_match]->GetIngridHit(i);
    ftrk[n_match] -> AddIngridHit( th );
  }

  //####

}

bool ReconTrackBasic::ReconTrack(bool cosmic=false){
  for(int iv=0; iv<2; iv++){
    vec_trk_point [iv] .clear();
    vec_uvetoistrk[iv] .clear();
    vec_ftplistrk [iv] .clear();
    vec_mu_pln    [iv] .clear();
  }


  bool fit0 = ReconTrackXY( FromX , cosmic);
  bool fit1 = ReconTrackXY( FromY , cosmic);
  if( !(fit0) || !(fit1) )
    return false;
#ifdef DEBUG111206
  cout << "1st tracking done..." << endl;
#endif


  //get 1st track
  for(int iv=0; iv<2; iv++){
    vec_trk_point [iv].push_back(trk_point [iv]);
    vec_uvetoistrk[iv].push_back(uvetoistrk[iv]);
    vec_ftplistrk [iv].push_back(ftplistrk [iv]);
    vec_mu_pln    [iv].push_back(mu_pln    [iv]);
  }
  //try to find 2 or more track
  for(int iv=0; iv<2; iv++){
    while( ReconTrackXY(iv, cosmic) ){
      vec_trk_point [iv].push_back(trk_point [iv]);
      vec_uvetoistrk[iv].push_back(uvetoistrk[iv]);
      vec_ftplistrk [iv].push_back(ftplistrk [iv]);
      vec_mu_pln    [iv].push_back(mu_pln    [iv]);
    }
  }

  //determine which track is used
  Int_t nact[2]={11,11};
  Int_t det[2];
  for(int iv=0; iv<2; iv++){
    Int_t ns = vec_trk_point[iv].size();
    for(int is=0; is<ns; is++){

      if( vec_mu_pln[iv][is] < nact[iv] ){
	det [iv] = is;
	nact[iv] = vec_mu_pln[iv][is];
	//nact[iv] = vec_trk_point[iv][is].size();
      }
      /*
      else if( vec_mu_pln[iv][is] == nact[iv] ){
     	if( !(vec_uvetoistrk[iv][is]) && !(vec_ftplistrk[iv][is]) ){
	  det[iv]=is;
	}
      }
      */
    }

    trk_point[iv] = vec_trk_point[iv][det[iv]];
  }

#ifdef DEBUG111206
  Print();
#endif
  /*
  n_already[FromX] = 0;
  n_already[FromY] = 0;
  while( ReconTrackXY( FromX, cosmic) ){
    n_already[FromX]++;
  }
  while( ReconTrackXY( FromY, cosmic) ){
    n_already[FromY]++;
  }

  if(n_already[FromX]==0 || n_already[FromY]==0)
    return false;
  */


  Point2Trk ( FromX );
  Point2Trk ( FromY );
  //TrkMatching( 0 ); //### Now Only 1 track is reconstructed
                      //### Not Ready 
  return true;
}

bool ReconTrackBasic::ReconTrackXY(int v, bool cosmic){
  int cpln = Last3TPL(v);

  //2011/

  if( cpln != -1 ){
    //cout << "view ;" << v << " Last3TPL success" << endl;
    bool gap=false;
    cpln--;
    while(1){
      bool miss = !( AddTrkPoint(v, cpln) );
      if( gap && miss )        //## current and last plane are missed
	break;
      else if( cpln == 0 )     //## last TPL
	break;
      else if( !gap && miss )  //## if current tpln is missed, gap is true
	gap = true;	  
      else if( !miss )
	gap = false;
      cpln--;
    }
 
    mu_pln[v] = cpln;


    //### if study of VETO efficiency with cosmic #####
    //#################################################
    cpln = md_pln[v] + 1;
    gap  = false;
    while(cosmic && cpln < nTPL){
      bool miss = !( AddDownStreamPoint(v, cpln) );
      if( gap && miss )        //## current and last plane are missed
	break;
      else if( cpln == 0 )     //## last TPL
	break;
      else if( !gap && miss )  //## if current tpln is missed, gap is true
	gap = true;	  
      else if( !miss )
	gap = false;
      cpln++;
    }
    AddVETOPoint(v);
    //### if study of VETO efficiency with cosmic #####
    //#################################################

    ftplistrk[v]   = FirstTPLisTrack(v);
    uvetohashit[v] = UVETOhasHit(v);
    return true;
  }

  return false;
}

IngridTrackSummary* ReconTrackBasic::GetTrack(int v, int i) const{
  if (i < n_trk[v] && i>=0 ) return ftrkxy[v][i];
  return 0;
}


void ReconTrackBasic::Print(){
  cout <<"----------------------------------------" <<endl;
  cout << "most down stream(side,view0)" << md_pln[0] << endl;
  cout << "most down stream(top,view1)"  << md_pln[1] << endl;
  cout << "-- hit information --" << endl;
  for( int i = 0 ; i< nhits; i++){
    cout << "view:"   << hit[i] -> view
	 << "\tpln:"  << hit[i] -> pln
	 << "\tch:"   << hit[i] -> ch
	 << "\tpe:"   << hit[i] -> pe
	 << "\tdummy:"<< hit[i] -> dummy
	 <<endl;
  }  
  for( int v = 0 ; v<2; v++){
    cout<<"-- view " << v << " track info. --" << endl; 
    for( int i = 0 ; i< trk_point[v].size(); i++){
      cout << "\tpln"   << trk_point[v][i].first 
	   << "\tch"    << trk_point[v][i].second 
	   <<endl;

    }  



  }
  

  /*
  cout <<"----------------------------------------" <<endl;
  cout << "most down stream(side)" << md_pln[0] << endl;
  cout << "# of hits in view0    " << nhits0   << endl;
  cout << "View:0 # of tracks = " << n_trk0 << endl;
  for(int i=0; i<n_trk0; i++){
    cout << " *track # " << i << " ------" << endl;
    for(int j=0; j<trk_ch0_z[i].size(); j++){
      cout << "pln:"  << trk_ch0_z[i][j] 
	   << "\tch:" << trk_ch0  [i][j] 
	   << endl;
    }
  }
  */

}
#endif

#ifndef _RECONTRACK_CXX
#define _RECONTRACK_CXX


#include"ReconTrack.hxx"

ReconTrack::ReconTrack(){
  ntrack_x = 0;
  ntrack_y = 0;
  char name[300];
  for(int n = 0; n < max_ntrack; n++){
    sprintf(name, "track_x_%d", n);
    fFtrack_x[n] = new TF1(name,  "pol1",  0,  120);
    sprintf(name, "track_y_%d", n);
    fFtrack_y[n] = new TF1(name,  "pol1",  0,  120);
    fFtrack_x[n] -> SetLineWidth(0.5);
    fFtrack_x[n] -> SetLineColor(2);
    fFtrack_y[n] -> SetLineWidth(0.5);
    fFtrack_y[n] -> SetLineColor(2);
  }
  for(int n = 0; n < THETA_RESOLUTION; n++){
    sn[n] = TMath::Sin(pai*n);
    cs[n] = TMath::Cos(pai*n);
  }

  fH2_theta_rho_x = new TH2F("fH2_theta_rho_x",
			     "fH2_theta_rho_x",
			     THETA_RESOLUTION, 0, THETA_RESOLUTION,
			     RHO_RESOLUTION_H,   0, RHO_RESOLUTION);
  fH2_theta_rho_y = new TH2F("fH2_theta_rho_y",
			     "fH2_theta_rho_y",
			     THETA_RESOLUTION, 0, THETA_RESOLUTION,
			     RHO_RESOLUTION_H,   0, RHO_RESOLUTION);
  for(int v=0; v<2; v++){
    for(int i=0; i<max_ntrack; i++){
      ftrk[v][i] = new IngridTrackSummary();
    }
  }
}

ReconTrack::~ReconTrack(){
  delete fFtrack_x;
  delete fFtrack_y;
}


//____________________________________________________________________
void ReconTrack::Reset(){
 if( fH2_theta_rho_x ){
    fH2_theta_rho_x -> Reset();
    delete fH2_theta_rho_x;
    fH2_theta_rho_x = new TH2F("fH2_theta_rho_x",
			       "fH2_theta_rho_x",
			       THETA_RESOLUTION, 0, THETA_RESOLUTION,
			       RHO_RESOLUTION,   0, RHO_RESOLUTION);

 }
 if( fH2_theta_rho_y ){
    fH2_theta_rho_y -> Reset();
    delete fH2_theta_rho_y;
    fH2_theta_rho_y = new TH2F("fH2_theta_rho_y",
			       "fH2_theta_rho_y",
			       THETA_RESOLUTION, 0, THETA_RESOLUTION,
			       RHO_RESOLUTION,   0, RHO_RESOLUTION);
  }
 for( int n=0; n<max_ntrack; n++ ){
   //if(fFtrack_x[n])fFtrack_x[n] -> Clear();
   //if(fFtrack_y[n])fFtrack_y[n] -> Clear();
 }
 ntrack_x = 0;
 ntrack_y = 0;
}
//____________________________________________________________________
bool ReconTrack::HoughTrans(IngridBasicReconSummary& reduc){
 
  //### Hough transformation        ###############################
  //### (transformation from x-y space to theta-rho space)  #######
  //###############################################################
  int nhit    = reduc.nhits;
  int trgbit  = reduc.inarowbit;
  for(int n=0; n<nhit; n++){
    IngridHitSummary* inghit = (IngridHitSummary*)reduc.GetIngridHit(n);
    float                 xy = inghit -> xy;
    float                  z = inghit -> z;
    int                    p = inghit -> pln;

    if(p>=11)continue;
    if( !( p & trgbit ) || p == 0 )continue;
    for(theta = 0; theta<THETA_RESOLUTION; theta++){
      rho = (int)( z * cs[theta] + xy * sn[theta] + 0.5);
      if(inghit->view ==1){
	fH2_theta_rho_x -> Fill( theta, rho + RHO_RESOLUTION_H );
      }

      if(inghit->view ==0){
	fH2_theta_rho_y -> Fill( theta, rho + RHO_RESOLUTION_H );
      }

    }//##theta
  }//nhit
}

bool ReconTrack::GetMaxInHoughSpace(){
 
  for(int v=0; v<2; v++){
    n_max[v]   = 0;
    ent_max[v].   clear(); 
    theta_max[v]. clear();
    rho_max[v].   clear();
  }

  int tm;
  int rm;
  int m;
  for(int v=0; v<2; v++){
    while(1){
      if(v==1){
	fH2_theta_rho_x     -> GetMaximumBin(tm, rm, m);
	m = fH2_theta_rho_x -> GetBinContent(tm, rm);
      }
      if(v==0){
	fH2_theta_rho_y     -> GetMaximumBin(tm, rm, m);
	m = fH2_theta_rho_y -> GetBinContent(tm, rm);
      }
      if( m < isReconTrackN )
	break;


      for(int r=-reset_region_rho; r<reset_region_rho; r++){
	for(int t=-reset_region_theta; t<reset_region_theta; t++){
	  int tcut = tm + t;
	  int rcut = rm + r;
	  if( tcut <= 0 ){
	    tcut += THETA_RESOLUTION;
	    rcut  = -rcut;
	  }
	  if( tcut > THETA_RESOLUTION ){
	    tcut -= THETA_RESOLUTION;
	    rcut  = -rcut;
	  }
	  if(v==1)
	    fH2_theta_rho_x -> SetBinContent(tcut, rcut, 0);
	  if(v==0)
	    fH2_theta_rho_y -> SetBinContent(tcut, rcut, 0);
	}
      }
      ent_max[v].  push_back(m);
      theta_max[v].push_back(tm);
      rho_max[v].  push_back(rm);
      n_max[v]++;

    }//while
  }//view
  /*
  for(int r=0; r<RHO_RESOLUTION; r++){
    for(int t=0; t<THETA_RESOLUTION; t++){
      int ent_bin_x = fH2_theta_rho_x -> GetCellContent(t, r);
      if( ent_bin_x >= isReconTrackN ){
	n_max[1]++;
	ent_max[1].   push_back(ent_bin_x);
	theta_max[1]. push_back(t);
	rho_max[1].   push_back(r);
      }
      int ent_bin_y = fH2_theta_rho_y -> GetCellContent(t, r);
      if( ent_bin_y >= isReconTrackN ){
	n_max[0]++;
	ent_max[0].   push_back(ent_bin_y);
	theta_max[0]. push_back(t);
	rho_max[0].   push_back(r);
      }
    }//theta
  }//rho
  */
}

bool ReconTrack::GetTrack(IngridBasicReconSummary& reduc){
  ntrack[0] = 0;
  ntrack[1] = 0;
  int nhit    = reduc.nhits;
  int trgbit  = reduc.inarowbit;

  for( int v=0; v<2; v++){
    for( int n = 0; n < n_max[v]; n++ ){
      //cout << theta_max[v][n] << "\t"  << rho_max[v][n] << endl;
      cout<<ntrack[v]<<endl;
      ftrk[v][ ntrack[v] ] -> Clear("C");
      for(int i=0; i<nhit; i++){
	IngridHitSummary* inghit = (IngridHitSummary*)reduc.GetIngridHit(i);
	float                 xy = inghit -> xy;
	float                  z = inghit -> z;
	int                 view = inghit -> view;
	int                    p = inghit -> pln;
	if( view != v )continue;
	if( !( p & trgbit ) || p == 0 )continue;

	for(int th = 0; th < THETA_RESOLUTION; th++){
	  int rh = (int)( z * cs[th] + xy * sn[th] + 0.5);
	  rh    += RHO_RESOLUTION_H;
	  if( fabs( rh - rho_max[v][n]   ) <= ReconTolerance &&
	      fabs( th - theta_max[v][n] ) <= ReconTolerance  ){
	    ftrk[v][ ntrack[v] ] -> AddIngridHit( inghit );
	    break;
	  }
	}//theta  
      }//nhit

      if( ftrk[v][ ntrack[v] ] -> NHits() < NHitofTrack )
	continue;
      if( NHitPln( *ftrk[v][ntrack[v]] )  < NHitLyr )
	continue;

      FitTrack( *ftrk[v][ntrack[v]] );
      ntrack[v]++;
    }//n_max
  }//View
  return false;

}

void ReconTrack::debug(){
  TCanvas* cdebug = new TCanvas("cdebug","cdebug",10,10,800,800);
  cdebug->cd();
  cdebug->Divide(1,2);
  cdebug->cd(1);
  fH2_theta_rho_y->Draw("colz");
  cdebug->cd(2);
  fH2_theta_rho_x->Draw("colz");
  cdebug->Update();

}


void ReconTrack::Draw_Hough(TCanvas& c1){
  c1.cd();
  c1.Divide(2,1);
  c1.cd(2);
  fH2_theta_rho_y -> Draw("colz");
  c1.cd(1);
  fH2_theta_rho_x -> Draw("colz");
}


int ReconTrack::NHitPln(IngridTrackSummary& trk){
  int trgbit = 0;
  int nhits = trk.NHits();
 
  for( Int_t i = 0 ; i < nhits ; i++){
    IngridHitSummary* hita = (IngridHitSummary*)trk.GetIngridHit(i);
    trgbit  =  trgbit | ((0x400)>>(10-hita->pln));
  }//i
 
  int nactpln=0;
  for(Int_t i=1;i<=1024;i=i<<1){
    if(trgbit & i){
      if(i!=1)nactpln++;
    }
  }
  return nactpln;
}
bool ReconTrack::FitTrack(IngridTrackSummary& trk){
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
  for( Int_t i = 0 ; i < nhits ; i++){
    IngridHitSummary* hita = (IngridHitSummary*)trk.GetIngridHit(i);
    float              xy  = hita -> xy;
    float               z  = hita ->  z;
    chi2 += pow( xy - ( a * z + b ), 2 ); 
  }//i
  trk. tx  = a;
  trk. etx = b;
  trk. chi2x = chi2;

}

#endif

//#define DEBUG
//#define DRAW
{
  TRandom3 * R = new TRandom3();
  //TFile * _flux = new TFile("flux/for_pm_v3/macro/flux_cov_full_BNJ.root");
  //TMatrixDSym * Cov = (TMatrixDSym*) _flux->Get("flux_cor");
  TFile * _flux = new TFile("Covariance_Flux13a.root");
  //TMatrixDSym * Cov = (TMatrixDSym*) _flux->Get("flux_cor");//The covariance matix
  TH2D * Cov = (TH2D*) _flux->Get("hcov_pmnd");//The covariance matix
  //TH2D * Cov = (TH2D*) _flux->Get("hcor_pmnd");//The covariance matix
  
  int NBins=Cov->GetNbinsX();//Number of INGRID module x energy bins.
  cout<<"Total number of bins="<<NBins<<endl;
  int NRand=1e3;
  int NBinsMod=20;//Number of energy bins.
  TMatrixDSym Cov2(NBinsMod);
  int NeutrinoMode=0;  

  //1. Fill the erro vector and the correlation matrix
  for(int k1=0;k1<NBins;k1++){//Run over all the bins: INGRID module and energy.
    //cout<<"k1="<<k1<<endl;
    int m1=(int) (k1/NBinsMod);//Provide the module number
    int x=k1-NBinsMod*NeutrinoMode;//Provide the energy bin numnber
    for(int k2=0;k2<NBins;k2++){
      //cout<<"k2="<<k2<<endl;
     int m2=(int) (k2/NBinsMod);
      if(m1!=NeutrinoMode || m2!=NeutrinoMode) continue;
      else{
	int y=k2-NBinsMod*NeutrinoMode;
	cout<<x<<", "<<y<<endl;
	//Cov2(0,0)=3;
	Cov2(x,y)=Cov->GetBinContent(k1+1,k2+1);
	cout<<Cov->GetBinContent(k1+1,k2+1)<<", validation="<<Cov2(x,y)<<endl;
      }
    }
  }
#ifdef DRAW
  TCanvas * c0 = new TCanvas();
  Cov2->Draw("colz");
#endif

  //////////////////////////////////////////////////////
  TDecompChol * Chol = new TDecompChol(Cov2);
  Chol->Decompose();
  TMatrixD * Chol2 = (TMatrixD*) Chol->GetU();//Invert the matrix
  TMatrixD * Chol3 = new TMatrixD(NBinsMod,NBinsMod);

  for(int k1=0;k1<NBinsMod;k1++){
    for(int k2=0;k2<NBinsMod;k2++){
      (*Chol3)(k2,k1)=(*Chol2)(k1,k2);
    }
  }
  //TMatrixD  * Chol2 = Chol->GetU();
#ifdef DRAW
  TCanvas * c1 = new TCanvas();
  Chol3->Draw("colz");  
#endif
  TVectorD * V = new TVectorD(NBinsMod);
  TVectorD * V2;
  cout<<"Number of bins="<<NBins<<endl;

  TVectorD * Var = new TVectorD(NBinsMod);
  int NMods=(int) NBins/NBinsMod;
  char Name[256];
  TFile *fwr = new TFile("ErrorFlux_mod3.root","RECREATE");

  
  //////////////////////////////////////////////////////////
  for(int k1=0;k1<NBinsMod;k1++){
    for(int k2=0;k2<NBinsMod;k2++){
      //cout<<"k1="<<k1<<", k2="<<k2<<", Choleski ="<<Chol2(k1,k2)<<", Cov="<<Cov(k1,k2)<</*", Test Choleski="<<Chol3(k1,k2)<<*/endl;
    }
  }


  //2. DRAW THE NUMBER OF EVENTS IN THE GAUSSIAN
  for(int r=0;r<NRand;r++){ //Loop over the number of toy experiment
    if(r%100==0) cout<<r<<endl;
    
    for(int k=0;k<NBinsMod;k++){//Loop over the bin in energy
      V(k)=R->Gaus();//A random number
    }

    V2 = new TVectorD(NBinsMod);
    *V2=(*Chol3)*(*V);//Allow to draw a number in the BinsMod dimensional gaussian.
    for(int i=0;i<NBinsMod;i++){
      //cout<<"Bin="<<i<<", Value="<<(*V2)(i);
      (*Var)(i)=(*V2)(i);
      //(*Var)(i)=(*V2)(i)*Cov2(i,i);
      //Var(i)=V2(i);
    }
#ifdef DEBUG
    cout<<"New Vector: "<<endl;
    for(int i=0;i<NBinsMod;i++){
      cout<<(*Var)(i)<<", ";
    }
    cout<<endl;
#endif
    sprintf(Name,"Var[%d]",r);
    Var.Write(Name);
    //cout<<endl;
  }
  //////////////////////////////////////////////

  
  fwr->Close();
}

{
  TRandom3 * R = new TRandom3();
  //TFile * _flux = new TFile("flux/for_pm_v3/macro/flux_cov_full_BNJ.root");
  //TMatrixDSym * Cov = (TMatrixDSym*) _flux->Get("flux_cor");
  TFile * _flux = new TFile("flux/for_pm_v3/macro/flux_cov_full_BNJ.root");
  TMatrixDSym * Cov = (TMatrixDSym*) _flux->Get("flux_cor");//The covariance matix
  
  int NBins=Cov->GetNrows();
  int NRand=10000;
  int NBinsMod=43;
  TMatrixDSym Cov2(NBinsMod);
  TVectorD * ErrorMod = (TVectorD*) _flux->Get("ErrorMod");
  TVectorD ErrorMod2(NBinsMod);

  for(int k1=0;k1<NBins;k1++){
    cout<<"k1="<<k1<<endl;
    int m1=(int) (k1/NBinsMod);
    int x=k1-NBinsMod*3;
    if(m1==3){
      ErrorMod2(x)=(*ErrorMod)(k1);
    }
    for(int k2=0;k2<NBins;k2++){
      //cout<<"k2="<<k2<<endl;
     int m2=(int) (k2/NBinsMod);
      if(m1!=3 || m2!=3) continue;
      else{
	int y=k2-NBinsMod*3;
	cout<<x<<", "<<y<<endl;
	//Cov2(0,0)=3;
	cout<<(*Cov)(k1,k2)<<endl;
	Cov2(x,y)=(*Cov)(k1,k2);
      }
    }
  }
  TCanvas * c0 = new TCanvas();
  Cov2->Draw("colz");

  //Double_t tol=0.0;
 
  /*  TDecompChol Chol;
  Chol.SetMatrix(*Cov);
  Chol.Decompose();
  //TMatrixD * Chol2 = (TMatrixD*) Chol->GetU();
  TMatrixD Chol2 = Chol.GetU();
*/
  
  TDecompChol * Chol = new TDecompChol(Cov2);
  Chol->Decompose();
  TMatrixD * Chol2 = (TMatrixD*) Chol->GetU();//Invert the matrix
  //Transpose(Chol2);
  //TMatrixD Chol3 = TMatrixT::Transpose((TMatrixD*)Chol2);
  TMatrixD * Chol3 = new TMatrixD(NBinsMod,NBinsMod);

  for(int k1=0;k1<NBinsMod;k1++){
    for(int k2=0;k2<NBinsMod;k2++){
      (*Chol3)(k2,k1)=(*Chol2)(k1,k2);
    }
  }
  //TMatrixD  * Chol2 = Chol->GetU();
  TCanvas * c1 = new TCanvas();
  Chol3->Draw("colz");  

  TVectorD * V = new TVectorD(NBinsMod);
  TVectorD * V2;
  cout<<"Number of bins="<<NBins<<endl;

  TVectorD * Var = new TVectorD(NBinsMod);
  int NMods=(int) NBins/NBinsMod;
  char Name[256];
  TFile *fwr = new TFile("ErrorVarFinal.root","RECREATE");

  //////////////////////////////////////////////////////////
  for(int k1=0;k1<NBinsMod;k1++){
    for(int k2=0;k2<NBinsMod;k2++){
      //cout<<"k1="<<k1<<", k2="<<k2<<", Choleski ="<<Chol2(k1,k2)<<", Cov="<<Cov(k1,k2)<</*", Test Choleski="<<Chol3(k1,k2)<<*/endl;
    }
  }
  for(int k=0;k<NBinsMod;k++){
    cout<<"Error Mod="<<ErrorMod2(k)<<endl;
  }


  //2. DRAW THE NUMBER OF EVENTS IN THE GAUSSIAN
  for(int r=0;r<NRand;r++){ //Loop over the number of toy experiment
    if(r%100==0) cout<<r<<endl;
    
    for(int k=0;k<NBinsMod;k++){//Loop over the bin in energy
      V(k)=R->Gaus();//A random number
    }
    /*
      for(int i=0;i<NBins;i++){
      for(int j=0;j<NBins;j++){
      //cout<<i<<", "<<j<<endl;
      //cout<<Chol2(i,j)<<endl;
      }
      }
    */
    //TVectorD V2(NBinsMod);

    V2 = new TVectorD(NBinsMod);
    *V2=(*Chol3)*(*V);//Allow to draw a number in the BinsMod dimensional gaussian.
    
    for(int i=0;i<NBinsMod;i++){
      //cout<<"Bin="<<i<<", Value="<<(*V2)(i);
      (*Var)(i)=(*V2)(i)*ErrorMod2(i);
      //Var(i)=V2(i);
      //cout<<", after Error added="<<V2(i)<<", Error Mod="<<ErrorMod(i)<<endl;
    }
    
    sprintf(Name,"Var[%d]",r);
    Var.Write(Name);
    //cout<<endl;
    }
  fwr->Close();
}

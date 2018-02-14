{
  char Name[256];
  //TFile * _flux = new TFile("ErrorVarFinal.root");
  TFile * _flux = new TFile("ErrorVarFinal_mod3.root");
  int NErrors=400;
  int NBinsMod=20;
  TH1D * Flux = new TH1D("Flux","",100,-5,5);
  TH1D * FluxBin[NBinsMod];
  double Cov[NBinsMod][NBinsMod];
  double Cor[NBinsMod][NBinsMod];
  double RMS[NBinsMod];
  /*
  for(int i=0;i<NBinsMod;i++){
    sprintf(Name,"FluxBin[%d]",i);
    FluxBin[i]= new TH1D(Name,"",100,-5,5);
  }
  for(int e=0;e<NErrors;e++){
    sprintf(Name,"Var[%d]",e);
    TVectorD * Var = (TVectorD*) _flux->Get(Name);
    //cout<<endl<<endl<<"Error ="<<e<<endl;
    double Total=0;
    for(int i=0;i<NBinsMod;i++){
      FluxBin[i]->Fill(Var(i));
      Total+=Var(i);
    }
    Flux->Fill(Total);
  }
*/
  for(int i=0;i<NBinsMod;i++){
    cout<<"Bin i="<<i<<endl;
    RMS[i]=0;
    for(int j=0;j<NBinsMod;j++){
      cout<<"Bin j="<<j<<endl;
      Cov[i][j]=0;
      Cor[i][j]=0;
      for(int e=0;e<NErrors;e++){
	sprintf(Name,"Var[%d]",e);
	TVectorD * Var = (TVectorD*) _flux->Get(Name);
	Cov[i][j]+=Var(i)*Var(j);
      }
      Cov[i][j]/=NErrors;
    }
  }
  for(int i=0;i<NBinsMod;i++) RMS[i]=TMath::Sqrt(Cov[i][i]);

  for(int i=0;i<NBinsMod;i++){
    for(int j=0;j<NBinsMod;j++){
      Cor[i][j]=Cov[i][j]/(RMS[i]*RMS[j]);
    }
  }
  /*
  TCanvas * c0 = new TCanvas();
  Flux->Draw();

  TCanvas * c1[NBinsMod];
  for(int i=0;i<NBinsMod;i++){
    c1[i]= new TCanvas();
    FluxBin[i]->Draw();
  }
  */
  TCanvas * c2 = new TCanvas();
  TH2D * hCor = new TH2D("hCor","",NBinsMod,0,NBinsMod,NBinsMod,0,NBinsMod);
  for(int i=0;i<NBinsMod;i++){
    for(int j=0;j<NBinsMod;j++){
      hCor->SetBinContent(i+1,j+1,Cor[i][j]);
    }
  }
  hCor->Draw("colz");
}

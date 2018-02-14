void rebin(){

  //from XS/inc/setup.h
  const int NBinsEnergyFlux=20;
  double BinningEnergyFlux[NBinsEnergyFlux+1];
  for(int i=0;i<=NBinsEnergyFlux;i++){
    if(i<=15) BinningEnergyFlux[i]=i*0.2;//in GeV
    else if(i<=16) BinningEnergyFlux[i]=BinningEnergyFlux[15]+(i-15)*1;
    else if(i<=19) BinningEnergyFlux[i]=BinningEnergyFlux[16]+(i-16)*2;
    else if(i<=20) BinningEnergyFlux[i]=30.;
  }
  //

  TFile * f = new TFile("tune_ingbg.root","read");
  TH1D* h=(TH1D*) f->Get("nd2_tune_numu");
  h->SetDirectory(0);
  f->Close();

  h->SetName("test");

  TH1D* wh=new TH1D("nd2_tune_numu","",NBinsEnergyFlux,BinningEnergyFlux);

  int bin=1;
  for(int wbin=1;wbin<=NBinsEnergyFlux;wbin++){
    double limit_up=wh->GetBinLowEdge(wbin+1);
    double Nneut=0.;
    while(bin<220 && h->GetBinLowEdge(bin+1)<limit_up){
      Nneut+=h->GetBinContent(bin);
      bin++;
    }
    wh->SetBinContent(wbin,Nneut);
  }

  cout<<"Original #neutrinos="<<h->Integral()<<endl;
  cout<<"After rebinning="<<wh->Integral()<<endl;


  TFile * wfile=new TFile("tune_ingbg_rebin.root","recreate");
  wfile->cd();
  wh->Write();
  wfile->Close();







}

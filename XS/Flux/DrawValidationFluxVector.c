//#define DEBUG
//#define DRAW
{
  //copied from inc/setup.h. Re-copy if changed
  const int NBinsEnergyFlux=20;
  double BinningEnergyFlux[NBinsEnergyFlux+1];
  for(int i=0;i<=NBinsEnergyFlux;i++){
    if(i<=15) BinningEnergyFlux[i]=i*0.2;//in GeV
    else if(i<=16) BinningEnergyFlux[i]=BinningEnergyFlux[15]+(i-15)*1;
    else if(i<=19) BinningEnergyFlux[i]=BinningEnergyFlux[16]+(i-16)*2; // ML corr 2017/10
    else if(i<=20) BinningEnergyFlux[i]=30.;
    else cout<<"Error in binning the flux in energy. Please look at setup.h"<<endl;
  }
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  
  TFile * f = new TFile("tune_ingbg.root");
  TH1D * NeutrinoFlux = (TH1D*) f->Get("nd2_tune_numu");
  TH1D * hFlux = new TH1D("hFlux","",1e3,1e11,1e14);
  //TFile *fwr = new TFile("ErrorFlux_mod3.root","READ");
  TFile *fwr = new TFile("ErrorVarFinal_mod3.root","READ");

  int NBinsMod=20;
  cout<<"Total number of neutrino ="<<NeutrinoFlux->Integral()<<endl;
  TVectorD * FluxVector;
  for(int n=0;n<1e3;n++){
    FluxVector = (TVectorD*) fwr->Get(Form("Var[%d]",n));
    double Total = 0;

    for(int i=1;i<=NeutrinoFlux->GetNbinsX();i++){
      double Nneut=NeutrinoFlux->GetBinContent(i);
      double Enu=NeutrinoFlux->GetBinCenter(i);//center energy of the bin
      //We check in which bin the center energy is
      for(int j=0;j<NBinsEnergyFlux;j++){
	if(Enu<BinningEnergyFlux[j+1]){
	  //cout<<"Enu="<<Enu<<", bin ["<<BinningEnergyFlux[j]<<","<<BinningEnergyFlux[j+1]<<endl;
	  Total += NeutrinoFlux->GetBinContent(i)*(1+FluxVector(j));
	  break;
	}
      }
    }
    hFlux->Fill(Total);
    if( n%100 == 0) cout<<"flux = "<<Total<<endl;  
  }
  cout<<hFlux->GetRMS()/hFlux->GetMean()<<endl;
  hFlux->Fit("gaus","0");
}

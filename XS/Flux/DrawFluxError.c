//#define DEBUG
//#define DRAW
{
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  
  TFile * _flux = new TFile("Covariance_Flux13a.root");
  TH2D * Cov = (TH2D*) _flux->Get("hcov_pmnd");//The covariance matix
  TH2D * Cor = (TH2D*) _flux->Get("hcor_pmnd");//The covariance matix

  int NeutrinoMode=0;  
  int NBins=Cov->GetNbinsX();//Number of INGRID module x energy bins.
  cout<<"Total number of bins="<<NBins<<endl;
  int NBinsMod=20;//Number of energy bins.
  double BinningEnergy[NBinsMod+1];
  for(int i=0;i<=NBinsMod;i++){
    if(i<=15) BinningEnergy[i] = i*0.2;
    else if(i==16) BinningEnergy[i] = BinningEnergy[i-1] + 1.;
    else if(i<=19) BinningEnergy[i] = BinningEnergy[i-1] + 2.;
    else BinningEnergy[i] = BinningEnergy[i-1] + 20.;
  }
  TH2D * Cov_Mod3 = new TH2D("Cov_Mod3","Covariance on the INGRID central module",NBinsMod,BinningEnergy,NBinsMod,BinningEnergy);
  TH2D * Cor_Mod3 = new TH2D("Cor_Mod3","Correlation on the INGRID central module",NBinsMod,BinningEnergy,NBinsMod,BinningEnergy);
  TH1D * Error_Mod3 = new TH1D("Error_Mod3","Error on the INGRID central module",NBinsMod,BinningEnergy);

  //1. Fill the erro vector and the correlation matrix
  for(int k1=0;k1<NBins;k1++){//Run over all the bins: INGRID module and energy.
    int m1=(int) (k1/NBinsMod);//Provide the module number
    int x=k1-NBinsMod*NeutrinoMode;//Provide the energy bin numnber
    double Error=0;
    cout<<"NeutrinoMode number ="<<m1<<endl;
    
    for(int k2=0;k2<NBins;k2++){
      //cout<<"k2="<<k2<<endl;
     int m2=(int) (k2/NBinsMod);
     if(k1 == k2) cout<<"Error = "<<TMath::Sqrt(Cov->GetBinContent(k1+1,k2+1))<<endl;
     if(m1!=NeutrinoMode || m2!=NeutrinoMode) continue;
     else{
       int y=k2-NBinsMod*NeutrinoMode;
       //cout<<x<<", "<<y<<endl;
       //Cov2(0,0)=3;
       Cov_Mod3->SetBinContent(x+1,y+1,Cov->GetBinContent(k1+1,k2+1));
       Cor_Mod3->SetBinContent(x+1,y+1,Cor->GetBinContent(k1+1,k2+1));
       if(k1 == k2){
	 //cout<<"Error = "<<TMath::Sqrt(Cov->GetBinContent(k1,k2))<<endl;
	 Error_Mod3->SetBinContent(x+1,TMath::Sqrt(Cov->GetBinContent(k1+1,k2+1)));
       }
     }
    }
  }

  TCanvas * c0 = new TCanvas();
  Error_Mod3->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  Error_Mod3->GetXaxis()->SetRangeUser(0,10);
  Error_Mod3->GetYaxis()->SetTitle("Fractional error");
  Error_Mod3->SetLineWidth(2);
  Error_Mod3->GetYaxis()->SetRangeUser(0,0.2);
  Error_Mod3->Draw();
  c0->SaveAs("../plots/FracError.pdf");
  
  TCanvas * c1 = new TCanvas();
  Cor_Mod3->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  Cor_Mod3->GetYaxis()->SetTitle("E_{#nu} (GeV)");
  Cor_Mod3->GetXaxis()->SetRangeUser(0,10);
  Cor_Mod3->GetYaxis()->SetRangeUser(0,10);
  Cor_Mod3->Draw("colz");
  c1->SaveAs("../plots/CorrelationFluxMatrix.pdf");
}

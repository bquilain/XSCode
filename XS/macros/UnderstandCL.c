{
  //The first parameter is a normalisation coefficient.
  //the second parameter is the most probable value.
  //the third parameter correspond to Lambda in most book
    
  TF1 * PDF = new TF1("PDF","landau",0.,100);
  PDF->SetParameter(0,1);
  PDF->SetParameter(1,5);
  PDF->SetParameter(2,1);

    TF1 * PDF2 = new TF1("PDF2","landau",0.,100);
  PDF2->SetParameter(0,1);
  PDF2->SetParameter(1,15);
  PDF2->SetParameter(2,1);
  TCanvas * c0 = new TCanvas();
  PDF->Draw();
  PDF2->Draw("same");

  TH1D * CumulativeInvert = new TH1D("CumulativeInvert","",1000,0,100);
  double Integral=PDF->Integral(0,100);
  for(int ibinx=1;ibinx<=CumulativeInvert->GetNbinsX();ibinx++){
    double Value=PDF->Integral(0,CumulativeInvert->GetBinLowEdge(ibinx));
    CumulativeInvert->SetBinContent(ibinx,1.-(Value/Integral));
  }
  TCanvas * c1 = new TCanvas();
  CumulativeInvert->Draw();

  TRandom3 * r = new TRandom3();
  TH1D * CL_Uniform = new TH1D("CL_Uniform","",100,0.,1.);
  TH1D * CL_Landau = new TH1D("CL_Landau","",100,0.,1.);
  TH1D * CL_Landau2 = new TH1D("CL_Landau2","",100,0.,1.);

  int NPoints=1000000;
  for(int i=0;i<NPoints;i++){
    double X=r->Uniform(0,100);
    int nbin=((int) (X*10))+1;
    CL_Uniform->Fill(CumulativeInvert->GetBinContent(nbin));
    //if(CumulativeInvert->GetBinContent(nbin)>0.9) cout<<X<<", nbin="<<nbin<<", value="<<CumulativeInvert->GetBinContent(nbin)<<endl;

    double X0=PDF->GetRandom();//Estimate the #pe
    int nbin0=((int) (X0*10))+1;
    //cout<<X0<<", nbin="<<nbin0<<", value="<<CumulativeInvert->GetBinContent(nbin0)<<endl;
    CL_Landau->Fill(CumulativeInvert->GetBinContent(nbin0));
    //cout<<

    double X1=PDF2->GetRandom();//Estimate the #pe
    int nbin1=((int) (X1*10))+1;
    //cout<<X1<<", nbin="<<nbin1<<", value="<<CumulativeInvert->GetBinContent(nbin1)<<endl;
    CL_Landau2->Fill(CumulativeInvert->GetBinContent(nbin1));

  }
  TCanvas * c2 = new TCanvas();
  CL_Landau->SetLineColor(2);
  CL_Landau->SetLineWidth(2);
  CL_Landau->Draw();
  CL_Uniform->SetLineColor(1);
  CL_Uniform->SetLineWidth(2);
  CL_Uniform->Draw("same");
  CL_Landau2->SetLineColor(4);
  CL_Landau2->Draw("same");
}

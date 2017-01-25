/**/
{
  TFile * f = new TFile("../src/mucl.root");
  
  TH1F * PDF_Sci = (TH1F*) f->Get("hsci");
  TH1F * Cumulative_Sci = (TH1F*) f->Get("lsci");
  TH1F * Distribution_Sci = new TH1F("Distribution_Sci","",100,0,1);

  TH1F * PDF_Ing = (TH1F*) f->Get("hing");
  TH1F * Cumulative_Ing = (TH1F*) f->Get("ling");
  TH1F * Distribution_Ing = new TH1F("Distribution_Ing","",100,0,1);

  TH1F * PDF_Sci_Distorted = (TH1F*) PDF_Sci->Clone("PDF_Sci_Distorted");
  //TH1F * PDF_Sci_Distorted = new TH1F("PDF_Sci_Distorted","",300,0,300);
  //Same cumulative as non distorted since we are using this in pratice
  TH1F * Distribution_Sci_Distorted = new TH1F("Distribution_Sci_Distorted","",100,0,1);//Let's assumme that Kikawa-san forgot to correct something when generating PSF_Sci. Let's assume this correction will create a PDF_Sci higher than PDF_Sci_Corrected, and that the data/MC follow the second one. Let's generate this second one distribution -> normally, should make CL closer to 1.

  double Max=PDF_Sci->GetMaximumBin();
  double RMS=PDF_Sci->GetRMS();
  for(int ibinx=1;ibinx<=PDF_Sci_Distorted->GetNbinsX();ibinx++){//let's shift everything by 5 bins (five pe)
    double Value;
    if(ibinx<(PDF_Sci_Distorted->GetNbinsX()-5)) Value=PDF_Sci->GetBinContent(ibinx+5);
    else Value=0;
    double Scaling=1;//TMath::Abs(ibinx-Max);
    if(Scaling!=0) Value/=Scaling;
    //else Value=PDF_Sci->GetBinContent(PDF_Sci->GetNbinsX());
    PDF_Sci_Distorted->SetBinContent(ibinx,Value);
  }
  
  TCanvas * c0 = new TCanvas();
  c0->Divide(2);
  c0->cd(1);
  PDF_Sci->SetLineColor(1);
  PDF_Sci->DrawNormalized();
  PDF_Ing->SetLineColor(2);
  PDF_Ing->DrawNormalized("same");
  PDF_Sci_Distorted->SetLineColor(4);
  PDF_Sci_Distorted->DrawNormalized("same"); 
  c0->cd(2);
  Cumulative_Sci->Draw();
  Cumulative_Ing->SetLineColor(2);
  Cumulative_Ing->Draw("same");


  int NPoints=100000;
  for(int i=0;i<NPoints;i++){
    double PE_Sci=PDF_Sci->GetRandom();
    int ibin_Sci=((int) (PE_Sci))+1;
    double P_Sci=Cumulative_Sci->GetBinContent(ibin_Sci);
    Distribution_Sci->Fill(P_Sci);

    double PE_Ing=PDF_Ing->GetRandom();
    int ibin_Ing=((int) (PE_Ing))+1;
    double P_Ing=Cumulative_Ing->GetBinContent(ibin_Ing);
    Distribution_Ing->Fill(P_Ing);

    double PE_Sci_Distorted=PDF_Sci_Distorted->GetRandom();
    int ibin_Sci_Distorted=((int) (PE_Sci_Distorted))+1;
    double P_Sci_Distorted=Cumulative_Sci->GetBinContent(ibin_Sci_Distorted);
    Distribution_Sci_Distorted->Fill(P_Sci_Distorted);
  }
  
  TCanvas * c1 = new TCanvas();
  Distribution_Sci->DrawNormalized();
  Distribution_Ing->SetLineColor(kRed);
  Distribution_Ing->DrawNormalized("same");
  Distribution_Sci_Distorted->SetLineColor(kBlue);
  Distribution_Sci_Distorted->DrawNormalized("same");
  
  /*
  TH1D * CumulativeInvert = new TH1D("CumulativeInvert","",1000,0,100);
  double Integral=PDF->Integral(0,100);
  for(int ibinx=1;ibinx<=CumulativeInvert->GetNbinsX();ibinx++){
    double Value=PDF->Integral(0,CumulativeInvert->GetBinLowEdge(ibinx));
    CumulativeInvert->SetBinContent(ibinx,1.-(Value/Integral));
  }

  TCanvas * c2 = new TCanvas();
  CL_Landau->SetLineColor(2);
  CL_Landau->SetLineWidth(2);
  CL_Landau->DrawNormalized();
  CL_Uniform->SetLineColor(1);
  CL_Uniform->SetLineWidth(2);
  CL_Uniform->DrawNormalized("same");
  CL_Landau2->SetLineColor(4);
  CL_Landau2->DrawNormalized("same");
  */
}

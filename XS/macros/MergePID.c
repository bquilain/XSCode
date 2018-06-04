{
  bool PM=true;
  bool PiCL=true;
  char suffix[3];sprintf(suffix,(PM?"":"_WM"));
  TFile * wfile = new TFile(Form("../PID/BuggedPDFMuCL_Likelihood%s.root",suffix),"recreate");

  TH1D * PDFMuCL_PMIng_Muon;TH1D * PDFMuCL_PMSci_Muon;TH1D * PDFMuCL_Ing_Muon;
  TH1D * PDFMuCL_PMIng_NotMuon;TH1D * PDFMuCL_PMSci_NotMuon;TH1D * PDFMuCL_Ing_NotMuon;
  TH1D * PDFPiCL_PMIng_Pion;TH1D * PDFPiCL_PMSci_Pion;TH1D * PDFPiCL_Ing_Pion;
  TH1D * PDFPiCL_PMIng_NotPion;TH1D * PDFPiCL_PMSci_NotPion;TH1D * PDFPiCL_Ing_NotPion;
  
  TH1D * PDFParticle;
  
  TFile * fInput;

  for(int i=0;i<100;i++){
    fInput = new TFile(Form("../PID/PDFMuCL_Likelihood%d.root",i));
    if(i == 0){
      PDFMuCL_PMIng_Muon = (TH1D*) ((TH1D*) fInput->Get("PDFMuCL_PMIng_Muon"))->Clone("PDFMuCL_PMIng_Muon");
      PDFMuCL_PMSci_Muon = (TH1D*) ((TH1D*) fInput->Get("PDFMuCL_PMSci_Muon"))->Clone("PDFMuCL_PMSci_Muon");
      PDFMuCL_Ing_Muon = (TH1D*) ((TH1D*) fInput->Get("PDFMuCL_Ing_Muon"))->Clone("PDFMuCL_Ing_Muon");

      PDFMuCL_PMIng_NotMuon = (TH1D*) ((TH1D*) fInput->Get("PDFMuCL_PMIng_NotMuon"))->Clone("PDFMuCL_PMIng_NotMuon");
      PDFMuCL_PMSci_NotMuon = (TH1D*) ((TH1D*) fInput->Get("PDFMuCL_PMSci_NotMuon"))->Clone("PDFMuCL_PMSci_NotMuon");
      PDFMuCL_Ing_NotMuon = (TH1D*) ((TH1D*) fInput->Get("PDFMuCL_Ing_NotMuon"))->Clone("PDFMuCL_Ing_NotMuon");

      if(PiCL){
	PDFPiCL_PMIng_Pion = (TH1D*) ((TH1D*) fInput->Get("PDFPiCL_PMIng_Pion"))->Clone("PDFPiCL_PMIng_Pion");
	PDFPiCL_PMSci_Pion = (TH1D*) ((TH1D*) fInput->Get("PDFPiCL_PMSci_Pion"))->Clone("PDFPiCL_PMSci_Pion");
	PDFPiCL_Ing_Pion = (TH1D*) ((TH1D*) fInput->Get("PDFPiCL_Ing_Pion"))->Clone("PDFPiCL_Ing_Pion");
	
	PDFPiCL_PMIng_NotPion = (TH1D*) ((TH1D*) fInput->Get("PDFPiCL_PMIng_NotPion"))->Clone("PDFPiCL_PMIng_NotPion");
	PDFPiCL_PMSci_NotPion = (TH1D*) ((TH1D*) fInput->Get("PDFPiCL_PMSci_NotPion"))->Clone("PDFPiCL_PMSci_NotPion");
	PDFPiCL_Ing_NotPion = (TH1D*) ((TH1D*) fInput->Get("PDFPiCL_Ing_NotPion"))->Clone("PDFPiCL_Ing_NotPion");
      }
     PDFParticle = (TH1D*) ((TH1D*) fInput->Get("PDFParticle"))->Clone("PDFParticle");
    }
    else{
      PDFMuCL_PMIng_Muon->Add((TH1D*) fInput->Get("PDFMuCL_PMIng_Muon"));
      PDFMuCL_PMSci_Muon->Add((TH1D*) fInput->Get("PDFMuCL_PMSci_Muon"));
      PDFMuCL_Ing_Muon->Add((TH1D*) fInput->Get("PDFMuCL_Ing_Muon"));

      PDFMuCL_PMIng_NotMuon->Add((TH1D*) fInput->Get("PDFMuCL_PMIng_NotMuon"));
      PDFMuCL_PMSci_NotMuon->Add((TH1D*) fInput->Get("PDFMuCL_PMSci_NotMuon"));
      PDFMuCL_Ing_NotMuon->Add((TH1D*) fInput->Get("PDFMuCL_Ing_NotMuon"));

      if(PiCL){
	PDFPiCL_PMIng_Pion->Add((TH1D*) fInput->Get("PDFPiCL_PMIng_Pion"));
	PDFPiCL_PMSci_Pion->Add((TH1D*) fInput->Get("PDFPiCL_PMSci_Pion"));
	PDFPiCL_Ing_Pion->Add((TH1D*) fInput->Get("PDFPiCL_Ing_Pion"));
	
	PDFPiCL_PMIng_NotPion->Add((TH1D*) fInput->Get("PDFPiCL_PMIng_NotPion"));
	PDFPiCL_PMSci_NotPion->Add((TH1D*) fInput->Get("PDFPiCL_PMSci_NotPion"));
	PDFPiCL_Ing_NotPion->Add((TH1D*) fInput->Get("PDFPiCL_Ing_NotPion"));
      }
      PDFParticle->Add((TH1D*) fInput->Get("PDFParticle"));
    }
  }
  fInput->Close();

  wfile->cd();  
  PDFMuCL_PMIng_Muon->Draw();
  if(PM){
  PDFMuCL_PMIng_Muon->Scale(1./PDFMuCL_PMIng_Muon->Integral());
  PDFMuCL_PMSci_Muon->Scale(1./PDFMuCL_PMSci_Muon->Integral());
  PDFMuCL_PMIng_Muon->Write();
  PDFMuCL_PMSci_Muon->Write();
  PDFMuCL_PMIng_NotMuon->Scale(1./PDFMuCL_PMIng_NotMuon->Integral());
  PDFMuCL_PMSci_NotMuon->Scale(1./PDFMuCL_PMSci_NotMuon->Integral());
  PDFMuCL_PMIng_NotMuon->Write();
  PDFMuCL_PMSci_NotMuon->Write();

  double x_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];double ex_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];double y_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];double ey_PMIng_Muon[PDFMuCL_PMIng_Muon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMIng_Muon->GetNbinsX();ibinx++){
    x_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinCenter(ibinx);
    ex_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinWidth(ibinx)/2;
    y_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinContent(ibinx);
    ey_PMIng_Muon[ibinx-1]=PDFMuCL_PMIng_Muon->GetBinError(ibinx);
  }
  TGraphErrors * g_PMIng_Muon = new TGraphErrors(PDFMuCL_PMIng_Muon->GetNbinsX(),x_PMIng_Muon,y_PMIng_Muon,ex_PMIng_Muon,ey_PMIng_Muon);
  g_PMIng_Muon->Write("g_PMIng_Muon");
  TSpline3* s_PMIng_Muon = new TSpline3("s_PMIng_Muon",g_PMIng_Muon);
  s_PMIng_Muon->Write("s_PMIng_Muon");
  
  double x_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];double ex_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];double y_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];double ey_PMSci_Muon[PDFMuCL_PMSci_Muon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMSci_Muon->GetNbinsX();ibinx++){
    x_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinCenter(ibinx);
    ex_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinWidth(ibinx)/2;
    y_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinContent(ibinx);
    ey_PMSci_Muon[ibinx-1]=PDFMuCL_PMSci_Muon->GetBinError(ibinx);
  }
  TGraphErrors * g_PMSci_Muon = new TGraphErrors(PDFMuCL_PMSci_Muon->GetNbinsX(),x_PMSci_Muon,y_PMSci_Muon,ex_PMSci_Muon,ey_PMSci_Muon);
  g_PMSci_Muon->Write("g_PMSci_Muon");
  TSpline3 * s_PMSci_Muon = new TSpline3("s_PMSci_Muon",g_PMSci_Muon);
  s_PMSci_Muon->Write("s_PMSci_Muon");


  double x_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];double ex_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];double y_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];double ey_PMIng_NotMuon[PDFMuCL_PMIng_NotMuon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMIng_NotMuon->GetNbinsX();ibinx++){
    x_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinCenter(ibinx);
    ex_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinWidth(ibinx)/2;
    y_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinContent(ibinx);
    ey_PMIng_NotMuon[ibinx-1]=PDFMuCL_PMIng_NotMuon->GetBinError(ibinx);
  }
  TGraphErrors * g_PMIng_NotMuon = new TGraphErrors(PDFMuCL_PMIng_NotMuon->GetNbinsX(),x_PMIng_NotMuon,y_PMIng_NotMuon,ex_PMIng_NotMuon,ey_PMIng_NotMuon);
  g_PMIng_NotMuon->Write("g_PMIng_NotMuon");
  TSpline3* s_PMIng_NotMuon = new TSpline3("s_PMIng_NotMuon",g_PMIng_NotMuon);
  s_PMIng_NotMuon->Write("s_PMIng_NotMuon");
  
  double x_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];double ex_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];double y_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];double ey_PMSci_NotMuon[PDFMuCL_PMSci_NotMuon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_PMSci_NotMuon->GetNbinsX();ibinx++){
    x_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinCenter(ibinx);
    ex_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinWidth(ibinx)/2;
    y_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinContent(ibinx);
    ey_PMSci_NotMuon[ibinx-1]=PDFMuCL_PMSci_NotMuon->GetBinError(ibinx);
  }
  TGraphErrors * g_PMSci_NotMuon = new TGraphErrors(PDFMuCL_PMSci_NotMuon->GetNbinsX(),x_PMSci_NotMuon,y_PMSci_NotMuon,ex_PMSci_NotMuon,ey_PMSci_NotMuon);
  g_PMSci_NotMuon->Write("g_PMSci_NotMuon");
  TSpline3 * s_PMSci_NotMuon = new TSpline3("s_PMSci_NotMuon",g_PMSci_NotMuon);
  s_PMSci_NotMuon->Write("s_PMSci_NotMuon");

  if(PiCL){
    PDFPiCL_PMIng_Pion->Scale(1./PDFPiCL_PMIng_Pion->Integral());
    PDFPiCL_PMSci_Pion->Scale(1./PDFPiCL_PMSci_Pion->Integral());
    PDFPiCL_PMIng_Pion->Write();
    PDFPiCL_PMSci_Pion->Write();
    PDFPiCL_PMIng_NotPion->Scale(1./PDFPiCL_PMIng_NotPion->Integral());
    PDFPiCL_PMSci_NotPion->Scale(1./PDFPiCL_PMSci_NotPion->Integral());
    PDFPiCL_PMIng_NotPion->Write();
    PDFPiCL_PMSci_NotPion->Write();

    double x_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];double ex_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];double y_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];double ey_PMIng_Pion[PDFPiCL_PMIng_Pion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMIng_Pion->GetNbinsX();ibinx++){
      x_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinCenter(ibinx);
      ex_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinWidth(ibinx)/2;
      y_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinContent(ibinx);
      ey_PMIng_Pion[ibinx-1]=PDFPiCL_PMIng_Pion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMIng_Pion = new TGraphErrors(PDFPiCL_PMIng_Pion->GetNbinsX(),x_PMIng_Pion,y_PMIng_Pion,ex_PMIng_Pion,ey_PMIng_Pion);
    g_PMIng_Pion->Write("g_PMIng_Pion");
    TSpline3* s_PMIng_Pion = new TSpline3("s_PMIng_Pion",g_PMIng_Pion);
    s_PMIng_Pion->Write("s_PMIng_Pion");

    double x_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];double ex_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];double y_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];double ey_PMSci_Pion[PDFPiCL_PMSci_Pion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMSci_Pion->GetNbinsX();ibinx++){
      x_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinCenter(ibinx);
      ex_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinWidth(ibinx)/2;
      y_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinContent(ibinx);
      ey_PMSci_Pion[ibinx-1]=PDFPiCL_PMSci_Pion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMSci_Pion = new TGraphErrors(PDFPiCL_PMSci_Pion->GetNbinsX(),x_PMSci_Pion,y_PMSci_Pion,ex_PMSci_Pion,ey_PMSci_Pion);
    g_PMSci_Pion->Write("g_PMSci_Pion");
    TSpline3 * s_PMSci_Pion = new TSpline3("s_PMSci_Pion",g_PMSci_Pion);
    s_PMSci_Pion->Write("s_PMSci_Pion");

    double x_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];double ex_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];double y_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];double ey_PMIng_NotPion[PDFPiCL_PMIng_NotPion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMIng_NotPion->GetNbinsX();ibinx++){
      x_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinCenter(ibinx);
      ex_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinWidth(ibinx)/2;
      y_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinContent(ibinx);
      ey_PMIng_NotPion[ibinx-1]=PDFPiCL_PMIng_NotPion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMIng_NotPion = new TGraphErrors(PDFPiCL_PMIng_NotPion->GetNbinsX(),x_PMIng_NotPion,y_PMIng_NotPion,ex_PMIng_NotPion,ey_PMIng_NotPion);
    g_PMIng_NotPion->Write("g_PMIng_NotPion");
    TSpline3* s_PMIng_NotPion = new TSpline3("s_PMIng_NotPion",g_PMIng_NotPion);
    s_PMIng_NotPion->Write("s_PMIng_NotPion");

    double x_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];double ex_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];double y_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];double ey_PMSci_NotPion[PDFPiCL_PMSci_NotPion->GetNbinsX()];
    for(int ibinx=1;ibinx<=PDFPiCL_PMSci_NotPion->GetNbinsX();ibinx++){
      x_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinCenter(ibinx);
      ex_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinWidth(ibinx)/2;
      y_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinContent(ibinx);
      ey_PMSci_NotPion[ibinx-1]=PDFPiCL_PMSci_NotPion->GetBinError(ibinx);
    }
    TGraphErrors * g_PMSci_NotPion = new TGraphErrors(PDFPiCL_PMSci_NotPion->GetNbinsX(),x_PMSci_NotPion,y_PMSci_NotPion,ex_PMSci_NotPion,ey_PMSci_NotPion);
    g_PMSci_NotPion->Write("g_PMSci_NotPion");
    TSpline3 * s_PMSci_NotPion = new TSpline3("s_PMSci_NotPion",g_PMSci_NotPion);
    s_PMSci_NotPion->Write("s_PMSci_NotPion");
  }
  }
  else{
    TGraphErrors * g_WM_Muon[3];
    TSpline3 * s_WM_Muon[3];
    TGraphErrors * g_WM_Pion[3];
    TSpline3 * s_WM_Pion[3];
    TGraphErrors * g_WM_NotMuon[3];
    TSpline3 * s_WM_NotMuon[3];
    TGraphErrors * g_WM_NotPion[3];
    TSpline3 * s_WM_NotPion[3];

    for(int i=0;i<3;i++){
      PDFMuCL_WM_Muon[i]->Scale(1./PDFMuCL_WM_Muon[i]->Integral());
      PDFMuCL_WM_Muon[i]->Write();

      double x_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];double ex_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];double y_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];double ey_WM_Muon[PDFMuCL_WM_Muon[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFMuCL_WM_Muon[i]->GetNbinsX();ibinx++){
	x_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinCenter(ibinx);
	ex_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinWidth(ibinx)/2;
	y_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinContent(ibinx);
	ey_WM_Muon[ibinx-1]=PDFMuCL_WM_Muon[i]->GetBinError(ibinx);
      }
       g_WM_Muon[i] = new TGraphErrors(PDFMuCL_WM_Muon[i]->GetNbinsX(),x_WM_Muon,y_WM_Muon,ex_WM_Muon,ey_WM_Muon);
       g_WM_Muon[i]->Write(Form("g_WM_Muon_%d",i));
       s_WM_Muon[i] = new TSpline3(Form("s_WM_Muon_%d",i),g_WM_Muon[i]);
       s_WM_Muon[i]->Write(Form("s_WM_Muon_%d",i));

       PDFMuCL_WM_NotMuon[i]->Scale(1./PDFMuCL_WM_NotMuon[i]->Integral());
       PDFMuCL_WM_NotMuon[i]->Write();
       
       double x_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];double ex_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];double y_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];double ey_WM_NotMuon[PDFMuCL_WM_NotMuon[i]->GetNbinsX()];
       for(int ibinx=1;ibinx<=PDFMuCL_WM_NotMuon[i]->GetNbinsX();ibinx++){
	 x_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinCenter(ibinx);
	 ex_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinWidth(ibinx)/2;
	 y_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinContent(ibinx);
	 ey_WM_NotMuon[ibinx-1]=PDFMuCL_WM_NotMuon[i]->GetBinError(ibinx);
       }
       g_WM_NotMuon[i] = new TGraphErrors(PDFMuCL_WM_NotMuon[i]->GetNbinsX(),x_WM_NotMuon,y_WM_NotMuon,ex_WM_NotMuon,ey_WM_NotMuon);
       g_WM_NotMuon[i]->Write(Form("g_WM_NotMuon_%d",i));
       s_WM_NotMuon[i] = new TSpline3(Form("s_WM_NotMuon_%d",i),g_WM_NotMuon[i]);
       s_WM_NotMuon[i]->Write(Form("s_WM_NotMuon_%d",i));
    
       if(PiCL){
      PDFPiCL_WM_Pion[i]->Scale(1./PDFPiCL_WM_Pion[i]->Integral());
      PDFPiCL_WM_Pion[i]->Write();

      double x_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];double ex_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];double y_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];double ey_WM_Pion[PDFPiCL_WM_Pion[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFPiCL_WM_Pion[i]->GetNbinsX();ibinx++){
	x_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinCenter(ibinx);
	ex_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinWidth(ibinx)/2;
	y_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinContent(ibinx);
	ey_WM_Pion[ibinx-1]=PDFPiCL_WM_Pion[i]->GetBinError(ibinx);
      }
      g_WM_Pion[i] = new TGraphErrors(PDFPiCL_WM_Pion[i]->GetNbinsX(),x_WM_Pion,y_WM_Pion,ex_WM_Pion,ey_WM_Pion);
      g_WM_Pion[i]->Write(Form("g_WM_Pion_%d",i));
      s_WM_Pion[i] = new TSpline3(Form("s_WM_Pion_%d",i),g_WM_Pion[i]);
      s_WM_Pion[i]->Write(Form("s_WM_Pion_%d",i));
      PDFPiCL_WM_NotPion[i]->Scale(1./PDFPiCL_WM_NotPion[i]->Integral());
      PDFPiCL_WM_NotPion[i]->Write();

      double x_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];double ex_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];double y_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];double ey_WM_NotPion[PDFPiCL_WM_NotPion[i]->GetNbinsX()];
      for(int ibinx=1;ibinx<=PDFPiCL_WM_NotPion[i]->GetNbinsX();ibinx++){
	x_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinCenter(ibinx);
	ex_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinWidth(ibinx)/2;
	y_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinContent(ibinx);
	ey_WM_NotPion[ibinx-1]=PDFPiCL_WM_NotPion[i]->GetBinError(ibinx);
      }
       g_WM_NotPion[i] = new TGraphErrors(PDFPiCL_WM_NotPion[i]->GetNbinsX(),x_WM_NotPion,y_WM_NotPion,ex_WM_NotPion,ey_WM_NotPion);
       g_WM_NotPion[i]->Write(Form("g_WM_NotPion_%d",i));
       s_WM_NotPion[i] = new TSpline3(Form("s_WM_NotPion_%d",i),g_WM_NotPion[i]);
       s_WM_NotPion[i]->Write(Form("s_WM_NotPion_%d",i));
       }
    }
  }

  PDFMuCL_Ing_Muon->Scale(1./PDFMuCL_Ing_Muon->Integral()); 
  PDFMuCL_Ing_Muon->Write();
  PDFMuCL_Ing_NotMuon->Scale(1./PDFMuCL_Ing_NotMuon->Integral());
  PDFMuCL_Ing_NotMuon->Write();

  double x_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];double ex_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];double y_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];double ey_Ing_Muon[PDFMuCL_Ing_Muon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_Ing_Muon->GetNbinsX();ibinx++){
    x_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinCenter(ibinx);
    ex_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinWidth(ibinx)/2;
    y_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinContent(ibinx);
    ey_Ing_Muon[ibinx-1]=PDFMuCL_Ing_Muon->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_Muon = new TGraphErrors(PDFMuCL_Ing_Muon->GetNbinsX(),x_Ing_Muon,y_Ing_Muon,ex_Ing_Muon,ey_Ing_Muon);
  g_Ing_Muon->Write("g_Ing_Muon");
  TSpline3 * s_Ing_Muon = new TSpline3("s_Ing_Muon",g_Ing_Muon);
  s_Ing_Muon->Write("s_Ing_Muon");

    double x_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];double ex_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];double y_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];double ey_Ing_NotMuon[PDFMuCL_Ing_NotMuon->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFMuCL_Ing_NotMuon->GetNbinsX();ibinx++){
    x_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinCenter(ibinx);
    ex_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinWidth(ibinx)/2;
    y_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinContent(ibinx);
    ey_Ing_NotMuon[ibinx-1]=PDFMuCL_Ing_NotMuon->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_NotMuon = new TGraphErrors(PDFMuCL_Ing_NotMuon->GetNbinsX(),x_Ing_NotMuon,y_Ing_NotMuon,ex_Ing_NotMuon,ey_Ing_NotMuon);
  g_Ing_NotMuon->Write("g_Ing_NotMuon");
  TSpline3 * s_Ing_NotMuon = new TSpline3("s_Ing_NotMuon",g_Ing_NotMuon);
  s_Ing_NotMuon->Write("s_Ing_NotMuon");

  if(PiCL){
  PDFPiCL_Ing_Pion->Scale(1./PDFPiCL_Ing_Pion->Integral()); 
  PDFPiCL_Ing_Pion->Write();

  double x_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];double ex_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];double y_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];double ey_Ing_Pion[PDFPiCL_Ing_Pion->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFPiCL_Ing_Pion->GetNbinsX();ibinx++){
    x_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinCenter(ibinx);
    ex_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinWidth(ibinx)/2;
    y_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinContent(ibinx);
    ey_Ing_Pion[ibinx-1]=PDFPiCL_Ing_Pion->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_Pion = new TGraphErrors(PDFPiCL_Ing_Pion->GetNbinsX(),x_Ing_Pion,y_Ing_Pion,ex_Ing_Pion,ey_Ing_Pion);
  g_Ing_Pion->Write("g_Ing_Pion");
  TSpline3 * s_Ing_Pion = new TSpline3("s_Ing_Pion",g_Ing_Pion);
  s_Ing_Pion->Write("s_Ing_Pion");

    PDFPiCL_Ing_NotPion->Scale(1./PDFPiCL_Ing_NotPion->Integral());
  PDFPiCL_Ing_NotPion->Write();
  
  double x_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];double ex_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];double y_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];double ey_Ing_NotPion[PDFPiCL_Ing_NotPion->GetNbinsX()];
  for(int ibinx=1;ibinx<=PDFPiCL_Ing_NotPion->GetNbinsX();ibinx++){
    x_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinCenter(ibinx);
    ex_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinWidth(ibinx)/2;
    y_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinContent(ibinx);
    ey_Ing_NotPion[ibinx-1]=PDFPiCL_Ing_NotPion->GetBinError(ibinx);
  }
  TGraphErrors * g_Ing_NotPion = new TGraphErrors(PDFPiCL_Ing_NotPion->GetNbinsX(),x_Ing_NotPion,y_Ing_NotPion,ex_Ing_NotPion,ey_Ing_NotPion);
  g_Ing_NotPion->Write("g_Ing_NotPion");
  TSpline3 * s_Ing_NotPion = new TSpline3("s_Ing_NotPion",g_Ing_NotPion);
  s_Ing_NotPion->Write("s_Ing_NotPion");
  }

  PDFParticle->Scale(1./PDFParticle->Integral());
  PDFParticle->Write();
  
  wfile->Write();
  wfile->Close();

}

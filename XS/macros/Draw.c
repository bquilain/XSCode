{
  TFile * fMC = new TFile("../plots/MCPlots_Plan.root");
  TH1D * MC_MuCL = (TH1D*) fMC->Get("Stack_MuCL");
  /*
  TFile * fData = new TFile("../plots/DataPlots.root");
  TH1D * Data_MuCL = (TH1D*) fData->Get("Data_MuCL");

  int NMCFiles=1000;
  cout<<"please enter the number of MC file processed to obtain the input root files:(default=1000)"<<endl;
  cin>>NMCFiles;

  double NPOTs=1e21
  cout<<"please enter the number of POT processed in the input root data files as ***e21:(default=1e21)"<<endl;
  cin>>NPOTs;
  NPOTs=NPOTs*pow(10,21);

  double ScaleFactor=(NPOTs/1e21)/NMCFiles;

  TCanvas * cMuCL = new TCanvas("c0");
  Data_MuCL->Draw();
  Data_MuCL->Draw("same");

  */
  TH1D * MuCL_TrueMuon = (TH1D*) fMC->Get("hMuCL_TrueMuon");
  TH1D * MuCL_TruePion = (TH1D*) fMC->Get("hMuCL_TruePion");
  TH1D * MuCL_TrueProton = (TH1D*) fMC->Get("hMuCL_TrueProton");
  MuCL_TrueMuon->Sumw2();
  MuCL_TruePion->Sumw2();
  MuCL_TrueProton->Sumw2();
  
  MuCL_TrueMuon->Scale(1./MuCL_TrueMuon->Integral());
  MuCL_TrueProton->Scale(1./MuCL_TrueProton->Integral());
  MuCL_TruePion->Scale(1./MuCL_TruePion->Integral());

  TCanvas * cCL = new TCanvas();
  MuCL_TrueMuon->SetLineColor(kBlue);
  MuCL_TruePion->SetLineColor(kGreen);
  MuCL_TrueProton->SetLineColor(kRed);
  MuCL_TrueMuon->Draw();
  MuCL_TruePion->Draw("same");
  MuCL_TrueProton->Draw("same");

  
  TCanvas * cHighCL = new TCanvas("cHighCL","Optimisation of the high confidence level");//goal is to remove most of the protons, in the case that there is a proton and a muon -> normalise the distribution in the same way, to show probability to obtain a muon or a proton
  TH1D * Muon_Purity = (TH1D*) MuCL_TrueMuon->Clone("Muon_Purity");//Purity of muon w/ confidence level, assuming same amount of protons and muons
  TH1D * Muon_Efficiency = (TH1D*) MuCL_TrueMuon->Clone("Muon_Efficiency");//Efficiency of muon w/ confidence level
  TH1D * Proton_Purity = (TH1D*) MuCL_TrueProton->Clone("Proton_Purity");//Purity of muon w/ confidence level, assuming same amount of protons and muons
  TH1D * Proton_Efficiency = (TH1D*) MuCL_TrueProton->Clone("Proton_Efficiency");//Efficiency of muon w/ confidence level
  TH1D * Pion_Purity = (TH1D*) MuCL_TruePion->Clone("Pion_Purity");//Purity of muon w/ confidence level, assuming same amount of protons and muons
  TH1D * Pion_Efficiency = (TH1D*) MuCL_TruePion->Clone("Pion_Efficiency");//Efficiency of muon w/ confidence level

  double NMuTotal=MuCL_TrueMuon->Integral(1,Muon_Efficiency->GetNbinsX());
  for(int ibinx=1;ibinx<=Muon_Efficiency->GetNbinsX();ibinx++){
    double NMu=MuCL_TrueMuon->Integral(ibinx,Muon_Efficiency->GetNbinsX());
    double NP=MuCL_TrueProton->Integral(ibinx,Proton_Efficiency->GetNbinsX());
    double EffMu=NMu;
    double PurMu=NMu;
    if(NMuTotal!=0) EffMu/=NMuTotal;
    if((NMu+NP)!=0) PurMu/=(NMu+NP);
    Muon_Efficiency->SetBinContent(ibinx,PurMu);
    Muon_Purity->SetBinContent(ibinx,EffMu);
  }
  //TH1D * MuonAndProton = (TH1D*) MuCL_TrueMuon->Clone("MuonAndProton");
  //MuonAndProton->Add(MuCL_TrueProton);
  //Muon_Purity->Divide(MuonAndProton);
  Muon_Purity->Draw();
  Muon_Efficiency->Draw("same");

  
  TCanvas * cLowCL = new TCanvas("cLowCL","Optimisation of the low confidence level");//goal is to remove most of the pions, in the case that there is a proton and a pion -> normalise the distribution in the same way, to show probability to obtain a pion or a proton

  double NPTotal=MuCL_TrueProton->Integral(1,Proton_Efficiency->GetNbinsX());
  for(int ibinx=1;ibinx<=Proton_Efficiency->GetNbinsX();ibinx++){
    double NP=MuCL_TrueProton->Integral(1,ibinx);
    double NPi=MuCL_TruePion->Integral(1,ibinx);
    double EffMu=NP; double PurMu=NP;
    if(NPTotal!=0) EffMu/=NPTotal;
    if((NP+NPi)!=0) PurMu/=(NP+NPi);
    Proton_Efficiency->SetBinContent(ibinx,PurMu);
    Proton_Purity->SetBinContent(ibinx,EffMu);
  }
  //TH1D * ProtonAndPion = (TH1D*) MuCL_TrueProton->Clone("ProtonAndPion");
  //ProtonAndPion->Add(MuCL_TruePion);
  //Proton_Purity->Divide(ProtonAndPion);
  Proton_Purity->Draw();
  Proton_Efficiency->Draw("same");

  int NFSIs=11;
  char Name[256];
  TH1D * hMuCL[NFSIs];
  TH1D * hMuCL_Lowest[NFSIs];

  TH1D * hMuCL_CC0pi = (TH1D*) fMC->Get("hMuCL0");
  TH1D * hMuCL_Lowest_CC0pi = (TH1D*) fMC->Get("hMuCL_Lowest0");
  TH1D * hMuCL_All = (TH1D*) hMuCL_CC0pi->Clone("hMuCL_All");
  TH1D * hMuCL_Lowest_All = (TH1D*) hMuCL_Lowest_CC0pi->Clone("hMuCL_Lowest_All");
  hMuCL_CC0pi->SetFillStyle(0);
  hMuCL_Lowest_CC0pi->SetFillStyle(0);
  hMuCL_All->SetFillStyle(0);
  hMuCL_Lowest_All->SetFillStyle(0);
 
  for(int fsi=0;fsi<NFSIs;fsi++){
    sprintf(Name,"hMuCL%d",fsi);
    hMuCL[fsi] = (TH1D*) fMC->Get(Name);
    sprintf(Name,"hMuCL_Lowest%d",fsi);
    hMuCL_Lowest[fsi] = (TH1D*) fMC->Get(Name);
    hMuCL[fsi]->SetFillStyle(0);
    hMuCL_Lowest[fsi]->SetFillStyle(0);
    
    if(fsi>0){
      if(fsi<3){
	hMuCL_CC0pi->Add(hMuCL[fsi]);
	hMuCL_Lowest_CC0pi->Add(hMuCL_Lowest[fsi]);
      }
      hMuCL_All->Add(hMuCL[fsi]);
      hMuCL_Lowest_All->Add(hMuCL_Lowest[fsi]);
    } 
  }

  TH1D * CC0pi_Purity_High = (TH1D*) hMuCL_CC0pi->Clone("CC0pi_Purity_High");
  TH1D * CC0pi_Efficiency_High = (TH1D*) hMuCL_CC0pi->Clone("CC0pi_Efficiency_High");

  TH1D * CC0pi_Purity_Low = (TH1D*) hMuCL_Lowest_CC0pi->Clone("CC0pi_Purity_Low");
  TH1D * CC0pi_Efficiency_Low = (TH1D*) hMuCL_Lowest_CC0pi->Clone("CC0pi_Efficiency_Low");

  double NCC0pi_Total_High=hMuCL_CC0pi->Integral(1,CC0pi_Efficiency_High->GetNbinsX());
  double NCC0pi_Total_Low=hMuCL_Lowest_CC0pi->Integral(1,CC0pi_Efficiency_Low->GetNbinsX());
  
  for(int ibinx=1;ibinx<=CC0pi_Efficiency_High->GetNbinsX();ibinx++){
    
    double NCC0pi_High=hMuCL_CC0pi->Integral(ibinx,CC0pi_Efficiency_High->GetNbinsX());
    double NAll_High=hMuCL_All->Integral(ibinx,CC0pi_Efficiency_High->GetNbinsX());
   
    double EffCC0pi_High=NCC0pi_High;
    double PurCC0pi_High=NCC0pi_High;
    if(NCC0pi_Total_High!=0) EffCC0pi_High/=NCC0pi_Total_High;
    if(NAll_High!=0) PurCC0pi_High/=NAll_High;

    CC0pi_Purity_High->SetBinContent(ibinx,PurCC0pi_High);
    CC0pi_Efficiency_High->SetBinContent(ibinx,EffCC0pi_High);

    double NCC0pi_Low=hMuCL_Lowest_CC0pi->Integral(1,ibinx);
    double NAll_Low=hMuCL_Lowest_All->Integral(1,ibinx);    

    double EffCC0pi_Low=NCC0pi_Low;
    double PurCC0pi_Low=NCC0pi_Low;
    if(NCC0pi_Total_Low!=0) EffCC0pi_Low/=NCC0pi_Total_Low;
    if(NAll_Low!=0) PurCC0pi_Low/=NAll_Low;

    CC0pi_Purity_Low->SetBinContent(ibinx,PurCC0pi_Low);
    CC0pi_Efficiency_Low->SetBinContent(ibinx,EffCC0pi_Low);

  }

  TCanvas * cEff_High = new TCanvas("cEff_High","Efficiency of CC0pi with MuCL cut value");
  CC0pi_Efficiency_High->Draw();
  CC0pi_Purity_High->Draw("same");

  TCanvas * cEff_Low = new TCanvas("cEff_Low","Efficiency of CC0pi with MuCL cut value");
  CC0pi_Efficiency_Low->Draw();
  CC0pi_Purity_Low->Draw("same");


}

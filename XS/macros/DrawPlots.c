//Macro helpful for efficiency / purity optimization
{
  gROOT->SetBatch(kTRUE);//Silent mode
  
  //TFile * fMC = new TFile("../files/MCSelected_PM_PG_Systematics0_0.root");
  TFile * fMC = new TFile("../files/MCSelected_PM_Systematics0_0.root");
  //TFile * fMC = new TFile("../files/MCSelected_Systematics0_0_clmuonplan.root");
  //TFile * fMC = new TFile("../plots/MCPlots_20160119_EfficiencyAdded.root");

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
  /////////////////////////////////////////////DRAW VARIABLES FOR INDIVIDUAL PARTICLES//////////////////////////////////////////////
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
  //MuCL_TrueMuon->SetFillColor(kBlue);
  MuCL_TruePion->SetLineColor(kGreen);
  //MuCL_TruePion->SetFillColor(kGreen);
  MuCL_TrueProton->SetLineColor(kRed);
  //MuCL_TrueProton->SetFillColor(kRed);
  MuCL_TrueMuon->Draw("");
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
  Muon_Efficiency->SetLineColor(kBlack);
  Muon_Purity->SetLineColor(kBlue);
  Pion_Purity->SetLineColor(kGreen);
  Proton_Purity->SetLineColor(kRed);
  Muon_Efficiency->Draw();
  Muon_Purity->Draw("same");
  Pion_Purity->Draw("same");
  Proton_Purity->Draw("same");
  TLegend * lpart = new TLegend(0.7,0.75,0.95,0.99);
  lpart->AddEntry(Muon_Efficiency,"#mu efficiency");
  lpart->AddEntry(Muon_Purity,"#mu purity");
  lpart->AddEntry(Pion_Purity,"#pi purity");
  lpart->AddEntry(Proton_Purity,"p purity");
  lpart->Draw("same");
  Muon_Efficiency->GetYaxis()->SetRangeUser(0,1.2);
  cHighCL->SetLogx();
  
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


  /////////////////////////////////////////////DRAW VARIABLES FOR INTERACTIONS//////////////////////////////////////////////
  int NFSIs=11;
  char Name[256];
#define CLMuon
#ifdef CLMuon
  int RebinValueCL=15;

  TH1D * hMuCL[NFSIs];
  TH2D * hMuCL_2tracks[NFSIs];
  TH1D * hMuCL_CC0pi = (TH1D*) fMC->Get("hMuCL0");
  TH2D * hMuCL_2tracks_CC0pi = (TH2D*) fMC->Get("hMuCL_2tracks0");
  TH1D * hMuCL_All = (TH1D*) hMuCL_CC0pi->Clone("hMuCL_All");
  TH2D * hMuCL_2tracks_All = (TH2D*) hMuCL_2tracks_CC0pi->Clone("hMuCL_2tracks_All");
  hMuCL_2tracks_CC0pi->Rebin2D(RebinValueCL,RebinValueCL);
  hMuCL_2tracks_All->Rebin2D(RebinValueCL,RebinValueCL);
  hMuCL_CC0pi->SetFillStyle(0);
  hMuCL_2tracks_CC0pi->SetFillStyle(0);
  hMuCL_All->SetFillStyle(0);
  hMuCL_2tracks_All->SetFillStyle(0);


  for(int fsi=1;fsi<NFSIs;fsi++){
    sprintf(Name,"hMuCL%d",fsi);
    hMuCL[fsi] = (TH1D*) fMC->Get(Name);
    hMuCL[fsi]->SetFillStyle(0);

    sprintf(Name,"hMuCL_2tracks%d",fsi);
    hMuCL_2tracks[fsi] = (TH2D*) fMC->Get(Name);
    hMuCL_2tracks[fsi]->SetFillStyle(0);
    hMuCL_2tracks[fsi]->Rebin2D(RebinValueCL,RebinValueCL);

    if(fsi>-1){
      if(fsi<3){
	hMuCL_CC0pi->Add(hMuCL[fsi]);
	hMuCL_2tracks_CC0pi->Add(hMuCL_2tracks[fsi]);
      }
      hMuCL_All->Add(hMuCL[fsi]);
      hMuCL_2tracks_All->Add(hMuCL_2tracks[fsi]);
    } 
  }


    
  TH1D * CC0pi_Purity_1track = (TH1D*) hMuCL_CC0pi->Clone("CC0pi_Purity_1track");
  TH1D * CC0pi_Efficiency_1track = (TH1D*) hMuCL_CC0pi->Clone("CC0pi_Efficiency_1track");

  TH2D * CC0pi_Purity_2tracks = (TH2D*) hMuCL_2tracks_CC0pi->Clone("CC0pi_Purity_2tracks");
  TH2D * CC0pi_Efficiency_2tracks = (TH2D*) hMuCL_2tracks_CC0pi->Clone("CC0pi_Efficiency_2tracks");

  double NBinsX=CC0pi_Efficiency_1track->GetNbinsX();
  double NBinsY=CC0pi_Efficiency_2tracks->GetNbinsY();
  double NCC0pi_Total_1track=hMuCL_CC0pi->Integral(1,NBinsX);
  double NCC0pi_Total_2tracks=hMuCL_2tracks_CC0pi->Integral(1,NBinsX,1,NBinsY);

  
  for(int ibinx=1;ibinx<=NBinsX;ibinx++){
      cout<<ibinx<<endl;    
    double NCC0pi_1track=hMuCL_CC0pi->Integral(ibinx,NBinsX);
    double NAll_1track=hMuCL_All->Integral(ibinx,NBinsX);
   
    double EffCC0pi_1track=NCC0pi_1track;
    double PurCC0pi_1track=NCC0pi_1track;
    if(NCC0pi_Total_1track!=0) EffCC0pi_1track/=NCC0pi_Total_1track;
    if(NAll_1track!=0) PurCC0pi_1track/=NAll_1track;

    CC0pi_Purity_1track->SetBinContent(ibinx,PurCC0pi_1track);
    CC0pi_Efficiency_1track->SetBinContent(ibinx,EffCC0pi_1track);

    for(int ibiny=1;ibiny<=NBinsY;ibiny++){
      double NCC0pi_2tracks=hMuCL_2tracks_CC0pi->Integral(ibinx,NBinsX,1,ibiny);
      double NAll_2tracks=hMuCL_2tracks_All->Integral(ibinx,NBinsX,1,ibiny);    
      
      double EffCC0pi_2tracks=NCC0pi_2tracks;
      double PurCC0pi_2tracks=NCC0pi_2tracks;

      if(NCC0pi_Total_2tracks!=0) EffCC0pi_2tracks/=NCC0pi_Total_2tracks;
      if(NAll_2tracks!=0) PurCC0pi_2tracks/=NAll_2tracks;
      //cout<<EffCC0pi_2tracks<<", "<<PurCC0pi_2tracks<<endl;
      CC0pi_Purity_2tracks->SetBinContent(ibinx,ibiny,PurCC0pi_2tracks);
      CC0pi_Efficiency_2tracks->SetBinContent(ibinx,ibiny,EffCC0pi_2tracks);
    }
    }

  TCanvas * cEff_1track = new TCanvas("cEff_1track","Tune mu-like cut: Efficiency/Purity of CC0pi with MuCL cut value");//Use the 1 track sample
  CC0pi_Efficiency_1track->SetLineColor(kRed);
  CC0pi_Efficiency_1track->Draw();
  CC0pi_Purity_1track->SetLineColor(kBlue);
  CC0pi_Purity_1track->Draw("same");
  cEff_1track->SaveAs("../plots/MVAOptimisation/CC0pi_1track_EfficiencyPurity.eps");
  
  gStyle->SetPaintTextFormat("2.2g");//Ste the text format in the option draw("colztext")
  TCanvas * cEff_2tracks = new TCanvas("cEff_2tracks","Tune p-like cut: Efficiency of CC0pi with MuCL cut value");
  //CC0pi_Efficiency_2tracks->Rebin2D(RebinValueCL,1);
  //CC0pi_Efficiency_2tracks->Scale(1./RebinValueCL);
  CC0pi_Efficiency_2tracks->Draw("colztext");
  CC0pi_Efficiency_2tracks->GetXaxis()->SetRangeUser(0,1);
  CC0pi_Efficiency_2tracks->GetYaxis()->SetRangeUser(0,1);
  cEff_2tracks->SaveAs("../plots/MVAOptimisation/CC0pi_2tracks_Efficiency.eps");
  
  TCanvas * cPur_2tracks = new TCanvas("cPur_2tracks","Tune p-like cut: Purity of CC0pi with MuCL cut value");
  //CC0pi_Purity_2tracks->Rebin2D(RebinValueCL,1);
  //CC0pi_Purity_2tracks->Scale(1./RebinValueCL);
  CC0pi_Purity_2tracks->Draw("colztext");
  CC0pi_Purity_2tracks->GetXaxis()->SetRangeUser(0,1);
  CC0pi_Purity_2tracks->GetYaxis()->SetRangeUser(0,1);
  cPur_2tracks->SaveAs("../plots/MVAOptimisation/CC0pi_2tracks_Purity.eps");
  
#endif
  
#define MVA  
#ifdef MVA
  int RebinValueMVA=5;

  TH1D * hMVAMuondiscriminant_1track[NFSIs];
  TH2D * hMVAMuonVSProtondiscriminant_2tracks[NFSIs];
  TH1D * hMVAMuondiscriminant_1track_CC0pi = (TH1D*) fMC->Get("hMVAMuondiscriminant_1track0");
  TH2D * hMVAMuonVSProtondiscriminant_2tracks_CC0pi = (TH2D*) fMC->Get("hMVAMuonVSProtondiscriminant_2tracks0");
  TH1D * hMVAMuondiscriminant_1track_All = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("hMVAdiscriminant_1track_All");
  TH2D * hMVAMuonVSProtondiscriminant_2tracks_All = (TH2D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("hMVAMuonVSProtondiscriminant_2tracks_All");

  hMVAMuondiscriminant_1track_CC0pi->SetFillStyle(0);
  hMVAMuonVSProtondiscriminant_2tracks_CC0pi->SetFillStyle(0);
  hMVAMuondiscriminant_1track_All->SetFillStyle(0);
  hMVAMuonVSProtondiscriminant_2tracks_All->SetFillStyle(0);
  
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track[NFSIs];
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi = (TH2D*) fMC->Get("hMVAMuondiscriminantVSMuonMomentum_1track0");
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track_All = (TH2D*) hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->Clone("hMVAMuondiscriminantVSMuonMomentum_1track_All");
  
  hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->SetFillStyle(0);
  hMVAMuondiscriminantVSMuonMomentum_1track_All->SetFillStyle(0);
  
  for(int fsi=1;fsi<NFSIs;fsi++){
    sprintf(Name,"hMVAMuondiscriminant_1track%d",fsi);
    hMVAMuondiscriminant_1track[fsi] = (TH1D*) fMC->Get(Name);
    hMVAMuondiscriminant_1track[fsi]->SetFillStyle(0);

    sprintf(Name,"hMVAMuonVSProtondiscriminant_2tracks%d",fsi);
    hMVAMuonVSProtondiscriminant_2tracks[fsi] = (TH2D*) fMC->Get(Name);
    hMVAMuonVSProtondiscriminant_2tracks[fsi]->SetFillStyle(0);

    sprintf(Name,"hMVAMuondiscriminantVSMuonMomentum_1track%d",fsi);
    hMVAMuondiscriminantVSMuonMomentum_1track[fsi] = (TH2D*) fMC->Get(Name);
    hMVAMuondiscriminantVSMuonMomentum_1track[fsi]->SetFillStyle(0);

    if(fsi>-1){
      if(fsi<3){
	hMVAMuondiscriminant_1track_CC0pi->Add(hMVAMuondiscriminant_1track[fsi]);
	hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Add(hMVAMuonVSProtondiscriminant_2tracks[fsi]); 
	hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->Add(hMVAMuondiscriminantVSMuonMomentum_1track[fsi]);
     }
      hMVAMuondiscriminant_1track_All->Add(hMVAMuondiscriminant_1track[fsi]);
      hMVAMuonVSProtondiscriminant_2tracks_All->Add(hMVAMuonVSProtondiscriminant_2tracks[fsi]);
      hMVAMuondiscriminantVSMuonMomentum_1track_All->Add(hMVAMuondiscriminantVSMuonMomentum_1track[fsi]);
    } 
  }

  hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Rebin2D(RebinValueMVA,RebinValueMVA);  
  hMVAMuonVSProtondiscriminant_2tracks_All->Rebin2D(RebinValueMVA,RebinValueMVA);  
  
  TH1D * mvaCC0pi_Purity_1track = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("mvaCC0pi_Purity_1track");
  TH1D * mvaCC0pi_Efficiency_1track = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("mvaCC0pi_Efficiency_1track");
  TH1D * mvaCC0pi_Purity_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("mvaCC0pi_Purity_2tracks");
  TH1D * mvaCC0pi_Efficiency_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("mvaCC0pi_Efficiency_2tracks");

  double mvaNBinsX=mvaCC0pi_Efficiency_1track->GetNbinsX();
  double mvaNBinsY=mvaCC0pi_Efficiency_1track->GetNbinsX();
  double mvaNCC0pi_Total_1track=hMVAMuondiscriminant_1track_CC0pi->Integral(1,mvaNBinsX);
  double mvaNCC0pi_Total_2tracks=hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Integral(1,mvaNBinsX,1,mvaNBinsY);
  
  for(int ibinx=1;ibinx<=mvaNBinsX;ibinx++){
      cout<<"MVA, "<<ibinx<<endl;    
      double NCC0pi_1track=hMVAMuondiscriminant_1track_CC0pi->Integral(ibinx,mvaNBinsX);
      double NAll_1track=hMVAMuondiscriminant_1track_All->Integral(ibinx,mvaNBinsX);
      double EffCC0pi_1track=NCC0pi_1track;
      double PurCC0pi_1track=NCC0pi_1track;
      if(mvaNCC0pi_Total_1track!=0) EffCC0pi_1track/=mvaNCC0pi_Total_1track;
      if(NAll_1track!=0) PurCC0pi_1track/=NAll_1track;
      mvaCC0pi_Purity_1track->SetBinContent(ibinx,PurCC0pi_1track);
      mvaCC0pi_Efficiency_1track->SetBinContent(ibinx,EffCC0pi_1track);

      for(int ibiny=1;ibiny<=mvaNBinsY;ibiny++){
	double NCC0pi_2tracks=hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Integral(ibinx,mvaNBinsX,ibiny,mvaNBinsY);
	double NAll_2tracks=hMVAMuonVSProtondiscriminant_2tracks_All->Integral(ibinx,mvaNBinsX,ibiny,mvaNBinsY);
	
	double EffCC0pi_2tracks=NCC0pi_2tracks;
	double PurCC0pi_2tracks=NCC0pi_2tracks;

	if(mvaNCC0pi_Total_2tracks!=0) EffCC0pi_2tracks/=mvaNCC0pi_Total_2tracks;
	if(NAll_2tracks!=0) PurCC0pi_2tracks/=NAll_2tracks;
	mvaCC0pi_Purity_2tracks->SetBinContent(ibinx,ibiny,PurCC0pi_2tracks);
	mvaCC0pi_Efficiency_2tracks->SetBinContent(ibinx,ibiny,EffCC0pi_2tracks);
      }
  }

  TCanvas * cmvaEfficiency_1track = new TCanvas("cmvaEfficiency_1track","1 track,Efficiency/Purity of CC0pi with MVA discriminant cut value");//Use the 1 track sample
  mvaCC0pi_Efficiency_1track->SetLineColor(kRed);
  mvaCC0pi_Purity_1track->SetLineColor(kBlue);
  mvaCC0pi_Efficiency_1track->Draw();
  mvaCC0pi_Purity_1track->Draw("same");
  cmvaEfficiency_1track->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_1track_EfficiencyPurity.eps");
  
  TCanvas * cmvaEfficiency_2tracks = new TCanvas("cmvaEfficiency_2tracks","2 tracks,Efficiency of CC0pi with MVA discriminant cut value");//Use the 2 tracks sample
  //mvaCC0pi_Efficiency_2tracks->SetLineColor(kRed);
  //mvaCC0pi_Purity_2tracks->SetLineColor(kBlue);
  mvaCC0pi_Efficiency_2tracks->Draw("COLZTEXT");
  cmvaEfficiency_2tracks->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_Efficiency.eps");

  TCanvas * cmvaPurity_2tracks = new TCanvas("cmvaPurity_2tracks","2 tracks,Purity of CC0pi with MVA discriminant cut value");//Use the 2 tracks sample
  mvaCC0pi_Purity_2tracks->Draw("COLZTEXT");
  cmvaPurity_2tracks->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_Purity.eps");

  TCanvas * cmvaVSMuonMomentum_1track = new TCanvas("cmvaVSMuonMomentum_1track","1 track,MVA discriminant cut value with true muon momentum");//Use the 1 track sample
  cmvaVSMuonMomentum_1track->Divide(1,2);
  cmvaVSMuonMomentum_1track->cd(1);
  hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->Draw("COLZ");
  cmvaVSMuonMomentum_1track->cd(2);
  hMVAMuondiscriminantVSMuonMomentum_1track_All->Draw("COLZ");
  
#endif


  TH1D * TrueParticleType_1track[NFSIs];
  TH2D * TrueParticleType_2tracks[NFSIs];
  TH1D * TrueParticleType_1track_CC0pi = (TH1D*) fMC->Get("TrueParticleType_1track0");
  TH2D * TrueParticleType_2tracks_CC0pi = (TH2D*) fMC->Get("TrueParticleType_2tracks0");
  TH1D * TrueParticleType_1track_All = (TH1D*) TrueParticleType_1track_CC0pi->Clone("TrueParticleType_1track_All");
  TH2D * TrueParticleType_2tracks_All = (TH2D*) TrueParticleType_2tracks_CC0pi->Clone("TrueParticleType_2tracks_All");
  TrueParticleType_1track_CC0pi->SetFillStyle(0);
  TrueParticleType_2tracks_CC0pi->SetFillStyle(0);
  TrueParticleType_1track_All->SetFillStyle(0);
  TrueParticleType_2tracks_All->SetFillStyle(0);

  
  for(int fsi=1;fsi<NFSIs;fsi++){
    sprintf(Name,"TrueParticleType_1track%d",fsi);
    TrueParticleType_1track[fsi] = (TH1D*) fMC->Get(Name);
    TrueParticleType_1track[fsi]->SetFillStyle(0);

    sprintf(Name,"TrueParticleType_2tracks%d",fsi);
    TrueParticleType_2tracks[fsi] = (TH2D*) fMC->Get(Name);
    TrueParticleType_2tracks[fsi]->SetFillStyle(0);

    if(fsi>-1){
      if(fsi<3){
	TrueParticleType_1track_CC0pi->Add(TrueParticleType_1track[fsi]);
	TrueParticleType_2tracks_CC0pi->Add(TrueParticleType_2tracks[fsi]);
      }
      TrueParticleType_1track_All->Add(TrueParticleType_1track[fsi]);
      TrueParticleType_2tracks_All->Add(TrueParticleType_2tracks[fsi]);
    } 
  }
  
  
  //1 track sample: good sample is assumed to be only true muon sample:
  double dCC0pi_Purity_1track = TrueParticleType_1track_CC0pi->GetBinContent(1)/TrueParticleType_1track_All->GetBinContent(1);
  double dCC0pi_Efficiency_1track = TrueParticleType_1track_CC0pi->GetBinContent(1)/TrueParticleType_1track_CC0pi->Integral();

  //2 tracks sample: good sample is assumed to be true muon + true proton sample:
  double dCC0pi_Purity_2tracks = (TrueParticleType_2tracks_CC0pi->GetBinContent(1,3)+TrueParticleType_2tracks_CC0pi->GetBinContent(3,1))/(TrueParticleType_2tracks_All->GetBinContent(1,3)+TrueParticleType_2tracks_All->GetBinContent(3,1));
  double dCC0pi_Efficiency_2tracks = (TrueParticleType_2tracks_CC0pi->GetBinContent(1,3)+TrueParticleType_2tracks_CC0pi->GetBinContent(3,1))/TrueParticleType_2tracks_CC0pi->Integral();

  cout<<setprecision(2)<<"1 Track sample:"<<endl;
  cout<<"Maximal purity = "<<dCC0pi_Purity_1track<<", reachable with an efficiency = "<<dCC0pi_Efficiency_1track<<endl;
  cout<<"2 Tracks sample:"<<endl;
  cout<<"Maximal purity = "<<dCC0pi_Purity_2tracks<<", reachable with an efficiency = "<<dCC0pi_Efficiency_2tracks<<endl; 
  cout<<" & \\mu & \\pi & p & Others \\\\"<<endl;
 cout<<"CC0\\pi 1 track & "<<TrueParticleType_1track_CC0pi->GetBinContent(1)/TrueParticleType_1track_CC0pi->Integral()<<" & "<<TrueParticleType_1track_CC0pi->GetBinContent(2)/TrueParticleType_1track_CC0pi->Integral()<<" & "<<TrueParticleType_1track_CC0pi->GetBinContent(3)/TrueParticleType_1track_CC0pi->Integral()<<" & "<<TrueParticleType_1track_CC0pi->GetBinContent(4)/TrueParticleType_1track_CC0pi->Integral()<<" \\ "<<endl;
  cout<<"CC0\\pi 2 tracks & \\mu & \\pi & p & Others \\\\"<<endl;
 cout<<"\\mu ";
 for(int i=1;i<5;i++) cout<<" & "<<TrueParticleType_2tracks_CC0pi->GetBinContent(1,i)/TrueParticleType_2tracks_CC0pi->Integral();
 cout<<" \\\\"<<endl;
 cout<<"\\pi ";
 for(int i=1;i<5;i++) cout<<" & "<<TrueParticleType_2tracks_CC0pi->GetBinContent(2,i)/TrueParticleType_2tracks_CC0pi->Integral();
 cout<<" \\\\"<<endl;
 cout<<"p ";
 for(int i=1;i<5;i++) cout<<" & "<<TrueParticleType_2tracks_CC0pi->GetBinContent(3,i)/TrueParticleType_2tracks_CC0pi->Integral();
 cout<<" \\\\"<<endl;
 cout<<"Others ";
 for(int i=1;i<5;i++) cout<<" & "<<TrueParticleType_2tracks_CC0pi->GetBinContent(4,i)/TrueParticleType_2tracks_CC0pi->Integral();
 cout<<" \\\\"<<endl;

#define DEDXPLOTS
#ifdef DEDXPLOTS

 int NSamples=6;
 TH1D * hEnergyDepositionLength_Muon_1D[NSamples];
 TH1D * hEnergyDepositionLength_Pion_1D[NSamples];
 TH1D * hEnergyDepositionLength_Proton_1D[NSamples];
 TH1D * hEnergyDepositionSplineLength_Muon_1D[NSamples];
 TH1D * hEnergyDepositionSplineLength_Pion_1D[NSamples];
 TH1D * hEnergyDepositionSplineLength_Proton_1D[NSamples];

 //Ratio plots
 TCanvas * cRatio_EnergyDepositionLength[NSamples];
 //Ratio
 TH1D * RatioPionMuon_EnergyDepositionLength[NSamples];
 TH1D * RatioProtonMuon_EnergyDepositionLength[NSamples];
 TH1D * RatioPionMuon_EnergyDepositionSplineLength[NSamples];
 TH1D * RatioProtonMuon_EnergyDepositionSplineLength[NSamples];
 //(Ratio -1)/Error
 TH1D * NumberOfSigmaPionMuon_EnergyDepositionLength[NSamples];
 TH1D * NumberOfSigmaProtonMuon_EnergyDepositionLength[NSamples];
 TH1D * NumberOfSigmaPionMuon_EnergyDepositionSplineLength[NSamples];
 TH1D * NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[NSamples];
 

 for(int is=0;is<NSamples;is++){
   hEnergyDepositionLength_Muon_1D[is] = (TH1D*) fMC->Get(Form("EnergyDepositionLength_Muon_%d_y",is));
   hEnergyDepositionLength_Pion_1D[is] = (TH1D*) fMC->Get(Form("EnergyDepositionLength_Pion_%d_y",is));
   hEnergyDepositionLength_Proton_1D[is] = (TH1D*) fMC->Get(Form("EnergyDepositionLength_Proton_%d_y",is));

   hEnergyDepositionSplineLength_Muon_1D[is] = (TH1D*) fMC->Get(Form("EnergyDepositionSplineLength_Muon_%d_y",is));
   hEnergyDepositionSplineLength_Pion_1D[is] = (TH1D*) fMC->Get(Form("EnergyDepositionSplineLength_Pion_%d_y",is));
   hEnergyDepositionSplineLength_Proton_1D[is] = (TH1D*) fMC->Get(Form("EnergyDepositionSplineLength_Proton_%d_y",is));
   
   RatioPionMuon_EnergyDepositionLength[is] = new TH1D(Form("RatioPionMuon_EnergyDepositionLength_%d",is),"",hEnergyDepositionLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioPionMuon_EnergyDepositionLength[is]->Sumw2();
   
   RatioProtonMuon_EnergyDepositionLength[is] = new TH1D(Form("RatioProtonMuon_EnergyDepositionLength_%d",is),"",hEnergyDepositionLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioProtonMuon_EnergyDepositionLength[is]->Sumw2();

   RatioPionMuon_EnergyDepositionSplineLength[is] = new TH1D(Form("RatioPionMuon_EnergyDepositionSplineLength_%d",is),"",hEnergyDepositionSplineLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioPionMuon_EnergyDepositionSplineLength[is]->Sumw2();
   
   RatioProtonMuon_EnergyDepositionSplineLength[is] = new TH1D(Form("RatioProtonMuon_EnergyDepositionSplineLength_%d",is),"",hEnergyDepositionSplineLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioProtonMuon_EnergyDepositionSplineLength[is]->Sumw2();

   //
   NumberOfSigmaPionMuon_EnergyDepositionLength[is] = new TH1D(Form("NumberOfSigmaPionMuon_EnergyDepositionLength_%d",is),"",hEnergyDepositionLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->Sumw2();
   
   NumberOfSigmaProtonMuon_EnergyDepositionLength[is] = new TH1D(Form("NumberOfSigmaProtonMuon_EnergyDepositionLength_%d",is),"",hEnergyDepositionLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->Sumw2();

   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is] = new TH1D(Form("NumberOfSigmaPionMuon_EnergyDepositionSplineLength_%d",is),"",hEnergyDepositionSplineLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->Sumw2();
   
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is] = new TH1D(Form("NumberOfSigmaProtonMuon_EnergyDepositionSplineLength_%d",is),"",hEnergyDepositionSplineLength_Muon_1D[is]->GetNbinsX(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmin(),hEnergyDepositionSplineLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->Sumw2();

   for(int ibinx=1;ibinx<=RatioPionMuon_EnergyDepositionLength[is]->GetNbinsX();ibinx++){

     double ValueMuon = hEnergyDepositionLength_Muon_1D[is]->GetBinContent(ibinx);
     double ErrorMuon = hEnergyDepositionLength_Muon_1D[is]->GetBinError(ibinx);
     double ValuePion = hEnergyDepositionLength_Pion_1D[is]->GetBinContent(ibinx);
     double ErrorPion = hEnergyDepositionLength_Pion_1D[is]->GetBinError(ibinx);
     double ValueProton = hEnergyDepositionLength_Proton_1D[is]->GetBinContent(ibinx);
     double ErrorProton = hEnergyDepositionLength_Proton_1D[is]->GetBinError(ibinx);

     double RatioPionMuon = ValuePion/ValueMuon ;
     double ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorPion/ValuePion),2.));
     double RatioProtonMuon = ValueProton/ValueMuon ;
     double ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorProton/ValueProton),2.));

     RatioPionMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_EnergyDepositionLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_EnergyDepositionLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     double NumberOfSigmaPionMuon = RatioPionMuon/ErrorPionMuon ;
     double NumberOfSigmaProtonMuon = RatioProtonMuon/ErrorProtonMuon ;

     NumberOfSigmaPionMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);
     
     //Spline
     ValueMuon = hEnergyDepositionSplineLength_Muon_1D[is]->GetBinContent(ibinx);
     ErrorMuon = hEnergyDepositionSplineLength_Muon_1D[is]->GetBinError(ibinx);
     ValuePion = hEnergyDepositionSplineLength_Pion_1D[is]->GetBinContent(ibinx);
     ErrorPion = hEnergyDepositionSplineLength_Pion_1D[is]->GetBinError(ibinx);
     ValueProton = hEnergyDepositionSplineLength_Proton_1D[is]->GetBinContent(ibinx);
     ErrorProton = hEnergyDepositionSplineLength_Proton_1D[is]->GetBinError(ibinx);
     
     RatioPionMuon = ValuePion/ValueMuon ;
     ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorPion/ValuePion),2.));
     RatioProtonMuon = ValueProton/ValueMuon ;
     ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorProton/ValueProton),2.));

     RatioPionMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_EnergyDepositionSplineLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_EnergyDepositionSplineLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     NumberOfSigmaPionMuon = RatioPionMuon/ErrorPionMuon ;
     NumberOfSigmaProtonMuon = RatioProtonMuon/ErrorProtonMuon ;

     NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);

     //cout<<"ibin="<<ibinx<<", Ratio="<<RatioPionMuon<<" +- "<<ErrorPionMuon<<endl;
     cout<<"ibin="<<ibinx<<": Ratio="<<RatioPionMuon_EnergyDepositionLength[is]->GetBinContent(ibinx)<<" +- "<<RatioPionMuon_EnergyDepositionLength[is]->GetBinError(ibinx)<<endl;
     cout<<"Error pion="<<ErrorPion<<", error muon="<<ErrorMuon<<endl;
     //cout<<" Spline, Ratio="<<RatioPionMuon_EnergyDepositionSplineLength[is]->GetBinContent(ibinx)<<" +- "<<RatioPionMuon_EnergyDepositionSplineLength[is]->GetBinError(ibinx)<<endl;

     
   }

   cRatio_EnergyDepositionLength[is] = new TCanvas(Form("cRatio_EnergyDepositionLength_%d",is),Form("Ratio of energy deposition with length of sample %d",is));
   cRatio_EnergyDepositionLength[is]->Divide(1,3);
   cRatio_EnergyDepositionLength[is]->cd(1);
   hEnergyDepositionLength_Muon_1D[is]->SetLineColor(kBlue);
   hEnergyDepositionLength_Pion_1D[is]->SetLineColor(kGreen+2);
   hEnergyDepositionLength_Proton_1D[is]->SetLineColor(kRed);
   hEnergyDepositionSplineLength_Muon_1D[is]->SetLineColor(kBlue);
   hEnergyDepositionSplineLength_Pion_1D[is]->SetLineColor(kGreen+2);
   hEnergyDepositionSplineLength_Proton_1D[is]->SetLineColor(kRed);
   hEnergyDepositionSplineLength_Muon_1D[is]->SetLineStyle(2);
   hEnergyDepositionSplineLength_Pion_1D[is]->SetLineStyle(2);
   hEnergyDepositionSplineLength_Proton_1D[is]->SetLineStyle(2);
   hEnergyDepositionLength_Muon_1D[is]->Draw("HIST");
   hEnergyDepositionLength_Muon_1D[is]->GetYaxis()->SetRangeUser(20,80);
   hEnergyDepositionLength_Pion_1D[is]->Draw("HISTsame");
   hEnergyDepositionLength_Proton_1D[is]->Draw("HISTsame");
   hEnergyDepositionSplineLength_Muon_1D[is]->Draw("HISTsame");
   hEnergyDepositionSplineLength_Pion_1D[is]->Draw("HISTsame");
   hEnergyDepositionSplineLength_Proton_1D[is]->Draw("HISTsame");

   cRatio_EnergyDepositionLength[is]->cd(2);
   RatioPionMuon_EnergyDepositionLength[is]->SetLineColor(kGreen+2);
   RatioPionMuon_EnergyDepositionSplineLength[is]->SetLineColor(kGreen+2);
   RatioProtonMuon_EnergyDepositionLength[is]->SetLineColor(kRed);
   RatioProtonMuon_EnergyDepositionSplineLength[is]->SetLineColor(kRed);

   RatioPionMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   RatioProtonMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   
   RatioPionMuon_EnergyDepositionLength[is]->Draw("HIST");   
   RatioPionMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");   
   RatioProtonMuon_EnergyDepositionLength[is]->Draw("HISTsame");   
   RatioProtonMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");   
   RatioPionMuon_EnergyDepositionLength[is]->GetYaxis()->SetRangeUser(0.8,1.6);
     
   cRatio_EnergyDepositionLength[is]->cd(3);
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->SetLineColor(kGreen+2);
   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->SetLineColor(kGreen+2);
   NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->SetLineColor(kRed);
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->SetLineColor(kRed);

   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->Draw("HIST");   
   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");   
   NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->Draw("HISTsame");   
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");   
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->GetYaxis()->SetRangeUser(1.,3.);
   cRatio_EnergyDepositionLength[is]->SaveAs(Form("../plots/MVAOptimisation/cRatio_EnergyDepositionLength_sample%d.eps",is));
 }				 
#endif

 
#define TRANSVERSEPLOTS
#ifdef TRANSVERSEPLOTS

 int NSamples=6;
 TH1D * hTransverseWidthLength_Muon_1D[NSamples];
 TH1D * hTransverseWidthLength_Pion_1D[NSamples];
 TH1D * hTransverseWidthLength_Proton_1D[NSamples];
 TH1D * hTransverseWidthNonIsolatedLength_Muon_1D[NSamples];
 TH1D * hTransverseWidthNonIsolatedLength_Pion_1D[NSamples];
 TH1D * hTransverseWidthNonIsolatedLength_Proton_1D[NSamples];

 //Ratio plots
 TCanvas * cRatio_TransverseWidthLength[NSamples];
 //Ratio
 TH1D * RatioPionMuon_TransverseWidthLength[NSamples];
 TH1D * RatioProtonMuon_TransverseWidthLength[NSamples];
 TH1D * RatioPionMuon_TransverseWidthNonIsolatedLength[NSamples];
 TH1D * RatioProtonMuon_TransverseWidthNonIsolatedLength[NSamples];
 //(Ratio -1)/Error
 TH1D * NumberOfSigmaPionMuon_TransverseWidthLength[NSamples];
 TH1D * NumberOfSigmaProtonMuon_TransverseWidthLength[NSamples];
 TH1D * NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[NSamples];
 TH1D * NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[NSamples];
 

 for(int is=0;is<NSamples;is++){
   hTransverseWidthLength_Muon_1D[is] = (TH1D*) fMC->Get(Form("TransverseWidthLength_Muon_%d_y",is));
   hTransverseWidthLength_Pion_1D[is] = (TH1D*) fMC->Get(Form("TransverseWidthLength_Pion_%d_y",is));
   hTransverseWidthLength_Proton_1D[is] = (TH1D*) fMC->Get(Form("TransverseWidthLength_Proton_%d_y",is));

   hTransverseWidthNonIsolatedLength_Muon_1D[is] = (TH1D*) fMC->Get(Form("TransverseWidthNonIsolatedLength_Muon_%d_y",is));
   hTransverseWidthNonIsolatedLength_Pion_1D[is] = (TH1D*) fMC->Get(Form("TransverseWidthNonIsolatedLength_Pion_%d_y",is));
   hTransverseWidthNonIsolatedLength_Proton_1D[is] = (TH1D*) fMC->Get(Form("TransverseWidthNonIsolatedLength_Proton_%d_y",is));
   
   RatioPionMuon_TransverseWidthLength[is] = new TH1D(Form("RatioPionMuon_TransverseWidthLength_%d",is),"",hTransverseWidthLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioPionMuon_TransverseWidthLength[is]->Sumw2();
   
   RatioProtonMuon_TransverseWidthLength[is] = new TH1D(Form("RatioProtonMuon_TransverseWidthLength_%d",is),"",hTransverseWidthLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioProtonMuon_TransverseWidthLength[is]->Sumw2();

   RatioPionMuon_TransverseWidthNonIsolatedLength[is] = new TH1D(Form("RatioPionMuon_TransverseWidthNonIsolatedLength_%d",is),"",hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioPionMuon_TransverseWidthNonIsolatedLength[is]->Sumw2();
   
   RatioProtonMuon_TransverseWidthNonIsolatedLength[is] = new TH1D(Form("RatioProtonMuon_TransverseWidthNonIsolatedLength_%d",is),"",hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmax());
   RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->Sumw2();

   //
   NumberOfSigmaPionMuon_TransverseWidthLength[is] = new TH1D(Form("NumberOfSigmaPionMuon_TransverseWidthLength_%d",is),"",hTransverseWidthLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaPionMuon_TransverseWidthLength[is]->Sumw2();
   
   NumberOfSigmaProtonMuon_TransverseWidthLength[is] = new TH1D(Form("NumberOfSigmaProtonMuon_TransverseWidthLength_%d",is),"",hTransverseWidthLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaProtonMuon_TransverseWidthLength[is]->Sumw2();

   NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[is] = new TH1D(Form("NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength_%d",is),"",hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[is]->Sumw2();
   
   NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[is] = new TH1D(Form("NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength_%d",is),"",hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetNbinsX(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmin(),hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetXaxis()->GetXmax());
   NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[is]->Sumw2();

   for(int ibinx=1;ibinx<=RatioPionMuon_TransverseWidthLength[is]->GetNbinsX();ibinx++){

     double ValueMuon = hTransverseWidthLength_Muon_1D[is]->GetBinContent(ibinx);
     double ErrorMuon = hTransverseWidthLength_Muon_1D[is]->GetBinError(ibinx);
     double ValuePion = hTransverseWidthLength_Pion_1D[is]->GetBinContent(ibinx);
     double ErrorPion = hTransverseWidthLength_Pion_1D[is]->GetBinError(ibinx);
     double ValueProton = hTransverseWidthLength_Proton_1D[is]->GetBinContent(ibinx);
     double ErrorProton = hTransverseWidthLength_Proton_1D[is]->GetBinError(ibinx);

     double RatioPionMuon = ValuePion/ValueMuon ;
     double ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorPion/ValuePion),2.));
     double RatioProtonMuon = ValueProton/ValueMuon ;
     double ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorProton/ValueProton),2.));

     RatioPionMuon_TransverseWidthLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_TransverseWidthLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_TransverseWidthLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_TransverseWidthLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     double NumberOfSigmaPionMuon = RatioPionMuon/ErrorPionMuon ;
     double NumberOfSigmaProtonMuon = RatioProtonMuon/ErrorProtonMuon ;

     NumberOfSigmaPionMuon_TransverseWidthLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_TransverseWidthLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);
     
     //NonIsolated
     ValueMuon = hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetBinContent(ibinx);
     ErrorMuon = hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetBinError(ibinx);
     ValuePion = hTransverseWidthNonIsolatedLength_Pion_1D[is]->GetBinContent(ibinx);
     ErrorPion = hTransverseWidthNonIsolatedLength_Pion_1D[is]->GetBinError(ibinx);
     ValueProton = hTransverseWidthNonIsolatedLength_Proton_1D[is]->GetBinContent(ibinx);
     ErrorProton = hTransverseWidthNonIsolatedLength_Proton_1D[is]->GetBinError(ibinx);
     
     RatioPionMuon = ValuePion/ValueMuon ;
     ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorPion/ValuePion),2.));
     RatioProtonMuon = ValueProton/ValueMuon ;
     ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ValueMuon),2.)+pow((ErrorProton/ValueProton),2.));

     RatioPionMuon_TransverseWidthNonIsolatedLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_TransverseWidthNonIsolatedLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     NumberOfSigmaPionMuon = RatioPionMuon/ErrorPionMuon ;
     NumberOfSigmaProtonMuon = RatioProtonMuon/ErrorProtonMuon ;

     NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);

     //cout<<"ibin="<<ibinx<<", Ratio="<<RatioPionMuon<<" +- "<<ErrorPionMuon<<endl;
     cout<<"ibin="<<ibinx<<": Ratio="<<RatioPionMuon_TransverseWidthLength[is]->GetBinContent(ibinx)<<" +- "<<RatioPionMuon_TransverseWidthLength[is]->GetBinError(ibinx)<<endl;
     cout<<"Error pion="<<ErrorPion<<", error muon="<<ErrorMuon<<endl;
     //cout<<" NonIsolated, Ratio="<<RatioPionMuon_TransverseWidthNonIsolatedLength[is]->GetBinContent(ibinx)<<" +- "<<RatioPionMuon_TransverseWidthNonIsolatedLength[is]->GetBinError(ibinx)<<endl;

     
   }

   cRatio_TransverseWidthLength[is] = new TCanvas(Form("cRatio_TransverseWidthLength_%d",is),Form("Ratio of energy deposition with length of sample %d",is));
   cRatio_TransverseWidthLength[is]->Divide(1,3);
   cRatio_TransverseWidthLength[is]->cd(1);
   hTransverseWidthLength_Muon_1D[is]->SetLineColor(kBlue);
   hTransverseWidthLength_Pion_1D[is]->SetLineColor(kGreen+2);
   hTransverseWidthLength_Proton_1D[is]->SetLineColor(kRed);
   hTransverseWidthNonIsolatedLength_Muon_1D[is]->SetLineColor(kBlue);
   hTransverseWidthNonIsolatedLength_Pion_1D[is]->SetLineColor(kGreen+2);
   hTransverseWidthNonIsolatedLength_Proton_1D[is]->SetLineColor(kRed);
   hTransverseWidthNonIsolatedLength_Muon_1D[is]->SetLineStyle(2);
   hTransverseWidthNonIsolatedLength_Pion_1D[is]->SetLineStyle(2);
   hTransverseWidthNonIsolatedLength_Proton_1D[is]->SetLineStyle(2);
   hTransverseWidthLength_Muon_1D[is]->Draw("HIST");
   hTransverseWidthLength_Muon_1D[is]->GetYaxis()->SetRangeUser(0,20.);
   hTransverseWidthLength_Pion_1D[is]->Draw("HISTsame");
   hTransverseWidthLength_Proton_1D[is]->Draw("HISTsame");
   hTransverseWidthNonIsolatedLength_Muon_1D[is]->Draw("HISTsame");
   hTransverseWidthNonIsolatedLength_Pion_1D[is]->Draw("HISTsame");
   hTransverseWidthNonIsolatedLength_Proton_1D[is]->Draw("HISTsame");

   cRatio_TransverseWidthLength[is]->cd(2);
   RatioPionMuon_TransverseWidthLength[is]->SetLineColor(kGreen+2);
   RatioPionMuon_TransverseWidthNonIsolatedLength[is]->SetLineColor(kGreen+2);
   RatioProtonMuon_TransverseWidthLength[is]->SetLineColor(kRed);
   RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->SetLineColor(kRed);

   RatioPionMuon_TransverseWidthNonIsolatedLength[is]->SetLineStyle(2);
   RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->SetLineStyle(2);
   
   RatioPionMuon_TransverseWidthLength[is]->Draw("HIST");   
   RatioPionMuon_TransverseWidthNonIsolatedLength[is]->Draw("HISTsame");   
   RatioProtonMuon_TransverseWidthLength[is]->Draw("HISTsame");   
   RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->Draw("HISTsame");   
   RatioPionMuon_TransverseWidthLength[is]->GetYaxis()->SetRangeUser(0.,3.);
     
   cRatio_TransverseWidthLength[is]->cd(3);
   NumberOfSigmaPionMuon_TransverseWidthLength[is]->SetLineColor(kGreen+2);
   NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[is]->SetLineColor(kGreen+2);
   NumberOfSigmaProtonMuon_TransverseWidthLength[is]->SetLineColor(kRed);
   NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[is]->SetLineColor(kRed);

   NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[is]->SetLineStyle(2);
   NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[is]->SetLineStyle(2);
   
   NumberOfSigmaPionMuon_TransverseWidthLength[is]->Draw("HIST");   
   NumberOfSigmaPionMuon_TransverseWidthNonIsolatedLength[is]->Draw("HISTsame");   
   NumberOfSigmaProtonMuon_TransverseWidthLength[is]->Draw("HISTsame");   
   NumberOfSigmaProtonMuon_TransverseWidthNonIsolatedLength[is]->Draw("HISTsame");   
   NumberOfSigmaPionMuon_TransverseWidthLength[is]->GetYaxis()->SetRangeUser(0.,3.);
   cRatio_TransverseWidthLength[is]->SaveAs(Form("../plots/MVAOptimisation/cRatio_TransverseWidthLength_sample%d.eps",is));
 }				 
#endif



 /*    

  cout<<"1 Track sample:"<<endl;
  cout<<"Maximal purity = "<<dCC0pi_Purity_1track<<", reachable with an efficiency = "<<dCC0pi_Efficiency_1track<<endl;
  cout<<"2 Tracks sample:"<<endl;
  cout<<"Maximal purity = "<<dCC0pi_Purity_2tracks<<", reachable with an efficiency = "<<dCC0pi_Efficiency_2tracks<<endl;
  */
#ifdef WRITE  
    TFile * wfile = new TFile("plotsclmuon.root","RECREATE");
    cHighCL->Write();
#ifdef CLMuon
    cEff_1track->Write();
    cEff_2tracks->Write();
    cPur_2tracks->Write();
#endif
#ifdef MVA
    cmvaEfficiency_1track->Write();
    cmvaEfficiency_2tracks->Write();
#endif
    wfile->Close();
#endif
  //CC0pi_Purity_2tracks->Draw("same");
  /*
  cout<<"Please choose a cut value for mu-like particles"<<endl;
  cin>>mulike;
cout<<"Please choose a cut value for 
*/
}

//Macro helpful for efficiency / purity optimization
{
  gROOT->SetBatch(kTRUE);//Silent mode
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetOptStat(kFALSE);

  double ErrorBkg=0.1;
  int NBinsTrueMom=4;
  int NBinsTrueAngle=4;
  
  //TFile * fMC = new TFile("../files/MCSelected_PM_PG_Systematics0_0PM.root");
  //TFile * fMC = new TFile("../files/MCSelected_PM_Systematics0_0.root");
  //TFile * fMC = new TFile("../files/MCSelected_WM_Systematics0_0_StoppedOnly.root");
  TFile * fMC = new TFile("../test.root.root");
  //TFile * fMC = new TFile("../files/MCSelected_PM_Systematics0_0PM_StoppedOnly.root");
  
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
  //TDirectory * PIDCutTuning = (TDirectory*) fMC->Get("PIDCutTuning");
  //PIDCutTuning->cd();
#define CLParticles
#ifdef CLParticles
  TH1D * MuCL_TrueMuon = (TH1D*) PIDCutTuning->Get("hMuCL_TrueMuon");
  TH1D * MuCL_TruePion = (TH1D*) PIDCutTuning->Get("hMuCL_TruePion");
  TH1D * MuCL_TrueProton = (TH1D*) PIDCutTuning->Get("hMuCL_TrueProton");
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
  cCL->SaveAs("../plots/MVAOptimisation/MuCL_particles.eps");
  
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

  /////////////////////////////////////////////////////////////////////
  int NSimplifiedPDG=4;
  TH2D * hMVAMuondiscriminantVSDistance_TrueParticle[NSimplifiedPDG];
  for(int i=0;i<NSimplifiedPDG;i++){
    hMVAMuondiscriminantVSDistance_TrueParticle[i] = (TH2D*) PIDCutTuning->Get(Form("hMVAMuondiscriminantVSDistance_TrueParticle%d",i));
  }
  TCanvas * cMVAMuondiscriminantVSDistance_TrueParticle = new TCanvas("cMVAMuondiscriminantVSDistance_TrueParticle","True particle vs distance and MVA");
  hMVAMuondiscriminantVSDistance_TrueParticle[0]->Draw("colz");
  cMVAMuondiscriminantVSDistance_TrueParticle->SaveAs("../plots/MVAOptimisation/MVAMuondiscriminantVSDistance_TrueParticle.eps");



  double ThresholdEfficiencyMuon=0.5;//Efficiency that we wish to keep constant
  TH2D * EfficiencyMuonVSDistance = (TH2D*) hMVAMuondiscriminantVSDistance_TrueParticle[0]->Clone("EfficiencyMuonVSDistance");
  EfficiencyMuonVSDistance->Reset();

  double distance[hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsX()];
  double mva[hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsX()];
  double errorx[hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsX()];//assume the error on distance is only the bin size in x
  double errory[hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsX()];//assume the error on mva is only the bin size in y

  for(int ibinx=1;ibinx<=hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsX();ibinx++){
    bool Found=false;//true when find the intersection between efficiency 2D plot and ThresholdEfficiencyMuon constant plane
    for(int ibiny=1;ibiny<=hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsY();ibiny++){
      double NMuon = hMVAMuondiscriminantVSDistance_TrueParticle[0]->Integral(ibinx,ibinx,ibiny,hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsY());
      double TotalMuon=hMVAMuondiscriminantVSDistance_TrueParticle[0]->Integral(ibinx,ibinx,1,hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsY());
      double Efficiency = NMuon / (TotalMuon !=0 ? TotalMuon : 1. );
      EfficiencyMuonVSDistance->SetBinContent(ibinx,ibiny,Efficiency);

      //Given the binning, we start from the largest to the lowest efficiency in a bin x
      if(!Found && Efficiency<ThresholdEfficiencyMuon){
	distance[ibinx-1]=hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetXaxis()->GetBinCenter(ibinx);
	mva[ibinx-1]=hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetYaxis()->GetBinCenter(ibiny);
	errorx[ibinx-1]=hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetXaxis()->GetBinWidth(ibinx)/2;
	errory[ibinx-1]=hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetYaxis()->GetBinWidth(ibiny);
	Found=true;
	//cout<<"For distance="<<hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetXaxis()->GetBinCenter(ibinx)<<", MVA value found is="<<hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetYaxis()->GetBinCenter(ibiny)<<endl;
	//cout<<"For distance="<<MVAVSDistance_ConstantEfficiency->GetXaxis()->GetBinCenter(ibinx)<<", MVA value found is="<<MVAVSDistance_ConstantEfficiency->GetBinContent(ibinx)<<endl;
	cout<<"Distance="<<distance[ibinx-1]<<", "<<mva[ibinx-1]<<endl;
      }
    }
  }
  TCanvas * cEfficiencyMuonVSDistance = new TCanvas("cEfficiencyMuonVSDistance","True particle vs distance and MVA");
  cEfficiencyMuonVSDistance->Divide(1,2);
  cEfficiencyMuonVSDistance->cd(1);
  EfficiencyMuonVSDistance->Draw("colz");
  cEfficiencyMuonVSDistance->cd(2);

  TGraphErrors * gMVAVSDistance_ConstantEfficiency = new TGraphErrors(hMVAMuondiscriminantVSDistance_TrueParticle[0]->GetNbinsX(),distance,mva,errorx,errory);
  TF1 * fpol1 = new TF1("fpol1","[0]+[1]*x",17.5,57.5);
  //fpol1->SetParameter(0,1.);
  //fpol1->SetParameter(1,1.);
  gMVAVSDistance_ConstantEfficiency->Fit("fpol1","R");
  gMVAVSDistance_ConstantEfficiency->Draw("AP");
  fpol1->Draw("same");

  cEfficiencyMuonVSDistance->SaveAs("../plots/MVAOptimisation/EfficiencyMuonVSDistance.eps");
  
#endif
  
  /////////////////////////////////////////////DRAW VARIABLES FOR INTERACTIONS//////////////////////////////////////////////
  int NFSIs=13;
  char Name[256];
#define CLMuon
#ifdef CLMuon
  

  int RebinValueCL=25;

  TH1D * hMuCL[NFSIs];
  TH2D * hMuCL_2tracks[NFSIs];
  TH1D * hMuCL_CC0pi = (TH1D*) PIDCutTuning->Get("hMuCL0");
  TH2D * hMuCL_2tracks_CC0pi = (TH2D*) PIDCutTuning->Get("hMuCL_2tracks0");
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
    hMuCL[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hMuCL[fsi]->SetFillStyle(0);

    sprintf(Name,"hMuCL_2tracks%d",fsi);
    hMuCL_2tracks[fsi] = (TH2D*) PIDCutTuning->Get(Name);
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
  TH1D * CC0pi_SignalBkg_1track = (TH1D*) hMuCL_CC0pi->Clone("CC0pi_SignalBkg_1track");
  
  TH2D * CC0pi_Purity_2tracks = (TH2D*) hMuCL_2tracks_CC0pi->Clone("CC0pi_Purity_2tracks");
  TH2D * CC0pi_Efficiency_2tracks = (TH2D*) hMuCL_2tracks_CC0pi->Clone("CC0pi_Efficiency_2tracks");
  TH2D * CC0pi_SignalBkg_2tracks = (TH2D*) hMuCL_2tracks_CC0pi->Clone("CC0pi_SignalBkg_2tracks");

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
    
    //double SqrtSB_1track = ErrorBkg*(NAll_1track - NCC0pi_1track) / pow(NCC0pi_1track,3/2);
    double SqrtSB_1track = sqrt( pow(sqrt(NAll_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_1track-NCC0pi_1track)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
    if(NCC0pi_1track!=0) SqrtSB_1track/=(NCC0pi_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
      CC0pi_SignalBkg_1track->SetBinContent(ibinx,ibiny,SqrtSB_1track);
    
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

      //double SqrtSB_2tracks = ErrorBkg*(NAll_2tracks - NCC0pi_2tracks) / pow(NCC0pi_2tracks,3/2);
    double SqrtSB_2tracks = sqrt( pow(sqrt(NAll_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_2tracks-NCC0pi_2tracks)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
    if(NCC0pi_2tracks!=0) SqrtSB_2tracks/=(NCC0pi_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
      CC0pi_SignalBkg_2tracks->SetBinContent(ibinx,ibiny,SqrtSB_2tracks);
      }
  }

  TCanvas * cEff_1track = new TCanvas("cEff_1track","Tune mu-like cut: Efficiency/Purity of CC0pi with MuCL cut value");//Use the 1 track sample
  CC0pi_Efficiency_1track->SetLineColor(kRed);
  CC0pi_Efficiency_1track->Draw();
  CC0pi_Purity_1track->SetLineColor(kBlue);
  CC0pi_Purity_1track->Draw("same");
  cEff_1track->SaveAs("../plots/MVAOptimisation/CC0pi_1track_EfficiencyPurity.eps");

  TCanvas * cSignalBkg_1track = new TCanvas("cSignalBkg_1track","Tune mu-like cut: SignalBkg of CC0pi with MuCL cut value");//Use the 1 track sample
  CC0pi_SignalBkg_1track->Draw("same");
  cSignalBkg_1track->SaveAs("../plots/MVAOptimisation/CC0pi_1track_SignalBkg.eps");
  
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

  TCanvas * cSignalBkg_2tracks = new TCanvas("cSignalBkg_2tracks","Tune p-like cut: SignalBkg of CC0pi with MuCL cut value");
  CC0pi_SignalBkg_2tracks->Draw("colztext");
  CC0pi_SignalBkg_2tracks->GetXaxis()->SetRangeUser(0,1);
  CC0pi_SignalBkg_2tracks->GetYaxis()->SetRangeUser(0,1);
  cSignalBkg_2tracks->SaveAs("../plots/MVAOptimisation/CC0pi_2tracks_SignalBkg.eps");
  
#endif
  
#define MVA  
#ifdef MVA
  double MVAZoomInf=-0.1;
  double MVAZoomSup=0.15;
  int RebinValueMVA=1;
  //For the LMVA, the unified discriminator
  double LMVAZoomInf=-1;
  double LMVAZoomSup=1;
  int RebinValueLMVA=5;

  TH1D * hMVAMuondiscriminant_1track[NFSIs];
  TH2D * hMVAMuonVSProtondiscriminant_2tracks[NFSIs];

  TH1D * hMVAMuondiscriminant_1track_CC0pi = (TH1D*) PIDCutTuning->Get("hMVAMuondiscriminant_1track0");
  TH2D * hMVAMuonVSProtondiscriminant_2tracks_CC0pi = (TH2D*) PIDCutTuning->Get("hMVAMuonVSProtondiscriminant_2tracks0");

  TH1D * hMVAMuondiscriminant_1track_All = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("hMVAdiscriminant_1track_All");
  TH2D * hMVAMuonVSProtondiscriminant_2tracks_All = (TH2D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("hMVAMuonVSProtondiscriminant_2tracks_All");
  
  hMVAMuondiscriminant_1track_CC0pi->SetFillStyle(0);
  hMVAMuonVSProtondiscriminant_2tracks_CC0pi->SetFillStyle(0);
  hMVAMuondiscriminant_1track_All->SetFillStyle(0);
  hMVAMuonVSProtondiscriminant_2tracks_All->SetFillStyle(0);
  
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track[NFSIs];
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi = (TH2D*) PIDCutTuning->Get("hMVAMuondiscriminantVSMuonMomentum_1track0");
  TH2D * hMVAMuondiscriminantVSMuonMomentum_1track_All = (TH2D*) hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->Clone("hMVAMuondiscriminantVSMuonMomentum_1track_All");
  hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->SetFillStyle(0);
  hMVAMuondiscriminantVSMuonMomentum_1track_All->SetFillStyle(0);

  TH2D * hMVAMuondiscriminantVSDistance_1track[NFSIs];
  TH2D * hMVAMuondiscriminantVSDistance_1track_CC0pi = (TH2D*) PIDCutTuning->Get("hMVAMuondiscriminantVSDistance_1track0");
  TH2D * hMVAMuondiscriminantVSDistance_1track_All = (TH2D*) hMVAMuondiscriminantVSDistance_1track_CC0pi->Clone("hMVAMuondiscriminantVSDistance_1track_All");
  hMVAMuondiscriminantVSDistance_1track_CC0pi->SetFillStyle(0);
  hMVAMuondiscriminantVSDistance_1track_All->SetFillStyle(0);

  TH2D * hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[NFSIs];
  TH2D * hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_CC0pi = (TH2D*) PIDCutTuning->Get("hMVAProtondiscriminantVSDistance_2tracks_LowestMVA0");
  TH2D * hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_All = (TH2D*) hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_CC0pi->Clone("hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_All");
  hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_CC0pi->SetFillStyle(0);
  hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_All->SetFillStyle(0);

  //New LMVA discriminator
  cout<<"New LMVA discriminator"<<endl;
  TH1D * hLMVA_1track[NFSIs];
  TH2D * hLMVA_2tracks[NFSIs];
  TH1D * hLMVA_LargestLMVA[NFSIs];

  TH1D * hLMVA_1track_CC0pi = (TH1D*) PIDCutTuning->Get("hLMVA_1track0");
  hLMVA_1track_CC0pi->SetFillStyle(0);
  TH1D * hLMVA_1track_All = (TH1D*) hLMVA_1track_CC0pi->Clone("hLMVA_1track_All");
  TH2D * hLMVA_2tracks_CC0pi = (TH2D*) PIDCutTuning->Get("hLMVA_2tracks0");
  TH2D * hLMVA_2tracks_All = (TH2D*) hLMVA_2tracks_CC0pi->Clone("hLMVA_2tracks_All");
  TH1D * hLMVA_LargestLMVA_CC0pi = (TH1D*) PIDCutTuning->Get("hLMVA_LargestLMVA0");
  hLMVA_LargestLMVA_CC0pi->SetFillStyle(0);
  TH1D * hLMVA_LargestLMVA_All = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("hLMVA_LargestLMVA_All");
  cout<<"done"<<endl;


  
  //For particle per particle testing
  TH2D * hMVAMuondiscriminantVSPDG_1track[NFSIs];
  TH2D * hMVAMuondiscriminantVSPDG_1track_All;
  TH2D * hMVAMuondiscriminantVSPDG_1track_All;

  TH2D * hMVAMuonVSProtondiscriminant_2tracks_CC1pi = (TH2D*) PIDCutTuning->Get("hMVAMuonVSProtondiscriminant_2tracks3");
  hMVAMuonVSProtondiscriminant_2tracks_CC1pi->SetFillStyle(0);

  for(int fsi=1;fsi<NFSIs;fsi++){
    sprintf(Name,"hMVAMuondiscriminant_1track%d",fsi);
    hMVAMuondiscriminant_1track[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hMVAMuondiscriminant_1track[fsi]->SetFillStyle(0);

    //
    sprintf(Name,"hMVAMuondiscriminantVSPDG_1track%d",fsi);
    hMVAMuondiscriminantVSPDG_1track[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hMVAMuondiscriminantVSPDG_1track[fsi]->SetFillStyle(0);
    
    if(fsi==1) hMVAMuondiscriminantVSPDG_1track_All = (TH2D*) hMVAMuondiscriminantVSPDG_1track[fsi]->Clone("hMVAMuondiscriminantVSPDG_1track_All");
    else hMVAMuondiscriminantVSPDG_1track_All->Add(hMVAMuondiscriminantVSPDG_1track[fsi]);
    
    
    sprintf(Name,"hMVAMuondiscriminantVSMuonMomentum_1track%d",fsi);
    hMVAMuondiscriminantVSMuonMomentum_1track[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hMVAMuondiscriminantVSMuonMomentum_1track[fsi]->SetFillStyle(0);
    
    sprintf(Name,"hMVAMuondiscriminantVSDistance_1track%d",fsi);
    hMVAMuondiscriminantVSDistance_1track[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hMVAMuondiscriminantVSDistance_1track[fsi]->SetFillStyle(0);
    
    sprintf(Name,"hMVAMuonVSProtondiscriminant_2tracks%d",fsi);
    hMVAMuonVSProtondiscriminant_2tracks[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hMVAMuonVSProtondiscriminant_2tracks[fsi]->SetFillStyle(0);

    sprintf(Name,"hMVAProtondiscriminantVSDistance_2tracks_LowestMVA%d",fsi);
    hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[fsi]->SetFillStyle(0);

    sprintf(Name,"hLMVA_1track%d",fsi);
    hLMVA_1track[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hLMVA_1track[fsi]->SetFillStyle(0);
    
    sprintf(Name,"hLMVA_2tracks%d",fsi);
    hLMVA_2tracks[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hLMVA_2tracks[fsi]->SetFillStyle(0);

    sprintf(Name,"hLMVA_LargestLMVA%d",fsi);
    hLMVA_LargestLMVA[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hLMVA_LargestLMVA[fsi]->SetFillStyle(0);
    

    if(fsi>-1){
      if(fsi<3){
	hMVAMuondiscriminant_1track_CC0pi->Add(hMVAMuondiscriminant_1track[fsi]);
	hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Add(hMVAMuonVSProtondiscriminant_2tracks[fsi]); 
	hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->Add(hMVAMuondiscriminantVSMuonMomentum_1track[fsi]);
	hMVAMuondiscriminantVSDistance_1track_CC0pi->Add(hMVAMuondiscriminantVSDistance_1track[fsi]);
	hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_CC0pi->Add(hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[fsi]);

	hLMVA_1track_CC0pi->Add(hLMVA_1track[fsi]);
	hLMVA_2tracks_CC0pi->Add(hLMVA_2tracks[fsi]);
	hLMVA_LargestLMVA_CC0pi->Add(hLMVA_LargestLMVA[fsi]);
      }

      hMVAMuondiscriminant_1track_All->Add(hMVAMuondiscriminant_1track[fsi]);
      hMVAMuonVSProtondiscriminant_2tracks_All->Add(hMVAMuonVSProtondiscriminant_2tracks[fsi]);
      hMVAMuondiscriminantVSMuonMomentum_1track_All->Add(hMVAMuondiscriminantVSMuonMomentum_1track[fsi]);
      hMVAMuondiscriminantVSDistance_1track_All->Add(hMVAMuondiscriminantVSDistance_1track[fsi]);
      hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_All->Add(hMVAProtondiscriminantVSDistance_2tracks_LowestMVA[fsi]);

      hLMVA_1track_All->Add(hLMVA_1track[fsi]);
      hLMVA_2tracks_All->Add(hLMVA_2tracks[fsi]);
      hLMVA_LargestLMVA_All->Add(hLMVA_LargestLMVA[fsi]);
    } 
  }

  hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Rebin2D(RebinValueMVA,RebinValueMVA);  
  hMVAMuonVSProtondiscriminant_2tracks_CC1pi->Rebin2D(RebinValueMVA,RebinValueMVA);  
  hMVAMuonVSProtondiscriminant_2tracks_All->Rebin2D(RebinValueMVA,RebinValueMVA);  

  hLMVA_1track_CC0pi->Rebin(RebinValueLMVA);
  hLMVA_1track_All->Rebin(RebinValueLMVA);
  hLMVA_LargestLMVA_CC0pi->Rebin(RebinValueLMVA);
  hLMVA_LargestLMVA_All->Rebin(RebinValueLMVA);
  hLMVA_2tracks_CC0pi->Rebin2D(RebinValueLMVA,RebinValueLMVA);
  hLMVA_2tracks_All->Rebin2D(RebinValueLMVA,RebinValueLMVA);

  TH1D * mvaCC0pi_Purity_1track = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("mvaCC0pi_Purity_1track");
  TH1D * mvaCC0pi_Efficiency_1track = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("mvaCC0pi_Efficiency_1track");
  TH1D * mvaCC0pi_SignalBkg_1track = (TH1D*) hMVAMuondiscriminant_1track_CC0pi->Clone("mvaCC0pi_SignalBkg_1track");
  TH1D * mvaCC0pi_Purity_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("mvaCC0pi_Purity_2tracks");
  TH1D * mvaCC0pi_Efficiency_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("mvaCC0pi_Efficiency_2tracks");
  TH1D * mvaCC0pi_SignalBkg_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Clone("mvaCC0pi_SignalBkg_2tracks");

  double mvaNBinsX=mvaCC0pi_Efficiency_1track->GetNbinsX();
  double mvaNBinsY=mvaCC0pi_Efficiency_1track->GetNbinsX();
  double mvaNCC0pi_Total_1track=hMVAMuondiscriminant_1track_CC0pi->Integral(1,mvaNBinsX);
  double mvaNCC0pi_Total_2tracks=hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Integral(1,mvaNBinsX,1,mvaNBinsY);
  
  
  TH1D * mvaCC1pi_Purity_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC1pi->Clone("mvaCC1pi_Purity_2tracks");
  TH1D * mvaCC1pi_Efficiency_2tracks = (TH1D*) hMVAMuonVSProtondiscriminant_2tracks_CC1pi->Clone("mvaCC1pi_Efficiency_2tracks");
  double mvaNCC1pi_Total_2tracks=hMVAMuonVSProtondiscriminant_2tracks_CC1pi->Integral(1,mvaNBinsX,1,mvaNBinsY);

  
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

      //double SqrtSB_1track = ErrorBkg*(NAll_1track - NCC0pi_1track) / pow(NCC0pi_1track,3/2);
      double SqrtSB_1track = sqrt( pow(sqrt(NAll_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_1track-NCC0pi_1track)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
      if(NCC0pi_1track!=0) SqrtSB_1track/=(NCC0pi_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
      mvaCC0pi_SignalBkg_1track->SetBinContent(ibinx,ibiny,SqrtSB_1track);
      
      
    //double SqrtSB_1track = TMath::Sqrt(NAll_1track)*(NAll_1track-NCC0pi_1track);//sqrt(signal+bkg)*bkg
    //if(NCC0pi_1track!=0) SqrtSB_1track/=pow(NCC0pi_1track,2.);//sqrt(signal+bkg)*bkg / signal*signal
     
      for(int ibiny=1;ibiny<=mvaNBinsY;ibiny++){
	//CC0pi
	double NCC0pi_2tracks=hMVAMuonVSProtondiscriminant_2tracks_CC0pi->Integral(ibinx,mvaNBinsX,ibiny,mvaNBinsY);
	double NAll_2tracks=hMVAMuonVSProtondiscriminant_2tracks_All->Integral(ibinx,mvaNBinsX,ibiny,mvaNBinsY);
	
	double EffCC0pi_2tracks=NCC0pi_2tracks;
	double PurCC0pi_2tracks=NCC0pi_2tracks;

	if(mvaNCC0pi_Total_2tracks!=0) EffCC0pi_2tracks/=mvaNCC0pi_Total_2tracks;
	if(NAll_2tracks!=0) PurCC0pi_2tracks/=NAll_2tracks;
	mvaCC0pi_Purity_2tracks->SetBinContent(ibinx,ibiny,PurCC0pi_2tracks);
	mvaCC0pi_Efficiency_2tracks->SetBinContent(ibinx,ibiny,EffCC0pi_2tracks);

	//double SqrtSB_2tracks = ErrorBkg*(NAll_2tracks - NCC0pi_2tracks) / pow(NCC0pi_2tracks,3/2);
	double SqrtSB_2tracks = sqrt( pow(sqrt(NAll_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_2tracks-NCC0pi_2tracks)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
	if(NCC0pi_2tracks!=0) SqrtSB_2tracks/=(NCC0pi_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
	mvaCC0pi_SignalBkg_2tracks->SetBinContent(ibinx,ibiny,SqrtSB_2tracks);

	//CC1pi
	double NCC1pi_2tracks=hMVAMuonVSProtondiscriminant_2tracks_CC1pi->Integral(ibinx,mvaNBinsX,1,ibiny-1);
	double NAll_2tracks2=hMVAMuonVSProtondiscriminant_2tracks_All->Integral(ibinx,mvaNBinsX,1,ibiny-1);
	
	double EffCC1pi_2tracks=NCC1pi_2tracks;
	double PurCC1pi_2tracks=NCC1pi_2tracks;

	if(mvaNCC1pi_Total_2tracks!=0) EffCC1pi_2tracks/=mvaNCC1pi_Total_2tracks;
	if(NAll_2tracks2!=0) PurCC1pi_2tracks/=NAll_2tracks2;
	mvaCC1pi_Purity_2tracks->SetBinContent(ibinx,ibiny,PurCC1pi_2tracks);
	mvaCC1pi_Efficiency_2tracks->SetBinContent(ibinx,ibiny,EffCC1pi_2tracks);

      }
  }

  TCanvas * cmvaEfficiency_1track = new TCanvas("cmvaEfficiency_1track","1 track,Efficiency/Purity of CC0pi with MVA discriminant cut value");//Use the 1 track sample
  mvaCC0pi_Efficiency_1track->SetLineColor(kRed);
  mvaCC0pi_Purity_1track->SetLineColor(kBlue);
  mvaCC0pi_Efficiency_1track->Draw();
  mvaCC0pi_Purity_1track->Draw("same");
  mvaCC0pi_Efficiency_1track->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaEfficiency_1track->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_1track_EfficiencyPurity.eps");

  TCanvas * cmvaCC0piSignalBkg_1track = new TCanvas("cmvaCC0piSignalBkg_1track","1 track,SignalBkg of CC0pi with MVA discriminant cut value");//Use the 1 track sample
  mvaCC0pi_SignalBkg_1track->SetLineColor(kBlack);
  mvaCC0pi_SignalBkg_1track->Draw();
  mvaCC0pi_SignalBkg_1track->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaCC0piSignalBkg_1track->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_1track_SignalBkgPurity.eps");

  ////////////////////////////////////////////

  TH1D * mvaMuon_1track = (TH1D*) hMVAMuondiscriminantVSPDG_1track_All->ProjectionY("mvaMuon_1track",1,1);
  TH1D * mvaMuon_1track = (TH1D*) hMVAMuondiscriminantVSPDG_1track_All->ProjectionY("mvaMuon_1track",1,1);
  TH1D * mvaAll_1track = (TH1D*) hMVAMuondiscriminantVSPDG_1track_All->ProjectionY("mvaAll_1track",1,hMVAMuondiscriminantVSPDG_1track_All->GetNbinsX());
  TH1D * mvaAll_1track = (TH1D*) hMVAMuondiscriminantVSPDG_1track_All->ProjectionY("mvaAll_1track",1,hMVAMuondiscriminantVSPDG_1track_All->GetNbinsX());
  
  mvaNBinsX=mvaMuon_1track->GetNbinsX();
  double mvaMuon_Total_1track=mvaMuon_1track->Integral();
  
  TH1D * mvaMuon_Purity_1track = (TH1D*) hMVAMuondiscriminantVSPDG_1track_All->ProjectionY("mvaMuon_Purity_1track",1,1);
  TH1D * mvaMuon_Efficiency_1track = (TH1D*) hMVAMuondiscriminantVSPDG_1track_All->ProjectionY("mvaMuon_Efficiency_1track",1,1);
  mvaMuon_Purity_1track->Reset();
  mvaMuon_Efficiency_1track->Reset();
  
  for(int ibinx=1;ibinx<=mvaNBinsX;ibinx++){
    cout<<"Number of muons="<<mvaMuon_1track->Integral(ibinx,mvaNBinsX)<<", and others="<<mvaAll_1track->Integral(ibinx,mvaNBinsX)<<endl;
    double EffMuon = mvaMuon_1track->Integral(ibinx,mvaNBinsX) / ( mvaMuon_Total_1track!=0 ?   mvaMuon_Total_1track : 1. );
    double PurMuon = mvaMuon_1track->Integral(ibinx,mvaNBinsX) / ( mvaAll_1track->Integral(ibinx,mvaNBinsX)!=0 ?   mvaAll_1track->Integral(ibinx,mvaNBinsX) : 1. );
    mvaMuon_Efficiency_1track->SetBinContent(ibinx,EffMuon);
    mvaMuon_Purity_1track->SetBinContent(ibinx,PurMuon);
  }
  TCanvas * cmvaMuonEfficiency_1track = new TCanvas("cmvaMuonEfficiency_1track","1 track,Efficiency/Purity of Muon with MVA discriminant cut value");//Use the 1 track sample
  mvaMuon_Efficiency_1track->SetLineColor(kRed);
  mvaMuon_Purity_1track->SetLineColor(kBlue);
  mvaMuon_Efficiency_1track->Draw();
  mvaMuon_Purity_1track->Draw("same");
  cmvaMuonEfficiency_1track->SaveAs("../plots/MVAOptimisation/MVA_Muon_1track_EfficiencyPurity.eps");

  //
  
  TCanvas * cmvaEfficiency_2tracks = new TCanvas("cmvaEfficiency_2tracks","2 tracks,Efficiency of CC0pi with MVA discriminant cut value");//Use the 2 tracks sample
  //mvaCC0pi_Efficiency_2tracks->SetLineColor(kRed);
  //mvaCC0pi_Purity_2tracks->SetLineColor(kBlue);
  mvaCC0pi_Efficiency_2tracks->Draw("COLZTEXT");
  mvaCC0pi_Efficiency_2tracks->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  mvaCC0pi_Efficiency_2tracks->GetYaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaEfficiency_2tracks->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_Efficiency.eps");

  TCanvas * cmvaPurity_2tracks = new TCanvas("cmvaPurity_2tracks","2 tracks,Purity of CC0pi with MVA discriminant cut value");//Use the 2 tracks sample
  mvaCC0pi_Purity_2tracks->Draw("COLZTEXT");
  mvaCC0pi_Purity_2tracks->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  mvaCC0pi_Purity_2tracks->GetYaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaPurity_2tracks->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_Purity.eps");

  TCanvas * cmvaSignalBkg_2tracks = new TCanvas("cmvaSignalBkg_2tracks","2 tracks,SignalBkg of CC0pi with MVA discriminant cut value");//Use the 2 tracks sample
  mvaCC0pi_SignalBkg_2tracks->Draw("COLZTEXT");
  mvaCC0pi_SignalBkg_2tracks->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  mvaCC0pi_SignalBkg_2tracks->GetYaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaSignalBkg_2tracks->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_SignalBkg.eps");

  TCanvas * cmvaVSMuonMomentum_1track = new TCanvas("cmvaVSMuonMomentum_1track","1 track,MVA discriminant cut value with true muon momentum");//Use the 1 track sample
  cmvaVSMuonMomentum_1track->Divide(1,2);
  cmvaVSMuonMomentum_1track->cd(1);
  hMVAMuondiscriminantVSMuonMomentum_1track_CC0pi->Draw("COLZ");
  cmvaVSMuonMomentum_1track->cd(2);
  hMVAMuondiscriminantVSMuonMomentum_1track_All->Draw("COLZ");

  TCanvas * cmvaVSDistance_1track = new TCanvas("cmvaVSDistance_1track","1 track,MVA discriminant cut value with true muon momentum");//Use the 1 track sample
  cmvaVSDistance_1track->Divide(1,2);
  cmvaVSDistance_1track->cd(1);
  hMVAMuondiscriminantVSDistance_1track_CC0pi->Draw("COLZ");
  cmvaVSDistance_1track->cd(2);
  hMVAMuondiscriminantVSDistance_1track_All->Draw("COLZ");
  
  TCanvas * cmvaVSDistance_1track_profile = new TCanvas("cmvaVSDistance_1track_profile","1 track,MVA discriminant cut value with true muon momentum");//Use the 1 track sample
  TProfile * pMVAMuondiscriminantVSDistance_1track_CC0pi = (TProfile*) hMVAMuondiscriminantVSDistance_1track_CC0pi->ProfileX("pMVAMuondiscriminantVSDistance_1track_CC0pi",0,hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsY());
  pMVAMuondiscriminantVSDistance_1track_CC0pi->SetLineColor(kRed);
  pMVAMuondiscriminantVSDistance_1track_CC0pi->Draw();
  TProfile * pMVAMuondiscriminantVSDistance_1track_All = (TProfile*) hMVAMuondiscriminantVSDistance_1track_All->ProfileX("pMVAMuondiscriminantVSDistance_1track_All",0,hMVAMuondiscriminantVSDistance_1track_All->GetNbinsY());
  pMVAMuondiscriminantVSDistance_1track_All->SetLineColor(kBlue);
  pMVAMuondiscriminantVSDistance_1track_All->Draw("same");
  pMVAMuondiscriminantVSDistance_1track_CC0pi->GetYaxis()->SetRangeUser(-0.2,0.2);
  cmvaVSDistance_1track_profile->SaveAs("../plots/MVAOptimisation/MVAVSDistance_CC0pi_1track_profile.eps");

  //CC1pi
  TCanvas * cmvaEfficiency_2tracks_CC1pi = new TCanvas("cmvaEfficiency_2tracks_CC1pi","2 tracks,Efficiency of CC0pi with MVA discriminant cut value For CC1pi");//Use the 2 tracks sample
  //mvaCC0pi_Efficiency_2tracks->SetLineColor(kRed);
  //mvaCC0pi_Purity_2tracks->SetLineColor(kBlue);
  mvaCC1pi_Efficiency_2tracks->Draw("COLZTEXT");
  mvaCC1pi_Efficiency_2tracks->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  mvaCC1pi_Efficiency_2tracks->GetYaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaEfficiency_2tracks_CC1pi->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_Efficiency_CC1pi.eps");

  TCanvas * cmvaPurity_2tracks_CC1pi = new TCanvas("cmvaPurity_2tracks_CC1pi","2 tracks,Purity of CC0pi with MVA discriminant cut value for CC1pi");//Use the 2 tracks sample
  //mvaCC0pi_Purity_2tracks->SetLineColor(kRed);
  //mvaCC0pi_Purity_2tracks->SetLineColor(kBlue);
  mvaCC1pi_Purity_2tracks->Draw("COLZTEXT");
  mvaCC1pi_Purity_2tracks->GetXaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  mvaCC1pi_Purity_2tracks->GetYaxis()->SetRangeUser(MVAZoomInf,MVAZoomSup);
  cmvaPurity_2tracks_CC1pi->SaveAs("../plots/MVAOptimisation/MVA_CC0pi_2tracks_Purity_CC1pi.eps");

  /////////
  TH2D * Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi = (TH2D*) hMVAMuondiscriminantVSDistance_1track_CC0pi->Clone("Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi");
  Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi->Reset();
  TH2D * Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi = (TH2D*) hMVAMuondiscriminantVSDistance_1track_CC0pi->Clone("Purity7_hMVAMuondiscriminantVSDistance_1track_CC0pi");
  Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi->Reset();

  TCanvas * cEffPurmvaVSDistance_1track_1D[hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsX()];
  TH1D * Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsX()];
  TH1D * Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsX()];
  double TF1_CC0pi_OverCut=0;
  double TF1_CC0pi=0;
  double TF1_All_OverCut=0;
  double TF1_All=0;
  
  for(int ibinx=1;ibinx<=hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsX();ibinx++){//Run over bins of distance
    double TotalCC0pi=0;
    double TotalAll=0;
    for(int ibiny=1;ibiny<=hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsY();ibiny++){
      TotalCC0pi+=hMVAMuondiscriminantVSDistance_1track_CC0pi->GetBinContent(ibinx,ibiny);
      TotalAll+=hMVAMuondiscriminantVSDistance_1track_All->GetBinContent(ibinx,ibiny);
    }
    //double CutValue=0.004*(hMVAMuondiscriminantVSDistance_2tracks_CC0pi->GetXaxis()->GetBinCenter(ibinx)-17.5);
    //if(hMVAMuondiscriminantVSDistance_2tracks_CC0pi->GetXaxis()->GetBinCenter(ibinx)>40) CutValue=0.1;
    
    for(int ibiny=1;ibiny<=hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsY();ibiny++){
      double CC0pi_OverCut=0;
      double All_OverCut=0;
      for(int ibiny2=ibiny;ibiny2<=hMVAMuondiscriminantVSDistance_1track_CC0pi->GetNbinsY();ibiny2++){
	CC0pi_OverCut+=hMVAMuondiscriminantVSDistance_1track_CC0pi->GetBinContent(ibinx,ibiny2);
	All_OverCut+=hMVAMuondiscriminantVSDistance_1track_All->GetBinContent(ibinx,ibiny2);
      }
      double Efficiency=CC0pi_OverCut / ( TotalCC0pi!=0 ?  TotalCC0pi : 1. );
      double Purity=CC0pi_OverCut / ( All_OverCut!=0 ?  All_OverCut : 1. );

      Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi->SetBinContent(ibinx,ibiny,Efficiency);
      Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi->SetBinContent(ibinx,ibiny,Purity);

      //double MVAValue_CC0pi=hMVAMuondiscriminantVSDistance_2tracks_CC0pi->GetYaxis()->GetBinCenter(ibiny);
      //if(MVAValue_CC0pi>CutValue) TF1_CC0pi_OverCut+=
    }
    cEffPurmvaVSDistance_1track_1D[ibinx] = new TCanvas(Form("cEffPurmvaVSDistance_1track_1D_%d",ibinx-1));
    Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1] = (TH1D*) Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi->ProjectionY(Form("Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D_%d",ibinx-1),ibinx,ibinx);
    Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1] = (TH1D*) Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi->ProjectionY(Form("Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D_%d",ibinx-1),ibinx,ibinx);
    Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1]->SetLineColor(kRed);
    Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1]->SetLineColor(kBlue);
    Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1]->SetLineWidth(2);
    Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1]->SetLineWidth(2);
    Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1]->Draw();
    Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi_1D[ibinx-1]->Draw("same");
    cEffPurmvaVSDistance_1track_1D[ibinx]->SaveAs(Form("cEffPurmvaVSDistance_1track_1D_Distance%3.0f",hMVAMuondiscriminantVSDistance_1track_CC0pi->GetXaxis()->GetBinCenter(ibinx)));
  }
  TCanvas * cEffPurmvaVSDistance_1track = new TCanvas("cEffPurmvaVSDistance_1track","1 track,MVA discriminant cut value with true muon momentum");//Use the 1 track sample
  cEffPurmvaVSDistance_1track->Divide(1,2);
  cEffPurmvaVSDistance_1track->cd(1);
  Efficiency_hMVAMuondiscriminantVSDistance_1track_CC0pi->Draw("COLZ");
  cEffPurmvaVSDistance_1track->cd(2);
  Purity_hMVAMuondiscriminantVSDistance_1track_CC0pi->Draw("COLZ");
  
  TCanvas * cmvaVSDistance_2tracks_LowestMVA = new TCanvas("cmvaVSDistance_2tracks_LowestMVA","1 track,MVA discriminant cut value with true muon momentum");//Use the 1 track sample
  cmvaVSDistance_2tracks_LowestMVA->Divide(1,2);
  cmvaVSDistance_2tracks_LowestMVA->cd(1);
  hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_CC0pi->Draw("COLZ");
  cmvaVSDistance_2tracks_LowestMVA->cd(2);
  hMVAProtondiscriminantVSDistance_2tracks_LowestMVA_All->Draw("COLZ");

  /////////////////////
  



  TH1D * lmvaCC0pi_Purity_1track = (TH1D*) hLMVA_1track_CC0pi->Clone("lmvaCC0pi_Purity_1track");
  TH1D * lmvaCC0pi_Efficiency_1track = (TH1D*) hLMVA_1track_CC0pi->Clone("lmvaCC0pi_Efficiency_1track");
  TH1D * lmvaCC0pi_SignalBkg_1track = (TH1D*) hLMVA_1track_CC0pi->Clone("lmvaCC0pi_SignalBkg_1track");

  TH1D * lmvaCC0pi_Purity_LargestLMVA = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("lmvaCC0pi_Purity_LargestLMVA");
  TH1D * lmvaCC0pi_Efficiency_LargestLMVA = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("lmvaCC0pi_Efficiency_LargestLMVA");
  TH1D * lmvaCC0pi_SignalBkg_LargestLMVA = (TH1D*) hLMVA_1track_CC0pi->Clone("lmvaCC0pi_SignalBkg_LargestLMVA");

  TH1D * lmvaCC0pi_Purity_2tracks = (TH1D*) hLMVA_2tracks_CC0pi->Clone("lmvaCC0pi_Purity_2tracks");
  TH1D * lmvaCC0pi_Efficiency_2tracks = (TH1D*) hLMVA_2tracks_CC0pi->Clone("lmvaCC0pi_Efficiency_2tracks");
  TH1D * lmvaCC0pi_SignalBkg_2tracks = (TH1D*) hLMVA_2tracks_CC0pi->Clone("lmvaCC0pi_SignalBkg_2tracks");

  double lmvaNBinsX=lmvaCC0pi_Efficiency_1track->GetNbinsX();
  cout<<"NBins LMVA = "<<lmvaNBinsX<<endl;
  double lmvaNBinsY=lmvaCC0pi_Efficiency_1track->GetNbinsX();
  double lmvaNCC0pi_Total_1track=hLMVA_1track_CC0pi->Integral(1,lmvaNBinsX);
  double lmvaNCC0pi_Total_LargestLMVA=hLMVA_LargestLMVA_CC0pi->Integral(1,lmvaNBinsX);
  double lmvaNCC0pi_Total_2tracks=hLMVA_2tracks_CC0pi->Integral(1,lmvaNBinsX,1,lmvaNBinsY);
  //TH1D * lmvaCC0pi_Purity_2tracks = (TH1D*) hLMVA_2tracks_CC0pi->Clone("lmvaCC0pi_Purity_2tracks");
  //TH1D * lmvaCC0pi_Efficiency_2tracks = (TH1D*) hLMVA_2tracks_CC0pi->Clone("lmvaCC0pi_Efficiency_2tracks");

  for(int ibinx=1;ibinx<=lmvaNBinsX;ibinx++){
      cout<<"LMVA, "<<ibinx<<endl;    
      double NCC0pi_1track=hLMVA_1track_CC0pi->Integral(ibinx,lmvaNBinsX);
      double NAll_1track=hLMVA_1track_All->Integral(ibinx,lmvaNBinsX);
      double EffCC0pi_1track=NCC0pi_1track;
      double PurCC0pi_1track=NCC0pi_1track;
      if(lmvaNCC0pi_Total_1track!=0) EffCC0pi_1track/=lmvaNCC0pi_Total_1track;
      if(NAll_1track!=0) PurCC0pi_1track/=NAll_1track;
      lmvaCC0pi_Purity_1track->SetBinContent(ibinx,PurCC0pi_1track);
      lmvaCC0pi_Efficiency_1track->SetBinContent(ibinx,EffCC0pi_1track);

      double SqrtSB_1track = sqrt( pow(sqrt(NAll_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_1track-NCC0pi_1track)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
	if(NCC0pi_1track!=0) SqrtSB_1track/=(NCC0pi_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
	//double SqrtSB_1track = ErrorBkg*(NAll_1track - NCC0pi_1track) / pow(NCC0pi_1track,3/2);
    //double SqrtSB_1track = sqrt( pow(sqrt(NAll_1track/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_1track-NCC0pi_1track)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
      lmvaCC0pi_SignalBkg_1track->SetBinContent(ibinx,ibiny,SqrtSB_1track);
      
      double NCC0pi_LargestLMVA=hLMVA_LargestLMVA_CC0pi->Integral(ibinx,lmvaNBinsX);
      double NAll_LargestLMVA=hLMVA_LargestLMVA_All->Integral(ibinx,lmvaNBinsX);
      double EffCC0pi_LargestLMVA=NCC0pi_LargestLMVA;
      double PurCC0pi_LargestLMVA=NCC0pi_LargestLMVA;
      if(lmvaNCC0pi_Total_LargestLMVA!=0) EffCC0pi_LargestLMVA/=lmvaNCC0pi_Total_LargestLMVA;
      if(NAll_LargestLMVA!=0) PurCC0pi_LargestLMVA/=NAll_LargestLMVA;
      lmvaCC0pi_Purity_LargestLMVA->SetBinContent(ibinx,PurCC0pi_LargestLMVA);
      lmvaCC0pi_Efficiency_LargestLMVA->SetBinContent(ibinx,EffCC0pi_LargestLMVA);

      double SqrtSB_LargestLMVA = sqrt( pow(sqrt(NAll_LargestLMVA/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_LargestLMVA-NCC0pi_LargestLMVA)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
	if(NCC0pi_LargestLMVA!=0) SqrtSB_LargestLMVA/=(NCC0pi_LargestLMVA/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
	//double SqrtSB_LargestLMVA = ErrorBkg*(NAll_LargestLMVA - NCC0pi_LargestLMVA) / pow(NCC0pi_LargestLMVA,3/2);
      //double SqrtSB_LargestLMVA = sqrt( pow(sqrt(NAll_LargestLMVA/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_LargestLMVA-NCC0pi_LargestLMVA)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
      lmvaCC0pi_SignalBkg_LargestLMVA->SetBinContent(ibinx,ibiny,SqrtSB_LargestLMVA);

      for(int ibiny=1;ibiny<=lmvaNBinsY;ibiny++){
	//CC0pi
	double NCC0pi_2tracks=hLMVA_2tracks_CC0pi->Integral(ibinx,lmvaNBinsX,1,ibiny);
	double NAll_2tracks=hLMVA_2tracks_All->Integral(ibinx,lmvaNBinsX,1,ibiny);
	
	double EffCC0pi_2tracks=NCC0pi_2tracks;
	double PurCC0pi_2tracks=NCC0pi_2tracks;

	if(lmvaNCC0pi_Total_2tracks!=0) EffCC0pi_2tracks/=lmvaNCC0pi_Total_2tracks;
	if(NAll_2tracks!=0) PurCC0pi_2tracks/=NAll_2tracks;
	lmvaCC0pi_Purity_2tracks->SetBinContent(ibinx,ibiny,PurCC0pi_2tracks);
	lmvaCC0pi_Efficiency_2tracks->SetBinContent(ibinx,ibiny,EffCC0pi_2tracks);

	double SqrtSB_2tracks = sqrt( pow(sqrt(NAll_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_2tracks-NCC0pi_2tracks)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
	if(NCC0pi_2tracks!=0) SqrtSB_2tracks/=(NCC0pi_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle));//sqrt(signal+bkg)*bkg / signal*signal
	//double SqrtSB_2tracks = ErrorBkg*(NAll_2tracks - NCC0pi_2tracks) / pow(NCC0pi_2tracks,3/2);
	//double SqrtSB_2tracks = sqrt( pow(sqrt(NAll_2tracks/TMath::Max(NBinsTrueMom,NBinsTrueAngle)),2) + pow(ErrorBkg*(NAll_2tracks-NCC0pi_2tracks)/TMath::Max(NBinsTrueMom,NBinsTrueAngle),2));//sqrt( sqrt(signal+bkg)**2+bkg**2). It is divided by the number of bins along momentum or angle in order to have an optimization for an average bin, not for the whole distribution
	lmvaCC0pi_SignalBkg_2tracks->SetBinContent(ibinx,ibiny,SqrtSB_2tracks);
      }
  }

  TCanvas * clmvaEfficiency_1track = new TCanvas("clmvaEfficiency_1track","1 track,Efficiency/Purity of CC0pi with LMVA discriminant cut value");//Use the 1 track sample
  lmvaCC0pi_Efficiency_1track->SetLineColor(kRed);
  lmvaCC0pi_Purity_1track->SetLineColor(kBlue);
  lmvaCC0pi_Efficiency_1track->Draw();
  lmvaCC0pi_Purity_1track->Draw("same");
  lmvaCC0pi_Efficiency_1track->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaEfficiency_1track->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_1track_EfficiencyPurity.eps");

  TCanvas * clmvaCC0piSignalBkg_1track = new TCanvas("clmvaCC0piSignalBkg_1track","1 track,SignalBkg of CC0pi with LMVA discriminant cut value");//Use the 1 track sample
  lmvaCC0pi_SignalBkg_1track->SetLineColor(kBlack);
  lmvaCC0pi_SignalBkg_1track->Draw();
  lmvaCC0pi_SignalBkg_1track->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaCC0piSignalBkg_1track->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_1track_SignalBkgPurity.eps");

  TCanvas * clmvaEfficiency_2tracks = new TCanvas("clmvaEfficiency_2tracks","2 tracks,Efficiency of CC0pi with LMVA discriminant cut value");//Use the 2 tracks sample
  //lmvaCC0pi_Efficiency_2tracks->SetLineColor(kRed);
  //lmvaCC0pi_Purity_2tracks->SetLineColor(kBlue);
  lmvaCC0pi_Efficiency_2tracks->Draw("COLZTEXT");
  lmvaCC0pi_Efficiency_2tracks->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  lmvaCC0pi_Efficiency_2tracks->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaEfficiency_2tracks->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_2tracks_Efficiency.eps");

  TCanvas * clmvaPurity_2tracks = new TCanvas("clmvaPurity_2tracks","2 tracks,Purity of CC0pi with LMVA discriminant cut value");//Use the 2 tracks sample
  lmvaCC0pi_Purity_2tracks->Draw("COLZTEXT");
  lmvaCC0pi_Purity_2tracks->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  lmvaCC0pi_Purity_2tracks->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaPurity_2tracks->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_2tracks_Purity.eps");

  TCanvas * clmvaSignalBkg_2tracks = new TCanvas("clmvaSignalBkg_2tracks","2 tracks,SignalBkg of CC0pi with LMVA discriminant cut value");//Use the 2 tracks sample
  lmvaCC0pi_SignalBkg_2tracks->Draw("COLZTEXT");
  lmvaCC0pi_SignalBkg_2tracks->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  lmvaCC0pi_SignalBkg_2tracks->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaSignalBkg_2tracks->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_2tracks_SignalBkg.eps");

  TCanvas * clmvaEfficiency_LargestLMVA = new TCanvas("clmvaEfficiency_LargestLMVA","Largest MVA,Efficiency/Purity of CC0pi with LMVA discriminant cut value");//Use the 1 track sample
  lmvaCC0pi_Efficiency_LargestLMVA->SetLineColor(kRed);
  lmvaCC0pi_Purity_LargestLMVA->SetLineColor(kBlue);
  lmvaCC0pi_Efficiency_LargestLMVA->Draw();
  lmvaCC0pi_Purity_LargestLMVA->Draw("same");
  lmvaCC0pi_Efficiency_LargestLMVA->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaEfficiency_LargestLMVA->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_LargestLMVA_EfficiencyPurity.eps");

  TCanvas * clmvaCC0piSignalBkg_LargestLMVA = new TCanvas("clmvaCC0piSignalBkg_LargestLMVA","1 track,SignalBkg of CC0pi with LMVA discriminant cut value");//Use the 1 track sample
  lmvaCC0pi_SignalBkg_LargestLMVA->SetLineColor(kBlack);
  lmvaCC0pi_SignalBkg_LargestLMVA->Draw();
  lmvaCC0pi_SignalBkg_LargestLMVA->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  clmvaCC0piSignalBkg_LargestLMVA->SaveAs("../plots/MVAOptimisation/LMVA_CC0pi_LargestMVA_SignalBkgPurity.eps");

  
#endif


  TH1D * TrueParticleType_1track[NFSIs];
  TH2D * TrueParticleType_2tracks[NFSIs];
  TH1D * TrueParticleType_1track_CC0pi = (TH1D*) PIDCutTuning->Get("TrueParticleType_1track0");
  TH2D * TrueParticleType_2tracks_CC0pi = (TH2D*) PIDCutTuning->Get("TrueParticleType_2tracks0");
  TH1D * TrueParticleType_1track_All = (TH1D*) TrueParticleType_1track_CC0pi->Clone("TrueParticleType_1track_All");
  TH2D * TrueParticleType_2tracks_All = (TH2D*) TrueParticleType_2tracks_CC0pi->Clone("TrueParticleType_2tracks_All");
  TrueParticleType_1track_CC0pi->SetFillStyle(0);
  TrueParticleType_2tracks_CC0pi->SetFillStyle(0);
  TrueParticleType_1track_All->SetFillStyle(0);
  TrueParticleType_2tracks_All->SetFillStyle(0);

  
  for(int fsi=1;fsi<NFSIs;fsi++){
    sprintf(Name,"TrueParticleType_1track%d",fsi);
    TrueParticleType_1track[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    TrueParticleType_1track[fsi]->SetFillStyle(0);

    sprintf(Name,"TrueParticleType_2tracks%d",fsi);
    TrueParticleType_2tracks[fsi] = (TH2D*) PIDCutTuning->Get(Name);
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

 TDirectory * PIDCutTuning = (TDirectory*) fMC->Get("PIDCutTuning");
 PIDCutTuning->cd();

 //#define DEDXPLOTS
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
   cout<<"Sample tested="<<is<<endl;
   hEnergyDepositionLength_Muon_1D[is] = (TH1D*) MVAInputVariables->Get(Form("EnergyDepositionLength_Muon_%d_y",is));
   hEnergyDepositionLength_Pion_1D[is] = (TH1D*) MVAInputVariables->Get(Form("EnergyDepositionLength_Pion_%d_y",is));
   hEnergyDepositionLength_Proton_1D[is] = (TH1D*) MVAInputVariables->Get(Form("EnergyDepositionLength_Proton_%d_y",is));

   hEnergyDepositionSplineLength_Muon_1D[is] = (TH1D*) MVAInputVariables->Get(Form("EnergyDepositionSplineLength_Muon_%d_y",is));
   hEnergyDepositionSplineLength_Pion_1D[is] = (TH1D*) MVAInputVariables->Get(Form("EnergyDepositionSplineLength_Pion_%d_y",is));
   hEnergyDepositionSplineLength_Proton_1D[is] = (TH1D*) MVAInputVariables->Get(Form("EnergyDepositionSplineLength_Proton_%d_y",is));
   
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

     double RatioPionMuon = ValuePion/ ( ValueMuon!=0 ?  ValueMuon : 1. );
     double ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorPion/ ( ValuePion!=0 ?  ValuePion : 1. )),2.));
     double RatioProtonMuon = ValueProton/ ( ValueMuon!=0 ?  ValueMuon : 1. );
     double ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorProton/ ( ValueProton!=0 ?  ValueProton : 1. )),2.));

     RatioPionMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_EnergyDepositionLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_EnergyDepositionLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     double NumberOfSigmaPionMuon = RatioPionMuon/ ( ErrorPionMuon!=0 ? ErrorPionMuon : 1. ) ;
     double NumberOfSigmaProtonMuon = RatioProtonMuon/ ( ErrorProtonMuon!=0 ? ErrorProtonMuon : 1. ) ;

     NumberOfSigmaPionMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);
     
     //Spline
     ValueMuon = hEnergyDepositionSplineLength_Muon_1D[is]->GetBinContent(ibinx);
     ErrorMuon = hEnergyDepositionSplineLength_Muon_1D[is]->GetBinError(ibinx);
     ValuePion = hEnergyDepositionSplineLength_Pion_1D[is]->GetBinContent(ibinx);
     ErrorPion = hEnergyDepositionSplineLength_Pion_1D[is]->GetBinError(ibinx);
     ValueProton = hEnergyDepositionSplineLength_Proton_1D[is]->GetBinContent(ibinx);
     ErrorProton = hEnergyDepositionSplineLength_Proton_1D[is]->GetBinError(ibinx);
     
     RatioPionMuon = ValuePion/ ( ValueMuon!=0 ?  ValueMuon : 1. ) ;
     ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorPion/ ( ValuePion!=0 ?  ValuePion : 1. )),2.));
     RatioProtonMuon = ValueProton/ ( ValueMuon!=0 ?  ValueMuon : 1. ) ;
     ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorProton/ ( ValueProton!=0 ?  ValueProton : 1. )),2.));

     RatioPionMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_EnergyDepositionSplineLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_EnergyDepositionSplineLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     NumberOfSigmaPionMuon = RatioPionMuon/ ( ErrorPionMuon!=0 ? ErrorPionMuon : 1. ) ;
     NumberOfSigmaProtonMuon = RatioProtonMuon/ ( ErrorProtonMuon!=0 ? ErrorProtonMuon : 1. ) ;

     NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);

     //cout<<"ibin="<<ibinx<<", Ratio="<<RatioPionMuon<<" +- "<<ErrorPionMuon<<endl;
     cout<<"ibin="<<ibinx<<": Ratio="<<RatioPionMuon_EnergyDepositionLength[is]->GetBinContent(ibinx)<<" +- "<<RatioPionMuon_EnergyDepositionLength[is]->GetBinError(ibinx)<<endl;
     cout<<"Error pion="<<ErrorPion<<", error muon="<<ErrorMuon<<endl;
     cout<<" Spline, Ratio="<<RatioPionMuon_EnergyDepositionSplineLength[is]->GetBinContent(ibinx)<<" +- "<<RatioPionMuon_EnergyDepositionSplineLength[is]->GetBinError(ibinx)<<endl;

     
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
   hEnergyDepositionLength_Muon_1D[is]->GetYaxis()->SetRangeUser(20,100);
   hEnergyDepositionLength_Pion_1D[is]->Draw("HISTsame");
   hEnergyDepositionLength_Proton_1D[is]->Draw("HISTsame");
#ifdef INTERPOLATION
   hEnergyDepositionSplineLength_Muon_1D[is]->Draw("HISTsame");
   hEnergyDepositionSplineLength_Pion_1D[is]->Draw("HISTsame");
   hEnergyDepositionSplineLength_Proton_1D[is]->Draw("HISTsame");
#endif

   cRatio_EnergyDepositionLength[is]->cd(2);
   RatioPionMuon_EnergyDepositionLength[is]->SetLineColor(kGreen+2);
   RatioPionMuon_EnergyDepositionSplineLength[is]->SetLineColor(kGreen+2);
   RatioProtonMuon_EnergyDepositionLength[is]->SetLineColor(kRed);
   RatioProtonMuon_EnergyDepositionSplineLength[is]->SetLineColor(kRed);

   RatioPionMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   RatioProtonMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   
   RatioPionMuon_EnergyDepositionLength[is]->Draw("HIST");
   RatioProtonMuon_EnergyDepositionLength[is]->Draw("HISTsame");   
   RatioPionMuon_EnergyDepositionLength[is]->GetYaxis()->SetRangeUser(0.8,3.);
   
#ifdef INTERPOLATION
   RatioPionMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");
   RatioProtonMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");   
#endif 
   
   cRatio_EnergyDepositionLength[is]->cd(3);
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->SetLineColor(kGreen+2);
   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->SetLineColor(kGreen+2);
   NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->SetLineColor(kRed);
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->SetLineColor(kRed);

   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->SetLineStyle(2);
   
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->Draw("HIST");   
   NumberOfSigmaProtonMuon_EnergyDepositionLength[is]->Draw("HISTsame");   
   NumberOfSigmaPionMuon_EnergyDepositionLength[is]->GetYaxis()->SetRangeUser(1.,4.);
#ifdef INTERPOLATION
   NumberOfSigmaPionMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");   
   NumberOfSigmaProtonMuon_EnergyDepositionSplineLength[is]->Draw("HISTsame");
#endif
   cRatio_EnergyDepositionLength[is]->SaveAs(Form("../plots/MVAOptimisation/cRatio_EnergyDepositionLength_sample%d.eps",is));
 }				 
#endif

 
 //#define TRANSVERSEPLOTS
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
   hTransverseWidthLength_Muon_1D[is] = (TH1D*) MVAInputVariables->Get(Form("TransverseWidthLength_Muon_%d_y",is));
   hTransverseWidthLength_Pion_1D[is] = (TH1D*) MVAInputVariables->Get(Form("TransverseWidthLength_Pion_%d_y",is));
   hTransverseWidthLength_Proton_1D[is] = (TH1D*) MVAInputVariables->Get(Form("TransverseWidthLength_Proton_%d_y",is));

   hTransverseWidthNonIsolatedLength_Muon_1D[is] = (TH1D*) MVAInputVariables->Get(Form("TransverseWidthNonIsolatedLength_Muon_%d_y",is));
   hTransverseWidthNonIsolatedLength_Pion_1D[is] = (TH1D*) MVAInputVariables->Get(Form("TransverseWidthNonIsolatedLength_Pion_%d_y",is));
   hTransverseWidthNonIsolatedLength_Proton_1D[is] = (TH1D*) MVAInputVariables->Get(Form("TransverseWidthNonIsolatedLength_Proton_%d_y",is));
   
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

     double RatioPionMuon = ValuePion/ ( ValueMuon!=0 ?  ValueMuon : 1. ) ;
     double ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorPion/ ( ValuePion!=0 ?  ValuePion : 1. )),2.));
     double RatioProtonMuon = ValueProton/ ( ValueMuon!=0 ?  ValueMuon : 1. ) ;
     double ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorProton/ ( ValueProton!=0 ?  ValueProton : 1. )),2.));

     RatioPionMuon_TransverseWidthLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_TransverseWidthLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_TransverseWidthLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_TransverseWidthLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     double NumberOfSigmaPionMuon = RatioPionMuon/ ( ErrorPionMuon!=0 ? ErrorPionMuon : 1. ) ;
     double NumberOfSigmaProtonMuon = RatioProtonMuon/ ( ErrorProtonMuon!=0 ? ErrorProtonMuon : 1. ) ;

     NumberOfSigmaPionMuon_TransverseWidthLength[is]->SetBinContent(ibinx,NumberOfSigmaPionMuon);
     NumberOfSigmaProtonMuon_TransverseWidthLength[is]->SetBinContent(ibinx,NumberOfSigmaProtonMuon);
     
     //NonIsolated
     ValueMuon = hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetBinContent(ibinx);
     ErrorMuon = hTransverseWidthNonIsolatedLength_Muon_1D[is]->GetBinError(ibinx);
     ValuePion = hTransverseWidthNonIsolatedLength_Pion_1D[is]->GetBinContent(ibinx);
     ErrorPion = hTransverseWidthNonIsolatedLength_Pion_1D[is]->GetBinError(ibinx);
     ValueProton = hTransverseWidthNonIsolatedLength_Proton_1D[is]->GetBinContent(ibinx);
     ErrorProton = hTransverseWidthNonIsolatedLength_Proton_1D[is]->GetBinError(ibinx);
     
     RatioPionMuon = ValuePion/ ( ValueMuon!=0 ?  ValueMuon : 1. ) ;
     ErrorPionMuon = RatioPionMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorPion/ ( ValuePion!=0 ?  ValuePion : 1. )),2.));
     RatioProtonMuon = ValueProton/ ( ValueMuon!=0 ?  ValueMuon : 1. ) ;
     ErrorProtonMuon = RatioProtonMuon*TMath::Sqrt(pow((ErrorMuon/ ( ValueMuon!=0 ?  ValueMuon : 1. )),2.)+pow((ErrorProton/ ( ValueProton!=0 ?  ValueProton : 1. )),2.));

     RatioPionMuon_TransverseWidthNonIsolatedLength[is]->SetBinContent(ibinx,RatioPionMuon);
     RatioPionMuon_TransverseWidthNonIsolatedLength[is]->SetBinError(ibinx,ErrorPionMuon);
     RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->SetBinContent(ibinx,RatioProtonMuon);
     RatioProtonMuon_TransverseWidthNonIsolatedLength[is]->SetBinError(ibinx,ErrorProtonMuon);

     NumberOfSigmaPionMuon = RatioPionMuon/ ( ErrorPionMuon!=0 ? ErrorPionMuon : 1. ) ;
     NumberOfSigmaProtonMuon = RatioProtonMuon/ ( ErrorProtonMuon!=0 ? ErrorProtonMuon : 1. ) ;

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

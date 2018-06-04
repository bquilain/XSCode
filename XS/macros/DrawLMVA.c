//For upper cut: Allow region > cut value
const int NSimplifiedPDG=4;
const int NSamples=6;
const int NFSIs=13;
double NBinsUnfolding = 1;//Divide the number of events per the number of bin to have the average number of event per bin. Useful to tune the stat error per bin and not on the whole sample!
double LMVAZoomInf=-1;
double LMVAZoomSup=1;
int RebinValueLMVA=4;
double ErrorBkg=0.15;
double ErrorFlux=0.;

void EfficiencyPurity1D(TH1D * hEfficiency_CC0pi, TH1D * hPurity_CC0pi, TH1D * hError_CC0pi, TH1D * CL_CC0pi, TH1D * CL_All, double ErrorBkg, bool UpCut = true){

  int NBinsX = CL_CC0pi->GetNbinsX();

  double Total_CC0pi = CL_CC0pi->Integral(1,NBinsX);
  double Total_All = CL_All->Integral(1,NBinsX);
  
  for(int ibinx=1;ibinx<=NBinsX;ibinx++){
    
    double Selected_CC0pi = UpCut? CL_CC0pi->Integral(ibinx,NBinsX):CL_CC0pi->Integral(1,ibinx-1);
    double Selected_All = UpCut? CL_All->Integral(ibinx,NBinsX):CL_All->Integral(1,ibinx-1);
      
    double Efficiency_CC0pi = (Total_CC0pi != 0)? Selected_CC0pi/Total_CC0pi:0;
    double Purity_CC0pi = (Selected_All != 0)? Selected_CC0pi/Selected_All:0;

    hEfficiency_CC0pi->SetBinContent(ibinx,Efficiency_CC0pi);
    hPurity_CC0pi->SetBinContent(ibinx,Purity_CC0pi);

    double Selected_All_PerBin = Selected_All/NBinsUnfolding;
    double TotalError;
    if(Selected_All_PerBin == 0) TotalError = 1;
    else TotalError = TMath::Sqrt( pow(TMath::Sqrt(Selected_All_PerBin)/Selected_All_PerBin,2) + pow(ErrorBkg*(Selected_All-Selected_CC0pi)/Selected_All,2) + pow(ErrorFlux,2));//sqrt( sqrt(signal+bkg)**2+bkg**2).
    hError_CC0pi->SetBinContent(ibinx,TotalError);

    //cout<<"MVA : "<<CL_CC0pi->GetBinCenter(ibinx)<<", signal selected : "<<Selected_CC0pi<<", efficiency : "<<Efficiency_CC0pi<<", purity: "<<Purity_CC0pi<<endl;
  }
  //hError_CC0pi->Scale(0.5*hEr);
}

//Assume a higher cut on x axis and lower on y axis
void EfficiencyPurity2D(TH2D * hEfficiency_CC0pi, TH2D * hPurity_CC0pi, TH2D * hError_CC0pi, TH2D * CL_CC0pi, TH2D * CL_All, double ErrorBkg, bool UpCutX = true, bool UpCutY = false){

  int NBinsX = CL_CC0pi->GetNbinsX();
  int NBinsY = CL_CC0pi->GetNbinsY();
  double Total_CC0pi = CL_CC0pi->Integral(1,NBinsX,1,NBinsY);
  double Total_All = CL_All->Integral(1,NBinsX,1,NBinsY);

  for(int ibinx=1;ibinx<=NBinsX;ibinx++){
    for(int ibiny=1;ibiny<=NBinsY;ibiny++){
    
      double Selected_CC0pi;
      if(UpCutX && UpCutY) Selected_CC0pi = CL_CC0pi->Integral(ibinx,NBinsX,ibiny,NBinsY); 
      else if(UpCutX && !UpCutY) Selected_CC0pi = CL_CC0pi->Integral(ibinx,NBinsX,1,ibiny-1);
      else if(!UpCutX && UpCutY) Selected_CC0pi = CL_CC0pi->Integral(1,ibinx-1,ibiny,NBinsY);
      else if(!UpCutX && !UpCutY) Selected_CC0pi = CL_CC0pi->Integral(1,ibinx-1,1,ibiny-1);

      double Selected_All;
      if(UpCutX && UpCutY) Selected_All = CL_All->Integral(ibinx,NBinsX,ibiny,NBinsY); 
      else if(UpCutX && !UpCutY) Selected_All = CL_All->Integral(ibinx,NBinsX,1,ibiny-1);
      else if(!UpCutX && UpCutY) Selected_All = CL_All->Integral(1,ibinx-1,ibiny,NBinsY);
      else if(!UpCutX && !UpCutY) Selected_All = CL_All->Integral(1,ibinx-1,1,ibiny-1);
            
      double Efficiency_CC0pi = (Total_CC0pi != 0)? Selected_CC0pi/Total_CC0pi:0;
      double Purity_CC0pi = (Selected_All != 0)? Selected_CC0pi/Selected_All:0;

      hEfficiency_CC0pi->SetBinContent(ibinx,ibiny,Efficiency_CC0pi);
      hPurity_CC0pi->SetBinContent(ibinx,ibiny,Purity_CC0pi);

    double Selected_All_PerBin = Selected_All/NBinsUnfolding;
    double TotalError;
    if(Selected_All_PerBin == 0) TotalError = 1;
    else TotalError = TMath::Sqrt( pow(TMath::Sqrt(Selected_All_PerBin)/Selected_All_PerBin,2) + pow(ErrorBkg*(Selected_All-Selected_CC0pi)/Selected_All,2) + pow(ErrorFlux,2));//sqrt( sqrt(signal+bkg)**2+bkg**2).
    hError_CC0pi->SetBinContent(ibinx,ibiny,TotalError);
    }
  }
}

void PrepareCanvas1D(TCanvas * c, TH1D * hEfficiency, TH1D * hPurity, TH1D * hError, TLegend * l){
  c->cd();
  hEfficiency->SetLineColor(kRed); hEfficiency->GetYaxis()->SetTitle("Efficiency");
  hPurity->SetLineColor(kBlue); hPurity->GetYaxis()->SetTitle("Purity");
  hError->SetLineColor(kBlack); hError->GetYaxis()->SetTitle("Error");
  hEfficiency->Draw(); hPurity->Draw("same");
  hEfficiency->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  c->Update();
  
    //scale hint1 to the pad coordinates
  hError->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  Float_t rightmax = 1.1*hError->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;
  hError->Scale(scale);
  hError->Draw("same");
  
  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
  //TGaxis *axis = new TGaxis(LMVAZoomSup,gPad->GetUymin(),LMVAZoomSup, gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.04);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.04);
  axis->SetLineColor(kBlack);
  axis->SetTextColor(kBlack);
  axis->SetTitle("Error");
  axis->Draw();
  //hEfficiency->GetYaxis()->SetRangeUser(0.05,1);
  //c->SetLogy();
  l->Draw("same");
}


void ProduceStackParticles(TH1D * hmu, TH1D * hpi, TH1D * hp, TH1D * hothers, THStack * hStack){

  hmu->GetYaxis()->SetTitleOffset(1.3);

  hmu->SetFillColor(kBlue);
  //hmu->SetFillStyle(3001);
  hpi->SetFillColor(kGreen+2);
  //hpi->SetFillStyle(3001);
  hp->SetFillColor(kRed);
  //hp->SetFillStyle(3001);
  hothers->SetFillColor(kGray);
  //hothers->SetFillStyle(3001);

  hmu->SetLineColor(1);
  hmu->SetLineWidth(2);
  hStack->Add(hmu);
  hpi->SetLineColor(1);
  hpi->SetLineWidth(2);
  hStack->Add(hpi);
  hp->SetLineColor(1);
  hp->SetLineWidth(2);
  hStack->Add(hp);
  hothers->SetLineColor(1);
  hothers->SetLineWidth(2);
  hStack->Add(hothers);

  hmu->SetTitle("#mu");
  hpi->SetTitle("#pi");
  hp->SetTitle("p");
  hothers->SetTitle("Others");

}


void DrawLMVA(){
  //gROOT->SetBatch(kTRUE);//Silent mode
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPaintTextFormat("2.2g");//Ste the text format in the option draw("colztext")
  gStyle->SetOptTitle(kFALSE);
 
  
  //TFile * fMC = new TFile("../files/MCSelected_PM_PG_Systematics0_0PM.root");
  //TFile * fMC = new TFile("../files/MCSelected_PM_Systematics0_0.root");
  //TFile * fMC = new TFile("../files/MCSelected_WM_Systematics0_0_StoppedOnly.root");
  TFile * fMC = new TFile("Test.root");
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
  //#define CLParticles
#ifdef CLParticles
  TH1D * MuCL_TrueMuon = (TH1D*) PIDCutTuning->Get("hMuCL_TrueMuon");MuCL_TrueMuon->Sumw2();MuCL_TrueMuon->Scale(1./MuCL_TrueMuon->Integral());
  TH1D * MuCL_TruePion = (TH1D*) PIDCutTuning->Get("hMuCL_TruePion");MuCL_TruePion->Sumw2();MuCL_TrueProton->Scale(1./MuCL_TrueProton->Integral());
  TH1D * MuCL_TrueProton = (TH1D*) PIDCutTuning->Get("hMuCL_TrueProton");MuCL_TrueProton->Sumw2();MuCL_TruePion->Scale(1./MuCL_TruePion->Integral());

  //1. Confidence level for true particles
  TCanvas * cCL = new TCanvas();
  MuCL_TrueMuon->SetLineColor(kBlue);
  MuCL_TruePion->SetLineColor(kGreen);
  MuCL_TrueProton->SetLineColor(kRed);
  MuCL_TrueMuon->Draw("");
  MuCL_TruePion->Draw("same");
  MuCL_TrueProton->Draw("same");
  cCL->SaveAs("../plots/MVAOptimisation/MuCL_particles.eps");
  
  //2. Efficiency and purity for true particles wrt Muon CL
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


  //3. Efficiency and purity for true particles wrt Proton CL
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
  Proton_Purity->Draw();
  Proton_Efficiency->Draw("same");


  //4. MVA vs penetration distance
  TH2D * hMVAMuondiscriminantVSDistance_TrueParticle[NSimplifiedPDG];
  for(int i=0;i<NSimplifiedPDG;i++){
    hMVAMuondiscriminantVSDistance_TrueParticle[i] = (TH2D*) PIDCutTuning->Get(Form("hMVAMuondiscriminantVSDistance_TrueParticle%d",i));
  }
  TCanvas * cMVAMuondiscriminantVSDistance_TrueParticle = new TCanvas("cMVAMuondiscriminantVSDistance_TrueParticle","True particle vs distance and MVA");
  hMVAMuondiscriminantVSDistance_TrueParticle[0]->Draw("colz");
  cMVAMuondiscriminantVSDistance_TrueParticle->SaveAs("../plots/MVAOptimisation/MVAMuondiscriminantVSDistance_TrueParticle.eps");

  
  //5. Purity and Efficiency of MVA wrt penetrating distance
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
  char Name[256];

#define MVA  
#ifdef MVA

  //New LMVA discriminator
  cout<<"New LMVA discriminator"<<endl;

  TH1D * hLMVA_1track[NFSIs];
  TH1D * hLMVA_1track_CC0pi = (TH1D*) PIDCutTuning->Get("hLMVA_1track0");hLMVA_1track_CC0pi->SetFillStyle(0);
  TH1D * hLMVA_1track_All = (TH1D*) hLMVA_1track_CC0pi->Clone("hLMVA_1track_All");
  
  TH2D * hLMVA_2tracks[NFSIs];
  TH2D * hLMVA_2tracks_CC0pi = (TH2D*) PIDCutTuning->Get("hLMVA_2tracks0");hLMVA_2tracks_CC0pi->SetFillStyle(0);
  TH2D * hLMVA_2tracks_All = (TH2D*) hLMVA_2tracks_CC0pi->Clone("hLMVA_2tracks_All");

  TH1D * hLMVA_LargestLMVA[NFSIs];
  TH1D * hLMVA_LargestLMVA_CC0pi = (TH1D*) PIDCutTuning->Get("hLMVA_LargestLMVA0");hLMVA_LargestLMVA_CC0pi->SetFillStyle(0);
  TH1D * hLMVA_LargestLMVA_All = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("hLMVA_LargestLMVA_All");

  TH1D * hLMVA_SmallestLMVA[NFSIs];
  TH1D * hLMVA_SmallestLMVA_CC0pi = (TH1D*) PIDCutTuning->Get("hLMVA_SmallestLMVA0");hLMVA_SmallestLMVA_CC0pi->SetFillStyle(0);
  TH1D * hLMVA_SmallestLMVA_All = (TH1D*) hLMVA_SmallestLMVA_CC0pi->Clone("hLMVA_SmallestLMVA_All");

  TH2D * hLMVA_2tracks_CC1pi = (TH2D*) PIDCutTuning->Get("hLMVA_2tracks3");hLMVA_2tracks_CC1pi->SetFillStyle(0);

  for(int fsi=1;fsi<NFSIs;fsi++){
    sprintf(Name,"hLMVA_1track%d",fsi);
    hLMVA_1track[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hLMVA_1track[fsi]->SetFillStyle(0);

    sprintf(Name,"hLMVA_2tracks%d",fsi);
    hLMVA_2tracks[fsi] = (TH2D*) PIDCutTuning->Get(Name);
    hLMVA_2tracks[fsi]->SetFillStyle(0);

    sprintf(Name,"hLMVA_LargestLMVA%d",fsi);
    hLMVA_LargestLMVA[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hLMVA_LargestLMVA[fsi]->SetFillStyle(0);

    sprintf(Name,"hLMVA_SmallestLMVA%d",fsi);
    hLMVA_SmallestLMVA[fsi] = (TH1D*) PIDCutTuning->Get(Name);
    hLMVA_SmallestLMVA[fsi]->SetFillStyle(0);
            
    if(fsi>-1){
      if(fsi<3){
	hLMVA_1track_CC0pi->Add(hLMVA_1track[fsi]);
	hLMVA_2tracks_CC0pi->Add(hLMVA_2tracks[fsi]);
	hLMVA_LargestLMVA_CC0pi->Add(hLMVA_LargestLMVA[fsi]);
	hLMVA_SmallestLMVA_CC0pi->Add(hLMVA_SmallestLMVA[fsi]);
      }
      hLMVA_1track_All->Add(hLMVA_1track[fsi]);
      hLMVA_2tracks_All->Add(hLMVA_2tracks[fsi]);
      hLMVA_LargestLMVA_All->Add(hLMVA_LargestLMVA[fsi]);
      hLMVA_SmallestLMVA_All->Add(hLMVA_SmallestLMVA[fsi]);
    } 
  }

  hLMVA_2tracks_CC0pi->Rebin2D(RebinValueLMVA,RebinValueLMVA);
  hLMVA_2tracks_CC1pi->Rebin2D(RebinValueLMVA,RebinValueLMVA);
  hLMVA_2tracks_All->Rebin2D(RebinValueLMVA,RebinValueLMVA);


  //Purity and efficiency histogram creation
  TH1D * hEfficiency_1track_CC0pi = (TH1D*) hLMVA_1track_CC0pi->Clone("hEfficiency_1track_CC0pi");hEfficiency_1track_CC0pi->Reset();
  TH1D * hPurity_1track_CC0pi = (TH1D*) hLMVA_1track_CC0pi->Clone("hPurity_1track_CC0pi");hPurity_1track_CC0pi->Reset();
  TH1D * hError_1track_CC0pi = (TH1D*) hLMVA_1track_CC0pi->Clone("hError_1track_CC0pi");hError_1track_CC0pi->Reset();
  EfficiencyPurity1D(hEfficiency_1track_CC0pi, hPurity_1track_CC0pi, hError_1track_CC0pi, hLMVA_1track_CC0pi, hLMVA_1track_All, ErrorBkg, true);
    
  TH2D * hEfficiency_2tracks_CC0pi = (TH2D*) hLMVA_2tracks_CC0pi->Clone("hEfficiency_2tracks_CC0pi");hEfficiency_2tracks_CC0pi->Reset();
  TH2D * hPurity_2tracks_CC0pi = (TH2D*) hLMVA_2tracks_CC0pi->Clone("hPurity_2tracks_CC0pi");hPurity_2tracks_CC0pi->Reset();
  TH2D * hError_2tracks_CC0pi = (TH2D*) hLMVA_2tracks_CC0pi->Clone("hError_2tracks_CC0pi");hError_2tracks_CC0pi->Reset();
  EfficiencyPurity2D(hEfficiency_2tracks_CC0pi, hPurity_2tracks_CC0pi, hError_2tracks_CC0pi, hLMVA_2tracks_CC0pi, hLMVA_2tracks_All, ErrorBkg, true, false);

  TH1D * hEfficiency_LargestLMVA_CC0pi = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("hEfficiency_LargestLMVA_CC0pi");hEfficiency_LargestLMVA_CC0pi->Reset();
  TH1D * hPurity_LargestLMVA_CC0pi = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("hPurity_LargestLMVA_CC0pi");hPurity_LargestLMVA_CC0pi->Reset();
  TH1D * hError_LargestLMVA_CC0pi = (TH1D*) hLMVA_LargestLMVA_CC0pi->Clone("hError_LargestLMVA_CC0pi");hError_LargestLMVA_CC0pi->Reset();
  EfficiencyPurity1D(hEfficiency_LargestLMVA_CC0pi, hPurity_LargestLMVA_CC0pi, hError_LargestLMVA_CC0pi, hLMVA_LargestLMVA_CC0pi, hLMVA_LargestLMVA_All, ErrorBkg, true);

  TH1D * hEfficiency_SmallestLMVA_CC0pi = (TH1D*) hLMVA_SmallestLMVA_CC0pi->Clone("hEfficiency_SmallestLMVA_CC0pi");hEfficiency_SmallestLMVA_CC0pi->Reset();
  TH1D * hPurity_SmallestLMVA_CC0pi = (TH1D*) hLMVA_SmallestLMVA_CC0pi->Clone("hPurity_SmallestLMVA_CC0pi");hPurity_SmallestLMVA_CC0pi->Reset();
  TH1D * hError_SmallestLMVA_CC0pi = (TH1D*) hLMVA_SmallestLMVA_CC0pi->Clone("hError_SmallestLMVA_CC0pi");hError_SmallestLMVA_CC0pi->Reset();
  EfficiencyPurity1D(hEfficiency_SmallestLMVA_CC0pi, hPurity_SmallestLMVA_CC0pi, hError_SmallestLMVA_CC0pi, hLMVA_SmallestLMVA_CC0pi, hLMVA_SmallestLMVA_All, ErrorBkg, false);

  TH2D * hEfficiency_2tracks_CC1pi = (TH2D*) hLMVA_2tracks_CC1pi->Clone("hEfficiency_2tracks_CC1pi");hEfficiency_2tracks_CC1pi->Reset();
  TH2D * hPurity_2tracks_CC1pi = (TH2D*) hLMVA_2tracks_CC1pi->Clone("hPurity_2tracks_CC1pi");hPurity_2tracks_CC1pi->Reset();
  TH2D * hError_2tracks_CC1pi = (TH2D*) hLMVA_2tracks_CC1pi->Clone("hError_2tracks_CC1pi");hError_2tracks_CC1pi->Reset();
  EfficiencyPurity2D(hEfficiency_2tracks_CC1pi, hPurity_2tracks_CC1pi, hError_2tracks_CC1pi, hLMVA_2tracks_CC1pi, hLMVA_2tracks_All, ErrorBkg, true, true);

  //Plots
  TLegend * l1D = new TLegend(0.11,0.11,0.35,0.35);
  l1D->SetFillColor(0);l1D->SetLineColor(0);
  l1D->AddEntry(hEfficiency_1track_CC0pi,"Efficiency");
  l1D->AddEntry(hPurity_1track_CC0pi,"Purity");
  l1D->AddEntry(hError_1track_CC0pi,"Total error");
  
  TCanvas * cLMVA_1track_CC0pi = new TCanvas("cLMVA_1track_CC0pi","Efficiency, purity and error of the CC0pi selected using the LMVA for 1 track sample");
  PrepareCanvas1D(cLMVA_1track_CC0pi, hEfficiency_1track_CC0pi, hPurity_1track_CC0pi, hError_1track_CC0pi, l1D);
  cLMVA_1track_CC0pi->SaveAs("../plots/MVAOptimisation/LMVA_1track_EfficiencyPurity.eps");

  TCanvas * cLMVA_LargestLMVA_CC0pi = new TCanvas("cLMVA_LargestLMVA_CC0pi","Efficiency, purity and error of the CC0pi selected using the LMVA for 1 track sample");
  PrepareCanvas1D(cLMVA_LargestLMVA_CC0pi, hEfficiency_LargestLMVA_CC0pi, hPurity_LargestLMVA_CC0pi, hError_LargestLMVA_CC0pi, l1D);
  cLMVA_LargestLMVA_CC0pi->SaveAs("../plots/MVAOptimisation/LMVA_LargestLMVA_EfficiencyPurity.eps");

  TCanvas * cLMVA_SmallestLMVA_CC0pi = new TCanvas("cLMVA_SmallestLMVA_CC0pi","Efficiency, purity and error of the CC0pi selected using the LMVA for 1 track sample");
  PrepareCanvas1D(cLMVA_SmallestLMVA_CC0pi, hEfficiency_SmallestLMVA_CC0pi, hPurity_SmallestLMVA_CC0pi, hError_SmallestLMVA_CC0pi, l1D);
  cLMVA_SmallestLMVA_CC0pi->SaveAs("../plots/MVAOptimisation/LMVA_SmallestLMVA_EfficiencyPurity.eps");

  TCanvas * cLMVAEfficiency_2tracks_CC0pi = new TCanvas("cLMVAEfficiency_2tracks_CC0pi","Efficiency of the CC0pi using the LMVA discriminant for 2 tracks sample");//Use the 2 tracks sample
  hEfficiency_2tracks_CC0pi->Draw("COLZTEXT");
  hEfficiency_2tracks_CC0pi->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  hEfficiency_2tracks_CC0pi->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  cLMVAEfficiency_2tracks_CC0pi->SaveAs("../plots/MVAOptimisation/LMVA_2tracks_CC0pi_Efficiency.eps");

  TCanvas * cLMVAPurity_2tracks_CC0pi = new TCanvas("cLMVAPurity_2tracks_CC0pi","Purity of the CC0pi using the LMVA discriminant for 2 tracks sample");//Use the 2 tracks sample
  hPurity_2tracks_CC0pi->Draw("COLZTEXT");
  hPurity_2tracks_CC0pi->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  hPurity_2tracks_CC0pi->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  cLMVAPurity_2tracks_CC0pi->SaveAs("../plots/MVAOptimisation/LMVA_2tracks_CC0pi_Purity.eps");

  TCanvas * cLMVAError_2tracks_CC0pi = new TCanvas("cLMVAError_2tracks_CC0pi","Error of the CC0pi using the LMVA discriminant for 2 tracks sample");//Use the 2 tracks sample
  hError_2tracks_CC0pi->Draw("COLZTEXT");
  hError_2tracks_CC0pi->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  hError_2tracks_CC0pi->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  cLMVAError_2tracks_CC0pi->SetLogz();
  cLMVAError_2tracks_CC0pi->SaveAs("../plots/MVAOptimisation/LMVA_2tracks_CC0pi_Error.eps");

  TCanvas * cLMVAEfficiency_2tracks_CC1pi = new TCanvas("cLMVAEfficiency_2tracks_CC1pi","Efficiency of the CC1pi using the LMVA discriminant for 2 tracks sample");//Use the 2 tracks sample
  hEfficiency_2tracks_CC1pi->Draw("COLZTEXT");
  hEfficiency_2tracks_CC1pi->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  hEfficiency_2tracks_CC1pi->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  cLMVAEfficiency_2tracks_CC1pi->SaveAs("../plots/MVAOptimisation/LMVA_2tracks_CC1pi_Efficiency.eps");

  TCanvas * cLMVAPurity_2tracks_CC1pi = new TCanvas("cLMVAPurity_2tracks_CC1pi","Purity of the CC1pi using the LMVA discriminant for 2 tracks sample");//Use the 2 tracks sample
  hPurity_2tracks_CC1pi->Draw("COLZTEXT");
  hPurity_2tracks_CC1pi->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  hPurity_2tracks_CC1pi->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  cLMVAPurity_2tracks_CC1pi->SaveAs("../plots/MVAOptimisation/LMVA_2tracks_CC1pi_Purity.eps");

  TCanvas * cLMVAError_2tracks_CC1pi = new TCanvas("cLMVAError_2tracks_CC1pi","Error of the CC1pi using the LMVA discriminant for 2 tracks sample");//Use the 2 tracks sample
  hError_2tracks_CC1pi->Draw("COLZTEXT");
  hError_2tracks_CC1pi->GetXaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  hError_2tracks_CC1pi->GetYaxis()->SetRangeUser(LMVAZoomInf,LMVAZoomSup);
  cLMVAError_2tracks_CC1pi->SetLogz();
  cLMVAError_2tracks_CC1pi->SaveAs("../plots/MVAOptimisation/LMVA_2tracks_CC1pi_Error.eps");
  ////////////////////////////////////////////


  TH1D * hLMVA_TrueParticle[NSimplifiedPDG][NSamples];

  TH1D * hLMVA_TrueParticle_Muon[NSamples];
  TH1D * hLMVA_TrueParticle_Proton[NSamples];
  TH1D * hLMVA_TrueParticle_All[NSamples];

  TH1D * hLMVA_TrueParticle_FullSample_Muon;
  TH1D * hLMVA_TrueParticle_FullSample_Proton;
  TH1D * hLMVA_TrueParticle_FullSample_All;
  
  for(int j=0;j<NSamples;j++){
    
    for(int i=0;i<NSimplifiedPDG;i++){
      hLMVA_TrueParticle[i][j] = (TH1D*) PIDCutTuning->Get(Form("hLMVA_TrueParticle%d_Sample%d",i,j));
      //hLMVA_TrueParticle[i][j]->SetFillStyle(0);

      if(i == 0){
	hLMVA_TrueParticle_Muon[j] = (TH1D*) hLMVA_TrueParticle[i][j]->Clone(Form("hLMVA_TrueParticle_Muon_Sample%d",j)); hLMVA_TrueParticle_Muon[j]->SetFillStyle(0);
	hLMVA_TrueParticle_All[j] = (TH1D*) hLMVA_TrueParticle[i][j]->Clone(Form("hLMVA_TrueParticle_All_Sample%d",j)); hLMVA_TrueParticle_All[j]->SetFillStyle(0);
      }	
      if(i == 2){ hLMVA_TrueParticle_Proton[j] = (TH1D*) hLMVA_TrueParticle[i][j]->Clone(Form("hLMVA_TrueParticle_Proton_Sample%d",j)); hLMVA_TrueParticle_Proton[j]->SetFillStyle(0);}
      
      hLMVA_TrueParticle_All[j]->Add(hLMVA_TrueParticle[i][j]);
    }


    if(j==0){
      hLMVA_TrueParticle_FullSample_Muon = (TH1D*) hLMVA_TrueParticle_Muon[j]->Clone(Form("hLMVA_TrueParticle_FullSample_Muon",j)); hLMVA_TrueParticle_FullSample_Muon[j]->SetFillStyle(0);
      hLMVA_TrueParticle_FullSample_Proton = (TH1D*) hLMVA_TrueParticle_Proton[j]->Clone(Form("hLMVA_TrueParticle_FullSample_Proton",j)); hLMVA_TrueParticle_FullSample_Proton[j]->SetFillStyle(0);
      hLMVA_TrueParticle_FullSample_All = (TH1D*) hLMVA_TrueParticle_All[j]->Clone(Form("hLMVA_TrueParticle_FullSample_All",j)); hLMVA_TrueParticle_FullSample_All[j]->SetFillStyle(0);
    }
    else{
      hLMVA_TrueParticle_FullSample_Muon->Add(hLMVA_TrueParticle_Muon[j]);
      hLMVA_TrueParticle_FullSample_Proton->Add(hLMVA_TrueParticle_Proton[j]);
      hLMVA_TrueParticle_FullSample_All->Add(hLMVA_TrueParticle_All[j]);
    }
  }    

  //To create: Stack plots of discriminant of muon, pion, proton, others for each sample
  //Efficiency/purity for an upper cut of muons for each sample and total
  //Efficiency/purity for a lower cut of protons for each sample and total

  TCanvas * cLMVA_TrueParticle[NSamples];
  THStack * sLMVA_TrueParticle[NSamples];
  TLegend * lPart = new TLegend(0.11,0.6,0.3,0.89);
  lPart->SetFillColor(0);
  lPart->SetLineColor(0);
  lPart->AddEntry(hLMVA_TrueParticle[0][0],"#mu");
  lPart->AddEntry(hLMVA_TrueParticle[1][0],"#pi");
  lPart->AddEntry(hLMVA_TrueParticle[2][0],"p");
  lPart->AddEntry(hLMVA_TrueParticle[3][0],"Others");
 
  for(int i=0;i<NSamples;i++){
    sLMVA_TrueParticle[i] = new THStack(Form("sLMVA_TrueParticle_Sample%d",i),"");
    ProduceStackParticles(hLMVA_TrueParticle[0][i],hLMVA_TrueParticle[1][i],hLMVA_TrueParticle[2][i],hLMVA_TrueParticle[3][i],sLMVA_TrueParticle[i]);
    cLMVA_TrueParticle[i] =  new TCanvas(Form("cLMVA_TrueParticle_Sample%d",i),"Stack of the CC0pi selected using the LMVA for 1 track sample");
    sLMVA_TrueParticle[i]->Draw();
    sLMVA_TrueParticle[i]->GetXaxis()->SetTitle(hLMVA_TrueParticle[0][i]->GetXaxis()->GetTitle());
    sLMVA_TrueParticle[i]->GetYaxis()->SetTitle(hLMVA_TrueParticle[0][i]->GetYaxis()->GetTitle());
    lPart->Draw("same");
  }
  
  TCanvas * cLMVA_TrueParticle_Muon[NSamples];
  TH1D * hEfficiency_TrueParticle_Muon[NSamples];
  TH1D * hPurity_TrueParticle_Muon[NSamples];
  TH1D * hError_TrueParticle_Muon[NSamples];
  
  for(int i=0;i<NSamples;i++){

    cout<<"Sample "<<i<<endl;
    cout<<"Integral signal: "<< hLMVA_TrueParticle_Muon[i] ->Integral() << ", all : " << hLMVA_TrueParticle_All[i]->Integral() << endl;

    hEfficiency_TrueParticle_Muon[i] = (TH1D*) hLMVA_TrueParticle_Muon[i]->Clone(Form("hEfficiency_TrueParticle_Muon_Sample%d",i)); hEfficiency_TrueParticle_Muon[i]->SetFillStyle(0);
    hPurity_TrueParticle_Muon[i] = (TH1D*) hLMVA_TrueParticle_Muon[i]->Clone(Form("hPurity_TrueParticle_Muon_Sample%d",i)); hPurity_TrueParticle_Muon[i]->SetFillStyle(0);
    hError_TrueParticle_Muon[i] = (TH1D*) hLMVA_TrueParticle_Muon[i]->Clone(Form("hError_TrueParticle_Muon_Sample%d",i)); hError_TrueParticle_Muon[i]->SetFillStyle(0);
    EfficiencyPurity1D(hEfficiency_TrueParticle_Muon[i], hPurity_TrueParticle_Muon[i], hError_TrueParticle_Muon[i], hLMVA_TrueParticle_Muon[i], hLMVA_TrueParticle_All[i], ErrorBkg, true);

    cLMVA_TrueParticle_Muon[i] = new TCanvas(Form("cLMVA_TrueParticle_Muon_Sample%d",i),"Efficiency, purity and error of the CC0pi selected using the LMVA for 1 track sample");
    PrepareCanvas1D(cLMVA_TrueParticle_Muon[i], hEfficiency_TrueParticle_Muon[i], hPurity_TrueParticle_Muon[i], hError_TrueParticle_Muon[i], l1D);
    cLMVA_TrueParticle_Muon[i]->SaveAs(Form("../plots/MVAOptimisation/LMVA_TrueParticle_Muon_Sample%d_EfficiencyPurity.eps",i));

  }

  

  ////////////////////////////////////////////
  

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

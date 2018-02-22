void ProduceStackParticles(TH1D * hmu, TH1D * hpi, TH1D * hp, TH1D * ho, THStack * hStack){

  hmu->GetYaxis()->SetTitleOffset(1.3);

  hmu->SetFillColor(kBlue);
  hmu->SetLineColor(1);
  hmu->SetLineWidth(2);

  hpi->SetFillColor(kGreen+2);
  hpi->SetLineColor(1);
  hpi->SetLineWidth(2);

  hp->SetFillColor(kRed);
  hp->SetLineColor(1);
  hp->SetLineWidth(2);

  ho->SetFillColor(kGray);
  ho->SetLineColor(1);
  ho->SetLineWidth(2);

  hStack->Add(hmu);
  hStack->Add(hpi);
  hStack->Add(hp);
  hStack->Add(ho);

  hmu->SetTitle("#mu");
  hpi->SetTitle("#pi");
  hp->SetTitle("p");
  ho->SetTitle("others");

}
void MVAOptimizer(){

  gROOT->SetBatch(kTRUE);//silent mode
  gROOT->SetStyle("Plain");
  gStyle->SetFillColor(10);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(kTRUE);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  //gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetPadLeftMargin(0.13);
  gStyle->SetOptTitle(kFALSE);

  //#################################################################################################
  //###############CHECK ALSO RIGHT TAIL CUT FOR PROTON. USE ANOTHER MVA, MVAProton##################
  bool MVAPion=false;
  bool MVAProton=true;
  //#################################################################################################

  
  const int NBinsMVA = 160; 
  double LowerBoundMVA = -0.5 ;
  double UpperBoundMVA = 0.5 ;
  const int NBinsCL = 100;
  double LowerBoundCL = 0;
  double UpperBoundCL = 1;

  //TFile * rfile = new TFile("../src/MVAparticle_MVAmuononly.root");
  //TFile * rfile = new TFile("../src/MVAparticle_MVAmuonpion.root");
  //TFile * rfile = new TFile("../src/MVAparticle.root");
  //TFile * rfile = new TFile("../src/MVAparticle.root");
  //TFile * rfile = new TFile("../src/MVAparticle_multiclassMVAmuonpion_10000trees.root");
  //TFile * rfile = new TFile("../src/MVAparticle_MVAmuonpion_bdtg.root");
 
  //TFile * rfile = new TFile("../src/MVAparticleMuon_1000trees_BugTrueParticleSolved.root");
  //TFile * rfile = new TFile("../src/MVAparticleMuon_1000trees.root");
  TFile * rfile = new TFile("../src/MVAWMparticleMuon_1000trees.root");

  //TFile * rfile = new TFile("../src/MVAparticleMuon_Forest1000trees_Input1000trees.root");
  //TFile * rfile = new TFile("../src/MVAparticleProton_1000trees.root");
  
  TTree * rtree = (TTree*) rfile->Get("TrainTree");
  //TTree * rtree = (TTree*) rfile->Get("TestTree");

  float  FSIInt;//Final state true information
  float  Num_Int;//Interaction number
  float  nTracks;//number of tracks reconstructed per vertex
  float weight;//weight of each event
  float IsFV;//True vertex is in FV
  float IsSand;
  float IsAnti;
  float IsBkgH;
  float IsBkgV;
  float IsNuE;
  float IsDetected;//Event is reconstructed (means >=1 vertex reconstructed)
  float SelectionFV;float SelectionOV;//Vertex is reconstructed within/out of FV
  float OpeningAngle;//Opening angle between reconstructed tracks - only relevant for 2 track samples
  float CoplanarityAngle;//Coplanarity angle between reconstructed tracks - only relevant for 2 track samples
  float Enu;//True neutrino energy
  float  nIngBasRec;//Number of vertices reconstructed per event
  float NewEvent;//true if new event, false if it is the same event than previous but only another reconstructed one.
  float Errorweight;
  float TrueMomentumMuon, TrueAngleMuon;//True muon kinematic variables
  float POT;//Number of POT

  float TrackAngle;//Reconstructed angle of the reconstructed track
  float  TypeOfTrack;//pdg of the reconstructed track. Careful that it is based on an algorithm defined in Reconstruction.cc
  float CLMuon_Likelihood;// = new vector<float> ;
  float TotalCharge;//Charge of the track per unit distance
  float Momentum;
  float EquivalentIronDistance;
  float  Sample;//Geometric properties of the track, defined in Reconstruction::SelectTrackSample
  float IsReconstructed;
  float BDT;
  
  TBranch * Br_TrackAngle;
  TBranch * Br_TypeOfTrack;
  TBranch * Br_IsReconstructed;
  TBranch * Br_Sample;
  TBranch * Br_CLMuon_Likelihood;
  TBranch * Br_EquivalentIronDistance;
  TBranch * Br_TotalCharge;
  TBranch * Br_Momentum;
  TBranch * Br_BDT;

  Br_TrackAngle = rtree->GetBranch(Form("TrackAngle"));
  Br_TrackAngle->SetAddress(&(TrackAngle));
  rtree->SetBranchAddress(Form("TrackAngle"),&(TrackAngle));
  
  Br_TypeOfTrack = rtree->GetBranch(Form("TypeOfTrack"));
  Br_TypeOfTrack->SetAddress(&(TypeOfTrack));
  rtree->SetBranchAddress(Form("TypeOfTrack"),&(TypeOfTrack));
  
  Br_IsReconstructed = rtree->GetBranch(Form("IsReconstructed"));
  Br_IsReconstructed->SetAddress(&(IsReconstructed));
  rtree->SetBranchAddress(Form("IsReconstructed"),&(IsReconstructed));
  
  Br_Sample = rtree->GetBranch(Form("Sample"));
  Br_Sample->SetAddress(&(Sample));
  rtree->SetBranchAddress(Form("Sample"),&(Sample));
  
  Br_CLMuon_Likelihood = rtree->GetBranch(Form("CLMuon_Likelihood"));
  Br_CLMuon_Likelihood->SetAddress(&(CLMuon_Likelihood));
  rtree->SetBranchAddress(Form("CLMuon_Likelihood"),&(CLMuon_Likelihood));
  
  Br_EquivalentIronDistance = rtree->GetBranch(Form("EquivalentIronDistance"));
  Br_EquivalentIronDistance->SetAddress(&(EquivalentIronDistance));
  rtree->SetBranchAddress(Form("EquivalentIronDistance"),&(EquivalentIronDistance));
  
  Br_TotalCharge = rtree->GetBranch(Form("TotalCharge"));
  Br_TotalCharge->SetAddress(&(TotalCharge));
  rtree->SetBranchAddress(Form("TotalCharge"),&(TotalCharge));
  
  Br_Momentum = rtree->GetBranch(Form("Momentum"));
  Br_Momentum->SetAddress(&(Momentum));
  rtree->SetBranchAddress(Form("Momentum"),&(Momentum));
  
  Br_BDT = rtree->GetBranch(Form("BDT"));
  Br_BDT->SetAddress(&(BDT));
  rtree->SetBranchAddress(Form("BDT"),&(BDT));

  Br_weight = rtree->GetBranch(Form("weight"));
  Br_weight->SetAddress(&(weight));
  rtree->SetBranchAddress(Form("weight"),&(weight));

  //Create the plot of the muon, pion and proton distribution.
  TH1D * MVAdiscriminant_TrueMuon[6];
  TH1D * MVAdiscriminant_TruePion[6];
  TH1D * MVAdiscriminant_TrueProton[6];
  TH1D * MVAdiscriminant_TrueOthers[6];
  THStack * StackMVAdiscriminant[6];

  for(int is=0;is<6;is++){
    MVAdiscriminant_TrueMuon[is] = new TH1D(Form("MVAdiscriminant_TrueMuon%d",is),"MVA discriminant of true muons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAdiscriminant_TruePion[is] = new TH1D(Form("MVAdiscriminant_TruePion%d",is),"MVA discriminant of true pions",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAdiscriminant_TrueProton[is] = new TH1D(Form("MVAdiscriminant_TrueProton%d",is),"MVA discriminant of true protons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAdiscriminant_TrueOthers[is] = new TH1D(Form("MVAdiscriminant_TrueOthers%d",is),"MVA discriminant of true others",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    StackMVAdiscriminant[is] = new THStack(Form("StackMVAdiscriminant%d",is),"");
  }
  TH1D * MVAdiscriminant_TrueMuon_AllSamples;
  TH1D * MVAdiscriminant_TruePion_AllSamples;
  TH1D * MVAdiscriminant_TrueProton_AllSamples;
  TH1D * MVAdiscriminant_TrueOthers_AllSamples;
  THStack * StackMVAdiscriminant_AllSamples = new THStack("StackMVAdiscriminant_AllSamples","");

  
  TH1D * EfficiencyMVA_TrueMuon = new TH1D("EfficiencyMVA_TrueMuon","Right tail cut efficiency of muons using MVA discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * PurityMVA_TrueMuon = new TH1D("PurityMVA_TrueMuon","Right tail cut efficiency of muons using MVA discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * EfficiencyMVA_TrueProton = new TH1D("EfficiencyMVA_TrueProton","Left tail cut efficiency of muons using MVA discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * PurityMVA_TrueProton = new TH1D("PurityMVA_TrueProton","Left tail cut efficiency of muons using MVA discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);

  TH1D * PurityMVA_TrueMuonVSPion = new TH1D("PurityMVA_TrueMuonVSPion","Right tail cut efficiency of muons vs pions using MVA discriminant: muon/pion separation power",NBinsMVA,LowerBoundMVA,UpperBoundMVA);

  TH1D * EfficiencyMVA_TrueMuon_INGRIDStopping = new TH1D("EfficiencyMVA_TrueMuon_INGRIDStopping","Right tail cut efficiency of muons using MVA discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * PurityMVA_TrueMuon_INGRIDStopping = new TH1D("PurityMVA_TrueMuon_INGRIDStopping","Right tail cut efficiency of muons using MVA discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);


  double PurityMuonMVA[NBinsMVA];
  double EfficiencyMuonMVA[NBinsMVA];
  double PurityProtonMVA[NBinsMVA];
  double EfficiencyProtonMVA[NBinsMVA];

  //
  TH1D * CLdiscriminant_TrueMuon = new TH1D("CLdiscriminant_TrueMuon","CL discriminant of true muons",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * CLdiscriminant_TruePion = new TH1D("CLdiscriminant_TruePion","CL discriminant of true pions",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * CLdiscriminant_TrueProton = new TH1D("CLdiscriminant_TrueProton","CL discriminant of true protons",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * CLdiscriminant_TrueOthers = new TH1D("CLdiscriminant_TrueOthers","CL discriminant of true others",NBinsCL,LowerBoundCL,UpperBoundCL);

  TH1D * EfficiencyCL_TrueMuon = new TH1D("EfficiencyCL_TrueMuon","Right tail cut efficiency of muons using CL discriminant",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * PurityCL_TrueMuon = new TH1D("PurityCL_TrueMuon","Right tail cut efficiency of muons using CL discriminant",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * EfficiencyCL_TrueProton = new TH1D("EfficiencyCL_TrueProton","Left tail cut efficiency of muons using CL discriminant",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * PurityCL_TrueProton = new TH1D("PurityCL_TrueProton","Left tail cut efficiency of muons using CL discriminant",NBinsCL,LowerBoundCL,UpperBoundCL);

  TH1D * PurityCL_TrueMuonVSPion = new TH1D("PurityCL_TrueMuonVSPion","Right tail cut efficiency of muons vs pions using CL discriminant: muon/pion separation power",NBinsCL,LowerBoundCL,UpperBoundCL);

  TH1D * EfficiencyCL_TrueMuon_INGRIDStopping = new TH1D("EfficiencyCL_TrueMuon_INGRIDStopping","Right tail cut efficiency of muons using CL discriminant",NBinsCL,LowerBoundCL,UpperBoundCL);
  TH1D * PurityCL_TrueMuon_INGRIDStopping = new TH1D("PurityCL_TrueMuon_INGRIDStopping","Right tail cut efficiency of muons using CL discriminant",NBinsCL,LowerBoundCL,UpperBoundCL);

  THStack * StackCLdiscriminant = new THStack("StackCLdiscriminant","");

  double PurityMuonCL[NBinsCL];
  double EfficiencyMuonCL[NBinsCL];
  double PurityProtonCL[NBinsCL];
  double EfficiencyProtonCL[NBinsCL];

  int nevt = rtree->GetEntries();
  bool preselection = (IsFV) && (SelectionFV) && (IsDetected) && (IsReconstructed);//Only MC in FV is used for training 

  
  for(int ievt=0;ievt<nevt;ievt++){
    if(ievt%10000==0) cout<<"event "<<ievt<<"/"<<nevt<<endl;
    rtree->GetEntry(ievt);

    if(1/*preselection*/){
      //cout<<TypeOfTrack<<", "<<BDT<<", "<<weight<<endl;
      int iSample=(int) Sample;
      if(TMath::Abs(TypeOfTrack)==13 /*&& (iSample>=3)*/){
	CLdiscriminant_TrueMuon->Fill(CLMuon_Likelihood,weight);
	MVAdiscriminant_TrueMuon[iSample]->Fill(BDT,weight); 
	//CLdiscriminant_TrueMuon[iSample]->Fill(CLMuon_Likelihood,weight);
      }
      else if(TMath::Abs(TypeOfTrack)==211){
	CLdiscriminant_TruePion->Fill(CLMuon_Likelihood,weight);
	MVAdiscriminant_TruePion[iSample]->Fill(BDT,weight); 
	//CLdiscriminant_TruePion[iSample]->Fill(CLMuon_Likelihood,weight);
      }
      else if(TMath::Abs(TypeOfTrack)==2212){
	CLdiscriminant_TrueProton->Fill(CLMuon_Likelihood,weight);
	MVAdiscriminant_TrueProton[iSample]->Fill(BDT,weight); 
	//CLdiscriminant_TrueProton[iSample]->Fill(CLMuon_Likelihood,weight);
      }
      else{
	CLdiscriminant_TrueOthers->Fill(CLMuon_Likelihood,weight);
	MVAdiscriminant_TrueOthers[iSample]->Fill(BDT,weight); 
	//CLdiscriminant_TrueOthers[iSample]->Fill(CLMuon_Likelihood,weight);
      }
    }
  }

  for(int is=0;is<6;is++){
    if(is==0){
      MVAdiscriminant_TrueMuon_AllSamples = (TH1D*) MVAdiscriminant_TrueMuon[is]->Clone("MVAdiscriminant_TrueMuon_AllSamples");
      MVAdiscriminant_TruePion_AllSamples = (TH1D*) MVAdiscriminant_TruePion[is]->Clone("MVAdiscriminant_TruePion_AllSamples");
      MVAdiscriminant_TrueProton_AllSamples = (TH1D*) MVAdiscriminant_TrueProton[is]->Clone("MVAdiscriminant_TrueProton_AllSamples");
      MVAdiscriminant_TrueOthers_AllSamples = (TH1D*) MVAdiscriminant_TrueOthers[is]->Clone("MVAdiscriminant_TrueOthers_AllSamples");
    }
    else{
      MVAdiscriminant_TrueMuon_AllSamples->Add(MVAdiscriminant_TrueMuon[is]);
      MVAdiscriminant_TruePion_AllSamples->Add(MVAdiscriminant_TruePion[is]);
      MVAdiscriminant_TrueProton_AllSamples->Add(MVAdiscriminant_TrueProton[is]);
      MVAdiscriminant_TrueOthers_AllSamples->Add(MVAdiscriminant_TrueOthers[is]);
    }
  }
  
  ////////////////////////// MVA based plots/////////////////////////////////////
  TCanvas * cMVA_AllSamples = new TCanvas("cMVA_AllSamples","MVA distributions of true particles");
  MVAdiscriminant_TrueMuon_AllSamples->SetLineColor(kBlue);
  MVAdiscriminant_TruePion_AllSamples->SetLineColor(kGreen+2);
  MVAdiscriminant_TrueProton_AllSamples->SetLineColor(kRed);
  MVAdiscriminant_TrueOthers_AllSamples->SetLineColor(kGray);
  MVAdiscriminant_TrueMuon_AllSamples->Draw();
  MVAdiscriminant_TrueMuon_AllSamples->GetXaxis()->SetTitle("#mu_{MVA}");
  MVAdiscriminant_TrueMuon_AllSamples->GetYaxis()->SetTitle("Number of events");
  MVAdiscriminant_TruePion_AllSamples->Draw("same");
  MVAdiscriminant_TrueProton_AllSamples->Draw("same");
  MVAdiscriminant_TrueOthers_AllSamples->Draw("same");
  
  
  TLegend * lParticles = new TLegend(0.1,0.7,0.3,0.9);
  lParticles->AddEntry(MVAdiscriminant_TrueMuon_AllSamples,"#mu");
  lParticles->AddEntry(MVAdiscriminant_TruePion_AllSamples,"#pi");
  lParticles->AddEntry(MVAdiscriminant_TrueProton_AllSamples,"p");
  lParticles->Draw("same");

  ProduceStackParticles(MVAdiscriminant_TrueMuon_AllSamples,MVAdiscriminant_TruePion_AllSamples,MVAdiscriminant_TrueProton_AllSamples,MVAdiscriminant_TrueOthers_AllSamples,StackMVAdiscriminant_AllSamples);
  TCanvas * cMVAstack_AllSamples = new TCanvas("cMVAstack_AllSamples","Stack of MVA distributions of true particles");
  StackMVAdiscriminant_AllSamples->Draw();
  lParticles->Draw("same");
  //

  TCanvas * cMVA[6];
  TCanvas * cMVAstack[6];
  for(int is=0;is<6;is++){
    cMVA[is] = new TCanvas(Form("cMVA%d",is),Form("MVA distributions of true particles for sample %d",is));
    MVAdiscriminant_TrueMuon[is]->SetLineColor(kBlue);
    MVAdiscriminant_TruePion[is]->SetLineColor(kGreen+2);
    MVAdiscriminant_TrueProton[is]->SetLineColor(kRed);
    MVAdiscriminant_TrueOthers[is]->SetLineColor(kGray);
    MVAdiscriminant_TrueMuon[is]->Draw();
    MVAdiscriminant_TrueMuon[is]->GetXaxis()->SetTitle("#mu_{MVA}");
    MVAdiscriminant_TrueMuon[is]->GetYaxis()->SetTitle("Number of events");
    MVAdiscriminant_TruePion[is]->Draw("same");
    MVAdiscriminant_TrueProton[is]->Draw("same");
    MVAdiscriminant_TrueOthers[is]->Draw("same");
    lParticles->Draw("same");

    ProduceStackParticles(MVAdiscriminant_TrueMuon[is],MVAdiscriminant_TruePion[is],MVAdiscriminant_TrueProton[is],MVAdiscriminant_TrueOthers[is],StackMVAdiscriminant[is]);
    cMVAstack[is] = new TCanvas(Form("cMVAstack%d",is),Form("Stack of MVA distributions of true particles for sample %d",is));
    StackMVAdiscriminant[is]->Draw();
    lParticles->Draw("same");
  //
  }
  //

  //Only INGRID stopping sample
  for(int ibinx=1;ibinx<=MVAdiscriminant_TrueMuon[3]->GetNbinsX();ibinx++){
    double RightTailMuon=MVAdiscriminant_TrueMuon[3]->Integral(ibinx,MVAdiscriminant_TrueMuon[3]->GetNbinsX());
    double RightTailPion=MVAdiscriminant_TruePion[3]->Integral(ibinx,MVAdiscriminant_TruePion[3]->GetNbinsX());
    double RightTailProton=MVAdiscriminant_TrueProton[3]->Integral(ibinx,MVAdiscriminant_TrueProton[3]->GetNbinsX());
    double RightTailOthers=MVAdiscriminant_TrueOthers[3]->Integral(ibinx,MVAdiscriminant_TrueOthers[3]->GetNbinsX());
    
    double TotalMuon=MVAdiscriminant_TrueMuon[3]->Integral();
    double TotalPion=MVAdiscriminant_TruePion[3]->Integral();
    double TotalProton=MVAdiscriminant_TrueProton[3]->Integral();
    double TotalOthers=MVAdiscriminant_TrueOthers[3]->Integral();

    double PurityMuon = RightTailMuon;
    double EfficiencyMuon = RightTailMuon;
    if(TotalMuon!=0) EfficiencyMuon /= TotalMuon;
    if(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers) PurityMuon/=(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers);
    EfficiencyMVA_TrueMuon_INGRIDStopping->SetBinContent(ibinx,EfficiencyMuon);
    PurityMVA_TrueMuon_INGRIDStopping->SetBinContent(ibinx,PurityMuon);
  }

    //All samples together
  for(int ibinx=1;ibinx<=MVAdiscriminant_TrueMuon_AllSamples->GetNbinsX();ibinx++){
    double RightTailMuon=MVAdiscriminant_TrueMuon_AllSamples->Integral(ibinx,MVAdiscriminant_TrueMuon_AllSamples->GetNbinsX());
    double RightTailPion=MVAdiscriminant_TruePion_AllSamples->Integral(ibinx,MVAdiscriminant_TruePion_AllSamples->GetNbinsX());
    double RightTailProton=MVAdiscriminant_TrueProton_AllSamples->Integral(ibinx,MVAdiscriminant_TrueProton_AllSamples->GetNbinsX());
    double RightTailOthers=MVAdiscriminant_TrueOthers_AllSamples->Integral(ibinx,MVAdiscriminant_TrueOthers_AllSamples->GetNbinsX());

    double TotalMuon=MVAdiscriminant_TrueMuon_AllSamples->Integral();
    double TotalPion=MVAdiscriminant_TruePion_AllSamples->Integral();
    double TotalProton=MVAdiscriminant_TrueProton_AllSamples->Integral();
    double TotalOthers=MVAdiscriminant_TrueOthers_AllSamples->Integral();

    double PurityMuon = RightTailMuon;
    double EfficiencyMuon = RightTailMuon;
    if(TotalMuon!=0) EfficiencyMuon /= TotalMuon;
    if(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers) PurityMuon/=(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers);
    EfficiencyMVA_TrueMuon->SetBinContent(ibinx,EfficiencyMuon);
    PurityMVA_TrueMuon->SetBinContent(ibinx,PurityMuon);
    
    double LeftTailMuon=MVAdiscriminant_TrueMuon_AllSamples->Integral(1,ibinx);
    double LeftTailPion=MVAdiscriminant_TruePion_AllSamples->Integral(1,ibinx);
    double LeftTailProton=MVAdiscriminant_TrueProton_AllSamples->Integral(1,ibinx);
    double LeftTailOthers=MVAdiscriminant_TrueOthers_AllSamples->Integral(1,ibinx);

    double PurityProton = LeftTailProton;
    double EfficiencyProton = LeftTailProton;
    if(TotalProton!=0) EfficiencyProton /= TotalProton;
    if(LeftTailMuon+LeftTailPion+LeftTailProton+LeftTailOthers) PurityProton/=(LeftTailMuon+LeftTailPion+LeftTailProton+LeftTailOthers);
    EfficiencyMVA_TrueProton->SetBinContent(ibinx,EfficiencyProton);
    PurityMVA_TrueProton->SetBinContent(ibinx,PurityProton);

    double PurityMuonVSPion = RightTailMuon;
    if(RightTailMuon+RightTailPion) PurityMuonVSPion/=(RightTailMuon+RightTailPion);
    PurityMVA_TrueMuonVSPion->SetBinContent(ibinx,PurityMuonVSPion);

    EfficiencyMuonMVA[ibinx-1]=EfficiencyMuon;
    PurityMuonMVA[ibinx-1]=PurityMuon;
    EfficiencyProtonMVA[ibinx-1]=EfficiencyProton;
    PurityProtonMVA[ibinx-1]=PurityProton;
    //cout<<"Efficiency="<<EfficiencyMuonMVA[ibinx-1]<<", purity="<<PurityMuonMVA[ibinx-1]<<endl;
  }

  TCanvas * cEfficiencyMVA_TrueMuon = new TCanvas("cEfficiencyMVA_TrueMuon","Efficiency and purity of muons using MVA");
  EfficiencyMVA_TrueMuon->SetLineColor(kBlue);
  PurityMVA_TrueMuon->SetLineColor(kRed);
  EfficiencyMVA_TrueMuon->Draw();
  EfficiencyMVA_TrueMuon->GetXaxis()->SetTitle("#mu_{MVA}");
  PurityMVA_TrueMuon->Draw("same");
  EfficiencyMVA_TrueMuon->GetYaxis()->SetRangeUser(0.,1.);
  PurityMVA_TrueMuonVSPion->SetLineColor(kGreen+2);
  PurityMVA_TrueMuonVSPion->Draw("same");
  TLegend * lLeftTail = new TLegend(0.7,0.1,0.9,0.3);
  lLeftTail->AddEntry(EfficiencyMVA_TrueMuon,"Eff. #mu");
  lLeftTail->AddEntry(PurityMVA_TrueMuon,"Pur. #mu");
  lLeftTail->AddEntry(PurityMVA_TrueMuonVSPion,"Pur. #mu vs #pi+#mu");
  lLeftTail->Draw("same");

  TCanvas * cEfficiencyMVA_TrueMuon_INGRIDStopping = new TCanvas("cEfficiencyMVA_TrueMuon_INGRIDStopping","Efficiency and purity of muons using MVA");
  EfficiencyMVA_TrueMuon_INGRIDStopping->SetLineColor(kBlue);
  PurityMVA_TrueMuon_INGRIDStopping->SetLineColor(kRed);
  EfficiencyMVA_TrueMuon_INGRIDStopping->Draw();
  EfficiencyMVA_TrueMuon_INGRIDStopping->GetXaxis()->SetTitle("#mu_{MVA}");
  PurityMVA_TrueMuon_INGRIDStopping->Draw("same");
  EfficiencyMVA_TrueMuon_INGRIDStopping->GetYaxis()->SetRangeUser(0.,1.);
  lLeftTail->Draw("same");
  
  TCanvas * cEfficiencyMVA_TrueProton = new TCanvas("cEfficiencyMVA_TrueProton","Efficiency and purity of protons using MVA");
  EfficiencyMVA_TrueProton->SetLineColor(kBlue);
  PurityMVA_TrueProton->SetLineColor(kRed);
  EfficiencyMVA_TrueProton->Draw();
  EfficiencyMVA_TrueProton->GetXaxis()->SetTitle("#mu_{MVA}");
  PurityMVA_TrueProton->Draw("same");
  EfficiencyMVA_TrueProton->GetYaxis()->SetRangeUser(0.,1.);
  TLegend * lRightTail = new TLegend(0.7,0.1,0.9,0.3);
  lRightTail->AddEntry(EfficiencyMVA_TrueProton,"Eff. p");
  lRightTail->AddEntry(PurityMVA_TrueProton,"Pur. p");
  lRightTail->Draw("same");
 
  ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////// CL based plots/////////////////////////////////////
  TCanvas * cCL = new TCanvas("cCL","CL distributions of true particles");
  CLdiscriminant_TrueMuon->SetLineColor(kBlue);
  CLdiscriminant_TruePion->SetLineColor(kGreen+2);
  CLdiscriminant_TrueProton->SetLineColor(kRed);
  CLdiscriminant_TrueOthers->SetLineColor(kGray);
  CLdiscriminant_TrueMuon->Draw();
  CLdiscriminant_TrueMuon->GetXaxis()->SetTitle("#mu_{CL}");
  CLdiscriminant_TruePion->Draw("same");
  CLdiscriminant_TrueProton->Draw("same");
  CLdiscriminant_TrueOthers->Draw("same");
  lParticles->Draw("same");
  //
  
  ProduceStackParticles(CLdiscriminant_TrueMuon,CLdiscriminant_TruePion,CLdiscriminant_TrueProton,CLdiscriminant_TrueOthers,StackCLdiscriminant);
  TCanvas * cCLstack = new TCanvas("cCLstack","Stack of CL distributions of true particles");
  StackCLdiscriminant->Draw();
  lParticles->Draw("same");

  /*
  //For stopping sample
  for(int ibinx=1;ibinx<=CLdiscriminant_TrueMuon->GetNbinsX();ibinx++){
    double RightTailMuon=CLdiscriminant_TrueMuon->Integral(ibinx,CLdiscriminant_TrueMuon->GetNbinsX());
    double RightTailPion=CLdiscriminant_TruePion->Integral(ibinx,CLdiscriminant_TruePion->GetNbinsX());
    double RightTailProton=CLdiscriminant_TrueProton->Integral(ibinx,CLdiscriminant_TrueProton->GetNbinsX());
    double RightTailOthers=CLdiscriminant_TrueOthers->Integral(ibinx,CLdiscriminant_TrueOthers->GetNbinsX());

    double TotalMuon=CLdiscriminant_TrueMuon->Integral();
    double TotalPion=CLdiscriminant_TruePion->Integral();
    double TotalProton=CLdiscriminant_TrueProton->Integral();
    double TotalOthers=CLdiscriminant_TrueOthers->Integral();

    double PurityMuon = RightTailMuon;
    double EfficiencyMuon = RightTailMuon;
    if(TotalMuon!=0) EfficiencyMuon /= TotalMuon;
    if(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers) PurityMuon/=(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers);
    EfficiencyCL_TrueMuon->SetBinContent(ibinx,EfficiencyMuon);
    PurityCL_TrueMuon->SetBinContent(ibinx,PurityMuon);
  }
  */
  //All samples together
  for(int ibinx=1;ibinx<=CLdiscriminant_TrueMuon->GetNbinsX();ibinx++){
    double RightTailMuon=CLdiscriminant_TrueMuon->Integral(ibinx,CLdiscriminant_TrueMuon->GetNbinsX());
    double RightTailPion=CLdiscriminant_TruePion->Integral(ibinx,CLdiscriminant_TruePion->GetNbinsX());
    double RightTailProton=CLdiscriminant_TrueProton->Integral(ibinx,CLdiscriminant_TrueProton->GetNbinsX());
    double RightTailOthers=CLdiscriminant_TrueOthers->Integral(ibinx,CLdiscriminant_TrueOthers->GetNbinsX());

    double TotalMuon=CLdiscriminant_TrueMuon->Integral();
    double TotalPion=CLdiscriminant_TruePion->Integral();
    double TotalProton=CLdiscriminant_TrueProton->Integral();
    double TotalOthers=CLdiscriminant_TrueOthers->Integral();

    double PurityMuon = RightTailMuon;
    double EfficiencyMuon = RightTailMuon;
    if(TotalMuon!=0) EfficiencyMuon /= TotalMuon;
    if(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers) PurityMuon/=(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers);
    EfficiencyCL_TrueMuon->SetBinContent(ibinx,EfficiencyMuon);
    PurityCL_TrueMuon->SetBinContent(ibinx,PurityMuon);

    double LeftTailMuon=CLdiscriminant_TrueMuon->Integral(1,ibinx);
    double LeftTailPion=CLdiscriminant_TruePion->Integral(1,ibinx);
    double LeftTailProton=CLdiscriminant_TrueProton->Integral(1,ibinx);
    double LeftTailOthers=CLdiscriminant_TrueOthers->Integral(1,ibinx);

    double PurityProton = LeftTailProton;
    double EfficiencyProton = LeftTailProton;
    if(TotalProton!=0) EfficiencyProton /= TotalProton;
    if(LeftTailMuon+LeftTailPion+LeftTailProton+LeftTailOthers) PurityProton/=(LeftTailMuon+LeftTailPion+LeftTailProton+LeftTailOthers);
    EfficiencyCL_TrueProton->SetBinContent(ibinx,EfficiencyProton);
    PurityCL_TrueProton->SetBinContent(ibinx,PurityProton);

    double PurityMuonVSPion = RightTailMuon;
    if(RightTailMuon+RightTailPion) PurityMuonVSPion/=(RightTailMuon+RightTailPion);
    PurityCL_TrueMuonVSPion->SetBinContent(ibinx,PurityMuonVSPion);

    EfficiencyMuonCL[ibinx-1]=EfficiencyMuon;
    PurityMuonCL[ibinx-1]=PurityMuon;
    EfficiencyProtonCL[ibinx-1]=EfficiencyProton;
    PurityProtonCL[ibinx-1]=PurityProton;
    //cout<<"Efficiency="<<EfficiencyMuonMVA[ibinx-1]<<", purity="<<PurityMuonMVA[ibinx-1]<<endl;
  }

  TCanvas * cEfficiencyCL_TrueMuon = new TCanvas("cEfficiencyCL_TrueMuon","Efficiency and purity of muons using CL");
  EfficiencyCL_TrueMuon->SetLineColor(kBlue);
  PurityCL_TrueMuon->SetLineColor(kRed);
  EfficiencyCL_TrueMuon->Draw();
  EfficiencyCL_TrueMuon->GetXaxis()->SetTitle("#mu_{CL}");
  PurityCL_TrueMuon->Draw("same");
  EfficiencyCL_TrueMuon->GetYaxis()->SetRangeUser(0.0,1);
  PurityCL_TrueMuonVSPion->SetLineColor(kGreen+2);
  PurityCL_TrueMuonVSPion->Draw("same");
  lLeftTail->Draw("same");
  
  TCanvas * cEfficiencyCL_TrueProton = new TCanvas("cEfficiencyCL_TrueProton","Efficiency and purity of protons using CL");
  EfficiencyCL_TrueProton->SetLineColor(kBlue);
  PurityCL_TrueProton->SetLineColor(kRed);
  EfficiencyCL_TrueProton->Draw();
  EfficiencyCL_TrueProton->GetXaxis()->SetTitle("#mu_{CL}");
  PurityCL_TrueProton->Draw("same");
  EfficiencyCL_TrueProton->GetYaxis()->SetRangeUser(0.0,1);
  lRightTail->Draw("same");
  

  if(MVAProton){
    //TFile * rfile_MVAProton = new TFile("../src/MVAparticleProton_1000trees.root");
    TFile * rfile_MVAProton = new TFile("../src/MVAWMparticleProton_1000trees.root");
    //TFile * rfile_MVAProton = new TFile("../src/MVAparticleProton_1000trees_BugTrueParticleSolved.root");
    //TFile * rfile_MVAProton = new TFile("../src/MVAparticleProton_500trees_BugTrueParticleSolved_90PerCentEventsUsedForTraining.root");

    //TFile * rfile_MVAProton = new TFile("../src/MVAparticleProton_Forest1000trees_Input500trees.root");

    TTree * rtree_MVAProton = (TTree*) rfile_MVAProton->Get("TrainTree");
    //TTree * rtree_MVAProton = (TTree*) rfile_MVAProton->Get("TestTree");

  float  FSIInt;//Final state true information
  float  Num_Int;//Interaction number
  float  nTracks;//number of tracks reconstructed per vertex
  float weight;//weight of each event
  float IsFV;//True vertex is in FV
  float IsSand;
  float IsAnti;
  float IsBkgH;
  float IsBkgV;
  float IsNuE;
  float IsDetected;//Event is reconstructed (means >=1 vertex reconstructed)
  float SelectionFV;float SelectionOV;//Vertex is reconstructed within/out of FV
  float OpeningAngle;//Opening angle between reconstructed tracks - only relevant for 2 track samples
  float CoplanarityAngle;//Coplanarity angle between reconstructed tracks - only relevant for 2 track samples
  float Enu;//True neutrino energy
  float  nIngBasRec;//Number of vertices reconstructed per event
  float NewEvent;//true if new event, false if it is the same event than previous but only another reconstructed one.
  float Errorweight;
  float TrueMomentumMuon, TrueAngleMuon;//True muon kinematic variables
  float POT;//Number of POT

  float TrackAngle;//Reconstructed angle of the reconstructed track
  float  TypeOfTrack;//pdg of the reconstructed track. Careful that it is based on an algorithm defined in Reconstruction.cc
  float CLMuon_Likelihood;// = new vector<float> ;
  float TotalCharge;//Charge of the track per unit distance
  float Momentum;
  float EquivalentIronDistance;
  float  Sample;//Geometric properties of the track, defined in Reconstruction::SelectTrackSample
  float IsReconstructed;
  float BDT;
  
  TBranch * Br_TrackAngle;
  TBranch * Br_TypeOfTrack;
  TBranch * Br_IsReconstructed;
  TBranch * Br_Sample;
  TBranch * Br_CLMuon_Likelihood;
  TBranch * Br_EquivalentIronDistance;
  TBranch * Br_TotalCharge;
  TBranch * Br_Momentum;
  TBranch * Br_BDT;

  Br_TrackAngle = rtree_MVAProton->GetBranch(Form("TrackAngle"));
  Br_TrackAngle->SetAddress(&(TrackAngle));
  rtree_MVAProton->SetBranchAddress(Form("TrackAngle"),&(TrackAngle));
  
  Br_TypeOfTrack = rtree_MVAProton->GetBranch(Form("TypeOfTrack"));
  Br_TypeOfTrack->SetAddress(&(TypeOfTrack));
  rtree_MVAProton->SetBranchAddress(Form("TypeOfTrack"),&(TypeOfTrack));
  
  Br_IsReconstructed = rtree_MVAProton->GetBranch(Form("IsReconstructed"));
  Br_IsReconstructed->SetAddress(&(IsReconstructed));
  rtree_MVAProton->SetBranchAddress(Form("IsReconstructed"),&(IsReconstructed));
  
  Br_Sample = rtree_MVAProton->GetBranch(Form("Sample"));
  Br_Sample->SetAddress(&(Sample));
  rtree_MVAProton->SetBranchAddress(Form("Sample"),&(Sample));
  
  Br_CLMuon_Likelihood = rtree_MVAProton->GetBranch(Form("CLMuon_Likelihood"));
  Br_CLMuon_Likelihood->SetAddress(&(CLMuon_Likelihood));
  rtree_MVAProton->SetBranchAddress(Form("CLMuon_Likelihood"),&(CLMuon_Likelihood));
  
  Br_EquivalentIronDistance = rtree_MVAProton->GetBranch(Form("EquivalentIronDistance"));
  Br_EquivalentIronDistance->SetAddress(&(EquivalentIronDistance));
  rtree_MVAProton->SetBranchAddress(Form("EquivalentIronDistance"),&(EquivalentIronDistance));
  
  Br_TotalCharge = rtree_MVAProton->GetBranch(Form("TotalCharge"));
  Br_TotalCharge->SetAddress(&(TotalCharge));
  rtree_MVAProton->SetBranchAddress(Form("TotalCharge"),&(TotalCharge));
  
  Br_Momentum = rtree_MVAProton->GetBranch(Form("Momentum"));
  Br_Momentum->SetAddress(&(Momentum));
  rtree_MVAProton->SetBranchAddress(Form("Momentum"),&(Momentum));
  
  Br_BDT = rtree_MVAProton->GetBranch(Form("BDT"));
  Br_BDT->SetAddress(&(BDT));
  rtree_MVAProton->SetBranchAddress(Form("BDT"),&(BDT));

  Br_weight = rtree_MVAProton->GetBranch(Form("weight"));
  Br_weight->SetAddress(&(weight));
  rtree_MVAProton->SetBranchAddress(Form("weight"),&(weight));

  //Create the plot of the muon, pion and proton distribution.
  TH1D * MVAProtondiscriminant_TrueMuon[6];
  TH1D * MVAProtondiscriminant_TruePion[6];
  TH1D * MVAProtondiscriminant_TrueProton[6];
  TH1D * MVAProtondiscriminant_TrueOthers[6];
  THStack * StackMVAProtondiscriminant[6];

  for(int is=0;is<6;is++){
    MVAProtondiscriminant_TrueMuon[is] = new TH1D(Form("MVAProtondiscriminant_TrueMuon%d",is),"MVAProton discriminant of true muons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAProtondiscriminant_TruePion[is] = new TH1D(Form("MVAProtondiscriminant_TruePion%d",is),"MVAProton discriminant of true pions",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAProtondiscriminant_TrueProton[is] = new TH1D(Form("MVAProtondiscriminant_TrueProton%d",is),"MVAProton discriminant of true protons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAProtondiscriminant_TrueOthers[is] = new TH1D(Form("MVAProtondiscriminant_TrueOthers%d",is),"MVAProton discriminant of true others",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    StackMVAProtondiscriminant[is] = new THStack(Form("StackMVAProtondiscriminant%d",is),"");
  }
  TH1D * MVAProtondiscriminant_TrueMuon_AllSamples;
  TH1D * MVAProtondiscriminant_TruePion_AllSamples;
  TH1D * MVAProtondiscriminant_TrueProton_AllSamples;
  TH1D * MVAProtondiscriminant_TrueOthers_AllSamples;
  THStack * StackMVAProtondiscriminant_AllSamples = new THStack("StackMVAProtondiscriminant_AllSamples","");
  /*
  TH1D * MVAProtondiscriminant_TrueMuon = new TH1D("MVAProtondiscriminant_TrueMuon","MVAProton discriminant of true muons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * MVAProtondiscriminant_TruePion = new TH1D("MVAProtondiscriminant_TruePion","MVAProton discriminant of true pions",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * MVAProtondiscriminant_TrueProton = new TH1D("MVAProtondiscriminant_TrueProton","MVAProton discriminant of true protons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * MVAProtondiscriminant_TrueOthers = new TH1D("MVAProtondiscriminant_TrueOthers","MVAProton discriminant of true others",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  */
  TH1D * EfficiencyMVAProton_TrueProton = new TH1D("EfficiencyMVAProton_TrueProton","Left tail cut efficiency of muons using MVAProton discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * PurityMVAProton_TrueProton = new TH1D("PurityMVAProton_TrueProton","Left tail cut efficiency of muons using MVAProton discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);

  TH1D * PurityMVAProton_TrueProtonVSPion = new TH1D("PurityMVAProton_TrueProtonVSPion","Right tail cut efficiency of muons vs pions using MVAProton discriminant: muon/pion separation power",NBinsMVA,LowerBoundMVA,UpperBoundMVA);

  //THStack * StackMVAProtondiscriminant = new THStack("StackMVAProtondiscriminant","");

  double PurityProtonMVAProton[NBinsMVA];
  double EfficiencyProtonMVAProton[NBinsMVA];

  nevt = rtree_MVAProton->GetEntries();

  
  for(int ievt=0;ievt<nevt;ievt++){
    if(ievt%10000==0) cout<<"event "<<ievt<<"/"<<nevt<<endl;
    rtree_MVAProton->GetEntry(ievt);

    if(1/*preselection*/){
      //cout<<TypeOfTrack<<", "<<BDT<<", "<<weight<<endl;
      int iSample=(int) Sample;
      if(TMath::Abs(TypeOfTrack)==13 /*&& (Sample>=3)*/){
	MVAProtondiscriminant_TrueMuon[iSample]->Fill(BDT,weight); 
      }
      else if(TMath::Abs(TypeOfTrack)==211){
	MVAProtondiscriminant_TruePion[iSample]->Fill(BDT,weight); 
      }
      else if(TMath::Abs(TypeOfTrack)==2212){
	MVAProtondiscriminant_TrueProton[iSample]->Fill(BDT,weight);
      }
      else{
	MVAProtondiscriminant_TrueOthers[iSample]->Fill(BDT,weight);
      }
    }
  }
  for(int is=0;is<6;is++){
    if(is==0){
      MVAProtondiscriminant_TrueMuon_AllSamples = (TH1D*) MVAProtondiscriminant_TrueMuon[is]->Clone("MVAProtondiscriminant_TrueMuon_AllSamples");
      MVAProtondiscriminant_TruePion_AllSamples = (TH1D*) MVAProtondiscriminant_TruePion[is]->Clone("MVAProtondiscriminant_TruePion_AllSamples");
      MVAProtondiscriminant_TrueProton_AllSamples = (TH1D*) MVAProtondiscriminant_TrueProton[is]->Clone("MVAProtondiscriminant_TrueProton_AllSamples");
      MVAProtondiscriminant_TrueOthers_AllSamples = (TH1D*) MVAProtondiscriminant_TrueOthers[is]->Clone("MVAProtondiscriminant_TrueOthers_AllSamples");
    }
    else{
      MVAProtondiscriminant_TrueMuon_AllSamples->Add(MVAProtondiscriminant_TrueMuon[is]);
      MVAProtondiscriminant_TruePion_AllSamples->Add(MVAProtondiscriminant_TruePion[is]);
      MVAProtondiscriminant_TrueProton_AllSamples->Add(MVAProtondiscriminant_TrueProton[is]);
      MVAProtondiscriminant_TrueOthers_AllSamples->Add(MVAProtondiscriminant_TrueOthers[is]);
    }
  }
  
  ////////////////////////// MVA based plots/////////////////////////////////////

  TCanvas * cMVAProton_AllSamples = new TCanvas("cMVAProton_AllSamples","MVAProton distributions of true particles");
  MVAProtondiscriminant_TrueMuon_AllSamples->SetLineColor(kBlue);
  MVAProtondiscriminant_TruePion_AllSamples->SetLineColor(kGreen+2);
  MVAProtondiscriminant_TrueProton_AllSamples->SetLineColor(kRed);
  MVAProtondiscriminant_TrueOthers_AllSamples->SetLineColor(kGray);
  MVAProtondiscriminant_TrueMuon_AllSamples->Draw();
  MVAProtondiscriminant_TrueMuon_AllSamples->GetXaxis()->SetTitle("#mu_{MVAProton}");
  MVAProtondiscriminant_TrueMuon_AllSamples->GetYaxis()->SetTitle("Number of events");
  MVAProtondiscriminant_TruePion_AllSamples->Draw("same");
  MVAProtondiscriminant_TrueProton_AllSamples->Draw("same");
  MVAProtondiscriminant_TrueOthers_AllSamples->Draw("same");

  ProduceStackParticles(MVAProtondiscriminant_TrueMuon_AllSamples,MVAProtondiscriminant_TruePion_AllSamples,MVAProtondiscriminant_TrueProton_AllSamples,MVAProtondiscriminant_TrueOthers_AllSamples,StackMVAProtondiscriminant_AllSamples);
  TCanvas * cMVAProtonstack_AllSamples = new TCanvas("cMVAProtonstack_AllSamples","Stack of MVAProton distributions of true particles");
  StackMVAProtondiscriminant_AllSamples->Draw();
  lParticles->Draw("same");
  //

  TCanvas * cMVAProton[6];
  TCanvas * cMVAProtonstack[6];
  for(int is=0;is<6;is++){
    cMVAProton[is] = new TCanvas(Form("cMVAProton%d",is),Form("MVAProton distributions of true particles for sample %d",is));
    MVAProtondiscriminant_TrueMuon[is]->SetLineColor(kBlue);
    MVAProtondiscriminant_TruePion[is]->SetLineColor(kGreen+2);
    MVAProtondiscriminant_TrueProton[is]->SetLineColor(kRed);
    MVAProtondiscriminant_TrueOthers[is]->SetLineColor(kGray);
    MVAProtondiscriminant_TrueMuon[is]->Draw();
    MVAProtondiscriminant_TrueMuon[is]->GetXaxis()->SetTitle("p_{MVA}");
    MVAProtondiscriminant_TrueMuon[is]->GetYaxis()->SetTitle("Number of events");
    MVAProtondiscriminant_TruePion[is]->Draw("same");
    MVAProtondiscriminant_TrueProton[is]->Draw("same");
    MVAProtondiscriminant_TrueOthers[is]->Draw("same");
    lParticles->Draw("same");

    ProduceStackParticles(MVAProtondiscriminant_TrueMuon[is],MVAProtondiscriminant_TruePion[is],MVAProtondiscriminant_TrueProton[is],MVAProtondiscriminant_TrueOthers[is],StackMVAProtondiscriminant[is]);
    cMVAProtonstack[is] = new TCanvas(Form("cMVAProtonstack%d",is),Form("Stack of MVAProton distributions of true particles for sample %d",is));
    StackMVAProtondiscriminant[is]->Draw();
    lParticles->Draw("same");
  //
  }
  //
  /*
  ////////////////////////// MVAProton based plots/////////////////////////////////////
  TCanvas * cMVAProton = new TCanvas("cMVAProton","MVAProton distributions of true particles");
  MVAProtondiscriminant_TrueMuon->SetLineColor(kBlue);
  MVAProtondiscriminant_TruePion->SetLineColor(kGreen+2);
  MVAProtondiscriminant_TrueProton->SetLineColor(kRed);
  MVAProtondiscriminant_TrueOthers->SetLineColor(kGray);
  MVAProtondiscriminant_TrueMuon->Draw();
  MVAProtondiscriminant_TrueMuon->GetXaxis()->SetTitle("p_{MVA}");
  MVAProtondiscriminant_TrueMuon->GetYaxis()->SetTitle("Number of events");
  MVAProtondiscriminant_TruePion->Draw("same");
  MVAProtondiscriminant_TrueProton->Draw("same");
  MVAProtondiscriminant_TrueOthers->Draw("same");
  lParticles->Draw("same");
  
  //
  
  ProduceStackParticles(MVAProtondiscriminant_TrueMuon,MVAProtondiscriminant_TruePion,MVAProtondiscriminant_TrueProton,MVAProtondiscriminant_TrueOthers,StackMVAProtondiscriminant);
  TCanvas * cMVAProtonstack = new TCanvas("cMVAProtonstack","Stack of MVAProton distributions of true particles");
  StackMVAProtondiscriminant->Draw();
  lParticles->Draw("same");
  */
  //

  for(int ibinx=1;ibinx<=MVAProtondiscriminant_TrueMuon_AllSamples->GetNbinsX();ibinx++){
    double RightTailMuon=MVAProtondiscriminant_TrueMuon_AllSamples->Integral(ibinx,MVAProtondiscriminant_TrueMuon_AllSamples->GetNbinsX());
    double RightTailPion=MVAProtondiscriminant_TruePion_AllSamples->Integral(ibinx,MVAProtondiscriminant_TruePion_AllSamples->GetNbinsX());
    double RightTailProton=MVAProtondiscriminant_TrueProton_AllSamples->Integral(ibinx,MVAProtondiscriminant_TrueProton_AllSamples->GetNbinsX());
    double RightTailOthers=MVAProtondiscriminant_TrueOthers_AllSamples->Integral(ibinx,MVAProtondiscriminant_TrueOthers_AllSamples->GetNbinsX());

    double TotalMuon=MVAProtondiscriminant_TrueMuon_AllSamples->Integral();
    double TotalPion=MVAProtondiscriminant_TruePion_AllSamples->Integral();
    double TotalProton=MVAProtondiscriminant_TrueProton_AllSamples->Integral();
    double TotalOthers=MVAProtondiscriminant_TrueOthers_AllSamples->Integral();

    double PurityProton = RightTailProton;
    double EfficiencyProton = RightTailProton;
    if(TotalProton!=0) EfficiencyProton /= TotalProton;
    if(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers) PurityProton/=(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers);
    EfficiencyMVAProton_TrueProton->SetBinContent(ibinx,EfficiencyProton);
    PurityMVAProton_TrueProton->SetBinContent(ibinx,PurityProton);

    double PurityProtonVSPion = RightTailProton;
    if(RightTailProton+RightTailPion) PurityProtonVSPion/=(RightTailProton+RightTailPion);
    PurityMVAProton_TrueProtonVSPion->SetBinContent(ibinx,PurityProtonVSPion);
    
    EfficiencyProtonMVAProton[ibinx-1]=EfficiencyProton;
    PurityProtonMVAProton[ibinx-1]=PurityProton;
    //cout<<"Efficiency="<<EfficiencyMuonMVAProton[ibinx-1]<<", purity="<<PurityMuonMVAProton[ibinx-1]<<endl;
  }
  
  TCanvas * cEfficiencyMVAProton_TrueProton = new TCanvas("cEfficiencyMVAProton_TrueProton","Efficiency and purity of muons using MVAProton");
  EfficiencyMVAProton_TrueProton->SetLineColor(kBlue);
  PurityMVAProton_TrueProton->SetLineColor(kRed);
  EfficiencyMVAProton_TrueProton->Draw();
  EfficiencyMVAProton_TrueProton->GetXaxis()->SetTitle("p_{MVA}");
  PurityMVAProton_TrueProton->Draw("same");
  EfficiencyMVAProton_TrueProton->GetYaxis()->SetRangeUser(0.0,1);
  PurityMVAProton_TrueProtonVSPion->SetLineColor(kGreen+2);
  PurityMVAProton_TrueProtonVSPion->Draw("same");

  TLegend * lLeftTailProton = new TLegend(0.7,0.1,0.9,0.3);
  lLeftTailProton->AddEntry(EfficiencyMVAProton_TrueProton,"Eff. p");
  lLeftTailProton->AddEntry(PurityMVAProton_TrueProton,"Pur. p");
  lLeftTailProton->AddEntry(PurityMVAProton_TrueProtonVSPion,"Pur. p vs #pi+p");
  lLeftTailProton->Draw("same");
  
  ////////////////////////////////////////////////////////////////////////////////
  }//End of MVA proton


  if(MVAPion){
  TFile * rfile_MVAPion = new TFile("../src/MVAparticleProton_1000trees.root");
  
  TTree * rtree_MVAPion = (TTree*) rfile_MVAPion->Get("TrainTree");

  float  FSIInt;//Final state true information
  float  Num_Int;//Interaction number
  float  nTracks;//number of tracks reconstructed per vertex
  float weight;//weight of each event
  float IsFV;//True vertex is in FV
  float IsSand;
  float IsAnti;
  float IsBkgH;
  float IsBkgV;
  float IsNuE;
  float IsDetected;//Event is reconstructed (means >=1 vertex reconstructed)
  float SelectionFV;float SelectionOV;//Vertex is reconstructed within/out of FV
  float OpeningAngle;//Opening angle between reconstructed tracks - only relevant for 2 track samples
  float CoplanarityAngle;//Coplanarity angle between reconstructed tracks - only relevant for 2 track samples
  float Enu;//True neutrino energy
  float  nIngBasRec;//Number of vertices reconstructed per event
  float NewEvent;//true if new event, false if it is the same event than previous but only another reconstructed one.
  float Errorweight;
  float TrueMomentumMuon, TrueAngleMuon;//True muon kinematic variables
  float POT;//Number of POT

  float TrackAngle;//Reconstructed angle of the reconstructed track
  float  TypeOfTrack;//pdg of the reconstructed track. Careful that it is based on an algorithm defined in Reconstruction.cc
  float CLMuon_Likelihood;// = new vector<float> ;
  float TotalCharge;//Charge of the track per unit distance
  float Momentum;
  float EquivalentIronDistance;
  float  Sample;//Geometric properties of the track, defined in Reconstruction::SelectTrackSample
  float IsReconstructed;
  float BDT;
  
  TBranch * Br_TrackAngle;
  TBranch * Br_TypeOfTrack;
  TBranch * Br_IsReconstructed;
  TBranch * Br_Sample;
  TBranch * Br_CLMuon_Likelihood;
  TBranch * Br_EquivalentIronDistance;
  TBranch * Br_TotalCharge;
  TBranch * Br_Momentum;
  TBranch * Br_BDT;

  Br_TrackAngle = rtree_MVAPion->GetBranch(Form("TrackAngle"));
  Br_TrackAngle->SetAddress(&(TrackAngle));
  rtree_MVAPion->SetBranchAddress(Form("TrackAngle"),&(TrackAngle));
  
  Br_TypeOfTrack = rtree_MVAPion->GetBranch(Form("TypeOfTrack"));
  Br_TypeOfTrack->SetAddress(&(TypeOfTrack));
  rtree_MVAPion->SetBranchAddress(Form("TypeOfTrack"),&(TypeOfTrack));
  
  Br_IsReconstructed = rtree_MVAPion->GetBranch(Form("IsReconstructed"));
  Br_IsReconstructed->SetAddress(&(IsReconstructed));
  rtree_MVAPion->SetBranchAddress(Form("IsReconstructed"),&(IsReconstructed));
  
  Br_Sample = rtree_MVAPion->GetBranch(Form("Sample"));
  Br_Sample->SetAddress(&(Sample));
  rtree_MVAPion->SetBranchAddress(Form("Sample"),&(Sample));
  
  Br_CLMuon_Likelihood = rtree_MVAPion->GetBranch(Form("CLMuon_Likelihood"));
  Br_CLMuon_Likelihood->SetAddress(&(CLMuon_Likelihood));
  rtree_MVAPion->SetBranchAddress(Form("CLMuon_Likelihood"),&(CLMuon_Likelihood));
  
  Br_EquivalentIronDistance = rtree_MVAPion->GetBranch(Form("EquivalentIronDistance"));
  Br_EquivalentIronDistance->SetAddress(&(EquivalentIronDistance));
  rtree_MVAPion->SetBranchAddress(Form("EquivalentIronDistance"),&(EquivalentIronDistance));
  
  Br_TotalCharge = rtree_MVAPion->GetBranch(Form("TotalCharge"));
  Br_TotalCharge->SetAddress(&(TotalCharge));
  rtree_MVAPion->SetBranchAddress(Form("TotalCharge"),&(TotalCharge));
  
  Br_Momentum = rtree_MVAPion->GetBranch(Form("Momentum"));
  Br_Momentum->SetAddress(&(Momentum));
  rtree_MVAPion->SetBranchAddress(Form("Momentum"),&(Momentum));
  
  Br_BDT = rtree_MVAPion->GetBranch(Form("BDT"));
  Br_BDT->SetAddress(&(BDT));
  rtree_MVAPion->SetBranchAddress(Form("BDT"),&(BDT));

  Br_weight = rtree_MVAPion->GetBranch(Form("weight"));
  Br_weight->SetAddress(&(weight));
  rtree_MVAPion->SetBranchAddress(Form("weight"),&(weight));

  //Create the plot of the muon, pion and proton distribution.
  TH1D * MVAPiondiscriminant_TrueMuon[6];
  TH1D * MVAPiondiscriminant_TruePion[6];
  TH1D * MVAPiondiscriminant_TrueProton[6];
  TH1D * MVAPiondiscriminant_TrueOthers[6];
  THStack * StackMVAPiondiscriminant[6];

  for(int is=0;is<6;is++){
    MVAPiondiscriminant_TrueMuon[is] = new TH1D(Form("MVAPiondiscriminant_TrueMuon%d",is),"MVAPion discriminant of true muons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAPiondiscriminant_TruePion[is] = new TH1D(Form("MVAPiondiscriminant_TruePion%d",is),"MVAPion discriminant of true pions",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAPiondiscriminant_TrueProton[is] = new TH1D(Form("MVAPiondiscriminant_TrueProton%d",is),"MVAPion discriminant of true protons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    MVAPiondiscriminant_TrueOthers[is] = new TH1D(Form("MVAPiondiscriminant_TrueOthers%d",is),"MVAPion discriminant of true others",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
    StackMVAPiondiscriminant[is] = new THStack(Form("StackMVAPiondiscriminant%d",is),"");
  }
  TH1D * MVAPiondiscriminant_TrueMuon_AllSamples;
  TH1D * MVAPiondiscriminant_TruePion_AllSamples;
  TH1D * MVAPiondiscriminant_TrueProton_AllSamples;
  TH1D * MVAPiondiscriminant_TrueOthers_AllSamples;
  THStack * StackMVAPiondiscriminant_AllSamples = new THStack("StackMVAPiondiscriminant_AllSamples","");
  /*
  TH1D * MVAPiondiscriminant_TrueMuon = new TH1D("MVAPiondiscriminant_TrueMuon","MVAPion discriminant of true muons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * MVAPiondiscriminant_TruePion = new TH1D("MVAPiondiscriminant_TruePion","MVAPion discriminant of true pions",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * MVAPiondiscriminant_TrueProton = new TH1D("MVAPiondiscriminant_TrueProton","MVAPion discriminant of true protons",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * MVAPiondiscriminant_TrueOthers = new TH1D("MVAPiondiscriminant_TrueOthers","MVAPion discriminant of true others",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  */
  TH1D * EfficiencyMVAPion_TrueProton = new TH1D("EfficiencyMVAPion_TrueProton","Left tail cut efficiency of muons using MVAPion discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);
  TH1D * PurityMVAPion_TrueProton = new TH1D("PurityMVAPion_TrueProton","Left tail cut efficiency of muons using MVAPion discriminant",NBinsMVA,LowerBoundMVA,UpperBoundMVA);

  TH1D * PurityMVAPion_TrueProtonVSPion = new TH1D("PurityMVAPion_TrueProtonVSPion","Right tail cut efficiency of muons vs pions using MVAPion discriminant: muon/pion separation power",NBinsMVA,LowerBoundMVA,UpperBoundMVA);

  //THStack * StackMVAPiondiscriminant = new THStack("StackMVAPiondiscriminant","");

  double PurityProtonMVAPion[NBinsMVA];
  double EfficiencyProtonMVAPion[NBinsMVA];

  nevt = rtree_MVAPion->GetEntries();

  
  for(int ievt=0;ievt<nevt;ievt++){
    if(ievt%10000==0) cout<<"event "<<ievt<<"/"<<nevt<<endl;
    rtree_MVAPion->GetEntry(ievt);

    if(1/*preselection*/){
      //cout<<TypeOfTrack<<", "<<BDT<<", "<<weight<<endl;
      int iSample=(int) Sample;
      if(TMath::Abs(TypeOfTrack)==13 /*&& (Sample>=3)*/){
	MVAPiondiscriminant_TrueMuon[iSample]->Fill(BDT,weight); 
      }
      else if(TMath::Abs(TypeOfTrack)==211){
	MVAPiondiscriminant_TruePion[iSample]->Fill(BDT,weight); 
      }
      else if(TMath::Abs(TypeOfTrack)==2212){
	MVAPiondiscriminant_TrueProton[iSample]->Fill(BDT,weight);
      }
      else{
	MVAPiondiscriminant_TrueOthers[iSample]->Fill(BDT,weight);
      }
    }
  }
  for(int is=0;is<6;is++){
    if(is==0){
      MVAPiondiscriminant_TrueMuon_AllSamples = (TH1D*) MVAPiondiscriminant_TrueMuon[is]->Clone("MVAPiondiscriminant_TrueMuon_AllSamples");
      MVAPiondiscriminant_TruePion_AllSamples = (TH1D*) MVAPiondiscriminant_TruePion[is]->Clone("MVAPiondiscriminant_TruePion_AllSamples");
      MVAPiondiscriminant_TrueProton_AllSamples = (TH1D*) MVAPiondiscriminant_TrueProton[is]->Clone("MVAPiondiscriminant_TrueProton_AllSamples");
      MVAPiondiscriminant_TrueOthers_AllSamples = (TH1D*) MVAPiondiscriminant_TrueOthers[is]->Clone("MVAPiondiscriminant_TrueOthers_AllSamples");
    }
    else{
      MVAPiondiscriminant_TrueMuon_AllSamples->Add(MVAPiondiscriminant_TrueMuon[is]);
      MVAPiondiscriminant_TruePion_AllSamples->Add(MVAPiondiscriminant_TruePion[is]);
      MVAPiondiscriminant_TrueProton_AllSamples->Add(MVAPiondiscriminant_TrueProton[is]);
      MVAPiondiscriminant_TrueOthers_AllSamples->Add(MVAPiondiscriminant_TrueOthers[is]);
    }
  }
  
  ////////////////////////// MVA based plots/////////////////////////////////////

  TCanvas * cMVAPion_AllSamples = new TCanvas("cMVAPion_AllSamples","MVAPion distributions of true particles");
  MVAPiondiscriminant_TrueMuon_AllSamples->SetLineColor(kBlue);
  MVAPiondiscriminant_TruePion_AllSamples->SetLineColor(kGreen+2);
  MVAPiondiscriminant_TrueProton_AllSamples->SetLineColor(kRed);
  MVAPiondiscriminant_TrueOthers_AllSamples->SetLineColor(kGray);
  MVAPiondiscriminant_TrueMuon_AllSamples->Draw();
  MVAPiondiscriminant_TrueMuon_AllSamples->GetXaxis()->SetTitle("#mu_{MVAPion}");
  MVAPiondiscriminant_TrueMuon_AllSamples->GetYaxis()->SetTitle("Number of events");
  MVAPiondiscriminant_TruePion_AllSamples->Draw("same");
  MVAPiondiscriminant_TrueProton_AllSamples->Draw("same");
  MVAPiondiscriminant_TrueOthers_AllSamples->Draw("same");

  ProduceStackParticles(MVAPiondiscriminant_TrueMuon_AllSamples,MVAPiondiscriminant_TruePion_AllSamples,MVAPiondiscriminant_TrueProton_AllSamples,MVAPiondiscriminant_TrueOthers_AllSamples,StackMVAPiondiscriminant_AllSamples);
  TCanvas * cMVAPionstack_AllSamples = new TCanvas("cMVAPionstack_AllSamples","Stack of MVAPion distributions of true particles");
  StackMVAPiondiscriminant_AllSamples->Draw();
  lParticles->Draw("same");
  //

  TCanvas * cMVAPion[6];
  TCanvas * cMVAPionstack[6];
  for(int is=0;is<6;is++){
    cMVAPion[is] = new TCanvas(Form("cMVAPion%d",is),Form("MVAPion distributions of true particles for sample %d",is));
    MVAPiondiscriminant_TrueMuon[is]->SetLineColor(kBlue);
    MVAPiondiscriminant_TruePion[is]->SetLineColor(kGreen+2);
    MVAPiondiscriminant_TrueProton[is]->SetLineColor(kRed);
    MVAPiondiscriminant_TrueOthers[is]->SetLineColor(kGray);
    MVAPiondiscriminant_TrueMuon[is]->Draw();
    MVAPiondiscriminant_TrueMuon[is]->GetXaxis()->SetTitle("p_{MVA}");
    MVAPiondiscriminant_TrueMuon[is]->GetYaxis()->SetTitle("Number of events");
    MVAPiondiscriminant_TruePion[is]->Draw("same");
    MVAPiondiscriminant_TrueProton[is]->Draw("same");
    MVAPiondiscriminant_TrueOthers[is]->Draw("same");
    lParticles->Draw("same");

    ProduceStackParticles(MVAPiondiscriminant_TrueMuon[is],MVAPiondiscriminant_TruePion[is],MVAPiondiscriminant_TrueProton[is],MVAPiondiscriminant_TrueOthers[is],StackMVAPiondiscriminant[is]);
    cMVAPionstack[is] = new TCanvas(Form("cMVAPionstack%d",is),Form("Stack of MVAPion distributions of true particles for sample %d",is));
    StackMVAPiondiscriminant[is]->Draw();
    lParticles->Draw("same");
  //
  }
  //
  /*
  ////////////////////////// MVAPion based plots/////////////////////////////////////
  TCanvas * cMVAPion = new TCanvas("cMVAPion","MVAPion distributions of true particles");
  MVAPiondiscriminant_TrueMuon->SetLineColor(kBlue);
  MVAPiondiscriminant_TruePion->SetLineColor(kGreen+2);
  MVAPiondiscriminant_TrueProton->SetLineColor(kRed);
  MVAPiondiscriminant_TrueOthers->SetLineColor(kGray);
  MVAPiondiscriminant_TrueMuon->Draw();
  MVAPiondiscriminant_TrueMuon->GetXaxis()->SetTitle("p_{MVA}");
  MVAPiondiscriminant_TrueMuon->GetYaxis()->SetTitle("Number of events");
  MVAPiondiscriminant_TruePion->Draw("same");
  MVAPiondiscriminant_TrueProton->Draw("same");
  MVAPiondiscriminant_TrueOthers->Draw("same");
  lParticles->Draw("same");
  
  //
  
  ProduceStackParticles(MVAPiondiscriminant_TrueMuon,MVAPiondiscriminant_TruePion,MVAPiondiscriminant_TrueProton,MVAPiondiscriminant_TrueOthers,StackMVAPiondiscriminant);
  TCanvas * cMVAPionstack = new TCanvas("cMVAPionstack","Stack of MVAPion distributions of true particles");
  StackMVAPiondiscriminant->Draw();
  lParticles->Draw("same");
  */
  //

  for(int ibinx=1;ibinx<=MVAPiondiscriminant_TrueMuon_AllSamples->GetNbinsX();ibinx++){
    double RightTailMuon=MVAPiondiscriminant_TrueMuon_AllSamples->Integral(ibinx,MVAPiondiscriminant_TrueMuon_AllSamples->GetNbinsX());
    double RightTailPion=MVAPiondiscriminant_TruePion_AllSamples->Integral(ibinx,MVAPiondiscriminant_TruePion_AllSamples->GetNbinsX());
    double RightTailProton=MVAPiondiscriminant_TrueProton_AllSamples->Integral(ibinx,MVAPiondiscriminant_TrueProton_AllSamples->GetNbinsX());
    double RightTailOthers=MVAPiondiscriminant_TrueOthers_AllSamples->Integral(ibinx,MVAPiondiscriminant_TrueOthers_AllSamples->GetNbinsX());

    double TotalMuon=MVAPiondiscriminant_TrueMuon_AllSamples->Integral();
    double TotalPion=MVAPiondiscriminant_TruePion_AllSamples->Integral();
    double TotalProton=MVAPiondiscriminant_TrueProton_AllSamples->Integral();
    double TotalOthers=MVAPiondiscriminant_TrueOthers_AllSamples->Integral();

    double PurityProton = RightTailProton;
    double EfficiencyProton = RightTailProton;
    if(TotalProton!=0) EfficiencyProton /= TotalProton;
    if(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers) PurityProton/=(RightTailMuon+RightTailPion+RightTailProton+RightTailOthers);
    EfficiencyMVAPion_TrueProton->SetBinContent(ibinx,EfficiencyProton);
    PurityMVAPion_TrueProton->SetBinContent(ibinx,PurityProton);

    double PurityProtonVSPion = RightTailProton;
    if(RightTailProton+RightTailPion) PurityProtonVSPion/=(RightTailProton+RightTailPion);
    PurityMVAPion_TrueProtonVSPion->SetBinContent(ibinx,PurityProtonVSPion);
    
    EfficiencyProtonMVAPion[ibinx-1]=EfficiencyProton;
    PurityProtonMVAPion[ibinx-1]=PurityProton;
    //cout<<"Efficiency="<<EfficiencyMuonMVAPion[ibinx-1]<<", purity="<<PurityMuonMVAPion[ibinx-1]<<endl;
  }

  TCanvas * cEfficiencyMVAPion_TrueProton = new TCanvas("cEfficiencyMVAPion_TrueProton","Efficiency and purity of muons using MVAPion");
  EfficiencyMVAPion_TrueProton->SetLineColor(kBlue);
  PurityMVAPion_TrueProton->SetLineColor(kRed);
  EfficiencyMVAPion_TrueProton->Draw();
  EfficiencyMVAPion_TrueProton->GetXaxis()->SetTitle("p_{MVA}");
  PurityMVAPion_TrueProton->Draw("same");
  EfficiencyMVAPion_TrueProton->GetYaxis()->SetRangeUser(0.0,1);
  PurityMVAPion_TrueProtonVSPion->SetLineColor(kGreen+2);
  PurityMVAPion_TrueProtonVSPion->Draw("same");

  TLegend * lLeftTailProton = new TLegend(0.7,0.1,0.9,0.3);
  lLeftTailProton->AddEntry(EfficiencyMVAPion_TrueProton,"Eff. p");
  lLeftTailProton->AddEntry(PurityMVAPion_TrueProton,"Pur. p");
  lLeftTailProton->AddEntry(PurityMVAPion_TrueProtonVSPion,"Pur. p vs #pi+p");
  lLeftTailProton->Draw("same");
  
  ////////////////////////////////////////////////////////////////////////////////
  }//End of MVA pion





  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  TGraph * EfficiencyVSPurityMVA_TrueMuon = new TGraph(NBinsMVA,EfficiencyMuonMVA,PurityMuonMVA);
  TGraph * EfficiencyVSPurityCL_TrueMuon = new TGraph(NBinsCL,EfficiencyMuonCL,PurityMuonCL);
  TGraph * EfficiencyVSPurityMVA_TrueProton = new TGraph(NBinsMVA,EfficiencyProtonMVA,PurityProtonMVA);
  TGraph * EfficiencyVSPurityCL_TrueProton = new TGraph(NBinsCL,EfficiencyProtonCL,PurityProtonCL);
  TGraph * EfficiencyVSPurityMVAProton_TrueProton;

  if(MVAProton){
    cout<<"You selected MVAProton. Therefore, the proton cut is done keeping the right tail of the distribution of MVA Proton"<<endl;
    EfficiencyVSPurityMVAProton_TrueProton = new TGraph(NBinsMVA,EfficiencyProtonMVAProton,PurityProtonMVAProton);
  }
    
  
  TCanvas * cEfficiencyVSPurity_MVA = new TCanvas("cEfficiencyVSPurityMVA_TrueMuon","Efficiency VS purity of muons and protons using MVA");
  EfficiencyVSPurityMVA_TrueMuon->SetLineColor(kBlue);
  EfficiencyVSPurityMVA_TrueMuon->SetLineWidth(2);
  EfficiencyVSPurityMVA_TrueMuon->SetFillColor(0);
  EfficiencyVSPurityMVA_TrueMuon->Draw("ALC");
  EfficiencyVSPurityMVA_TrueMuon->GetXaxis()->SetTitle("Efficiency");
  EfficiencyVSPurityMVA_TrueMuon->GetYaxis()->SetTitle("Purity");

  EfficiencyVSPurityMVA_TrueProton->SetLineColor(kOrange);
  EfficiencyVSPurityMVA_TrueProton->SetLineWidth(2);
  EfficiencyVSPurityMVA_TrueProton->SetFillColor(0);
  EfficiencyVSPurityMVA_TrueProton->Draw("Csame");

  if(MVAProton){
    EfficiencyVSPurityMVAProton_TrueProton->SetLineColor(kRed);
    EfficiencyVSPurityMVAProton_TrueProton->SetLineWidth(2);
    EfficiencyVSPurityMVAProton_TrueProton->SetFillColor(0);
    EfficiencyVSPurityMVAProton_TrueProton->Draw("Csame");
  }
  EfficiencyVSPurityMVA_TrueMuon->GetXaxis()->SetTitle("Efficiency");
  EfficiencyVSPurityMVA_TrueMuon->GetYaxis()->SetTitle("Purity");
  EfficiencyVSPurityMVA_TrueMuon->GetXaxis()->SetRangeUser(0.,1.);
  EfficiencyVSPurityMVA_TrueMuon->GetYaxis()->SetRangeUser(0.,1.);

  TLegend * lEffficiencyPurity = new TLegend(0.3,0.1,0.55,0.3);
  lEffficiencyPurity->AddEntry(EfficiencyVSPurityMVA_TrueMuon,"#mu, #mu_{MVA}");
  lEffficiencyPurity->AddEntry(EfficiencyVSPurityMVA_TrueProton,"p, #mu_{MVA}");
  if(MVAProton) lEffficiencyPurity->AddEntry(EfficiencyVSPurityMVAProton_TrueProton,"p, p_{MVA}");
  lEffficiencyPurity->Draw("same");

  TCanvas * cEfficiencyVSPurity_CL = new TCanvas("cEfficiencyVSPurityCL_TrueMuon","Efficiency VS purity of muons and protons using CL");
  EfficiencyVSPurityCL_TrueMuon->SetLineColor(kBlue);
  EfficiencyVSPurityCL_TrueProton->SetLineColor(kRed);
  EfficiencyVSPurityCL_TrueMuon->SetLineWidth(2);
  EfficiencyVSPurityCL_TrueProton->SetLineWidth(2);
  EfficiencyVSPurityCL_TrueMuon->SetFillColor(0);
  EfficiencyVSPurityCL_TrueProton->SetFillColor(0);
  EfficiencyVSPurityCL_TrueMuon->Draw("ALC");
  EfficiencyVSPurityCL_TrueMuon->GetXaxis()->SetTitle("Efficiency");
  EfficiencyVSPurityCL_TrueMuon->GetYaxis()->SetTitle("Purity");

  EfficiencyVSPurityCL_TrueProton->Draw("Csame");
  EfficiencyVSPurityCL_TrueMuon->GetXaxis()->SetTitle("Efficiency");
  EfficiencyVSPurityCL_TrueMuon->GetYaxis()->SetTitle("Purity");
  EfficiencyVSPurityCL_TrueMuon->GetXaxis()->SetRangeUser(0.,1.);
  EfficiencyVSPurityCL_TrueMuon->GetYaxis()->SetRangeUser(0.,1.);
  lEffficiencyPurity->Draw("same");

  TCanvas * cEfficiencyVSPurity_MVACL = new TCanvas("cEfficiencyVSPurityMVACL_TrueMuon","Efficiency VS purity of muons and protons using both MVA and CL");
  EfficiencyVSPurityMVA_TrueMuon->Draw("ALC");
  EfficiencyVSPurityMVA_TrueProton->Draw("Csame");
  if(MVAProton) EfficiencyVSPurityMVAProton_TrueProton->Draw("Csame");
  EfficiencyVSPurityMVA_TrueMuon->GetXaxis()->SetTitle("Efficiency");
  EfficiencyVSPurityMVA_TrueMuon->GetYaxis()->SetTitle("Purity");
  EfficiencyVSPurityMVA_TrueMuon->GetXaxis()->SetRangeUser(0.4,1.);
  EfficiencyVSPurityMVA_TrueMuon->GetYaxis()->SetRangeUser(0.,1.);
  TGraph * EfficiencyVSPurityCL_TrueMuon_clone = (TGraph*) EfficiencyVSPurityCL_TrueMuon->Clone("EfficiencyVSPurityCL_TrueMuon_clone");
  TGraph * EfficiencyVSPurityCL_TrueProton_clone = (TGraph*) EfficiencyVSPurityCL_TrueProton->Clone("EfficiencyVSPurityCL_TrueProton_clone");
  EfficiencyVSPurityCL_TrueMuon_clone->SetLineColor(kCyan);
  EfficiencyVSPurityCL_TrueProton_clone->SetLineColor(kMagenta);
  EfficiencyVSPurityCL_TrueMuon_clone->Draw("Csame");
  EfficiencyVSPurityCL_TrueProton_clone->Draw("Csame");
  
  TLegend * lEffficiencyPurity_MVACL = new TLegend(0.3,0.1,0.55,0.45);
  lEffficiencyPurity_MVACL->AddEntry(EfficiencyVSPurityCL_TrueMuon_clone,"#mu, #mu_{CL}");
  lEffficiencyPurity_MVACL->AddEntry(EfficiencyVSPurityCL_TrueProton_clone,"p, #mu_{CL}");
  lEffficiencyPurity_MVACL->AddEntry(EfficiencyVSPurityMVA_TrueMuon,"#mu, #mu_{MVA}");
  lEffficiencyPurity_MVACL->AddEntry(EfficiencyVSPurityMVA_TrueProton,"p, #mu_{MVA}");
  if(MVAProton) lEffficiencyPurity_MVACL->AddEntry(EfficiencyVSPurityMVAProton_TrueProton,"p, p_{MVA}");
  lEffficiencyPurity_MVACL->Draw("same");

  /*
  TFile * output = new TFile("../plots/MVAOptimisation.root","RECREATE");
  cEfficiencyCL_TrueMuon->Write();
  cEfficiencyCL_TrueProton->Write();
  cEfficiencyVSPurity_CL->Write();
  */  
  /*  
  TCanvas * cEfficiencyVSPurity_TrueMuon = new TCanvas("cEfficiencyVSPurity_TrueMuon","Efficiency VS purity of muons");
  EfficiencyVSPurityMVA_TrueMuon->SetLineColor(kBlue);
  EfficiencyVSPurityMVA_TrueMuon->Draw("ACP");
  EfficiencyVSPurityCL_TrueMuon->SetLineColor(kRed);
  EfficiencyVSPurityCL_TrueMuon->Draw("CPsame");
  EfficiencyVSPurityMVA_TrueMuon->GetXaxis()->SetRangeUser(0,1);
  */
  
  /////////////////////////////////////////////////////////////////////////////////////////
  for(is=0;is<6;is++){
    cMVAstack[is]->SaveAs(Form("../plots/MVAOptimisation/MVADistributionParticles_stack_sample%d.png",is));
    if(MVAPion) cMVAPionstack[is]->SaveAs(Form("../plots/MVAOptimisation/MVAPionDistributionParticles_stack_sample%d.png",is));
    if(MVAProton) cMVAProtonstack[is]->SaveAs(Form("../plots/MVAOptimisation/MVAProtonDistributionParticles_stack_sample%d.png",is));
  }
  cMVAstack_AllSamples->SaveAs("../plots/MVAOptimisation/MVADistributionParticles_stack.png");
  if(MVAPion) cMVAPionstack_AllSamples->SaveAs("../plots/MVAOptimisation/MVAPionDistributionParticles_stack.png");
  if(MVAProton) cMVAProtonstack_AllSamples->SaveAs("../plots/MVAOptimisation/MVAProtonDistributionParticles_stack.png");
  cEfficiencyMVA_TrueMuon->SaveAs("../plots/MVAOptimisation/EfficiencyMVA_TrueMuon.png");
  cEfficiencyMVA_TrueMuon_INGRIDStopping->SaveAs("../plots/MVAOptimisation/EfficiencyMVA_TrueMuon_INGRIDStopping.png");
  //cEfficiencyMVA_TruePion->SaveAs("../plots/MVAOptimisation/EfficiencyMVA_TruePion.png");
  cEfficiencyMVA_TrueProton->SaveAs("../plots/MVAOptimisation/EfficiencyMVA_TrueProton.png");
  cEfficiencyVSPurity_MVA->SaveAs("../plots/MVAOptimisation/EfficiencyVSPurity_MVA.png");
  cCLstack->SaveAs("../plots/MVAOptimisation/CLDistributionParticles_stack.png");
  cEfficiencyCL_TrueMuon->SaveAs("../plots/MVAOptimisation/EfficiencyCL_TrueMuon.png");
  //cEfficiencyCL_TrueMuon_INGRIDStopping->SaveAs("../plots/MVAOptimisation/EfficiencyCL_TrueMuon_INGRIDStopping.png");
  //cEfficiencyCL_TruePion->SaveAs("../plots/MVAOptimisation/EfficiencyCL_TruePion.png");
  cEfficiencyCL_TrueProton->SaveAs("../plots/MVAOptimisation/EfficiencyCL_TrueProton.png");
  cEfficiencyVSPurity_CL->SaveAs("../plots/MVAOptimisation/EfficiencyVSPurity_CL.png");
  
  if(MVAProton) cEfficiencyMVAProton_TrueProton->SaveAs("../plots/MVAOptimisation/EfficiencyMVAProton_TrueProton.png");
  cEfficiencyVSPurity_MVACL->SaveAs("../plots/MVAOptimisation/EfficiencyVSPurity_MVACL.png");
}


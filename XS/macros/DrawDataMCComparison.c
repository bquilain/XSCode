{
  TFile * fMC = new TFile("../files/MCSelected_PM_Systematics0_0.root","READ");
  THStack * Stack_NTracks = (THStack*) fMC->Get("Stack_NTracks");
  THStack * Stack_RecMom = (THStack*) fMC->Get("Stack_RecMom");
  THStack * Stack_RecAngle = (THStack*) fMC->Get("Stack_RecAngle");
  THStack * Stack_MuCL = (THStack*) fMC->Get("Stack_MuCL");
  THStack * Stack_MVAMuondiscriminant_1track = (THStack*) fMC->Get("Stack_MVAMuondiscriminant_1track");
  THStack * Stack_MVAProtondiscriminant_2tracks = (THStack*) fMC->Get("Stack_MVAProtondiscriminant_2tracks");

  //double ScalingData = 5.8e20 / 1.2956e19;
  double ScalingData = 1;//5.8e20 / 1.2956e19;

  TFile * fData = new TFile("../files/DataSelected_PM_Systematics0_0.root","READ");
  TH1D * Data_NTracks = (TH1D*) fData->Get("Data_NTracks"); Data_NTracks->Scale(ScalingData);
  TH1D * Data_RecMom = (TH1D*) fData->Get("Data_RecMom"); Data_RecMom->Scale(ScalingData);
  TH1D * Data_RecAngle = (TH1D*) fData->Get("Data_RecAngle"); Data_RecAngle->Scale(ScalingData);
  TH1D * Data_MuCL = (TH1D*) fData->Get("Data_MuCL"); Data_MuCL->Scale(ScalingData);
  TH1D * Data_MVAMuondiscriminant_1track = (TH1D*) fData->Get("Data_MVAMuondiscriminant_1track"); Data_MVAMuondiscriminant_1track->Scale(ScalingData);
  TH1D * Data_MVAProtondiscriminant_2tracks = (TH1D*) fData->Get("Data_MVAProtondiscriminant_2tracks"); Data_MVAProtondiscriminant_2tracks->Scale(ScalingData);

  TCanvas * cNTracks = new TCanvas(Form("c%s",Data_NTracks->GetName()));
  Stack_NTracks->Draw();
  Data_NTracks->Draw("E1same");
  //cNTracks->SaveAs(Form("%s.png",Data_NTracks->GetName()));

  TCanvas * cRecMom = new TCanvas(Form("c%s",Data_RecMom->GetName()));
  Stack_RecMom->Draw();
  Data_RecMom->Draw("E1same");
  //cRecMom->SaveAs(Form("%s.png",Data_RecMom->GetName()));

  TCanvas * cRecAngle = new TCanvas(Form("c%s",Data_RecAngle->GetName()));
  Stack_RecAngle->Draw();
  Data_RecAngle->Draw("E1same");
  //cRecAngle->SaveAs(Form("%s.png",Data_RecAngle->GetName()));
  
  TCanvas * cMuCL = new TCanvas(Form("c%s",Data_MuCL->GetName()));
  Stack_MuCL->Draw();
  Data_MuCL->Draw("E1same");
  //cMuCL->SaveAs(Form("%s.png",Data_MuCL->GetName()));
  
  TCanvas * cMVAMuondiscriminant_1track = new TCanvas(Form("c%s",Data_MVAMuondiscriminant_1track->GetName()));
  Stack_MVAMuondiscriminant_1track->Draw();
  Data_MVAMuondiscriminant_1track->Draw("E1same");
  //cMVAMuondiscriminant_1track->SaveAs(Form("%s.png",Data_MVAMuondiscriminant_1track->GetName()));
  
  TCanvas * cMVAProtondiscriminant_2tracks = new TCanvas(Form("c%s",Data_MVAProtondiscriminant_2tracks->GetName()));
  Stack_MVAProtondiscriminant_2tracks->Draw();
  Data_MVAProtondiscriminant_2tracks->Draw("E1same");
  //cMVAProtondiscriminant_2tracks->SaveAs(Form("%s.png",Data_MVAProtondiscriminant_2tracks->GetName()));
  
  
}

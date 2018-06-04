{

  TFile * fInput;
  TH1D * Output_PMIng_Muon;
  TH1D * Output_PMSci_Muon;
  TH1D * Output_Ing_Muon;
  TH1D * Output_PMIng;
  TH1D * Output_PMSci;
  TH1D * Output_Ing;

  for(int i=0;i<100;i++){
    fInput = new TFile(Form("../root_input/XSFormat_PM_Run1_%d_Systematics0_0_Plan.root",i));
    if(i == 0){
      Output_PMIng_Muon = (TH1D*) ((TH1D*) fInput->Get("hTest_PMIng_Muon"))->Clone("Output_PMIng_Muon");
      Output_PMSci_Muon = (TH1D*) ((TH1D*) fInput->Get("hTest_PMSci"))->Clone("Output_PMSci");
      Output_Ing_Muon = (TH1D*) ((TH1D*) fInput->Get("hTest_Ing_Muon"))->Clone("Output_Ing_Muon");
      Output_PMIng = (TH1D*) ((TH1D*) fInput->Get("hTest_PMIng"))->Clone("Output_PMIng");
      Output_PMSci = (TH1D*) ((TH1D*) fInput->Get("hTest_PMSci"))->Clone("Output_PMSci");
      Output_Ing = (TH1D*) ((TH1D*) fInput->Get("hTest_Ing"))->Clone("Output_Ing");
    }
    else{
      Output_PMIng_Muon->Add((TH1D*) fInput->Get("hTest_PMIng_Muon"));
      Output_PMSci_Muon->Add((TH1D*) fInput->Get("hTest_PMSci_Muon"));
      Output_Ing_Muon->Add((TH1D*) fInput->Get("hTest_Ing_Muon"));
      Output_PMIng->Add((TH1D*) fInput->Get("hTest_PMIng"));
      Output_PMSci->Add((TH1D*) fInput->Get("hTest_PMSci"));
      Output_Ing->Add((TH1D*) fInput->Get("hTest_Ing"));
    }
  }
  //TFile * fPDF = new TFile("../src/PDFMuCL_Likelihood.root");
  TFile * fPDF = new TFile("../src/TestPreFinalPDFMuCL_Likelihood.root");
  TH1D * PDF_PMIng_Muon = (TH1D*) fPDF->Get("PDFMuCL_PMIng_Muon");
  TH1D * PDF_PMSci_Muon = (TH1D*) fPDF->Get("PDFMuCL_PMSci_Muon");
  TH1D * PDF_Ing_Muon = (TH1D*) fPDF->Get("PDFMuCL_Ing_Muon");

  TCanvas * cPMIng_Muon = new TCanvas();
  Output_PMIng_Muon->DrawNormalized();
  PDF_PMIng_Muon->SetLineColor(2);
  PDF_PMIng_Muon->Draw("same");
  TCanvas * cPMSci_Muon = new TCanvas();
  Output_PMSci_Muon->DrawNormalized();
  PDF_PMSci_Muon->SetLineColor(2);
  PDF_PMSci_Muon->Draw("same");
  TCanvas * cIng_Muon = new TCanvas();
  Output_Ing_Muon->DrawNormalized();
  PDF_Ing_Muon->SetLineColor(2);
  PDF_Ing_Muon->Draw("same");

  TH1D * PDF_PMIng = (TH1D*) fPDF->Get("PDFMuCL_PMIng");
  TH1D * PDF_PMSci = (TH1D*) fPDF->Get("PDFMuCL_PMSci");
  TH1D * PDF_Ing = (TH1D*) fPDF->Get("PDFMuCL_Ing");

  TCanvas * cPMIng = new TCanvas();
  Output_PMIng->DrawNormalized();
  PDF_PMIng->SetLineColor(2);
  PDF_PMIng->Draw("same");
  TCanvas * cPMSci = new TCanvas();
  Output_PMSci->DrawNormalized();
  PDF_PMSci->SetLineColor(2);
  PDF_PMSci->Draw("same");
  TCanvas * cIng = new TCanvas();
  Output_Ing->DrawNormalized();
  PDF_Ing->SetLineColor(2);
  PDF_Ing->Draw("same");
  
}

{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Jan 10 05:08:59 2012) by ROOT version5.26/00c
   TCanvas *c1 = new TCanvas("c1", "c1",0,36,600,400);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   
   TH1F *hrat = new TH1F("hrat","hrat",70,0,1.4e+13);
   hrat->SetBinContent(6,2.897569e-14);
   hrat->SetBinContent(8,1.637573e-14);
   hrat->SetBinContent(9,1.493377e-14);
   hrat->SetBinContent(10,1.584563e-14);
   hrat->SetBinContent(11,1.520555e-14);
   hrat->SetBinContent(12,1.517352e-14);
   hrat->SetBinContent(15,1.866239e-14);
   hrat->SetBinContent(16,1.605333e-14);
   hrat->SetBinContent(17,1.515253e-14);
   hrat->SetBinContent(18,1.484844e-14);
   hrat->SetBinContent(19,1.542421e-14);
   hrat->SetBinContent(20,1.528746e-14);
   hrat->SetBinContent(21,1.4778e-14);
   hrat->SetBinContent(23,7.545561e-15);
   hrat->SetBinContent(24,1.566394e-14);
   hrat->SetBinContent(25,1.644795e-14);
   hrat->SetBinContent(26,1.7037e-14);
   hrat->SetBinError(6,9.658563e-15);
   hrat->SetBinError(8,7.594061e-16);
   hrat->SetBinError(9,6.505252e-16);
   hrat->SetBinError(10,1.584563e-14);
   hrat->SetBinError(11,2.165821e-16);
   hrat->SetBinError(12,3.066766e-16);
   hrat->SetBinError(15,1.866239e-14);
   hrat->SetBinError(16,5.15176e-16);
   hrat->SetBinError(17,2.113301e-16);
   hrat->SetBinError(18,4.409346e-16);
   hrat->SetBinError(19,1.177729e-16);
   hrat->SetBinError(20,1.112146e-16);
   hrat->SetBinError(21,9.043987e-16);
   hrat->SetBinError(23,4.356432e-15);
   hrat->SetBinError(24,2.813327e-15);
   hrat->SetBinError(25,1.084546e-15);
   hrat->SetBinError(26,1.249214e-15);
   hrat->SetEntries(2.734107e-13);
   hrat->Draw("E");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}

{

  TChain * h1 = new TChain("h1","h1");
  h1->Add("test1.root");
  h1->Add("test2.root");
  
  TChain * tree = new TChain("tree","tree");
  tree->Add("test1.root");
  tree->Add("test2.root");

  TFile * f = new TFile("f.root","recreate");
  h1->Write();
  tree->Write();
  f->Close();
  
  //test
  TFile * f1 = new TFile("f.root","read");
  TTree * tree2 = (TTree*) f1->Get("tree");
  cout<<tree2->GetEntries()<<endl;

}

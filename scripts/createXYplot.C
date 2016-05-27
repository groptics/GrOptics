
void createXYplot() {
  TFile *f = new TFile("photonLocation.root");
  TTree *t1 = (TTree*)f->Get("T1");
  TCanvas *c1 = new TCanvas("c1");
  t1->Draw("photonY:photonX","","", 1, 0);

  TCanvas *c2 = new TCanvas("c2");
  TTree *t2 = (TTree*)f->Get("T2");
  t2->Draw("photonY:photonX","","", 1, 0);


}

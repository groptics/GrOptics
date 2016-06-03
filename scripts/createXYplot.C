// script for plotting photonX vs photonY for four telescope images
// requires single shower in photonLocation.root from grOptics
//  used for testing grOptics 

{
  TCanvas *c1 = new TCanvas("c1");
  TFile *f = new TFile("photonLocation.root");

  TTree *t1 = (TTree*)f->Get("T1");
  if (t1 != 0) {
    t1->Draw("photonY:photonX","","", 1, 0);
  }
  
  TTree *t2 = (TTree*)f->Get("T2");
  if (t2 != 0) {
    t2->Draw("photonY:photonX","","same", 1, 0);
  }
  
  TTree *t3 = (TTree*)f->Get("T3");
  if (t3 != 0) {
    t3->Draw("photonY:photonX","","same", 1, 0);
  }
  
  TTree *t4 = (TTree*)f->Get("T4");
  if (t4 != 0) {
    t4->Draw("photonY:photonX","","same", 1, 0);
  }
  
  TTree *t5 = (TTree*)f->Get("T5");
  cout << "t5 " << t5 << endl;
  if (t5 != 0) {
    t5->Draw("photonY:photonX","","same", 1, 0);
  }
}

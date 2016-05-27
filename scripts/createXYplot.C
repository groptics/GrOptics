// script for plotting photonX vs photonY for four telescope images
// requires single shower in photonLocation.root from grOptics
//  used for testing grOptics 

{
  TCanvas *c1 = new TCanvas("c1");
  TFile *f = new TFile("photonLocation.root");
 
  TTree *t1 = (TTree*)f->Get("T1");
  t1->Draw("photonY:photonX","","", 1, 0);

  TTree *t2 = (TTree*)f->Get("T2");
  t2->Draw("photonY:photonX","","same", 1, 0);

  TTree *t3 = (TTree*)f->Get("T3");
  t3->Draw("photonY:photonX","","same", 1, 0);

  TTree *t4 = (TTree*)f->Get("T4");
  t4->Draw("photonY:photonX","","same", 1, 0);

}

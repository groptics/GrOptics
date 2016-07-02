// createXYplot.C
/*
VERSION4.0 
30May2016
      C. Duke
  Grinnell College

 Usage: execute  "root scripts/createXYplot.C"
 useful for testing grOptics

 script for plotting photonX vs photonY for up to five telescope images.
 requires single shower in photonLocation.root from grOptics
 can easily be extended to more than five telescopes.
  
 NOTE: The units are degrees on the sky.  I've made no attempt to fix up the
 axis labels and the axis ranges.  Maybe later.
*/
{
  Double_t plateSFacSC = 97.5;
  Double_t plateSFacDC = 209.5;
  
  int w = 600;
  int h = 600;
  TCanvas *c1 = new TCanvas("c1","c1",w,h);
 
  //c1->SetCanvasSize(w,h);
  
  TFile *f = new TFile("photonLocation.root");

  vector<double> *telPlateSFacVector;
  vector<int> *telTypeVector;
  
  TTree *allT = (TTree*)f->Get("allT");
  allT->SetBranchAddress("telPlateSFacVector",&telPlateSFacVector);
  allT->SetBranchAddress("telTypeVector",&telTypeVector);
  allT->GetEntry(0);
  
  int numTel = (*telPlateSFacVector).size();
  cout << "Number of telescopes " << numTel << endl;
  cout << " Tel   TelType   TelPlateScaleFactor " << endl;
  for (int i = 0;i<4;i++) {
    cout << "  " << i+1 << "      " << (*telTypeVector)[i] << "          " <<  (*telPlateSFacVector)[i] << endl;
  }

  ostringstream convert;
  for (int i= 0;i<numTel;i++) {
    string telNumS;
    convert.str("");
    convert << i+1;
    telNumS = "T" + convert.str();

    TTree *tt = (TTree*)f->Get(telNumS.c_str());
    string ss = "same";
    if (i == 0) {
      ss = "";
    }
    double plf = (*telPlateSFacVector)[i];
    // CONVERT TO STRING
    convert.str("");
    convert << "photonY/"<< plf << ":photonX/" << plf;
    cout << "string in TTree->Draw " << convert.str() << endl;
    string fString = convert.str();
    tt->Draw(fString.c_str(),"",ss.c_str(), 1, 0);
  }
}
  
 

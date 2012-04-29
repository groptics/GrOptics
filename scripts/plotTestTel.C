// script to produce plots from testTel.root
// assumes spot plots for field angle every 1.0 degrees
//   C. Duke
//   Grinnell College
//   Feb. 13 2012

void plotTestTel() {
  cout << " -- starting plotTestTel " << endl;

  // base image file name and file type
  string baseFilename = "psfHistSCCameraWT2.0IR1.5";
  string imageFileType = "png";

  TFile *fi = new TFile("testTel.root");
  
  TH2D *histSpot[4];
  TH2D *histSpot[0] = (TH2D*)fi->Get("hist0"); 
  TH2D *histSpot[1] = (TH2D*)fi->Get("hist2"); 
  TH2D *histSpot[2] = (TH2D*)fi->Get("hist4"); 
  TH2D *histSpot[3] = (TH2D*)fi->Get("hist6"); 

  TH2D *histTime[4];
  TH2D *histTime[0] = (TH2D*)fi->Get("histT0"); 
  TH2D *histTime[1] = (TH2D*)fi->Get("histT2"); 
  TH2D *histTime[2] = (TH2D*)fi->Get("histT4"); 
  TH2D *histTime[3] = (TH2D*)fi->Get("histT6"); 

  double rmsX[4];
  double rmsY[4];
  double deg[4];
  for (int i = 0;i<4;i++) {
    deg[i] = (double)i ;
    cout << "deg " << i << "  " << deg[i] << endl;
  }

  TGraph* graRMS = new TGraph;
  graRMS->SetTitle(";Field angle (deg);2 #times max{RMS_{sagital}, RMS_{tangential}} (arcmin)");

  // make spot plots for 0.0, 1.0, 2.0, and 3.0 degrees
  TCanvas* canSpot = new TCanvas("canSpot", "canSpot", 800, 200);
  canSpot->Divide(4, 1);

  //TCanvas* canTime = new TCanvas("canTime", "canTime", 600, 600);
  //canTime->Divide(2, 2);
  TString str0 = "#scale[1.5]{";
  TString str2 = "}";

  for (int i = 0;i<4;i++) {
    canSpot->cd(i+1);
    histSpot[i]->SetStats(0);
    histSpot[i]->SetMarkerColor(4);

    rmsX[i] = histSpot[i]->GetRMS(1);
    rmsY[i] = histSpot[i]->GetRMS(2);
    graRMS->SetPoint(graRMS->GetN(),deg[i], (rmsX[i] > rmsY[i] ? rmsX[i]: rmsY[i])*2);

    histSpot[i]->GetXaxis()->SetTitleSize(0.05);
    histSpot[i]->GetXaxis()->SetTitleOffset(0.8);
    histSpot[i]->GetYaxis()->SetTitleSize(0.05);
    histSpot[i]->GetYaxis()->SetTitleOffset(0.8);

    TString str1 = Form("#it{#theta} = %3.1f (deg)",deg[i]);
    TString str = str0 + str1 + str2;
    cout << "str " << str << endl;
    //gStyle->SetTitle
      //TString title = "new title";
    //title->SetTextSize(2.0);
    histSpot[i]->SetTitle(str);
    //histSpot[i]->SetTitle(#scale[1.2]{'Form("#it{#theta} = %3.1f (deg)",deg[i])'} );
    //histSpot[i]->SetTitle(#scale[1.2]{'Form("#it{#theta} = %3.1f (deg)",deg[i])}'  ;
    gPad->SetGrid();
    
    cout << " deg rmsX rmsY " << deg[i] << " " << rmsX[i] 
	 << "  " << rmsY[i] << endl;
    histSpot[i]->Draw("col");
  }
  string imageFilename = baseFilename + "_PsfHist." + imageFileType;
  canSpot->SaveAs(imageFilename.c_str());

  //for (int i = 0;i<4;i++) {
    // space here for time histograms, add later.
  //}


  double rangeGraphPsfX = 4.0;
  double rangeGraphPsfY = 6.0;
  TCanvas* canGrPSF = new TCanvas("canGrPSF", "canGrPSF", 600, 300);
  gPad->SetGrid();
  graRMS->SetMarkerStyle(25);
  graRMS->GetXaxis()->SetLimits(0, rangeGraphPsfX);
  graRMS->GetYaxis()->SetRangeUser(0, rangeGraphPsfY);
  graRMS->GetXaxis()->SetTitleSize(0.05);
  graRMS->GetXaxis()->SetTitleOffset(0.8);
  graRMS->GetYaxis()->SetTitleSize(0.05);
  graRMS->GetYaxis()->SetTitleOffset(0.8);
  
  graRMS->Draw("apl");
  string graphFile = baseFilename + "_graph." + imageFileType;
  canGrPSF->SaveAs(graphFile.c_str());
  // TBrowser *b = new TBrowser();

}

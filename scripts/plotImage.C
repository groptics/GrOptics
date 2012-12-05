/* root script to superimpose photonx vs -photony camera plots for a single
   shower in the photonLocation.root output file of grOptics. The camera
   units are mm. Note the "-photonY" to plot in telescope coor. system.

   Used to test images produced from GrOptics. To test, set the wobble offset
   in opticsSimuation.pilot and turn on the DEBUGBRANCHES option. This option
   provides core location and telescope location information for the telescope
   and shower planes in each output telescope tree.  So, you can have a look
   at the core location compared to the telescope location, plot the images
   using this script and see if the images point to the source point.

   Note that the Az/Zn values of the primary and telescope determine the 
   orientation of the image on the camera.

   Note that the plots are of (photonX,-photonY) to produce plots in 
   telescope coordinates rather than in camera coordinates (easier to 
   compare image orientation with core location which is also in 
   telescope coordinates.

   This script assumes "photonLocation.root" for the name of the GrOptics 
   output file.

   To run this script:
   First, set the title of the plot, I use, e.g.  
   string plotTitle = "Images, tel coor: (wE,wN,lat) = (0.0,0.0,90.0";
   so the wobble offsets and latitudes are given.

   Then, set the name of the output file, e.g. image.png. The code will 
   produce the image on the screen and save the canvas to this file.
   Then, execute:
        root plotImage.C

   Charlie Duke
   Grinnell College
   December 6, 2012
 */
{
  cout << " entering plotImage" << endl;
  gROOT->Reset();

  string plotTitle ="Camera Images, tel coor: (wE,wN,lat) = (0.5,0.0,70.0)"; 
  string CanvasFilename = "cameraImage_5.png";

  TCanvas *c1 = new TCanvas("c1","camera Image",600,600);
  c1->SetTitle("Camera Images - tel.coor.");
  gPad->SetGrid();
  gStyle->SetOptStat(0);

  TFile *fi = new TFile("photonLocation.root");

  TTree *t1 = (TTree*)fi->Get("T1");
  TTree *t2 = (TTree*)fi->Get("T2");
  TTree *t3 = (TTree*)fi->Get("T3");
  TTree *t4 = (TTree*)fi->Get("T4");

  TH2D *h1 = new TH2D("h1","h1",600,-300,300,600,-300,300);
  TH2D *h2 = new TH2D("h2","h2",600,-300,300,600,-300,300);
  TH2D *h3 = new TH2D("h3","h3",600,-300,300,600,-300,300);
  TH2D *h4 = new TH2D("h4","h4",600,-300,300,600,-300,300);

  h1->SetTitle(plotTitle.c_str());
  h1->GetXaxis()->SetTitle("camera X axis (mm)");
  h1->GetYaxis()->SetTitle("camera Y axis (mm):tel.coorSys");
  h1->GetXaxis()->SetLabelSize(0.025);
  h1->GetYaxis()->SetLabelSize(0.025);

  h1->SetMarkerColor(1);
  h2->SetMarkerColor(2);
  h3->SetMarkerColor(3);
  h4->SetMarkerColor(4);

  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);

  h1->SetLineWidth(4);
  h2->SetLineWidth(4);
  h3->SetLineWidth(4);
  h4->SetLineWidth(4);
  
  t1->Draw("-photonY:photonX >> h1","","", 1, 0);
  t2->Draw("-photonY:photonX >> h2","","", 1, 0);
  t3->Draw("-photonY:photonX >> h3","","", 1, 0);
  t4->Draw("-photonY:photonX >> h4","","", 1, 0);

  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");

  TLegend *leg = new TLegend(0.68,0.80,0.95,0.93);
  leg->AddEntry(h1,"Telescope 1","L");
  leg->AddEntry(h2,"Telescope 2","L");
  leg->AddEntry(h3,"Telescope 3","L");
  leg->AddEntry(h4,"Telescope 4","L");
  leg->Draw();

  c1->SaveAs(CanvasFilename.c_str()); 
}

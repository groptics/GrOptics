/*  script to test adding roughness to mirrors
    Start with a given unit vector 
    and add a linear error or a Gaussian error.

    resulting in a unit vector with a new direction.

    see comments at the end of this file for more details.

    Feb. 27, 2012
    C. Duke
    Grinnell College
 */

void addError() {

  // select linear or Gaussian offset
  bool doLinear = false;  // if false, use Gaussian roughness

  // select output root file (histograms in this root file)
  string rootfileOut = "addError.root";

  // choose direction cosines for the unit vector
  Double_t xcos0 = 0.3;
  Double_t ycos0 = 0.2;
  Double_t zcos0 = sqrt( 1 - xcos0*xcos0 - ycos0*ycos0);
  ROOT::Math::XYZVector vec(xcos0,ycos0,zcos0);

  // choose a maximum offset angle
  Double_t maxAngle = 0.9;  // degrees
  Double_t maxAngleR = maxAngle*(TMath::DegToRad());

  // choose a rms sigma for the Gaussian
  Double_t sigma = 0.2; // degrees
  Double_t sigmaR = sigma*(TMath::DegToRad());

  // set the number of points
  Int_t numPts = 100000;  // number of histogram points
  TFile *fi = new TFile(rootfileOut.c_str(),"RECREATE"); // fileout.root

  /////////////////// Set up histograms /////////////////
  // histogram binning
  Int_t nbinsxy = 50;
  // abs(upper and lower limits) > max.offset angle
  Double_t xylow = -1.0;
  Double_t xyUp  = 1.0;
  TH2D *his = new TH2D("his","offset x/y (degrees)",
		       nbinsxy,xylow,xyUp,nbinsxy,xylow,xyUp);
  his->GetXaxis()->SetTitle("x offset (deg)");
  his->GetYaxis()->SetTitle("y offset (deg)");

  // set binning for annular area histogram, 
  // two histograms will be produced.
  //      annular area histogram
  //      histogram of radius angle
  // binsHH > maxAngle
  Int_t nbinsH = 50;
  Double_t binsHL = 0.0;
  Double_t binsHH = 1.0;
  TH1D *h1 = new TH1D("h1","offset density vs r",nbinsH,binsHL,binsHH);
  h1->GetXaxis()->SetTitle("annulus radius r");
  h1->GetYaxis()->SetTitle("density within annulus (r and r + binW)");
  h1->SetStats(kFALSE);
  
  TH1D *harea = new TH1D("harea","areaHist",nbinsH,binsHL,binsHH);
  Double_t binW = h1->GetBinWidth(3);

  // fill area histogram
  for (int i = 0;i < nbinsH+1;i++) {
    Double_t rcenter = h1->GetBinCenter(i);
    Double_t area = (TMath::TwoPi())*rcenter*binW;
    harea->SetBinContent(i,area);
  }

  // some input details to the screen  
  cout << "  -- addError script";
  if (doLinear) {
    cout << "  with linear error offsets " << endl;;
  }
  else {
    cout << "  with Gaussian error offsets" << endl;
  }
  cout << "       Direction of initial unit vector (dirCosx/y/z)  " ;
  cout << xcos0 << "  " << ycos0 << "  " << zcos0 << endl;
  
  cout << "       maxAngle (Deg)       " << maxAngle << endl;
  
  if (!doLinear) {
    cout << "       Gaussian sigma (Deg) " << sigma << endl;
  }
  /////////////////// start loop ////////////////////  
  TRandom3 TR3;
  TR3.Rndm();

  for (int i = 0;i<numPts;i++) {
    
    Double_t r;
    Double_t theta;
    // set to 1 to add Gaussian error
    if (doLinear) {
      
      // linear error,see notes at end of this script for method
      r = maxAngleR*sqrt(TR3.Rndm());
      theta   = ( TMath::TwoPi())* TR3.Rndm() - (TMath::Pi());
    }
    else {
      // Gaussian error,see notes at end of this script for method
      Double_t n = TR3.Rndm();
      Double_t R = maxAngleR;
      Double_t s = sigmaR;
      Double_t a = 1/(1 - exp(-R*R/(2.*sigmaR*sigmaR))  );
      
      Double_t tmp = 1 - n/a;
      r = s * sqrt( -2.* log(tmp) );
      
      theta   = ( TMath::TwoPi())* TR3.Rndm() - (TMath::Pi());
      
    }  // doLinear

    // get direction cosines of offset vector
    Double_t xcosZNew = sin(r)*cos(theta);
    Double_t ycosZNew = sin(r)*sin(theta);
    Double_t zcosZNew = cos(r);

    // Get x/y and radius in x/y plane in Degrees
    Double_t xcosZNewD = sin(r)*cos(theta)*TMath::RadToDeg();
    Double_t ycosZNewD = sin(r)*sin(theta)*TMath::RadToDeg();;
    Double_t RadN = sqrt(xcosZNewD*xcosZNewD + ycosZNewD*ycosZNewD);

    // fill his with x/y components of offset vector
    his->Fill(xcosZNewD,ycosZNewD);

    // fill h1 with radial distance in degrees in x/y plane
    h1->Fill(RadN);

    // n1 is the new offset unit vector
    ROOT::Math::XYZVector n1(xcosZNew,ycosZNew,zcosZNew);
    
    // get the phi,theta components of the original unit vector
    Double_t phiR = vec.Phi();
    Double_t thetaR = vec.Theta();
    
    // set up rotation matrix
    ROOT::Math::RotationZ rz(-((TMath::PiOver2()) -phiR));
    ROOT::Math::RotationX rx(-thetaR);
    ROOT::Math::Rotation3D rotM = rz*rx;
    
    // change coordinate system, rotated n1 replaces the original
    // given unit vector.
    vec = rotM*n1;
    
    // new direction cosines, print this if necessary.
    if (0) {
      cout << "   new direction cosines " << vec.X() << "  " 
	   << vec.Y() << "  " << vec.Z() << endl;
    }
  }  // loop finished
 ///////////////////////////////////////////////////////////
  // make density histogram by dividing each radius bin by the
  // corresponding annular area.
  h1->Divide(harea);

  TCanvas *cn = new TCanvas("cn","addError", 800,400);
  cn->Divide(2,1);
  cn->cd(1);
  
  if (!doLinear) {
    // make Gaussian fit, results printed to screen (not on graph).
    // fit will show as red line on graph
    h1->Fit("gaus","R","",0.0,.9);
  }
  else {
    // fit a horizontal line, results printed to screen
    // fit will show as red line on graph
    h1->Fit("pol0","R","",0.0,.9);
  }

  cn->cd(2);
  if (doLinear) {
    // the stats are not useful in the linear case
    // but they are in the Gaussian case
    his->SetStats(kFALSE);
  }
  his->Draw();

  // write histograms to root file and we're done
  h1->Write();
  his->Write();
  
}

/*
  Linear deviation

  it's easiest to place a random point with a circle using polar coordinates.
  In this case, find using r = rmax*sqrt(uniform rand.num between 0 and 1).
  Can derive expression from probability transformation, given that we want
  points having constant density within the circle bounded by rmax.

  Alternative method: random numbers to determine an x and a y. Repeat if
  radius is greater than rmax.  The above method is faster; the point is 
  always within rmax and only one random number is required.
 */

/*
  Gaussian deviation
  Ensure a Gaussian distribution in two dimensions by finding x < rmax and 
  y < rmax from the usual Gaussian deviate (both having same rmax and same
  variance). But, doing so requires testing if 
  x*x + y*y < rmax*rmax and repeating if not.  Using a 
  modified Rayleigh distribution that always produces a point 0 < r < rmax 
  uses 1/2 the computing time.  

  Here's an outline of the derivation.

  (A standard Gaussian deviate for r doesn't apply since the circumference
  of the circle with radius r increases with r, i.e. the area of the 
  annulus between r and r + dr increases with r. The standard Gaussian 
  deviate does not compensate for this increase in area while the 
  Rayleigh distribution can be exactly derived using the annulus area
  element. 

  n = uniform random deviate
  R = rmax 
  s = sigma of Gaussian distribution
  a = normalization constant

  pdf = a * ( 1/(2*s*s) )* r * exp(r*r/(2*s*s))
  
  normalize pdf with limits between 0 and R to get "a"

  a = 1/(1 - exp(-R*R/(2*s*s)) )
  
  do probability transformation with uniform deviate pdf and
  solve for r.
  
  r = s* sqrt( -2*ln(1 - n/a) )

  note that the usual expression for R infinitely large (see
  gsl code for example) is
      r = s*sqrt(-2*ln(n) ) with n > 0,
  but we want r to be less than R.
 
  In the specific case of reflection from a mirror, if we assume that 
  the deviation of the normal to the mirror is small and follows two 
  independent Gaussian deviations in the normal plane, then the above
  Rayleigh distribution applies. 
 */

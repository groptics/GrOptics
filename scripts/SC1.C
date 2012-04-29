// $Id$
// Author: Akira Okumura 2010/11/28
// revised for os8 sc structure
// C. Duke 6Feb2011
/******************************************************************************
 * Copyright (C) 2006-, Akira Okumura                                         *
 * All rights reserved.                                                       *
 *****************************************************************************/

// An example of Schwarzschild-Couder optical system
// Parameters are taken from the OS2 configuration in
// Vassiliev, et al. (2007) Astropart. Phys. 28 10-27

// define useful units
static const Double_t cm = AOpticsManager::cm();
static const Double_t mm = AOpticsManager::mm();
static const Double_t um = AOpticsManager::um();
static const Double_t nm = AOpticsManager::nm();
static const Double_t  m = AOpticsManager::m();

void SC1()
{

  TGeoManager *manager = new AOpticsManager("manager","My 3D Project");
  
  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  vacuum->SetTransparency(80);   
  TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);

  /////////////////////////// top volume ///////////////////////////////////
  double fTX = 15.0;  //!< top volume dimensions
  double fTY = 15.0;
  double fTZ; = 15.0
  TGeoVolume *top = manager->MakeBox("top",Air,fTX*m,fTY*m,fTZ*m);
  manager->SetTopVolume(top);
  manager->SetTopVisible(0); 
 
   ////////////// telescope parameters ////////////////////
   // Note: The center of the primary mirror is located at (X, Y, Z) = (0, 0, 0)
 
  double kFp      = 5.5863*m;   // Fp
  double kDp      = 4.8319*2*m; // Dp
  double kDpinner = 2.1933*2*m; // Dpinner
  double kZp      = 0*m;        // Primary mirror's position
  double kZs      = 3./2.*kFp;  // Secondary mirror's position
  double kDs      = 2.7083*2*m; // Ds
  double kDsinner = 0.3945*2*m; // Dsinner
  double kZf      = (3./2. - 1./3.)*kFp;  // Focal plane position
  double kk1      = -19.75*mm; // kappa1
  double kk2      = 0.503*mm; // kappa2

  // The actual dead area should be considered carefully
  // See the spec sheet of H8500
  // const double kPixelSize = 5.8*mm; // width of a MAPMT pixel
  // We do not care the dead area for now

  double kPixelSize   = 6.08*mm; // width of an MAPMT pixel
  double kPixelPitch  = 6.08*mm; // pitch of MAPMT channels
  double kMAPMTWidth  = 52.0*mm;
  double kMAPMTLength = 32.7*mm; // between input window and anode pins
  double kInputWindowThickness = 1.5*mm;

  double kMAPMTAngularSize = 4.00/60.; // (deg)

   // 14th order polynomials
   // fitted with gnuplot by Akira
   const int kNPar = 8;
   
   double kzp[kNPar] = {TMath::Power(kFp,  1)*8.4374e-06,
   TMath::Power(kFp, -1)*0.110917,
   TMath::Power(kFp, -3)*-0.00511208,
   TMath::Power(kFp, -5)*-0.0118961,
   TMath::Power(kFp, -7)*0.0253067,
   TMath::Power(kFp, -9)*-0.0460152,
   TMath::Power(kFp,-11)*0.0413689,
   TMath::Power(kFp,-13)*-0.0172745};

   double kzs[kNPar] = {TMath::Power(kFp,  1)*6.45608e-08,
			      TMath::Power(kFp, -1)*-0.416688,
			      TMath::Power(kFp, -3)*-0.144035,
			      TMath::Power(kFp, -5)*0.647955,
			      TMath::Power(kFp, -7)*-2.96087,
			      TMath::Power(kFp, -9)*9.39256,
			      TMath::Power(kFp,-11)*-18.0811,
			      TMath::Power(kFp,-13)*15.3711};

 //////////////// Make the primary mirror ////////////////////////
   AGeoAsphericDisk* primaryV = 0;
   primaryV = new AGeoAsphericDisk("primaryV", kZp + kzp[0] - 1*um, 
				   0, kZp + kzp[
   0] , 0, kDp/2., kDpinner/2.); // curvatures are set to infinity (0 = 1./inf)
   primaryV->SetPolynomials(kNPar - 1, &kzp[1], kNPar - 1, &kzp[1]);
   

   // need to set reflection graph here, after each mirror.
   // have default 1.0
   AMirror* primaryMirror = new AMirror("primaryMirror", primaryV);
   if (gPrimRefl!=0) {
     primaryMirror->SetReflectivity(gPrimRefl);
   }

   top->AddNode(primaryMirror, 1);

   ////////////////// make secondary mirror /////////////////////
   AGeoAsphericDisk* secondaryV = 0;
   secondaryV = new AGeoAsphericDisk("secondaryV", 
				     kZs + kzs[0], 0, 
				     kZs + kzs[0]  + 1*um, 
				     0, kDs/2., kDsinner/2.); 
   // curvatures are set to infinity (0 = 1./inf)
   secondaryV->SetPolynomials(kNPar - 1, &kzs[1], kNPar - 1, &kzs[1]);
   
   AMirror* secondaryMirror = new AMirror("secondaryMirror", secondaryV); 
   if (gSeconRefl!=0) {
     secondaryMirror->SetReflectivity(gSeconRefl);
   }
   //secondaryMirror->SetLineColor(kBlue);
   //secondaryMirror->SetFillColor(kBlue);
  top->AddNode(secondaryMirror, 1);

  /////////////////////////////////////////////////////
  ////// make Dummy transparent material //////////////
  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(0); // needed in OpenGL view
  TGeoMedium* med = new TGeoMedium("med", 1, mat);

  //////////////////////////////////////////////////////////////////
  ///// Make a dummy obscuration before the secondary mirror ////////
  AGeoAsphericDisk* secondaryObsV = 0;

  secondaryObsV = new AGeoAsphericDisk("secondaryObsV", kZs + kzs[0] + 10*um, 
				       0, kZs + kzs[0] + 1*um + 10*um, 
				       0, kDs/2., 0);

  secondaryObsV->SetPolynomials(kNPar - 1, &kzs[1], kNPar - 1, &kzs[1]);

  AObscuration* secondaryObs = new AObscuration("secondaryObs", secondaryObsV, med);

  top->AddNode(secondaryObs, 1);

  bool makeMAPMTCamera = false;

  if (makeMAPMTCamera ) {
    //////////// Make a single MAPMT photocathode /////////////////////
    //a very thin box
    TGeoBBox* mapmtCathodeV = new TGeoBBox("mapmtCathodeV", kPixelSize/2., 
					   kPixelSize/2., 100*um); 

    AFocalSurface* mapmtCathode = new AFocalSurface("mapmtCathode", 
						    mapmtCathodeV);

    //////////////// Make a single MAPMT //////////////////////////
    TGeoBBox* mapmtV = new TGeoBBox("mapmtV", kMAPMTWidth/2., kMAPMTWidth/2.,
				    kMAPMTLength/2.);
    AOpticalComponent* mapmt = new AOpticalComponent("mapmt", mapmtV);
    TGeoBBox* mapmtInputWindowV = new TGeoBBox("mapmtInputWindowV",
					       kMAPMTWidth/2., kMAPMTWidth/2.,
					       kInputWindowThickness/2.);
    
    ALens* mapmtInputWindow = new ALens("mapmtInputWindow", 
					mapmtInputWindowV, med);
    
    // The refractive index of the input window must be checked
    ARefractiveIndex* bk7 = AGlassCatalog::GetRefractiveIndex("N-BK7");
    mapmtInputWindow->SetRefractiveIndex(bk7);
    mapmt->AddNode(mapmtInputWindow, 1, new TGeoTranslation(0, 
							    0, 
		     kMAPMTLength/2. - kInputWindowThickness/2.));

    int n = 1;
    for(int i = 0; i < 8; i++){
      double dx = (i - 3.5)*kPixelPitch;
      for(int j = 0; j < 8; j++){
	double dy = (j - 3.5)*kPixelPitch;
	mapmt->AddNode(mapmtCathode, n, new TGeoTranslation(dx, dy, kMAPMTLength/2. - kInputWindowThickness - 100*um));
	n++;
      } // j
    } // i
    
    TGeoBBox* mapmtBackObsV = new TGeoBBox("mapmtBackObsV",
					   kMAPMTWidth/2., kMAPMTWidth/2.,
					   15*mm);
    AObscuration* mapmtBackObs = new AObscuration("mapmtBackObs", mapmtBackObsV);
    mapmt->AddNode(mapmtBackObs, 1, new TGeoTranslation(0, 
							0, 
							-kMAPMTLength/2. + 15*mm));

    /////////////////// Make the focal plane ///////////////////////
    n = 1;
    ofstream ofile("mapmtLocations.txt");
    bool prtMAPMT = false;
    if (prtMAPMT) {
      ofile << " -- GSCTelescope::buildTelescope, prtMAPMT " << endl;
      ofile << "               all distances in centimeters" << endl;
      ofile << "        kMAPMTWidth / kMAPMTLength " << kMAPMTWidth 
	    << "  /  " << kMAPMTLength << endl;
      ofile << "        i   j     dx        dy      dz    r" << endl;
    }
    for(int i = -7; i <= 7; i++){
      double dx = i*kMAPMTWidth;
      for(int j = -7; j <= 7; j++){
	if((TMath::Abs(i) + TMath::Abs(j) >= 11) || (TMath::Abs(i)*TMath::Abs(j) == 21)){
	  continue;
	} // if
	double dy = j*kMAPMTWidth;
	double r2 = (i/7.)*(i/7.) + (j/7.)*(j/7.);
	double dz = kk1*r2 + kk2*r2*r2;

	if (prtMAPMT) {
	  ofile << "      " << i << "   " << j  << "      " << dx
		<< "     " << dy << "     " << dz << "     "
		<< sqrt(dx*dx + dy*dy) << endl;
	}
	
	top->AddNode(mapmt, n, new TGeoTranslation(dx, dy, kZf - kMAPMTLength/2. + dz));
	n++;
      } // y
    } // x
  }
  else {
    *oLog << " ready to make focal plane with no camera " << endl;
    *oLog << " kZf " << kZf << endl;
    *oLog << " kZf - 1*mm " << kZf-1*mm << endl;
    *oLog << " 38.*cm " << 38.*cm << endl;
    *oLog << " kMAPMTWidth " << kMAPMTWidth << endl;
   *oLog << " kk1  kk2 " << kk1 << "  " << kk2 << endl;
    Double_t tmp = 1/(7.*kMAPMTWidth);
    Double_t tmp2 = tmp*tmp;
    Double_t beta[2];
    beta[0] = tmp2*kk1;
    beta[1] = tmp2*tmp2*kk2;
    *oLog << "  beta[0]  beta[1] " << beta[0] << "  " << beta[1] << endl;
    //double planeRad = 38.*cm;
    double planeRad = 38.*cm;
    *oLog << " planeRad " << planeRad << endl;

    AGeoAsphericDisk* focalV = new AGeoAsphericDisk("focalV", kZf + kzs[0] - 1*mm, 0, kZf + kzs[0], 0, planeRad, 0.); // curvatures are set to infinity (0 = 1./inf)
    focalV->SetPolynomials(2, beta, 2, beta);
    AFocalSurface* focalPlane = new AFocalSurface("focalPlane", focalV);
    top->AddNode(focalPlane, 1);

    // make an obscuration disk before the focal plane
    AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV",
                        kZf + kzs[0] - 2*mm,0,kZf + kzs[0]-1*mm, 0, planeRad, 0.);
    focalObsV->SetPolynomials(2, beta, 2, beta);

    AObscuration* focalObs = new AObscuration("focalObs", focalObsV);
    /*
    AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV", focalV->CalcF1(10*cm*7) - 1*cm, 0, focalV->CalcF1(10*cm*7)    
  AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV", focalV->CalcF1(10*cm*7) - 1*cm, 0, focalV->CalcF1(10*cm*7), 0, 10*cm*7.3, 0.);
    */
    top->AddNode(focalObs, 1);
    
    *oLog << " focal plane positions: " << kZf + kzs[0] - 1*mm << "  " << kZf + kzs[0] << endl;
    *oLog << " obs plane positions:   " << kZf + kzs[0] - 2*mm << "  " << kZf + kzs[0]-1*mm << endl;

 
  manager->CloseGeometry();

  TCanvas* canGeometry = new TCanvas("canGeometry", "canGeometry", 800, 800);
  top->Draw();

  // Start ray-tracing
  const int kN = 8; // 0 to 7 [deg]
  TH2D* hist[kN];
  TH1D* histT[kN];
  TGraph* graAeff = new TGraph;
  graAeff->SetTitle(";Field angle (deg);Effective Area (m^{2})");
  TGraph* graRMS = new TGraph;
  graRMS->SetTitle(";Field angle (deg);2 #times max{RMS_{sagital}, RMS_{tangential}} (arcmin)");
  TGraph* graT = new TGraph;
  graT->SetTitle(";Field angle (deg);Photon propagation time spread (RMS) (ns)");

  for(Int_t n = 0; n < kN; n++){
    hist[n] = new TH2D(Form("hist%d", n), Form("#it{#theta} = %d (deg);X (arcmin);Y (arcmin)", n), 1000, -20, 20, 1000, -20, 20);
    histT[n]= new TH1D(Form("histT%d",n), Form("#it{#theta} = %d (deg);Propagation delay (ns);Entries", n), 120, -6, 6);

    Double_t deg = n;

    TGeoTranslation* raytr = new TGeoTranslation("raytr", -1.2*kZs*TMath::Sin(deg*TMath::DegToRad()), 0, 1.2*kZs*TMath::Cos(deg*TMath::DegToRad()));

    TVector3 dir;
    dir.SetMagThetaPhi(1, TMath::Pi() - deg*TMath::DegToRad(), 0);
    Double_t lambda = 400*nm; // does not affect the results because we have no lens
    // 1 photon per 0.0025 m^2
    ARayArray* array = ARayShooter::Square(lambda, 20*m, 401, 0, raytr, &dir);
    manager->TraceNonSequential(*array);
    TObjArray* focused = array->GetFocused();

    Double_t Aeff = 0.;
    for(Int_t j = 0; j <= focused->GetLast(); j++){
      ARay* ray = (ARay*)(*focused)[j];
      if(!ray) continue;

      // Calculate the effective area from the number of focused photons
      Aeff += 0.0025; // 0.0025 (m^2)

      Double_t p[4];    
      ray->GetLastPoint(p);
      ray->SetLineWidth(1);
      /* uncomment here if you want to draw all photon trajectories
      TPolyLine3D* pol = ray->MakePolyLine3D();
      pol->Draw();
      */
      Double_t x = deg*10*cm;
      hist[n]->Fill((p[0] - x)/(10*cm)*60, p[1]/(10*cm)*60);
      histT[n]->Fill((p[3] - (3.2*kZs - kZf)/(TMath::C()*m))/1e-9); // ns
    } // j

    graAeff->SetPoint(graAeff->GetN(), deg, Aeff);

    Double_t rmsx = hist[n]->GetRMS(1);
    Double_t rmsy = hist[n]->GetRMS(2);
    
    graRMS->SetPoint(graRMS->GetN(), deg, (rmsx > rmsy ? rmsx : rmsy)*2);

    graT->SetPoint(graT->GetN(), deg, histT[n]->GetRMS());

    delete array;
    delete raytr;
  } // n

  TCanvas* canSpot = new TCanvas("canSpot", "canSpot", 1200, 600);
  canSpot->Divide(4, 2);

  TCanvas* canTime = new TCanvas("canTime", "canTime", 1200, 600);
  canTime->Divide(4, 2);

  for(Int_t i = 0; i < kN; i++){
    canSpot->cd(i + 1);
    hist[i]->Draw("colz");

    canTime->cd(i + 1);
    histT[i]->Draw();
  } // i

  // Figure 5 in the paper
  TCanvas* canFig5 = new TCanvas("canFig5", "canFig5", 1200, 600);
  canFig5->Divide(2, 1);
  canFig5->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graAeff->Draw("apl");
  graAeff->SetMarkerStyle(25);
  graAeff->GetXaxis()->SetLimits(0, 7);
  graAeff->GetYaxis()->SetRangeUser(0, 60);

  // PSF is not consistent with the original paper, but the spot diagram at
  // 5 (deg) is consistent with each other by eye comparison. There may be a
  // difference between calculations of RMS in my code and the paper
  canFig5->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  graRMS->Draw("apl");
  graRMS->SetMarkerStyle(25);
  graRMS->GetXaxis()->SetLimits(0, 7);
  graRMS->GetYaxis()->SetRangeUser(0, 10);

  // Figure 10 in the paper
  // Time spread is 2 times larger in the original paper. I believe the paper
  // is wrong. You can roughly calculate the spread width by
  // Dp * sin(angle)/c ~ 2.5 (ns)
  TCanvas* canFig10 = new TCanvas("canFig10", "canFig10", 1200, 600);
  canFig10->Divide(2, 1);
  canFig10->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graT->Draw("apl");
  graT->SetMarkerStyle(25);
  graT->GetXaxis()->SetLimits(0, 7);
  graT->GetYaxis()->SetRangeUser(0, 1.8);

  canFig10->cd(2);
  histT[5]->Draw();
}

{
  gROOT->Reset();

  const   double cm = 1.;
  const   double mm = 1e-1;
  const   double um = 1e-4;
  const   double nm = 1e-7;
  const double  m = 100.;

  const double kk1      = -19.75*mm;
  const double kk2      = 0.503*mm;

   const double kZf = 0.0;
  TGeoManager *geom = new TGeoManager();

  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);

  TGeoVolume *top = geom->MakeBox("top",Air,300,300,300);
  geom->SetTopVolume(top);

  geom->SetTopVisible(0); 
  
  const double kPixelSize   = 6.08*mm; // width of an MAPMT pixel
  const double kPixelPitch  = 6.08*mm; // pitch of MAPMT channels
  const double kMAPMTWidth  = 52.0*mm;
  const double kMAPMTLength = 32.7*mm; // between input window and anode pins
  const double kInputWindowThickness = 1.5*mm;

  TGeoVolume *mapmtCathode = geom->MakeBox("mapmtCathode",Air,
					   kPixelSize/2., kPixelSize/2., 100*um);
  mapmtCathode->SetLineColor(kBlue);
  //mapmtCathode->SetFillColor(kRed);
  //mapmtCathode->Draw("ogle");
  TGeoVolume *mapmt = geom->MakeBox("mapmt",Air,kMAPMTWidth/2., kMAPMTWidth/2.,
                                  kMAPMTLength/2.);
  //mapmt->SetLineColor(2);

  TGeoVolume *mapmtInputWindow = geom->MakeBox("mapmtInputWindow",Air,
					       kMAPMTWidth/2., kMAPMTWidth/2.,
					       kInputWindowThickness/2.);
  mapmtInputWindow->SetLineColor(2);
  mapmt->AddNode(mapmtInputWindow, 1, new TGeoTranslation(0, 0, kMAPMTLength/2. - kInputWindowThickness/2.));

 int n = 1;
  for(int i = 0; i < 8; i++){
    double dx = (i - 3.5)*kPixelPitch;
    for(int j = 0; j < 8; j++){
      double dy = (j - 3.5)*kPixelPitch;
      mapmt->AddNode(mapmtCathode, n, new TGeoTranslation(dx, dy, kMAPMTLength/2. - kInputWindowThickness - 100*um));
      n++;
    } // j
  } // i

  TGeoVolume *mapmtBackObs = geom->MakeBox("mapmtBackObs",Air,kMAPMTWidth/2., kMAPMTWidth/2.,
					   15*mm);
  mapmtBackObs->SetLineColor(3);
  mapmt->AddNode(mapmtBackObs, 1, new TGeoTranslation(0, 0, -kMAPMTLength/2. + 15*mm));
  //mapmt->SetTransparency(60);
  // Make the focal plane
  n = 1;
  for(int i = -7; i <= 7; i++){
    double dx = i*kMAPMTWidth;
    for(int j = -7; j <= 7; j++){
      if((TMath::Abs(i) + TMath::Abs(j) >= 11) || (TMath::Abs(i)*TMath::Abs(j) == 21)){
        continue;
      } // if
      double dy = j*kMAPMTWidth;
      double r2 = (i/7.)*(i/7.) + (j/7.)*(j/7.);
      double dz = kk1*r2 + kk2*r2*r2;
      top->AddNode(mapmt, n, new TGeoTranslation(dx, dy, kZf - kMAPMTLength/2. + dz));
      n++;
    } // y
  } // x

  //top->AddNode(mapmt,1,new TGeoTranslation(0,0,0));
  //mapmt->Draw("ogle");
  geom->CloseGeometry();
  geom->SetVisLevel(4);
  TCanvas *cCamera = new TCanvas("cCamera","cCamera");
  geom->GetMasterVolume()->Draw("gl");
  mapmt->SetLineColor(4);
  TCanvas *cMapmt = new TCanvas("cMapmt","cMapmt");
  mapmt->Draw("gl");
}

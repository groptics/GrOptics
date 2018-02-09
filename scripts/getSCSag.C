{
  gROOT->Reset();

  // calculate mirror sag at its edge, e.g. r = mirror radius

  const Double_t cm = 1.;
  const   Double_t mm = 1e-1;
  const   Double_t um = 1e-4;
  const   Double_t nm = 1e-7;
  const Double_t  m = 100.;

  const Double_t kFp = 5.5863*m;  // focal length (parameter in coeff.)
  const Double_t fZ1 = 0.0;   // Z location of center of mirror
  const double kDp   =  4.8319*2*m;  // mirror diameter
  const double kDpinner = 2.1933*2*m; // Dpinner

  cout << "kFp " << kFp << endl;
  cout << "kDp " << kDp << endl;


  const Double_t fCurve1 = 0.0; // set curvature to zero.

  // polynomial parameters, note array starts from 1; see
  // GSCTelescope.cpp
  // SetPolynomials(kNPar - 1, &kzp[1], kNPar - 1, &kzp[1]);
  const int kNPar = 8 -1;
  const double kzp[kNPar] = { // TMath::Power(kFp,  1)*8.4374e-06,
			     TMath::Power(kFp, -1)*0.110917,
			     TMath::Power(kFp, -3)*-0.00511208,
			     TMath::Power(kFp, -5)*-0.0118961,
			     TMath::Power(kFp, -7)*0.0253067,
			     TMath::Power(kFp, -9)*-0.0460152,
			     TMath::Power(kFp,-11)*0.0413689,
			     TMath::Power(kFp,-13)*-0.0172745};

  // Calculate z value of surface at edge of mirror
  Double_t r = kDp/2.0;
  //Double_t r = kDpinner/2.0;
  //Double_t r = 100.0;
  cout << "  get sag at this distance from center " << r << endl;

  Int_t fNPol1 = kNPar;

  Double_t p = r*r*fCurve1;
  
  Double_t ret = fZ1 + p/(1 + sqrt(1 - p*fCurve1));
  cout << " retinit " << ret << endl;
  
  for(Int_t i = 0; i < fNPol1; i++){
    ret += kzp[i]*pow(r, 2*(i+1));
    cout << i << " " << kzp[i] << "   " << pow(r, 2*(i+1)) << "  " <<  ret << endl;
  } // i
  cout << " z value (cm) at edge of mirror " << ret << endl;
  //////////////////////////////

  // calculate radius of enclosing sphere (or almost enclosing sphere).
  // radius = distance from center to edge.

  Double_t R = sqrt(r*r + ret*ret);
  cout << "radius of enclosing (almost) sphere " << R << endl;
  
}

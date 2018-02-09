
void getPolyCoeffsFocalPlane() {

  Double_t km() { return 1e3*m();};
  Double_t  m() { return 1e2*cm();};
  Double_t cm() { return 1;};
  Double_t mm() { return 1e-3*m();};
  Double_t um() { return 1e-6*m();};
  Double_t nm() { return 1e-9*m();};

  cout << "  starting getPolyCoeffsFocalPlane " << endl;

  Double_t kMAPMTWidth = 52.0 * mm();
  Double_t kk1 = -19.75*mm(); // mm
  Double_t kk2 = 0.503*mm(); // mm

  cout << "  kMAPMTWidth " << kMAPMTWidth << endl;
  cout << "  kk1         " << kk1 << endl;
  cout << "  kk2         " << kk2 << endl << endl;

  Int_t i = 4;
  Int_t j = 3;
  Double_t dx = i*kMAPMTWidth;
  Double_t dy = j*kMAPMTWidth;
 
  cout << " i  j " << i << "    " << j << endl;
  cout << " dx  dy  " << dx << "   " << dy 
       << "   " <<  endl << endl;

  // calculate sag from GSCTelescope, this agrees with GDCTelescope calc.
  Double_t r2 = (i/7.)*(i/7.) + (j/7.)*(j/7.);
  Double_t dz =  kk1*r2 + kk2*r2*r2;
  cout << "       sag from GSCTelescope r2  dz " << r2 
       << "   " << dz << endl;

  // recalculate z based on expansion in dx2+dy2
  Double_t tmp = 1/(7.*kMAPMTWidth);
  Double_t tmp2 = tmp*tmp;
  Double_t beta1 = tmp2*kk1;
  Double_t beta2 = tmp2*tmp2*kk2;
  Double_t rr2 = dx*dx + dy*dy;
  Double_t zv = beta1*rr2 + beta2*rr2*rr2;
  cout << "  beta1  beta2 " << beta1 << "  " << beta2 << endl;
  cout << "    zv  " << zv << endl;
  cout << "      compare with dz above " << dz << endl;
 
}

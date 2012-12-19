// New Schwartschild-Couder optical system with segmented mirrors and the
// telescope structure.
// Nov/21/2012 Akira Okumura <oxon@mac.com>

#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TVector3.h"

#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include "TGeoTrd1.h"
#include "TGeoTube.h"
#include "TGeoXtru.h"
#include "TPolyLine3D.h"
#include "TRandom.h"
#include "TROOT.h"

#include "ABorderSurfaceCondition.h"
#include "AGeoAsphericDisk.h"
#include "AGlassCatalog.h"
#include "ALens.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ARay.h"
#include "ARayArray.h"
#include "ARayShooter.h"

#include <iostream>

// define useful units
static const Double_t cm = AOpticsManager::cm();
static const Double_t mm = AOpticsManager::mm();
static const Double_t um = AOpticsManager::um();
static const Double_t nm = AOpticsManager::nm();
static const Double_t  m = AOpticsManager::m();

class SegmentedMirror
{
  // The base class for segmented mirrors. This class holds some Parameters
  // to define the shape and alignment errors of a mirror. But the surfaces of
  // the primary and secondary mirrors are given in NewSCT class.

private:

protected:
  Double_t fRmin; // The innner radius of a sector
  Double_t fRmax; // The outer radius of a sector
  Double_t fPhimin; // The minimum azimuthal angle
  Double_t fPhimax; // The maximum azimuthal angle

  Double_t fMargin; // Margin around a mirror
  Double_t fPositionErrorX; // Alignment error along X axis
  Double_t fPositionErrorY; // Alignment error along Y axis
  Double_t fPositionErrorZ; // Alignment error along Z axis
  Double_t fRotationErrorXY; // Rotation error in X-Y plane
  Double_t fRotationErrorSagittal; // Rotation error in sagittal plane
  Double_t fRotationErrorTangential; // Rotation error in tangential plane
  Double_t fRoughness; // Mirror roughness (deg)

public:
  SegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  virtual ~SegmentedMirror() {}

  virtual AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary) = 0;
  TGeoCombiTrans* BuildMirrorCombiTrans(AGeoAsphericDisk* disk, Bool_t isPrimary);

  Double_t GetRoughness() const {return fRoughness;}

  void SetMargin(Double_t margin);
  void SetPositionErrors(Double_t x, Double_t y, Double_t z);
  void SetRotationErrors(Double_t xy, Double_t sag, Double_t tan);
  void SetRougness(Double_t roughness);
};

//______________________________________________________________________________
SegmentedMirror::SegmentedMirror(Double_t rmin, Double_t rmax,
                                 Double_t phimin, Double_t phimax)
  : fRmin(rmin), fRmax(rmax), fPhimin(phimin), fPhimax(phimax),
    fMargin(0), fPositionErrorX(0), fPositionErrorY(0), fPositionErrorZ(0),
    fRotationErrorXY(0), fRotationErrorSagittal(0), fRotationErrorTangential(0),
    fRoughness(0)
{
}

//______________________________________________________________________________
TGeoCombiTrans* SegmentedMirror::BuildMirrorCombiTrans(AGeoAsphericDisk* disk, Bool_t isPrimary)
{
  Double_t cphi = (fPhimax + fPhimin)/2.;
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  TGeoRotation rot1("", -90. + cphi + fRotationErrorXY, fRotationErrorSagittal, 0);
  TGeoRotation rot2("", -cphi, 0., 0.);
  TGeoRotation rot3("", 0., fRotationErrorTangential, 0.);
  TGeoRotation rot4("", cphi, 0., 0.);

  TGeoTranslation tr4(cr*TMath::Cos(cphi*TMath::DegToRad()) + fPositionErrorX, cr*TMath::Sin(cphi*TMath::DegToRad()) + fPositionErrorY, cz + fPositionErrorZ);
 
  TGeoCombiTrans* combi = new TGeoCombiTrans(tr4, ((rot4*rot3)*rot2)*rot1);

  return combi;
}

//______________________________________________________________________________
void SegmentedMirror::SetMargin(Double_t margin)
{
  // Set the margin around the mirror
  fMargin = margin;
}

//______________________________________________________________________________
void SegmentedMirror::SetPositionErrors(Double_t x, Double_t y, Double_t z)
{
  // Set the position errors along X/Y/Z axes
  fPositionErrorX = x;
  fPositionErrorY = y;
  fPositionErrorZ = z;
}

//______________________________________________________________________________
void SegmentedMirror::SetRotationErrors(Double_t xy, Double_t sag, Double_t tan)
{
  // Set the rotation errors in XY/sagittal/tangential planes
  fRotationErrorXY = xy;
  fRotationErrorSagittal = sag;
  fRotationErrorTangential = tan;
}

//______________________________________________________________________________
void SegmentedMirror::SetRougness(Double_t roughness)
{
  // Set the mirror roughness in unit of degree.
  fRoughness = roughness;
}

//______________________________________________________________________________
class SectorSegmentedMirror : public SegmentedMirror
{
private:

public:
  SectorSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~SectorSegmentedMirror(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};

//______________________________________________________________________________
SectorSegmentedMirror::SectorSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax)
  : SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

//______________________________________________________________________________
AMirror* SectorSegmentedMirror::BuildMirror(const char* name,
                                            AGeoAsphericDisk* disk,
                                            Bool_t isPrimary)
{
  Double_t zmin = disk->GetOrigin()[2] - disk->GetDZ();
  Double_t zmax = disk->GetOrigin()[2] + disk->GetDZ();
  Double_t dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // Define the base of the segment shape
  TGeoTubeSeg* seg1 = new TGeoTubeSeg(Form("%s_seg1", name),
                                      fRmin + fMargin, fRmax - fMargin,
                                      (zmax - zmin)/2., -dphi + 90., dphi + 90.);

  // To be used for the margin cut of the side edges
  TGeoTubeSeg* seg2 = new TGeoTubeSeg(Form("%s_seg2", name),
                                      0, fRmax - fMargin,
                                      (zmax - zmin)/2., -dphi + 90., dphi + 90.);

  TGeoTranslation* tr1
    = new TGeoTranslation(Form("%s_tr1", name), 0, -cr, (zmax + zmin)/2. - cz);
  tr1->RegisterYourself();

  Double_t d = fMargin/TMath::Sin(dphi*TMath::DegToRad());
  TGeoTranslation* tr2
    = new TGeoTranslation(Form("%s_tr2", name), 0, d - cr, (zmax + zmin)/2. - cz);
  tr2->RegisterYourself();

  TGeoTranslation* tr3 = new TGeoTranslation(Form("%s_tr3", name), 0, -cr, -cz);
  tr3->RegisterYourself();

  TGeoCompositeShape* comp
    = new TGeoCompositeShape(Form("%s_comp", name),
                             Form("%s:%s*%s:%s*%s:%s",
                                  seg1->GetName(), tr1->GetName(),
                                  seg2->GetName(), tr2->GetName(),
                                  disk->GetName(), tr3->GetName()));

  AMirror* mirror = new AMirror(Form("%s_mirror", name), comp);
  
  return mirror;
}

//______________________________________________________________________________
class TetragonSegmentedMirror : public SegmentedMirror
{
private:

public:
  TetragonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~TetragonSegmentedMirror(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};

//______________________________________________________________________________
TetragonSegmentedMirror::TetragonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax) :
  SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

//______________________________________________________________________________
AMirror* TetragonSegmentedMirror::BuildMirror(const char* name,
                                              AGeoAsphericDisk* disk,
                                              Bool_t isPrimary)
{
  Double_t zmin = disk->GetOrigin()[2] - disk->GetDZ();
  Double_t zmax = disk->GetOrigin()[2] + disk->GetDZ();
  Double_t dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // Define the base of the segment shape
  // ____
  //     ----____
  // ------- B    Rmax  
  // ----p1/
  //    / /
  // --p0/
  // ----____
  //     A   ---- Rmin
  //
  Double_t d2r = TMath::DegToRad();
  Double_t ax = fRmin*TMath::Cos((90. - dphi)*d2r);
  Double_t ay = fRmin*TMath::Sin((90. - dphi)*d2r);
  Double_t bx = fRmax*TMath::Cos((90. - dphi)*d2r);
  Double_t by = fRmax*TMath::Sin((90. - dphi)*d2r);
  Double_t p0x = ax - fMargin*(1./TMath::Sin(dphi*d2r) - 1.)/TMath::Tan((90. - dphi)*d2r);
  Double_t p0y = ay + fMargin;
  Double_t p1x = ax + (-ay + by - fMargin*(1./TMath::Sin(dphi*d2r) + 1.))/TMath::Tan((90. - dphi)*d2r);
  Double_t p1y = by - fMargin;
  TGeoTrd1* trd1 = new TGeoTrd1(Form("%s_trd1", name), p0x, p1x, (zmax - zmin)/2., (p1y - p0y)/2.);

  TGeoRotation* rot1 = new TGeoRotation(Form("%s_rot1", name), 0., -90., 0.);
  rot1->RegisterYourself();

  TGeoCombiTrans* combi1 = new TGeoCombiTrans(Form("%s_combi1", name), 0, (p0y + p1y)/2. - cr, (zmax + zmin)/2. - cz, rot1);
  combi1->RegisterYourself();

  TGeoTranslation* tr2 = new TGeoTranslation(Form("%s_tr2", name), 0, -cr, -cz);
  tr2->RegisterYourself();
  
  TGeoCompositeShape* comp
    = new TGeoCompositeShape(Form("%s_comp", name),
                             Form("%s:%s*%s:%s",
                                  trd1->GetName(), combi1->GetName(),
                                  disk->GetName(), tr2->GetName()));
  
  AMirror* mirror = new AMirror(Form("%s_mirror", name), comp);
  
  return mirror;
}

//______________________________________________________________________________
class PentagonSegmentedMirror : public SegmentedMirror
{
private:

public:
  PentagonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~PentagonSegmentedMirror(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};

//______________________________________________________________________________
PentagonSegmentedMirror::PentagonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax) :
  SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

//______________________________________________________________________________
AMirror* PentagonSegmentedMirror::BuildMirror(const char* name,
                                              AGeoAsphericDisk* disk,
                                              Bool_t isPrimary)
{
  Double_t zmin = disk->GetOrigin()[2] - disk->GetDZ();
  Double_t zmax = disk->GetOrigin()[2] + disk->GetDZ();
  Double_t dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // Define the base of the segment shape
  // ____
  // ___ ----____
  // C  ---- B    Rmax  
  //p2---p1/
  //    / /
  // --p0/
  // ----____
  //     A   ---- Rmin
  //
  Double_t d2r = TMath::DegToRad();
  Double_t ax = fRmin*TMath::Cos((90. - dphi)*d2r);
  Double_t ay = fRmin*TMath::Sin((90. - dphi)*d2r);
  Double_t bx = fRmax*TMath::Cos((90. - dphi)*d2r);
  Double_t by = fRmax*TMath::Sin((90. - dphi)*d2r);
  Double_t cy = fRmax;

  Double_t p0x = ax - fMargin*(1./TMath::Sin(dphi*d2r) - 1.)/TMath::Tan((90. - dphi)*d2r);
  Double_t p0y = ay + fMargin;

  Double_t t1 = TMath::Tan((90. - dphi)*d2r);
  Double_t t2 = TMath::Tan((90. - dphi/2.)*d2r);
  Double_t s1 = TMath::Sin(dphi*d2r);
  Double_t s2 = TMath::Sin((90. - dphi/2.)*d2r);
  Double_t p1x = 1/(t1 + 1/t2)*(cy - fMargin/s2 + t1*ax - ay - fMargin/s1);
  Double_t p1y = t1*(p1x - ax) + ay + fMargin/s1;
  Double_t p2x = 0;
  Double_t p2y = cy - fMargin/s2;

  TGeoXtru* xtru1 = new TGeoXtru(2);
  xtru1->SetName(Form("%s_xtru1", name));
  Double_t x[5] = {p0x, p1x, p2x, -p1x, -p0x};
  Double_t y[5] = {p0y, p1y, p2y, p1y, p0y};
  xtru1->DefinePolygon(5, x, y);
  xtru1->DefineSection(0, zmin);
  xtru1->DefineSection(1, zmax);

  TGeoTranslation* tr1 = new TGeoTranslation(Form("%s_tr1", name), 0, -cr, -cz);
  tr1->RegisterYourself();
  
  TGeoCompositeShape* comp
    = new TGeoCompositeShape(Form("%s_comp", name),
                             Form("%s:%s*%s:%s",
                                  disk->GetName(), tr1->GetName(),
                                  xtru1->GetName(), tr1->GetName()));
  
  AMirror* mirror = new AMirror(Form("%s_mirror", name), comp);
  
  return mirror;
}

//______________________________________________________________________________
class NewSCT
{
private:
  AOpticsManager*   fManager;
  AGeoAsphericDisk* fPrimaryV;
  AGeoAsphericDisk* fSecondaryV;

  // Optics Parameters
  Double_t fF;     // Focal length
  Double_t fAlpha; // \alpha
  Double_t fQ;     // q
  Double_t fRpMax; // Primary radius max
  Double_t fRpMin; // Primary radius min
  Double_t fRsMax; // Secondary radius max
  Double_t fRsMin; // Secondary radius min
  Double_t fKappa1;// Focal plane sag constant
  Double_t fKappa2;// Focal plane sag constant

  Int_t     fNp;   // Number of coefficients for the primary
  Int_t     fNs;   // Number of coefficients for the secondary
  Double_t* fP;   // Polynomial coefficients (p0, p1 ...)
  Double_t* fS;   // Polynomial coefficients (s0, s1 ...)

  // Focal plane
  Double_t fRf;

  // MAPMT Parameters
  Double_t fPixelSize;            // Width of a MAPMT pixel
  Double_t fPixelPitch;           // Pitch of MAPMT pixels
  Double_t fMAPMTWidth;           // Housing width
  Double_t fMAPMTLength;          // between input window and anode pins
  Double_t fInputWindowThickness; // Thickness of the input window

  // For TestPerformance()
  TCanvas* fCanvasSpot; // canvas for spot diagram
  TCanvas* fCanvasTime; // canvas for time gradient
  TCanvas* fCanvasEffArea; // canvas for effective area
  TCanvas* fCanvasSigma; // canvas for sigma vs incident angle

public:
  NewSCT();
  ~NewSCT();

  void AddIdealFocalPlane();
  void AddMAPMTFocalPlane();
  void AddPrimaryMirror(const char* name, SegmentedMirror* mirror);
  void AddSecondaryMirror(const char* name, SegmentedMirror* mirror);
  void AddSecondaryObscuration();
  void CloseGeometry() {fManager->CloseGeometry();}
  AOpticsManager* GetManager() const {return fManager;}
  void InitParameters();
  void TestPerformance() const;
};

//______________________________________________________________________________
NewSCT::NewSCT()
{
  InitParameters();

  // This method must be called after adding segmented mirrors
  fManager = new AOpticsManager("manager", "The optics manager of SCT");
  fManager->SetVisLevel(5);
  fManager->SetNsegments(50);

  // Make dummy material
  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(70); // needed in OpenGL view
  new TGeoMedium("med", 1, mat);

  // Make the world
  TGeoBBox* worldbox = new TGeoBBox("worldbox", 30*m, 30*m, 30*m);
  AOpticalComponent* world = new AOpticalComponent("world", worldbox);
  fManager->SetTopVolume(world);

  const Double_t kZp = 0.*m;
  const Double_t kZs = fF/fQ;

  // Make the ideal volume of the primary mirror
  fPrimaryV = new AGeoAsphericDisk("primaryV",
                                   kZp + fP[0] - 1*um, 0, kZp + fP[0] , 0, fRpMax, fRpMin);
  fPrimaryV->SetPolynomials(fNp - 1, &fP[1], fNp - 1, &fP[1]);

  // Make the ideal volume of the secondary mirror
  fSecondaryV = new AGeoAsphericDisk("secondaryV",
                                     kZs + fS[0], 0, kZs + fS[0]  + 1*um, 0, fRsMax, fRsMin);
  fSecondaryV->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);
}

//______________________________________________________________________________
NewSCT::~NewSCT()
{
  // Note that all the geometry components are automatically deleted by ROOT
  delete fManager;

  delete [] fP;
  delete [] fS;

  delete fCanvasSpot;
  fCanvasSpot = 0;

  delete fCanvasTime;
  fCanvasTime = 0;

  delete fCanvasEffArea;
  fCanvasEffArea = 0;

  delete fCanvasSigma;
  fCanvasSigma = 0;
}

//______________________________________________________________________________
void NewSCT::AddIdealFocalPlane()
{
  const Double_t kZs = fF/fQ;
  const Double_t kZf = kZs - (1 - fAlpha)*fF;

  AGeoAsphericDisk* idealCameraV = new AGeoAsphericDisk("idealCameraV", kZf - 1*um, 0, kZf, 0, fRf, 0);
  Double_t sagPar[2] = {fKappa1*TMath::Power(fF, -1),
                        fKappa2*TMath::Power(fF, -3)};
  idealCameraV->SetPolynomials(2, sagPar, 2, sagPar);
  AFocalSurface* idealCamera = new AFocalSurface("idealCamera", idealCameraV);
  AObscuration* idealCameraObs = new AObscuration("idealCameraObs", idealCameraV);
  fManager->GetTopVolume()->AddNode(idealCamera, 1);
  fManager->GetTopVolume()->AddNode(idealCameraObs, 1, new TGeoTranslation(0, 0, -1*um));
}

//______________________________________________________________________________
void NewSCT::AddMAPMTFocalPlane()
{
  // Make MAPMT photocathode without pixel structure
  TGeoBBox* mapmtCathodeV = new TGeoBBox("mapmtCathodeV", fPixelSize*4, fPixelSize*4, 100*um); // very thin box
  AFocalSurface* mapmtCathode = new AFocalSurface("mapmtCathode", mapmtCathodeV);

  // Make a single MAPMT
  TGeoBBox* mapmtV = new TGeoBBox("mapmtV", fMAPMTWidth/2., fMAPMTWidth/2.,
                                  fMAPMTLength/2.);
  AOpticalComponent* mapmt = new AOpticalComponent("mapmt", mapmtV);
  TGeoBBox* mapmtInputWindowV = new TGeoBBox("mapmtInputWindowV",
                                             fMAPMTWidth/2., fMAPMTWidth/2.,
                                             fInputWindowThickness/2.);

  TGeoMedium* med = fManager->GetMedium("med");
  ALens* mapmtInputWindow = new ALens("mapmtInputWindow", mapmtInputWindowV, med);
  ARefractiveIndex* bk7 = AGlassCatalog::GetRefractiveIndex("N-BK7");
  mapmtInputWindow->SetRefractiveIndex(bk7);
  mapmt->AddNodeOverlap(mapmtInputWindow, 1, new TGeoTranslation(0, 0, fMAPMTLength/2. - fInputWindowThickness/2.));

  mapmt->AddNode(mapmtCathode, 1, new TGeoTranslation(0, 0, fMAPMTLength/2. - fInputWindowThickness - 100*um));

  TGeoBBox* mapmtBackObsV = new TGeoBBox("mapmtBackObsV",
                                         fMAPMTWidth/2., fMAPMTWidth/2.,
                                         15*mm);
  AObscuration* mapmtBackObs = new AObscuration("mapmtBackObs", mapmtBackObsV);
  mapmt->AddNode(mapmtBackObs, 1, new TGeoTranslation(0, 0, -fMAPMTLength/2. + 15*mm));

  const Double_t kZs = fF/fQ;
  const Double_t kZf = kZs - (1 - fAlpha)*fF;

  // Make the focal plane
  Int_t n = 1;
  for(Int_t i = -7; i <= 7; i++){
    Double_t dx = i*fMAPMTWidth;
    for(Int_t j = -7; j <= 7; j++){
      if((TMath::Abs(i) + TMath::Abs(j) >= 11) || (TMath::Abs(i)*TMath::Abs(j) == 21)){
        continue;
      } // if
      Double_t dy = j*fMAPMTWidth;
      Double_t r2 = (i*i + j*j)*fMAPMTWidth*fMAPMTWidth;
      Double_t dz = fKappa1*TMath::Power(fF, -1)*r2 + fKappa2*TMath::Power(fF, -3)*r2*r2;
      fManager->GetTopVolume()->AddNode(mapmt, n, new TGeoTranslation(dx, dy, kZf - fMAPMTLength/2. + dz));
      n++;
    } // y
  } // x
}

//______________________________________________________________________________
void NewSCT::AddPrimaryMirror(const char* name, SegmentedMirror* mirror)
{
  AMirror* mir = mirror->BuildMirror(name, fPrimaryV, kTRUE);
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fPrimaryV, kTRUE);

  ABorderSurfaceCondition * condition
    = new ABorderSurfaceCondition(fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  fManager->GetTopVolume()->AddNode(mir, 1, combi);
}

//______________________________________________________________________________
void NewSCT::AddSecondaryMirror(const char* name, SegmentedMirror* mirror)
{
  AMirror* mir = mirror->BuildMirror(name, fSecondaryV, kFALSE);
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fSecondaryV, kFALSE);

  ABorderSurfaceCondition * condition
    = new ABorderSurfaceCondition(fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());

  fManager->GetTopVolume()->AddNode(mir, 1, combi);
}

//______________________________________________________________________________
void NewSCT::AddSecondaryObscuration()
{
  const Double_t kZs = fF/fQ;
  AGeoAsphericDisk* disk
    = new AGeoAsphericDisk("secondaryObsV", kZs + fS[0] + 1.*cm, 0, kZs + fS[0]  + 1*cm, 0, fRsMax, 0);
  disk->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  TGeoMedium* med = fManager->GetMedium("med");
  AObscuration* secondaryObs = new AObscuration("secondaryObs", disk, med);
  fManager->GetTopVolume()->AddNode(secondaryObs, 1);
}

//______________________________________________________________________________
void NewSCT::InitParameters()
{
  // SCT class has many parameters which should be given by setter methods as in
  // SegmentedMirror class. But I do not want to many setters here.

  // These values are retrieved from SCT-OPTMO/121108 ver. 2.00
  // "Optical System of 9.5m Schwarzschild-Couder Telescope for CTA"
  // Note that some parameter definitions have been changed since the time of
  // the first document. GrOptics config must be updated.
  fF      = 5.5863*m;
  fAlpha  = 2./3.;
  fQ      = 2./3.;
  fRpMax  = 4.8319*m;
  fRpMin  = 2.1933*m;
  fRsMax  = 2.7083*m;
  fRsMin  = 0.3945*m;
  fKappa1 = -0.8327;
  fKappa2 = 4.9950;

  fRf = 0.39*m;

  // These coefficients were calculated from polynomial fittings by A.O.
  // Note that the above documet uses different values analliticaly obtained
  // from Taylor expansions

  fNp = 8;
  fP = new Double_t[fNp];
  fP[0] = TMath::Power(fF,  1)* 8.4374e-06;
  fP[1] = TMath::Power(fF, -1)* 0.110917;
  fP[2] = TMath::Power(fF, -3)*-0.00511208;
  fP[3] = TMath::Power(fF, -5)*-0.0118961;
  fP[4] = TMath::Power(fF, -7)* 0.0253067;
  fP[5] = TMath::Power(fF, -9)*-0.0460152;
  fP[6] = TMath::Power(fF,-11)* 0.0413689;
  fP[7] = TMath::Power(fF,-13)*-0.0172745;

  fNs = 8;
  fS = new Double_t[fNs];
  fS[0] = TMath::Power(fF,  1)* 6.45608e-08;
  fS[1] = TMath::Power(fF, -1)*-0.416688;
  fS[2] = TMath::Power(fF, -3)*-0.144035;
  fS[3] = TMath::Power(fF, -5)* 0.647955;
  fS[4] = TMath::Power(fF, -7)*-2.96087;
  fS[5] = TMath::Power(fF, -9)* 9.39256;
  fS[6] = TMath::Power(fF,-11)*-18.0811;
  fS[7] = TMath::Power(fF,-13)* 15.3711;

  // MAPMT parameters 
  fPixelSize            = 6.08*mm;
  fPixelPitch           = 6.08*mm;
  fMAPMTWidth           = 52.0*mm;
  fMAPMTLength          = 32.7*mm;
  fInputWindowThickness = 1.5*mm;

  fCanvasSpot = 0;
  fCanvasTime = 0;
  fCanvasEffArea = 0;
  fCanvasSigma = 0;
}

//______________________________________________________________________________
void NewSCT::TestPerformance() const
{
  const Int_t kN = 42; // 0 to 4.1 [deg]
  const Double_t kDegStep = 0.1;
  TH2D* histSpot[kN];
  TH1D* histTime[kN];

  TGraph* graEffArea = new TGraph;
  graEffArea->SetTitle(";Field angle (deg);Effective Area (m^{2})");

  TGraph* graSigma[2];
  for(Int_t i = 0; i < 2; i++){
    graSigma[i] = new TGraph;
    graSigma[i]->SetTitle(";Field angle (deg);2 #times max{#sigma_{sagital}, #sigma_{tangential}} (arcmin)");
    graSigma[i]->SetLineStyle(i + 1);
  } // i

  TGraph* graTime = new TGraph;
  graTime->SetTitle(";Field angle (deg);Photon propagation time spread (#sigma) (ns)");

  for(Int_t i = 0; i < kN; i++){
    Double_t deg = i*kDegStep;

    TObject* obj;
    
    obj = gROOT->Get(Form("histSpot%d", i));
    if(obj){
      delete obj;
      obj = 0;
    } // if
      
    histSpot[i] = new TH2D(Form("histSpot%d", i), Form("#it{#theta} = %.1f (deg);X (arcmin);Y (arcmin)", deg), 1000, -10, 10, 1000, -10, 10);

    obj = gROOT->Get(Form("histTime%d", i));
    if(obj){
      delete obj;
      obj = 0;
    } // if
    
    histTime[i]= new TH1D(Form("histTime%d",i), Form("#it{#theta} = %.1f (deg);Propagation delay (ns);Entries", deg), 120, -6, 6);

    const Double_t kZs = fF/fQ;
    const Double_t kZf = kZs - (1 - fAlpha)*fF;

    TGeoTranslation raytr("raytr", -2*kZs*TMath::Sin(deg*TMath::DegToRad()), 0, 2*kZs*TMath::Cos(deg*TMath::DegToRad()));

    TVector3 dir;
    dir.SetMagThetaPhi(1, TMath::Pi() - deg*TMath::DegToRad(), 0);
    Double_t lambda = 400*nm;
    ARayArray* array = ARayShooter::Square(lambda, 2*fRpMax, 401, 0, &raytr, &dir);
    //    ARayArray* array = ARayShooter::Square(lambda, 2*fRpMax, 11, 0, &raytr, &dir);
    GetManager()->TraceNonSequential(*array);
    TObjArray* focused = array->GetFocused();

    Double_t Aeff = 0.;
    for(Int_t j = 0; j <= focused->GetLast(); j++){
      ARay* ray = (ARay*)(*focused)[j];
      if(!ray) continue;
      
      // Calculate the effective area from the number of focused photons
      Aeff += 2*fRpMax*2*fRpMax/400./400./m/m;
      //      Aeff += 2*fRpMax*2*fRpMax/10./10./m/m;
      
      Double_t p[4];
      ray->GetLastPoint(p);
      ray->SetLineWidth(1);
      if(deg == 0 && gRandom->Uniform(1) < 0.01){
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineColor(4);
        pol->SetLineStyle(2);
        pol->Draw();
      } // if

      Double_t deg2dist = 1.62499*mm*60.;
      Double_t x = deg*deg2dist;
      histSpot[i]->Fill((p[0] - x)/deg2dist*60, p[1]/deg2dist*60);
      histTime[i]->Fill((p[3] - (4*kZs - kZf)/(TMath::C()*m))/1e-9); // ns
    } // j

    graEffArea->SetPoint(graEffArea->GetN(), deg, Aeff);

    Double_t rmsx = histSpot[i]->GetRMS(1);
    Double_t rmsy = histSpot[i]->GetRMS(2);

    graSigma[0]->SetPoint(graSigma[0]->GetN(), deg, 2*rmsx);
    graSigma[1]->SetPoint(graSigma[1]->GetN(), deg, 2*rmsy);

    graTime->SetPoint(graTime->GetN(), deg, histTime[i]->GetRMS());

    delete array;
  } // n

  TCanvas* canSpot = new TCanvas("canSpot", "canSpot", 900, 900);
  canSpot->Divide(3, 3, 1e-10, 1e-10);

  //  TCanvas* canTime = ExactSizeCanvas("canTime", "canTime", 900, 900);
  TCanvas* canTime = new TCanvas("canTime", "canTime", 900, 900);
  canTime->Divide(3, 3, 1e-10, 1e-10);

  for(Int_t i = 0; i < kN; i += 5){
    canSpot->cd(i/5 + 1);
    histSpot[i]->Draw("colz");

    canTime->cd(i/5 + 1);
    histTime[i]->Draw();
  } // i

  // Figure 5 in the paper
  TCanvas* canFig5 = new TCanvas("canFig5", "canFig5", 1200, 600);
  canFig5->Divide(2, 1);
  canFig5->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graEffArea->Draw("apl");
  graEffArea->SetMarkerStyle(25);
  graEffArea->GetXaxis()->SetLimits(0, 5);
  graEffArea->GetYaxis()->SetRangeUser(0, 60);

  // PSF is not consistent with the original paper, but the spot diagram at
  // 5 (deg) is consistent with each other by eye comparison. There may be a
  // difference between calculations of RMS in my code and the paper
  canFig5->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  graSigma[0]->Draw("apl");
  graSigma[1]->Draw("pl same");
  graSigma[0]->SetMarkerStyle(25);
  graSigma[1]->SetMarkerStyle(20);
  graSigma[0]->GetXaxis()->SetLimits(0, 5);
  graSigma[0]->GetYaxis()->SetRangeUser(0, 10);

  // Figure 10 in the paper
  // Time spread is 2 times larger in the original paper. I believe the paper
  // is wrong. You can roughly calculate the spread width by
  // Dp * sin(angle)/c ~ 2.5 (ns)
  TCanvas* canFig10 = new TCanvas("canFig10", "canFig10", 1200, 600);
  canFig10->Divide(2, 1);
  canFig10->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graTime->Draw("apl");
  graTime->SetMarkerStyle(25);
  graTime->GetXaxis()->SetLimits(0, 5);
  graTime->GetYaxis()->SetRangeUser(0, 1.8);

  canFig10->cd(2);
  histTime[5]->Draw();

}

//______________________________________________________________________________
void AddPrimaryF(NewSCT* sct, Double_t margin = 3.5*mm)
{
  Int_t count = 0;

  // P1 mirrors
  for(Int_t i = 0; i < 16; i++){
    Double_t rmin = 2.19350*m - margin/TMath::Cos(11.25/2.*TMath::DegToRad());
    Double_t rmax = 3.40000*m;

    Double_t phimin = 22.5*i;
    Double_t phimax = 22.5*(i + 1);
    PentagonSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
    mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
    mirror.SetRotationErrors(0., 0., 0.);
    mirror.SetRougness(0.);
    mirror.SetMargin(margin);
    sct->AddPrimaryMirror(Form("primary%d", count), &mirror);
    count++;
  } // i
 
  // P2 mirrors
  for(Int_t i = 0; i < 32; i++){
    Double_t rmin = 3.40000*m;
    Double_t rmax = 4.831875*m + margin/TMath::Cos(11.25/2.*TMath::DegToRad());

    Double_t phimin = 11.25*i;
    Double_t phimax = 11.25*(i + 1);
    TetragonSegmentedMirror mirror(rmin, rmax + margin, phimin, phimax);
    mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
    mirror.SetRotationErrors(0., 0., 0.);
    mirror.SetRougness(0.);
    mirror.SetMargin(margin);
    sct->AddPrimaryMirror(Form("primary%d", count), &mirror);
    count++;
  } // i
}

//______________________________________________________________________________
void AddSecondaryJ(NewSCT* sct, Double_t margin = 3.5*mm)
{
  Int_t count = 0;

  // S1 mirrors
  for(Int_t i = 0; i < 8; i++){
    Double_t rmin = 0.3945*m;
    Double_t rmax = 1.5965*m;

    Double_t phimin = 45.*i;
    Double_t phimax = 45.*(i + 1);
    SectorSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
    mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
    mirror.SetRotationErrors(0., 0., 0.);
    mirror.SetRougness(0.);
    mirror.SetMargin(margin);
    sct->AddSecondaryMirror(Form("secondary%d", count), &mirror);
    count++;
  } // i

  // S2 mirrors
  for(Int_t i = 0; i < 16; i++){
    Double_t rmin = 1.5965*m;
    Double_t rmax = 2.7083*m;

    Double_t phimin = 22.5*i;
    Double_t phimax = 22.5*(i + 1);
    SectorSegmentedMirror mirror(rmin, rmax + margin, phimin, phimax);
    mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
    mirror.SetRotationErrors(0., 0., 0.);
    mirror.SetRougness(0.);
    mirror.SetMargin(margin);
    sct->AddSecondaryMirror(Form("secondary%d", count), &mirror);
    count++;
  } // i
}

//______________________________________________________________________________
void sct_test()
{
  NewSCT* sct = new NewSCT;

  sct->GetManager()->DisableFresnelReflection(1);

  Double_t margin = 0.001*mm;
  //  Double_t margin = 3.5*mm;

  AddPrimaryF(sct, margin);
  AddSecondaryJ(sct, margin);

  sct->AddSecondaryObscuration();
  //sct->AddIdealFocalPlane();
  sct->AddMAPMTFocalPlane();

  sct->CloseGeometry();

  // under development/b
  // sct->GetManager()->SetMultiThread(kTRUE);
  // sct->GetManager()->SetMaxThreads(8);

  sct->TestPerformance();

  sct->GetManager()->GetTopVolume()->Draw("ogl");
}

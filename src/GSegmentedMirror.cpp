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
#include "GSegmentedMirror.h"

#include <iostream>

// define useful units
static const Double_t cm = AOpticsManager::cm();
static const Double_t mm = AOpticsManager::mm();
static const Double_t um = AOpticsManager::um();
static const Double_t nm = AOpticsManager::nm();
static const Double_t  m = AOpticsManager::m();

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

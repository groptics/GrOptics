/*
VERSION3.1
2March2015
*/
/*  GSegSCTelescope.cpp

    Akira Okumura 

        and

    Charlie Duke
    Grinnell College
    May 2011
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>

using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TRandom3.h"
#include "TPolyLine3D.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
//#include "TGeoRotation.h"
#include "TGeoMatrix.h"
//#include "TGeoCombiTrans.h"
#include "AGeoAsphericDisk.h"
#include "AGlassCatalog.h"
#include "TView.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVector3.h"
#include "TGeoPgon.h"
#include "TVersionCheck.h"
#include "TGraph.h"
#include "TBrowser.h"

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

#include "GUtilityFuncts.h"
#include "GDefinition.h"

#include "GSegmentedMirror.h"
#include "GSegmentedObscuration.h"
#include "GTelescope.h"
#include "GSegSCTelescope.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGW(x) *oLog << "         " << #x << " = " << x << endl

// define useful units
static const double cm = AOpticsManager::cm();
static const double mm = AOpticsManager::mm();
static const double um = AOpticsManager::um();
static const double nm = AOpticsManager::nm();
static const double  m = AOpticsManager::m();

GSegSCTelescope::GSegSCTelescope() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::GSegSCTelescope() " << endl;
  }
  initialize();
  };
/********************** end of GSegSCTelescope *****************/
GSegSCTelescope::~GSegSCTelescope() {

  bool debug=false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::~GSegSCTelescope " << endl;
  }
  if (fManager !=0) {
    gGeoManager = fManager;
    SafeDelete(fManager);
  }

  if (hisF != 0) SafeDelete(hisF);
 
  map<int, TGraph *>::iterator itmGRefl; 
  for (itmGRefl=mGRefl->begin();
       itmGRefl!=mGRefl->end(); itmGRefl++) {
    SafeDelete(itmGRefl->second ); 
  }
  SafeDelete(mGRefl);
  SafeDelete(ray);

  if (fS != 0) delete[] fS;
  if (fP != 0) delete[] fP;
};
/********************** end of ~GSegSCTelescope *****************/

void GSegSCTelescope::buildTelescope(bool os8)
{
  // fix units
  fFMeters = fF;
  fF = fF*m;
  fZp = fZp*m;

  fTelRadius = fRpMax;
  fRpMax = fRpMax*m;
  fRpMin = fRpMin*m;
  fRsMax = fRsMax*m;
  fRsMin = fRsMin*m;
  //fTX = fTX*m;
  //fTY = fTY*m;
  //fTZ = fTZ*m;

  fPixelSize = fPixelSize*mm;
  fMAPMTWidth = fMAPMTWidth*mm;
  fMAPMTLength = fMAPMTLength*mm;
  fInputWindowThickness = fInputWindowThickness*mm;
  fMAPMTOffset = fMAPMTOffset*mm;
  fMAPMTGap = fMAPMTGap*mm;

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::buildTelescope" << endl;
  }
  gGeoManager = 0;
  fManager = new AOpticsManager("manager","The optics manager of SEGSC");
  //fManager->SetVisLevel(5);// should be 0 or 1
  //fManager->SetNsegments(50);
  fManager->DisableFresnelReflection(1);

  // Make dummy material
  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(70); // needed in OpenGL view, 70
  new TGeoMedium("med", 1, mat);

  // Make the world
  TGeoBBox* worldbox = new TGeoBBox("worldbox", fTX*m, fTY*m, fTZ*m);
  AOpticalComponent* world = new AOpticalComponent("world", worldbox);
  fManager->SetTopVolume(world);
  fManager->SetTopVisible(0); 
  //fManager->SetTopVisible(.9); 

  makePrimarySecondaryDisks();
  addPrimaryF();
  addSecondaryJ();
  addSecondaryObscuration();

  if (bCameraFlag) {
    addMAPMTFocalPlane();
  }
  else {
    addIdealFocalPlane();
  }

  closeGeometry();

  *oLog << "ssssssssssssssssssiPrtMode " << iPrtMode << endl;
  printTelescope();
  //testPerformance();

  return;
};
/*************************************************************************************/
void GSegSCTelescope::makePrimarySecondaryDisks() {

  bool debug = true;

  if (debug) {
    *oLog << "  --  GSegSCTelescope::makePrimarySecondaryDisks" << endl;
  }

  const Double_t kZp = fZp;
  const Double_t kZs = (fF)*fZs;
  
  if (0) {
    *oLog << "   fF " << fF << endl;
    *oLog << " primary" << endl;
    *oLog << "   kZp+fP[0]      " << kZp + fP[0] << endl;
    *oLog << "   kZp+fP[0]-1*um " << kZp + fP[0] - 1*um << endl;
    *oLog << "   fRpMax, fRpMin " << fRpMax << " " << fRpMin << endl;
    *oLog << " secondary" << endl;
    *oLog << "   kZs + fS[0]    " << kZs + fS[0] << endl;
    *oLog << "   kZs + fS[0]  + 1*um " << kZs + fS[0]  + 1*um << endl;
    *oLog << "   fRsMax, fRsMin " << fRsMax << "  " << fRsMin << endl;
  }
  fPrimaryV = new AGeoAsphericDisk("primaryV",
                                   kZp + fP[0] - 1*um, 0,
                                   kZp + fP[0] , 0, 
                                    fRpMax, fRpMin);

  if (0) {
    *oLog << "   primary coefficients " << endl;
    *oLog << "   fNp " << fNp << endl;
    for (int i = 0;i<fNp; i++ ) {
      *oLog << "   i  fp[i] " << i << "   " << fP[i] << endl;
      *oLog << "     fzp[i]    " << fzp[i] << endl;
    }
  }
  fPrimaryV->SetPolynomials(fNp - 1, &fP[1], fNp - 1, &fP[1]);

  // Make the ideal volume of the secondary mirror
  if (0) {
    *oLog << "   secondary coefficients " << endl;
    *oLog << "   fNs " << fNs << endl;
    for (int i = 0;i<fNs; i++ ) {
      *oLog << "   i  fs " << i << "   " << fS[i] << endl;
      *oLog << "      fzs[i] " << fzs[i] << endl;
    }
  }
  fSecondaryV = new AGeoAsphericDisk("secondaryV",
                                     kZs + fS[0] , 0, 
                                     kZs + fS[0]  + 1*um , 
                                     0, fRsMax, fRsMin);
  fSecondaryV->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  // make volume for secondary obscurations
  // not currently used: needed if make obscurations for individual secondary segments
  fSecondaryObsV = new AGeoAsphericDisk("secondaryObsV",
                                     kZs + fS[0] + 5.*cm, 0, 
                                     kZs + fS[0] + 5.*cm + 1*um, 
                                     0, fRsMax, fRsMin);

  fSecondaryObsV->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);
    
  return;
};
/*************************************************************************************/

void GSegSCTelescope::addPrimaryF() {
  bool debug = true;

  if (debug) {
    *oLog << "  --  GSegSCTelescope::addPrimaryF" << endl;
   }
  Int_t count = 0;

  // segment not created if reflect = 0
  // P1 mirrors
  for (Int_t i = 0; i < iNumP1Mirrors; i++) {
    if ( ( (*(vSegP1.at(i))).reflect) > 0) {
      Double_t rmin = (*(vSegP1.at(i))).rmin;
      Double_t rmax = (*(vSegP1.at(i))).rmax;
      Double_t margin = ( (*(vSegP1.at(i))).margin )*mm;
      rmin = rmin*m - margin/TMath::Cos(11.25/2.*TMath::DegToRad());
      rmax = rmax*m;
      
      Double_t phimin = ( (*(vSegP1.at(i))).delPhi)*i;
      Double_t phimax = ( (*(vSegP1.at(i))).delPhi)*(i+1);
      if (0) {
        *oLog << "    P1 pentagon segmented mirror number " << i << endl;
        *oLog << "        rmin/rmax " << rmin << " " << rmax << endl;
        *oLog << "        margin    " << margin << endl;
        *oLog << "        phimin/max " << phimin << "  " << phimax << endl;
      }
    PentagonSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
    Double_t posErrorX = (*(vSegP1.at(i))).posErrorX;
    Double_t posErrorY = (*(vSegP1.at(i))).posErrorY;
    Double_t posErrorZ = (*(vSegP1.at(i))).posErrorZ;
    mirror.SetPositionErrors(posErrorX*mm, posErrorY*mm, posErrorZ*mm);
 
    Double_t rotErrorPhi = (*(vSegP1.at(i))).rotErrorPhi;
    Double_t rotErrorTheta = (*(vSegP1.at(i))).rotErrorTheta;
    Double_t rotErrorPsi = (*(vSegP1.at(i))).rotErrorPsi;
    mirror.SetRotationErrors(rotErrorPhi, rotErrorTheta,
                             rotErrorPsi);
 
    Double_t roughness = (*(vSegP1.at(i))).roughness; 
    mirror.SetRougness(roughness);

    mirror.SetMargin(margin);
    iReflect = (*(vSegP1.at(i))).reflect;
   // add mirror segment
    addPrimaryMirror(Form("primary%d", count), &mirror);
    count++;
    } 
  }

  // P2 mirrors
  for (Int_t i = 0; i < iNumP2Mirrors; i++) {
    if ( ( (*(vSegP2.at(i))).reflect) > 0) {
      Double_t rmin = (*(vSegP2.at(i))).rmin;
      Double_t rmax = (*(vSegP2.at(i))).rmax;
      Double_t margin = ( (*(vSegP2.at(i))).margin )*mm;
      rmin = rmin*m - margin/TMath::Cos(11.25/2.*TMath::DegToRad());
      rmax = rmax*m;
      
      Double_t phimin = ( (*(vSegP2.at(i))).delPhi)*i;
      Double_t phimax = ( (*(vSegP2.at(i))).delPhi)*(i+1);
      if (0) {
        *oLog << "    P1 pentagon segmented P2 mirror number " << i << endl;
        *oLog << "        rmin/rmax " << rmin << " " << rmax << endl;
        *oLog << "        margin    " << margin << endl;
        *oLog << "        phimin/max " << phimin << "  " << phimax << endl;
      }
    PentagonSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
    Double_t posErrorX = (*(vSegP2.at(i))).posErrorX;
    Double_t posErrorY = (*(vSegP2.at(i))).posErrorY;
    Double_t posErrorZ = (*(vSegP2.at(i))).posErrorZ;
    mirror.SetPositionErrors(posErrorX*mm, posErrorY*mm, posErrorZ*mm);
 
    Double_t rotErrorPhi = (*(vSegP2.at(i))).rotErrorPhi;
    Double_t rotErrorTheta = (*(vSegP2.at(i))).rotErrorTheta;
    Double_t rotErrorPsi = (*(vSegP2.at(i))).rotErrorPsi;
    mirror.SetRotationErrors(rotErrorPhi, rotErrorTheta,
                             rotErrorPsi);
 
    Double_t roughness = (*(vSegP2.at(i))).roughness; 
    mirror.SetRougness(roughness);

    mirror.SetMargin(margin);
    iReflect = (*(vSegP2.at(i))).reflect;
    // add mirror segment
    addPrimaryMirror(Form("primary%d", count), &mirror);
    count++;
    } 
  }
};
/*******************************************************************/
void GSegSCTelescope::addPrimaryMirror(const char*name,
                                       SegmentedMirror *mirror) {

  gGeoManager = fManager;
  
  bool debug = false;
  if (debug) {
    *oLog << "  --  GSegSCTelescope::addPrimaryMirror" << endl;
  }
  
  AMirror* mir = mirror->BuildMirror(name, fPrimaryV, kTRUE);
  mir->SetLineColor(iPrimaryColor);
  // get and add TGraph for reflectivity mir->SetReflectivity(TGraph *)

  TGraph * graph = makeReflectivityGraph(iReflect);
  mir->SetReflectivity(graph); // graph owned by AMirror (and deleted)
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fPrimaryV, kTRUE);

  ABorderSurfaceCondition * condition
    = new ABorderSurfaceCondition(fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  fManager->GetTopVolume()->AddNode(mir, 1, combi);
  
};
/*******************************************************************/

void GSegSCTelescope::addSecondaryJ() {

  bool debug = true;
  
  if (debug) {
    *oLog << "  --  GSegSCTelescope::addSecondaryJ" << endl;
  }
  Int_t count = 0;
  Int_t count1 = 0;

  // S1 mirrors
  for (Int_t i = 0; i < iNumS1Mirrors; i++) {
    if ( ( (*(vSegS1.at(i))).reflect) > 0) {
      Double_t rmin = (*(vSegS1.at(i))).rmin;
      Double_t rmax = (*(vSegS1.at(i))).rmax;
      Double_t margin = ( (*(vSegS1.at(i))).margin )*mm;
      rmin = rmin*m - margin/TMath::Cos(11.25/2.*TMath::DegToRad());
      rmax = rmax*m;
      
      Double_t phimin = ( (*(vSegS1.at(i))).delPhi)*i;
      Double_t phimax = ( (*(vSegS1.at(i))).delPhi)*(i+1);
      SectorSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
      Double_t posErrorX = (*(vSegS1.at(i))).posErrorX;
      Double_t posErrorY = (*(vSegS1.at(i))).posErrorY;
      Double_t posErrorZ = (*(vSegS1.at(i))).posErrorZ;
      mirror.SetPositionErrors(posErrorX*mm, posErrorY*mm, posErrorZ*mm);
      
      Double_t rotErrorPhi = (*(vSegS1.at(i))).rotErrorPhi;
      Double_t rotErrorTheta = (*(vSegS1.at(i))).rotErrorTheta;
      Double_t rotErrorPsi = (*(vSegS1.at(i))).rotErrorPsi;
      mirror.SetRotationErrors(rotErrorPhi, rotErrorTheta,
                               rotErrorPsi);
      
      Double_t roughness = (*(vSegS1.at(i))).roughness; 
      mirror.SetRougness(roughness);
      
      mirror.SetMargin(margin);
      iReflect = (*(vSegS1.at(i))).reflect;
      // add mirror segment
      addSecondaryMirror(Form("secondary%d", count1), &mirror);    
      
      //SectorSegmentedObscuration  obscuration(rmin + margin, rmax,phimin,phimax);
      //obscuration.SetPositionErrors(posErrorX*mm, posErrorY*mm, 
      //                            posErrorZ*mm);
      //obscuration.SetRotationErrors(rotErrorPhi, rotErrorTheta,
      //                               rotErrorPsi);
      //obscuration.SetMargin(margin);
      //addSecondaryObscurationSeg(Form("secondaryObs%d", count1), &obscuration);    
      count++;
      count1++;
    }
  }

  // S2 mirrors
  for (Int_t i = 0; i < iNumS2Mirrors; i++) {
    if ( ( (*(vSegS2.at(i))).reflect) > 0) {
      Double_t rmin = (*(vSegS2.at(i))).rmin;
      Double_t rmax = (*(vSegS2.at(i))).rmax;
      Double_t margin = ( (*(vSegS2.at(i))).margin )*mm;
      rmin = rmin*m - margin/TMath::Cos(11.25/2.*TMath::DegToRad());
      rmax = rmax*m;
      
      Double_t phimin = ( (*(vSegS2.at(i))).delPhi)*i;
      Double_t phimax = ( (*(vSegS2.at(i))).delPhi)*(i+1);
      SectorSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
      Double_t posErrorX = (*(vSegS2.at(i))).posErrorX;
      Double_t posErrorY = (*(vSegS2.at(i))).posErrorY;
      Double_t posErrorZ = (*(vSegS2.at(i))).posErrorZ;
      mirror.SetPositionErrors(posErrorX*mm, posErrorY*mm, posErrorZ*mm);
      
      Double_t rotErrorPhi = (*(vSegS2.at(i))).rotErrorPhi;
      Double_t rotErrorTheta = (*(vSegS2.at(i))).rotErrorTheta;
      Double_t rotErrorPsi = (*(vSegS2.at(i))).rotErrorPsi;
      mirror.SetRotationErrors(rotErrorPhi, rotErrorTheta,
                               rotErrorPsi);
      
      Double_t roughness = (*(vSegS2.at(i))).roughness; 
      mirror.SetRougness(roughness);
      
      mirror.SetMargin(margin);
      iReflect = (*(vSegS2.at(i))).reflect;
      // add mirror segment
      addSecondaryMirror(Form("secondary%d", count), &mirror);  
      
      //SectorSegmentedObscuration  obscuration(rmin + margin, rmax, phimin, phimax);
      //obscuration.SetPositionErrors(posErrorX*mm, posErrorY*mm, posErrorZ*mm);
      //obscuration.SetRotationErrors(rotErrorPhi, rotErrorTheta,
      //                            rotErrorPsi);
      //obscuration.SetMargin(margin);
      //addSecondaryObscurationSeg(Form("secondaryObs%d", count1), &obscuration);    
      count++;
      count1++;
      
    }
  }

};
/*******************************************************************/

void GSegSCTelescope::addSecondaryObscuration() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::AddSecondaryObscuration " << endl;
  }

  //const Double_t kZs = fF/fQ;
  const Double_t kZs = fF*fZs;
  AGeoAsphericDisk* disk
    = new AGeoAsphericDisk("secondaryObsV", 
                           kZs + fS[0] + 1.*cm, 0, 
                           kZs + fS[0]  + 1.*cm + 1.*um, 0, fRsMax, 0);
  disk->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  TGeoMedium* med = fManager->GetMedium("med");
  AObscuration* secondaryObs = new AObscuration("secondaryObs", disk, med);
  secondaryObs->SetLineColor(iSecondaryObscurationColor);
  fManager->GetTopVolume()->AddNode(secondaryObs, 1);

};
/*******************************************************************/
void GSegSCTelescope::addSecondaryObscurationSeg(const char*name, 
                               SegmentedObscuration *obscuration) {

  gGeoManager = fManager;
  
  bool debug = false;
  if (debug) {
    *oLog << "  --  GSegSCTelescope::addSecondaryObscurationSeg" << endl;
  }

  AObscuration* obs = obscuration->BuildObscuration(name, fSecondaryObsV, kFALSE);
  obs->SetLineColor(iSecondaryObscurationColor);
  TGeoCombiTrans* combi = obscuration->BuildObscurationCombiTrans(fSecondaryObsV, kFALSE);


  fManager->GetTopVolume()->AddNode(obs, 1, combi);
  
};

/*******************************************************************/
void GSegSCTelescope::addSecondaryMirror(const char*name, SegmentedMirror *mirror) {
  gGeoManager = fManager;
  
  bool debug = false;
  if (debug) {
    *oLog << "  --  GSegSCTelescope::addSecondaryMirror" << endl;
  }

  AMirror* mir = mirror->BuildMirror(name, fSecondaryV, kFALSE);
  mir->SetLineColor(iSecondaryColor);
  TGraph * graph = makeReflectivityGraph(iReflect);
  mir->SetReflectivity(graph);

  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fSecondaryV, kFALSE);

  ABorderSurfaceCondition * condition
    = new ABorderSurfaceCondition(fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());

  fManager->GetTopVolume()->AddNode(mir, 1, combi);
  
};
/*************************************************************************************/
void GSegSCTelescope::addIdealFocalPlane()  {
  bool debug = true;

  if (debug) {
    *oLog << "  --  GSegSCTelescope::addIdealFocalPlane" << endl;
    *oLog << "        fF fQ fAlpha " << fF << " " << fQ << " " << fAlpha
          << endl;
  }

  //const Double_t kZs = (fF)/fQ;
  //const Double_t kZf = kZs - (1 - fAlpha)*fF;
  const Double_t kZf = fF * fZf;

  //Double_t focalPlaneHalfThickness = 1*um;
  Double_t focalPlaneHalfThickness = 10*cm;
  AGeoAsphericDisk* idealCameraV = new AGeoAsphericDisk("idealCameraV", kZf - 1*um, 0, kZf, 0, fRf*m, 0);
  Double_t sagPar[2] = {fKappa1*TMath::Power(fF, -1),
                        fKappa2*TMath::Power(fF, -3)};
  idealCameraV->SetPolynomials(2, sagPar, 2, sagPar);
  AFocalSurface* idealCamera = new AFocalSurface("idealCamera", idealCameraV);
  idealCamera->SetLineColor(iMAPMTCathodeColor);

  AObscuration* idealCameraObs = new AObscuration("idealCameraObs", idealCameraV);
  idealCameraObs->SetLineColor(iMAPMTObscurationColor);
  fManager->GetTopVolume()->AddNode(idealCamera, 1);

  Double_t obscurationOffset = -30*cm;
  //Double_t obscurationOffset = -100.*um;

  fManager->GetTopVolume()->AddNode(idealCameraObs, 1, new TGeoTranslation(0, 0, obscurationOffset));

};
/*************************************************************************************/

void GSegSCTelescope::addMAPMTFocalPlane()  {
  bool debug = true;

  if (debug) {
    *oLog << "  --  GSegSCTelescope::addMAPMTFocalPlane" << endl;
  }
  Double_t fWidthBox = 50.0*cm;
  Double_t fHeightBox = 10.*cm;
  TGeoMedium* med = fManager->GetMedium("med");

 // make a new volume for the camera
  // size adequately covers os8 camera/focal surface 
  TGeoVolume *focVol = gGeoManager->MakeBox("focVol",med,fWidthBox,
                                            fWidthBox,fHeightBox);
  ////////////////////////////////////////////////////////////////////////
  // Make MAPMT photocathode without pixel structure 
  Double_t cathodeHalfThick = 100*um;
  //Double_t cathodeHalfThick = 2.0*mm;
  TGeoBBox* mapmtCathodeV = new TGeoBBox("mapmtCathodeV", fPixelSize*4, 
                                         fPixelSize*4, cathodeHalfThick); // very thin box
  AFocalSurface* mapmtCathode = new AFocalSurface("mapmtCathode", mapmtCathodeV);
  mapmtCathode->SetLineColor(iMAPMTCathodeColor);
  if (debug) *oLog << "cathodeHalfThick " << cathodeHalfThick << endl;

  //////////////////////////////////////////////////////////////////////
  // Make a single MAPMT
  TGeoBBox* mapmtV = new TGeoBBox("mapmtV", fMAPMTWidth/2., fMAPMTWidth/2.,
                                  fMAPMTLength/2.);
  AOpticalComponent* mapmt = new AOpticalComponent("mapmt", mapmtV);

  ///////////////////////////////////////////////////////////
  // make input window
  TGeoBBox* mapmtInputWindowV = new TGeoBBox("mapmtInputWindowV",
                                             fMAPMTWidth/2., fMAPMTWidth/2.,
                                             fInputWindowThickness/2.);
  if (debug) *oLog << " fInputWindowThickness/2. " << fInputWindowThickness/2. << endl;

  ALens* mapmtInputWindow = new ALens("mapmtInputWindow", mapmtInputWindowV, med);
  mapmtInputWindow->SetLineColor(iMAPMTWindowColor);
  ARefractiveIndex* bk7 = AGlassCatalog::GetRefractiveIndex("N-BK7");
  mapmtInputWindow->SetRefractiveIndex(bk7);
  mapmt->AddNodeOverlap(mapmtInputWindow, 
                        1, new TGeoTranslation(0, 0, fMAPMTLength/2. - fInputWindowThickness/2.));

  if (debug) *oLog << "fMAPMTLength/2. - fInputWindowThickness/2. " 
                   << fMAPMTLength/2. - fInputWindowThickness/2. << endl;

  Double_t fWindowBottomRelToMapmtCenter = 
    (fMAPMTLength/2. - fInputWindowThickness/2.) - fInputWindowThickness/2.; // rel. to mapmt center

  Double_t cathodePosition = fMAPMTLength/2. - fInputWindowThickness - fMAPMTGap - cathodeHalfThick;
  mapmt->AddNode(mapmtCathode, 1, new TGeoTranslation(0, 0, cathodePosition));
  if (debug) *oLog << "cathodePosition " << cathodePosition << endl;

  Double_t fCathodeTopRelToMapmtCenter = cathodePosition + cathodeHalfThick;
  
  Double_t backObsThickness = 1*mm;
  TGeoBBox* mapmtBackObsV = new TGeoBBox("mapmtBackObsV",
                                         fMAPMTWidth/2., fMAPMTWidth/2.,
                                         backObsThickness);
  
  AObscuration* mapmtBackObs = new AObscuration("mapmtBackObs", mapmtBackObsV);
  mapmtBackObs->SetLineColor(iMAPMTObscurationColor);
  Double_t backObsPosition = -fMAPMTLength/2. + backObsThickness;
  mapmt->AddNode(mapmtBackObs, 1, new TGeoTranslation(0, 0,backObsPosition ));
  Double_t backObsTopPositionRelToMapmtCenter = backObsPosition + backObsThickness;
  
  *oLog << " fWindowBottomRelToMapmtCenter  " << fWindowBottomRelToMapmtCenter << endl;
  *oLog << " fCathodeTopRelToMapmtCenter  " << fCathodeTopRelToMapmtCenter << endl;
  *oLog << " backObsTopPositionRelToMapmtCenter  " << backObsTopPositionRelToMapmtCenter << endl;

  fCathodeTopRelToFocalSurface =  fCathodeTopRelToMapmtCenter + fMAPMTOffset -
    fCathodeTopRelToMapmtCenter;
  fWindowBottomRelToFocalSurface = fWindowBottomRelToMapmtCenter + fMAPMTOffset -
    fCathodeTopRelToMapmtCenter;
  fMAPOscurationTopRelToFocalSurface = backObsTopPositionRelToMapmtCenter + fMAPMTOffset -
    fCathodeTopRelToMapmtCenter;
  
  // sanity check
  Double_t fCathodeBottomRelToFocalSurface = fCathodeTopRelToFocalSurface - cathodeHalfThick*2.0; 
  fCathodeBottomRelToOscurationTop = fCathodeBottomRelToFocalSurface - fMAPOscurationTopRelToFocalSurface;

  if (fCathodeBottomRelToOscurationTop < 0.0) {
    *oLog << " fCathodeBottomRelToOscurationTop, " << fCathodeBottomRelToOscurationTop 
          << ",  is less than zero. Cathode is below the obscuration. Need to increase MAPMT length"
          << "    stopping code " << endl;
    exit(0);
  }
  
  //const Double_t kZs = fF/fQ;
  //const Double_t kZf = kZs - (1 - fAlpha)*fF;
  const Double_t kZf = fF * fZf;
  //*oLog << "  =================== kZf " << kZf << endl;
  // Make the focal plane
  Double_t mapmtPositionReltoFocalSurface = - fMAPMTLength/2. + fInputWindowThickness + fMAPMTGap;
  //*oLog << " xxxxxxxxxxx  mapmtPositionReltoFocalSurface " << mapmtPositionReltoFocalSurface << endl;
  Int_t n = 1;
  // loop from -iNum to +iNum
  Int_t iNum = 7; // set to 1, make gl plot and use CheckPoint to see location of center module
  if (bSingleMAPMTmodule == false) {
    for(Int_t i = -iNum; i <= iNum; i++){
      Double_t dx = i*fMAPMTWidth;
      for(Int_t j = -iNum; j <= iNum; j++){
        if((TMath::Abs(i) + TMath::Abs(j) >= 11) || (TMath::Abs(i)*TMath::Abs(j) == 21)){
          continue;
        } // if
        Double_t dy = j*fMAPMTWidth;
        Double_t r2 = (i*i + j*j)*fMAPMTWidth*fMAPMTWidth;
        Double_t dz = fKappa1*TMath::Power(fF, -1)*r2 + fKappa2*TMath::Power(fF, -3)*r2*r2;
        focVol->AddNode(mapmt, n, new TGeoTranslation(dx, dy, 
                                                      mapmtPositionReltoFocalSurface +
                                                      + fMAPMTOffset + dz));
        n++;
      } // y
    } // x
  }
  else {
    Double_t dx = 0.0;
    Double_t dy = 0.0;
    Double_t dz = 0.0;
    focVol->AddNode(mapmt, 1, new TGeoTranslation(dx, dy, 
                                                  mapmtPositionReltoFocalSurface +
                                                  + fMAPMTOffset + dz));

  }
  fManager->GetTopVolume()->AddNode(focVol,1,new TGeoCombiTrans("cFocS",
                                                               0.0,
                                                               0.0,
                                                               kZf,
                                                               new TGeoRotation("rFocS",
                                                                                0.0,
                                                                                0.0,
                                                                                0.0)));
                                                               
};
/*************************************************************************************/

void GSegSCTelescope::closeGeometry()  {
  gGeoManager = fManager;
  fManager->CloseGeometry();
};
/*************************************************************************************/

void GSegSCTelescope::testFocalPlane() {

  // inject photons 

};
/****************************************************************************************/

void GSegSCTelescope::injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                                const ROOT::Math::XYZVector &photonDirT,
				const double &photWaveLgt) {
  gGeoManager = fManager;

  bool debug = false;
  if (debug) {
    *oLog << " -- GSegSCTelescope::injectPhoton " << endl;
    *oLog << "       photonLocT  ";
    GUtilityFuncts::printGenVector(photonLocT); *oLog << endl;
    *oLog << "       photonDirT  ";
    GUtilityFuncts::printGenVector(photonDirT); *oLog << endl;
    *oLog << "      bPhotonHistoryFlag " << bPhotonHistoryFlag << endl;
  } 
  // for debugging and testing only

  /*
  fphotonInjectLoc[0] = 3.0;
  fphotonInjectLoc[1] = 0.0;
  fphotonInjectLoc[2] = 0.0;

  fphotonInjectDir[0] = 0.0;
  fphotonInjectDir[1] = 0.0;
  fphotonInjectDir[2] = -1.0;
  */
  //*oLog << "        TESTING WITH THESE VALUES" << endl;
  //*oLog << "             specified location and direction " << endl;
  //for (int i = 0;i<3;i++) {
  //*oLog << i << "  " << fphotonInjectLoc[i] << "  " 
  //	  << fphotonInjectDir[i] << endl;
  //}

  photonLocT.GetCoordinates(fInitialInjectLoc);
  photonLocT.GetCoordinates(fphotonInjectLoc);
  photonDirT.GetCoordinates(fphotonInjectDir); 

  // convert to cm as required for robast
  fphotWaveLgt = photWaveLgt*nm;

  // move base of injection location vector to center of primary mirror
  // (add rotation offset)
  fphotonInjectLoc[2] = fphotonInjectLoc[2] + fRotationOffset;

  movePositionToTopOfTopVol();
 
  // initialize photonHistory parameters if necessary
  if (bPhotonHistoryFlag) {
    initializePhotonHistoryParms();    
  }
  
  // Assuming that three arguments are given in units of (m), (m), (nm)
  double t = 0;
  double x  = fphotonInjectLoc[0];
  double y  = fphotonInjectLoc[1];
  double z  = fphotonInjectLoc[2];
  double dx = fphotonInjectDir[0];
  double dy = fphotonInjectDir[1];
  double dz = fphotonInjectDir[2];

  SafeDelete(ray);
  ray = new ARay(0, fphotWaveLgt, x*m, y*m, z*m, t, dx, dy, dz);

  gGeoManager = fManager;

  fManager->TraceNonSequential(*ray);
 
  // Here you can get the traced result
  double new_xyzt[4];
  ray->GetLastPoint(new_xyzt);
  fNPoints = ray->GetNpoints();
  double new_x = new_xyzt[0];
  double new_y = new_xyzt[1];
  double new_z = new_xyzt[2];
  double new_t = new_xyzt[3];
  if (debug) {
    *oLog << " from GetLastPoint " << new_x << " " << new_y 
	  << " " << new_z << " " << new_t << endl;
  }

};
/********************** end of injectPhoton *****************/
void GSegSCTelescope::movePositionToTopOfTopVol() {

  gGeoManager = fManager;

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::movePositionToTopOfTopVol " << endl;
    *oLog << "        position prior to move to top ";
    *oLog << fphotonInjectLoc[0] << "  " << fphotonInjectLoc[1] << "  " 
	  << fphotonInjectLoc[2] << endl;
  }

  Double_t rfx = fphotonInjectLoc[0];
  Double_t rfy = fphotonInjectLoc[1];
  Double_t rfz = fphotonInjectLoc[2];

  //Double_t Z = fTZ; // just inside top volume
  Double_t Z = fTZ; // just inside top volume

  Double_t dl = fphotonInjectDir[0];
  Double_t dm = fphotonInjectDir[1];
  Double_t dn = fphotonInjectDir[2];

  Double_t Rx = rfx - (rfz - Z)*dl/dn;
  Double_t Ry = rfy - (rfz - Z)*dm/dn;
  Double_t Rz = Z;

  fphotonInjectLoc[0] = Rx;
  fphotonInjectLoc[1] = Ry;
  fphotonInjectLoc[2] = Rz;

  // distance traveled from top to inject location
  Double_t dist = (rfz - Z)/dn;
  fphotonToTopVolTime = - dist/(TMath::C());
  if (debug) {
    *oLog << "        TopVolPos in focal point coor.  ";
    for (int i = 0;i<3;i++) {
      *oLog << fphotonInjectLoc[i] << " ";
    }
    *oLog << endl;
    *oLog << "        distance from top to inject loc " << dist << endl;
    *oLog << "        fphotonToTopVolTime " << fphotonToTopVolTime << endl;
    *oLog << endl;
  }
  
};
//****************************************************
bool GSegSCTelescope::getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                                           ROOT::Math::XYZVector *photonDcos,
                                           double *photonTime) {
  
  gGeoManager = fManager;

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::getCameraPhotonLocation " << endl;
  }

  photonLoc->SetCoordinates(0.0,0.0,0.0);
  photonDcos->SetCoordinates(0.0,0.0,-1.0);
  *photonTime = 10.0;
  
  //return true;
  double x[4],dir[3];
  for (int i = 0;i<3;i++) {
    x[i] = 0.0;
    dir[i] = 0.0;
  }
  x[3] = 0.0;

  ray->GetLastPoint(x);
  ray->GetDirection(dir);
  fStatusLast = -1;

  // convert distances to mm from cm.
  for (int i = 0;i<3;i++) {
    x[i] = x[i] * 10.0;
  }
  photonLoc->SetCoordinates(x[0],x[1],x[2]);
  photonDcos->SetCoordinates(dir[0],dir[1],dir[2]);
  *photonTime = (x[3] + fphotonToTopVolTime)*1.0e09;
 
  enum {kRun, kStop, kExit, kFocus, kSuspend, kAbsorb};
  
  if (ray->IsExited()) fStatusLast = kExit;
  else if (ray->IsFocused()) {
    fStatusLast = kFocus;
  }
  else if (ray->IsStopped()) fStatusLast = kStop;
  else if (ray->IsSuspended()) fStatusLast = kSuspend;
  else if (ray->IsAbsorbed() ) fStatusLast = kAbsorb;
  
  if (debug) {
    *oLog << "       GetLastPOint " << endl;
    for (int i = 0; i < 3; i++) {
      *oLog << "            " << i << "  loc " << x[i] << endl;
      *oLog << "            " << i << "  dir " << dir[i] << endl;
    }
    *oLog <<   "            time " << x[3] << endl;
    *oLog <<   "            fStatus " << fStatusLast << endl;
    *oLog <<   "            enum for fStatus {kRun, kStop, kExit, kFocus, kSuspend, kAbsorb} " << endl;
    *oLog <<   "            *photonTime " << *photonTime << endl;
  }

  if (bPhotonHistoryFlag) {
    for (int i = 0;i<3;i++) {
      fLocLast[i] = x[i];
      fDirLast[i] = dir[i];
    }
    
    fTimeLast = (x[3] + fphotonToTopVolTime)*1.0e09;
    fillPhotonHistory();
  }

  // the following will draw polylines that end up on the camera/focal plane.
  // you have to instantiate app and do a app.run in grOptics.cpp by 
  // uncommenting these lines in grOptics (shortly after the start of main
  // shortly before the end of main.
  // This also works with the testtel option; however, you have to set
  // nPhotons to a small number, e.g. 10 or less in GArrayTel.cpp (in 
  // the testTelescope method. The TPolyLine3D::Print("all") will print
  // the start of the line, the intermediate points, and the end of the 
  // line.  A good way to test the code.

  // fStatusLast can be 0 through 5.

  if (bRayPlotModeFlag) {
    // draw the telescope only once.
    static int idraw = 1;
    if (idraw) {
      //fManager->GetTopVolume()->Draw("ogl");
      drawTelescope();
      idraw = 0;
    }
    // do we draw the ray.
    if ( ( (eRayPlotType == FOCUSONLY) && (fStatusLast == 3) ) ||
         (eRayPlotType == ALLSURFACES) ) {
      *oLog << " ready to draw polyline" << endl;
      TPolyLine3D *pol = ray->MakePolyLine3D();
      pol->Print("all");
      cout << " fStatusLast " << fStatusLast << endl;
      pol->SetLineColor(2);
      pol->Draw();
    }
  }

  return ray->IsFocused();
};
/********************** end of getCameraPhotonLocation *****************/

//void GSegSCTelescope::setLogFile(const ofstream &logFile) {
//bool debug = false;
//if (debug) {
//  *oLog << " -- GSegSCTelescope::setLogFile " << logFile << endl;
//}
//
//};
/********************** end of setLogFile *****************/

void GSegSCTelescope::printTelescope() {

  bool debug = true;
  if (debug) {
    *oLog << " -- GSegSCTelescope::printTelescope" << endl;
    *oLog << "      iPrtMode " << iPrtMode << endl;
  }
  if (iPrtMode > 0) {
    *oLog << "      fF " << fF << endl;
    *oLog << "      fPlateScaleFactor " << fPlateScaleFactor << endl;
    *oLog << "      fAvgTransitTime   " << fAvgTransitTime << endl;
    *oLog << "      fRotationOffset   " << fRotationOffset << endl;
    *oLog << "      fZp               " << fZp << endl;
    *oLog << "      fZs               " << fZs << endl;
    *oLog << "      number of P1 segments " << vSegP1.size() << endl;
    *oLog << "      number of P2 segments " << vSegP2.size() << endl;
    *oLog << "      number of S1 segments " << vSegS1.size() << endl;
    *oLog << "      number of S2 segments " << vSegS2.size() << endl;
    *oLog << "      fZf                   " << fZf << endl;
    *oLog << "      fRf                   " << fRf << endl;
    if (bCameraFlag) {
      *oLog << "      build MAPMT camera " << endl;
    }
    else {
      *oLog << "      build ideal focal plane " << endl;
    }
    *oLog << "        focalSurfaceOffsets from FOCALSURFACEOFFSET record" << endl;
    *oLog << "           x/y/z " << fFocalSurfaceXOffset << " " 
          << fFocalSurfaceYOffset << " " << fFocalSurfaceZOffset << endl;
    *oLog << "           phi/theta/psi " << fFocalSurfacePhiOffset << "  "
          << fFocalSurfaceThetaOffset << " " << fFocalSurfacePsiOffset << endl;
    
  }

  if ( iPrtMode == 2) {
    // print P1 segment details
    *oLog << "       Number of P1 segments " << iNumP1Mirrors << endl;
    GUtilityFuncts::printSegVector(vSegP1);  
    // print P2 segment details
    *oLog << "       Number of P2 segments " << iNumP2Mirrors << endl;
    GUtilityFuncts::printSegVector(vSegP2);  
    // print S1 segment details
    *oLog << "       Number of S1 segments " << iNumS1Mirrors << endl;
    GUtilityFuncts::printSegVector(vSegS1);  
    // print P2 segment details
    *oLog << "       Number of S2 segments " << iNumS2Mirrors << endl;
    GUtilityFuncts::printSegVector(vSegS2);  
    
  }
    
  if (iPrtMode >0 ) {  

    *oLog << "      ******* Focal Surface " << endl;
    *oLog << "      focal surface " << endl;
    *oLog << "        fKappa1 fKappa2 " << fKappa1 << " " << fKappa2 << " " 
          << fRf << endl;
    *oLog << "        fZf fRf " << fRf << endl;
    *oLog << "      ******* Camera " << endl;
    *oLog << "        fPixelSize  " << fPixelSize << endl;
    *oLog << "        fMAPMTWidth  " << fMAPMTWidth  << endl;
    *oLog << "        fMAPMTLength  " << fMAPMTLength << endl;
    *oLog << "        fInputWindowThickness  " << fInputWindowThickness << endl;
    *oLog << "        fMAPMTOffset  " << fMAPMTOffset << endl;
    *oLog << "        fMAPMTGap  " << fMAPMTGap << endl;
    *oLog << "        fMAPMTRefIndex  " << fMAPMTRefIndex << endl << endl;
    *oLog << "        fCathodeTopRelToFocalSurface "
          << fCathodeTopRelToFocalSurface<< endl;
    *oLog << "        fWindowBottomRelToFocalSurface "
          << fWindowBottomRelToFocalSurface << endl;
    *oLog << "        fMAPOscurationTopRelToFocalSurface  " 
          << fMAPOscurationTopRelToFocalSurface << endl;
    *oLog << "        fCathodeBottomRelToOscurationTop " 
          << fCathodeBottomRelToOscurationTop << endl << endl;
  }
  if (debug) *oLog << "exiting GSegSCTelescope::printTelescope" << endl; 
};
/********************** end of printTelescope *****************/

void GSegSCTelescope::drawTelescope(const int &option) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::drawTelescope" << endl; 
    *oLog << "       option: "  << option << endl;
    if (option > 2) {
      *oLog << "valid options are 0, 1, or 2 " << endl;
    }
  }
  gGeoManager = fManager;

  if ( (option == 0) || (option == 2) ){
    TCanvas *cTelescope = new TCanvas("cTelescope","cTelescope",300,300);
    if (debug) *oLog << "   ready to draw: option " << option << endl;
    fManager->GetTopVolume()->Draw("ogl");
    if (debug) *oLog << "   finished drawing " << endl;
    //gGeoManager->GetTopVolume()->Draw("x3d");
  }
  if (bCameraFlag) {
    if ( (option == 1) || (option == 2)) {
      TCanvas * cCamera = new TCanvas("cCamera","cCamera",300,300);
      gGeoManager->GetVolume("focVol")->Draw("ogl");
      
      if (bSingleMAPMTmodule == false) {
        TCanvas * cMAPMT = new TCanvas("cMAPMT","cMAPMT",300,300);
        gGeoManager->GetVolume("focVol")->GetNode(1)->GetVolume()->Draw("ogl");;
        //gGeoManager->GetVolume("focVol")->GetNode(1)->InspectNode();      
      }
    }
  }

};
/********************** end of drawTelescope *****************/

void GSegSCTelescope::setPrintMode(ostream &oStr,const int prtMode) {
  bool debug = false;
  iPrtMode = prtMode;
  if (debug) {
    *oLog << " -- GSegSCTelescope::setPrintMode " << iPrtMode << endl;
  } 
};
/********************** end of setPrintMode *****************/

void GSegSCTelescope::setPhotonHistory(const string &rootFile,
                                       const string &treeName,
                                       const int &option) {
  bool debug = true;

  
  bPhotonHistoryFlag = true;

  historyFileName = rootFile;
  historyTreeName = treeName;
  iHistoryOption = option;
  
  if (debug) {
    *oLog << "  -- setPhotonHistory " << endl;
    *oLog << "        historyFileName / rootTree " <<  historyFileName 
          << "  " << historyTreeName << " /option " << option 
          << endl;
  }
  // add telID to historyFileName
  //   make string to insert
  stringstream osTmp;
  osTmp << "_Tel" << iTelID;
  string strInsert = osTmp.str();

  // find insertion point, append if . isn't found
  size_t idx = historyFileName.rfind('.');
  
  if (idx != string::npos) {
    historyFileName.insert(idx,strInsert);
    if (debug) *oLog << historyFileName << endl;
  }
  else {
    historyFileName = historyFileName + strInsert;
  }

  if (debug) {
    *oLog << "     opening photon history file / tree:  " 
          << historyFileName << " / " << historyTreeName << endl;
  }
  hisF = new TFile(historyFileName.c_str(),"RECREATE");
  hisT = new TTree(historyTreeName.c_str(),historyTreeName.c_str());

  makePhotonHistoryBranches();
  initializePhotonHistoryParms();
  
};
/********************** end of setPhotonHistory *****************/

void GSegSCTelescope::makePhotonHistoryBranches() {

  hisF->cd();
  bool debug = false;
  if (debug) {
    *oLog << "  -- makePhotonHistoryBranches " << endl;
  }

  hisT->Branch("status",&fStatusLast,"status/I");
  hisT->Branch("nPoints",&fNPoints,"nPoints/I");
  hisT->Branch("injectX",&fInitialInjectLoc[0],"injectX/D");
  hisT->Branch("injectY",&fInitialInjectLoc[1],"injectY/D");
  hisT->Branch("injectZ",&fInitialInjectLoc[2],"injectZ/D");

  hisT->Branch("xLast",&fLocLast[0],"xLast/D");
  hisT->Branch("yLast",&fLocLast[1],"yLast/D");
  hisT->Branch("zLast",&fLocLast[2],"zLast/D");
  hisT->Branch("xLastDir",&fDirLast[0],"xLastDir/D");
  hisT->Branch("yLastDir",&fDirLast[1],"yLastDir/D");
  hisT->Branch("zLastDir",&fDirLast[2],"zLastDir/D");
  hisT->Branch("timeLast",&fTimeLast,"timeLast/D");
};
/************************* end of makePhotonHistoryBranches *****/

void GSegSCTelescope::fillPhotonHistory() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- fillPhotonHistory " << endl;
  }
  hisF->cd();
  hisT->Fill();

};
/************************* end of fillPhotonHistory *****/

void GSegSCTelescope::initializePhotonHistoryParms() {

};
/************************* end of initializePhotonHistoryParms *****/

void GSegSCTelescope::writePhotonHistory() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- in GSegSCTelescope::writePhotonHistory " << endl;
  }
  if (hisF != 0) {
    hisF->cd();  
    hisT->Write();
    // deleting the file will close the file;
    // the file owns the tree and will delete it
    SafeDelete(hisF);
    hisF = 0;
 }
 
};
/************************* end of writePhotonHistory *****/
void GSegSCTelescope::setReflCoeffMap(map<int, TGraph *> *mGr) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::setReflCoeffMap" << endl;
    *oLog << "       mGr->size() " << mGr->size() << endl;
  }
  mGRefl = new map<int, TGraph *>;

  // make a copy of the map
  map<int, TGraph *>::iterator iter;
  for (iter = mGr->begin();iter != mGr->end(); iter++) {
    Int_t id = iter->first;
    if (debug) {
      *oLog << " building new graph " << iter->first << endl;
    }
    TGraph *tmpOld = iter->second;
    TGraph *tmpNew = new TGraph(*tmpOld);
    (*mGRefl)[id] = tmpNew;
  }
  if (debug) {
    *oLog << "      size of mGRefl " << mGRefl->size() << endl;
  }
};
/************************* end of setReflCoeffMap *****/

TGraph * GSegSCTelescope::makeReflectivityGraph(const Int_t &irefl) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::makeReflectivityGraph" << endl;
    *oLog << "    size of mGRefl " << mGRefl->size() << endl;
  }
  Int_t id = irefl;
  TGraph *tmpNew;
  tmpNew = 0;

  map<int, TGraph *>::iterator iter;
  if ( (iter = mGRefl->find(irefl) ) != mGRefl->end() ) {
    TGraph *tmpOld = iter->second;
    tmpNew = new TGraph(*tmpOld);
    if (debug) {
      *oLog << "     ready to print TGraph: reflect curve = " << id << endl;
      tmpNew->Print();
    }

  }
  return tmpNew;
};
/************************* end of makeReflectivityGraph *****/

void GSegSCTelescope::initialize() {
  bool debug = true;

  if (debug) {
    *oLog << "  -- GSegSCTelescope::initialization " << endl;
  }

  fManager = 0;
  fPrimaryV      = 0;
  fSecondaryV    = 0;
  fSecondaryObsV = 0;
  iPrtMode = 0;

  ray = 0;
  hisF = 0;
  hisT = 0;

  iTelID = 0;
  iStdID = 0;

  fTX = 15.0;  // set to 15 later, after confirming code
  fTY = 15.0;
  fTZ = 15.0;

  fTelRadius = 0.0;
  fFMeters = 0.0;

  fAvgTransitTime = 0.0;
  fPlateScaleFactor = 0.0;
  fphotWaveLgt = 0.0;
  fphotonToTopVolTime = 0.0;
  fInjectTime = 0.0; 
  fInjectLambda = 0.0;
  fF = 0.0;
  fAlpha = 0.0;
  fQ = 0.0;
  fRpMax = 0.0; 
  fRpMin = 0.0; 
  fZp = 0.0;
  fP = 0;
  fRsMax = 0.0;
  fRsMin = 0.0;
  fZs = 0.0;
  fNs = 0;
  iNumP1Mirrors = 0;
  iNumP2Mirrors = 0;
  iNParS = 0;
  iNumS1Mirrors = 0;
  iNumS2Mirrors = 0;

  historyFileName = "";
  fTimeLast = 0.0;
  iHistoryOption = 0;
  fStatusLast = 0;
  fNPoints = 0;
  fRotationOffset = 0.0;
  fKappa1 = 0.0;
  fKappa2 = 0.0;
  fRf = 0.0;
  fZf = 0.0;

  bCameraFlag = false;
  bPhotonHistoryFlag = false;
  fPixelSize   = 0.0;
  fMAPMTWidth  = 0.0;
  fMAPMTLength = 0.0;
  fInputWindowThickness = 0.0;
  fMAPMTOffset = 0.0;
  fMAPMTGap    = 0.0;
  fMAPMTRefIndex        = 0.0;
  bSingleMAPMTmodule = false;

  fCathodeTopRelToFocalSurface       = 0.0;
  fWindowBottomRelToFocalSurface     = 0.0;
  fMAPOscurationTopRelToFocalSurface = 0.0;
  fCathodeBottomRelToOscurationTop   = 0.0;

  fFocalSurfaceXOffset     = 0.0;
  fFocalSurfaceYOffset     = 0.0;
  fFocalSurfaceZOffset     = 0.0;
  fFocalSurfacePhiOffset   = 0.0;
  fFocalSurfaceThetaOffset = 0.0;
  fFocalSurfacePsiOffset   = 0.0;  

  bRayPlotModeFlag = false;
  eRayPlotType = FOCUSONLY;  

  eTelType = SEGSC;

  // set gl picture colors: black/1, red/2/, green/3,blue/4, brown/28/49
  iPrimaryColor =38;
  iSecondaryColor = iPrimaryColor;
  iSecondaryObscurationColor = 1;
  iMAPMTObscurationColor = iSecondaryObscurationColor;
  iMAPMTCathodeColor = 2;
  iMAPMTWindowColor = 3; 

  for (int i = 0;i<3;i++) {
    fphotonInjectLoc[i] = 0.0;
    fphotonInjectDir[i] = 0.0;
    fInjectLoc[i] = 0.0;
    fInjectDir[i] = 0.0;
    fLocLast[i] = 0.0;
    fDirLast[i] = 0.0;
    fInitialInjectLoc[i] = 0.0;
  }
};
/************************* end of writePhotonHistory *****/

void GSegSCTelescope::CloseGeometry() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::CloseGeometry" << endl;
  }
  gGeoManager = fManager;
  fManager->CloseGeometry();
}; 
/************************* end of CloseGeometry *****/

void GSegSCTelescope::testPerformance() {

  // unused - from Akira's NewSCT.C, may not work here as is.
  /*
  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::testPerformance " << endl;
  }

  gGeoManager = fManager;

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
    fManager->TraceNonSequential(*array);
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
  */
};
/************************* end of testPerformance *****/

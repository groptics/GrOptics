/*
VERSION2.7aBeta
17Dec2012
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

#include "AOpticsManager.h"
#include "ABorderSurfaceCondition.h"

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

#include "GUtilityFuncts.h"
#include "GDefinition.h"

#include "GSegmentedMirror.h"
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

  bool debug=true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::~GSegSCTelescope " << endl;
  }
  if (fManager !=0) {
    gGeoManager = fManager;
    SafeDelete(fManager);
  }
 
};
/********************** end of ~GSegSCTelescope *****************/

void GSegSCTelescope::buildTelescope(bool os8)
{

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::buildTelescope" << endl;
  }
  gGeoManager = 0;
  fManager = new AOpticsManager("manager","The optics manager of SEGSC");
  fManager->SetVisLevel(5);
  fManager->SetNsegments(50);

  // Make dummy material
  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(70); // needed in OpenGL view
  new TGeoMedium("med", 1, mat);

  // Make the world
  TGeoBBox* worldbox = new TGeoBBox("worldbox", fTX*m, fTY*m, fTZ*m);
  AOpticalComponent* world = new AOpticalComponent("world", worldbox);
  fManager->SetTopVolume(world);
  //manager->SetTopVisible(0); 
  //fManager->SetTopVisible(.5); 

  makePrimarySecondaryDisks();
  addPrimaryF();
  
  closeGeometry();
  return;
};
/*************************************************************************************/
void GSegSCTelescope::makePrimarySecondaryDisks() {

  bool debug = true;

  if (debug) {
    *oLog << "  --  GSegSCTelescope::makePrimarySecondaryDisks" << endl;
  }

  const Double_t kZp = fZp;
  const Double_t kZs = (fF)/fQ;
  if (debug) {
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

  if (debug) {
    *oLog << "   primary coefficients " << endl;
    *oLog << "   fNp " << fNp << endl;
    for (int i = 0;i<fNp; i++ ) {
      *oLog << "   i/fNp " << i << " " << fP[i] << endl;
    }
  }
  fPrimaryV->SetPolynomials(fNp - 1, &fP[1], fNp - 1, &fP[1]);

  // Make the ideal volume of the secondary mirror
  if (debug) {
    *oLog << "   secondary coefficients " << endl;
    *oLog << "   fNs " << fNs << endl;
    for (int i = 0;i<fNs; i++ ) {
      *oLog << "   i/fNs " << i << " " << fS[i] << endl;
    }
  }
  fSecondaryV = new AGeoAsphericDisk("secondaryV",
                                     kZs + fS[0], 0, 
                                     kZs + fS[0]  + 1*um, 
                                     0, fRsMax, fRsMin);
  fSecondaryV->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  return;
};
/*************************************************************************************/

void GSegSCTelescope::addPrimaryF() {
  bool debug = true;

  if (debug) {
    *oLog << "  --  GSegSCTelescope::addPrimaryF" << endl;
    /*
    *oLog << "     primaryF parameters " << endl;
    *oLog << "       rmin/rmax " <<  (*(vSegP1.at(0))).rmin << "  "
          << (*(vSegP1.at(0))).rmax << endl;
    
    *oLog << "       delPhi " << (*(vSegP1.at(0))).delPhi << endl;
    *oLog << "       margin " << (*(vSegP1.at(0))).margin << endl;
    *oLog << "       pos.errors " << (*(vSegP1.at(0))).posErrorX << "  "
          << (*(vSegP1.at(0))).posErrorY << "  " << (*(vSegP1.at(0))).posErrorZ
          << endl;
    *oLog << "       rot.errors " << (*(vSegP1.at(0))).rotErrorX << "  "
          << (*(vSegP1.at(0))).rotErrorY << "  " << (*(vSegP1.at(0))).rotErrorZ
          << endl;
    *oLog << "       reflect curve " << (*(vSegP1.at(0))).reflect << endl;
    *oLog << "       iNumP1Mirrors " << iNumP1Mirrors << endl;
    */
   }
  Int_t count = 0;

  // segment not created if reflect = 0
  // P1 mirrors
  for (Int_t i = 0; i < iNumP1Mirrors; i++) {
    if ( ( (*(vSegP1.at(i))).reflect) > 0) {
      Double_t rmin = (*(vSegP1.at(0))).rmin;
      Double_t rmax = (*(vSegP1.at(0))).rmax;
      Double_t margin = (*(vSegP1.at(0))).margin;
      rmin = rmin*m - margin/TMath::Cos(11.25/2.*TMath::DegToRad());
      rmax = rmax*m;
      
      Double_t phimin = ( (*(vSegP1.at(0))).delPhi)*i;
      Double_t phimax = ( (*(vSegP1.at(0))).delPhi)*(i+1);
      if (debug) {
        *oLog << "    P1 pentagon segmented mirror number " << i << endl;
        *oLog << "        rmin/rmax " << rmin << " " << rmax << endl;
        *oLog << "        margin    " << margin << endl;
        *oLog << "        phimin/max " << phimin << "  " << phimax << endl;
      }
    PentagonSegmentedMirror mirror(rmin + margin, rmax, phimin, phimax);
    mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
    mirror.SetRotationErrors(0., 0., 0.);
    mirror.SetRougness(0.);
    mirror.SetMargin(margin);
    // add mirror segment
    addPrimaryMirror(Form("primary%d", count), &mirror);
    count++;
    } 
  }
  /*
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
  */
};
/*******************************************************************/
void GSegSCTelescope::addPrimaryMirror(const char*name,
                                       SegmentedMirror *mirror) {

  gGeoManager = fManager;
  
  bool debug = true;
  if (debug) {
    *oLog << "  --  GSegSCTelescope::addPrimaryMirror" << endl;
  }
  
  AMirror* mir = mirror->BuildMirror(name, fPrimaryV, kTRUE);
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fPrimaryV, kTRUE);

  ABorderSurfaceCondition * condition
    = new ABorderSurfaceCondition(fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  fManager->GetTopVolume()->AddNode(mir, 1, combi);
  
};
/*******************************************************************/

void GSegSCTelescope::addSecondaryJ() {
  Int_t count = 0;
  /*
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
  */
};
/*******************************************************************/

void GSegSCTelescope::AddSecondaryObscuration() {
  /*
  const Double_t kZs = fF/fQ;
  AGeoAsphericDisk* disk
    = new AGeoAsphericDisk("secondaryObsV", kZs + fS[0] + 1.*cm, 0, kZs + fS[0]  + 1*cm, 0, fRsMax, 0);
  disk->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  TGeoMedium* med = fManager->GetMedium("med");
  AObscuration* secondaryObs = new AObscuration("secondaryObs", disk, med);
  fManager->GetTopVolume()->AddNode(secondaryObs, 1);
  */
};

/*******************************************************************/
void GSegSCTelescope::addSecondaryMirror() {
  /*
  AMirror* mir = mirror->BuildMirror(name, fSecondaryV, kFALSE);
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fSecondaryV, kFALSE);

  ABorderSurfaceCondition * condition
    = new ABorderSurfaceCondition(fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());

  fManager->GetTopVolume()->AddNode(mir, 1, combi);
  */
};
/*************************************************************************************/
void GSegSCTelescope::addIdealFocalPlane()  {
  /*
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
  */
};
/*************************************************************************************/

void GSegSCTelescope::addMAPMTFocalPlane()  {
  /*
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
*/
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
  //gGeoManager = manager;

  bool debug = true;
  if (debug) {
    *oLog << " -- GSegSCTelescope::injectPhoton " << endl;
    *oLog << "       photonLocT  ";
    GUtilityFuncts::printGenVector(photonLocT); *oLog << endl;
    *oLog << "       photonDirT  ";
    GUtilityFuncts::printGenVector(photonDirT); *oLog << endl;
    *oLog << "      bPhotonHistoryFlag " << bPhotonHistoryFlag << endl;
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
  
  //gGeoManager = manager;

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::getCameraPhotonLocation " << endl;
  }
  // just return zeros
  photonLoc->SetCoordinates(0.0,0.0,0.0);
  photonDcos->SetCoordinates(0.0,0.0,-1.0);
  *photonTime = 10.0;
  
  return true;

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
    *oPrtStrm << " -- GSegSCTelescope::printTelescope" << endl;
  }
  *oLog << "    fF " << fF << endl;
};
/********************** end of printTelescope *****************/

void GSegSCTelescope::drawTelescope() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::drawTelescope" << endl; 
  }
  gGeoManager = fManager;
  //gGeoManager->GetMasterVolume()->Draw("ogl");
  gGeoManager->GetTopVolume()->Draw("ogl");
};
/********************** end of drawTelescope *****************/

void GSegSCTelescope::setPrintMode(ostream &oStr,const int prtMode) {
  bool debug = true;
  if (debug) {
    *oLog << " -- GSegSCTelescope::setPrintMode" << endl;
  } 
};
/********************** end of setPrintMode *****************/

void GSegSCTelescope::setPhotonHistory(const string &rootFile,
                                       const string &treeName,
                                       const int &option) {
  bool debug = true;

};
/********************** end of setPhotonHistory *****************/

void GSegSCTelescope::makePhotonHistoryBranches() {

  //hisF->cd();
  bool debug = true;
  if (debug) {
    *oLog << "  -- makePhotonHistoryBranches " << endl;
  }
 
};
/************************* end of makePhotonHistoryBranches *****/

void GSegSCTelescope::fillPhotonHistory() {

  // careful have to move to correct hisF ?, do a changedirectory?

  bool debug = true;
  if (debug) {
    *oLog << "  -- fillPhotonHistory " << endl;
  }

};
/************************* end of fillPhotonHistory *****/

void GSegSCTelescope::initializePhotonHistoryParms() {

};
/************************* end of initializePhotonHistoryParms *****/

void GSegSCTelescope::writePhotonHistory() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- in GSegSCTelescope::writePhotonHistory " << endl;
  }
 
};
/************************* end of writePhotonHistory *****/
void GSegSCTelescope::setReflCoeffMap(map<int, TGraph *> *mGRefl) {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::setReflCoeffMap" << endl;
  }
  // make a copy of the map
  map<int, TGraph *>::iterator iter;
  for (iter = mGRefl->begin();iter != mGRefl->end(); iter++) {
    Int_t id = iter->first;
    TGraph *tmpOld = iter->second;
    TGraph *tmpNew = new TGraph(*tmpOld);
    (*mGRefl)[id] = tmpNew;
  }
  if (debug) {
    *oLog << "      size of mGRefl " << mGRefl->size() << endl;
  }
};
/************************* end of writePhotonHistory *****/

void GSegSCTelescope::initialize() {
  bool debug = true;

  if (debug) {
    *oLog << "  -- GSegSCTelescope::initialization " << endl;
  }

  fManager = 0;
  ray = 0;
  hisF = 0;
  hisT = 0;

  iTelID = 0;
  iStdID = 0;

  fTX = 30.0;  // set to 15 later, after confirming code
  fTY = 30.0;
  fTZ = 30.0;

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
  iNParP = 0;
  iNumP1Mirrors = 0;
  iNumP2Mirrors = 0;
  iNParS = 0;
  iNumS1Mirrors = 0;
  iNumS2Mirrors = 0;
  gPrimRefl = 0;
  fPrimMaxLmda = 0.0;
  fPrimMinLmda = 0.0;
  gSeconRefl = 0;
  fSeconMaxLmda = 0.0;
  fSeconMinLmda = 0.0;
  historyFileName = "";
  fTimeLast = 0.0;
  iHistoryOption = 0;
  fStatusLast = 0;
  fNPoints = 0;
  fRotationOffset = 0.0;

  eTelType = SEGSC;

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

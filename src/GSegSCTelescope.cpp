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
  TGeoBBox* worldbox = new TGeoBBox("worldbox", 30*m, 30*m, 30*m);
  AOpticalComponent* world = new AOpticalComponent("world", worldbox);
  fManager->SetTopVolume(world);

  const Double_t kZp = fZp*m;
  const Double_t kZs = fF/fQ;
  fPrimaryV = new AGeoAsphericDisk("primaryV",
                                   kZp + fP[0] - 1*um, 0,
                                   kZp + fP[0] , 0, 
                                   fRpMax, fRpMin);
  fPrimaryV->SetPolynomials(fNp - 1, &fP[1], fNp - 1, &fP[1]);

  // Make the ideal volume of the secondary mirror
  fSecondaryV = new AGeoAsphericDisk("secondaryV",
                                     kZs + fS[0], 0, 
                                     kZs + fS[0]  + 1*um, 
                                     0, fRsMax, fRsMin);
  fSecondaryV->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  return;
};
/*************************************************************************************/
void GSegSCTelescope::testFocalPlane() {

  // inject photons 

};

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

  //gGeoManager = manager;

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::movePositionToTopOfTopVol " << endl;
    *oLog << "        position prior to move to top ";
    *oLog << fphotonInjectLoc[0] << "  " << fphotonInjectLoc[1] << "  " 
	  << fphotonInjectLoc[2] << endl;
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
  gGeoManager->GetMasterVolume()->Draw("ogl");
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

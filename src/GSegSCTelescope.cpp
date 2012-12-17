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

  };
/********************** end of GSegSCTelescope *****************/
GSegSCTelescope::~GSegSCTelescope() {

  bool debug=true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::~GSegSCTelescope " << endl;
  }

 
};
/********************** end of ~GSegSCTelescope *****************/

void GSegSCTelescope::buildTelescope(bool os8)
{

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::buildTelescope" << endl;
  }
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

};
/********************** end of printTelescope *****************/

void GSegSCTelescope::drawTelescope() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::drawTelescope" << endl; 
  }

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
void GSegSCTelescope::setReflCv(const int &primID, const int &secondID,
			     map<int, TGraph *> *mGRefl) {
  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescope::setReflCv" << endl;
  }
  
};
/************************* end of writePhotonHistory *****/

void GSegSCTelescope::initialization() {
  bool debug = true;

  if (debug) {
    *oLog << "  -- GSegSCTelescope::initialization " << endl;
  }

};

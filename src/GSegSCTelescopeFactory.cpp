/*
VERSION3.1
2March2015
*/
/*!  GSegSCTelescopeFactory.cpp

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
#include <assert.h>

#include "TROOT.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"

#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "AOpticsManager.h" 
#include "GTelescope.h"
#include "GSegSCTelescope.h"
#include "GSCTelescope.h"
#include "GTelescopeFactory.h"
#include "GSCTelescopeFactory.h"
#include "GSegSCTelescopeFactory.h"

#include "GReadSegSCStd.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

// define useful units
static const double cm = AOpticsManager::cm();
static const double mm = AOpticsManager::mm();
static const double um = AOpticsManager::um();
static const double nm = AOpticsManager::nm();
static const double  m = AOpticsManager::m();

/*!  /brief SegSCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
SegSCStdOptics::SegSCStdOptics() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::SegSCStdOptics " << endl;
  }

  iPrimReflID = 0;
  iSecReflID = 0;
  iPrtMode = 0;
  oStr     = oLog;

  stdType = SEGSC; 
  stdNum = 0;
  fAvgTransitTime = 0.0;
  fRotationOffset = 0.0;
  fPlateScaleFactor = 0.0;
  
  fF = 0.0;     // Focal length
  fAlpha = 2./3.; // \alpha
  fQ = 2./3.;     // q
  fRpMax = 0.0; // Primary radius max
  fRpMin = 0.0; // Primary radius min
  fRsMax = 0.0; // Secondary radius max
  fRsMin = 0.0; // Secondary radius min
  fKappa1 = 0.0;// Focal plane sag constant
  fKappa2 = 0.0;// Focal plane sag constant
  fRf = 0.0;

  fZp = 0.0;
  fZs = 1.5;
  fZf = (1./fQ) - (1 - fAlpha);
  
  fNp = 0;   // Number of coefficients for the primary
  fNs = 0;   // Number of coefficients for the secondary
  fP = 0;   // Polynomial coefficients (p0, p1 ...)
  fS = 0;   // Polynomial coefficients (s0, s1 ...)
    
  // MAPMT Parameters
  bCameraFlag = false;
  fPixelSize = 0.0;            // Width of a MAPMT pixel
  fMAPMTWidth = 0.0;           // Housing width
  fMAPMTLength = 0.0;          // between input window and anode pins
  fInputWindowThickness = 0.0; // Thickness of the input window
  fMAPMTGap = 0.0;
  fMAPMTRefIndex = 0.0;
  fMAPMTOffset = 0.0;
  bSingleMAPMTmodule = false;

  iNParP = 0;
  iNParS = 0;

  iNumP1Mirrors = 0;
  iNumP2Mirrors = 0;
  iNumS1Mirrors = 0;
  iNumS2Mirrors = 0;  

  fFocalSurfaceXOffset     = 0.0;
  fFocalSurfaceYOffset     = 0.0;
  fFocalSurfaceZOffset     = 0.0;
  fFocalSurfacePhiOffset   = 0.0;
  fFocalSurfaceThetaOffset = 0.0;
  fFocalSurfacePsiOffset   = 0.0;  
};
/************** end of SegSCStdOptics ***********************/

SegSCStdOptics::SegSCStdOptics(const SegSCStdOptics &sco) {

  bool debug = true;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::SegSCStdOptics " << endl;
  }
  

};

/************** end of SegSCStdOptics ***********************/
SegSCStdOptics::~SegSCStdOptics() {
  bool debug = true;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::~SegSCStdOptics" << endl;
  }

  for (int unsigned i = 0;i< vSegP1.size();i++) {
    SafeDelete(vSegP1.at(i));
  }
  for (int unsigned i = 0;i< vSegP2.size();i++) {
    SafeDelete(vSegP2.at(i));
  }
  for (int unsigned i = 0;i< vSegS1.size();i++) {
    SafeDelete(vSegS1.at(i));
  }
  for (int unsigned i = 0;i< vSegS2.size();i++) {
    SafeDelete(vSegS2.at(i));
  }
};
/********************* end of ~SegSCStdOptics *****************/

void SegSCStdOptics::printSegSCStdOptics() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::printSegSCStdOptics " << endl;
  }
  if (iPrtMode == 0) return;
  *oLog << "    SegSCStdOptics::printSegSCStdOptics() " << endl;
  *oLog << "        stdType " << stdType << endl;
  *oLog << "        telType " << getTelType(stdType) << endl;
  *oLog << "        stdNum  " << stdNum << endl;
  *oLog << "        fF      " << fF << endl;
  *oLog << "        fAlpha      " << fAlpha << endl;
  *oLog << "        fQ      " << fQ << endl;
  *oLog << "        fPlateScaleFactor " <<  fPlateScaleFactor << endl;
  *oLog << "        fAvgTransitTime " << fAvgTransitTime << endl;
  *oLog << "        fRotationOffset " << fRotationOffset << endl;
  *oLog << "        Primary  fRpMax/Min fZp " << fRpMax << " " 
        << fRpMin  << " " << fZp << endl;
  *oLog << "        Secondary fRsMax/Min fZs " << fRsMax << " " 
        << fRsMin  << " " << fZs << endl;
  *oLog << "        iNParP " << iNParP << endl;
  if (iNParP > 0) {
    *oLog << "          primPolyCoeff " << endl;
    for (int i = 0;i<iNParP;i++) {
      *oLog << "               " << i << "    "  << fzp[i] << endl;
    }
  }
  
  *oLog << "        iNParS " << iNParS << endl;
  if (iNParS > 0) {
    *oLog << "          secondaryPolyCoeff " << endl;
    for (int i = 0;i<iNParS;i++) {
      *oLog << "               " << i << "    " << fzs[i] << endl;
    }
  }
  *oLog << "        Focal Surface fKappa1 fKappa2 fRf " 
        << fKappa1 << " " << fKappa2 << " " 
        << fRf << endl;
  
  *oLog << "        bCameraFlag  " << bCameraFlag << endl;
  *oLog << "        fPixelSize   " << fPixelSize << endl;
  *oLog << "        fMAPMTWidth   " << fMAPMTWidth << endl;
  *oLog << "        fMAPMTLength   " << fMAPMTLength << endl;
  *oLog << "        fInputWindowThickness   " << fInputWindowThickness << endl;
  *oLog << "        fMAPMTOffset        " << fMAPMTOffset << endl;
  *oLog << "        fMAPMTGap        " << fMAPMTGap << endl;
  *oLog << "        fMAPMTRefIndex        " << fMAPMTRefIndex << endl;
  *oLog << "        bSingleMAPMTmodule " << bSingleMAPMTmodule << endl;
  *oLog << "        focalSurfaceOffsets from FOCALSURFACEOFFSET record" << endl;
  *oLog << "           x/y/z " << fFocalSurfaceXOffset << " " 
        << fFocalSurfaceYOffset << " " << fFocalSurfaceZOffset << endl;
  *oLog << "           phi/theta/psi " << fFocalSurfacePhiOffset << "  "
        << fFocalSurfaceThetaOffset << " " << fFocalSurfacePsiOffset << endl;

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

  *oLog << "       iPrtMode " << iPrtMode << endl;
    
};
/************** end of printSegSCStdOptics ***********************/
void SegSCStdOptics::printSegVector (const vector<mirrorSegmentDetails *> &vec) {
  *oLog << "        num rmin    rmax marg dPhi refl rough errXYZ rotErrABC"
        << endl;
  int numElem = vec.size();
  for (int i = 0;i<numElem; i++) {
    mirrorSegmentDetails *t = vec[i];
    *oLog << "         " << i+1 << "  " 
          << t->rmin << "  " << t->rmax << "  " << t->margin << "  " 
          << t->delPhi << "  " << t->reflect << "   " << t->roughness
          << "      " << t->posErrorX
          << " " << t->posErrorY << " " << t->posErrorZ << "    " 
          << t->rotErrorPhi << " " << t->rotErrorTheta << " " << t->rotErrorPsi
          << endl;
  }

};
//////////////////////////////////////////////////////////////////////////

GSegSCTelescopeFactory::
GSegSCTelescopeFactory(GReadSegSCStd &scReader,
		    const string &editPilotFile) {

  bool debug = true;
  sPilotEdit = editPilotFile;
  if (debug) {
    *oLog << "  -- GSegSCTelescopeFactory Constructor:  " << sPilotEdit << endl;
  }
  iNumSCTelMade = 0;
  readSegSC = 0;
  pi = 0;
  sPilotEdit = "";
  mGRefl = 0;
  opt = 0;
  SCTel = 0;

  // make the reflectivity map 
  mGRefl = new map<int, TGraph *>;

  readSegSC = &(scReader);
  readSegSC->setSegSCTelescopeFactory(this);

  sPilotEdit = editPilotFile;

  pi = new GPilot(sPilotEdit);
 
  if (debug) {
    *oLog << "    -- end of GSCTelescopeFactory constructor" << endl;
  }
  
};

/************** end of GSegSCTelescopeFactory ***********************/

GSegSCTelescopeFactory::~GSegSCTelescopeFactory() {
 
  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescopeFactory::~GSegSCTelescopeFactory" << endl;
  }

  for (itmStdOp=mStdOptics.begin();
     itmStdOp!=mStdOptics.end(); itmStdOp++) {
    // WHY CAN'T I DELETE THIS CLASS? thinks I have already deleted the class?
    // was created in reader: SCFac->mStdOptics[i] = new SegSCStdOptics();
    //SegSCStdOptics *tmpSeg;
    //tmpSeg = itmStdOp->second;
    //*oLog << "tmpSeg->stdNum " <<tmpSeg->stdNum << endl;
    //delete tmpSeg;
    SafeDelete( itmStdOp->second );
    //*oLog << " (itmStdOp->second)->stdNum " <<  (itmStdOp->second)->stdNum << endl;
    //*oLog << "itmStdOp->second " << itmStdOp->second << endl;
  }
   for (itmGRefl=mGRefl->begin();
    itmGRefl!=mGRefl->end(); itmGRefl++) {
     SafeDelete(itmGRefl->second ); 
   }
  SafeDelete(mGRefl);
  SafeDelete(pi);
  SafeDelete(readSegSC);
 
};
/************** end of ~GSegSCTelescopeFactory ***********************/

GSegSCTelescope* GSegSCTelescopeFactory::makeTelescope(const int &id,
                                                     const int &std) {
  
  int debug = true;
  if (debug) {
    *oLog << " -- GSegSCTelescopeFactory::makeTelescope" << endl;
    //*oLog << "      telID  = " << id << endl;
    //*oLog << "      telStd = " << std << endl;
  }

  Int_t idTel = id;
  Int_t iStdID = std;

  // get parameters for this telescope
  itmStdOp = mStdOptics.find(iStdID);
  assert(itmStdOp != mStdOptics.end());
  // make pointer to working stdoptics structure, easier typing
  opt = mStdOptics[iStdID];

  // make the telescope
  SCTel = new GSegSCTelescope;

  SCTel->iPrtMode = opt->iPrtMode;

  // move over all reflection coefficients (the entire map)
  SCTel->setReflCoeffMap(mGRefl);
  
  SCTel->setTelID(idTel);
  SCTel->setStdID(iStdID);
  SCTel->eTelType = opt->stdType;
  SCTel->fAvgTransitTime = opt->fAvgTransitTime;
  SCTel->fRotationOffset = opt->fRotationOffset;
  SCTel->fPlateScaleFactor = opt->fPlateScaleFactor;

  // general telescope parameters
  SCTel->fF = (opt->fF);
  SCTel->fAlpha = opt->fAlpha;
  SCTel->fQ = opt->fQ;

  // primary parameters
  SCTel->fRpMax = (opt->fRpMax);
  SCTel->fRpMin = (opt->fRpMin);
  SCTel->fZp = (opt->fZp);
  //SCTel->fNp = opt->fNp;
  SCTel->fNp = opt->iNParP;
  SCTel->fzp = opt->fzp;
  // make primary poly coefficients
  SCTel->fP = new Double_t[SCTel->fNp];
  Double_t fF = (SCTel->fF)*m;
  
  (SCTel->fP)[0] = TMath::Power(fF,  1)* ((opt->fzp[0]));

  bool printParameters = true;
  if (printParameters) {
    *oLog << endl;
    *oLog << "  calculation of primary parameters " << endl;
    *oLog << "           fzp[i]: input parameters from configuration file " << endl;
    *oLog << "           fP[i]  =  TMath::Power(fF,  iPowerF)* fzp[i] " << endl;  
    *oLog << "  i     iPowerF     fzp[i]           fP[i] " << endl;
    *oLog <<"  0" << "        " << "1" << "       " << opt->fzp[0] << "       " 
          << (SCTel->fP)[0] << endl;   
  }
  int iPowerF = -1;
  for (int i = 1;i < SCTel->fNp; i++) {
    (SCTel->fP)[i] = TMath::Power(fF,  iPowerF)* ((opt->fzp[i]));  
    if (printParameters) {
      *oLog << "  " <<  i << "       " << iPowerF
            << "       " << opt->fzp[i] << "       " << (SCTel->fP)[i] << endl;
    }
    iPowerF = iPowerF -2;
  }
  if(printParameters) *oLog << endl;
  
  // secondary parameters
  SCTel->fRsMax = (opt->fRsMax);
  SCTel->fRsMin = (opt->fRsMin);
  SCTel->fZs = (opt->fZs);
  SCTel->fNs = opt->fNs;
  SCTel->fS = opt->fS;
  SCTel->iNParS = opt->iNParS;
  SCTel->fNs = opt->iNParS;
  SCTel->fzs = opt->fzs;
  // make secondary poly coefficients
  SCTel->fS = new Double_t[SCTel->iNParS];
  fF = (SCTel->fF)*m;

  (SCTel->fS)[0] = TMath::Power(fF,  1)* ((opt->fzs[0]));

  if (printParameters) {
    *oLog << "  calculation of secondary parameters " << endl;
    *oLog << "           fzs[i]: input parameters from configuration file " << endl;
    *oLog << "           fS[i]  =  TMath::Power(fF,  iPowerF)* fzs[i] " << endl;  
    *oLog << "  i     iPowerF     fzs[i]           fS[i] " << endl;
    *oLog <<"  0" << "        " << "1" << "       " << opt->fzs[0] << "       " 
          << (SCTel->fS)[0] << endl;   
  }

  iPowerF = -1;
  for (int i = 1;i < SCTel->fNs; i++) {
    (SCTel->fS)[i] = TMath::Power(fF,  iPowerF)* ((opt->fzs[i]));  
    if (printParameters) {
      *oLog << "  " <<  i << "       " << iPowerF
            << "       " << opt->fzs[i] << "       " << (SCTel->fS)[i] << endl;
    }
    iPowerF = iPowerF -2;
  }
  if (printParameters) *oLog << endl;

  // primary segment details
  SCTel->iNumP1Mirrors = opt->iNumP1Mirrors;
  SCTel->iNumP2Mirrors = opt->iNumP2Mirrors;
  SCTel->vSegP1 = opt->vSegP1;
  SCTel->vSegP2 = opt->vSegP2;

 // secondary segment S1details
  SCTel->iNumS1Mirrors = opt->iNumS1Mirrors;
  SCTel->iNumS2Mirrors = opt->iNumS2Mirrors;
  SCTel->vSegS1 = opt->vSegS1;
  SCTel->vSegS2 = opt->vSegS2;

  // focal plane
  SCTel->fKappa1 = opt->fKappa1;
  SCTel->fKappa2 = opt->fKappa2;
  SCTel->fRf     = opt->fRf;
  SCTel->fZf     = opt->fZf;

  SCTel->bCameraFlag = opt->bCameraFlag;

  // camera
  SCTel->fPixelSize   = opt->fPixelSize;
  SCTel->fMAPMTWidth  = opt->fMAPMTWidth;
  SCTel->fMAPMTLength = opt->fMAPMTLength;
  SCTel->fInputWindowThickness   = opt->fInputWindowThickness;
  SCTel->fMAPMTOffset = opt->fMAPMTOffset;
  SCTel->fMAPMTGap    = opt->fMAPMTGap;
  SCTel->fMAPMTRefIndex   = opt->fMAPMTRefIndex;
  SCTel->bSingleMAPMTmodule = opt->bSingleMAPMTmodule;

  SCTel->fFocalSurfaceXOffset     = opt->fFocalSurfaceXOffset;
  SCTel->fFocalSurfaceYOffset     = opt->fFocalSurfaceYOffset;
  SCTel->fFocalSurfaceZOffset     = opt->fFocalSurfaceZOffset;
  SCTel->fFocalSurfacePhiOffset   = opt->fFocalSurfacePhiOffset;
  SCTel->fFocalSurfaceThetaOffset = opt->fFocalSurfaceThetaOffset;
  SCTel->fFocalSurfacePsiOffset   = opt->fFocalSurfacePsiOffset;
 
  SCTel->buildTelescope();

  return SCTel; 
};
/************** end of makeTelescope ***********************/

void GSegSCTelescopeFactory::editWorkingTelescope(GSegSCTelescope *SCTel1) {
  bool debug = true;
  if (debug) {
    *oLog << " -- GSegSCTelescopeFactory::editWorkingTelescope" << endl;
  }

 
};
/************** end of editWorkingTelescope ***********************/

void GSegSCTelescopeFactory::printStdTelescope(const int &iStd, 
                                               const int &mode,ostream &oStr) {

  // DO NOT USE.
  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GSegSCTelescopeFactory::printStdTelescope" << endl;
  }
  
};
/************** end of :printStdTelescope ***********************/

void GSegSCTelescopeFactory::setPrintMode(ostream &oStr,
                                          const int prtMode) {

  int debug = 1;
  if (debug > 0) {
    *oLog << " -- GSegSCTelescopeFactory::setPrintMode" << endl;
  }
 
}; 
/************** end of setPrintMode  ***********************/




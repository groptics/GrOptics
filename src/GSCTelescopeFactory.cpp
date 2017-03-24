/*
VERSION4.0
30May2016
*/
/*!  GSCTelescopeFactory.cpp

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

#include "TGraph.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"

#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "AOpticsManager.h" 
#include "ATelescope.h"
#include "GSCTelescope.h"
#include "ATelescopeFactory.h"
#include "GSCTelescopeFactory.h"

#include "GReadSCStd.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

/*!  /brief SCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
SCStdOptics::SCStdOptics() {

  iPrimReflID = 0;
  iSecReflID  = 0;
  stdNum      = 0;
  fAvgTransitTime = 0;
  
  oStr     = oLog;
  iPrtMode = 0;
  stdType  = SC; 
  fRotationOffset = 0.0;
  fFocLgt = 0.0;
  fPlateScaleFactor = 0.0;

  // primary
  fDp = 0.0;
  fDpinner = 0.0;
  fZp = 0.0;
  fPrimXoffset = 0.0;
  fPrimYoffset = 0.0;
  fPrimZoffset = 0.0;
  fPrimThetaOffset = 0.0;
  fPrimPhiAngle    = 0.0;
  fPrimRoughSigma  = 0.0;
  fPrimRoughMax    = 0.0;

  // secondary
   fZs = 0.0;
   fDs = 0.0;
   fDsinner = 0.0;
   fSecondXoffset     = 0.0;
   fSecondYoffset     = 0.0;
   fSecondZoffset     = 0.0;
   fSecondThetaOffset = 0.0;
   fSecondPhiAngle    = 0.0;
   fSecondRoughSigma  = 0.0;
   fSecondRoughMax    = 0.0;
  
   // focal plane
   fZf = 0.0;
   fk1 = 0.0;
   fk2 = 0.0;
   fFocalPlXoffset     = 0.0;
   fFocalPlYoffset     = 0.0;
   fFocalPlZoffset     = 0.0;
   fFocalPlThetaOffset = 0.0;
   fFocalPlPhiAngle    = 0.0;
   
   // camera
   fPixelSize = 0.0;
   fPixelPitch = 0.0;
   fMAPMTWidth = 0.0;
   fMAPMTLength = 0.0;
   fInputWindowThickness = 0.0;
   fMAPMTAngularSize = 0.0;
   fMAPMTOffset = 0.0;
   fMAPMTGap = 0.0;
   fMAPMTRefIndex = 0.0; 
   bCameraFlag = true;

   iNParP = 0;
   iNParS = 0;

};
/************** end of SCStdOptics ***********************/

SCStdOptics::SCStdOptics(const SCStdOptics &sco) {
 
  oStr = sco.oStr;
  iPrtMode = sco.iPrtMode;
  stdType = sco.stdType;
  fAvgTransitTime = sco.fAvgTransitTime;
  fRotationOffset = sco.fRotationOffset;
};

/************** end of SCStdOptics ***********************/
SCStdOptics::~SCStdOptics() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- SCStdOptics::~SCStdOptics" << endl;
  }
};
/********************* end of ~SCStdOptics *****************/

void SCStdOptics::printSCStdOptics() {

 if (iPrtMode == 0) return;
 *oLog << "    SCStdOptics::printSCStdOptics() " << endl;
 *oLog << "        telType " << getTelType(stdType) << endl;
 *oLog << "        stdNum  " << stdNum << endl;
 *oLog << "        fFocLgt    " << fFocLgt << endl;
 *oLog << "        fPlateScaleFactor " <<  fPlateScaleFactor << endl;
 *oLog << "        fAvgTransitTime " << fAvgTransitTime << endl;
 *oLog << "        Primary/Secondary mirror refl.ID " 
       << iPrimReflID << "  /  " << iSecReflID << endl;
 *oLog << "        fRotationOffset " << fRotationOffset << endl;
 *oLog << "        Primary   fDp fDpinner fZp " << fDp << " " 
       << fDpinner << " " << fZp << endl;
 *oLog << "        Primary Translation Offsetsx/y/z " << fPrimXoffset
       << "   " << fPrimYoffset << "  " << fPrimZoffset << endl;
 *oLog << "        Primary rotation Theta Phi offsets " << fPrimThetaOffset 
       << "   " << fPrimPhiAngle << endl;
 *oLog << "        Primary Surface roughness sigma " << fPrimRoughSigma 
       << endl;
 *oLog << "        Secondary fDs fDsinner fZs " << fDs << " "
       << fDsinner << " " << fZs << endl;
 *oLog << "        Secondary Translation Offsetsx/y/z " << fSecondXoffset
       << "   " << fSecondYoffset << "  " << fSecondZoffset << endl;
 *oLog << "        Secondary rotation Theta Phi offsets " << fSecondThetaOffset 
       << "   " << fSecondPhiAngle << endl;
 *oLog << "        Secondary Surface roughness sigma " << fSecondRoughSigma 
       << endl;
 *oLog << "        Focal Surface fk1 fk2 fZf " << fk1 << " " 
       << fk2 << " " << fZf << endl;
 *oLog << "        Focal Surface Translation Offsetsx/y/z " << fFocalPlXoffset
       << "   " << fFocalPlYoffset << "  " << fFocalPlZoffset << endl;
 *oLog << "        Focal Surface rotation Theta Phi offsets " 
       << fFocalPlThetaOffset 
       << "   " << fFocalPlPhiAngle << endl;
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
 *oLog << "        bCameraFlag  " << bCameraFlag << endl;
 *oLog << "        fPixelSize   " << fPixelSize << endl;
 *oLog << "        fPixelPitch  " << fPixelPitch  << endl;
 *oLog << "        fMAPMTWidth   " << fMAPMTWidth << endl;
 *oLog << "        fMAPMTLength   " << fMAPMTLength << endl;
 *oLog << "        fInputWindowThickness   " << fInputWindowThickness << endl;
 *oLog << "        fMAPMTAngularSize   " << fMAPMTAngularSize << endl;
 *oLog << "        fMAPMTOffset        " << fMAPMTOffset << endl;
 *oLog << "        fMAPMTGap        " << fMAPMTGap << endl;
 *oLog << "        fMAPMTRefIndex        " << fMAPMTRefIndex << endl;

};
/************** end of printSCStdOptics ***********************/

//////////////////////////////////////////////////////////////////////////

GSCTelescopeFactory::
GSCTelescopeFactory(GReadSCStd &scReader,
		    const string &editPilotFile) {

  bool debug = false;
  sPilotEdit = editPilotFile;
  if (debug) {
    *oLog << "  -- SCTelescopeFactory Constructor:  " << sPilotEdit << endl;
  }
  iNumSCTelMade = 0;
  readSC = 0;
  pi = 0;
  sPilotEdit = "";
  mGRefl = 0;
  opt = 0;
  SCTel = 0;

  // make the reflectivity map 
  mGRefl = new map<int, TGraph *>;

  readSC = &(scReader);
  readSC->setSCTelescopeFactory(this);

  sPilotEdit = editPilotFile;

  pi = new GPilot(sPilotEdit);
 
  if (debug) {
    *oLog << "    -- end of GSCTelescopeFactory constructor" << endl;
  }
};

/************** end of GSCTelescopeFactory ***********************/

GSCTelescopeFactory::~GSCTelescopeFactory() {
 
  bool debug = false;
  if (debug) {
    *oLog << "  -- GSCTelescopeFactory::~GSCTelescopeFactory" << endl;
  }
  SafeDelete(pi);
  SafeDelete(readSC);

  for (itmStdOp=mStdOptics.begin();
       itmStdOp!=mStdOptics.end(); itmStdOp++) {
    SafeDelete(itmStdOp->second);
  }

  for (itmGRefl=mGRefl->begin();
       itmGRefl!=mGRefl->end(); itmGRefl++) {
    SafeDelete(itmGRefl->second ); 
  }
  SafeDelete(mGRefl);
};
/************** end of ~GSCTelescopeFactory ***********************/

GSCTelescope* GSCTelescopeFactory::makeTelescope(const int &id,
                                                     const int &std) {
  
  int debug = 1;
  if (debug>0) {
    *oLog << " -- GSCTelescopeFactory::makeTelescope" << endl;
    *oLog << "      telID  = " << id << endl;
    *oLog << "      telStd = " << std << endl;
  }

  int idTel = id;
  int iStdID  = std;

  // get parameters for this telescope  
  itmStdOp = mStdOptics.find(iStdID);
  assert(itmStdOp != mStdOptics.end());

  // make the telescope
  GSCTelescope *SCTel1 = 0;
  SCTel1 = new GSCTelescope;

  // opt is pointer to working SCStdOptics structure
  opt = mStdOptics[std];

  // move over all reflection coeffients (the refl.coeff map)
  int primID = opt->iPrimReflID;    // will shortly come from pilot file
  int secondID = opt->iSecReflID;
  //int primID = 1;
  //int secondID = 1;
 
  SCTel1->setReflCv(primID,secondID,mGRefl);

  // get the id and stdOpticsID for the telescope
  SCTel1->setTelID(idTel);
  SCTel1->setStdID(iStdID);

  // set the eTelType and the avg. transition time
  SCTel1->eTelType = opt->stdType;  
  SCTel1->fAvgTransitTime = opt->fAvgTransitTime;
  SCTel1->fRotationOffset = opt->fRotationOffset;

  SCTel1->fFocLgt = opt->fFocLgt;
  SCTel1->fPlateScaleFactor = opt->fPlateScaleFactor;
  SCTel1->fDp = opt->fDp;
  SCTel1-> fDpinner= opt->fDpinner;
  SCTel1->fZp = opt->fZp;

  SCTel1->fPrimXoffset = opt->fPrimXoffset;
  SCTel1->fPrimYoffset = opt->fPrimYoffset;
  SCTel1->fPrimZoffset = opt->fPrimZoffset;
  SCTel1->fPrimThetaOffset = opt->fPrimThetaOffset;
  SCTel1->fPrimPhiAngle = opt->fPrimPhiAngle;
  SCTel1->fPrimRoughSigma = opt->fPrimRoughSigma;
  SCTel1->fPrimRoughMax = opt->fPrimRoughMax;

  SCTel1->fZs = opt->fZs;
  SCTel1->fDs = opt->fDs;
  SCTel1->fDsinner = opt->fDsinner;

  SCTel1->fSecondXoffset = opt->fSecondXoffset;
  SCTel1->fSecondYoffset = opt->fSecondYoffset;
  SCTel1->fSecondZoffset = opt->fSecondZoffset;
  SCTel1->fSecondThetaOffset = opt->fSecondThetaOffset;
  SCTel1->fSecondPhiAngle = opt->fSecondPhiAngle;
  SCTel1->fSecondRoughSigma = opt->fSecondRoughSigma;
  SCTel1->fSecondRoughMax = opt->fSecondRoughMax;

  SCTel1->fZf = opt->fZf;
  SCTel1->fk1 = opt->fk1;
  SCTel1->fk2 = opt->fk2;

  SCTel1->fFocalPlXoffset = opt->fFocalPlXoffset;
  SCTel1->fFocalPlYoffset = opt->fFocalPlYoffset;
  SCTel1->fFocalPlZoffset = opt->fFocalPlZoffset;
  SCTel1->fFocalPlThetaOffset = opt->fFocalPlThetaOffset;
  SCTel1->fFocalPlPhiAngle = opt->fFocalPlPhiAngle;


  SCTel1->fPixelSize = opt->fPixelSize;
  SCTel1->fPixelPitch = opt->fPixelPitch;
  SCTel1->fMAPMTWidth = opt->fMAPMTWidth;
  SCTel1->fMAPMTLength = opt->fMAPMTLength;
  SCTel1->fInputWindowThickness = opt->fInputWindowThickness;
  SCTel1->fMAPMTAngularSize = opt->fMAPMTAngularSize;
  SCTel1->fMAPMTOffset = opt->fMAPMTOffset;
  SCTel1->fMAPMTGap = opt->fMAPMTGap;
  SCTel1->fMAPMTRefIndex = opt->fMAPMTRefIndex;
  SCTel1->bCameraFlag = opt->bCameraFlag;
  SCTel1->iNParP = opt->iNParP;
  SCTel1->fzp = opt->fzp;
  SCTel1->iNParS = opt->iNParS;
  SCTel1->fzs = opt->fzs;

  // now can edit SCTel1
  editWorkingTelescope(SCTel1);
  // can now build the telescope
  SCTel1->buildTelescope();

  iNumSCTelMade++; // increment number of SC telescopes made
  return SCTel1;
};
/************** end of makeTelescope ***********************/

void GSCTelescopeFactory::editWorkingTelescope(GSCTelescope *SCTel1) {
  int iTelID = SCTel1->iTelID;
  (void) iTelID;  // unused
  
  int debug = 0;
  if (debug>0) {
    *oLog << " -- GSCTelescopeFactory::editWorkingTelescope" << endl;
  }

  /*
  int iTelID = SCTel->iTelID;

  string flag = "EDITSCTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    vector<int> listTel;    
    vector<int>::iterator itv;
    GUtilityFuncts::decodeMatlabString(tokens[0],listTel);
    
    // are there edits for this telescope
    itv = find(listTel.begin(),listTel.end(),iTelID);
    if (itv != listTel.end() ) {
      // there are edits for this telescope
      if (tokens[1] == "DUMMY") {
	double tmpfDummy = atof(tokens[2].c_str() );
	int tmpiDummy = atoi(tokens[3].c_str() );
	SCTel->fDummy = tmpfDummy;
	SCTel->iDummy = tmpiDummy;
	
      }  // if tokens = DUMMY
    }  // if edits for this telescope
    
  } // while loop
  */  
};
/************** end of editWorkingTelescope ***********************/

void GSCTelescopeFactory::printStdTelescope(const int &iStd, 
                                               const int &mode,ostream &oStr) {
  (void) mode;  // unused
  oStr << "in not working printStdTelescope" << endl;
  
  // DO NOT USE.
  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GSCTelescopeFactory::printStdTelescope" << endl;
  }

  int stdNum = iStd;

  bool allFlag = false;
  if (iStd < 0) {
    allFlag = true;  // print all the standard telescopes
    stdNum = -iStd;
  }

  
  
};
/************** end of :printStdTelescope ***********************/

void GSCTelescopeFactory::setPrintMode(ostream &oStr,
                                          const int prtMode) {

  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GSCTelescopeFactory::setPrintMode" << endl;
  }
  oPrtStrm = &oStr;
  iPrtMode = prtMode;
}; 
/************** end of setPrintMode  ***********************/




/*
VERSION3.1
2March2015
*/
/*  GDCTelescope.cpp

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

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGW(x) *oLog << "         " << #x << " = " << x << endl

#include "GUtilityFuncts.h"
#include "GDefinition.h"
#include "GRootDCNavigator.h"
#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GGeometryBase.h"
#include "GDCGeometry.h"
#include "GRayTracerBase.h"
#include "GDCRayTracer.h"

DCStdFacet::DCStdFacet() {
  type = 0;
  radius = 0.0;
  curv = 0.0;
  xm = 0.0;
  ym = 0.0;
  zm = 0.0;
  mis_align = 0.0;
  ftprot = 0.0;
  ffprot = 0.0;
  roughness = 0.0;
  reflect = 0.0;
  rflctid = 0;
};
/**************** end of DCStdFacet::DCStdFacet()***********/

DCStdFacet::DCStdFacet(const DCStdFacet &dcf) {

  facNum = dcf.facNum;

  type = dcf.type;     
  sides = dcf.sides;
  radius = dcf.radius;   
  curv = dcf.curv; 
  xm = dcf.xm;         
  ym = dcf.ym;
  zm = dcf.zm;
  vFacLoc = dcf.vFacLoc;
  vFacCentrCurv = dcf.vFacCentrCurv;
  vFacPlLoc = dcf.vFacPlLoc;
  vUnitFacPlToCC = dcf.vUnitFacPlToCC;
  rotTelToFP = dcf.rotTelToFP;
  mis_align = dcf.mis_align;
  ftprot = dcf.ftprot;    
  ffprot = dcf.ffprot;
  roughness = dcf.roughness;
  reflect = dcf.reflect;
  rflctid = dcf.rflctid;
};
/******************* end of DCStdFacet::DCStdFacet **********/

void DCStdFacet::findFacetCurvatureCenter(const double &focLgt) {
 
  bool debug = false;

  if (facNum==7) {
    debug = false;
  }
  else {
    debug = false;
  }

  if (debug) {
    *oLog << "  -- DCStdFacet::findFacetCurvatureCenters()" << endl;
    *oLog << "        facNum " << facNum << endl;
  }

  // make vector from focal pt. to z = foc.length, x=y=0
  ROOT::Math::XYZVector rzfl(0.0,0.0,focLgt);

  if (debug) {
    *oLog << "        vFacLoc ";
    GUtilityFuncts::printGenVector(vFacLoc); *oLog << endl;

    *oLog << "        rzfl ";
    GUtilityFuncts::printGenVector(rzfl); *oLog << endl;
  }

  ROOT::Math::XYZVector tmp = rzfl - vFacLoc;
  vUnitFacPlToCC = tmp.Unit();
      
  if (debug) {
  
    *oLog << "        ivFacLoc ";
    GUtilityFuncts::printGenVector(vFacLoc); *oLog << endl;
    *oLog << "        vector from facetToC ";
    GUtilityFuncts::printGenVector(tmp); *oLog << endl;
    *oLog << "        unit vector from facetToC ";
    GUtilityFuncts::printGenVector(vUnitFacPlToCC); 
    *oLog << endl;
  }
  // add error 
  GUtilityFuncts::addErrorRoot(&vUnitFacPlToCC,mis_align);

  // find center of curvature vector
  vFacCentrCurv = vFacLoc + (curv)*vUnitFacPlToCC;

  if (debug) {
    *oLog << "        adding error mis_align " 
          << mis_align*(TMath::RadToDeg()) << " degrees" << endl;
    *oLog << "        unit vector from facetToC ";
    GUtilityFuncts::printGenVector(vUnitFacPlToCC); *oLog << endl;
    *oLog << "        center of curv. from foc. pt ";
    GUtilityFuncts::printGenVector(vFacCentrCurv);  *oLog << endl;
    *oLog << endl;

  }

};
/******************* end of DCStdFacet::findFacetCurvatureCenters **********/
void DCStdFacet::findFacetPlaneRotationMatrix() {
  
  double tmp1 = 0.0;
  double alpha = 0.0;
  double t[3][3];
  for (int i = 0;i<3;i++) {
    for (int j = 0;j<3;j++) {
      t[i][j] = 0.0;
    }
  }

  double dlfp = vUnitFacPlToCC.X(); 
  double dmfp = vUnitFacPlToCC.Y(); 
  double dnfp = vUnitFacPlToCC.Z(); 

  tmp1 = dmfp/dnfp;
  alpha = 1/sqrt(1 + tmp1*tmp1 );

  t[0][0] = alpha*(dnfp +(tmp1*dmfp));
  t[0][1] = -alpha*tmp1*dlfp;
  t[0][2] = -alpha*dlfp;

  t[0][0] = 0.0;
  t[0][1] = alpha;
  t[0][2] = -alpha*tmp1;
  t[0][0] = dlfp;
  t[0][1] = dmfp;
  t[0][2] = dnfp;

  rotTelToFP.SetComponents(t[0][0],t[0][1],t[0][2],
                           t[1][0],t[1][1],t[1][2],
                           t[2][0],t[2][1],t[2][2]);

};
/******************* end of DCStdFacet::findFacetPlaneRotationMatrix *******/
  
void DCStdFacet::findFacetPlane() {

  //vFacLoc vFacCentrCurv vFacPlLoc vUnitFacPlToC
  // unit vector perpendicular to facet plane
  vUnitFacPlToCC = vFacCentrCurv - vFacLoc;
  vUnitFacPlToCC = vUnitFacPlToCC.Unit();

  // location of center of facet plane
  double tmp = radius / curv;
  double d = curv*( 1 - sqrt(1 - tmp*tmp) );
  vFacPlLoc = vFacLoc + (d*vUnitFacPlToCC);

};
/******************* end of DCStdFacet::findFacetPlane **********/

void DCStdFacet::printDCStdFacet(ostream &oStr) {

  oStr<<  "  -- DCStdFacet::printDCStdFacet " << endl;
  oStr << "      facetNumber (Starting from zero) " << facNum << endl;
  oStr << "      radius curv " << radius << "  " << curv << endl;
  oStr << "      xm  ym  zm " << xm << "  " << ym << "  " 
       << zm << endl;
  oStr<<  "      mis_align ftprot ffprot " 
      << mis_align*(TMath::RadToDeg()) << "  " 
      << ftprot  << "  " << ffprot << endl;
  oStr<<  "      roughness reflect reflctid  " 
      << roughness*(TMath::RadToDeg()) << " "  << reflect 
      << "  " << rflctid << endl;

  oStr << "      vFacLoc ";
  GUtilityFuncts::printGenVector(vFacLoc); oStr << endl;

  oStr << "        vFacCentrCurv  ";
  GUtilityFuncts::printGenVector(vFacCentrCurv); oStr << endl;

  oStr << "        vFacPlLoc  ";
  GUtilityFuncts::printGenVector(vFacPlLoc); oStr << endl;

  oStr << "        vUnitFacPlToCC  ";
  GUtilityFuncts::printGenVector(vUnitFacPlToCC); oStr << endl;
  
  oStr << "   end of DCStdFacet::printDCStdFacet" << endl;

};
/********************** end of DCStdFacet::printDCStdFacet ***********/

GDCTelescope::GDCTelescope() {

  // initialize parameters
  dRotationOffset = 0.0;
  dAvgTransitTime = 0.0;
  dPointOffsetX   = 0.0;
  dPointOffsetY   = 0.0;
  dRadius         = 0.0;
  dFocLgt         = 0.0;
  dCamRad         = 0.0;
  dPlateScaleFactor = 0.0;
  dFocError       = 0.0;

  iTelID          = 0;
  iStdID          = 0;
  eTelType     = DC;
  eGeoType     = DCVERITAS;
  bDoGeoStruct = false;
  geoStruct    = 0;
  //mReflC = 0;
  mVReflWaveLgts = 0;
  mVCoeffs = 0;

  rayTracer = 0;
  //  nbinsx = 0;
  //nbinsy = 0;
  
  bEditAlignFlag    = false;    
  bEditReflectFlag  = false;  
  eRayTracerType = RTDC;

  historyFileName = "";
  historyTreeName = "";
  historyOption   = 0;
  hisF = 0;
  hisT = 0;

  vPhotonCameraLoc.SetCoordinates(0.0,0.0,0.0);
  vPhotonCameraDcos.SetCoordinates(0.0,0.0,0.0);
  bPhotonOnCameraFlag = false;
  fTransitTime = 0.0;
  
  bPhotonHistoryFlag = false;

  oPrtStrm = oLog;
  iPrtMode = 0;
};
/********************** end of GDCTelescope *****************/
GDCTelescope::~GDCTelescope() {

  bool debug = false;

  if (debug) {
    *oLog << "  -- GDCTelescope::~GDCTelescope: iTelID = " 
	  << iTelID << endl;
  }
  SafeDelete(rayTracer);
  /*
    we've already closed the file in writePhotonHistory()
    if (hisF != 0) {
    hisF->cd();  
    hisT->Write();
    SafeDelete(hisF);
  }
  */
};
/********************** end of ~GDCTelescope *****************/
 
void GDCTelescope::injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                                const ROOT::Math::XYZVector &photonDirT,
                                const double &photWaveLgt) {

  bool debug = false;
  if (iTelID==1) debug = false;
  if (debug) {
    *oLog << " -- GDCTelescope::injectPhoton " << endl;
    *oLog << "       photonLocT  ";
    GUtilityFuncts::printGenVector(photonLocT); *oLog << endl;
    *oLog << "       photonDirT  ";
    GUtilityFuncts::printGenVector(photonDirT); *oLog << endl;
  } 

  // initialize photonHistory parameters if necessary
  if (historyFileName!="") {
    initializePhotonHistoryParms();    
    injectXTC = photonLocT.X();
    injectYTC = photonLocT.Y();
    injectZTC = photonLocT.Z();
    injectDcosXTC = photonDirT.X();
    injectDcosYTC = photonDirT.Y();
    injectDcosZTC = photonDirT.Z();
  }

  // initialize rayTracing flags
  bPhotonOnTelescopeFlag = false; // photon hasn't hit telescope yet
  bPhotonOnFacetFlag     = false; // no facet hit yet
  bFacetReflectFlag      = false; // no facet reflection yet
  bPhotonOnCameraFlag    = false; // no camera hit yet

  // set inject location and Dcos's: telescope coordinates
  vPhotonInjectLocRT = photonLocT;
  vPhotonInjectDcosT = photonDirT;
  fPhotWaveLgt = photWaveLgt;

  bPhotonOnCameraFlag = rayTracer->getPhotonLocCamera(&vPhotonCameraLoc,
                                                  &vPhotonCameraDcos,
                                                  &fRayTracerTime);
  //fNetTimeToCamera = fRayTracerTime - dAvgTransitTime;
  fNetTimeToCamera = fRayTracerTime;
  if (debug) {
    *oLog << "       fRayTracerTime    " << fRayTracerTime << endl;
    *oLog << "       dAvgTransitTime  " << dAvgTransitTime << endl;
    *oLog << "       fNetTimeToCamera  " << fNetTimeToCamera << endl;
  }

  if (bPhotonHistoryFlag) {
    fillPhotonHistory();
  }
};
/********************** end of injectPhoton *****************/
bool GDCTelescope::getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                                           ROOT::Math::XYZVector *photonDcos,
                                           double *photonTime) {

  *photonLoc = vPhotonCameraLoc;
  *photonDcos = vPhotonCameraDcos;

  *photonTime = fNetTimeToCamera;

  bool debug = false;
  if ( debug && bPhotonOnCameraFlag) {
    *oLog << "  -- GDCTelescope::getCameraPhotonLocation " << endl;
    *oLog << "        iTelID *photonTime " << iTelID << " " 
	  << *photonTime << " " << fNetTimeToCamera << endl;
  }

  return bPhotonOnCameraFlag;
};
/********************** end of getCameraPhotonLocation *****************/

void GDCTelescope::printTelescope() {

  bool debug = false;
  if (debug) {
    *oPrtStrm << " -- GDCTelescope::printTelescope" << endl;
  }
  if (iPrtMode == 0)  return;

  *oPrtStrm << "       GDCTelescope::printTelescope" << endl;
  *oPrtStrm << "         Telescope Number: " << iTelID << endl;
  *oPrtStrm << "         constructed from standard DC telescope: " << iStdID << endl;
  
  *oPrtStrm<< "         dRotationOffset:  " << dRotationOffset << endl;
  *oPrtStrm<< "         dPointOffsetX/Y:  " << dPointOffsetX << "  " 
	   << dPointOffsetY << endl;
  *oPrtStrm<< "         dRadius:          " << dRadius << endl;
  *oPrtStrm<< "         dFocLgt:           " << dFocLgt << endl;    
  *oPrtStrm<< "         dCamRad:         " << dCamRad << endl;          
  *oPrtStrm<< "         dPlateScaleFactor:   " << dPlateScaleFactor << endl;
  *oPrtStrm<< "         dFocError:         " << dFocError << endl;          
  *oPrtStrm<< "         eTelType:          " << eTelType << endl;
  *oPrtStrm<< "         eRayTracerType:    " << eRayTracerType << endl;
  //*oPrtStrm<< "         nbinsx/y:          " << nbinsx << "   "
  //  << nbinsy << endl;
  *oPrtStrm<< "         eGeoType:          " << eGeoType << endl;
  *oPrtStrm<< "         facet.size()       " << facet.size() << endl;
  
  if (mVReflWaveLgts == 0) {
    *oPrtStrm << "         mVReflWaveLgts pointer zero" << endl;
  } 
  else {
    *oPrtStrm << "         mVReflWaveLgts pointer nonzero" << endl;
    map<int, vector<double> >::iterator itr;
    for (itr = mVReflWaveLgts->begin();
	 itr != mVReflWaveLgts->end(); itr++) {
      *oPrtStrm << "           reflNum numPts " 
		<< (itr->first) << "  "
		<< (itr->second).size() << endl;
    }

  }
  if (mVCoeffs == 0) {
    *oPrtStrm << "       mVCoeffs pointer zero" << endl;
  } 
  else {
    *oPrtStrm << "         mVReflCoeffs pointer nonzero" << endl;
    map<int, vector<double> >::iterator itr;
    for (itr = mVCoeffs->begin();
	 itr != mVCoeffs->end(); itr++) {
      *oPrtStrm << "           reflCoeffNum numPts " 
		<< (itr->first) << "  "
		<< (itr->second).size() << endl;
    }

  }
  
  if (geoStruct == 0) {
    *oPrtStrm << "         geoStruct pointer zero" << endl;
      
  }
  else {
    *oPrtStrm << "         geoStruct pointer nonzero " << endl;
  }
  
  if (rayTracer == 0)  {
    *oPrtStrm << "         rayTracer pointer zero " << endl;
  }
  else {
    *oPrtStrm << "         rayTracer pointer nonzero " << endl;
  }
  *oPrtStrm << "         avgTransitTime " << dAvgTransitTime << endl;
  

  if ( (iPrtMode > 1) && (geoStruct !=0)) {
    geoStruct->printGeometry1(*oPrtStrm,iPrtMode);
  }

  if  (iPrtMode > 2 ) { 
    for (unsigned i = 0;i<facet.size();i++) {
      facet[i].printDCStdFacet(*oPrtStrm);
      *oPrtStrm << endl;
    }
  }
  if (iPrtMode > 3) {
    // print reflection coefficients
    map<int,vector< double > >::iterator itvw;
    map<int,vector< double > >::iterator itvr;
    itvr= mVReflWaveLgts->begin();
    
    for (itvw=mVCoeffs->begin(); itvw !=mVCoeffs->end(); ++itvw) {
      *oPrtStrm << "     printing reflection coefficents curve number: " 
	    << (*itvw).first << "  with " << (*itvw).second.size() 
	    << " points" << endl;
      
      for (unsigned i = 0;i< (*itvw).second.size();i++) {
        double waveln =  (*itvr).second[i];
        double refl =  (*itvw).second[i];
	
        *oPrtStrm << i << "  " << waveln << "  " << refl << endl;
      }    
      itvr++;
    }
  }
  *oPrtStrm << endl;
};
/********************** end of printTelescope *****************/

void GDCTelescope::drawTelescope(const int &option) {
  if (eRayTracerType == RTDCROOT ) {
    rayTracer->printRayTracer();
    *oLog << " ready to draw the telescope" << endl;
  }
  else {
    *oLog << " can not draw telescope without a root raytracer"
	  << endl;
  }
};
/********************** end of drawTelescope *****************/

void GDCTelescope::setPrintMode(ostream &oStr,const int prtMode) {
  bool debug = false;
  if (debug) {
    *oLog << " -- GDCTelescope::setPrintMode" << endl;
  }

  oPrtStrm = &oStr;
  iPrtMode = prtMode;
};
/********************** end of setPrintMode *****************/

void GDCTelescope::setPhotonHistory(const string &rootFile,
                                       const string &treeName,
                                       const int &option) {
  bool debug = false;

  bPhotonHistoryFlag = true;

  historyFileName = rootFile;
  historyTreeName = treeName;
  historyOption = option;
  
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

void GDCTelescope::makePhotonHistoryBranches() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- makePhotonHistoryBranches " << endl;
  }

  hisT->Branch("injectXTC",&injectXTC,"injectXTC/D");
  hisT->Branch("injectYTC",&injectYTC,"injectYTC/D");
  hisT->Branch("injectZTC",&injectZTC,"injectZTC/D");
  hisT->Branch("injectDcosXTC",&injectDcosXTC,"injectDcosXTC/D");
  hisT->Branch("injectDcosYTC",&injectDcosYTC,"injectDcosYTC/D");
  hisT->Branch("injectDcosZTC",&injectDcosZTC,"injectDcosZTC/D");

  hisT->Branch("onTelFlag",&onTelFlag,"onTelFlag/I");
  hisT->Branch("telX",&telX,"telX/D");
  hisT->Branch("telY",&telY,"telY/D");
  hisT->Branch("telZ",&telZ,"telZ/D");

  hisT->Branch("onFacetFlag",&onFacetFlag,"onFacetFlag/I");
  hisT->Branch("facetNum",&facetNum,"facetNum/I");
  hisT->Branch("facetX",&facetX,"facetX/D");
  hisT->Branch("facetY",&facetY,"facetY/D");
  hisT->Branch("facetZ",&facetZ,"facetZ/D");
  hisT->Branch("timeTelFac",&fTimeOnTelToFacet,"timeTelFac/D");
  hisT->Branch("reflectFlag",&reflectFlag,"reflectFlag/I");
  hisT->Branch("onCameraFlag",&onCameraFlag,"onCameraFlag/I");
  hisT->Branch("cameraX",&cameraX,"cameraX/D");
  hisT->Branch("cameraY",&cameraY,"cameraY/D");
  hisT->Branch("cameraZ",&cameraZ,"cameraZ/D");
  hisT->Branch("cameraDcosX",&cameraDcosXTC,"cameraDcosX/D");
  hisT->Branch("cameraDcosY",&cameraDcosYTC,"cameraDcosY/D");
  hisT->Branch("cameraDcosZ",&cameraDcosZTC,"cameraDcosZ/D");
  hisT->Branch("timeFacetToCmra",&fTimeFacetToCamera,"timeFacetToCmra/D");
  hisT->Branch("RayTracerTime",&fRayTracerTime,"RayTracerTime/D");
  hisT->Branch("NetTimeToCamera",&fNetTimeToCamera,"NetTimeToCamera/D");
               
 
};
/************************* end of makePhotonHistoryBranches *****/

void GDCTelescope::fillPhotonHistory() {

  if (bFacetReflectFlag) {
    reflectFlag = 1;
  }
  else {
    reflectFlag = 0;
  }

  if (bPhotonOnCameraFlag) {
    onCameraFlag = 1;
  }
  else {
    onCameraFlag = 0;
  }

  bool debug = false;
  if (debug) {
    *oLog << "  -- fillPhotonHistory " << endl;
  }
  hisF->cd();
  hisT->Fill();
};
/************************* end of fillPhotonHistory *****/

void GDCTelescope::initializePhotonHistoryParms() {

  double defaultF = -9999.9;
  int defaultI    = -9999;
  injectXTC = defaultF;
  injectYTC = defaultF;
  injectZTC = defaultF;

  injectDcosXTC = defaultF;
  injectDcosYTC = defaultF;
  injectDcosZTC = defaultF;

  onTelFlag = 0;
  telX = defaultF;
  telY = defaultF;
  telZ = defaultF;
  fTimeToOnTelescope = defaultF;

  onFacetFlag = 0;
  facetNum = defaultI;
  facetX = defaultF;
  facetY = defaultF;
  facetZ = defaultF;
  reflectFlag = 0;

  cameraDcosXTC = defaultF;
  cameraDcosYTC = defaultF;
  cameraDcosZTC = defaultF;

  onCameraFlag = 0;
  cameraX = defaultF;
  cameraY = defaultF;
  cameraZ = defaultF;

};
/************************* end of initializePhotonHistoryParms *****/

void GDCTelescope::writePhotonHistory() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- in GDCTelescope::writePhotonHistory " << endl;
  }
  // set hisF to zero after deleting so can also do the writing in
  // the destructor.
  if (hisF != 0) {
    hisF->cd();  
    hisT->Write();
    hisF->Close();
    // do I need to delete hisF after Close()? the hisF pointer is not set
    // to zero by its Close() method.
    //SafeDelete(hisF);
  }
};
/************************* end of writePhotonHistory *****/

/*
VERSION3.1
2March2015
*/
/*!  GDCTelescopeFactory.cpp

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

#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"

#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GTelescopeFactory.h"
#include "GDCTelescopeFactory.h"

#include "GGeometryBase.h"
#include "GDCGeometry.h"

#include "GRayTracerBase.h"
#include "GDCRayTracer.h"

#include "GReadDCStdBase.h"
#include "GReadDCStdGrISU.h"
#include "GOrderedGrid.h"
#include "GRootDCNavigator.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       " << #x << " = " << x << endl

DCStdOptics::DCStdOptics() {

  //geoStruct = 0; 

  igeoStd = 0;
  stdType         = DC;
  defoc           = 0.0;
  foclength       = 0.0;
  camRadius       = 0.0;
  plateScaleFactor = 0.0;
  radius          = 0.0;
  rotation_offset = 0.0;
  xoff            = 0.0;
  yoff            = 0.0;
  nbr_mir         = 0;
  nbinsx          = 0;
  nbinsy          = 0;
  gridOption      = 0;
  gridFile        = "";
  grid            = 0;
  rayTracerEnum   = RTDC;
  fAvgTransitTime = 0.0;
  
  cameraRadius = 0.0;

  oStr = oLog;
  iPrtMode = 0;

}
/************** end of DCStdOptics ***********************/

DCStdOptics::~DCStdOptics() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- DCStdOptics::~DCStdOptics " << endl;
  }
  SafeDelete(grid);
  for (unsigned i = 0;i<facet.size();i++) {
    SafeDelete(facet[i]);
  }
}
/************** end of ~DCStdOptics ***********************/

void DCStdOptics::printDCStdOptics() {

  bool debug = false;
  if (debug) {
    *oStr << "  -- DCStdOptics::printDCStdOptics" << endl;
    *oStr << "       iPrtMode " << iPrtMode << endl;
  }
  if (iPrtMode == 0) return;
  *oStr << " -- GDCTelescopeFactory  DCStdOptics::printDCStdOptics" << endl;
  *oStr << "       stdType:        " << getTelType(stdType) << endl;
  *oStr << "       defoc           " << defoc << endl;
  *oStr << "       foclength       " << foclength << endl;
  *oStr << "       camRadius       " << camRadius << endl;
  *oStr << " plateScaleFactor      " << plateScaleFactor << endl;
  *oStr << "       radius          " << radius << endl;
  *oStr << "       rotation_offset " << rotation_offset << endl;
  *oStr << "       xoff            " << xoff << endl;
  *oStr << "       yoff            " << yoff << endl;
  *oStr << "       nbr_mir         " << nbr_mir << endl;
  *oStr << "       nbinsx          " << nbinsx  << endl;
  *oStr << "       nbinsy          " << nbinsy  << endl;
  *oStr << "       gridOption      " << gridOption << endl;
  *oStr << "       gridFile        " << gridFile << endl; 
  *oStr << "       rayTracer       " << getRayTracerType(rayTracerEnum)
	<< endl;
  *oStr << "       facet.size()    " << facet.size() << endl;
  *oStr << "       avgTransitTime  " << fAvgTransitTime <<endl;
  if (grid == 0) {
    *oStr << "      grid not instantiated " << endl;
  }

  /*
  if (mVReflWaveLgts == 0 ) {
    *oStr << "       mVReflWaveLgts pointer zero" << endl;
  }
  else {
    *oStr << "       mVReflWaveLgts.size()   " << mVReflWaveLgts->size() 
         << endl;
  }

  if (mVCoeffs == 0 ) {
    *oStr << "       mVCoeffs pointer zero" << endl;
  }
  else {
    *oStr << "       mVCoeffs.size()   " << mVCoeffs->size() 
         << endl;
  }
  */
  //if (geoStruct == 0) {
  //*oStr << "       geoStruct not instantiated" << endl;
  //}
  //else {
  //*oStr << "       geoStruct instantiated " << endl;
  //}

  //if ( (iPrtMode > 1) && (geoStruct !=0) ) {
  //geoStruct->printGeometry1(*oStr,iPrtMode);   
  //}

  if ( iPrtMode > 2 ) { 
    for (unsigned i=0;i<facet.size(); i++) {
      facet[i]->printDCStdFacet(*oStr);
      *oStr << endl;
    }
  }

  if (iPrtMode > 3 ) {

  }
  *oStr << endl;
};
/************** end of printDCStdOptics ***********************/

//////////////////////////////////////////////////////////////////////////

GDCTelescopeFactory::
GDCTelescopeFactory(GReadDCStdBase &dcReader,const string &editPilotFile) {

  bool debug = false;
  sPilotEdit = editPilotFile;
  if (debug) {
    *oLog << "  -- DCTelescopeFactory Constructor:  " << sPilotEdit << endl;
  }
  opt = 0;
  mVReflWaveLgts = 0;
  mVCoeffs = 0;

  // number of telescopes made by the factory
  iNumTelMake = 0;

  readDC = &(dcReader);
  readDC->setDCTelescopeFactory(this);

  sPilotEdit = editPilotFile;

  pi = new GPilot(sPilotEdit);
 
  if (debug) {
    *oLog << "    -- end of GDCTelescopeFactory constructor" << endl;
  }
};
/************** end of GDCTelescopeFactory ***********************/

GDCTelescopeFactory::~GDCTelescopeFactory() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GDCTelescopeFactory::~GDCTelescopeFactory" << endl;
  }
  SafeDelete(pi);
  SafeDelete(readDC);

  // reflection coeffs owned by the factory so delete here
  SafeDelete(mVReflWaveLgts);
  SafeDelete(mVCoeffs);

  for (itmStdOp=mStdOptics.begin();
       itmStdOp!=mStdOptics.end();
       itmStdOp++) {
    if (itmStdOp->second) SafeDelete(itmStdOp->second);
  }
 
  map<int,GGeometryBase*>::iterator iter;
  for (iter = mDCGeo.begin(); iter != mDCGeo.end();iter++) {
    SafeDelete(iter->second);
  }
 
};
/************** end of ~GDCTelescopeFactory ***********************/

GDCTelescope* GDCTelescopeFactory::makeTelescope(const int &id,
                                                     const int &std) {
  
  int debug = 1;
  if (debug>0) {
    *oLog << " -- GDCTelescopeFactory::makeTelescope" << endl;
    *oLog << "      telID  = " << id << endl;
    *oLog << "      telStd = " << std << endl;
  }

  int idTel = id;
  int iStd  = std;

  // get parameters for this telescope  
  itmStdOp = mStdOptics.find(iStd);
  if (itmStdOp == mStdOptics.end()) {
    cerr << "in GDCTelescopeFactory:  telescope standard " << iStd 
         << " does not exist" << endl;
    cerr << "stopping code " << endl;
    exit(0);
  }

  //opt is pointer to working DCStdOptics
  opt = mStdOptics[std];

    // make the telescope
  GDCTelescope *DCTel;
  DCTel = new GDCTelescope;

  // fill telescope parameters
  DCTel->dAvgTransitTime = opt->fAvgTransitTime;  // from input file
  //*oLog << "dTimeOffset " << DCTel->dTimeOffset << endl;
  //exit(0);
  DCTel->dRotationOffset = opt->rotation_offset;
  DCTel->vRotationOffsetT.SetCoordinates(0.0,0.0,DCTel->dRotationOffset);
  DCTel->dPointOffsetX = opt->xoff;
  DCTel->dPointOffsetY = opt->yoff;
  DCTel->dRadius = opt->radius;
  DCTel->dFocLgt = opt->foclength;
  DCTel->dCamRad = opt->camRadius;
  DCTel->dPlateScaleFactor = opt->plateScaleFactor;
  DCTel->dFocError = opt->defoc;
  DCTel->iTelID = idTel;
  DCTel->iStdID = iStd;
  DCTel->eTelType = opt->stdType;

  //  new geoStruct map in factory
  int geostd = opt->igeoStd;
  DCTel->geoStruct = mDCGeo[geostd];
  DCTel->eGeoType = DCTel->geoStruct->type;

  //DCTel->geoStruct = opt->geoStruct;
  //DCTel->eGeoType = opt->geoStruct->type; 

  if (opt->gridOption > 0) {
    DCTel->bGridOption = true;
  }
  else {
    DCTel->bGridOption = false;
  }

  DCTel->eRayTracerType = opt->rayTracerEnum;

  //DCTel->mReflC = opt->mReflC;
  DCTel->mVReflWaveLgts = mVReflWaveLgts;
  DCTel->mVCoeffs =  mVCoeffs;

  DCTel->facet.resize(opt->nbr_mir);
  for (unsigned i =0;i<opt->facet.size(); i++) {
    DCTel->facet[i] = *(opt->facet[i]);
  }

  // now can edit DCTel
  editWorkingTelescope(DCTel);

  for (unsigned i =0;i<opt->facet.size(); i++) {
    DCTel->facet[i].findFacetCurvatureCenter(DCTel->dFocLgt);

    if ( (DCTel->eRayTracerType == RTDC) ||
         (DCTel->eRayTracerType == RTDCROOT) ) {
      // find location and normal to facet plane
      DCTel->facet[i].findFacetPlane();

      // find rotation matrix from telescope coor to facet coor
      DCTel->facet[i].findFacetPlaneRotationMatrix();
    }

  }

  // can now make the rayTracer, set to zero for the moment
  DCTel->rayTracer = 0;
  makeRayTracer(DCTel,opt);

  //DCTel->findFacetCurvatureCenters();
  for (unsigned i =0;i<opt->facet.size(); i++) {
    DCTel->facet[i].findFacetCurvatureCenter(DCTel->dFocLgt);
  }

  // set reflection coefficient pointers to map with vector
  DCTel->mVReflWaveLgts = mVReflWaveLgts;
  DCTel->mVCoeffs = mVCoeffs;

  iNumTelMake += 1;  // increment number of telescopes make
  return DCTel;
};
/************** end of makeTelescope ***********************/

void GDCTelescopeFactory::editWorkingTelescope(GDCTelescope *DCTel) {
  int debug = 0;
  if (debug>0) {
    *oLog << " -- GDCTelescopeFactory::editWorkingTelescope" << endl;
  }

  int iTelID1 = DCTel->iTelID;

  string flag = "EDITDCTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    vector<int> listTel;    
    vector<int>::iterator itv;
    GUtilityFuncts::decodeMatlabString(tokens[0],listTel);
    
    // are there edits for this telescope
    itv = find(listTel.begin(),listTel.end(),iTelID1);
    if (itv != listTel.end() ) {
      // there are edits for this telescope
      if (tokens[1] == "FACET") {
        vector<int> listFacet;    
        vector<int>::iterator itfv;
        GUtilityFuncts::decodeMatlabString(tokens[2],listFacet);
        for (unsigned i = 0;i<listFacet.size();i++) {
          unsigned facNum = listFacet[i] -1 ;
          if (facNum < DCTel->facet.size() ) {
            if (tokens[3]=="align") {
             DCTel->facet[facNum].mis_align = atof(tokens[4].c_str())
		*(TMath::DegToRad());
	     DCTel->bEditAlignFlag = true;
           }
            else if (tokens[3]=="reflect") {
              DCTel->facet[facNum].roughness = atof(tokens[4].c_str())
		*(TMath::DegToRad());
	      if (tokens.size() > 5 ) {
                DCTel->facet[facNum].reflect   = atof(tokens[5].c_str());
	      }
	      if (tokens.size() > 6 ) {	      
		DCTel->facet[facNum]. rflctid  = atoi(tokens[6].c_str());
	      }
              DCTel->bEditReflectFlag = true;
            }
          }
        }
      }  // if tokens = FACET
    }  // if edits for this telescope
    
  } // while loop
  
};
/************** end of editWorkingTelescope ***********************/

void GDCTelescopeFactory::printStdTelescope(const int &iStd, 
                                               const int &mode,ostream &oStr) {

  // this prints all standard telescopes, not currently used.

  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GDCTelescopeFactory::printStdTelescope" << endl;
  }

  int stdNum = iStd;

  bool allFlag = false;
  if (iStd < 0) {
    allFlag = true;  // print all the standard telescopes
    stdNum = -iStd;
  }

  // get pointer to standard telescope and print details
  map<int,DCStdOptics*>::iterator it;
  map<int,DCStdOptics*>::iterator itm;

  if (!allFlag) {
    // set iterator to iStd telescope if necessary
    it = mStdOptics.find(iStd);

    if (it == mStdOptics.end() ) {
      oStr << "mStdOptics key: " << iStd << " does not exist" << endl;
      oStr << "     returning from GDCTelescopeFactory::printStdTelescope";
      oStr << " without printing" << endl;

      return;
    }
  }
  
  exit(0);
  oStr << "        printing from GDCTelescopeFactory::printStdTelescope" 
       << endl;
  for (itm = mStdOptics.begin(); itm!=mStdOptics.end(); ++itm) {
    if (!allFlag) {
      if (it != itm) continue;
    }
    DCStdOptics *opt1 = (*itm).second;
    
    oStr << "          printing standard telescope number " << (*itm).first 
         << endl;
    opt1->printDCStdOptics();
  }
};
/************** end of :printStdTelescope ***********************/

void GDCTelescopeFactory::setPrintMode(ostream &oStr,
                                          const int prtMode) {

  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GDCTelescopeFactory::setPrintMode" << endl;
  }
  oPrtStrm = &oStr;
  iPrtMode = prtMode;
}; 
/************** end of setPrintMode  ***********************/

void GDCTelescopeFactory::makeRayTracer(GDCTelescope *DCTel,DCStdOptics *opt1) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GDCTelescopeFactory::makeRayTracer " << endl;
    *oLog << "        DCTel->eRayTracerType " << DCTel->eRayTracerType << endl;
  }

  // instantiate raytracer
  if  (DCTel->eRayTracerType == RTDC) {
    DCTel->bDoGeoStruct = false;
  }

  if ( (DCTel->eRayTracerType == RTDC) || 
       (DCTel->eRayTracerType == RTDCROOT ) ) { 

    // determine whether or not to turn off struct root geom
    if  (DCTel->eRayTracerType == RTDC) {
      DCTel->bDoGeoStruct = false;
    }
    else {
      DCTel->bDoGeoStruct = true;
    }
    
    // make a  tmp pointer and use  methods,
    // then store as a GRayTracerBase pointer 
    GDCRayTracer *tmp = new GDCRayTracer(DCTel);
    
    if (DCTel->bGridOption) {
      tmp->setFacetGrid(opt->grid);
    }
    
    DCTel->rayTracer = tmp;
  }
  else if (DCTel->eRayTracerType == RTDCGEO2) {
    *oLog << " there is no Geo2 raytracer " << endl;
    DCTel->rayTracer = 0;
    
  }

};
/************** end of makeRayTracer  ***********************/



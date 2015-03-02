/*
VERSION3.1
2March2015
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
#include "GUtilityFuncts.h"
#include "GRayTracerBase.h"
#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GOrderedGrid.h"
#include "GRootDCNavigator.h"
#include "GGeometryBase.h"
#include "GDCRayTracer.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "      " << #x << " = " << x << endl

GDCRayTracer::GDCRayTracer(GDCTelescope *dcTel) 
  : DCTel(dcTel) {

  // make GRootDCNavigator root ray tracer driver for RTDCROOT 
  // ray tracer. Otherwise, don't need it
  fTimeOnTel = 0.0;
  bOnTelescopeFlag = false; 
  bOnFacetFlag = false;
  bReflFacetFlag = false;
  bOnCameraFlag = false;
  bOnQuadArmFlagIn = false;
  bOnQuadArmFlagOut = false;
  bOnCrossArmFlagIn = false;
  bOnCrossArmFlagOut = false;
  bOnFocusBoxIn = false;
  bOnEdgeBoxFlagIn = false;
  bOnEdgeBoxFlagOut = false;
  bOnShutterIn = false;
  bOnFocusBoxOut = false;
  bOnTopIn = false;
  bOnTopOut = false;
  // default no root ray tracing, no need to create GRootDCNavigator
  geoT = 0;
  bDoGeoStruct = false;

  facetGrid = 0;
  iFacet = 0;
  fTimeFacetToCamera = 0.0;
  fTopVolOrigin2FocBox = 0.0;

  // make the navigator if necessary
  if ( (dcTel->eRayTracerType == RTDCROOT) && 
       ( dcTel->geoStruct->type != NOSTRUCT ) ){
    geoT = new GRootDCNavigator(DCTel);
    // add small amount since geometry step is across boundary
    fTopVolOrigin2FocBox = geoT->getFocalBoxZBottomTopVolCoor() + 0.001;
 
    bDoGeoStruct = true;
  }

};
/*************** end of GDCRayTracer ****************/

GDCRayTracer::~GDCRayTracer() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GDCRayTracer::~GDCRayTracer " << endl;
  }
  SafeDelete(geoT);

};
/*************** end of ~GDCRayTracer ****************/

bool GDCRayTracer::getPhotonOnTelescope() {

  /*  determine if photon strikes telescope sphere; if so
      determine the location RELATIVE TO THE ROTATIONOFFSET
      POINT, NORMALLY THE FOCAL POINT

   */
  DCTel->vPhotonInjectLocTRO = DCTel->vPhotonInjectLocRT -
    DCTel->vRotationOffsetT;

  bool debug = false;
  if (DCTel->iTelID==1) {
    debug = false;
  }
  if (debug) {
    *oLog << "  -- GDCRayTracer::getPhotonOnTelescope " << endl;
    *oLog << "        photon inject relative location ";
    GUtilityFuncts::printGenVector(DCTel->vPhotonInjectLocRT); *oLog << endl;

    *oLog << "        photon inject Dcos ";
    GUtilityFuncts::printGenVector(DCTel->vPhotonInjectDcosT); *oLog << endl;

    *oLog << "        tel.focal length " << DCTel->dFocLgt << endl;

    *oLog << "        vRotationOffsetT ";
    GUtilityFuncts::printGenVector(DCTel->vRotationOffsetT); *oLog << endl;

    *oLog << "        vPhotonInjectLocTRO ";
    GUtilityFuncts::printGenVector(DCTel->vPhotonInjectLocTRO); *oLog << endl;
  }

  bOnTelescopeFlag = false;
  double plDotpUnit = DCTel->vPhotonInjectDcosT.Dot(DCTel->vPhotonInjectLocTRO);
  double rPho2 = DCTel->vPhotonInjectLocTRO.Dot(DCTel->vPhotonInjectLocTRO);
 
  double focLgt2 = DCTel->dFocLgt*(DCTel->dFocLgt);

  double d2 = rPho2 - (plDotpUnit*plDotpUnit);
  double beta2 = focLgt2 - d2;

 if (debug) {
    *oLog << "        plDotpUnit  " << plDotpUnit << endl;
    *oLog << "        rPho2       " << rPho2 << endl;
    *oLog << "        focLgt2     " << focLgt2 << endl;
    *oLog << "        d2          " << d2 << endl;
    *oLog << "        beta2       " << beta2 << endl;
  }
 
  if (beta2 > 0.0) {
    double beta = sqrt(beta2);
    fTimeOnTel = (beta - plDotpUnit);
    DCTel->fTimeToOnTelescope = (beta - plDotpUnit);
    vPhotonOnTelT = DCTel->vPhotonInjectLocTRO +
      fTimeOnTel*(DCTel->vPhotonInjectDcosT);
    vPhotonDirT = DCTel->vPhotonInjectDcosT;
    bOnTelescopeFlag = true;

    if (DCTel->bPhotonHistoryFlag) {
      DCTel->onTelFlag = 1;
      DCTel->telX = vPhotonOnTelT.X();
      DCTel->telY = vPhotonOnTelT.Y();
      DCTel->telZ = vPhotonOnTelT.Z();
    }

    if (debug) {
      *oLog << "            photon on telescope OK " << endl;
      *oLog << "            beta        " << beta << endl;
      *oLog << "            vPhotonOnTelescope ";
      GUtilityFuncts::printGenVector(vPhotonOnTelT);
      *oLog << endl;
      *oLog << "            fTimeOnTel " << fTimeOnTel << endl;
      *oLog << "        bOnTelescopeFlag  " 
	    << bOnTelescopeFlag << endl;
    }
  }
 
  return bOnTelescopeFlag;
};
/************* end of getPhotonOnTelescope *******************/

bool GDCRayTracer::findFacet() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GDCRayTracer::findFacet " << endl;
    *oLog << "       DCTel->bGridOption  " << DCTel->bGridOption
          << endl;
  }
  
  bOnFacetFlag = false;  
  // onTelescope hit location
  double x = vPhotonOnTelT.X();
  double y = vPhotonOnTelT.Y();

  // grid facet parameters, will use if grid flag is true
  list<GridFacet> *gridList = 0;
  list<GridFacet>::iterator gridIter;
  bool facetFlag;

  // maximum number of facets for looping, either all facets, or grid facets
  int maxfac;

  if (DCTel->bGridOption) {
  
    facetFlag = facetGrid->getGridBinList(x,y,gridList);

    if (debug) {
      *oLog << "         facetFlag " << facetFlag << endl;
      if (facetFlag) {
        *oLog << "         number of grid elements " << gridList->size() << endl;
        *oLog << "         facet num ";
        for (gridIter = gridList->begin();
             gridIter != gridList->end(); gridIter++) {
          *oLog << (*gridIter).facetNum << " ";
        }
        *oLog << endl;
      }
    }

  }
  // we now have a group of grid elements to loop over, built in 
  // loop over all facets for testing.

  if (!DCTel->bGridOption) {
    // will loop over all facets
    maxfac = DCTel->facet.size();  // maxfac = number of facets
  }
  else {
    if (facetFlag) {
      maxfac = gridList->size();   // maxfac = number of facets in the gridbin
      gridIter = gridList->begin();
    }
    else {
      //*oLog << "points outside the grid" << endl; // no gridbin available
      return false;  // points are outside the grid
    }
  }

  // loop over selected facets or all facets
  int fNum; 
  for (int k1 = 0;k1 < maxfac;k1++) {

    // set up facet number for this value of k1
    if (!DCTel->bGridOption) {
      fNum = k1;    // loop over all facets
    }
    else {
      // get the k1th grinbin list element
      GridFacet gFacet = *gridIter;
      fNum = gFacet.facetNum;
       if (debug) {
	 DCTel->facet[gFacet.facetNum].printDCStdFacet(*oLog);
        *oLog << "       gFacet " << fNum << "  "
              << gFacet.dist << endl;
        *oLog << "       photonOnT ";
        GUtilityFuncts::printGenVector(vPhotonOnTelT);*oLog << endl;
      }
      gridIter++;
    }

    //DCTel->facet[k1].vFacPlLoc facet plane center
    //vPhotonDirT  photon direction cosines
    //DCTel->facet[k1].vUnitFacPlToCC unit vector
    //vPhotonOnTelT photon telescope impact 
    if (debug) {
      *oLog << "       DCTel->facet[fNum].vFacPlLoc  ";
      *oLog << " fNum " << fNum << endl;
    }

    ROOT::Math::XYZVector vDelta = vPhotonOnTelT 
      - DCTel->facet[fNum].vFacPlLoc;
    double dot1 = vDelta.Dot(DCTel->facet[fNum].vUnitFacPlToCC);
    double dot2 =DCTel->facet[fNum].vUnitFacPlToCC.Dot(vPhotonDirT);
    
    double dd = dot1/dot2;
    
    // get location of photon on the facet plane
    ROOT::Math::XYZVector vPhotonFP = vDelta - dd*vPhotonDirT;
    
    // get distance from facet center and see if too large
    double r = sqrt(vPhotonFP.Dot(vPhotonFP));
    if (debug) {
      *oLog << "       print during find photon loc. on facet plane"
            << endl;
      *oLog << "         facet number k1 facNum " 
            << fNum << "  " << DCTel->facet[fNum].facNum << endl;;
      
      *oLog << "         vPhotonOnTelT  ";
      GUtilityFuncts::printGenVector(vPhotonOnTelT); *oLog << endl;
      
      *oLog << "         vPhotonDirT ";
      GUtilityFuncts::printGenVector(vPhotonDirT); *oLog << endl;
      
      *oLog << "         DCTel->facet[k1].vFacPlLoc ";
      GUtilityFuncts::printGenVector(DCTel->facet[fNum].vFacPlLoc); 
      *oLog << endl;
      
      *oLog << "         vDelta  ";
      GUtilityFuncts::printGenVector(vDelta); *oLog << endl;
      *oLog << "         DCTel->facet[k1].vUnitFacPlToCC  ";
      GUtilityFuncts::printGenVector(DCTel->facet[fNum].vUnitFacPlToCC); 
      *oLog << endl;
      
      *oLog << "         dot1 " << dot1 << endl;
      *oLog << "         dot2 " << dot2 << endl;
      *oLog << "         dd   " << dd << endl;
      
      *oLog << "         vPhotonFP   ";
      GUtilityFuncts::printGenVector(vPhotonFP); *oLog << endl;
      *oLog << "       distance from facet center " << r << endl;
      
    }
    
    if (r < DCTel->facet[fNum].radius) {
      // found a facet within reach
      // but is the location within the facet polygon
      double alpha = DCTel->facet[fNum].ftprot;
      double radius = DCTel->facet[fNum].radius;
      int sides = DCTel->facet[fNum].sides;

      double xfp = vPhotonFP.X();
      double yfp = vPhotonFP.Y();
      if (debug) {
        *oLog << "   Found a facet within reach: number " 
              << fNum << endl;
        DEBUGS(alpha);DEBUGS(radius);DEBUGS(sides);
        DEBUGS(xfp); DEBUGS(yfp);
        *oLog << "   ready to run polyInside " << endl;
        
      }
      
      bool polyTest = GUtilityFuncts::polyInside(sides,alpha,
                                                 radius,xfp,yfp);
      if (debug) {
        *oLog << "    polyTest " << polyTest << endl;

      }

      if (polyTest) {
        // find intersection of photon with facet.
        // same method as finding intersection of photon 
        // with telescope sphere. Use center of curvature
        // as origin

        ROOT::Math::XYZVector vPhotonCC = vPhotonOnTelT - 
          DCTel->facet[fNum].vFacCentrCurv;
        
        if (debug) {
          *oLog << "    vPhotonOnTelT " << endl;
          GUtilityFuncts::printGenVector(vPhotonOnTelT); *oLog << endl;
 
          *oLog << "    vFacCentrCurv ";
          GUtilityFuncts::printGenVector(DCTel->facet[fNum].vFacCentrCurv); 
          *oLog << endl;

          *oLog << "    vPhotonCC ";
         GUtilityFuncts::printGenVector(vPhotonCC); *oLog << endl;

          *oLog << "    vPhotonDirT ";
         GUtilityFuncts::printGenVector(vPhotonDirT); *oLog << endl;
        }
        
        double rLTotal = vPhotonCC.Dot(vPhotonDirT);
        double rp2 = vPhotonCC.Dot(vPhotonCC);
        double d2 = rp2 - (rLTotal*rLTotal);
        double curvR = DCTel->facet[fNum].curv;

        double l2 = curvR*curvR - d2;
        l2 = sqrt(l2);

        double del = rLTotal - l2;  
        ROOT::Math::XYZVector vPhotonOnFacCC = vPhotonCC 
          - del*vPhotonDirT;  // photon on facet wrt center of curv.

        vPhotonOnFacT = vPhotonOnFacCC + DCTel->facet[fNum].vFacCentrCurv;
        fTimeOnFacet = -del;
	DCTel->fTimeOnTelToFacet = -del;

        bOnFacetFlag = true;  
        iFacet = fNum;
        
        if (debug) {
          DEBUGS(rLTotal);DEBUGS(rp2);DEBUGS(d2);
          DEBUGS(curvR);DEBUGS(l2);DEBUGS(del);
          DEBUGS(DCTel->fTimeOnTelToFacet);
          
          *oLog << "    vPhotonOnFacCC " << endl;
          GUtilityFuncts::printGenVector(vPhotonOnFacCC); *oLog << endl;
          *oLog << "    vPhotonOnFacT " << endl;
          GUtilityFuncts::printGenVector(vPhotonOnFacT); *oLog << endl;
        }
        if (DCTel->bPhotonHistoryFlag) {
          DCTel->onFacetFlag = 1;
          DCTel->facetNum = iFacet;
          DCTel->facetX = vPhotonOnFacT.X();
          DCTel->facetY = vPhotonOnFacT.Y();
          DCTel->facetZ = vPhotonOnFacT.Z();
        }
        break;
      }
    }
  }

  return bOnFacetFlag;
};
/******************* end of getPhotonOnFacet ******************/
bool GDCRayTracer::getPhotonReflect() {
 
  // determine if the photon will successfully reflect from the facet
 bool debug = false;
  if (debug) {
    *oLog << " -- GDCRayTracer::getPhotonReflect" << endl;
    *oLog << "     using GUtilityFuncts::photonReflect " << endl;
  }

  DCTel->bFacetReflectFlag = false;

  // set variablesa and parameters
  int reflctid = DCTel->facet[iFacet].rflctid;
  double reflect = DCTel->facet[iFacet].reflect;
  double wavelgt = DCTel->fPhotWaveLgt;

  // use Utility function
  vector<double> vRefWaveLgt = (*DCTel->mVReflWaveLgts)[reflctid];
  vector<double> vRefCoeffs    = (*DCTel->mVCoeffs)[reflctid];
  
  DCTel->bFacetReflectFlag = GUtilityFuncts::photonReflect(vRefCoeffs,
                                                           vRefWaveLgt,
                                                           reflect, 
                                                           wavelgt);
  
  if (debug) {
    *oLog << "     DCTel->bFacetReflectFlag " << DCTel->bFacetReflectFlag 
          << endl;
  }

  return DCTel->bFacetReflectFlag;
  
};
/******************* end of getPhotonReflect ******************/

bool GDCRayTracer::getPhotonOnCamera() {

  // reflect the photon and get its direction cosines
 bool debug = false;
  if (debug) {
    *oLog << " -- GDCRayTracer::getPhotonOnCamera" << endl;
    GUtilityFuncts::printGenVector(vPhotonOnFacT); *oLog << endl;
    *oLog << "        DCTel->vPhotonCameraDcos ";
    GUtilityFuncts::printGenVector(DCTel->vPhotonCameraDcos); *oLog << endl;
    *oLog << "        DCTel->dFocError " << DCTel->dFocError << endl;;    
  }

  // get camera location
  double deFocusing = DCTel->dFocError; // in meters, >0 camera too far

  double ct = (deFocusing - vPhotonOnFacT.Z() ) / (DCTel->vPhotonCameraDcos).Z();

  DCTel->vPhotonCameraLoc = (vPhotonOnFacT + ct*DCTel->vPhotonCameraDcos)*1000.0;
  
  (DCTel->vPhotonCameraLoc).SetZ(deFocusing);

  fTimeFacetToCamera = ct;
  DCTel->fTimeFacetToCamera = fTimeFacetToCamera; 

  if (DCTel->iTelID==1) {
    //*oLog << "DCTel->fTimeToOnTelescope " << DCTel->fTimeToOnTelescope << endl;
    //*oLog << "DCTel->fTimeOnTelToFacet  " << DCTel->fTimeOnTelToFacet  << endl;
    //*oLog << "DCTel->fTimeFacetToCamera " << DCTel->fTimeFacetToCamera << endl;
  }

  DCTel->fRayTracerTime = (DCTel->fTimeToOnTelescope +
			   DCTel->fTimeOnTelToFacet +
			   DCTel->fTimeFacetToCamera ) *
    (1.0E09) / (TMath::C() ) ;
  if (DCTel->iTelID==1) {
    //*oLog << "DCTel->fRayTracerTime " << DCTel->fRayTracerTime << endl;
    //*oLog << endl;
  }

  // just a note, with normal photon incidence, 
  //DCTel->fTransitTime = (fTimeOnTel + fTimeOnFacet + fTimeFacetToCamera) *
  //(1.0E09) / (TMath::C());

  if (debug) {
    *oLog << "        ct               " << ct << endl;
    *oLog << "        DCTel->vPhotonCameraLoc";
    GUtilityFuncts::printGenVector(DCTel->vPhotonCameraLoc); *oLog << endl;
  }

  if (DCTel->bPhotonHistoryFlag) {
    DCTel->cameraDcosXTC = DCTel->vPhotonCameraDcos.X();
    DCTel->cameraDcosYTC = DCTel->vPhotonCameraDcos.Y();
    DCTel->cameraDcosZTC = DCTel->vPhotonCameraDcos.Z();
    DCTel->cameraX = DCTel->vPhotonCameraLoc.X();
    DCTel->cameraY = DCTel->vPhotonCameraLoc.Y();
    DCTel->cameraZ = DCTel->vPhotonCameraLoc.Z();
  }

  // do a camera radius check if necessary
  if (DCTel->dCamRad > 0.1) {
    double xdist = (DCTel->vPhotonCameraLoc).X();
    double ydist = (DCTel->vPhotonCameraLoc).Y();
    double dist = sqrt(xdist*xdist + ydist*ydist);
    if (dist > DCTel->dCamRad) return false;
  }
  return true;
};
/******************* end of getPhotonOnCamera ******************/

void GDCRayTracer::initializeRayTracer() {

  // initialize flags prior to start of ray tracing for any photon
  bOnTelescopeFlag     = false; 
  bOnFacetFlag         = false;
  bReflFacetFlag       = false;
  bOnCameraFlag        = false;
  bOnQuadArmFlagIn     = false;  
  bOnQuadArmFlagOut    = false; 
  bOnCrossArmFlagIn    = false; 
  bOnCrossArmFlagOut   = false;
  bOnFocusBoxIn        = false; 
  bOnEdgeBoxFlagIn     = false;
  bOnEdgeBoxFlagOut    = false;
  bOnShutterIn         = false;      

  bOnFocusBoxOut = false;    
  bOnTopIn       = false;          
  bOnTopOut      = false;         

  // reset time parameters , can eliminate fTimeOnTel later
  fTimeOnTel = 0.0;  // initialize time from injectLoc to TelLoc (meters)
  DCTel->fTimeToOnTelescope = 0.0;
  DCTel->fTimeOnTelToFacet = 0.0;
  DCTel->fTimeFacetToCamera = 0.0;

  DCTel->fRayTracerTime = 0.0;
};
/******************* end of initializeRayTracer ******************/

bool GDCRayTracer::getPhotonLocCamera(ROOT::Math::XYZVector *pCam,
                                          ROOT::Math::XYZVector *pDcos,
                                          double *ptime) {

  // will probably enter results directly in DCTel object, easier.

  bool debug = false;
  if (debug) {
    *oLog << " -- GDCRayTracer::getPhotonLocCamera" << endl;
  }
  // do raytracing

  initializeRayTracer();

  bool returnFlag = false;

  //vPhotonInjectLocTRO
  // getPhotonOnTelescope moves photon to telescope ground coordinates with
  // coor. system origin at the focal point (depending on the rotation offset param.

  double pos[3];
  double dir[3];
  for (int i = 0;i<3;i++) {
    pos[i] = 0.0;
    dir[i] = 0.0;
  }

  if ( (bOnTelescopeFlag = getPhotonOnTelescope()) ) {
    if (bDoGeoStruct) {
      pos[0] = DCTel->vPhotonInjectLocTRO.X();
      pos[1] = DCTel->vPhotonInjectLocTRO.Y();
      pos[2] = DCTel->vPhotonInjectLocTRO.Z();
    
      dir[0] = DCTel->vPhotonInjectDcosT.X();
      dir[1] = DCTel->vPhotonInjectDcosT.Y();
      dir[2] = DCTel->vPhotonInjectDcosT.Z();
      
      // can do tests here with particular values of pos and dir
      /*       
	       pos[0] = 1.5;
	       pos[1] = 0.0;
	       pos[2] = -12.0;
	       dir[0] = 0.0;
	       dir[1] = 0.0;
	       dir[2] = -1.0;
      */
      // place the photon at the top of fTopVol_1 if dir[2] < 0.0 
      // and make current position and direction.
      
      geoT->setTrackingDebug(false);
      geoT->setPositionDirection(pos,dir);
      string sNodeName = geoT->getNextNodeName();
      if (debug) *oLog << "sNodeName In Check" << sNodeName << endl;
      
      if (sNodeName != "fTopVol_1" ){
	if (sNodeName.find("FocBox") ) {
	  bOnFocusBoxIn = true;
	}
	else if (sNodeName.find("Quad") ) {
	  bOnQuadArmFlagIn = true;
	}
	else if (sNodeName.find("Cross") ) {
	  bOnCrossArmFlagIn = true;
	}
	else if (sNodeName.find("EdgeBox") ) {
	  bOnEdgeBoxFlagIn = true;
	}
	else if (sNodeName.find("Shutter") ) {
	  bOnShutterIn = true;
	}
	returnFlag = false;
	return returnFlag;
      }
    }  // bDoGeoStruct
    // continue on, no shadowing for incoming photon

    if ( (bOnFacetFlag = findFacet() ) ) {

      if ( (bReflFacetFlag = getPhotonReflect() ) ) {

	ROOT::Math::XYZVector vnorm = DCTel->facet[iFacet].vFacCentrCurv 
	  - vPhotonOnFacT;
	ROOT::Math::XYZVector vnormUnit = vnorm.Unit();
	ROOT::Math::XYZVector photonToCameraDcos;
	
	GUtilityFuncts::reflectDirection(vnormUnit,vPhotonDirT,
					 &photonToCameraDcos,
					 DCTel->facet[iFacet].roughness);
        
	DCTel->vPhotonCameraDcos =  photonToCameraDcos;
	// see if outgoing photon is shadowed by quad or cross arms
	if (bDoGeoStruct) {
	  pos[0] = vPhotonOnFacT.X();
	  pos[1] = vPhotonOnFacT.Y();
	  pos[2] = vPhotonOnFacT.Z();
	  
	  dir[0] = DCTel->vPhotonCameraDcos.X();
	  dir[1] = DCTel->vPhotonCameraDcos.Y();
	  dir[2] = DCTel->vPhotonCameraDcos.Z();
	  
       
          // comment out for testing 
          //////////////////////////////////////////////////
          // can track location of focus box, for example
          // position of edges of focus box are as specified
          // C. Duke, May 3, 2012
          /*
          pos[0] = 0.4;
          pos[1] = -0.4;
          pos[2] = -2.0;

          dir[0] = dir[1] = 0.0;
          dir[2] = 1.0;
            
          for (int i = 0;i<3;i++) {
            cout << "i  pos dir " << i << "  " << pos[i] << " " 
                 << dir[i] << endl;
          }
    
          /////////////////////////////////////////////////
          for (int j = 0;j<20; j++) {
            pos[0] = pos[0] + 0.01;
            pos[1] = pos[1] - 0.01;
            for (int i = 0;i<3;i++) {
              cout << "i  pos dir " << pos[i] << " " << dir[i] << endl;
            }
            geoT->setPositionDirection(pos,dir);
            string sNodeName3 = geoT->getNextNodeName();
            
            /////////// check for shadowing on the way out////////////
            // get current position after stepping to next node
            double pos4[3],dir4[3];
            geoT->getPositionDirection(pos4,dir4);
            
            *oLog << "       nextNodeName " << j << " " << sNodeName3 << endl;
            for (int i = 0;i<3;i++) {
              *oLog << "     " << i << "  " << pos4[i] << endl;
            }
            *oLog << endl;
         }
          exit(0);
          */
          ////////////////////////////////////////////////////////	  

          // position is a facet hit location
          geoT->setPositionDirection(pos,dir);
	  string sNodeName = geoT->getNextNodeName();
	  if (sNodeName != "fFocBoxVol_1") {
	    return false;
	  }
          // get current position after stepping to next node
          double pos2[3],dir2[3];
          geoT->getPositionDirection(pos2,dir2);
 
          // this removes photons striking side of focus box
          if (pos2[2] > fTopVolOrigin2FocBox ) {
            return false;
          }
    	  
	}  // bDoGeoStruct

	if ( (bOnCameraFlag = getPhotonOnCamera() )) {
	  *pCam = DCTel->vPhotonCameraLoc;
	  *pDcos = DCTel->vPhotonCameraDcos;
	  *ptime = DCTel->fRayTracerTime;
	  returnFlag = true;
        }  // getPhotonOnCamera
      } //  bReflFacetFlag
    } // bOnFacetFlag
  }
  
  if (debug) {
    *oLog << "      flags after ray tracing" << endl;
    *oLog << "        bOnTelescopeFlag " << bOnTelescopeFlag << endl;
    *oLog << "        bOnFacetFlag " << bOnFacetFlag << endl;
    *oLog << "        bReflFacetFlag " << bReflFacetFlag << endl;
    *oLog << "        bOnCameraFlag " << bOnCameraFlag << endl;
    *oLog << "        returnFlag    " << returnFlag << endl << endl;
  }
  return returnFlag;
};
/************* end of getPhotonLocCamera *************/
 
void GDCRayTracer::printRayTracer() {

  if (DCTel->eRayTracerType == RTDCROOT) {
    geoT->drawTelescope();
  }

};
/************* end of printRayTracer ***********************/

void GDCRayTracer::setPrintMode(ostream &oStr,
                                          const int prtMode) {
  bool debug = false;
  if (debug) {
    *oLog << " -- GDCRayTracer::setPrintMode" << endl;
  }
  oPrtStrm = &oStr;
  iPrtMode = prtMode;
}; 
/************* end of setPrintMode ***********************/


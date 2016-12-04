/*
VERSION4.0
30May2016
*/
/*!  testUtilities.cpp
     Charlie Duke
     Grinnell College

This code uses the coordinate transformation functions in GUtilityFuncts.cpp using direction parameters
that you set in the code.  It is useful for testing these functions and as a tutorial on the 
coordinate systems associated with GrOptics.  Understanding the use of these functions is crucial for
understanding the GrOptics code as the location and direction of the Cherenkov photons are passed to the
telescope ray tracing classes.  The coordinate systems are described in a README directory file and also
in the following (and I'll avoid using astronomical descriptions):

Ground-based system: x to the East, y to the North, z to the zenith

Telescope system: on the unit sphere at the pointing location of the telescope, y axis tangent to the sphere 
                  and pointed down. z axis pointed outward, perpendicular to sphere.  x axis perpendicular 
                  to y and z axes where xyz forms a right-handed system.
                  The xy plane forms the "tangent plane" centered at the telescope location on the unit sphere.
                  
Telescope location: By a unit vector from the origin to the telescope's position on the unit sphere.  
                    Location also given with respect to ground coordinates by telescope azimuth and 
                    elevation or zenith angles. The zenith angle: angular distance from the zenith to the telescope.  
                    The azimuth in the ground xy plane from the North (y axis) toward the East (x axis).

Tangent plane: See above definition of the telescope system.  Also, given a point S on the unit sphere and the 
               location of Polaris from the latitude, a second set of tangent plane axes given by: construct
               celestial sphere with Polaris at the top.  The y axis is tangent to the great circle in the
               S to P direction. The z axis points inward. The x axis is the East; a better statement is to 
               choose the x axis such that xyz is a right-handed system.  Now, rotate this sphere through 
	       an angle equal to the complement of the latitiude so that its vertical axis points to the 
	       zenith and you'll see that the two sets of axes for the tangent plane differ only differ by
	       a rotation within the tangent plane. (see the PrintFieldRot function in this file. Determining
	       this angle requires standard spherical trig functions).

Camera coordinates: This gets complicated. You would think that the camera coordinates of the image is just 
                    the reflection of the telescope system coordinates of the photon on the camera, x to -x
		    and y to -y.  However, for historical reasons, the y coordinate on the camera is upward 
		    on the celestial sphere toward the zenith so the y axis has already been reflected.  
		    Thus, letting x go to -x and y remain the same gives the camera photon location 
		    to correctly produce the sky image.
		 
 */


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <deque>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <ctime>
#include <limits>

using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "TRint.h"
#include "TROOT.h"
#include "TRandom3.h"
//#include "TString.h"
#include "TSystem.h"
#include "TVector3.h"

#include "TMatrixD.h"
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/EulerAngles.h>
#include <Math/VectorUtil.h>


#include "GDefinition.h"
#include "GUtilityFuncts.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl

ostream *oLog;
TRandom3 TR3;

void printFieldRot(double az, double zn, double latitude);
void printAzZnToXYcos(double az, double zn);
void printXYcosToAzZn(double xcos, double ycos);

void testAzZnToRotMat(double az, double zn, double x, double y, double z);

void testOffsetXYToAzZn(double offsetX, double offsetY,
                        double az, double zn);

void testtelescopeAzZn(const double &primAz,
                       const double &primZn,
                       const double &wobbleN,
                       const double &wobbleE,
                       const double &telOffsetX,
                       const double &telOffsetY,
                       const double &latitude);

void testtelescopeAzZnRot(const double &primAz,
                       const double &primZn,
                       const double &wobbleN,
                       const double &wobbleE,
                       const double &telOffsetX,
                       const double &telOffsetY,
                       const double &latitude);

//int main(int argc, char *argv[]) {
int main() {

  // TRandom used in GUtilityFuncts
  UInt_t seed = 123456;
  TR3.SetSeed(seed);
  // log file output stream, initialize
  oLog = &cerr;

  *oLog << " ---------------entering testUtilities-------------- " << endl;
  
  if (0) {  
    *oLog << "===================== test  GUtilityFunct::fieldRot ===================" << endl << endl;
    printFieldRot(350.0,  // az
		  0.01,   // zn
		  45.0);  // latitude
   
    *oLog << "=================  end of tests for fieldRot ================" << endl;
  }
  if (0) {
    *oLog << "============ tests for AzZnToXYcos ===============" << endl;
    // note: if zn = 0.0, then xcos = ycos = 0.0 for all az.
    printAzZnToXYcos(90.0,   // az
		     0.0);   // zn

    *oLog << "================ end of tests for AzZnToXYcos ==============" << endl;
  }
  if (0) {
    *oLog << "====================tests fo XYcosToAzZn ======================" << endl;
    printXYcosToAzZn(0.1,   // az
		     0.1);  // zn
    *oLog << "=================== end of tests for XYcosToAzZn =============== " << endl;
  }
  if (1) {
    *oLog << "================= tests for telescopeAzZn ====================" << endl;
    testtelescopeAzZn(30.0,  //azimuth
		      60.0, //zenith
		      0.0,  // wobbleN
		      0.0,  // wobbleE
		      0.0,  // offsetX
		      3.0,  // offsetY
		      45.0); // latitude
    *oLog << "================== end of tests for telescopeAzZn =============" <<endl;
  }
  if(1) {
    *oLog << "================= tests for telescopeAzZnRot ==================" << endl;
    testtelescopeAzZnRot(30.0,  //azimuth
		      60.0, //zenith
		      0.0,  // wobbleN
		      0.0,  // wobbleE
		      0.0,  // offsetX
		      3.0,  // offsetY
		      45.0); // latitude
  }
  
  if (0) {
    *oLog << "=============== tests for GUtilityFunct::AzZnToRotMat ============" << endl;
    // checks out ok. November 2016 C.Duke
    testAzZnToRotMat(90.0,  // az
		     30.0, // zn
		     0.0,  // x  in ground coordinates
  		     1.0,  // y
		     0.0); // z
    *oLog << "=================== end of tests for AzZnToRotMat =============== " << endl;
  }

  if (0) {
 
    *oLog << "================= tests for GUtilityFunct::offsetXYToAzZn " << endl;
    testOffsetXYToAzZn(0.0,   // offsetX
		       0.0,   // offsetY
		       0.0,   // az
		       20.0); // zn
    *oLog << "=================end of tests for offsetXYToAzZn ==================" << endl;
 
  }
};

void printFieldRot(double az, double zn, double latitude){
  
  double azr = az*TMath::DegToRad();
  double znr = zn*TMath::DegToRad();
  double latituder = latitude*TMath::DegToRad();

  double rot = GUtilityFuncts::fieldRot(znr,azr, latituder);
  cout << "az       zn    latitude    fieldRot:  " << endl; 
  cout << az << "       " << zn << "       " << latitude << "          "
       << rot*TMath::RadToDeg() << endl << endl;
};
void printAzZnToXYcos(double az, double zn) {
    double azRad = (TMath::DegToRad())*az;
    double znRad = (TMath::DegToRad())*zn;
    double xcos, ycos;
    *oLog << "      given az " << az << "  and zn = " << zn << endl;
    
    GUtilityFuncts::AzZnToXYcos(azRad,znRad,&xcos,&ycos);
    *oLog << "      xcos and ycos from AzZnToXYcos " << endl;
    *oLog << "         xcos = "<< xcos << "  ycos = " << ycos << endl << endl;
    
    *oLog << "      And returning through XYcosToAzZn " << endl;
    GUtilityFuncts::XYcosToAzZn(xcos,ycos,&azRad,&znRad);
    *oLog << "         az = " << azRad*(TMath::RadToDeg()) << "  zn = "<< znRad*(TMath::RadToDeg()) << endl;
}

void printXYcosToAzZn(double xcos, double ycos) {
    *oLog << "      given xcos " << xcos << "  and ycos = " << ycos << endl;

    double az = 0.0, zn = 0.0;
    GUtilityFuncts::XYcosToAzZn(xcos,ycos,&az,&zn);
    *oLog << "      Az and Zn from XYcosToAzZn " << endl;
    *oLog << "         Az = "<< az*(TMath::RadToDeg()) << "  Zn = " << zn*(TMath::RadToDeg()) << endl << endl;
    
    *oLog << "      And returning through AzZnToXYcos " << endl;
    GUtilityFuncts::AzZnToXYcos(az,zn,&xcos,&ycos);
    *oLog << "         xcos = " << xcos << "  ycos = "<< ycos << endl;
}

void testAzZnToRotMat(double az, double zn, double x, double y, double z) {
  double azr = (TMath::DegToRad())*az;
  double znr = (TMath::DegToRad())*zn;
  
  *oLog << " az, zn " << az << " " << zn << "  " << endl;
  //*oLog << " x, y, z " << x << "   " << y << "   " << z << endl;
  *oLog << "   first check individual rotations, az about z, then zn about x" << endl;
  *oLog << "   then call AzZnToRotMat(az,zn,rotMat) " << endl;
    
  // get rotation matrix
  ROOT::Math::Rotation3D *rotMat = new ROOT::Math::Rotation3D;
  GUtilityFuncts::AzZnToRotMat(azr,znr,rotMat);

  ROOT::Math::XYZVector vVec(x,y,z);
  
  GUtilityFuncts::printGenVector(vVec); *oLog << endl;
  
  *oLog << "       rotated vector using rotation matrix " << endl;
  ROOT::Math::XYZVector vVecRot = (*rotMat)*vVec;
  GUtilityFuncts::printGenVector(vVecRot); *oLog << endl;

  *oLog << "       rotated vector using inverse matrix " << endl;
  ROOT::Math::XYZVector vVecRot2 = (rotMat->Inverse())*vVecRot;
  GUtilityFuncts::printGenVector(vVecRot2);
  *oLog << endl;
};

void testOffsetXYToAzZn(double offsetX, double offsetY,
                        double az, double zn) {
  *oLog << " testOffsetXYToAzZn " << endl;
  //void GUtilityFuncts::offsetXYToAzZn(const double &offsetX, const double &offsetY,
  //                const double &az, const double &zn,
  //                double *azOffset, double *znOffset) {

  double offsetXr = (TMath::DegToRad()) * offsetX; 
  double offsetYr = (TMath::DegToRad()) * offsetY; 
  double azr = (TMath::DegToRad()) * az; 
  double znr = (TMath::DegToRad()) * zn;
    

  double azOffset= 0.0; double znOffset=0.0;
  GUtilityFuncts::offsetXYToAzZn(offsetXr, offsetYr, azr,znr, &azOffset,
                                 &znOffset);
  
  //  *oLog << " offsetX offsetY az  zn azOffset znOffset " << offsetX << "  "
  //    << offsetY << "  " << az << "  " << zn << "     " 
  //    << azOffset << "  " << znOffset << endl;

};

void testtelescopeAzZn(const double &primAz,
                       const double &primZn,
                       const double &wobbleN,
                       const double &wobbleE,
                       const double &telOffsetX,
                       const double &telOffsetY,
                       const double &latitude) {

  double telAzr = 0.0;
  double telZnr = 0.0;
  double sourceXr = 0.0;
  double sourceYr = 0.0;

  double primAzr = primAz * (TMath::DegToRad()); 
  double primZnr = primZn * (TMath::DegToRad()); 
  double wobbleNr = wobbleN * (TMath::DegToRad()); 
  double wobbleEr = wobbleE * (TMath::DegToRad()); 
  double telOffsetXr = telOffsetX * (TMath::DegToRad()); 
  double telOffsetYr = telOffsetY * (TMath::DegToRad()); 
  double latituder = latitude * (TMath::DegToRad()); 

  *oLog << "primAz primZn " << primAz << "  " << primZn << endl;
  *oLog << "t1.1elOffsetX telOffsetY " << telOffsetX << "  " << telOffsetY << endl;
  *oLog << "wobbleN wobbleE " << wobbleN << "  " << wobbleE << endl;
  *oLog << "latitude " << latitude << endl;

  //USE EITHER NEW OR NON-NEW
  GUtilityFuncts::telescopeAzZnNew(primAzr,primZnr,wobbleNr,wobbleEr,
                                telOffsetXr,telOffsetYr,latituder,
                                &telAzr,&telZnr,&sourceXr,&sourceYr);

  DEBUG(telAzr*TMath::RadToDeg()); DEBUG(telZnr*TMath::RadToDeg());
  return;
  //==========================================
  double cameraX = -(sourceXr) * (TMath::RadToDeg());
  double cameraY =  (sourceYr) * (TMath::RadToDeg());

  *oLog << "   telescope az zn " << telAzr* (TMath::RadToDeg())
	<< "   " << telZnr * (TMath::RadToDeg()) << endl;
  
  *oLog << "   source X/Y in telescope coordinates " << (sourceXr) * (TMath::RadToDeg())
	<< "    " << (sourceYr) * (TMath::RadToDeg()) << endl;
  
  *oLog << endl << "   source X/Y on camera " << cameraX << "  " 
        << cameraY << endl << endl;
  double totalSource = sqrt(cameraX*cameraX + cameraY*cameraY); 
  *oLog << "    total source distance " << totalSource << endl;
  double diffZn, diffAz, diffTheta;
  diffZn = telZnr - primZnr;
  diffAz = telAzr - primAzr;
  diffTheta = diffAz*sin(primZnr);
  *oLog << "   diffAz diffZn (tel - prim)   diffTheta " 
        << diffAz * (TMath::RadToDeg()) << "  "
        << diffZn  * (TMath::RadToDeg()) << "       "
        << diffTheta * (TMath::RadToDeg()) << endl;

  if (0) {
    // now get the az/zn of the telescope via rotations on the unit sphere
    // first do the teloffset rotations
    
    // az1, zn1 are the offset az and zn in equatorial coordinates
    *oLog << "=============== all rotations ===============" << endl;
    *oLog << " primAz primZn " << primAzr*TMath::RadToDeg()
	  << "   " << primZnr*TMath::RadToDeg() << endl;
    
    double zn1 = primZnr + telOffsetYr;
    double az1;
    if (primZnr > 0.00001) {
      az1= primAzr + telOffsetXr/sin(primZnr);
    }
    else {
      az1 = primAz;
    }
    *oLog << " -----  make rotations on the sphere ------- " << endl;
    *oLog << "        first make the teloffset adjustments ------ " << endl;
    *oLog << "     az1 zn1 " << az1*(TMath::RadToDeg())
          << " " << zn1*(TMath::RadToDeg()) << endl;

    *oLog << "        now make the wobble adjustments using the az1/zn1 from teloffset " << endl;
    // get unit vector pointing toward telescope on unit sphere
    double xcos1, ycos1;
    GUtilityFuncts::AzZnToXYcos(az1,zn1,&xcos1,&ycos1);
    double zcos1 = sqrt(1 - xcos1*xcos1 - ycos1*ycos1);
    ROOT::Math::RotationX rx(TMath::PiOver2() - latituder);
    ROOT::Math::XYZVector vVec1(xcos1,ycos1,zcos1);
    GUtilityFuncts::printGenVector(vVec1); *oLog << "   vVec1 " << endl;

    // get tel location vector in equatorial coor. system 
    ROOT::Math::XYZVector vVec2 = rx*vVec1;
    GUtilityFuncts::printGenVector(vVec2); *oLog <<  "   vVec2" << endl;

    // first get az and zn within this coor. system
    double az2 = 0.0,zn2 = 0.0;

    GUtilityFuncts::XYcosToAzZn(vVec2.X(), vVec2.Y(),&az2,&zn2);
    *oLog << "     az2 zn2 " << az2*(TMath::RadToDeg())
          << " " << zn2*(TMath::RadToDeg()) << endl;

    ROOT::Math::XYZVector vVec3 = (rx.Inverse()) * vVec2;
    GUtilityFuncts::printGenVector(vVec3); *oLog << "   vVec3  " << endl;

    double az3 = 0.0,zn3 = 0.0;
    GUtilityFuncts::XYcosToAzZn(vVec3.X(), vVec3.Y(),&az3,&zn3);
    *oLog << "  az3  zn3 " <<  az3*(TMath::RadToDeg())
          << " " << zn3*(TMath::RadToDeg()) << endl;
  }
  
  
};
void testtelescopeAzZnRot(const double &primAz,
			  const double &primZn,
			  const double &wobbleN,
			  const double &wobbleE,
			  const double &telOffsetX,
			  const double &telOffsetY,
			  const double &latitude) {
  double telAzr = 0.0;
  double telZnr = 0.0;
  double sourceXr = 0.0;
  double sourceYr = 0.0;

  double primAzr = primAz * (TMath::DegToRad()); 
  double primZnr = primZn * (TMath::DegToRad()); 
  double wobbleNr = wobbleN * (TMath::DegToRad()); 
  double wobbleEr = wobbleE * (TMath::DegToRad()); 
  double telOffsetXr = telOffsetX * (TMath::DegToRad()); 
  double telOffsetYr = telOffsetY * (TMath::DegToRad()); 
  double latituder = latitude * (TMath::DegToRad()); 

  *oLog << "primAz primZn " << primAz << "  " << primZn << endl;
  *oLog << "t1.1elOffsetX telOffsetY " << telOffsetX << "  " << telOffsetY << endl;
  *oLog << "wobbleN wobbleE " << wobbleN << "  " << wobbleE << endl;
  *oLog << "latitude " << latitude << endl;
  
  GUtilityFuncts::telescopeAzZnRot(primAzr,primZnr,wobbleNr,wobbleEr,
                                telOffsetXr,telOffsetYr,latituder,
                                &telAzr,&telZnr,&sourceXr,&sourceYr);

  DEBUG(telAzr*TMath::RadToDeg()); DEBUG(telZnr*TMath::RadToDeg());
  return;
   
  
};



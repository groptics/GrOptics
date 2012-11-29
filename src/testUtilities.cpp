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
#include "TRint.h"
#include "TROOT.h"
#include "TRandom3.h"
//#include "TString.h"
#include "TSystem.h"
#include "TVector3.h"

#include "GDefinition.h"
#include "GUtilityFuncts.h"

ostream *oLog;
TRandom3 TR3;

void printFieldRot(double az, double zn, double latitude); 
void printWobbleToAzZn(double wobbleN, double wobbleE,double latitude,
                       double az, double zn);
void testAzZnToRotMat(double az, double zn, double latitude, double wobbleN,
                      double wobbleE);

void testOffsetXYToAzZn(double offsetX, double offsetY,
                        double az, double zn);

void testtelescopeAzZn(const double &primAz,
                       const double &primZn,
                       const double &wobbleN,
                       const double &wobbleE,
                       const double &telOffsetX,
                       const double &telOffsetY,
                       const double &latitude);

int main(int argc, char *argv[]) {
  double xcos = 0.0;
  double ycos = 0.0;
  double zcos = 0.0;
  double az = 0.0;
  double zn = 0.0;

  // TRandom used in GUtilityFuncts
  UInt_t seed = 123456;
  TR3.SetSeed(seed);
 // log file output stream, initialize
  oLog = &cerr;

  *oLog << " entering testUtilities " << endl;

  /*
  *oLog << endl 
        << " =============== tests of transformation functions ============" 
        << endl << endl;
  *oLog << " -- GUtilityFunct::XYcosToAzZn " << endl;
  *oLog << "      Expected results " << endl;
  *oLog << "            xcos = ycos = 0.0      should give Az = Zn = 0.0 " 
        << endl;
  GUtilityFuncts::XYcosToAzZn(xcos,ycos,&az,&zn);
  *oLog << "      calculated az zn " << az*(TMath::RadToDeg())
        << " " << zn*(TMath::RadToDeg()) << endl << endl;

  xcos = 0.3; 
  ycos = 0.2;
  zcos = sqrt(1 - xcos*xcos - ycos*ycos);
  *oLog << "      xcos  ycos zcos " << xcos << "  " << ycos 
        << "  " << zcos << endl;
  *oLog << "      expected results for Az Zn: 56.3099 21.1343 " << endl;
  GUtilityFuncts::XYcosToAzZn(xcos,ycos,&az,&zn);
  *oLog << "      calculated az zn " << az*(TMath::RadToDeg())
        << " " << zn*(TMath::RadToDeg()) << endl << endl;
 
  double az1 = az;  
  double zn1 = zn;
 
  xcos = -0.3;
  ycos = 0.2;
  zcos = sqrt(1 - xcos*xcos - ycos*ycos);
  *oLog << "      xcos  ycos zcos " << xcos << "  " << ycos 
        << "  " << zcos << endl;
  
  GUtilityFuncts::XYcosToAzZn(xcos,ycos,&az,&zn);
  *oLog << "      calculated az zn " << az*(TMath::RadToDeg())
        << " " << zn*(TMath::RadToDeg()) << endl << endl;
  *oLog << "      Note: azimuth range 0 to 2*pi " << endl;

  double az2 = az;
  double zn2 = zn;

  *oLog << "----------------------------------------------" << endl;
  *oLog << "  -- GUtilityFunct::AZZnToXYcos " << endl;
  *oLog << "      take values for az,zn from test for XYcosToAzZn above "
        << " and reverse " << endl;
  az = 0.0;
  zn = 0.0;
  GUtilityFuncts::AzZnToXYcos(az,zn,&xcos,&ycos);
  cout << "      az, zn: xcos ycos " << az << "  " << zn << ":     " 
       << xcos << "  " << ycos << endl << endl;

  GUtilityFuncts::AzZnToXYcos(az1,zn1,&xcos,&ycos);
  cout << "      az, zn: xcos ycos " << az1*(TMath::RadToDeg()) 
       << "  " << zn1*(TMath::RadToDeg()) << ":     " 
       << xcos << "  " << ycos << endl << endl;
  
  GUtilityFuncts::AzZnToXYcos(az2,zn2,&xcos,&ycos);
  cout << "      az, zn: xcos ycos " << az2*(TMath::RadToDeg()) 
       << "  " << zn2*(TMath::RadToDeg()) << ":     " 
       << xcos << "  " << ycos << endl << endl;

  *oLog << "----------------------------------------------" << endl;
  *oLog << "  -- GUtilityFunct::fieldRot " << endl;
  printFieldRot(0.0, 30.0, 90.0);
  printFieldRot(45.0,90.0, 90.0 );
  printFieldRot(45.0,90.0, 0.0 );
  printFieldRot(0.0,30.0, 45.0 );
  printFieldRot(0.0,60.0, 45.0 );

  printFieldRot(30.0,60.0, 45.0 );
  printFieldRot(30.0,50.0, 45.0 );
  printFieldRot(30.0,40.0, 45.0 );
  printFieldRot(30.0,30.0, 45.0 );

  printFieldRot(135.0,60.0, 45.0 );
  printFieldRot(135.0,50.0, 45.0 );
  printFieldRot(135.0,40.0, 45.0 );
  printFieldRot(135.0,30.0, 45.0 );

  printFieldRot(180.0,60.0, 45.0 );

  printFieldRot(225.0,60.0, 45.0 );
  printFieldRot(225.0,50.0, 45.0 );
  printFieldRot(225.0,40.0, 45.0 );
  printFieldRot(225.0,30.0, 45.0 );

  printFieldRot(315.0,60.0, 45.0 );
  printFieldRot(315.0,50.0, 45.0 );
  printFieldRot(315.0,40.0, 45.0 );
  printFieldRot(315.0,30.0, 45.0 );

  printFieldRot(350.0,60.0, 45.0 );
  printFieldRot(350.0,50.0, 45.0 );
  printFieldRot(350.0,40.0, 45.0 );
  printFieldRot(350.0,30.0, 45.0 );

  printFieldRot(359.0,60.0, 45.0 );
  printFieldRot(359.0,50.0, 45.0 );
  printFieldRot(359.0,40.0, 45.0 );
  printFieldRot(359.0,30.0, 45.0 );

  printFieldRot(0.0,60.0, 45.0 );
  printFieldRot(0.0,50.0, 45.0 );
  printFieldRot(0.0,40.0, 45.0 );
  printFieldRot(0.0,30.0, 45.0 );

  *oLog << "----------------------------------------------" << endl;
  *oLog << "  -- GUtilityFunct::wobbleToAzZn " << endl;
  printWobbleToAzZn(0.0,0.0,90.0,0.0,30.0);
  printWobbleToAzZn(0.5,0.0,90.0,0.0,30.0);
  printWobbleToAzZn(0.0,0.5,90.0,0.0,30.0);
  printWobbleToAzZn(0.0,0.5,90.0,0.0,60.0);
  printWobbleToAzZn(0.0,0.5,90.0,30.0,30.0);
  printWobbleToAzZn(0.0,0.5,90.0,30.0,60.0);

  printWobbleToAzZn(0.5,0.0,45.0,0.0,30.0);
  printWobbleToAzZn(0.5,0.0,45.0,0.0,60.0);

  //printWobbleToAzZn(double wobbleN, double wobbleE,double latitude,
  //                     double az, double zn)

  *oLog << "----------------------------------------------" << endl;
  *oLog << "  -- GUtilityFunct::testAzZnToRotMat " << endl;

  testAzZnToRotMat(0.0,20.0,90.0, 0.5, 0.0);

  //testAzZnToRotMat(double az, double zn, double latitude, double wobbleN,
  //                  double wobbleE)
  *oLog << "----------------------------------------------" << endl;
  *oLog << "  -- GUtilityFunct::testOffsetXYToAzZn " << endl;
  testOffsetXYToAzZn(0.0,1.0,0.0,20.0);
  //void testOffsetXYToAzZn(double offsetX, double offsetY,
  //                    double az, double zn)
  */
  *oLog << "----------------------------------------------" << endl;
  *oLog << "  -- testtelescopeAzZn " << endl;
  testtelescopeAzZn(0.0,  //azimuth
                    30.0, //zenith
                    0.0,  // wobbleN
                    -0.5,  // wobbleE
                    0.0,  // offsetX
                    0.0,  // offsetY
                    90.0); // latitude


};

void printFieldRot(double az, double zn, double latitude){
  
  double azr = az*TMath::DegToRad();
  double znr = zn*TMath::DegToRad();
  double latituder = latitude*TMath::DegToRad();

  double rot = GUtilityFuncts::fieldRot(znr,azr, latituder);
  cout << " az, zn, latitude, fieldRot:  " 
       << az << "  " 
       << zn << "  "
       << latitude << "  "
       << rot*TMath::RadToDeg() << endl;
};

void printWobbleToAzZn(double wobbleN, double wobbleE,double latitude,
                       double az, double zn) {

  double wobbleNr = wobbleN*TMath::DegToRad();
  double wobbleEr = wobbleE*TMath::DegToRad(); 
  double latituder = latitude*TMath::DegToRad();
  double azr = az*TMath::DegToRad();
  double znr = zn*TMath::DegToRad();

  double azNew = 0.0; double znNew = 0.0;

  GUtilityFuncts::wobbleToAzZn(wobbleNr,wobbleEr,latituder,
                               azr,znr,&azNew,&znNew);
                               
  cout << "  wobbleN wobbleE latitude az zn azNew znNew "
       << wobbleN << "  " << wobbleE << "  " << latitude << "  " 
       << az << "  " << zn << "     " << azNew*TMath::RadToDeg()
       << "   " << znNew*TMath::RadToDeg() << endl;
  ////////////////////////////////////////////////

  
};
void testAzZnToRotMat(double az, double zn, double latitude, double wobbleN,
                      double wobbleE) {
  double azr = (TMath::DegToRad())*az;
  double znr = (TMath::DegToRad())*zn;
  double latituder = (TMath::DegToRad())*latitude;
  double wobbleNr = (TMath::DegToRad())*wobbleN;
  double wobbleEr = (TMath::DegToRad())*wobbleE;

  *oLog << " az,zn, latitude, wobbleN, wobbleS: " << az << " " << zn
        << "  " << latitude << " " << wobbleN << " " << wobbleE << endl;

  // test az to rotation matrix first
  ROOT::Math::Rotation3D *rotPrim = new ROOT::Math::Rotation3D;
  GUtilityFuncts::AzZnToRotMat(azr,znr,rotPrim);

  double dl = 0.0;
  double dm = sin(20.0*(TMath::DegToRad()) );
  double dn = sqrt(1 - dl*dl - dm*dm);

  //ROOT::Math::XYZVector vVec(dl,dm,dn);
  ROOT::Math::XYZVector vVec(0.0,0.0,1.0);
  GUtilityFuncts::printGenVector(vVec); *oLog << endl;
  ROOT::Math::XYZVector vVecRot = (rotPrim->Inverse())*vVec;
  GUtilityFuncts::printGenVector(vVecRot); *oLog << endl;
                 

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

  GUtilityFuncts::telescopeAzZn(primAzr,primZnr,wobbleNr,wobbleEr,
                                telOffsetXr,telOffsetYr,latituder,
                                &telAzr,&telZnr,&sourceXr,&sourceYr);

  double cameraX = -(sourceXr) * (TMath::RadToDeg());
  double cameraY =  (sourceYr) * (TMath::RadToDeg());
  *oLog << endl << "   sourceX/Y on camera " << cameraX << "  " 
        << cameraY << endl << endl;
  double totalSource = sqrt(cameraX*cameraX + cameraY*cameraY); 
  *oLog << "    total source distance " << totalSource << endl;
  double diffZn, diffAz, diffTheta;
  diffZn = telZnr - primZnr;
  diffAz = telAzr - primAzr;
  diffTheta = diffAz*sin(primZnr);
  *oLog << "   diffAz diffZn   diffTheta " 
        << diffAz * (TMath::RadToDeg()) << "  "
        << diffZn  * (TMath::RadToDeg()) << "       "
        << diffTheta * (TMath::RadToDeg()) << endl;

  

};



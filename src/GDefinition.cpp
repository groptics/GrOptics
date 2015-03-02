/*
VERSION3.1
2March2015
*/
/*!  GDefinition.cpp

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
#include "TRandom3.h"

using namespace std;

#include "GDefinition.h"

string getRdType(const enum RdType &rdType) {

  string rdtypeStr("");

  if (rdType==GRISU) {
    rdtypeStr = "GRISU";
  }
  else if (rdType==CORSIKA) {
    rdtypeStr = "CORSIKA";
  }

  return rdtypeStr;

};
enum RdType getRdTypeEnum(const string &rdTypeStr) {

  enum RdType rdType;

  string teststr = rdTypeStr;

  transform(teststr.begin(),teststr.end(),
            teststr.begin(),::toupper);
  
  if (teststr=="GRISU") {
    rdType = GRISU;
  }
  else if (teststr=="CORSIKA") {
    rdType = CORSIKA;
  }
  else {
    *oLog << " in getRdTypeEnum " << endl;
    *oLog << "     string for RdType enum not defined: ";
    *oLog << teststr << endl;
    *oLog << "         STOPPING CODE " << endl;
    exit(0);
  }

  return rdType;
};
/************* end of getRdType **********************/

string getOfType(const enum OfType &rdType) {
  string rdtypeStr("");

  if (rdType==ASCI) {
    rdtypeStr = "ASCI";
  }
  else if (rdType==ROOTLOC) {
    rdtypeStr = "ROOTLOC";
  }
  else if (rdType==ROOTPIX) {
    rdtypeStr = "ROOTPIX";
  }

  return rdtypeStr;

};
/************* end of getOfType **********************/

enum OfType getOfTypeEnum(const string &ofTypeStr) {

  string tmp = ofTypeStr;
  enum OfType tenum;

  string teststr = ofTypeStr;

  transform(teststr.begin(),teststr.end(),
            teststr.begin(),::toupper);
  
  if (teststr=="ASCI") {
    tenum = ASCI;
  }
  else if (teststr=="ROOTLOC") {
    tenum = ROOTLOC;
  }
  else if (teststr=="ROOTPIX") {
    tenum = ROOTPIX;
  }
  else {
    *oLog << " in getOfTypeEnum " << endl;
    *oLog << "     string for OfType enum not defined: ";
    *oLog << teststr << endl;
    *oLog << "         STOPPING CODE " << endl;
    exit(0);
  }

  return tenum;

}
/************* end of getOfType  **********************/

string getTelType(const enum TelType &rdType) {
    string rdtypeStr("");

  if (rdType==DC) {
    rdtypeStr = "DC";
  }
  else if (rdType==SC) {
    rdtypeStr = "SC";
  }
  else if (rdType==SEGSC) {
    rdtypeStr = "SEGSC";
  }

  return rdtypeStr;

};
/************* end of getTelType **********************/

enum TelType getTelTypeEnum(const string &telTypeStr) {


  enum TelType tenum;
  string teststr = telTypeStr;

  transform(teststr.begin(),teststr.end(),
            teststr.begin(),::toupper);
  
  if (teststr=="DC") {
    tenum = DC;
  }
  else if (teststr=="SC") {
    tenum = SC;
  }
  else if (teststr=="SEGSC") {
    tenum = SEGSC;
  }
  else {
    *oLog << " in getTelTypeEnum " << endl;
    *oLog << "     string for TelType enum not defined: ";
    *oLog << teststr << endl;
    *oLog << "         STOPPING CODE " << endl;
    exit(0);
  }

  return tenum;

}
/************* end of getTelTypeEnum **********************/

string getFacetShape(const enum FacetShape &rdType) {
    string rdtypeStr("");

  if (rdType==CIR) {
    rdtypeStr = "CIR";
  }
  else if (rdType==HEX) {
    rdtypeStr = "HEX";
  }
  else if (rdType==SQR) {
    rdtypeStr = "SQR";
  }

  return rdtypeStr;

};
/************* end of getFacetShape **********************/

enum FacetShape getFacetShapeEnum(const string &facetShapeStr) {

  enum FacetShape tenum;
  string teststr = facetShapeStr;

  transform(teststr.begin(),teststr.end(),
            teststr.begin(),::toupper);
  
  if (teststr=="CIR") {
    tenum = CIR;
  }
  else if (teststr=="HEX") {
    tenum = HEX;
  }
  else if (teststr=="SQR") {
    tenum = SQR;
  }
  else {
    *oLog << " in getFacetShapeEnum " << endl;
    *oLog << "     string for FacetShape enum not defined: ";
    *oLog << teststr << endl;
    *oLog << "         STOPPING CODE " << endl;
    exit(0);
  }

  return tenum;

}
/************* end of getFacetShapeEnum **********************/

string getRayTracerType(const enum RayTracerType &rdType) {
  string rdtypeStr("");
  
  if (rdType==RTDC) {
    rdtypeStr = "RTDC";
  }
  else if (rdType==RTDCROOT) {
    rdtypeStr = "RTDCROOT";
  }
  else if (rdType==RTDCGEO2) {
    rdtypeStr = "RTDCGEO2";
  }

  return rdtypeStr;
  
};
/************* end of getRayTracerType **********************/

enum GeoType getGeoTypeEnum(const string &geoTypeStr) {

  string teststr = geoTypeStr;
  enum GeoType tenum;

  transform(teststr.begin(),teststr.end(),
            teststr.begin(),::toupper);
  
  if (teststr=="DCVERITAS") {
    tenum = DCVERITAS ;
  }
  else if (teststr=="DCxxx") {
    tenum = DCxxx;
  }
  else if (teststr=="SCxxx") {
    tenum = SCxxx;
  }
  else if (teststr=="NOSTRUCT") {
    tenum = NOSTRUCT;
  }
  else {
    *oLog << " in getGeoTypeEnum " << endl;
    *oLog << "     string for GeoType enum not defined: ";
    *oLog << teststr << endl;
    *oLog << "         STOPPING CODE " << endl;
    exit(0);
  }

  return tenum;

}
/************* end of getGeoTypeEnum **********************/

string getGeoType(const enum GeoType &geoType) {

  string geoTypeStr("");
  if (geoType==DCVERITAS) {
    geoTypeStr = "DCVERITAS";
  }
  else if (geoType==DCxxx) {
    geoTypeStr = "DCxxx";
  }
  else if (geoType==SCxxx) {
    geoTypeStr = "SCxxx";
  }
  else if (geoType==NOSTRUCT) {
    geoTypeStr = "NOSTRUCT";
  }
  
  return geoTypeStr;
};
/************* end of getGeoType **********************/

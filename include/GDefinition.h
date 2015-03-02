/*
VERSION3.1
2March2015
*/
/*! \brief GDefinition Functions: general definitions, enums and associated
  functions to return enum strings.      

 */

#ifndef GDEFINITION
#define GDEFINITION

// version number
const string VERSN("3.0");

// log file
extern ostream *oLog;

// may want to put this in a namespace later

/*! telescope type enum
    DC == Davis Cotton, SC == Swartzshield/Coudee
*/
enum TelType {DC,SC,SEGSC};

/*! facet shape enum
     CIR == circle, HEX == hexagon, SQR == square
*/
enum FacetShape {CIR=1,HEX,SQR};

/*! ray tracing method enum
   RTDC == manual(no ROOT), RTDCROOT == manual+Geo struct, 
   RTDCGEO2 == geo for all
*/
enum RayTracerType {RTDC, RTDCROOT, RTDCGEO2}; 

/*!  telescope structure enum (more than one DC type possible, etc)
   DCVERITAS == veritas quad-arm,cross-arm, focus box type structure
   DCxxx and SCxxx just place holders (not yet used)
*/
enum GeoType {DCVERITAS,DCxxx,SCxxx,NOSTRUCT};

/*!  reader type enum
     GRISU == GrISU file format (from cherenkf7 or from corsikaIOreader
     CORSIKA not yet used
 */
enum RdType {GRISU,CORSIKA};

/*!// output file type enum, ASCI, ROOTLOC, ROOTPIX
 */
enum OfType {ASCI, ROOTLOC,ROOTPIX};

enum MirSeg {P1,P2,S1,S2};

/*!// record type for GrISU input files,S,P,R,H,eof lines
 */
// CD:2Mar2015 added C record CREC
enum GrISURecType {SREC,PREC,RREC,HREC,CREC,EOFREC};

/*!// type for rayplots 
 */
enum RayPlotType {FOCUSONLY, ALLSURFACES  };

// set up common random number generator for all classes, move the 
// include to .cpp files later.  maybe have a separate include file for 
// globals, GGlobal.h, fix this later.

// global variables
#include "TRandom3.h"  // fix this up later
extern TRandom3 TR3; //!< random number global class
extern string Versn;  //!< version number, set in main

/*! returns string with name of RdType enum parameter
    returns "" if no match
 */
string getRdType(const enum RdType &rdType);

/*! returns string with name of RdType enum parameter
    returns "" if no match
 */
enum RdType getRdTypeEnum(const string &rdTypeStr);

/*! returns string with name of RdType enum parameter
    returns "" if no match
 */
string getOfType(const enum OfType &rdType);

/*! returns string with name of OfType enum parameter
    returns "" if no match
 */
enum OfType getOfTypeEnum(const string &ofTypeStr);

/*! returns string with name of TelType enum parameter
    returns "" if no match
 */
string getTelType(const enum TelType &rdType);

/*! returns string with name of OfType enum parameter
    returns "" if no match
 */
enum TelType getTelTypeEnum(const string &telTypeStr);

/*! returns string with name of FacetShape enum parameter
    returns "" if no match
 */
string getFacetShape(const enum FacetShape &rdType);

/*! returns string with name of FacetShape enum parameter
    returns "" if no match
 */
enum FacetShape getFacetShapeEnum(const string &facetShapeStr);

/*! returns string with name of RayTracerType enum parameter
    returns "" if no match
 */
string getRayTracerType(const enum RayTracerType &rdType);

/*! returns string with name of FacetShape enum parameter
    returns "" if no match
 */
enum RayTracerType getRayTracerTypeEnum(const string &rayTracerTypeStr);

/*! returns string with name of Type enum parameter
    returns "" if no match
 */
string getGeoType(const enum GeoType &geoType);

/*! returns string with name of FacetShape enum parameter
    returns "" if no match
 */
enum GeoType getGeoTypeEnum(const string &geoTypeStr);

extern std::string TESTSR;

struct mirrorSegmentDetails {
  Double_t rmin;
  Double_t rmax;
  Double_t margin;
  Double_t delPhi;
  Int_t reflect;
  Double_t posErrorX;
  Double_t posErrorY;
  Double_t posErrorZ;
  Double_t rotErrorPhi;
  Double_t rotErrorTheta;
  Double_t rotErrorPsi;
  Double_t roughness;
  Int_t bRead; // if 0, set from BASIC; if 1, set from CHANGE
};
#endif

/*
VERSION3.1
2March2015
*/
/*!  GReadPhotonGrISU.cpp
     Charlie Duke
     Grinnell College
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
#include "TMatrixD.h"
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/EulerAngles.h>
#include <Math/VectorUtil.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "GUtilityFuncts.h"
#include "GDefinition.h"

#include "GReadPhotonBase.h"
#include "GReadPhotonGrISU.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "      " << #x << " = " << x << endl

GReadPhotonGrISU::GReadPhotonGrISU() {
  bool debugG = false;
  if (debugG) {
    *oLog << "  -- GReadPhotonGrISU::GReadPhotonGrISU" << endl;
  }
  pInStream = 0;
  sInFileStr = "";
  sInFileHeader = "";
  fObsHgt = 0.0;
  fGlobalEffic = 0.0;
  iParticleType = 0;
  sFileLine = "";
  eRecType = SREC;
  fSEnergy = 0.0;
  fSAz = 0.0;
  fSZn = 0.0;
  for (int i = 0;i<3;i++) {
    iSSeed[i] = 0;
  }
  fPHgtEmiss = 0.0;
  fPTime = 0.0;
  fPWaveLgt = 0.0;
  iPType = 0;
  iPTel = 0;
  fPAz = 0.0;
  fPZn = 0.0;
  fWobbleE = 0.0;
  fWobbleN = 0.0;
  fWobbleR = 0.0;
  fLatitude = 0.0;

};

/*****************end of GReadPhotonGrISU ********************************/
GReadPhotonGrISU::~GReadPhotonGrISU() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadPhotonGrISU::~GReadPhotonGrISU" << endl;
  }
  if ( (sInFileStr != "") && (pInStream != 0) ) { 
    SafeDelete(pInStream);
  }
};
/*****************end of ~GReadPhotonGrISU ********************************/
bool GReadPhotonGrISU::setInputFile(const string &infile) {

  bool debugL = false;

  if (debugL) {
    *oLog << "  -- GReadPhotonGrISU::setInputFile " << endl;
  }

  sInFileStr = infile;  // file name from input parameter

  // open istream from filename or set istream to cin
  if (sInFileStr=="") {
    pInStream = &cin;   // set to cin if ""
  }
  else {
    pInStream = new ifstream(sInFileStr.c_str());
    if ( (pInStream->rdstate() & ifstream::failbit) ) {
      cerr << "    -- GReadPhotonGrISU::setInputFile " << endl;
      cerr << " could not open " << sInFileStr << endl;
      cerr << "      STOPPING CODE " << endl;
      exit(0);  //<! stop code if file cannot be opened
    }
  }

  // read header input records
  string headerStart("* HEADF");  // start of header flag
  string headerEnd("* DATAF");    // end of header flag
  string::size_type idx,idpt;
  string fileline = "";

  int recordCt = 0;     // set record counter
  sInFileHeader = "";   // set header string

  while (getline(*pInStream,fileline,'\n')) {  
    recordCt++;

    if (fileline=="") continue;  // check for blank lines

    // look for start of header
    idx = fileline.find(headerStart);
    if (idx != string::npos) {
      // start of header line
      //sInFileHeader = fileline + '\n';
      sInFileHeader = "";

      // read lines, looking for end of header
      while(getline(*pInStream,fileline,'\n') ) {
	if (fileline == "") continue;  // check for blank lines
       
        recordCt++;
        
        // look for end of header
        idx = fileline.find(headerEnd);
        if  (idx == string::npos) {
	  // we haven't reached the end of the header yet
	  if (recordCt > 1000) {  // stop after reading 1000 records
            *oLog << "header record count exceeds 1000 records" << endl;
            *oLog << "  check input file for end of header flag" << endl;
            *oLog << "   start of header flag '* HEADF'" << endl;
            *oLog << "   end of header flag   '* DATAF'" << endl;
            *oLog << "stopping code" << endl;
            exit(0);
          }
	  sInFileHeader += fileline + '\n';
	  // pull out the particle type (iParticleType)
	  if ( (idpt = fileline.find("PTYPE: ") ) != string::npos) {
	    string tmp = fileline.substr(idpt+6);
	    iParticleType = atoi(tmp.c_str());
	  }
        }
        else {
          break;   // break out of header while loop
        }
        
      }
    }   // header reader while  
    break;       // break out of main while loop
    
  }     // end of main while loop for reading header

  // now get R line, first check for blank lines
  // R line must be present
  while (getline(*pInStream,fileline,'\n')) {
    recordCt++;
    if (fileline!="") break; 
  }
  // read R line fGlobalEffic
  istringstream is(fileline);
  char cRecType;

  if (fileline[0]=='R') {
    is >> cRecType >> fGlobalEffic;  // get global efficiency from R line
    is.str("");
    is.clear();

  }

  else {
    *oLog << "R record out of place or nonexistant in input file" << endl;
    *oLog << "fileline read " << fileline << endl;
    *oLog << "    STOPPING CODE " << endl;
    exit(0);
  }

  // get H line, H line must be present, no blank line check
  getline(*pInStream,fileline,'\n');
  if (fileline[0] == 'H') {

    is.clear();
    is.str(fileline);
    is >> cRecType >> fObsHgt; // get obs.height from H line

    is.str("");
    is.clear();
  }
  else {
    *oLog << "H record out of place or nonexistant in input file" << endl;
    *oLog << "fileline read " << fileline << endl;
    *oLog << "    STOPPING CODE " << endl;
    exit(0);
  }
  
  getLine();  // fill the current record string for later use
  return true;
};
/*****************end of setInputFile ********************************/

bool GReadPhotonGrISU::getLine() {

  bool debugL = false;

  bool okLine = false;  // if false, eof indicated

  if (debugL) { 
    *oLog << "  -- GReadPhotonGrISU::getLine" << endl;
  }

  if (getline(*pInStream,sFileLine,'\n')) {
    char c = sFileLine[0];
    okLine = true;     // new record has been read

    if (c =='S') {       // primary record type
      eRecType = SREC;
    }
    else if (c == 'P') { // photon record type
      eRecType = PREC;
    }
    else if (c =='C') {   // C record
      eRecType = CREC;
    }
    else {
      cerr << "  unknown record type in input file" << endl;
      cerr << "  record:  " << sFileLine << endl;
      cerr << "   STOPING CODE " << endl;
    }
  }
  else {
    eRecType = EOFREC;  // record type is EOF
  }

  if (debugL) {
    DEBUGS(eRecType);
    DEBUGS(sFileLine);
    *oLog << "  end of GReadPhotonGrISU::getLine" << endl;
  }
  
  return okLine;
};
/*****************end of getLine ********************************/

bool GReadPhotonGrISU::getPrimary(ROOT::Math::XYZVector *pCore, 
				  ROOT::Math::XYZVector *pDCos,double *Az,
				  double *Zn, double *energy, 
				  unsigned int *particleType,
                                  double *firstIntHgt, double *firstIntDpt,
                                  unsigned int *showerid) {
  *particleType = iParticleType; // from header

  bool debugS = false;

  if (debugS) {
    *oLog << "  -- GReadPhotonGrISU::getPrimary" << endl;
    *oLog << "      " << sFileLine << endl;
  }
  
  bool okReadS = false;  // initialize to no S record

  if (eRecType == SREC) {
    okReadS = true;   // S record in current record
    char c;
    double tmp = 0.0;
    istringstream is(sFileLine);
    double fSXcore = 0.0;
    double fSYcore = 0.0;
    double fSZcore = 0.0;
    double fSXcos  = 0.0;
    double fSYcos  = 0.0;
    double fSZcos  = 0.0;

    //CD:2Mar2015 default values
    fFirstIntHgt = -999.9;
    fFirstIntDpt = -9999.9;
    iShowerID     = 99;

    // read primary parameters from string stream 
    is >> c >> fSEnergy >> fSXcore >> fSYcore >> fSXcos >> fSYcos
       >> tmp >> iSSeed[0] >> iSSeed[1] >> iSSeed[2];

    // the input coordinate system has x-axis East, y-axis South,
    // and z-axis Down. Our ground coor. system has x-axis East,
    // y-axis North, and z-axis Up. Thus, we have to reflect the y axis
    // after reading the shower ycore and ydir values. 
    fSYcore = -fSYcore;
    fSYcos  = -fSYcos;
    fSZcos  = -sqrt(1 - fSXcos*fSXcos - fSYcos*fSYcos);  
    
    fSZcore = 0.0;  // z location is relative to observatory height
    
    if (debugS) {
      DEBUGS(fSEnergy);
      DEBUGS(fSXcore);DEBUGS(fSYcore);DEBUGS(fSZcore); 
      DEBUGS(fSXcos);DEBUGS(fSYcos);
      DEBUGS(fSZcos); DEBUGS(tmp); DEBUGS( iSSeed[0]); DEBUGS(iSSeed[1]);
      DEBUGS(iSSeed[2]); 
      DEBUGS(iParticleType);
    }
    
    vSCore.SetXYZ(fSXcore,fSYcore,fSZcore);
    vSDcos.SetXYZ(fSXcos,fSYcos,fSZcos);

    if (debugS) {
      *oLog << "      Primary Core ";
      GUtilityFuncts::printGenVector(vSCore); *oLog << endl;
      *oLog << "      Primary DirCosines ";
      GUtilityFuncts::printGenVector(vSDcos); *oLog << endl;
    }
    *pCore = vSCore;  // set primary parameters
    *pDCos = vSDcos;
    *energy=fSEnergy;

    //CD:new parameters
    *firstIntHgt = fFirstIntHgt;
    *firstIntDpt = fFirstIntDpt;
    *showerid    = iShowerID;
   
    // get zenith angle and azimuth without adding wobble

    GUtilityFuncts::XYcosToAzZn(-fSXcos, -fSYcos,&fSAz,&fSZn);

    *Az = fSAz;
    *Zn = fSZn;

    if (debugS) {
      DEBUGS(fSZn*(TMath::RadToDeg()));
      DEBUGS(fSAz*(TMath::RadToDeg()));
    }

    if (debugS) *oLog << "        ready to get next line " << endl;

    // get C line if it's there
    getLine();
    istringstream is1(sFileLine);
    char c1;

    if (eRecType == CREC) {
      is1 >> c1 >> fFirstIntHgt >> fFirstIntDpt >> iShowerID;  // C line parms
      *firstIntHgt = fFirstIntHgt;
      *firstIntDpt = fFirstIntDpt;
      *showerid    = iShowerID;
      is1.str("");
      is1.clear();
      getLine();
    }
    //if C record is not present, set to default values

  }
  if (debugS) {
    *oLog << "   end of GReadPhotonGrISU::getPrimary" << endl;
  }

  return okReadS;
  
};
/*****************end of getPrimary ********************************/

bool GReadPhotonGrISU::getPhoton(ROOT::Math::XYZVector *pGrd,
                                   ROOT::Math::XYZVector *pDcos,
                                   double *pAz,double *pZn,
                                   double *pHgtEmiss,double *pTime,
                                   double *pWaveLgt, int *pType,
                                   int *pTel) {

  bool debugP = false;
  if (debugP) {
    *oLog << "  -- GReadPhotonGrISU::getPhoton" << endl;
    *oLog << "      " << sFileLine << endl;
  }
  
  bool okReadP = false; // initialize to no P record

  if (eRecType == PREC) {
    okReadP = true;  // we have a P record
    double xGrd = 0.0;
    double yGrd = 0.0;
    double zGrd = 0.0; // z at observatory height
    double xcos = 0.0;
    double ycos = 0.0;
    double zcos = 0.0;
    
    // read photon parameters from stringstream
    char c;
    istringstream is(sFileLine);
    is >> c >> xGrd >> yGrd >> xcos >> ycos >> fPHgtEmiss
       >> fPTime >> fPWaveLgt >> iPType >> iPTel;  
    
    // convert to ground coordinate system from kascade system
    // input system has z down, y South rather than z up and y North
    yGrd = -yGrd;
    ycos = -ycos;
    zcos = -sqrt(1 - xcos*xcos - ycos*ycos);

    if (debugP) {
      DEBUGS(xGrd);DEBUGS(yGrd);DEBUGS(zGrd);DEBUGS(xcos);
      DEBUGS(ycos);DEBUGS(zcos);DEBUGS(fPHgtEmiss);DEBUGS(fPTime);
      DEBUGS(fPWaveLgt);DEBUGS(iPType);DEBUGS(iPTel);
    }
    vPGrd.SetXYZ(xGrd,yGrd,zGrd);
    vPDcos.SetXYZ(xcos,ycos,zcos);

    if (debugP) {
    *oLog << "      Photon ground ";
    GUtilityFuncts::printGenVector(vPGrd); *oLog << endl;
    *oLog << "      Photon DirCosines ";
    GUtilityFuncts::printGenVector(vPDcos); *oLog << endl;
    }

    GUtilityFuncts::XYcosToAzZn(-xcos,-ycos,
                   &fPAz, &fPZn);
   
    if (debugP) {
      DEBUGS(fPZn*(TMath::RadToDeg()));
      DEBUGS(fPAz*(TMath::RadToDeg()));
    }
    *pGrd = vPGrd;
    *pDcos = vPDcos;
    *pAz = fPAz;
    *pZn = fPZn;
    *pHgtEmiss = fPHgtEmiss;
    *pTime = fPTime;
    *pWaveLgt = fPWaveLgt;
    *pType = iPType;
    *pTel = iPTel;
 
    getLine();  // get next line for use later
  }
  if (debugP) {
    *oLog << "   end of GReadPhotonGrISU::getPhoton" << endl;
  }

  return okReadP;
};

/*****************end of getPhoton ********************************/

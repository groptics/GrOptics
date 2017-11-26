/*
VERSION4.0
30May2016
*/
/*!  GReadSegSCStd.cpp
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

#include "TGraph.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"
#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "ATelescope.h"
#include "GSegSCTelescope.h"

#include "ATelescopeFactory.h"
#include "GSegSCTelescopeFactory.h"

#include "GReadSegSCStd.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

GReadSegSCStd::GReadSegSCStd(const string &pilotfile,GSegSCTelescopeFactory *SegSCFactory  ) {

  // initialize variables
  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::GReadSegSCStd" << endl;
  }

  SegSCFac = 0;
  spilotfile = "";
  pi = 0;
  flagline = "";
  iStdNum = 0;
  opt = 0;
  
  SegSCFac = SegSCFactory;
  spilotfile = pilotfile;
  *oLog << "spilotfile " << spilotfile << endl;
  setupSegSCFactory();

};
/****************** end of GReadSegSCStd **********/

GReadSegSCStd::GReadSegSCStd(const string &pilotfile) {
  spilotfile = pilotfile;

  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::GReadSegSCStd(const string pilotfile) " 
          << pilotfile << endl;
  }

 // initialize variables
  SegSCFac = 0;
  pi = 0;
  flagline = "";
  iStdNum = 0;
  opt = 0;
 
  
};
/****************** end of GReadSegSCStd **********/

GReadSegSCStd::~GReadSegSCStd() { 
  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::~GReadSegSCStd " << endl;
  }
  //SafeDelete(pi); 

};
/****************** end of ~GReadSegSCStd **********/

void GReadSegSCStd::
setSegSCTelescopeFactory(GSegSCTelescopeFactory *SegSCFactory) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::setSegSCTelescopeFactory" << endl;
  }

  if (SegSCFac != 0 ) {
    cerr << "GReadSegSCStd::setSCTelescopeFactory " << endl;
    cerr << "   ERROR: SegSCFac pointer to SegSCFactory previously set" << endl;
    cerr << "   stopping code " << endl;
    exit(0);
  }

  SegSCFac = SegSCFactory;
  setupSegSCFactory();
  
};

/****************** end of setSCTelescopeFactory **********/
void GReadSegSCStd::setupSegSCFactory() {

  bool debug = false;

  if (debug) {
    *oLog << "  -- GReadSegSCStd::setupSegSCFactory" << endl;
    //*oLog << "       spilotfile " << spilotfile << endl;
  }

  if (SegSCFac == 0) {
    cerr << " GReadSegSCStd doesn't have a SegSCFactory object" << endl;
    cerr << "exiting code" << endl;
    exit(0);
  }

  getReflCoeff();

  pi = new GPilot(spilotfile);

  // check for additional pilot files to append to this file in the pilotreader
  string flag = "PILOTF";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    string newPilotFile = tokens.at(0);
    pi->addPilotFile(newPilotFile);
  }

  flag = "TELSEGSCSTD";
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {

    int iStdOptNum = atoi(tokens.at(0).c_str());
    SegSCFac->mStdOptics[iStdOptNum] =  new SegSCStdOptics();
    opt = SegSCFac->mStdOptics[iStdOptNum];

    // set pointer in SCStdOptics for reflection coefficients.
    // maps to coefficients are now in all SCStdOptics objects.
    //opt->mVReflWaveLgts = SegSCFac->mVReflWaveLgts;
    //opt->mVCoeffs = SegSCFac->mVCoeffs;

    // set standard number (also available as map key)
    opt->stdNum = iStdOptNum;

    double fFocLgt = atof(tokens.at(1).c_str());
    opt->fF = fFocLgt;

    double fAvgTransitTime = atof(tokens.at(2).c_str());
    opt->fAvgTransitTime = fAvgTransitTime;

    double rotationOffset = atof(tokens.at(3).c_str());
    opt->fRotationOffset = rotationOffset;

    if (tokens.size() >= 2) {
      int iprintmode = atoi(tokens.at(4).c_str());
      opt->iPrtMode = iprintmode;
    }
  }
  // read PLATESCALE record ///////////////
  flag = "PLATESCALE"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    opt->fPlateScaleFactor = atof(tokens.at(1).c_str() );
  }

  // read PRIMARY record ///////////////
  flag = "PRIMARY"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
 
    opt->fRpMax = atof(tokens.at(1).c_str() );
    opt->fRpMin = atof(tokens.at(2).c_str());
    if (tokens.size() == 4) {
      opt->fZp = atof(tokens.at(3).c_str());
    }
  }

  // read SECONDARY record ///////////////// 
  flag = "SECONDARY"; 
   pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    opt->fRsMax = atof(tokens.at(1).c_str() );
    opt->fRsMin = atof(tokens.at(2).c_str());
    // initialized to 1.5 meters in opt
    if (tokens.size() == 4) {
      opt->fZs = atof(tokens.at(3).c_str());
    }
  }

  // read PRIMSEGP1, both BASIC and CHANGE
  flag = "PRIMSEGP1"; 
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   

    if (tokens[1] == "BASIC") {
      readBasicRecord(tokens,P1,opt);
    }
    else if (tokens[1] == "CHANGE") {
      readChangeRecord(tokens,P1,opt);
    }
  }
  // read PRIMSEGP2, both BASIC and CHANGE
  flag = "PRIMSEGP2"; 
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   

    if (tokens[1] == "BASIC") {
      readBasicRecord(tokens,P2,opt);
    }
    else if (tokens[1] == "CHANGE") {
      readChangeRecord(tokens,P2,opt);
    }
  }
  // read SECONDSEGS1, both BASIC and CHANGE
  flag = "SECONDSEGS1"; 
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   

    if (tokens[1] == "BASIC") {
      readBasicRecord(tokens,S1,opt);
    }
    else if (tokens[1] == "CHANGE") {
      readChangeRecord(tokens,S1,opt);
    }
  }
  // read SECONDSEGS2, both BASIC and CHANGE
  flag = "SECONDSEGS2"; 
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   

    if (tokens[1] == "BASIC") {
      readBasicRecord(tokens,S2,opt);
    }
    else if (tokens[1] == "CHANGE") {
      readChangeRecord(tokens,S2,opt);
    }
  }
    
  flag = "FOCALSURFACE"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    opt->fKappa1 = atof(tokens.at(1).c_str() );
    opt->fKappa2 = atof(tokens.at(2).c_str());
    opt->fRf = atof(tokens.at(3).c_str());
    // see initialization in opt 
    if ( tokens.size() == 5) {
      opt->fZf = atof(tokens.at(3).c_str());
    }
  }

  flag = "FOCALSURFACEOFFSET"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    opt->fFocalSurfaceXOffset = atof(tokens.at(1).c_str() );
    opt->fFocalSurfaceYOffset = atof(tokens.at(2).c_str() );
    opt->fFocalSurfaceZOffset = atof(tokens.at(3).c_str() );
    opt->fFocalSurfacePhiOffset = atof(tokens.at(4).c_str() );
    opt->fFocalSurfaceThetaOffset = atof(tokens.at(5).c_str() );
    opt->fFocalSurfacePsiOffset = atof(tokens.at(6).c_str() );
  }

  flag = "CAMERA"; //CAMERA 1 1 6.75 0.0 54.0 32.7 0.1 0.0 0.0 0.0 0.0 1.55 0 2
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    int cflag = atoi(tokens.at(1).c_str() );
    if (cflag) {
      opt->bCameraFlag = true;
    }
    else {
      opt->bCameraFlag = false;
    }
    opt->fPixelSize = atof(tokens.at(2).c_str() );
    //opt->fPixelPitch = atof(tokens.at(3).c_str());
    opt->fMAPMTWidth = atof(tokens.at(4).c_str());
    opt->fMAPMTLength = atof(tokens.at(5).c_str());
    opt->fInputWindowThickness = atof(tokens.at(6).c_str());
    //opt->fMAPMTAngularSize = atof(tokens.at(7).c_str());
    opt->fMAPMTOffset_x =  atof(tokens.at(7).c_str());
    opt->fMAPMTOffset_y =  atof(tokens.at(8).c_str());
    opt->fMAPMTOffset =  atof(tokens.at(9).c_str()); //this is z axis offset
    opt->fMAPMTRoll =  atof(tokens.at(10).c_str());
    opt->fMAPMTPitch =  atof(tokens.at(11).c_str());
    opt->fMAPMTGap =  atof(tokens.at(12).c_str());
    if (tokens.size() > 12) {
      opt->fMAPMTRefIndex =  atof(tokens.at(13).c_str());
    }
    if (tokens.size() > 13) {
      Int_t tmpi = atoi(tokens.at(14).c_str());
      if (tmpi > 0) 
        opt->bSingleMAPMTmodule = true;
      else {
        opt->bSingleMAPMTmodule = false;
      }
    }
    if (tokens.size() > 14) {
      opt->fSubCells =  atof(tokens.at(15).c_str());
    }
  }

  flag = "WINDOW"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    int cflag = atoi(tokens.at(1).c_str() );
    if (cflag) {
      opt->bEntranceWindowFlag = true;
    }
    else {
      opt->bEntranceWindowFlag = false;
    }
    opt->fEntranceWindowThickness = atof(tokens.at(2).c_str() );
      opt->iEntranceTrCurveIndex = atof(tokens.at(3).c_str()) ;
      if(atof(tokens.at(3).c_str()) > 0 ) {
          getTranCurve();
      }
    opt->fEntranceWindowN = atof(tokens.at(4).c_str() );
    opt->fEntranceWindowOffset = atof(tokens.at(5).c_str() );
    int aflag = atof(tokens.at(6).c_str() );
    if (aflag) {
      opt->bEntranceWindowAbsFlag = true;
    }
    else {
      opt->bEntranceWindowAbsFlag = false;
    }
    opt->fEntranceWindowAbsLength = atof(tokens.at(7).c_str() );
  }

  flag = "PRIMFRAME";
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str());
    opt = SegSCFac->mStdOptics[iStdOptNum];
    int cflag = atoi(tokens.at(1).c_str());
    if (cflag) {
      opt->frameFlag = true;
    } else {
      opt->frameFlag = false;
    }
  }

  flag = "PRIMBAFFLE"; 
  pi->set_flag(flag);


  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    int cflag = atoi(tokens.at(1).c_str() );
    if (cflag) {
      opt->bpBaffleFlag = true;
    }
    else {
      opt->bpBaffleFlag = false;
    }
    opt->fpBLen = atof(tokens.at(2).c_str() );
    opt->fpBTilt = atof(tokens.at(3).c_str() );
    opt->fpBRadOffset = atof(tokens.at(4).c_str() );
    opt->fpBZOffset = atof(tokens.at(5).c_str() );
  }  

  flag = "SECONDBAFFLE"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    int cflag = atoi(tokens.at(1).c_str() );
    if (cflag) {
      opt->bsBaffleFlag = true;
    }
    else {
      opt->bsBaffleFlag = false;
    }
    opt->fsBLen = atof(tokens.at(2).c_str() );
    opt->fsBTilt = atof(tokens.at(3).c_str() );
    opt->fsBRadOffset = atof(tokens.at(4).c_str() );
    opt->fsBZOffset = atof(tokens.at(5).c_str() );
  }  

  /*
  flag = "REFLID"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SegSCFac->mStdOptics[iStdOptNum];   
    opt->iPrimReflID = atoi(tokens.at(1).c_str() );
    opt->iSecReflID = atoi(tokens.at(2).c_str());
  }
  */

  // calculate plate scale factor if no record present.
  if ( (opt->fPlateScaleFactor) < 0.00001) {
    opt->fPlateScaleFactor = tan( 1.0*(TMath::DegToRad())) *
      opt->fF * 100.0;
  }  
 
  SafeDelete(pi);

  getPolyCoeffs();
  // 

  //if (iPrtMode > 0) {
  if (debug) {
    *oLog << endl;
    *oLog << "       printing all SEGSC standard optics objects" << endl;
    for (SegSCFac->itmStdOp = SegSCFac->mStdOptics.begin();
    	 SegSCFac->itmStdOp != SegSCFac->mStdOptics.end();
    	 SegSCFac->itmStdOp++) {
      (*(SegSCFac->itmStdOp)).second->printSegSCStdOptics();
    }
  }

  if (debug) {
    *oLog << "  -- end of setupSegSCFactory" << endl;
  }

};

/******************** end of setupSegSCFactory ****************/
 void GReadSegSCStd::readBasicRecord(const vector<string> &tokens1,
                                     const MirSeg &eMirPS,
                                     SegSCStdOptics *opt1) {
   mirrorSegmentDetails segTmp;
   Int_t numMirrors = atoi(tokens1.at(2).c_str());  

   //opt1->iNumP1Mirrors = atoi(tokens1.at(2).c_str()); 
   segTmp.rmin =  atof(tokens1.at(3).c_str());
   segTmp.rmax =  atof(tokens1.at(4).c_str());
   segTmp.margin =  atof(tokens1.at(5).c_str());
   segTmp.delPhi =  atof(tokens1.at(6).c_str());
   segTmp.reflect = atoi(tokens1.at(7).c_str() );
   segTmp.roughness =  atof(tokens1.at(8).c_str());
   segTmp.posErrorX =  atof(tokens1.at(9).c_str());
   segTmp.posErrorY =  atof(tokens1.at(10).c_str());
   segTmp.posErrorZ =  atof(tokens1.at(11).c_str());
   segTmp.rotErrorPhi =  atof(tokens1.at(12).c_str());
   segTmp.rotErrorTheta =  atof(tokens1.at(13).c_str());
   segTmp.rotErrorPsi =  atof(tokens1.at(14).c_str());
   segTmp.bRead = 0;

   vector<mirrorSegmentDetails *> vSeg;
   for (int i = 0;i<numMirrors;i++) {
     vSeg.push_back(new mirrorSegmentDetails);
     *vSeg[i] = segTmp; 
     //vSeg.push_back(segTmp);
     //vSeg.push_back(tmpvv);
   }

   if (eMirPS == P1) {
     opt1->iNumP1Mirrors = numMirrors;
     opt1->vSegP1 = vSeg;
   }
   else if (eMirPS == P2) {
     opt1->iNumP2Mirrors = numMirrors;
     opt1->vSegP2 = vSeg;

   }
   else if (eMirPS == S1) {
     opt1->iNumS1Mirrors = numMirrors;
     opt1->vSegS1 = vSeg;

   }

  else if (eMirPS == S2) {
     opt1->iNumS2Mirrors = numMirrors;
     opt1->vSegS2 = vSeg;

   }  

 };

/******************** end of readBasicRecord ****************/
void GReadSegSCStd::readChangeRecord(const vector<string> &tokens1,
                       const MirSeg &eMirPS,
                                     SegSCStdOptics *opt1) {
  mirrorSegmentDetails segTmp;
  vector<int> vListSeg;    
  vector<int>::iterator itv;
  GUtilityFuncts::decodeMatlabString(tokens1[2],vListSeg);
  *oLog << "VLISTSEG " << endl;
  GUtilityFuncts::printVector(vListSeg);

  segTmp.rmin =  atof(tokens1.at(3).c_str());
  segTmp.rmax =  atof(tokens1.at(4).c_str());
  segTmp.margin =  atof(tokens1.at(5).c_str());
  segTmp.delPhi =  atof(tokens1.at(6).c_str());
  segTmp.reflect = atoi(tokens1.at(7).c_str() );
  segTmp.roughness =  atof(tokens1.at(8).c_str());
  segTmp.posErrorX =  atof(tokens1.at(9).c_str());
  segTmp.posErrorY =  atof(tokens1.at(10).c_str());
  segTmp.posErrorZ =  atof(tokens1.at(11).c_str());
  segTmp.rotErrorPhi =  atof(tokens1.at(12).c_str());
  segTmp.rotErrorTheta =  atof(tokens1.at(13).c_str());
  segTmp.rotErrorPsi =  atof(tokens1.at(14).c_str());
  segTmp.bRead = 1;

  for (unsigned int i = 0;i<vListSeg.size();i++) {
 
    int ix = vListSeg[i];
 
    if (eMirPS == P1) {
      // update element
      *( (opt1->vSegP1)[ix -1] ) = segTmp;
   }
    else if (eMirPS == P2) {
      // update element
      *( (opt1->vSegP2)[ix -1] ) = segTmp;
    }
    else if (eMirPS == S1) {
      // update element
      *( (opt1->vSegS1)[ix -1] ) = segTmp;
    }
    else if (eMirPS == S2) {
      // update element
      *( (opt1->vSegS2)[ix -1] ) = segTmp;
    }
  } 

};
/******************** end of readChangeRecord ****************/

void GReadSegSCStd::getPolyCoeffs() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::getPolyCoeffs " << endl;
  }

  ifstream inFile(spilotfile.c_str(),ios::in);
  if (! inFile) {
    cerr << "  -- GReadSegSCStd::getPolyCoeffs " << endl;
    cerr << "    could not open file: " << spilotfile << endl;
    exit(0);
  }

  string pilotline;

  // get primary polynomial coefficients
  while( getline(inFile,pilotline,'\n')) {
    vector<string> tokens1;
    GUtilityFuncts::tokenizer(pilotline,tokens1);
    if ( (tokens1.size() > 0) &&
         (tokens1.at(0) == "*")  &&
         (tokens1.at(1) == "PRIMCOEFF") )  {
      int index = atoi(tokens1.at(2).c_str());
      int number = atoi(tokens1.at(3).c_str());
      opt = SegSCFac->mStdOptics[index];
      opt->iNParP = number;
      // read and store the coefficients
      for (int i = 0;i<number;i++) {
        double polyc     = 0.0;
        inFile >> polyc;
	opt->fzp.push_back(polyc);
      }
    }
  }       

  // somewhat klugy, but that's ok for now
  ifstream inFile1(spilotfile.c_str(),ios::in);

  while( getline(inFile1,pilotline,'\n')) {
    vector<string> tokens1;
    GUtilityFuncts::tokenizer(pilotline,tokens1);
    if ( (tokens1.size() > 0) &&
         (tokens1.at(0) == "*")  &&
         (tokens1.at(1) == "SECONDCOEFF") )  {
      int index = atoi(tokens1.at(2).c_str());
      int number = atoi(tokens1.at(3).c_str());
      opt = SegSCFac->mStdOptics[index];
      opt->iNParS = number;
      // read and store the coefficients
      for (int i = 0;i<number;i++) {
        double polyc     = 0.0;
        inFile1 >> polyc;
	opt->fzs.push_back(polyc);
      }
    }
  }       
  
};
/******************** end of getPolyCoeffs ****************/

void GReadSegSCStd::getReflCoeff() {

  // wavelengths in the config. file have nm units
  // here convert to cm as required robast usage
  bool debug = false;

  if (debug) {
    *oLog << "  -- GReadSegSCStd::getReflCoeff() " << endl;
    *oLog << "       spilotfile " << spilotfile << endl;
  }

  ifstream inFile(spilotfile.c_str(),ios::in);
  if (! inFile) {
    cerr << "  -- GReadSegSCStd::getReflCoeff " << endl;
    cerr << "    could not open file: " << spilotfile << endl;
    exit(0);
  }

  string pilotline;
  while( getline(inFile,pilotline,'\n')) {
    vector<string> tokens1;
    GUtilityFuncts::tokenizer(pilotline,tokens1);
    if ( (tokens1.size() > 0) &&
         (tokens1.at(0) == "*")  &&
         (tokens1.at(1) == "RFCRV") )  {

      int index = atoi(tokens1.at(2).c_str());
      int number = atoi(tokens1.at(3).c_str());
      map<int, TGraph *>::iterator itGr;
      itGr =  SegSCFac->mGRefl->find(index);

      if (itGr!=SegSCFac->mGRefl->end()) {
        cerr << "trying to load reflection curve " << index << endl;
        cerr << "but curve already in the map: check pilot files " << endl;
        exit(0);
      }

      (*(SegSCFac->mGRefl))[index] = new TGraph(number);
      TGraph *grTmp = (*(SegSCFac->mGRefl))[index];

      // get ready to read and load reflection coefficients

      // read and store the reflection coefficients
      for (int i = 0;i<number;i++) {
        double wavel     = 0.0;
        double reflcoeff = 0.0;
        inFile >> wavel >> reflcoeff;
	// convert to cm from nm.
	wavel = wavel*1.0e-07;
        grTmp->SetPoint(i,wavel,reflcoeff);
      }
      getline(inFile,pilotline);// finish the end of the line
   
      if (debug) {
	*oLog << "       reading in reflectivity curve " << index 
	      << "   with " << number << " points" << endl;
	*oLog << "             converting from nm to cm for robast" << endl;
      }
      
      if (debug) {
        *oLog << "      print reflectance vectors for index: " << index 
             << "  number of points " << number << endl;
	*oLog << "             wavelengths stored in cm, not nm" << endl;
        for (int i=0;i<number;i++) {
	  double x = 0.0;
	  double y = 0.0;
	  grTmp->GetPoint(i,x,y);
	  *oLog << "          " << setw(4) << i << "   " 
		<< x/1.0e-07 << "     "
		<<  y << endl;
        }
        *oLog << endl;
      }        
    }
  }

  if (SegSCFac->mGRefl->size() == 0) {
    *oLog << "no reflection coeffient table found in config file" << endl;
    *oLog << "    STOPPING CODE" << endl;
    exit(0);
  }
  
};

/******************** end of getReflCoeff ****************/

void GReadSegSCStd::getTranCurve() {

    // wavelengths in the config. file have nm units
    // here convert to cm as required robast usage
    bool debug = false;

    if (debug) {
        *oLog << "  -- GReadSegSCStd::getTranCurve() " << endl;
        *oLog << "       spilotfile " << spilotfile << endl;
    }

    ifstream inFile(spilotfile.c_str(),ios::in);
    if (! inFile) {
        cerr << "  -- GReadSegSCStd::getTranCurve " << endl;
        cerr << "    could not open file: " << spilotfile << endl;
        exit(0);
    }

    string pilotline;
    while( getline(inFile,pilotline,'\n')) {
        vector<string> tokens1;
        GUtilityFuncts::tokenizer(pilotline,tokens1);
        if ( (tokens1.size() > 0) &&
             (tokens1.at(0) == "*")  &&
             (tokens1.at(1) == "TRCRV") )  {

            int index = atoi(tokens1.at(2).c_str());
            int number = atoi(tokens1.at(3).c_str());

            map<int, TGraph *>::iterator itGrAbs;
            itGrAbs =  SegSCFac->mGTranAbsLength->find(index);

            if (itGrAbs!=SegSCFac->mGTranAbsLength->end()) {
                cerr << "trying to load transmittance curve " << index << endl;
                cerr << "but curve already in the map: check pilot files " << endl;
                exit(0);
            }

            map<int, TGraph *>::iterator itGrN;
            itGrAbs =  SegSCFac->mGTranN->find(index);

            if (itGrAbs!=SegSCFac->mGTranN->end()) {
                cerr << "trying to load transmittance curve " << index << endl;
                cerr << "but curve already in the map: check pilot files " << endl;
                exit(0);
            }

            (*(SegSCFac->mGTranAbsLength))[index] = new TGraph(number);
            TGraph *grTmpAbs = (*(SegSCFac->mGTranAbsLength))[index];

            (*(SegSCFac->mGTranN))[index] = new TGraph(number);
            TGraph *grTmpN = (*(SegSCFac->mGTranN))[index];

            // get ready to read and load transmittance coefficients

            // read and store the transmittance coefficients
            for (int i = 0;i<number;i++) {
                double wavel     = 0.0;
                double refindex = 0.0;
                double absl = 0.0;
                inFile >> wavel >> refindex >> absl;
                // convert to cm from nm.
                wavel = wavel*1.0e-07;
                grTmpN->SetPoint(i,wavel,refindex);
                grTmpAbs->SetPoint(i,wavel,absl);
            }
            getline(inFile,pilotline);// finish the end of the line

            if (debug) {
                *oLog << "       reading in transmittance curve " << index
                      << "   with " << number << " points" << endl;
                *oLog << "             converting from nm to cm for robast" << endl;
            }

            if (debug) {
                *oLog << "      print transmittance vectors for index: " << index
                      << "  number of points " << number << endl;
                *oLog << "             wavelengths stored in cm, not nm" << endl;
                for (int i=0;i<number;i++) {
                    double x = 0.0;
                    double y = 0.0;
                    double z = 0.0;
                    grTmpN->GetPoint(i,x,y);
                    grTmpAbs->GetPoint(i,x,z);
                    *oLog << "          " << setw(4) << i << "   "
                          << x/1.0e-07 << "     "
                          <<  y <<  z << endl;
                }
                *oLog << endl;
            }
        }
    }

    if (SegSCFac->mGTranN->size() == 0) {
        *oLog << "no transmittance coeffient table found in config file" << endl;
        *oLog << "    STOPPING CODE" << endl;
        exit(0);
    }

};


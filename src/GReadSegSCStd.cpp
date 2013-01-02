/*
VERSION2.7aBeta
17Dec2012
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

#include "GTelescope.h"
#include "GSegSCTelescope.h"

#include "GTelescopeFactory.h"
#include "GSegSCTelescopeFactory.h"

#include "GReadSegSCStd.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

GReadSegSCStd::GReadSegSCStd(const string &pilotfile,GSegSCTelescopeFactory *SCFactory  ) {

  // initialize variables
  bool debug = true;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::GReadSegSCStd" << endl;
  }

  SCFac = 0;
  spilotfile = "";
  pi = 0;
  flagline = "";
  iStdNum = 0;
  opt = 0;
  
  SCFac = SCFactory;
  spilotfile = pilotfile;
  *oLog << "spilotfile " << spilotfile << endl;
  setupSCFactory();

};
/****************** end of GReadSegSCStd **********/

GReadSegSCStd::GReadSegSCStd(const string &pilotfile) {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::GReadSegSCStd" << endl;
  }
  /*
 // initialize variables
  SCFac = 0;
  spilotfile = "";
  pi = 0;
  flagline = "";
  iStdNum = 0;
  opt = 0;
 
  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadSegSCStd(const string pilotfile) " << pilotfile << endl;
  }

  SCFac = 0;
  */ 
  spilotfile = pilotfile;
  
};
/****************** end of GReadSegSCStd **********/

GReadSegSCStd::~GReadSegSCStd() { 
  bool debug = true;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::~GReadSegSCStd " << endl;
  }
  //SafeDelete(pi); 

};
/****************** end of ~GReadSegSCStd **********/

void GReadSegSCStd::
setSegSCTelescopeFactory(GSegSCTelescopeFactory *SCFactory) {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::setSegSCTelescopeFactory" << endl;
  }

  if (SCFac != 0 ) {
    cerr << "GReadSegSCStd::setSCTelescopeFactory " << endl;
    cerr << "   ERROR: SCFac pointer to SCFactory previously set" << endl;
    cerr << "   stopping code " << endl;
    exit(0);
  }

  SCFac = SCFactory;
  setupSCFactory();
  
};

/****************** end of setSCTelescopeFactory **********/
void GReadSegSCStd::setupSCFactory() {

  bool debug = true;

  if (debug) {
    *oLog << "  -- GReadSegSCStd::setupSCFactory" << endl;
    //*oLog << "       spilotfile " << spilotfile << endl;
  }

  if (SCFac == 0) {
    cerr << " GReadSegSCStd doesn't have a SCFactory object" << endl;
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
    SCFac->mStdOptics[iStdOptNum] =  new SegSCStdOptics();
    opt = SCFac->mStdOptics[iStdOptNum];

    // set pointer in SCStdOptics for reflection coefficients.
    // maps to coefficients are now in all SCStdOptics objects.
    //opt->mVReflWaveLgts = SCFac->mVReflWaveLgts;
    //opt->mVCoeffs = SCFac->mVCoeffs;

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
  // read PRIMARY record ///////////////
  flag = "PRIMARY"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
 
    opt->fRpMax = atof(tokens.at(1).c_str() );
    opt->fRpMin = atof(tokens.at(2).c_str());
    opt->fZp = atof(tokens.at(3).c_str());
  }

  // read SECONDARY record ///////////////// 
  flag = "SECONDARY"; 
   pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    opt->fRsMax = atof(tokens.at(1).c_str() );
    opt->fRsMin = atof(tokens.at(2).c_str());
    opt->fZs = atof(tokens.at(3).c_str());
  }

  // read PRIMSEGP1, both BASIC and CHANGE
  flag = "PRIMSEGP1"; 
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   

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
    opt = SCFac->mStdOptics[iStdOptNum];   

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
    opt = SCFac->mStdOptics[iStdOptNum];   

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
    opt = SCFac->mStdOptics[iStdOptNum];   

    if (tokens[1] == "BASIC") {
      readBasicRecord(tokens,S2,opt);
    }
    else if (tokens[1] == "CHANGE") {
      readChangeRecord(tokens,S2,opt);
    }
  }
    
  /*
  flag = "PRIMARYOFFSET"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    opt->fPrimXoffset = atof(tokens.at(1).c_str() );
    opt->fPrimYoffset = atof(tokens.at(2).c_str() );
    opt->fPrimZoffset = atof(tokens.at(3).c_str() );
    opt->fPrimThetaOffset = atof(tokens.at(4).c_str() );
    opt->fPrimPhiAngle = atof(tokens.at(5).c_str() );
    if (tokens.size() > 6) 
      opt->fPrimRoughSigma = atof(tokens.at(6).c_str() );
  }
  */

   /*
  flag = "SECONDARYOFFSET"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    opt->fSecondXoffset = atof(tokens.at(1).c_str() );
    opt->fSecondYoffset = atof(tokens.at(2).c_str() );
    opt->fSecondZoffset = atof(tokens.at(3).c_str() );
    opt->fSecondThetaOffset = atof(tokens.at(4).c_str() );
    opt->fSecondPhiAngle = atof(tokens.at(5).c_str() );
    if (tokens.size() > 6) 
      opt->fSecondRoughSigma = atof(tokens.at(6).c_str() );
  }
  */

  flag = "FOCALPLANE"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    opt->fKappa1 = atof(tokens.at(1).c_str() );
    opt->fKappa2 = atof(tokens.at(2).c_str());
    opt->fRf = atof(tokens.at(3).c_str());
  }
  /*
  flag = "FOCALPLANEOFFSET"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    opt->fFocalPlXoffset = atof(tokens.at(1).c_str() );
    opt->fFocalPlYoffset = atof(tokens.at(2).c_str() );
    opt->fFocalPlZoffset = atof(tokens.at(3).c_str() );
    opt->fFocalPlThetaOffset = atof(tokens.at(4).c_str() );
    opt->fFocalPlPhiAngle = atof(tokens.at(5).c_str() );
  }
  */
  flag = "CAMERA"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    int cflag = atoi(tokens.at(1).c_str() );
    if (cflag) {
      opt->bCameraFlag = true;
    }
    else {
      opt->bCameraFlag = false;
    }
    opt->fPixelSize = atof(tokens.at(2).c_str() );
    opt->fPixelPitch = atof(tokens.at(3).c_str());
    opt->fMAPMTWidth = atof(tokens.at(4).c_str());
    opt->fMAPMTLength = atof(tokens.at(5).c_str());
    opt->fInputWindowThickness = atof(tokens.at(6).c_str());
    opt->fMAPMTAngularSize = atof(tokens.at(7).c_str());
    opt->fMAPMTOffset =  atof(tokens.at(8).c_str());
    opt->fMAPMTGap =  atof(tokens.at(9).c_str());
    if (tokens.size() == 11) {
      opt->fMAPMTRefIndex =  atof(tokens.at(10).c_str());
    }
  }
  /*
  flag = "REFLID"; 
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens.at(0).c_str() );
    opt = SCFac->mStdOptics[iStdOptNum];   
    opt->iPrimReflID = atoi(tokens.at(1).c_str() );
    opt->iSecReflID = atoi(tokens.at(2).c_str());
  }
  */
  SafeDelete(pi);

  getPolyCoeffs();
  // 

  //if (iPrtMode > 0) {

  if (1) {
    *oLog << endl;
    *oLog << "       printing all SC standard optics objects" << endl;
    for (SCFac->itmStdOp = SCFac->mStdOptics.begin();
    	 SCFac->itmStdOp != SCFac->mStdOptics.end();
    	 SCFac->itmStdOp++) {
      (*(SCFac->itmStdOp)).second->printSegSCStdOptics();
    }
  }

  if (debug) {
    *oLog << "  -- end of setupSCFactory" << endl;
  }

};

/******************** end of setupSCFactory ****************/
 void GReadSegSCStd::readBasicRecord(const vector<string> &tokens,
                                     const MirSeg &eMirPS,
                                     SegSCStdOptics *opt) {
   mirrorSegmentDetails * segTmp = new mirrorSegmentDetails();
   Int_t numMirrors = atoi(tokens.at(2).c_str());  

   //opt->iNumP1Mirrors = atoi(tokens.at(2).c_str()); 
   segTmp->rmin =  atof(tokens.at(3).c_str());
   segTmp->rmax =  atof(tokens.at(4).c_str());
   segTmp->margin =  atof(tokens.at(5).c_str());
   segTmp->delPhi =  atof(tokens.at(6).c_str());
   segTmp->reflect = atoi(tokens.at(7).c_str() );
   segTmp->roughness =  atof(tokens.at(8).c_str());
   segTmp->posErrorX =  atof(tokens.at(9).c_str());
   segTmp->posErrorY =  atof(tokens.at(10).c_str());
   segTmp->posErrorZ =  atof(tokens.at(11).c_str());
   segTmp->rotErrorX =  atof(tokens.at(12).c_str());
   segTmp->rotErrorY =  atof(tokens.at(13).c_str());
   segTmp->rotErrorZ =  atof(tokens.at(14).c_str());
   segTmp->bRead = 0;

   vector<mirrorSegmentDetails *> vSeg;
   for (int i = 0;i<numMirrors;i++) {
     vSeg.push_back(segTmp);
   }

   if (eMirPS == P1) {
     opt->iNumP1Mirrors = numMirrors;
     opt->vSegP1 = vSeg;
   }
   else if (eMirPS == P2) {
     opt->iNumP2Mirrors = numMirrors;
     opt->vSegP2 = vSeg;

   }
   else if (eMirPS == S1) {
     opt->iNumS1Mirrors = numMirrors;
     opt->vSegS1 = vSeg;

   }

  else if (eMirPS == S2) {
     opt->iNumS2Mirrors = numMirrors;
     opt->vSegS2 = vSeg;

   }  
 };

/******************** end of readBasicRecord ****************/
void GReadSegSCStd::readChangeRecord(const vector<string> &tokens,
                       const MirSeg &eMirPS,
                                     SegSCStdOptics *opt) {
  mirrorSegmentDetails *segTmp = new mirrorSegmentDetails();
  vector<int> vListSeg;    
  vector<int>::iterator itv;
  GUtilityFuncts::decodeMatlabString(tokens[2],vListSeg);
  *oLog << "VLISTSEG " << endl;
  GUtilityFuncts::printVector(vListSeg);
  for (int i = 0;i<vListSeg.size();i++) {
    int ix = vListSeg[i];
    segTmp->rmin =  atof(tokens.at(3).c_str());
    segTmp->rmax =  atof(tokens.at(4).c_str());
    segTmp->margin =  atof(tokens.at(5).c_str());
    segTmp->delPhi =  atof(tokens.at(6).c_str());
    segTmp->reflect = atoi(tokens.at(7).c_str() );
    segTmp->roughness =  atof(tokens.at(8).c_str());
    segTmp->posErrorX =  atof(tokens.at(9).c_str());
    segTmp->posErrorY =  atof(tokens.at(10).c_str());
    segTmp->posErrorZ =  atof(tokens.at(11).c_str());
    segTmp->rotErrorX =  atof(tokens.at(12).c_str());
    segTmp->rotErrorY =  atof(tokens.at(13).c_str());
    segTmp->rotErrorZ =  atof(tokens.at(14).c_str());
    segTmp->bRead = 1;
  
    if (eMirPS == P1) {
      (opt->vSegP1)[ix -1] = segTmp;
    }
    else if (eMirPS == P2) {
      (opt->vSegP2)[ix -1] = segTmp;
    }
    else if (eMirPS == S1) {
      (opt->vSegS1)[ix -1] = segTmp;
    }
    else if (eMirPS == S2) {
      (opt->vSegS2)[ix -1] = segTmp;
    }
  } 

};

/******************** end of readChangeRecord ****************/
void GReadSegSCStd::getPolyCoeffs() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- GReadSegSCStd::getPolyCoeffs " << endl;
  }

  ifstream inFile(spilotfile.c_str(),ios::in);
  if (! inFile) {
    cerr << "  -- GReadSegSCStd::getReflCoeff " << endl;
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
      opt = SCFac->mStdOptics[index];
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
      opt = SCFac->mStdOptics[index];
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
      itGr =  SCFac->mGRefl->find(index);

      if (itGr!=SCFac->mGRefl->end()) {
        cerr << "trying to load reflection curve " << index << endl;
        cerr << "but curve already in the map: check pilot files " << endl;
        exit(0);
      }

      (*(SCFac->mGRefl))[index] = new TGraph(number);
      TGraph *grTmp = (*(SCFac->mGRefl))[index];

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
   
      *oLog << "       reading in reflectivity curve " << index 
	    << "   with " << number << " points" << endl;
      *oLog << "             converting from nm to cm for robast" << endl;
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
  if (SCFac->mGRefl->size() == 0) {
    *oLog << "no reflection coeffient table found in config file" << endl;
    *oLog << "    STOPPING CODE" << endl;
    exit(0);
  }
  
};
/*
VERSION3.1
2March2015
VERSION with reordered calls to open output file and rootwriter instances
*/
/*!  gropt.cpp
          optics simulations main code
           Charlie Duke
           Grinnell College
           April,2011
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
#include "GPilot.h"

//#include "GTelescopeFactory.h"
//#include "GDCTelescopeFactory.h"
//#include "GSCTelescopeFactory.h"
//#include "GSegSCTelescopeFactory.h"

#include "GSegmentedMirror.h"
#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GSCTelescope.h"
#include "GSegSCTelescope.h"

#include "GTelescopeFactory.h"
#include "GDCTelescopeFactory.h"
#include "GSCTelescopeFactory.h"
#include "GSegSCTelescopeFactory.h"

#include "GReadDCStdBase.h"
#include "GReadDCStdGrISU.h"
#include "GReadSCStd.h"
#include "GReadSegSCStd.h"

#include "GReadPhotonBase.h"
#include "GReadPhotonGrISU.h"
#include "GArrayTel.h"
#include "GSimulateOptics.h"
#include "GRootWriter.h"

TRandom3 TR3;

/*! \brief structure to hold command line entries
 */
struct Cline {
  string pilotfile;    //!< name of pilot file
  enum RdType inType;  //!< input file type enum 
  bool inTypeFlag;     //!< true if inType on command line
  string inFileName;   //!< name of cherenkov input file
  enum OfType outType; //!< output file type enum
  bool outTypeFlag;    //!< true if outType on command line
  string outFileName;  //!< name of output file
  vector<double> vWobble;    //!< wobble entries
};

/*! \brief structure to hold pilot file entries
 */
struct Pilot {
  string pilotfile;     //!< name of pilot file
  enum RdType inType;  //!< input file type enum 
  string inFileName;  //!< name of cherenkov input file
  enum OfType outType;  //!< output file type enum
  string outFileName;  //!< name of output file
  string outFileHeaderTree; //!< name of header tree
  string outFileTelTreeName; //!< base name of telescope tree
  bool outFileDCos;     //!< add DCos branches to outFileTelTree 
  string arrayConfigFile;  //!< name of array configuration file
  string logFileName;      //!< name of log file

  int printStdTelNumber;   //!< print this standard telescope
  int printStdTelMode;     //!< use this print mode for std. tels.

  double wobble[3];  //!< wobblex/y/r all in radians
  double latitude;  //!<  latitude in radians

  int nShower;  //!< 
  int nPhoton;  //!< 

  string photonHistoryFile;  //!< 
  string photonHistoryTree;  //!< 
  UInt_t seed;  //!< 

  int telToDraw;
  int testTel;  //!< telescope number for test graphs (>0). if zero no test produced
  string testTelFile; //!< base filename for test output
  bool debugBranchesFlag; //!< if true, create debug branches in output root file
  unsigned iNInitEvents;
};

/*! \brief structure to hold telescope factory parameters
 */
struct TelFactory {
  TelType telType;  //!< telescope type (DC/SC)
  RdType rdType;    //!< reader type GRISU or CORSIKA
  string configFile; //!< configuration filename
  string editFile;  //!< edit filename
};

struct TelDetails {
  int telID;       //!< number of telescope in the array
  int telStd;      //!< telescope standard number
  TelType telType; //!< telescope type enum
  ROOT::Math::XYZVector telLocGrd;  //!< telescope loc.grdCoordinates

  double telOffSetX; //!< telescope pointing offset, X dir (radians)
  double telOffSetY; //!< telescope pointing offset, Y dir (radians)
  int printMode;     //!< print mode for this telescope
};

#define DEBUGS(x) *oLog << "      " << #x << " = " << x << endl

/*! /brief function to show command-line help
 */
int commandLineHelp();

/*! /brief function to read command-line entries
 */
int readCommandLine(int argc, char *argv[], Cline *cline);

/*! /brief function to print command-line entries
 */
void printCommandLine( Cline *cline);

/*! /brief function to print pilot structure
 */
int pilotPrint(const Pilot &pilot);

/*!  /brief function to read pilot file
 */
int readPilot(Pilot *pilot);

/*!  /brief function to reconcile Cline and Pilot structures
 */
int updatePilot(const Cline &cline,Pilot *pilot);

/*!  /brief function to read arrayConfigFile factory records
 */
int getTelescopeFactoryDetails(vector<TelFactory *> *vTelFac,
                                 map<int,TelDetails *> *mTelDetails,
                                 const string &arrayConfigFile);

int closeRootFiles();

ostream *oLog;

/**************************   main   *****************************/
int main(int argc, char *argv[]) {

  //bool bDrawTel = false;
  bool bDrawTel = true;
  // log file output stream, initialize
  oLog = &cerr;

  string readerFile = "./Config/stdSegSCTelescopes.cfg";
  *oLog << "readerFile " << readerFile << endl;

  *oLog << "making reader " << endl;
  GReadSegSCStd *readerSegSC = new GReadSegSCStd(readerFile);
  *oLog << "making factory " << endl;
  string editFile = "./Config/arrayConfig.cfg";

  TRint *app = 0;
  Int_t pseudo_argc = 1;
  if (bDrawTel) {
    *oLog << "creating app " << endl;
    app = new TRint("app",&pseudo_argc, argv,0,0,kFALSE );
    //runApp = true;
  }
  GSegSCTelescopeFactory *SegSCFac = new GSegSCTelescopeFactory
    (*readerSegSC,editFile);

  int id = 1;
  int std = 1;
  *oLog << "READY to build TELESCOPE" << endl;
  GSegSCTelescope *segTel = SegSCFac->makeTelescope(id,std);
  *oLog << "TELESCOPE built" << endl;
  if (bDrawTel) {
    //segTel->drawTelescope();
    //segTel->getManager()->GetTopVolume()->Draw("ogl");
    app->Run(); 
    //return 0;
  }
 
  SafeDelete(readerSegSC);
  SafeDelete(SegSCFac);
  SafeDelete(segTel);
  return 0;
};
/********************** end of main ***************************************/

int readCommandLine(int argc, char *argv[], Cline *cline) {

  // move arguments into strings, easier to work with strings
  vector<string> clArg;
  for (int i = 0;i<argc;i++) {
    clArg.push_back(argv[i]);
  }

  for (int i = 1;i<argc;i++) {
    *oLog << "COMMANDLINE " << clArg.at(i) << endl;
    if (clArg.at(i)=="-p" ) {
      cline->pilotfile = clArg.at(i+1);
      i++;
    }
    else if (clArg.at(i)=="-it" ) {
      cline->inType = getRdTypeEnum(clArg.at(i+1));
      cline->inTypeFlag = true;
      i++; 
    }
    else if (clArg.at(i)=="-if" ) {
      cline->inFileName = clArg.at(i+1);
      i++;
    } 
    else if (clArg.at(i)=="-p" ) {
      cline->pilotfile = clArg.at(i+1);
      i++;
    }
    else if (clArg.at(i)=="-ot" ) {
      cline->outType = getOfTypeEnum(clArg.at(i+1));
      cline->outTypeFlag = true;
      i++;
    }
    else if (clArg.at(i)=="-of" ) {
      cline->outFileName = clArg.at(i+1);
      i++;
    }
    else if (clArg.at(i)=="-wb" ) {
      for (int j = 0;j<3;j++) {
	double w1 = atof(clArg.at(i+1).c_str());
	w1 = w1*(TMath::DegToRad());
	cline->vWobble.push_back(w1);
	i++;
      }
    }
    else if (clArg.at(i)=="-h" ) {
      commandLineHelp();
      exit(0);
    }
    else {
      *oLog << "   INCORRECT COMMMANDLINE SWITCHES" << endl;
      commandLineHelp();
      exit(0);
    }
  }
  return 1;
};
/********************** end of readCommandLine ************************/

void printCommandLine( Cline *cline ) {

  *oLog << "    -- command line entries; includes defaults" << endl;
  *oLog << "         pilotfile     " << cline->pilotfile << endl;
  *oLog << "         inType/file   " << getRdType(cline->inType) 
       << "   " << cline->inFileName << endl;
  *oLog << "         outType/file  " << getOfType(cline->outType)
       << "   " << cline->outFileName << endl;

  if (cline->vWobble.size() > 0) {
    *oLog << "         vWobble ";
    for ( unsigned i=0; i < cline->vWobble.size(); i++ )  {
      *oLog << cline->vWobble.at(i)*(TMath::RadToDeg()) << "   ";
    }
    *oLog << endl;
  }
};
/******************** end of printCommandLine ***********/
int commandLineHelp() {

  *oLog << "grOptics:  commandline options with default values" << endl;
  *oLog << "    -it <inputFileType = grisu>" << endl;
  *oLog << "        possible types: grisu, corsika, only grisu implemented" 
	<< endl;
  *oLog << "    -if <inputFileName> " << endl;
  *oLog << "    -ot <outputFileType> " << endl;
  *oLog << "        possible types: asci, rootloc, rootpix, not yet implemented " 
	<< endl;
  *oLog << "    -of <outputFileName> " << endl;
  *oLog << "    -p  <pilotFileName = ./Config/opticsSimulation.pilot>  " 
	<< endl;
  *oLog << "    -wb <wobbleX> <wobbleY> <source extension>" << endl;
  return 1;
};
/********************** end of commandLineHelp **************************************/

int pilotPrint(const Pilot &pilot) {

  *oLog << "    -- pilotPrint: defaults included " << endl;
  *oLog << "         pilot.pilotfile    " << pilot.pilotfile << endl;
  *oLog << "         pilot.inType       " << getRdType(pilot.inType) << endl;
  *oLog << "         pilot.inFileName   " << pilot.inFileName << endl;
  *oLog << "         pilot.outType      " << getOfType(pilot.outType) << endl;
  *oLog << "         pilot.outFileName  " << pilot.outFileName << endl;
  *oLog << "         pilot.outFileHeaderTree " << pilot.outFileHeaderTree << endl;
  *oLog << "         pilot.outFileTelTreeName " << pilot.outFileTelTreeName << endl;
  *oLog << "         pilot.outFileDCos  " << pilot.outFileDCos << endl;
  *oLog << "         arrayConfigFile    " << pilot.arrayConfigFile << endl;
  *oLog << "         logFileName        " << pilot.logFileName << endl;
  //*oLog << "         printStdTelNumber  " << pilot.printStdTelNumber << endl;
  //*oLog << "         printStdTelMode    " << pilot.printStdTelMode << endl;
  *oLog << "         nShower / nPhoton  " << pilot.nShower
       << " / " << pilot.nPhoton<< endl;

  *oLog << "         wobble x/y/r       " 
        << pilot.wobble[0]*(TMath::RadToDeg()) << " " 
        << pilot.wobble[1]*(TMath::RadToDeg()) << " "
        << pilot.wobble[2]*(TMath::RadToDeg()) << endl;
  *oLog << "         latitude           " 
        << pilot.latitude*(TMath::RadToDeg()) << endl;
  *oLog << "         photonHistoryFile/Tree " << pilot.photonHistoryFile
       << " / " << pilot.photonHistoryTree << endl;
  *oLog << "         random number seed " << pilot.seed << endl;
  *oLog << "         vector capacities  " << pilot.iNInitEvents << endl;
  *oLog << "         telToDraw " << pilot.telToDraw << endl;
  *oLog << "         testTel   " << pilot.testTel << endl;
  *oLog << "         testTel base filename   " << pilot.testTelFile << endl;
  *oLog << "         debugBranchesFlag       " << pilot.debugBranchesFlag << endl;
  *oLog << endl;

  return 1;
};
/*********************** end of pilotPrint ***************************************/

int readPilot(Pilot *pilot) {

  pilot->inType = GRISU;
  pilot->inFileName = "";
  pilot->outType    = ROOTLOC;
  pilot->outFileName = "";
  pilot->outFileHeaderTree = "";
  pilot->outFileTelTreeName = "";
  pilot->outFileDCos = false;
  pilot->arrayConfigFile = "./Config/arrayConfigSC.cfg";
  pilot->printStdTelNumber = 0;
  pilot->printStdTelMode = 0;
  pilot->nShower = -1;
  pilot->nPhoton = -1;

  for (int i = 0;i<3;++i) {
    pilot->wobble[i] = 0.0;
  }
  pilot->latitude = 0.0;
  pilot->photonHistoryFile = "";
  pilot->photonHistoryTree = "";
  pilot->logFileName       = "";

  pilot->seed = 0;
  pilot->telToDraw = 0;
  pilot->testTel = 0;
  pilot->testTelFile = "";
  pilot->debugBranchesFlag = false;
  pilot->iNInitEvents = 100;
  vector<string> tokens;
  string spilotfile = pilot->pilotfile;

  GPilot *pi = new GPilot(spilotfile);

  string flag = "FILEIN";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
        pilot->inFileName = tokens.at(0);
   if (tokens.size() == 2) {
      pilot->inType = getRdTypeEnum(tokens.at(1));
    }

  }

  flag = "FILEOUT";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->outFileName = tokens.at(0);
    // rootloc type is only implementation for now
    //if (tokens.size() == 2) {
    //pilot->outType = getOfTypeEnum(tokens.at(1));
    //}
    pilot->outFileHeaderTree = tokens.at(1);
    pilot->outFileTelTreeName = tokens.at(2);
    int fileDCos = atoi(tokens.at(3).c_str() );
    if (fileDCos > 0) pilot->outFileDCos = true;
  }
  flag = "ARRAYCONFIG";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->arrayConfigFile = tokens.at(0);
  }
  flag = "PRINTSTDTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->printStdTelNumber = atoi(tokens.at(0).c_str());
    pilot->printStdTelMode   = atoi(tokens.at(1).c_str());
    if (pilot->printStdTelNumber == 0) {
      *oLog << "CANNOT PRINT STANDARD TELESCOPE 0 "
           << " STANDARD TELESCOPE NUMBERING STARTS AT 1"
           << " PROBLEM IN PILOT FILE, PRINTSTDTEL RECORD" << endl;
      *oLog << "    ending code now" << endl;
    }
  }

  flag = "WOBBLE";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->wobble[0] = atof(tokens.at(0).c_str())*(TMath::DegToRad());
    pilot->wobble[1] = atof(tokens.at(1).c_str())*(TMath::DegToRad());
    pilot->wobble[2] = atof(tokens.at(2).c_str())*(TMath::DegToRad());
    pilot->latitude  = atof(tokens.at(3).c_str())*(TMath::DegToRad());
  }  
  flag = "PHOTONHISTORY";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->photonHistoryFile = tokens.at(0);
    if (tokens.size() == 2) {
      pilot->photonHistoryTree = tokens.at(1);
    }
  }  
  flag = "NSHOWER";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->nShower = atoi(tokens.at(0).c_str());
    if (tokens.size() == 2) {
      pilot->nPhoton = atoi(tokens.at(1).c_str());
    }
  }  
  flag = "LOGFILE";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->logFileName = tokens.at(0);
  }  
  flag = "SEED";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->seed = (UInt_t)atoi(tokens.at(0).c_str());
  }
  flag = "VECCAPACITY";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->iNInitEvents = (UInt_t)atoi(tokens.at(0).c_str());
  }
  flag = "DEBUGBRANCHES";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int tmp = (Int_t)atoi(tokens.at(0).c_str());
    if (tmp > 0) pilot->debugBranchesFlag = true;
  }  
  flag = "DRAWTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->telToDraw = atoi(tokens.at(0).c_str());
  }  
  flag = "TESTTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->testTel = atoi(tokens.at(0).c_str());
    if (tokens.size() == 2 ) {
      pilot->testTelFile = tokens.at(1);
    }
  }  
    
  delete pi;
  return 1;
};
/******************************* end of readPilot *****************************/
int updatePilot(const Cline &cline,Pilot *pilot) {

  if (cline.inTypeFlag) {
    pilot->inType = cline.inType;
  }
  if (cline.inFileName != "") {
    pilot->inFileName = cline.inFileName;
  }
  if (cline.outTypeFlag) {
    pilot->outType = cline.outType;
  }
  if (cline.outFileName != "") {
    pilot->outFileName = cline.outFileName;
  }
  if (cline.vWobble.size() > 0 ) {
    if (cline.vWobble.size() != 3) {
      *oLog << "    INCORRECT NUMBER OF WOBBLE ENTRIES ON COMMAND LINE "
	    << endl;
      commandLineHelp();
      exit(0);
    }
    for (unsigned j = 0;j<3;j++) {
      pilot->wobble[j] = cline.vWobble.at(j);
    } 
  } 
  return 1;
};
/****************************** end of updatePilot ******************************/

int getTelescopeFactoryDetails(vector<TelFactory *> *vTelFac,
                                 map<int,TelDetails *> *mTelDetails,
                                 const string &arrayConfigFile) {

  // normally leave debug = true for logfile printout
  bool debug = false;
  if (debug) {
    *oLog << "    -- getTelescopeFactoryDetails" << endl;
  }
  string spilotfile = arrayConfigFile;
  vector<string> tokens;
  
  if (debug) {
    *oLog << "          reading  " << spilotfile  << endl;
  }
  GPilot *pi = new GPilot(spilotfile);
  unsigned ic = 0;

  string flag = "TELFAC";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    vTelFac->push_back(new TelFactory);
    (*vTelFac).at(ic)->telType = getTelTypeEnum(tokens.at(0));
    (*vTelFac).at(ic)->rdType  = getRdTypeEnum(tokens.at(1)); 
    (*vTelFac).at(ic)->configFile = tokens.at(2);
    (*vTelFac).at(ic)->editFile   = tokens.at(3);
    ic++;
  }
  flag = "ARRAYTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int telID = atoi(tokens.at(0).c_str());
    (*mTelDetails)[telID] = new TelDetails();
    (*mTelDetails)[telID]->telID = telID;
    (*mTelDetails)[telID]->telStd = atoi(tokens.at(1).c_str());
    (*mTelDetails)[telID]->telType = getTelTypeEnum(tokens.at(2));
    double xl,yl,zl;
    xl =  atof(tokens.at(3).c_str());
    yl =  atof(tokens.at(4).c_str());
    zl =  atof(tokens.at(5).c_str());
    (*mTelDetails)[telID]->telLocGrd.SetCoordinates(xl,yl,zl);
    double xoffset,yoffset;
    xoffset = atof(tokens.at(6).c_str())*(TMath::DegToRad());
    yoffset = atof(tokens.at(7).c_str())*(TMath::DegToRad());
    (*mTelDetails)[telID]->telOffSetX = xoffset;
    (*mTelDetails)[telID]->telOffSetY = yoffset;
    (*mTelDetails)[telID]->printMode = atoi(tokens[8].c_str());
  }
  delete pi;

  if (debug) {
    *oLog << "       Telescope factories to instantiate" << endl;
    for (unsigned i = 0;i < vTelFac->size();++i) {
      *oLog << "         factory number " << i+1 << endl;
      *oLog << "         factory type/reader type  " 
            << getTelType((*vTelFac).at(i)->telType)
            << "  " << getRdType((*vTelFac).at(i)->rdType) << endl;
      *oLog << "         configuration file        " 
            << (*vTelFac).at(i)->configFile << endl;
      *oLog << "         edit file                 " 
            << (*vTelFac).at(i)->editFile << endl<< endl;
    }
 
    // print telescope details structure
    *oLog << "       array telescope details " << endl;
    *oLog << "                    ";
    *oLog << "type std   locX     locY    locZ  osX  osY  prtMode" << endl;
    map<int,TelDetails *>::iterator mIter;    
    
    for (mIter=mTelDetails->begin();mIter!=mTelDetails->end();mIter++) {
      *oLog << "         telNbr= " << mIter->first << "   ";
      *oLog  << getTelType(mIter->second->telType) << "   ";
      *oLog  << mIter->second->telStd << "   " 
	     << mIter->second->telLocGrd.X() 
	     << "    " << mIter->second->telLocGrd.Y() <<  "   " 
	     << mIter->second->telLocGrd.Z();
      *oLog  << "   " 
	     << mIter->second->telOffSetX*(TMath::RadToDeg()) 
	     << "     " 
	     << mIter->second->telOffSetY*(TMath::RadToDeg()) 
	     << "   " 
	     << mIter->second->printMode << endl;
    }
  }
  return 1;
};
/******************************** end of getTelescopeFactoryDetails *****/

int closeRootFiles() {
  // procedure to close all open root files, useful for debugging
  *oLog << "        closeRootFiles " << endl;
  int icount=0;
  TIter next( gROOT->GetListOfFiles() );

  TFile *fi = 0;
  while ( (fi = (TFile *)next() ) ) {
    icount++;
    *oLog << "            closing root file if non-zero: " 
	  << fi->GetName() << endl;
    if (!fi) {
      fi->Close();
    }
  }

  return icount;
}
/******************************** end of closeRootFiles *****/

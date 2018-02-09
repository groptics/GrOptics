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

#include "ADefinition.h"
#include "GUtilityFuncts.h"
#include "GPilot.h"

#include "ATelescope.h"
#include "GDCTelescope.h"
#include "GSCTelescope.h"

#include "ATelescopeFactory.h"
#include "GDCTelescopeFactory.h"
#include "GSCTelescopeFactory.h"

#include "GReadDCStdBase.h"
#include "GReadDCStdGrISU.h"
#include "GReadSCStd.h"

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

  // log file output stream, initialize
  oLog = &cerr;

  // TRint for possible use later with 
  // drawing a telescope. Instantiate later only if necessary
  int pseudo_argc = 1;
  TRint *app = 0;
  app = new TRint("app",&pseudo_argc, argv,0,0,kFALSE );

  time_t startTime = time(NULL);
  clock_t startClock = clock();

  // write start date/time to stdout
  ostringstream ostr;
  ostr << "============================================" << endl;
  ostr << endl;
  ostr << "  -- grOptics main program: version " << VERSN << endl;
  ostr << "       " << ctime(&startTime) << endl;
  ostr << "============================================" << endl;
 
  *oLog << ostr.str();
  
  ////////////////////////////////////////////////////////////
  /////// read command line and pilot file 
  /////// set random number generator seed
  /////// open and write to log file
  ////////////////////////////////////////////////////////////
  
  Pilot pilot;
  Cline cline;

  //set command line defaults
  cline.pilotfile = "./Config/opticsSimulation.pilot"; // default value
  //cline.pilotfile = "./Config/opticsSimulation.pilot"; // default value
  cline.inFileName    = "";
  cline.inType        = GRISU;
  cline.inTypeFlag    = false;
  cline.outFileName   = "";
  cline.outType       = ROOTLOC;
  cline.outTypeFlag   = false;
  
  // command line printed to cerr in this function
  readCommandLine(argc,argv,&cline);

  // set pilotfile name, read pilot and update from command line
  // structure
  pilot.pilotfile = cline.pilotfile;
  readPilot(&pilot);
  updatePilot(cline,&pilot);

  // set seed so can print out seed if from machine clock
  TR3.SetSeed(pilot.seed); 

  if (pilot.logFileName == "") {
    oLog = &cerr;
  }
  else {

    oLog = new ofstream( (pilot.logFileName).c_str() );
    // print heading from above
    *oLog << ostr.str();
 
   // redirect root output to log file.
    gSystem->RedirectOutput( (pilot.logFileName).c_str() ); 
  }

  // print command line entries to logfile
  printCommandLine(&cline);
  // print pilot file details to logfile
  *oLog << "     pilot structure after updating from command line entries" 
	<< endl;
  pilotPrint(pilot);
  if (pilot.seed == 0) {
    UInt_t tr3Seed = TR3.GetSeed();
    *oLog << "         seed set by machineClock to " << tr3Seed << endl;
  }
  
  *oLog << endl 
	<< "============= finished read command line and pilot file" 
	<< endl << endl;

  ////////////////////////////////////////////////////////////
  /////// fill map (mTelDetails) of TelDetails structure
  /////// fill vector (vTelFac) of GTelescopeFactory's 
  /////// use factories to fill map (mArrayTel) of array telescopes
  ////////////////////////////////////////////////////////////
 
  // get listings of required telescope factories from 
  // arrayConfiguration file. 
  // get map of telescope details

  if (pilot.telToDraw > 0) {
    app = new TRint("app",&pseudo_argc, argv,0,0 );
  }
  map<int,TelDetails *> mTelDetails;
  map<int,TelDetails *>::iterator mIter; 

  vector< TelFactory* > vTelFac;

  getTelescopeFactoryDetails(&vTelFac,&mTelDetails,pilot.arrayConfigFile);

  // set up factories
  ATelescopeFactory *DCFac = 0;
  ATelescopeFactory *SCFac = 0;
  GReadDCStdBase *readerDC = 0;
  GReadSCStd *readerSC = 0;

  for (unsigned i = 0;i< vTelFac.size(); ++i) {
    if (vTelFac.at(i)->telType == DC) {
      if (vTelFac.at(i)->rdType == GRISU) {
        // set up GRISU reader, make this a switch statement later
        readerDC = new GReadDCStdGrISU(vTelFac.at(i)->configFile);
	// the factory now owns the readerDC
        DCFac = new GDCTelescopeFactory(*readerDC,vTelFac.at(i)->editFile);
      }
    }
    else if ( vTelFac.at(i)->telType == SC) {
      readerSC = new GReadSCStd(vTelFac.at(i)->configFile);
      SCFac = new GSCTelescopeFactory(*readerSC,vTelFac.at(i)->editFile);
    }
  }

  // make arrayTelescope map<telid, ArrayTel*>
  map<int, GArrayTel *> mArrayTel;
  map<int, GArrayTel *>::iterator mArrayTelIter;

  int telnumm = 0;

  for (mIter=mTelDetails.begin();mIter!=mTelDetails.end();mIter++) {
    telnumm++;
    int telId = mIter->first;
    int telStd = mIter->second->telStd;
    TelType telType = mIter->second->telType;
    int printMode = mIter->second->printMode;

    ROOT::Math::XYZVector telLocGrd;
    telLocGrd = mIter->second->telLocGrd;

    double xoffsettel, yoffsettel;
    xoffsettel = mIter->second->telOffSetX;
    yoffsettel = mIter->second->telOffSetY;
    if (telType==DC) {
      ATelescope *tel = DCFac->makeTelescope(telId,telStd);
      if (pilot.telToDraw == telId) {
	app = new TRint("app",&pseudo_argc, argv,0,0,kFALSE );
	tel->drawTelescope();
	app->Run(); 
	return 0;
      }
      mArrayTel[telId] = new GArrayTel(telLocGrd,xoffsettel,
                                       yoffsettel,telType,
                                       telId,telStd,printMode,tel);

      if (pilot.photonHistoryFile != "") {
        mArrayTel[telId]->setPhotonHistory(pilot.photonHistoryFile,
                                           pilot.photonHistoryTree);
      }                                               
    }
    else if (telType==SC) {

      ATelescope *tel = SCFac->makeTelescope(telId,telStd);
      if (pilot.telToDraw == telId) {
	app = new TRint("app",&pseudo_argc, argv,0,0,kFALSE );
	tel->drawTelescope();
	app->Run(); 
	return 0;
      }
      mArrayTel[telId] = new GArrayTel(telLocGrd,xoffsettel,
                                       yoffsettel,telType,
                                       telId,telStd,printMode,tel);

      if (pilot.photonHistoryFile != "") {
        mArrayTel[telId]->setPhotonHistory(pilot.photonHistoryFile,
                                           pilot.photonHistoryTree);
      }

    }
  }
  *oLog << endl
	<< "============= finished factories and telescope making" 
	<< endl << endl;

  ////////////////////////////////////////////////////////////
  /////// make a photon reader (GReadPhotonBase)
  /////// create an instance of GSimulateOptics for the raytracing
  ////////////////////////////////////////////////////////////
  
  // make a photon reader ////////////////////////////
   
  GReadPhotonBase *readP = 0;
  if (pilot.inType==GRISU) {
    readP = new GReadPhotonGrISU();    
    readP->setInputFile(pilot.inFileName);
  }
  else {
    *oLog << " can't open reader type: " << pilot.inType << endl;
    *oLog << " stopping code, check pilot file" << endl;
    exit(0);
  }

  // open root output file and write header as a TString to the file,
  // trees are filled later in GSimulateOptics
  
  string outFileName = pilot.outFileName;  //!< name of output file
  
  TFile *fO = new TFile(outFileName.c_str(),"RECREATE");
  if (fO->IsZombie() ) {
    *oLog << "error opening root output file: " << outFileName << endl;
    *oLog << "...exiting" << endl;
    exit(-1);
  }
  
  // make map of photon writers, one per telescope, using telID
  // as key
  map<int,GRootWriter *> mRootWriter;
  map<int,GRootWriter *>::iterator mRootWriterIter;
  for (mIter=mTelDetails.begin();mIter!=mTelDetails.end();mIter++) {
    int telID = mIter->first;
    mRootWriter[telID] = new GRootWriter(fO,telID,pilot.outFileTelTreeName,
					 pilot.outFileDCos);
  }
  
  ////////////////////////////////////////
  /////// instantiate the simulate optics class (GSimulateOptics)

  ////////////////////////////////////////////////////////////
  
  GSimulateOptics *siO = new GSimulateOptics(readP,&mArrayTel,
					     &mRootWriter,pilot.outFileHeaderTree);
  siO->setWobble(pilot.wobble[0],pilot.wobble[1],
		 pilot.wobble[2],pilot.latitude);
  
  /////////////////////////////////////////////////////////////
  /////// do the simulations (where do we create the output class).
  ////////////////////////////////////////////////////////////
  *oLog << "================= starting simulations " << endl << endl;
  
  siO->startSimulations(pilot.nShower,pilot.nPhoton);
  
  *oLog << "======================================== " << endl;
  
  // have to loop thru the telescopes here since only one
  // telescope selected by last photon in the file.
  // Could write the history file in the telescope
  // destructor.
  for (mArrayTelIter=mArrayTel.begin();
       mArrayTelIter!=mArrayTel.end();
       mArrayTelIter++) {
    if (pilot.photonHistoryFile != "") {
      mArrayTelIter->second->writePhotonHistory();
    }
    //delete mArrayTelIter->second; 
  }
  
  // Close() probably deletes the TFile; use SafeDelete anyway
  fO->cd();
  fO->Close();
  SafeDelete(fO);
  
  clock_t endClock = clock();
  double elapsedTime = (double)((endClock - startClock)/ 
				(double)CLOCKS_PER_SEC);
  fprintf(stderr,"elapsed time %12.6f\n",elapsedTime);
  
  time_t endTime = time(NULL);
  fprintf(stderr,"---- RUN TIME %12.6f seconds\n\n",
          difftime(endTime,startTime));
  
  *oLog << "**** ALL DONE:  " << ctime(&endTime) << endl;
  
  // uncomment this Run statement if wish to remain in ROOT 
  // after end of execution.
  app->Run(); 

  // clean up
  if (readP) SafeDelete(readP);  // delete photon reader

  for (mIter = mTelDetails.begin(); mIter!=mTelDetails.end();
       mIter++) {
    SafeDelete(mIter->second);
  }
  for (unsigned ij = 0;ij < vTelFac.size();ij++) {
    SafeDelete(vTelFac[ij]);
  }
  if (DCFac) SafeDelete(DCFac);
  if (SCFac) SafeDelete(SCFac);

  // this deletes the telescopes
  for (mArrayTelIter = mArrayTel.begin();
       mArrayTelIter != mArrayTel.end();
       mArrayTelIter++ ) {
    if (mArrayTelIter->second) SafeDelete(mArrayTelIter->second);
  }
  // delete root writer map entries
  for (mRootWriterIter = mRootWriter.begin();
       mRootWriterIter != mRootWriter.end();
       mRootWriterIter++) {
    if (mRootWriterIter->second) SafeDelete(mRootWriterIter->second);
  } 
  if (siO) SafeDelete(siO);

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

};
/******************** end of printCommandLine ***********/
int commandLineHelp() {

  *oLog << "gropt:  commandline options with default values" << endl;
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
  *oLog << "         telToDraw " << pilot.telToDraw << endl;
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
      exit(0);
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
  flag = "DRAWTEL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->telToDraw = atoi(tokens.at(0).c_str());
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

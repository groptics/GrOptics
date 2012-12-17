/*
VERSION2.7aBeta
17Dec2012
*/
/*!  GSegSCTelescopeFactory.cpp

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
#include <assert.h>

#include "TROOT.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"

#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "AOpticsManager.h" 
#include "GTelescope.h"
#include "GSegSCTelescope.h"
#include "GSCTelescope.h"
#include "GTelescopeFactory.h"
#include "GSCTelescopeFactory.h"
#include "GSegSCTelescopeFactory.h"

#include "GReadSegSCStd.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

/*!  /brief SegSCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
SegSCStdOptics::SegSCStdOptics() {

  
};
/************** end of SegSCStdOptics ***********************/

SegSCStdOptics::SegSCStdOptics(const SegSCStdOptics &sco) {

  bool debug = true;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::SegSCStdOptics " << endl;
  }

};

/************** end of SegSCStdOptics ***********************/
SegSCStdOptics::~SegSCStdOptics() {
  bool debug = true;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::~SegSCStdOptics" << endl;
  }
};
/********************* end of ~SegSCStdOptics *****************/

void SegSCStdOptics::printSegSCStdOptics() {

  bool debug = true;
  if (debug) {
    *oLog << "  -- SegSCStdOptics::printSegSCStdOptics " << endl;
  }
 
};
/************** end of printSegSCStdOptics ***********************/

//////////////////////////////////////////////////////////////////////////

GSegSCTelescopeFactory::
GSegSCTelescopeFactory(GReadSegSCStd &scReader,
		    const string &editPilotFile) {

  bool debug = true;
  sPilotEdit = editPilotFile;
  if (debug) {
    *oLog << "  -- SCTelescopeFactory Constructor:  " << sPilotEdit << endl;
  }
  
};

/************** end of GSegSCTelescopeFactory ***********************/

GSegSCTelescopeFactory::~GSegSCTelescopeFactory() {
 
  bool debug = true;
  if (debug) {
    *oLog << "  -- GSegSCTelescopeFactory::~GSegSCTelescopeFactory" << endl;
  }
  
};
/************** end of ~GSegSCTelescopeFactory ***********************/

GSegSCTelescope* GSegSCTelescopeFactory::makeTelescope(const int &id,
                                                     const int &std) {
  
  int debug = 1;
  if (debug>0) {
    *oLog << " -- GSegSCTelescopeFactory::makeTelescope" << endl;
    //*oLog << "      telID  = " << id << endl;
    //*oLog << "      telStd = " << std << endl;
  }
  SCTel = new GSegSCTelescope;
  return SCTel;
 
};
/************** end of makeTelescope ***********************/

void GSegSCTelescopeFactory::editWorkingTelescope(GSegSCTelescope *SCTel1) {
  int debug = 1;
  if (debug>0) {
    *oLog << " -- GSegSCTelescopeFactory::editWorkingTelescope" << endl;
  }

 
};
/************** end of editWorkingTelescope ***********************/

void GSegSCTelescopeFactory::printStdTelescope(const int &iStd, 
                                               const int &mode,ostream &oStr) {

  // DO NOT USE.
  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GSegSCTelescopeFactory::printStdTelescope" << endl;
  }
  
};
/************** end of :printStdTelescope ***********************/

void GSegSCTelescopeFactory::setPrintMode(ostream &oStr,
                                          const int prtMode) {

  int debug = 1;
  if (debug > 0) {
    *oLog << " -- GSegSCTelescopeFactory::setPrintMode" << endl;
  }
 
}; 
/************** end of setPrintMode  ***********************/




/*
VERSION4.0
30May2016
*/
/* ATelescopeFactory.cpp

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

using namespace std;
// C.Duke
//#include "ADefinition.h"
#include "ATelescope.h"
#include "ATelescopeFactory.h"

//ClassImp(ATelescopeFactory)

ATelescopeFactory::ATelescopeFactory() {

  // initialize parameters
  // C.Duke
  //oPrtStrm = oLog;
  iPrtMode = 0;

};
/**************** end of ATelescopeFactory ****************/


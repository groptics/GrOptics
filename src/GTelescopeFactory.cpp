/*
VERSION3.1
2March2015
*/
/* GTelescopeFactory.cpp

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

#include "GDefinition.h"
#include "GTelescope.h"
#include "GTelescopeFactory.h"

GTelescopeFactory::GTelescopeFactory() {

  // initialize parameters
  oPrtStrm = oLog;
  iPrtMode = 0;

};
/**************** end of GTelescopeFactory ****************/


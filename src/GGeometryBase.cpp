/*
VERSION3.1
2March2015
*/
/*    GGeometryBase.cpp

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

#include "GGeometryBase.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl

GGeometryBase::GGeometryBase() {
  oPrtStrm = &*oLog;
  iPrtMode = 0;

};
/***************** end of GGeometryBase ********************/

GGeometryBase::~GGeometryBase() {

};
/***************** end of ~GGeometryBase ********************/




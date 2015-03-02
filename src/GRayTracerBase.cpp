/*
VERSION3.1
2March2015
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

#include "GRayTracerBase.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl

GRayTracerBase::~GRayTracerBase() {

  oPrtStrm = oLog;
  iPrtMode = 0;

};

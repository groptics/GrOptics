/*
VERSION4.0
30May2016
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
#include "GUtilityFuncts.h"
#include "GPilot.h"

#include "ATelescope.h"
#include "GDCTelescope.h"

#include "ATelescopeFactory.h"
#include "GDCTelescopeFactory.h"

#include "GReadDCStdBase.h"
#include "GReadDCStdGrISU.h"

#include "GRayTracerBase.h"

GReadDCStdBase::GReadDCStdBase() {

  oPrtStrm = oLog;
  iPrtMode = 0;

  DCFac = 0;
};



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

using namespace std;

#include "AOpticsManager.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TView.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVector3.h"
#include "TGeoPgon.h"
#include "TRandom3.h"

#include "ADefinition.h"
#include "ATelescope.h"

//#include "GGeometryBase.h"
///#include "GDCGeometry.h"

//#include "GRayTracerBase.h"
//ClassImp(ATelescope)
ATelescope::ATelescope() {
  oPrtStrm = oLog;
  iPrtMode = 0;


};


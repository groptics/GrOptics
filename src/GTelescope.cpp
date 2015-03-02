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

#include "GDefinition.h"
#include "GTelescope.h"

#include "GGeometryBase.h"
#include "GDCGeometry.h"

#include "GRayTracerBase.h"

GTelescope::GTelescope() {
  oPrtStrm = oLog;
  iPrtMode = 0;


};


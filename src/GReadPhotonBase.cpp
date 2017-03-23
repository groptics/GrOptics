/*
VERSION4.0
30May2016
*/
/*!  GReadPhotonBase.cpp
     Charlie Duke
     Grinnell College
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

#include "ADefinition.h"
#include "GReadPhotonBase.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl

GReadPhotonBase::GReadPhotonBase() {
  bool debugB = false;
  if (debugB) {
    *oLog << "  -- GReadPhotonBase::GReadPhotonBase" << endl;
  }
};
/******************* GReadPhotonBase ****************************/

GReadPhotonBase::~GReadPhotonBase() {

};

/*
VERSION3.1
2March2015
*/
/* GDCGeometry.cpp

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
#include "GDCGeometry.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGW(x) *oLog << "         " << #x << " = " << x << endl
#define DEBUGST(x)  *oLog << "   " << #x << " = " << x 

GDCGeometry::GDCGeometry() {

  type = DCVERITAS;

  for (int i=0; i<3; i++) {
    focBoxDim[i] = 0.0;
    focBoxRot[i] = 0.0;
  }

  for (int i = 0;i<4;i++) {
    edgeX[i] = 0.0;
    edgeY[i] = 0.0;
    edgeZ[i] = 0.0;
    edgeOffset[i] = 0.0;
    edgeRot1[i] = 0.0;
    edgeRot2[i] = 0.0;
    edgeRot3[i] = 0.0;
    quadArmBottomX[i] = 0.0;
    quadArmBottomY[i] = 0.0;
  }
  shutterX = 0.0;
  shutterZ = 0.0;
  shutterRot1 = 0.0;
  shutterRot2 = 0.0;
  shutterRot3 = 0.0;

  // quad arms
  quadArmX = 0.0;
  quadArmY = 0.0;
  quadArmOffset = 0.0;

  for (int i=0;i<2;i++) {
    crossBarX[i] = 0.0;
    crossBarY[i] = 0.0;
    crossBarDistBelowFocBox[i] = 0.0;
  }
  // other details
  cameraRadius = 0.0;

};
/************************* end of GDCGeometry **********************/

GDCGeometry::~GDCGeometry() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GDCGeometry::~GDCGeometry " << endl;
  }
};
/********************** end of ~GDCGeometry ****************/

GDCGeometry::GDCGeometry(const GDCGeometry &c) {

  epsil = c.epsil;

  for (int i = 0;i<4;i++) {
    edgeX[i] = c.edgeX[i];
    edgeY[i] = c.edgeY[i];
    edgeZ[i] = c.edgeZ[i];
    edgeOffset[i] = c.edgeOffset[i];
    edgeRot1[i] = c.edgeRot1[i];
    edgeRot2[i] = c.edgeRot2[i];
    edgeRot3[i] = c.edgeRot3[i];
    quadArmBottomX[i] = c.quadArmBottomX[i];
    quadArmBottomY[i] = c.quadArmBottomY[i];
  }

  shutterX = c.shutterX;
  shutterZ = c.shutterZ;
  shutterRot1 = c.shutterRot1;
  shutterRot2 = c.shutterRot2;
  shutterRot3 = c.shutterRot3;

  // quad arms
  quadArmX = c.quadArmX;
  quadArmY = c.quadArmY;
  quadArmOffset = c.quadArmOffset;

 for (int i=0;i<2;i++) {
    crossBarX[i] = c. crossBarX[i];
    crossBarY[i] = c. crossBarY[i];
    crossBarDistBelowFocBox[i] = c.crossBarDistBelowFocBox[i];
 }
 for (int i=0; i<3; i++) {
   focBoxDim[i] = c.focBoxDim[i];
   focBoxRot[i] = c.focBoxRot[i];
 }
 
 type = c.type;
 cameraRadius = c.cameraRadius;

};
/************************* end of GDCGeometry **********************/

GeoType GDCGeometry::getType() {

  return type;
};
/************************* end of GDCGeometry **********************/

void GDCGeometry::printGeometry() {
  int debug = 1;
  if (debug > 0) {
    *oLog << " -- GDCGeometry::printGeometry" << endl;
  }
  if (type == NOSTRUCT) {
    *oLog << "       geometry is NOSTRUCT type"<< endl;
    *oLog << "       nothing to print " << endl;
    return;
  }

  for (int i = 0;i<4;i++) {
    *oLog << "         " << i;
    DEBUGST(edgeX[i]);DEBUGST(edgeY[i]);
    DEBUGST(edgeZ[i]); 
    *oLog << endl;
  }
  for (int i = 0;i<4;i++) {
    *oLog << "         " << i ; 
    DEBUGST(edgeOffset[i]); 
    *oLog << endl;

  }
  for (int i = 0;i<4;i++) {
    *oLog << "         " << i ;
    DEBUGST(edgeRot1[i]);DEBUGST(edgeRot2[i]);
    DEBUGST(edgeRot3[i]);
    *oLog << endl;

  }
  for (int i = 0;i<4;i++) {
    *oLog << "         " << i ;
    DEBUGST(quadArmBottomX[i]);DEBUGST(quadArmBottomY[i]);
    *oLog << endl;
  }
    DEBUGW(shutterRot1);DEBUGW(shutterRot2);
    DEBUGW(shutterRot3);
    DEBUGW(quadArmX);DEBUGW(quadArmY);
    DEBUGW(quadArmOffset);
    DEBUGW(cameraRadius);

    *oLog << " -- End of printGeometry" << endl << endl;
};
/************************* end of printGeometry **********************/

void GDCGeometry::setPrintMode(ostream &oStr,
                                          const int prtMode) {
  int debug = 1;
  if (debug > 0) {
    *oLog << " -- GDCTelescopeFactory::setPrintMode" << endl;
  }
  oPrtStrm = &oStr;
  iPrtMode = prtMode;
}; 
/************************* end of setPrintMode **********************/

void GDCGeometry::printGeometry1(ostream &oStr,const int &prtMode) {
  oStr << " -- GDCGeometry::printGeometry1" << endl;
  oStr << "       edgeX/Y/Z [4] " << endl;
  for (int i = 0;i<4;i++) {
    oStr << "         " << i;
    oStr << "            " << edgeX[i] << " " << edgeY[i]
         << " " << edgeZ[i] << endl;
  }
  oStr << "       edgeOffset [4] " << endl;
  for (int i = 0;i<4;i++) {
    oStr << "         " << i << "  " << edgeOffset[i] << endl; 
  }
  oStr << "       edgeRot1/2/3 [4] " << endl;
  for (int i = 0;i<4;i++) {
    oStr << "         " << i << " " 
         << edgeRot1[i] << " " << edgeRot2[i] << " " 
         << edgeRot3[i] << endl;

  }
  oStr << "       quadArmBottomX/Y [4] " << endl;
  for (int i = 0;i<4;i++) {
    oStr << "         " << i << " " << quadArmBottomX[i] << " "
         << " " << quadArmBottomY[i] << endl;
  }

  oStr << "      crossArms, x,y,distance below focus box " << endl;
  for (int i = 0;i<2;i++) {
    oStr << "         " << i << "  " << crossBarX[i] << "  " 
	 << crossBarY[i] << "  " << crossBarDistBelowFocBox[i] << endl;
  }
  oStr << "       shutterRot1/2/3 " << endl;
  oStr << "         " << shutterRot1 << " " 
       << shutterRot2 << " " << shutterRot3 << endl;
  oStr << "       quadArmX/Y" << endl;
  oStr << "         " << quadArmX << " " << quadArmY << endl;
  oStr << "       quadArmOffset " << quadArmOffset << endl;
  oStr << "       end of printGeometry" << endl;
  oStr << endl;

};
/************************* end of printGeometry1 **********************/


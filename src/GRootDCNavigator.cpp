/*
VERSION3.1
2March2015
*/
// standard includes
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

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TMatrixD.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TView.h"
#include "TPad.h"
#include "TVector3.h"
#include "TGeoPgon.h"

//#include ".h"
#include "GDefinition.h"
#include "GUtilityFuncts.h"
#include "GGeometryBase.h"
#include "GDCGeometry.h"

#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GRootDCNavigator.h"

ClassImp(GRootDCNavigator);

// debug macro
#define DEBUG(x) *oLog << #x << " = " << x << endl

//****************************************************

GRootDCNavigator:: GRootDCNavigator(GDCTelescope *dcTel,
                                    int makeFacetsFlag, bool debugT) {

  DCTel = dcTel;
  fDebugT = debugT;  
  fMakeFacetsFlag = makeFacetsFlag;

  if (fDebugT) {
    *oLog << "  -- GRootDCNavigator::GRootDCNavigator" << endl;
    *oLog << "       fMakeFacetsFlag: " << fMakeFacetsFlag << endl;
  }
  // set pointers to zero
  gDC = 0;
  fGeom = 0;
  fMatVacuum = 0;
  fMatAl = 0;
  fVacuum = 0;
  fAl = 0;
  fTopPosV = 0;
  fTopVol = 0;

  // initialize single variables
  fFL = 0.0;
  fnodeNum = 0;
  fepsil = 0.0; 
  fFocBoxZTop = 0.0;
  fMoveToTop = false;
  fNodeC = 0;
  fNodeN = 0;

  // initialize arrays
  for (int i=0;i<3;i++) {
    fTopDim[i] = 0.0;
    fFocusBoxDim[i] = 0.0;
    fFocusBoxRot[i] = 0.0;
    fQuadArmPosV[i] = 0;
    fQuadArmR2R1V[i] = 0;
  }
  fQuadArmPosV[3] = 0;
  fQuadArmR2R1V[3] = 0;

  initialize();
  setupGeometry();

  // turn on via setTrackingDebug
  fDebugTr = false;

  // take this out at some point
  if (0) {
  // ============ set up default veritas telescope
    fFL = 12.0;            // focal length
    
    // default focus box dimensions (edge to edge, meters)
    // have to set focus box dimensions prior to setting
    // top dimensions
    // normally 1.5
    fFocusBoxDim[0] = 1.5;    // set focus box x dimension
    fFocusBoxDim[1] = 1.5;    // set focus box y dimension
    fFocusBoxDim[2] = 1.02;    // set focus box z dimension
    
    // focus box rotation (thru 45 degrees about orig. z axis)
    //fFocBoxRot[0] = 45.0;
    //fFocBoxRot[1] =  0.0;
    //fFocBoxRot[2] =  0.0;
    
    fepsil = 1.00;
    // Top volume dimensions (+-x,+-y,+-z) meters
    fTopDim[0] = 10.0;  // 6.5 covers all facets. set at 10.0 for now
    fTopDim[1] = 10.0;
    fTopDim[2] = (fFL + fFocusBoxDim[2] + fepsil)/2.0; 
    

    // location of Top center, in telescope coordinates
    //Double_t topx = 0.0;
    //Double_t topy = 0.0;
    //Double_t topz = (fFL + fFocusBoxDim[2] + fepsil)/2.0;
    //fTopPosV = new TVector3(topx,topy,topz);
  }

  if (fDebugT) {
    *oLog << "       fepsil: " << fepsil << endl;
    *oLog << "       FocusBoxDim, side length: " << fFocusBoxDim[0] << " " 
          << fFocusBoxDim[1] << " " << fFocusBoxDim[2] << endl;
    *oLog << "       Focal Length; " << fFL << endl;
    
    *oLog << "       TopVolume half-side length : " << fTopDim[0]  << " " 
          << fTopDim[1] << " " << fTopDim[2] << endl;
    *oLog << "       fTopPosV: top volume vector, TopCenter tele.coor." 
          << endl;
    fTopPosV->Print();
    double tmp1 = (*fTopPosV)[2] - fepsil - fFocusBoxDim[2];
    *oLog << "       Focal Point to TopVol center: " << tmp1 
          << endl; 
    *oLog << "       tele.base center to TopVol center: " 
          << (*fTopPosV)[2] << endl;
    *oLog << " *************End of constructor printing****************" 
          << endl;
    
  }
  // this creates a new TGeoManager for each telescope. 
  // If facets are added to the Geometry, then individual
  // geometries are required because of differing blur radii.  
  // However, at this time, root doesn't permit multiple
  // geometries. So, I'll leave the code like this for now.
  // but only the last geometry created will be open and will
  // serve for all telescopes. This will work fine for veritas
  // as long as the telescopes are identical.

  ostringstream os;
  os << "geomTel" << DCTel->iTelID;
  //*oLog << "geomTel " << os.str() << endl;
  string nameGeom = os.str();

  gGeoManager = 0;
  //*oLog << endl << " ********** new GeoManager ********" << endl << endl;
  //*oLog << "numGeoMan before " << gROOT->GetListOfGeometries()->GetSize() << endl;
  fGeom = new TGeoManager(nameGeom.c_str(), nameGeom.c_str());
  //*oLog << "numGeoMan after " << gROOT->GetListOfGeometries()->GetSize() << endl;

  gGeoManager = fGeom;

  makeMaterialMedia();
  makeTop();
  makeFocusBox();
  makeEdgeBoxes();
  makeQuadArms();
  makeCrossArms();
  makeShutter();
  if (fMakeFacetsFlag) {
    makeFacets();
  }
  //fGeom->SetVisLevel(4);
  
  //--- close the geometry
  fGeom->CloseGeometry();
  //fTopVol->Draw("ogl");
   
   //drawTelescope();
};
//****************************************************

GRootDCNavigator::  ~GRootDCNavigator() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootDCNavigator::~GRootDCNavigator" << endl;
  }
  // fGeom owns all the geometry classes (from root)
  if (fGeom) {
    gGeoManager = fGeom;
    SafeDelete(fGeom);
    // gGeoManager is set to zero ok

  }
  
  SafeDelete(fTopPosV);
  
  for (int i = 0; i<3;i++) {
    SafeDelete(fQuadArmR2R1V[i]);
    SafeDelete(fQuadArmPosV[i]);
  }

};
//****************************************************

void GRootDCNavigator::setupGeometry() {

  bool debug = fDebugT;
  if (debug) {
    *oLog << "  -- GRootDCNavigator::setupGeometry" << endl;
  }

  fFL = DCTel->dFocLgt;

  gDC = dynamic_cast<GDCGeometry*>(DCTel->geoStruct);

  // move structural parameters from gDC to here
  fepsil = gDC->epsil;

  // get focus box dimensions
  for (int i = 0;i<3;i++) {
    fFocusBoxDim[i] = gDC->focBoxDim[i];
    fFocusBoxRot[i] = gDC->focBoxRot[i];
    if (debug) *oLog << "        i fFocusBoxDim  " << i << endl;
  }   
  for (int i = 0;i<3;i++) {
    if (debug) *oLog << "   i fFocusBoxRot " << i << " " << fFocusBoxRot[i] << endl;
  }
 
  // calculate top volume
  //  add 5.0 meters to x/y dimensions to ensure that can move
  //  all photons to top of topVol, else get seg. fault in tracker
  fTopDim[0] = ( (fFL + fepsil)/2.0 ) + 5.0;
  fTopDim[1] = ( (fFL + fepsil)/2.0 ) + 5.0;
  fTopDim[2] = (fFL + fFocusBoxDim[2] + fepsil)/2.0; 

  for (int i = 0;i<3;i++) {
    if (debug) *oLog << "        i fTopDim  " << i << " " << fTopDim[i] << endl;
  }   

  Double_t topx = 0.0;
  Double_t topy = 0.0;
  Double_t topz = (fFL + fFocusBoxDim[2] + fepsil)/2.0;
  fTopPosV = new TVector3(topx,topy,topz);
  for (int i = 0;i<3;i++) {
    if (debug) *oLog << "        i fTopPosV  " << i << " " << (*fTopPosV)[i] << endl;
  }   
  if (debug) {
    *oLog << "      fFL " << fFL << endl;
    *oLog << "      fepsil " << fepsil << endl;
    *oLog << "  LEAVING SETUP GEOMETRY" << endl;
   }
  
};
//****************************************************

void GRootDCNavigator::initialize() {

  fnodeNum = 1;
  // initialize photon tracking variables
  for (int i = 0;i<3;i++) {
    fPosC[i] = 0.0;
    fDirC[i] = 0.0;
    fPosN[i] = 0.0;
    fDirN[i] = 0.0;
    fNodeC   = 0;
    fNodeN   = 0;
  }

};
//****************************************************

void GRootDCNavigator::makeMaterialMedia() {

   //--- define some materials
   fMatVacuum = new TGeoMaterial("fMatVacuum", 0,0,0);
   fMatVacuum->SetTransparency(80); 
   fMatAl = new TGeoMaterial("fMatAl", 26.98,13,2.7);

   //   //--- define some media
   fVacuum = new TGeoMedium("fVacuum",0, fMatVacuum);
   fAl = new TGeoMedium("fAl",1, fMatAl);

};
//****************************************************

void GRootDCNavigator::makeTop() {

  //--- make the top container volume
  //     +- x and y = (focal Length + 2 meters);
  //     +-z = (focal length + box z)/2.0
  //--- center of top at x=0,y=0,z = (focal length + box z)/2

  Double_t topx = fTopDim[0];
  Double_t topy = fTopDim[1];
  Double_t topz = fTopDim[2];

  // half side length of topbox is topx,topy,topz
  fTopVol = fGeom->MakeBox("fTopVol",fVacuum,topx,topy,topz);
  fGeom->SetTopVolume(fTopVol);
  //fGeom->SetTopVisible(1);
  fGeom->SetTopVisible(0);

};
//****************************************************

void GRootDCNavigator::makeFocusBox() {

  TGeoCombiTrans *combiFocBox;
  TGeoVolume *fFocBoxVol;

  // dimensions of box, to be created at center to TOP
  Double_t boxX = fFocusBoxDim[0]/2.0;
  Double_t boxY = fFocusBoxDim[1]/2.0;
  Double_t boxZ = fFocusBoxDim[2]/2.0;

  // focus box position relative to TOP center
  Double_t boxPosX = 0.0;
  Double_t boxPosY = 0.0;
  Double_t boxPosZ = (*fTopPosV)[2] - fepsil - boxZ;
  fFocBoxZTop = boxPosZ;

  //Double_t boxR1 = 45.0;
  //Double_t boxR2 = 0.0;
  //Double_t boxR3 = 0.0;
  Double_t boxR1 = fFocusBoxRot[0];
  Double_t boxR2 = fFocusBoxRot[1];
  Double_t boxR3 = fFocusBoxRot[2]; 

  if (fDebugT) {
    *oLog << "********** focus box translation & rotation *************" << endl;
    *oLog << "   size+- " 
         << boxX << " " << boxY << " " << boxZ << endl;
    *oLog << "   position wrt TopVol:  "
         << boxPosX << " " << boxPosY << " " << boxPosZ << endl;
    *oLog << "   rot.A  "
         << boxR1 << " " << boxR2 << " " << boxR3 << endl;
    *oLog << "*****************************************" << endl;
  }

  fFocBoxVol = fGeom->MakeBox("fFocBoxVol", fAl, boxX, boxY, boxZ);
  fFocBoxVol->SetLineColor(kRed);
  fFocBoxVol->SetVisibility(kTRUE);

  fFocBoxVol = fGeom->MakeBox("fFocBoxVol", fAl, boxX, boxY, boxZ);
  fFocBoxVol->SetLineColor(kRed);
  fFocBoxVol->SetVisibility(kTRUE);
  
  combiFocBox = new TGeoCombiTrans("combiFocBox",boxPosX,boxPosY,boxPosZ,
                         new TGeoRotation("rFocBox",boxR1,boxR2,boxR3)); 

  fTopVol->AddNode(fFocBoxVol, fnodeNum,combiFocBox);
  fnodeNum++;
};
//****************************************************************

void GRootDCNavigator::makeEdgeBoxes() {

  bool debug = false ;
  if (debug) {
    *oLog << "  -- GRootDCNavigator::makeEdgeBoxes" << endl;
    for (int i = 0;i<4;i++) {
      *oLog << " i gDC->edgeX/Y/Z " << gDC->edgeX[i] << "  " 
            << gDC->edgeY[i] << " " << gDC->edgeZ[i] << endl;
    }
  }

  // load edgebox size, offset, and rotation vectors
  vector<double> vsizex;
  vector<double> vsizey;
  vector<double> vsizez;
  vector<double> voffset;
  vector<double>  vdiagh;
  vector<double> vrot1;
  vector<double> vrot2;
  vector<double> vrot3;

  for (int i = 0;i<4;i++) {
    vsizex.push_back(gDC->edgeX[i]/2.0);
    vsizey.push_back(gDC->edgeY[i]/2.0);
    vsizez.push_back(gDC->edgeZ[i]/2.0);
    voffset.push_back(gDC->edgeOffset[i]);

    vrot1.push_back(gDC->edgeRot1[i]);
    vrot2.push_back(gDC->edgeRot2[i]);
    vrot3.push_back(gDC->edgeRot3[i]);

    vdiagh.push_back(( (fFocusBoxDim[0]/2.0)*sqrt(2.0) ) + voffset[i]);
   }

  // create edgebox volume vectors
  vector<TGeoVolume *> vEdgeBoxVol;
  for (int i = 0;i<4;i++) { 
    ostringstream os;
    os << "fEdgeBoxVol" << i+1;
    string edgebStr = os.str();
    vEdgeBoxVol.push_back(fGeom->MakeBox(edgebStr.c_str(),fAl,vsizex[i],vsizey[i],vsizez[i]));
    vEdgeBoxVol[i]->SetLineColor(kRed);
    vEdgeBoxVol[i]->SetVisibility(kTRUE);
  }

  // edgebox locations
  double aXbox[4] = {vdiagh[0],- vdiagh[1],0.0,0.0};
  double aYbox[4] = {0.0,0.0,vdiagh[2],- vdiagh[3]};
  double aZbox[4] = {fFocBoxZTop,fFocBoxZTop,fFocBoxZTop,fFocBoxZTop,};

  for (int i = 0;i<4;i++) {
    ostringstream os;
    ostringstream os1;
    os << "combiEdgeBox" << i+1;
    os1 << "rEdgeBox" << i+1;
    string combi = os.str();
    string redge = os1.str();
    fTopVol->AddNode(vEdgeBoxVol[1], fnodeNum,new TGeoCombiTrans(combi.c_str(),
								 aXbox[i],aYbox[i],aZbox[i],
								 new TGeoRotation(redge.c_str(),vrot1[1],vrot2[1],vrot3[1])) );
 fnodeNum++;
  }

};
//****************************************************

void GRootDCNavigator::makeShutter() {

  Double_t fShutterDim[3];
  // set dimensions of shutter
  fShutterDim[0] = gDC->shutterX;   // 9 inches
  //fShutterDim[0] = 0.2286;   // 9 inches
  fShutterDim[1] = fFocusBoxDim[0]/2.;   // width of focus box
  //fShutterDim[2] = .01;   // make it thin
  fShutterDim[2] = gDC->shutterZ;   // make it thin

  Double_t xs = fShutterDim[0];
  Double_t ys = fShutterDim[1];
  Double_t zs = fShutterDim[2];


  // Set up translation of shutter from Top center 
  Double_t xt = - (fFocusBoxDim[0] + fShutterDim[0])*sqrt(2.0)/4.0;
  Double_t yt = - xt;
  Double_t zt = (fFL - fepsil - fFocusBoxDim[2])/2;

  TGeoVolume *fShutterVol = fGeom->MakeBox("fShutterVol",fAl,xs,ys,zs);
  fShutterVol->SetLineColor(kGreen);
  //Double_t rotate = -45.0;
  Double_t rotate = gDC->shutterRot1;
  Double_t rotate2 = gDC->shutterRot2;
  Double_t rotate3 = gDC->shutterRot3;
  

  if (fDebugT) {
    *oLog << "************* shutter ************************" << endl;
    *oLog << "   shutter size: " << 2.*xs << " " << 2.*ys << " " << 2.*zs << endl;
    *oLog << "   translation:  " << xt << " " << yt << " " << zt << endl;
    *oLog << "   rotation:     " << rotate << "  " << rotate2 << "  " 
	  << rotate3 << endl;
    *oLog << "******************* end shutter ********************" << endl;
  }

  fTopVol->AddNode(fShutterVol,fnodeNum,
                   new TGeoCombiTrans("combiShutter",xt,yt,zt,
                                      new TGeoRotation("rShutter",
                                                       rotate,rotate2,rotate3)) ); 
  fnodeNum++;

};
//***************************************************

void GRootDCNavigator::makeQuadArms() {

  TGeoVolume *fQuadArmVol[4];

  Double_t fQuadArmDim[2];      // quad arm dimension
  Double_t offset;
  //fQuadArmDim[0] = .1778; // meters, 7 inches
  //fQuadArmDim[1] = .1778; // meters, 7 inches
  //fQuadArmDim[2] = 0.0;   // to be determined
  fQuadArmDim[0] = gDC->quadArmX;
  fQuadArmDim[1] = gDC->quadArmY;
  offset = gDC->quadArmOffset;

  if (fDebugT) {
    *oLog << "*************** quad arms **********************" << endl;
    *oLog << "    quad arm x and y dimensions (7 inches): "
	  << fQuadArmDim[0] << " " << fQuadArmDim[1] << "  " << offset << endl;
  }

  for (int i = 0;i<4;i++) {
    fQuadArmR2R1V[i] = new TVector3();
    fQuadArmPosV[i]  = new TVector3();
  }
  // locate top of arm, telescope coor.
  Double_t diagh = ( (fFocusBoxDim[0]/2.0)*sqrt(2.0) ) + offset;

  // ============== QuadArm0
 
  TVector3 r2(0.0,diagh,fFL);
 
  // locate bottom of arm, telescope coor
  //Double_t xArmM = 0.0; 
  //Double_t yArmM = 4.436;
  Double_t xArmM = gDC->quadArmBottomX[0]; 
  Double_t yArmM = gDC->quadArmBottomY[0];
  Double_t zArmM = fFL - sqrt(fFL*fFL - xArmM*xArmM - yArmM*yArmM);

  TVector3 r1(xArmM,yArmM,zArmM);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R = r2 - r1;
  Double_t Mag = R.Mag();
  TVector3 Unit = R.Unit();
  *fQuadArmR2R1V[0] = R;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc = r1 + r2; 
  Double_t tmp1 = (Rloc.Mag() )/2.0;
  Rloc.SetMag(tmp1);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop = Rloc - (*fTopPosV);
  Double_t xarmT = RlocTop(0);
  Double_t yarmT = RlocTop(1);
  Double_t zarmT = RlocTop(2);

  *fQuadArmPosV[0] = RlocTop;

  // set up rotation for quadarm1
  Double_t theta = R.Theta()*(TMath::RadToDeg());
  Double_t phi   = R.Phi()*(TMath::RadToDeg());

  Double_t euler0 = -(90.0 - phi);
  Double_t euler1 = -theta;

  Double_t xsize = fQuadArmDim[0]/2.;
  Double_t ysize = fQuadArmDim[1]/2.;
  Double_t zsize = Mag/2.;

  if (fDebugT) {
    *oLog << "      Quad Arm 1 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 1" << endl;
    r2.Print();
    *oLog << "   Bottom of arm 1" << endl;
    r1.Print();
    *oLog << "   length of quad arm: " << 2.*zsize << endl;
    *oLog << "   translate (TopVol coor): " << xarmT << " " << yarmT 
         << " " << zarmT << endl;
    *oLog << "   rotation: " << euler0 << " " << euler1 << "  0.0" 
         << endl<< endl;
  }


  fQuadArmVol[0] = fGeom->MakeBox("fQuadArm0Vol",fAl,xsize,ysize,zsize);
  fTopVol->AddNode(fQuadArmVol[0],fnodeNum,
                   new TGeoCombiTrans("combiQArm0",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm0",
                                                       euler0,euler1,0.0) ) );
  fnodeNum++;


  // ============== QuadArm1
  // locate top of arm, telescope coor.
 
  TVector3 r21(0.0,-diagh,fFL);
  // locate bottom of arm, telescope coor
  //Double_t xArmM1 = 0.0; 
  //Double_t yArmM1 = -4.436;
  Double_t xArmM1 = gDC->quadArmBottomX[1]; 
  Double_t yArmM1 = gDC->quadArmBottomY[1];
  Double_t zArmM1 = fFL - sqrt(fFL*fFL - xArmM1*xArmM1 - yArmM1*yArmM1);
  //Double_t zArmM1 = 0.0;
 
  TVector3 r11(xArmM1,yArmM1,zArmM1);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R1 = r21 - r11;
  Double_t Mag1 = R1.Mag();
  TVector3 Unit1 = R1.Unit();
  *fQuadArmR2R1V[1] = R1;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc1 = r11 + r21; 
  Double_t tmp11 = (Rloc1.Mag() )/2.0;
  Rloc1.SetMag(tmp11);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop1 = Rloc1 - (*fTopPosV);
  Double_t xarmT1 = RlocTop1(0);
  Double_t yarmT1 = RlocTop1(1);
  Double_t zarmT1 = RlocTop1(2);

  *fQuadArmPosV[1] = RlocTop1;

  // set up rotation for quadarm1
  Double_t theta1 = R1.Theta()*(TMath::RadToDeg());
  Double_t phi1   = R1.Phi()*(TMath::RadToDeg());

  Double_t euler01 = -(90.0 - phi1);
  Double_t euler11 = -theta1;
  Double_t xsize1 = fQuadArmDim[0]/2.;
  Double_t ysize1 = fQuadArmDim[1]/2.;
  Double_t zsize1 = Mag1/2.;

  if (fDebugT) {
    *oLog << "      Quad Arm 2 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 2" << endl;
    r21.Print();
    *oLog << "   Bottom of arm 2" << endl;
    r11.Print();
    *oLog << "   length of quad arm: " << 2.*zsize1 << endl;
    *oLog << "   translate (TopVol coor): " << xarmT1 << " " << yarmT1 
         << " " << zarmT1 << endl;
    *oLog << "   rotation: " << euler01 << " " << euler11 << "  0.0" 
         << endl<< endl;
  }

  fQuadArmVol[1] = fGeom->MakeBox("fQuadArm1Vol",fAl,xsize1,ysize1,zsize1);
  fTopVol->AddNode(fQuadArmVol[1],fnodeNum,
                   new TGeoCombiTrans("combiQArm1",xarmT1,yarmT1,zarmT1,
                                      new TGeoRotation("rArm1",
                                                       euler01,euler11,0.0)) );
  fnodeNum++;
  
  // ============== QuadArm2
  // locate top of arm, telescope coor.
 
  TVector3 r22(diagh,0.0,fFL);
  // locate bottom of arm, telescope coor
  //Double_t xArmM2 = 4.18; 
  //Double_t yArmM2 = 0.0;
  Double_t xArmM2 = gDC->quadArmBottomX[2]; 
  Double_t yArmM2 = gDC->quadArmBottomY[2];
  Double_t zArmM2 = fFL - sqrt(fFL*fFL - xArmM2*xArmM2 - yArmM2*yArmM2);
 
  TVector3 r12(xArmM2,yArmM2,zArmM2);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R2 = r22 - r12;
  Double_t Mag2 = R2.Mag();
  TVector3 Unit2 = R2.Unit();
  *fQuadArmR2R1V[2] = R2;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc2 = r12 + r22; 
  Double_t tmp12 = (Rloc2.Mag() )/2.0;
  Rloc2.SetMag(tmp12);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop2 = Rloc2 - (*fTopPosV);
  Double_t xarmT2 = RlocTop2(0);
  Double_t yarmT2 = RlocTop2(1);
  Double_t zarmT2 = RlocTop2(2);

  *fQuadArmPosV[2] = RlocTop2;

  // set up rotation for quadarm2
  Double_t theta2 = R2.Theta()*(TMath::RadToDeg());
  Double_t phi2   = R2.Phi()*(TMath::RadToDeg());
  
  Double_t euler02 = -(90.0 - phi2);
  Double_t euler12 = -theta2;

  Double_t xsize2 = fQuadArmDim[0]/2.;
  Double_t ysize2 = fQuadArmDim[1]/2.;
  Double_t zsize2 = Mag2/2.;

  if (fDebugT) {
    *oLog << "      Quad Arm 3 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 3" << endl;
    r22.Print();
    *oLog << "   Bottom of arm 3" << endl;
    r12.Print();
    *oLog << "   length of quad arm: " << 2.*zsize2 << endl;
    *oLog << "   translate (TopVol coor): " << xarmT2 << " " << yarmT2 
         << " " << zarmT2 << endl;
    *oLog << "   rotation: " << euler02 << " " << euler12 << "  0.0" 
         << endl<< endl;
  }

  fQuadArmVol[2] = fGeom->MakeBox("fQuadArm2Vol",fAl,xsize2,ysize2,zsize2);
  fTopVol->AddNode(fQuadArmVol[2],fnodeNum,
                   new TGeoCombiTrans("combiQArm2",xarmT2,yarmT2,zarmT2,
                                      new TGeoRotation("rArm2",
                                                       euler02,euler12,0.0)) );
  fnodeNum++;

  // ============== QuadArm3
  // locate top of arm, telescope coor.
 
  TVector3 r23(-diagh,0.0,fFL);
  // locate bottom of arm, telescope coor
  //Double_t xArmM3 = -4.18; 
  //Double_t yArmM3 = 0.0;
  Double_t xArmM3 = gDC->quadArmBottomX[3]; 
  Double_t yArmM3 = gDC->quadArmBottomY[3];
  Double_t zArmM3 = fFL - sqrt(fFL*fFL - xArmM3*xArmM3 - yArmM3*yArmM3);
 
  TVector3 r13(xArmM3,yArmM3,zArmM3);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R3 = r23 - r13;
  Double_t Mag3 = R3.Mag();
  TVector3 Unit3 = R3.Unit();
  *fQuadArmR2R1V[3] = R3;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc3 = r13 + r23; 
  Double_t tmp13 = (Rloc3.Mag() )/2.0;
  Rloc3.SetMag(tmp13);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop3 = Rloc3 - (*fTopPosV);
  Double_t xarmT3 = RlocTop3(0);
  Double_t yarmT3 = RlocTop3(1);
  Double_t zarmT3 = RlocTop3(2);

  *fQuadArmPosV[3] = RlocTop3;

  // set up rotation for quadarm3
  Double_t theta3 = R3.Theta()*(TMath::RadToDeg());
  Double_t phi3   = R3.Phi()*(TMath::RadToDeg());
  
  Double_t euler03 = -(90.0 - phi3);
  Double_t euler13 = -theta3;

  Double_t xsize3 = fQuadArmDim[0]/2.;
  Double_t ysize3 = fQuadArmDim[1]/2.;
  Double_t zsize3 = Mag3/2.;
 
  if (fDebugT) {
    *oLog << "      Quad Arm 4 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 4" << endl;
    r23.Print();
    *oLog << "   Bottom of arm 4" << endl;
    r13.Print();
    *oLog << "   length of quad arm: " << 2.*zsize3 << endl;
    *oLog << "   translate (TopVol coor): " << xarmT3 << " " << yarmT3 
         << " " << zarmT3 << endl;
    *oLog << "   rotation: " << euler03 << " " << euler13 << "  0.0" 
         << endl<< endl;
  }


  fQuadArmVol[3] = fGeom->MakeBox("fQuadArm3Vol",fAl,xsize3,ysize3,zsize3);
  fTopVol->AddNode(fQuadArmVol[3],fnodeNum,
                   new TGeoCombiTrans("combiQArm3",xarmT3,yarmT3,zarmT3,
                                      new TGeoRotation("rArm3",
                                                       euler03,euler13,0.0)) );
  fnodeNum++;
  
};
//********************************************************

void GRootDCNavigator::makeFacets() {
  // all vectors in telescope coordinates (origin at base)
  int ifirst = 0;  // flag for first time through facet loop
  TGeoVolume *facetVol = 0;

  TVector3 fpV(0.0,0.0,2*fFL);  // vector to 2FL point

  //open facet file
  ifstream infile("facet_loc_blrad",ios::in);
  if (! infile) {
    cerr << "    --  GRootDCNavigator::makeFacets " << endl;
    cerr << "      could not open file: " << "facet_loc_blrad" << endl;
    exit(0);
  }

  int facNum = 0;
  Double_t facRad = 0.0;
  Double_t fFX = 0.0, fFY=0.0, fFZ=0.0, fUse=0.0;
  
  while (infile >> facNum >> facRad >> fFX >> fFY >> fUse ) {

    if (fUse > 0.0) {

      // create a facet volume if this is first time thru loop
      // replicate this volume after the first time.
      if (ifirst == 0) {
        ifirst = 1;
        Double_t phi = 90.0;
        Double_t dphi = 360.0;
        Int_t nedges = 6;
        Int_t nz = 2;
        facetVol = fGeom->MakePgon("facet0",fAl,phi,dphi,
                                   nedges,nz);
        facetVol->SetLineColor(kBlue);
        Double_t zSect1 = -.1;
        Double_t zSect2 =  0.0;
        Double_t rmin = facRad - 0.001; 
        Double_t rmax = facRad;
        
        TGeoPgon *pgon = (TGeoPgon *)(facetVol->GetShape());
        pgon->DefineSection(0,zSect1,rmin,rmax);
        pgon->DefineSection(1,zSect2,rmin,rmax);
      }      
      
      fFZ = fFL - sqrt( (fFL*fFL) - (fFX*fFX) - (fFY*fFY) );

      TVector3 fLocV(fFX,fFY,fFZ);  // vector to facet center

     // vector from facet location to 2*fFL point
      TVector3 fcV = fpV - fLocV;
  
      Double_t thetaf = fcV.Theta()* (TMath::RadToDeg());
      Double_t phif   = fcV.Phi()* (TMath::RadToDeg());
      
      // Euler angle rotations 1,2,3
      Double_t r1 = -(90.0 - phif);
      Double_t r2 = - thetaf;
      Double_t r3 = - r1;
      
      // get translation components, facet loc wrt top center
      TVector3 fLocFacTopV = fLocV - (*fTopPosV);
      Double_t xtop = fLocFacTopV[0];
      Double_t ytop = fLocFacTopV[1];
      Double_t ztop = fLocFacTopV[2];

      fTopVol->AddNode(facetVol,fnodeNum,
                       new TGeoCombiTrans("combiFacet",
                                          xtop,ytop,ztop,
                                          new TGeoRotation("rFacet",
                                                           r1,r2,r3)) );
      fnodeNum++;
    }
  }
};
//************************************************

void GRootDCNavigator::makeCrossArms() {

  //fQuadArmR2R1V[3]
  //fQuadArmPosV[4]; // position of quadarms, TOP coordinates
  //fQuadArmR2R1V[4]; // r2-r1 vector for each quad arm

  /*  cross arms attach to the quad arms at a specificied distance
      below the focus box (the arms attach to the mirror at slightly
      different radial distances, so start at the focus box).  
      Assume the length of the quad arm is 12 meters (close), and that 
      the cross arms attach 0.4 *12 meters along the quad arm below the
      focus box, or 4.8 meters below.  
      Set this up in terms of the quad arm length, the position of
      the quad arm in TOP coordinates, and find the position of the cross
      arms on each quad arm, etc. using these vectors.

      fQuadArmR2R1V[3]; 
      fQuadArmPosV[4];  position of quadarms, TOP coordinates

  */
  // cross bar dimensions, length comes later
  // use two arrays since lengths may be different
  Double_t fCrossBar1Dim[3];  // cross bar dimensions
  Double_t fCrossBar2Dim[3];

  //fCrossBar1Dim[0] = fCrossBar2Dim[0] = 0.127;  // meters 
  //fCrossBar1Dim[1] = fCrossBar2Dim[1] = 0.127; 
  fCrossBar1Dim[0] = gDC->crossBarX[0]; 
  fCrossBar2Dim[0] = gDC->crossBarX[1];   // meters 
  fCrossBar1Dim[1] = gDC->crossBarY[0];
  fCrossBar2Dim[1] = gDC->crossBarY[1];

  fCrossBar1Dim[2] = fCrossBar2Dim[2] = 0.0;  // calc later

  //Double_t dcross = 4.8; // cross arms attach dcross from focus box
  Double_t dcross = gDC->crossBarDistBelowFocBox[0];
  //Double_t dcross2 =  gDC->crossBarDistBelowFocBox[1];

  // calculate distance of cross arm attachment from fQuadArmPos
  Double_t dc[4];
  for (int i = 0;i<3;i++) {
    dc[i] = 0.0;
  }

  TVector3 crossEndV[4];  // location vector of end of xarm,1,1,2,2
  
  // cross arm end locations, first two indices, arm1; then arm2
  if (fDebugT) {
    *oLog << "************ crossarm construction ********* " << endl;
    *oLog << "    cross arm distances below focus box " 
	  << dcross << endl;
  }

  for (int i = 0;i<4;i++) {
    Double_t armLength = fQuadArmR2R1V[i]->Mag();
    dc[i] = (armLength/2.) - dcross;
    TVector3 dcV = fQuadArmR2R1V[i]->Unit();
    dcV.SetMag(dc[i]);
    crossEndV[i] = (*fQuadArmPosV[i]) + dcV;

    if (fDebugT) {
        *oLog << "       1/2 length of quad arm " << i+1 << " " 
             << armLength/2 << endl;
        *oLog << "       attach above midpoint of quad: " 
             << i+1 << " " << dc[i] << endl;
      
    }
  }

  // vector parallel to crossarm from end to end
  TVector3 R = crossEndV[1] - crossEndV[0];

  // find cross arm location and set up directions
  // *********ARM 1 NEXT.***************
  TVector3 RlocTop = crossEndV[0] + crossEndV[1];
  Double_t tmp1 = (RlocTop.Mag() )/2.0;
  RlocTop.SetMag(tmp1);

  Double_t xarmT = RlocTop(0);
  Double_t yarmT = RlocTop(1);
  Double_t zarmT = RlocTop(2);
 
  // set up rotation for quadarm1
  Double_t theta = R.Theta()*(TMath::RadToDeg());
  Double_t phi   = R.Phi()*(TMath::RadToDeg());
  Double_t psi   = 0.0;

  Double_t euler0 = -(45.0 - phi);
  Double_t euler1 = -theta;

  TGeoVolume *fCrossArmVol[2];

  Double_t xsize = fCrossBar1Dim[0]/2.0;
  Double_t ysize = fCrossBar1Dim[1]/2.0;
  Double_t zsize = R.Mag()/2.;

  if (fDebugT) {
    *oLog << "      cross arm 1 dimensions: " << 2.*xsize << " " 
         << 2.*ysize << " " << 2.*zsize << endl;
    *oLog << "      cross arm 1 translation: " << xarmT << " " 
         << yarmT << " " << zarmT << endl;
    *oLog << "      cross arm 1 rotation: " << "0.0  90.0  0.0" << endl;
  }

  fCrossArmVol[0] = fGeom->MakeBox("fCrossArm0Vol",fAl,
                                   xsize,ysize,zsize);
  
  fTopVol->AddNode(fCrossArmVol[0],fnodeNum,
                   new TGeoCombiTrans("combiQArm0",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm0",
                                                       0.0,90.0,0.0)) );
  fnodeNum++;
  //****************** SECOND ARM ***********
  // vector parallel to crossarm from end to end
  R = crossEndV[3] - crossEndV[2];

  RlocTop = crossEndV[2] + crossEndV[3];
  tmp1 = (RlocTop.Mag() )/2.0;
  RlocTop.SetMag(tmp1);

  xarmT = RlocTop(0);
  yarmT = RlocTop(1);
  zarmT = RlocTop(2);
 
  // set up rotation for quadarm1
  theta = R.Theta()*(TMath::RadToDeg());
  phi   = R.Phi()*(TMath::RadToDeg());
  psi   = 0.0;

  euler0 = -(45.0 - phi);
  euler1 = -theta;

  xsize = fCrossBar2Dim[0]/2.0;
  ysize = fCrossBar2Dim[1]/2.0;
  zsize = R.Mag()/2.;
  
  if (fDebugT) {
    *oLog << "      cross arm 2 dimensions: " << 2.*xsize << " " 
         << 2.*ysize << " " << 2.*zsize << endl;
    *oLog << "      cross arm 2 translation: " << xarmT << " " 
         << yarmT << " " << zarmT << endl;
    *oLog << "      cross arm 2 rotation: " << "90.0  90.0  0.0" << endl;
    *oLog << " *********** end of cross arm construction *********" << endl;
  }

  fCrossArmVol[1] = fGeom->MakeBox("fCrossArm1Vol",fAl,xsize,ysize,zsize);
  
  fTopVol->AddNode(fCrossArmVol[1],fnodeNum,
                   new TGeoCombiTrans("combiQArm1",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm0",
                                                       90.0,90.0,0.0) ) );
  fnodeNum++;

};
//****************************************************

void GRootDCNavigator::movePositionToTopOfTopVol() {

  if (fDebugTr) {
    *oLog << "  -- GRootDCNavigator::movePositionToTopOfTopVol " << endl;
    *oLog << "        position prior to move to top ";
    *oLog << fPosC[0] << "  " << fPosC[1] << "  " << fPosC[2] << endl;
  }

  double rfx = fPosC[0];
  double rfy = fPosC[1];
  double rfz = fPosC[2];

  double Z = fepsil + fFocusBoxDim[2];

  double dl = fDirC[0];
  double dm = fDirC[1];
  double dn = fDirC[2];

  double Rx = rfx - (rfz - Z)*dl/dn;
  double Ry = rfy - (rfz - Z)*dm/dn;
  double Rz = Z;

  fPosC[0] = Rx;
  fPosC[1] = Ry;
  fPosC[2] = Rz;

  if (fDebugTr) {
    *oLog << "        TopVolPos in focal point coor.  ";
    for (int i = 0;i<3;i++) {
      *oLog << fPosC[i] << " ";
    }
    *oLog << endl;
  }
  
};
//****************************************************

const char * GRootDCNavigator::setPositionDirection( double *pos, double *dir) {

  gGeoManager = fGeom;
  // position and direction: in telescope coordinates
  // translate to TOP coordinates
  for (int i = 0;i< 3;i++) {
    fPosC[i] = pos[i];
    fDirC[i] = dir[i];
  }
  if (fDebugTr) {
    *oLog << "  -- GRootDCNavigator::setPositionDirection " << endl;
    *oLog << "       initial position, focal point coor: " << pos[0] 
          << " " << pos[1] << " " << pos[2] << endl;
    *oLog << "       initial direction: " << dir[0] 
          << " " << dir[1] << " " << dir[2] << endl;
  }
  // move position to intersection with topVol top surface
  if (dir[2] < 0) {
    fMoveToTop = true;
    //*oLog << "moving to top of topVol" << endl;
    movePositionToTopOfTopVol();
  } 

  // this is z component of position relative to the topVol Center
  fPosC[2] = fPosC[2] +(*fTopPosV)[2] - fFocusBoxDim[2] - fepsil; 

  fGeom->InitTrack(fPosC,fDirC);
  TGeoNode *cnode = fGeom->GetCurrentNode();
  fNodeC = cnode->GetName();

  if (fDebugTr) {
    *oLog << "           movePosition to top if dir[2] < 0 " << endl;
    *oLog << "       initial pos. TopVol coor " << fPosC[0] << " " 
          << fPosC[1] << " " << fPosC[2] << endl;

    //Double_t zOffset = (*fTopPosV)[2];
    Double_t zOffset = (*fTopPosV)[2] - fFocusBoxDim[2] - fepsil;

    *oLog << "       ftopPosV[2]: " << (*fTopPosV)[2] << endl;
    
    const Double_t *cpos = fGeom->GetCurrentPoint();
    *oLog << "       current position, tele. coor.: " << cpos[0] << " " 
         << cpos[1] << " " << cpos[2] - zOffset << endl;
    
    const Double_t *cdir = fGeom->GetCurrentDirection();
    *oLog << "       current direction: " << cdir[0] << " " 
         << cdir[1] << " " << cdir[2] << endl;
    *oLog << "       InitialNode: " << fNodeC << endl;
  }
  return fNodeC;
};

//*****************************************************
const char * GRootDCNavigator:: getNextNodeName() {

  gGeoManager = fGeom;

  TGeoNode *nextNode = fGeom->FindNextBoundary();
  //nextNode = fGeom->Step();
  fNodeN = nextNode->GetName();
  nextNode = fGeom->Step();

  if (fDebugTr) {
    //Double_t zOffset = (*fTopPosV)[2];
    Double_t zOffset = (*fTopPosV)[2] - fFocusBoxDim[2] - fepsil;
    //Double_t zOffset = 0.0;
    *oLog << "============ DebugT: First Boundary ===========" << endl;
    *oLog << "ftopPosV[2]: " << (*fTopPosV)[2] << endl;
    
    const Double_t *cpos = fGeom->GetCurrentPoint();
    *oLog << "current position, tele. coor.: " << cpos[0] << " " 
         << cpos[1] << " " << cpos[2] - zOffset << endl;
    
    const Double_t *cdir = fGeom->GetCurrentDirection();
    *oLog << "current direction: " << cdir[0] << " " 
         << cdir[1] << " " << cdir[2] << endl;
    *oLog << "InitialNode: " << fNodeN << endl << endl;
  }
  //exit(0);
  return fNodeN;
};

//************************************************************
const char * GRootDCNavigator::getNextNodeName(const double *pos, 
                                              const double *dir) {
  gGeoManager = fGeom;
  TGeoNode *nextNode = fGeom->FindNextBoundary();
  fNodeN = nextNode->GetName();
  fGeom->Step();

  TGeoNode *cnode = fGeom->GetCurrentNode();
  fNodeC = cnode->GetName();

  pos = fGeom->GetCurrentPoint();
  dir = fGeom->GetCurrentDirection();

  return fNodeC;
};

//****************************************************
const char * GRootDCNavigator::getPositionDirection(double *pos1, 
                                                    double *dir1) {
  gGeoManager = fGeom;

  const double *p,*d;
  
  p = fGeom->GetCurrentPoint();
  d = fGeom->GetCurrentDirection();
  for (int i = 0;i<3;i++) {
    pos1[i] = p[i];
    dir1[i] = d[i];
  }

  TGeoNode *cnode = fGeom->GetCurrentNode();
  const char * currnode = cnode->GetName();


  return currnode;
};

//****************************************************
void GRootDCNavigator::drawTelescope(const int &option) {
  gGeoManager = fGeom;
  gGeoManager->GetMasterVolume()->Draw("ogl");

  //   Draw("pad");        // use defaults in Draw if no openGL graphics
};

//****************************************************
void GRootDCNavigator::setTrackingDebug(bool setDebug) {
  fDebugTr = setDebug;

};
double GRootDCNavigator::getFocalBoxZBottomTopVolCoor() {
  double tmp = (*fTopPosV)[2] - fFocusBoxDim[2] - fepsil;
  return tmp;
};

//****************************************************
void GRootDCNavigator::printVariables() {

  *oLog << "============== printVariables ===========" << endl;
  for (int i = 0;i < 3;i++) {
    DEBUG(fPosC[i]);
  }
  *oLog << endl;
  for (int i = 0;i < 3;i++) {
    DEBUG(fDirC[i]);
  }
  *oLog << endl;

  DEBUG(fNodeC);

  for (int i = 0;i < 3;i++) {
    DEBUG(fPosN[i]);
  }
  *oLog << endl;

  for (int i = 0;i < 3;i++) {
    DEBUG(fDirN[i]);
  }
  *oLog << endl;

  DEBUG(fNodeN);
  *oLog << "============= finished printVariables ===" << endl;

};


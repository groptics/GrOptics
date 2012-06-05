/*
VERSION2.2
10MAY2012
*/
#ifndef GGEOTELESCOPE
#define GGEOTELESCOPE

#include "TObject.h"

class TGeoManager;
class TGeoManager;
class TGeoMaterial;
class TGeoMedium;
class TGeoVolume;
class TGeoCombiTrans;
class TGeoHMatrix;
class TGeoRotation;
class TVector3;
class GDCTelescope;
class GDCGeometry;

class GRootDCNavigator {

 private:

  GDCTelescope *DCTel;
  GDCGeometry *gDC;

  double fFL;           // focal length
  int fnodeNum;          // running node number, inc. on add

  // turn on to print details of telescope volumes (not facets)
  bool fDebugT;            // prints Top details if true

  // turn on via setTrackingDebug method
  bool fDebugTr;            // tracking details

  TGeoManager *fGeom;        // geomanager for telescope
  TGeoMaterial *fMatVacuum;  // vacuum material
  TGeoMaterial *fMatAl;      // aluminum material
  TGeoMedium *fVacuum;       // vacuum medium
  TGeoMedium *fAl;           // aluminum medium

  double fepsil;
  double fTopDim[3];    // TOP dimensions, side to side 
  TVector3 *fTopPosV;     // location of center of TOP, tel.coor.
  TGeoVolume *fTopVol;

  double fFocusBoxDim[3];   // focus box dimensions, side to side
  double fFocusBoxRot[3];   // focus box rotation (degrees)
  double fFocBoxZTop;
  
  TVector3 *fQuadArmPosV[4]; // position of quadarms, TOP coordinates
  TVector3 *fQuadArmR2R1V[4]; // r2-r1 vector for each quad arm
             
  bool fMakeFacetsFlag;
  // ========== photon location, direction, node variables

  double fPosC[3];  // current position
  double fDirC[3];  // current direction

  double fPosN[3];  // next position
  double fDirN[3];  // next location

  bool fMoveToTop;

  const char * fNodeC;      // current node
  const char * fNodeN;      // next node
 
  void initialize();         // initialization
  void setupGeometry();

  void makeMaterialMedia();  // make all media
  void makeTop();            // make top volume
  void makeFocusBox();       // make and place focus box
  void makeEdgeBoxes();      // make and place edge boxes
  void makeShutter();        // make and place shutter
  void makeQuadArms();       // make quad arms
  void makeFacets();
  void makeCrossArms();      // make and place cross arms

 public:

  // default constructor
  GRootDCNavigator(GDCTelescope *dcTel,int makeFacetsFlag = 0, bool debugT = false);

  // destructor
  ~GRootDCNavigator();

  // set initial position and direction for tracking
  // returns name of current node
  const char * setPositionDirection( double *pos, double *dir);

  // get node number after tracking from set position and direction
  const char * getNextNodeName();

  // get node number after tracking from position and direction set
  // by the parameters in this method
  const char * getNextNodeName(const double *pos, const double *dir);

  // get current position and direction from parameter arrays
  // return node number occupied by current position
  const char * getPositionDirection(double *pos1, double *dir1);

  void movePositionToTopOfTopVol();

  // draw the telescope volumes
  void drawTelescope();
  
  // print out position, direction, and node name in telescope coordinates
  void setTrackingDebug(bool setDebug);

  // get z location of bottom of focal box in topvol coordinates
  double getFocalBoxZBottomTopVolCoor();

  // debug print: print all variables, not totally implemented, don't use
  void printVariables();

  ClassDef(GRootDCNavigator,1);
};

#endif

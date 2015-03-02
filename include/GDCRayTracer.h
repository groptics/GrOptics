/*
VERSION3.1
2March2015
*/
/*!  GDCRayTracer concrete class for implementing DC raytracer with
  using ROOT for shadowing only, as in grisudet

     Version 0.0
     C. Duke
     Grinnell College
     May, 2011

*/

#ifndef GDCGEO1RAYTRACER
#define GDCGEO1RAYTRACER

class GDCTelescope;
class GOrderedGrid;
class GRootDCNavigator;

#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"

class GDCRayTracer : public GRayTracerBase {

  GDCTelescope *DCTel;

  ROOT::Math::XYZVector vPhotonLocT;  //!< photon Grd.Loc.Tel.coords
  ROOT::Math::XYZVector vPhotonDirT;  //!< photon Grd.DCos.Tel.coords

  ROOT::Math::XYZVector vPhotonOnTelT;  //!< photon locaton on tel.sphere,
  double fTimeOnTel;    //!< time (meters) from ground loc. to tel.Sphere

  bool bOnTelescopeFlag; 
  bool bOnFacetFlag;
  bool bReflFacetFlag;
  bool bOnCameraFlag;

  bool bDoGeoStruct;

  bool bOnQuadArmFlagIn;  //!< true if incoming photon strikes quad arm
  bool bOnQuadArmFlagOut; //!< true if outgoing photon strikes quad arm
  bool bOnCrossArmFlagIn; //!< true if incoming photon strikes cross arm
  bool bOnCrossArmFlagOut;//!< true if outgoing photon strikes cross arm
  bool bOnFocusBoxIn;     //!< true if incoming photon strikes focus box
  bool bOnEdgeBoxFlagIn;
  bool bOnEdgeBoxFlagOut;
  bool bOnShutterIn;      //!< true if incoming photon strikes shutter

  bool bOnFocusBoxOut;    //!< true if outgoing photon strikes focus box
  bool bOnTopIn;          //!< true if incoming photon strikes top(bottom)
  bool bOnTopOut;         //!< true if outgoing photon strikes top(top)

  GOrderedGrid *facetGrid;

  int iFacet;  //!< facet number of photon hit
  ROOT::Math::XYZVector vPhotonOnFacT; //!< photon location on facet
  double fTimeOnFacet;  //!< time change from tel.hit to facet hit(meters)

  double fTimeFacetToCamera;

  GRootDCNavigator *geoT; //!< root geometry ray tracer from GrISU  

  double fTopVolOrigin2FocBox;  //!< distance from topvolOrigin to focBox
  /*!<  \brief  determine the location of the photon on the tel.sphere

   */
  bool getPhotonOnTelescope();

  void initializeRayTracer();

  bool findFacet();

  bool getPhotonReflect();
  
  bool getPhotonOnCamera();

 public:

  GDCRayTracer(GDCTelescope *dcTel);
  
  ~GDCRayTracer();

  void setFacetGrid(GOrderedGrid *grid) {
    facetGrid = grid;
  };

  bool getPhotonLocCamera(ROOT::Math::XYZVector *pCam,
                          ROOT::Math::XYZVector *pDcos,
                          double *ptime);

  void printRayTracer();

  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
};


#endif

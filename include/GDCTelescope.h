/*
VERSION4.0
30May2016
*/

#ifndef GDCTELESCOPE
#define GDCTELESCOPE

// forward declarations
class TFile;
class TTree;
#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"

/*! DCStdFacet structure contains parameters for one 
         facet and a print method for printing these parameters
 */
struct DCStdFacet {
  
  int facNum;      /*!< facet number, starting from 0 */

  int type;        /*!< 1=Circular, 2=Hexagonal*/
  int sides;        /*!< number of sides in facet polygon (1 = circle) */
  double radius;    /*!< External radius*/
  double curv;      /*!< Spherical curvature radius*/
  double xm;        /*!< x location on the dish*/
  double ym;        /*!< y location on the dish*/
  double zm;        /*!< z location on the dish*/
  ROOT::Math::XYZVector vFacLoc; //*! facet location on dish from foc.pt. 
  ROOT::Math::XYZVector vFacCentrCurv; /*!< center of curvature, from foc.pt.*/
  ROOT::Math::XYZVector vFacPlLoc;    /*!< location of facet plane origin */
  ROOT::Math::XYZVector vUnitFacPlToCC; /*!< unit vector: facet plane to cur.cent.*/
  ROOT::Math::Rotation3D rotTelToFP;  /*! rotation matrix from tel. plane to facet plane*/
  
  double mis_align; /*!< Maximum misalignment of element axis*/
  
  double ftprot;    /*!< rotation angle of facet dir. in tel.plane */
  double ffprot;    /*!< rotation angle of facet in facet plane */
  /* currently, there is no distinction between
     these angles and ftprot=ffprot=0.0. Soon,
     these values can be read from the config.file
  */
  double roughness; /*!< Maximum angle between the real reflecting 
                      surface and the ideal sphere.*/
  double reflect;   /*!< Mirror element reflectivity degradation factor*/
  int rflctid;   /*!< Mirror element reflectivity curve identifier*/
 
  /*! \brief constructor: initializes parameters to zero
   */
  DCStdFacet();

  /*!  \brief copy constructor
   */
  DCStdFacet(const DCStdFacet &dcf);

  /*!  \brief determine location of facet curvature center
              with added mis_align error.  

      \param focLgt telescope focal length
   */
  void findFacetCurvatureCenter(const double &focLgt);

  /*! determine rotation matrix for coordinate system
      rotation from telescope frame to facet surface frame.

      The facet frame has its z axis perpendicular to facet surface
      and y axis in the telescope y-z plane. Thus, the z axis of the
      facet frame is toward the curvature center, the j axis of
      the facet frame is the intersection line of the two planes
      defined by the vector toward the curvature center and by
      the x axis of the telescope frame.  That is, the dot product 
      of the y axis of the facet frame with the x axis of the 
      telescope frame is 0.0.
   */
  void findFacetPlaneRotationMatrix();
  
  /*!  Determine the unit vector normal to the facet plane and location
       vector of the center of the facet plane. The facet plane covers the
       front of the facet. The intersection of the facet with its plane is the
       facet polygon. The polygon external radius is the external radius of
       the facet.
   */
  void findFacetPlane();

  /*! printDCStdFacet prints facet parameters.

      \param oStr output stream
   */
  void printDCStdFacet(ostream &oStr=cout);

 
};
/**************************end of DCStdFacet ******************************/

// forward declarations
class GGeometryBase;
class GRayTracerBase;
enum TelType;
enum GeoType;
enum RayTracerType;

/*! GDCTelescope concrete class for ACT telescopes
  inherits from GTelescope.

GDCTelescope provides the concrete class for Davies-Cotton
telescopes.

*/
class GDCTelescope : public GTelescope {

  friend class GRayTracerBase;
  friend class GDCRayTracer;
  friend class GDCTelescopeFactory;
  friend class GRootDCNavigator;

  void makePhotonHistoryBranches();
  //public:
  
  // parameters common to all concrete classes
  //ROOT::Math::XYZVector vLocation;  //!< telescope location (meters)
  double dRotationOffset;  //!< telescope rotation offset (meters)
  double dAvgTransitTime;      //!<  time offset (nanoseconds)
  double dPointOffsetX;    //!< telescope pointing offsetX (radians)
  double dPointOffsetY;    //!< telescope pointing offsetY (radians)
  double dRadius;          //!< telescope primary radius (meters)
  double dFocLgt;          //!< telescope focal length (meters)
  double dCamRad;          //!< camera radius (mm)
  double dPlateScaleFactor; //!< plate scale factor (cm/deg)
  double dPlateScaleFactor10; //!< plate scale factor (mm/deg)
  double dFocError;        //!< focusing error (meters) >0 camera too far
  
  double dCameraRadius;    //!< camera radius

  int iTelID;  //!< telescope id in the array.
  int iStdID;  //!< telescope standard number from the telescope factory

  string historyFileName; /*!< name of photon history file, 
                            if "", no history written */
  bool bPhotonHistoryFlag; //*< if false, no photon history file written

  string historyTreeName;
  int historyOption;
  TFile *hisF;
  TTree *hisT;

  TelType eTelType; //!< telescope type enum

  GeoType eGeoType; //!< geometry type enum
  bool bDoGeoStruct; //!< do/don't do structural ray tracing

  RayTracerType eRayTracerType; //!< ray tracer type enum

  GGeometryBase *geoStruct; //!< pointer to geometry 

  GRayTracerBase *rayTracer; //!< pointer to ray tracer

  // reflectivity map for wavelength and reflectivity
  //map<int,vector< pair<double,double> > > *mReflC; //!< reflectivity map<tab#,pair<lambda,reflectivity>>
  // map of reflection coeff. wavelgts and coeffs.
  map<int, vector<double> > *mVReflWaveLgts;
  map<int, vector<double> > *mVCoeffs;

  bool bGridOption; //!< if false, will loop through all facets, no grid search.
  //int nbinsx;  //!< number of grid bins in x direction for ordered hash table used in selecting a facet
  //int nbinsy;  //!< number of grid bins in y direction for ordered hash table used in selecting a facet
  
  vector<DCStdFacet> facet; //!< mirror element array  
  bool bEditAlignFlag;    // true if facet misalign parameter differs from standard misalignment.
  bool bEditReflectFlag;  // true if facet refl details differ from standard parameter.

  ROOT::Math::XYZVector vPhotonInjectLocRT;  //!< photon injection loc. Relative T coords
  ROOT::Math::XYZVector vPhotonInjectDcosT; //!< photon injection Dcos. T coords
  double fPhotWaveLgt;                      //!< photon wavelength (nm)

  ROOT::Math::XYZVector vPhotonInjectLocTRO; //!< photon loc.Tcoor. relative to rotation offset (focal point)
  ROOT::Math::XYZVector vRotationOffsetT;  //!< vector from tel.loc to rotation offset pt (focal Pt.)
  bool bPhotonOnTelescopeFlag;
  bool bPhotonOnFacetFlag;
  bool bFacetReflectFlag;
  ROOT::Math::XYZVector vPhotonCameraLoc;  //!< photon location on camera
  ROOT::Math::XYZVector vPhotonCameraDcos; //!< photon Dcos on camera

double fTransitTime;                //!< photon total transit time (from inject location)
  bool bPhotonOnCameraFlag;  //!< if true, photon impacts camera ok.

  ///////////////////////////////////////////////////////
  // parameters for storing in photon history root file

  double injectXTC,injectYTC,injectZTC;  //!< inject location telescope coor.
  double injectDcosXTC,injectDcosYTC,injectDcosZTC;  //!< inject Dcos telescope coor.

  int onTelFlag;  //<! if 1, photon strikes the telescope
  double telX,telY,telZ; //<! photon location on telescope sphere
  double fTimeToOnTelescope;

  int onFacetFlag; //<! if 1, photon strikes a facet
  int facetNum;  //!< number of facet hit (start counting at 1, not 0)
  double facetX,facetY,facetZ;  //<! photon facet hit location
  double fTimeOnTelToFacet;

  int reflectFlag; //!< if 1, photon reflection successful
  
  double cameraDcosXTC, cameraDcosYTC,cameraDcosZTC; //!< reflected photon dir.cos. Tel.coor.
  int onCameraFlag; //!< if 1, photon strikes camera (may be camera radius limit)
  double cameraX, cameraY, cameraZ; //<! camera hit location
  double fTimeFacetToCamera;
  
  double fRayTracerTime;
  double fNetTimeToCamera;  //<! ray tracer time - avgTransitTime

  //////////////////////////////////
  void fillPhotonHistory();

  void initializePhotonHistoryParms();

 public:

  /*! GDCTelescope constructor.
   */
  GDCTelescope();

  /*! GDCTelescope destructor.
   */
  ~GDCTelescope();

  /*! Inject Photon for ray tracing through telescope.

    \param photonLoc  photon ground location relative to telescope (tel.coor)
    \param photonDir  photon dirCosines in telescope coor. 
    \param photWavegt photon wavelength
  */
  void injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                    const ROOT::Math::XYZVector &photonDirT,
                    const double &photWaveLgt);

  /*! Obtain camera location following ray tracing. 

         \param photonLoc pointer to vector with photon camera location
         \param photonDcos pointer to vector with photon dir.cosines prior to camera hit
         \param photonTime pointer to photon transit time from injection point to camera location
	 \return true if photon strikes the camera, false otherwise.
  */
  bool getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                               ROOT::Math::XYZVector *photonDcos,
                               double *photonTime);
  
  /*! Print telescope details.
    
  */
  void printTelescope();

  /*! Draws the telescope in opengl drawing.

    /param option currently unused, no need to include in call.

   */
  void drawTelescope(const int &option = 0);

  /*! Sets PrintMode used for setting printing details
    
    \param oStr output stream
    \param prtMode print mode, see printTelescope method
  */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
    
  /*! setup PhotonHistory details.

    \param rootFile root filename for photonhistory tree
    \param treeName root history tree name 
    \param option not used.
   */ 
  void setPhotonHistory(const string &rootFile,const string &treeName,
                        const int &option = 0);

  /*! write photon history file.*/
  void writePhotonHistory();

  /*! Return average transit time. */
  double getAvgTransitTime() {
    return dAvgTransitTime;
  }

  /*! Return telescope radius. */
  double getTelescopeRadius() {
    return dRadius;
  }

  /*! Return focal length. */
  double getFocalLength() {
    return dFocLgt;
  }

  /*! Return platescale factor. */
  double getPlateScaleFactor() {
    return dPlateScaleFactor;
  }

  /*! Return camera radius */
  double getCameraRadius() {
    return dCamRad;
  }

  /*! Return ideal transit time */
  double getIdealTransitTime() {
    double tm = (dFocLgt/TMath::C() ) * 1.0e09;
    return tm;
  }

  /*! Set rayplot mode.

    \param eRayplot ray plot mode, note enum RayPlotType.
   */
  void setRayPlotMode(const enum RayPlotType &eRayPlot) {
    (void) eRayPlot;  // unused
  };
};


#endif

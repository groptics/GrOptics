/*
VERSION3.1
2March2015
*/
/*! \brief GDCTelescope concrete class for ACT telescopes
  inherits from GTelescope

GDCTelescope provides the concrete class for Davies-Cotton
telescopes

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

/*! \brief DCStdFacet structure contains parameters for one 
         facet and a print method for printing these parameters
 */
struct DCStdFacet {
  
  /*! \brief constructor: initializes parameters to zero
   */
  DCStdFacet();

  /*!  \brief copy constructor
   */
  DCStdFacet(const DCStdFacet &dcf);

  /*!  \brief determent location of facet curvature center
              with added mis_align error 

      \param focLgt telescope focal length
   */
  void findFacetCurvatureCenter(const double &focLgt);

  /*! \brief determine rotation matrix for coordinate system
             rotation from telescope frame to facet surface frame.

             Facet frame has z axis perpendicular to facet surface
             and y axis in the telescope y-z plane.
   */
  void findFacetPlaneRotationMatrix();
  
  /*!  \brief determine unit vector normal to the facet plane and location
              vector of the center of the facet plane. The facet plane covers the
              front of the facet. The intersection of the facet with its plane is the
              facet polygon. The polygon external radius is the external radius of
              the facet.
   */
  void findFacetPlane();

  /*! \brief printDCStdFacet prints facet parameters
      \param oStr output stream
   */
  void printDCStdFacet(ostream &oStr=cout);

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

};
/**************************end of DCStdFacet ******************************/

/*! \brief GDCTelescope contains the full description of a single
       telescope, including a ray-tracing method. Class inherits from
       GTelescope
*/

// forward declarations
class GGeometryBase;
class GRayTracerBase;
enum TelType;
enum GeoType;
enum RayTracerType;


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

  /*! \brief GDCTelescope constructor
   */
  GDCTelescope();

  /*! \brief GDCTelescope destructor 
   */
  ~GDCTelescope();



  /*! \brief injectPhoton for ray tracing through telescope (may need additional parameters later
    e.g. time?? 

    \param photonLoc  photon ground location relative to telescope (tel.coor)
    \param photonDir  photon dirCosines in telescope coor. 
  */
  void injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                    const ROOT::Math::XYZVector &photonDirT,
                    const double &photWaveLgt);

  /*! \brief getCameraPhotonLocation gets camera location following ray tracing. 
         SHOULD RETURN BOOL IN CASE PHOTON DOESN'T REACH CAMERA

         \param x  photon camera x location (CHANGE TO POINTER)
         \param y  photon camera y location (CHANGE TO POINTER)
         \param z  photon camera z location (CHANGE TO POINTER)
  */
  //void getCameraPhotonLocation(double &x,double &y,double &z);
  bool getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                               ROOT::Math::XYZVector *photonDcos,
                               double *photonTime);

/*! \brief setLogFile  used for logging photons, document later
 */
  //void setLogFile(const ofstream &logFile);
  
  /*! \brief printTelescope prints various details, document later
    
    \param oStr output stream
    \param prtMode print mode, to be documented later
  */
  void printTelescope();

  void drawTelescope(const int &option = 0);

  /*! \brief setPrintMode used for setting printing details
    
    \param oStr output stream
    \param prtMode pring mode, fix up later
  */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
    
  /*! \brief 
    
   */ 
  void setPhotonHistory(const string &rootFile,const string &treeName,
                        const int &option = 0);

  void writePhotonHistory();

  double getAvgTransitTime() {
    return dAvgTransitTime;
  }

  double getTelescopeRadius() {
    return dRadius;
  }

  double getFocalLength() {
    return dFocLgt;
  }

  double getPlateScaleFactor() {
    return dPlateScaleFactor;
  }

  double getCameraRadius() {
    return dCamRad;
  }

  double getIdealTransitTime() {
    double tm = (dFocLgt/TMath::C() ) * 1.0e09;
    return tm;
  }

  void setRayPlotMode(const enum RayPlotType &eRayPlot) {

  };
};


#endif

/*
VERSION3.1
2March2015
*/
/*! \brief GSCTelescope concrete class for ACT telescopes
  inherits from GTelescope

GSCTelescope provides the concrete class for a Schwarzschild-Couder
telescope.  Inherits from GTelescope base class. This class is based on 
Akira Okamura's NewSCT.C root script included in the GrOptics/scripts/SegSC
directory

*/
#ifndef GSEGSCTELESCOPE
#define GSEGSCTELESCOPE

#include "TGraph.h"

#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"

#include "AOpticsManager.h"
//class AOpticalComponent;
/*! \brief GSegSCTelescope contains the full description of a single
       SegSC telescope. Class inherits from GTelescope
*/

// forward declarations
enum TelType;
class mirrorSegmentDetails;
class TFile;
class TTree;
class ARay;
class TGraph;
class AGeoAsphericDisk;
class SegmentedMirror;
class SegmentedObscuration;

class GSegSCTelescope : public GTelescope {

  friend class GSegSCTelescopeFactory;
 
  AOpticsManager* fManager;
  AGeoAsphericDisk * fPrimaryV;
  AGeoAsphericDisk * fSecondaryV;
  AGeoAsphericDisk * fSecondaryObsV;

  ARay *ray;
  TFile *hisF;
  TTree *hisT;

  //!< top volume dimensions
  Double_t fTX;  // set at 30.0 for now, in SCTelescope = 15.0. make graphs first
  Double_t fTY;
  Double_t fTZ;

  int iTelID;  //!< telescope id in the array.
  int iStdID;  //!< telescope standard number from the telescope factory

  double fAvgTransitTime;
  Double_t fPlateScaleFactor;

  double fphotonInjectLoc[3];
  double fphotonInjectDir[3];
  double fphotWaveLgt;
  Double_t fphotonToTopVolTime;

  Double_t fInjectLoc[3];
  Double_t fInjectDir[3];
  Double_t fInjectTime;  // only sent to history file
  Double_t fInjectLambda;
 
  // general telescope parameters
  Double_t fF;     // Focal length
  Double_t fFMeters; // focal length in meters
  Double_t fAlpha; // \alpha
  Double_t fQ;     // q
  Double_t fTelRadius;

  // primary parameters
  Double_t fRpMax; // Primary radius max
  Double_t fRpMin; // Primary radius min
  Double_t fZp;  // primary position on z axis
  Int_t     fNp;   // Number of coefficients for the primary
  Double_t* fP;   // Polynomial coefficients (p0, p1 ...)
  vector<Double_t> fzp;

  // colors of elements in gl picture (set in initialize method)
  Int_t iPrimaryColor;
  Int_t iSecondaryColor;
  Int_t iSecondaryObscurationColor;
  Int_t iMAPMTCathodeColor;
  Int_t iMAPMTWindowColor;
  Int_t iMAPMTObscurationColor;
  Bool_t bSingleMAPMTmodule;

  // secondary parameters
  Double_t fRsMax; // Secondary radius max
  Double_t fRsMin; // Secondary radius min
  Double_t fZs;  // secondary position on z axis
  Int_t     fNs;   // Number of coefficients for the secondary
  Double_t* fS;   // Polynomial coefficients (s0, s1 ...)
  Int_t iNParS;
  vector<Double_t> fzs;

  // primary segment details
  Int_t iNumP1Mirrors;
  Int_t iNumP2Mirrors;
  vector<mirrorSegmentDetails *> vSegP1;
  vector<mirrorSegmentDetails *> vSegP2;

  // secondary segment details
  Int_t iNumS1Mirrors;
  Int_t iNumS2Mirrors;
  vector<mirrorSegmentDetails *> vSegS1;
  vector<mirrorSegmentDetails *> vSegS2;

  map<int, TGraph *> *mGRefl;
  Int_t iReflect;

  string historyFileName; /*!< name of photon history file, 
                            if "", no history written */
  string historyTreeName;
  bool bPhotonHistoryFlag; //*< if false, no photon history file written

  Double_t fLocLast[3];
  Double_t fDirLast[3];
  Double_t fTimeLast;
  int iHistoryOption;
  Int_t fStatusLast;
  Int_t fNPoints;
  Double_t fInitialInjectLoc[3];

  Double_t fRotationOffset;

  // camera and focal surface
  Double_t fKappa1;
  Double_t fKappa2;
  Double_t fRf;
  Double_t fZf;

  // for either ideal focal surface or camera
  Double_t fFocalSurfaceXOffset;     // xOffset, entered in mm
  Double_t fFocalSurfaceYOffset;     // yOffset, entered in mm
  Double_t fFocalSurfaceZOffset;     // zOffset, entered in mm 
  Double_t fFocalSurfacePhiOffset;   // Euler angle offset, degrees
  Double_t fFocalSurfaceThetaOffset; // Euler angle offset, degrees
  Double_t fFocalSurfacePsiOffset;   // Euler angle offset, degrees


  bool bCameraFlag; //*< if true, use MAPMT camera; false: ideal focal surface
  Double_t fPixelSize; 
  Double_t fMAPMTWidth; 
  Double_t fMAPMTLength; 
  Double_t fInputWindowThickness; 
  Double_t fMAPMTOffset; 
  Double_t fMAPMTGap; 
  Double_t fMAPMTRefIndex;  

  Double_t fCathodeTopRelToFocalSurface;
  Double_t fWindowBottomRelToFocalSurface;
  Double_t fMAPOscurationTopRelToFocalSurface;
  Double_t fCathodeBottomRelToOscurationTop;

  bool bRayPlotModeFlag;
  enum RayPlotType eRayPlotType;

  TelType eTelType; //!< telescope type enum (here will be SC)

  void makePhotonHistoryBranches();

  void fillPhotonHistory();

  void initializePhotonHistoryParms();

  void movePositionToTopOfTopVol();

  void initialize();

  void makePrimarySecondaryDisks();

  void addPrimaryF();

  void addPrimaryMirror(const char*name, SegmentedMirror *mirror);

  void addSecondaryJ();

  void addSecondaryMirror(const char*name, SegmentedMirror *mirror);

  void addSecondaryObscuration();

  void addSecondaryObscurationSeg(const char*name, 
                                  SegmentedObscuration *obscuration);

 public:

  /*! \brief GSCTelescope constructor
   */
  GSegSCTelescope();

  /*! \brief GSCTelescope destructor 
   */
  ~GSegSCTelescope();

  /*! \brief Build an ideal SC telescope with an FOV of 8 deg diameter. More
      realistic design and parameter reader should be implemented.
      See http://jelley.wustl.edu/actwiki/images/0/05/VVV-OSdefinitions_v2.pdf

      \param os8 If true, OS8 parameters in the PDF will be used, or else OS10
   */
  void buildTelescope(bool os8 = true);

  void addIdealFocalPlane();

  void addMAPMTFocalPlane();

  void closeGeometry();

  /*! \brief injectPhoton for ray tracing through telescope

    \param photonLoc  photon ground location relative to telescope (tel.coor)
    \param photonDir  photon dirCosines in telescope coor. 
  */
  void injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                    const ROOT::Math::XYZVector &photonDirT,
                    const double &photWaveLgt);

  /*! \brief getCameraPhotonLocation gets camera location following ray tracing. 
         RETURN FALSE IN CASE PHOTON DOESN'T REACH CAMERA

         \param x  photon camera x location (CHANGE TO POINTER)
         \param y  photon camera y location (CHANGE TO POINTER)
         \param z  photon camera z location (CHANGE TO POINTER)
	 \return true if photon reaches focal surface
  */
  bool getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                               ROOT::Math::XYZVector *photonDcos,
                               double *photonTime);

/*! \brief setLogFile  used for logging photons, document later
  not currently used,  *oLog global set as default. oPrtStrm is
  set in GTelescope constructor to *oLog.
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
    return fAvgTransitTime;
  }

  void setTelID(const int &telID) {
    iTelID = telID; };
  
  void setStdID(const int &stdID) {
    iStdID = stdID; };
 
  void setReflCoeffMap(map<int, TGraph *> *mGr);

  void testFocalPlane();

  double getTelescopeRadius() {
    return fTelRadius;  // max radius of primary
  }

  double getFocalLength() {
    return fFMeters;
  }

  double getPlateScaleFactor() {
    return fPlateScaleFactor;
  }

  double getIdealTransitTime() {
    //double tm = ( ( 2*fZs - fZp -fZf )* fFocLgt / TMath::C()) * 1.0e09;
    //return tm;
    return fAvgTransitTime;
  }

  TGraph * makeReflectivityGraph(const Int_t &irefl);

  void testPerformance();

  void CloseGeometry(); 

  AOpticsManager *getManager() const { return fManager;};

  void setRayPlotMode(const enum RayPlotType &eRayPlot) {
    bRayPlotModeFlag = true;
    eRayPlotType = eRayPlot;    
    bool debug = true;
    if (debug) {
      *oLog << "  -- GSegSCTelescope::setRayPlotMode " << endl;
      *oLog << "bRayPlotModeFlag eRayPlotType " 
            <<  bRayPlotModeFlag << "  " << eRayPlotType  << endl;
    }
    
  };
};


#endif

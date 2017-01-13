/*
VERSION4.0
30May2016
*/

#ifndef GSCTELESCOPE
#define GSCTELESCOPE

#include "TGraph.h"

// forward declarations
class TFile;
class TTree;
class ARay;
class TGraph;

#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"

#include "AOpticsManager.h"

/*! GSCTelescope models a single SC telescope; inherits from GTelescope
*/

// forward declarations
enum TelType;

/*! GSCTelescope class modeling a Schwarzschild-Couder telescope; 
  inherits from GTelescope 
*/
class GSCTelescope : public GTelescope {

  friend class GSCTelescopeFactory;
 
  AOpticsManager* manager;
  ARay *ray;
  TFile *hisF;
  TTree *hisT;

  double fTX;  //!< top volume dimensions
  double fTY;
  double fTZ;

  int iTelID;  //!< telescope id in the array.
  int iStdID;  //!< telescope standard number from the telescope factory

  double fAvgTransitTime;

  double fphotonInjectLoc[3];
  double fphotonInjectDir[3];
  double fphotWaveLgt;
  Double_t fphotonToTopVolTime;

  TGraph *gPrimRefl;
  double fPrimMaxLmda;
  double fPrimMinLmda;

  TGraph *gSeconRefl;
  double fSeconMaxLmda;
  double fSeconMinLmda;

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

  double fFocLgt;
  double fPlateScaleFactor;
  // primary
  double fDp;
  double fDpinner;
  double fZp;
  Double_t fPrimXoffset;
  Double_t fPrimYoffset;
  Double_t fPrimZoffset;
  Double_t fPrimThetaOffset;
  Double_t fPrimPhiAngle;
  Double_t fPrimRoughSigma;
  Double_t fPrimRoughMax;

  // secondary
  double fZs;
  double fDs;
  double fDsinner;
  Double_t fSecondXoffset;
  Double_t fSecondYoffset;
  Double_t fSecondZoffset;
  Double_t fSecondThetaOffset;
  Double_t fSecondPhiAngle;
  Double_t fSecondRoughSigma;
  Double_t fSecondRoughMax;

  // focal plane
  double fZf;
  double fk1;
  double fk2;
  Double_t fFocalPlXoffset;
  Double_t fFocalPlYoffset;
  Double_t fFocalPlZoffset;
  Double_t fFocalPlThetaOffset;
  Double_t fFocalPlPhiAngle;

  // camera
  double fPixelSize;
  double fPixelPitch;
  double fMAPMTWidth;
  double fMAPMTLength;
  double fInputWindowThickness;
  double fMAPMTAngularSize;
  double fMAPMTOffset;
  double fMAPMTGap;
  double fMAPMTRefIndex;
  bool bCameraFlag;

  bool bRayPlotModeFlag;
  enum RayPlotType eRayPlotType;

  int iNParP;
  vector<double> fzp;
  int iNParS;
  vector<double> fzs;


  TelType eTelType; //!< telescope type enum (here will be SC)

  void makePhotonHistoryBranches();

  void fillPhotonHistory();

  void initializePhotonHistoryParms();

  void movePositionToTopOfTopVol();

  void initialization();

 public:

  /*! GSCTelescope constructor.
   */
  GSCTelescope();

  /*! GSCTelescope destructor. 
   */
  ~GSCTelescope();

  /*! Build an ideal SC telescope.  More
      realistic design and parameter reader should be implemented.
      See http://jelley.wustl.edu/actwiki/images/0/05/VVV-OSdefinitions_v2.pdf;
   */
  void buildTelescope();

  /*! inject photon into the telescope in telescope coordinate system.

    \param photonLocT 3D vector giving photon location, telescope coordinates
    \param photonDirT 3D vector giving photon dir.cosines in telescope coordinates
    \param photWaveLgt photon wavelength in nanometers
  */
  void injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                    const ROOT::Math::XYZVector &photonDirT,
                    const double &photWaveLgt);

  /*!  get photon location on telescope camera after ray tracing.

    \param photonLoc pointer to 3D vector giving camera location of photon
    \param photonDcos pointer to 3D vector giving dir.cosines of photon prior to 
    striking camera.
    \param photonTime transit time of photon through telescope.
    \return true if photon reaches camera
  */
  bool getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                               ROOT::Math::XYZVector *photonDcos,
                               double *photonTime);
  
  /*! print telescope details to *oLog.
    
    \param oStr output stream
    \param prtMode print mode, see pilot file for documentation
  */
  void printTelescope();

  void drawTelescope(const int &option = 0);

  /*! select print mode
    
    \param oStr output stream
    \param prtMode print mode
  */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
   
  /*! 
    
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
 
  void setReflCv(const int &primID, const int &secondID,
		 map<int, TGraph *> *mGRefl);
  

  void testFocalPlane();

  double getTelescopeRadius() {
    return fDp;  // this is the radius not the diameter
  }

  double getFocalLength() {
    return fFocLgt;
  }

  double getPlateScaleFactor() {
    return fPlateScaleFactor;
  }

  double getIdealTransitTime() {
    double tm = ( ( 2*fZs - fZp -fZf )* fFocLgt / TMath::C()) * 1.0e09;
    return tm;
  }

  void setRayPlotMode(const enum RayPlotType &eRayPlot) {
    bRayPlotModeFlag = true;
    eRayPlotType = eRayPlot;    
    bool debug = true;
    if (debug) {
      *oLog << "  -- GSCTelescope::setRayPlotMode " << endl;
      *oLog << "bRayPlotModeFlag eRayPlotType " 
            <<  bRayPlotModeFlag << "  " << eRayPlotType  << endl;
    }
    
  };

};


#endif

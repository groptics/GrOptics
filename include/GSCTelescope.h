/*
VERSION3.1
2March2015
*/
/*! \brief GSCTelescope concrete class for ACT telescopes
  inherits from GTelescope

GSCTelescope provides the concrete class for a Schwarzschild-Couder
telescope.  Inherits from GTelescope base class.

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
//class AOpticalComponent;
/*! \brief GSCTelescope contains the full description of a single
       SC telescope. Class inherits from GTelescope
*/

// forward declarations
enum TelType;

// stubs for all methods
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

  Double_t fInjectLoc[3];
  Double_t fInjectDir[3];
  Double_t fInjectTime;  // only sent to history file
  Double_t fInjectLambda;
   
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

  /*! \brief GSCTelescope constructor
   */
  GSCTelescope();

  /*! \brief GSCTelescope destructor 
   */
  ~GSCTelescope();

  /*! \brief Build an ideal SC telescope with an FOV of 8 deg diameter. More
      realistic design and parameter reader should be implemented.
      See http://jelley.wustl.edu/actwiki/images/0/05/VVV-OSdefinitions_v2.pdf

      \param os8 If true, OS8 parameters in the PDF will be used, or else OS10
   */
  void buildTelescope(bool os8 = true);

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

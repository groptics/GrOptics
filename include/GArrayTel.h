/*
VERSION4.0
30May2016
*/
/*!  /brief GArrayTel class contains all telescope details 
            including a pointer to a GTelescope instance
 */

#ifndef GARRAYTEL
#define GARRAYTEL

#include "ATelescopeFactory.h"
#include "ATelescope.h"

// forward declarations
class ATelescope;
class GDCTelescope;
class TH2D;
class TH1D;
class TGraph;
class TCanvas;
#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"

class GArrayTel {

  ROOT::Math::XYZVector telLocGrdGC;  //!< tel.loc.ground coor.
  ROOT::Math::XYZVector telLocGrdTC; //!< tel.loc.telescope coor.
  double fWobbleN;  //!< wobble offset celestial N (radians)
  double fWobbleE;  //!< wobble offset celestial E (radians)
  double fLatitude;
  double fpointingOffsetX; //!< tel pointing offset X (radians)
  double fpointingOffsetY; //!< tel pointing offset Y (radians)

  TelType telType;    //!< telescope type enum
  int telID;          //!< telescope array ID
  int telStd;         //!< telescope standard number
  double fPlateScaleFactor;

  // primary parameters
  ROOT::Math::XYZVector vSCoreGC;  //!< primary coreHit ground coor.
  ROOT::Math::XYZVector vSDcosGC;  //!< primary dirCos. ground coor.
  ROOT::Math::XYZVector vSCoreTC;  //!< primary coreHit telescope coor.
  ROOT::Math::XYZVector vSDcosTC;  //!< primary dirCos. telescope coor.
  ROOT::Math::XYZVector vSCoreSC;  //!< primary coreHit primary(shower) coor.
  ROOT::Math::XYZVector vSDcosSC;  //!< primary dirCos. primary(shower) coor.
  double fAzPrim;                  //!< primary azimuth
  double fZnPrim;                  //!< primary zenith angle
  double fEnergy;                  //!< primary energy
  double fDelay;                   //!< timing delay

  ROOT::Math::Rotation3D rotGrdToTel;  //!< rotation matrix, grd to tel coor.
  ROOT::Math::Rotation3D rotGrdToShower;  //!< rot. matrix, grd to shower coor.
  double fAzTel;          //!< azimuthal angle (radians) for telescope
  double fZnTel;          //!< zenith angle (radians) for telescope
  ROOT::Math::XYZVector vTelDcosGC; //!< tel.dirCos. ground coor.  
  double fSrcRelToTelescopeX;
  double fSrcRelToTelescopeY;
  double fSrcRelToCameraX;
  double fSrcRelToCameraY;

    int iPrintMode;//**  moved from 91 to handle Wreorder -k

  bool   bFixedPointing   ;
  double fFixedPointingAz ;
  double fFixedPointingZn ;

  // photon parameters
  ROOT::Math::XYZVector vPhotonGrdLocGC; //!< photon grd.loc.ground coors.
  ROOT::Math::XYZVector vPhotonDcosGC;   //!< photon dir.cos.ground coors.
  double fAzPhotGr;
  double fZnPhotGr;
  double fPhotHgtEmiss;
  double fPhotGrdTime;
  double fPhotWaveLgt;
  int iPhotType;
  int iPhotTelHitNum;

  ROOT::Math::XYZVector vPhotonGrdLocTC;  //!< photon grd.loc.Tel.Coor
  ROOT::Math::XYZVector vPhotonRelLocTC;  /*< photon grd.loc.rel.to tel
                                          in tel.coor. */

  ROOT::Math::XYZVector vPhotonDcosTC;  //!< photon dir.cos.telescope coors.  

  bool bOnCamera;  //!< if true, photon on camera
  ROOT::Math::XYZVector vPhotonCameraLoc;
  ROOT::Math::XYZVector vPhotonCameraDcos;
  double fTotalPhotonTime;

  // have to define storage structures, need time, two deques are ok.
  deque<float> *dStoreLoc; //!< photon camera location deque storage
  deque<float> *dStorePix; //!< pixel location deque storage



  // variables for testTelescope
  vector<TH2D*> vHist;
  vector<TH1D*> vHistT;

  ATelescope *tel;  //!< pointer to associated telescope

 protected:

 public:

  /*! default constructor
   */
  GArrayTel();

  /*!  constructor
    \param telLocGrd telescope location in ground coordinates
    \param telOffsetX telescope x-offset on the sky, telescope coordinates
    \param telOffsetY telescope y-offset on the sky, telescope coordinates
    \param teltype enum designating telescope type
    \param telid telescope id number
    \param telstd standard telescope id for this telescope
    \param printMode printMode, see ARRAYTEL input record for details.
    \param ctel pointer to telescope instance produced by telescope factory
   */
  // C.Duke 2Feb2017, check if teltype and telstd are set in telescope instance earlier
  GArrayTel(const ROOT::Math::XYZVector telLocGrd,
            const double &telOffsetX,
            const double &telOffsetY,
            const TelType &teltype,
	    const int &telid,
            const int &telstd,
	    const int &printMode,
            const bool fixPtFlag,
            const double azFix,
            const double znFix,
	    ATelescope *ctel);

  /*! destructor
   */
  ~GArrayTel();

  /*!  \ brief set all primary parameters
    \param vSCorec primary ground loc. vector, ground coor. system
    \param vSDcosGdc primary dir.cos. vector, ground coor. system
    \param AzPrimc primary azimuthal angle
    \param ZnPrimc primary zenith angle
    \param energyc primary energy
    \param WobbleTNc  wobble offset, celestial North (radians)
    \param WobbleTEc wobble offset, celestial East (radians)
    \param Latitude observatory latitude
   */
  void setPrimary(const ROOT::Math::XYZVector &vSCorec,
                  const ROOT::Math::XYZVector &vSDcosGdc,
                  const double &AzPrimc,const double &ZnPrimc,
                  const double &energyc,const double &WobbleTNc,
                  const double &WobbleTEc,
                  const double &Latitude);
  
  /*!   set all photon parameters

    \param pGrd photon ground location vector, ground coordinates
    \param pDcos photon dir.cos. vector, ground coordinates
    \param pAz photon azimuthal angle (radians)
    \param pZn photon zenith angle (radians)
    \param pHgtEmiss photon height of emission (meters)
    \param pTime photon ground arrival time (nanoseconds)
    \param pWaveLgt photon wavelength (nanometers)
    \param pType particle type id that emitted the photon
    \param pTel telescope number for photon
   */ 
  void setPhoton(const ROOT::Math::XYZVector &pGrd,
                 const ROOT::Math::XYZVector &pDcos,
                 const double &pAz,const double &pZn,
                 const double &pHgtEmiss,const double &pTime,
                 const double &pWaveLgt, const int &pType,
                 const int &pTel);
  
  /*!  prints some ArrayTel telescope parameters, useful in debugging.
    Not currently used.
   */
  void printArrayTel();

  /*! sets photon history parameters 
    \param rootFile root filename for photon history tree
    \param treeName name of root tree containing history details
    \param option option parameter (unused at the moment?)
   */
  void setPhotonHistory(const string &rootFile,
                        const string &treeName,
                        const int &option = 0) {

    tel->setPhotonHistory(rootFile,treeName,option);
  };
  
  void writePhotonHistory() {

    tel->writePhotonHistory();
  };

  bool getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
					ROOT::Math::XYZVector *photonDcos,
			       double *photonTime);

  double getTimingDelay() {
    return fDelay;
  };

  int getTelescopeID() {
    return telID;
  };

  int getTelescopeStdNum() {
    return telStd;
  };

  TelType getTelescopeType() {
    return telType;
  };

  double getPlateScaleFactor() {
    return fPlateScaleFactor;
  }

  double getAvgTransitTime() {
    return tel->getAvgTransitTime();
  };

  void makeTelescopeTest(const string& testfile);

  void getTelLocTC(ROOT::Math::XYZVector *telLocTC) {
    *telLocTC = telLocGrdTC;
  };

  void getTelLocGC(double *telLocX, double *telLocY, double *telLocZ) {
    *telLocX = telLocGrdGC.X();
    *telLocY = telLocGrdGC.Y();
    *telLocZ = telLocGrdGC.Z();
  };

  void getAzZnTelescope(double *AzTel, double *ZnTel) {
    // get telescope Azimuth and Zenith angle, includes telescope offsets
    *AzTel = fAzTel;
    *ZnTel = fZnTel;
    return;
  };

  void getSrcRelativeToTelescope(double *srcTelX, double *srcTelY) {
    // get source location on sky relative to telescope, 
    *srcTelX = fSrcRelToTelescopeX;
    *srcTelY = fSrcRelToTelescopeY;
    return;
  };

  void getSrcRelativeToCamera(double *srcTelX, double *srcTelY) {
    // get source location on sky relative to telescope, 
    *srcTelX = fSrcRelToCameraX;
    *srcTelY = fSrcRelToCameraY;
    return;
  };

  void getCoreLocDCosTC(ROOT::Math::XYZVector *coreLocTC,
                        ROOT::Math::XYZVector *coreDcosTC) {
    *coreLocTC = vSCoreTC;
    *coreDcosTC = vSDcosTC;
  };
 
 void getCoreLocDCosSC(ROOT::Math::XYZVector *coreLocSC,
                        ROOT::Math::XYZVector *coreDcosSC) {
    *coreLocSC = vSCoreSC;
    *coreDcosSC = vSDcosSC;
  };
  
};

#endif

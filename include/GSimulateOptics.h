/*
VERSION3.1
2March2015
*/
/*!  /brief GSimulateOptics class directs optical simulations 
            
 */

#ifndef GSIMULATEOPTICS
#define GSIMULATEOPTICS

// forward declarations
class GReadPhotonBase;
class GArrayTel;
class GRootWriter;

#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"

class GSimulateOptics {

  GReadPhotonBase *reader; //!<
  GArrayTel *arrayTel;
  GRootWriter *rootWriter;

  //string sHeaderTree;  //!< header tree name
  string sFileHeader;  //!< input file header string
  string sVersion;
  double fObsHgt;      //!< observatory height from input record
  double fGlobalEffic;   //!< global efficiency from input record

  int iNShowers; //!< number of showers, < 0 no limit
  int iNPhotons; //!< number of photons/shower, < 0 no limit

  map<int, GArrayTel *> *mArrayTel;
  map<int, GArrayTel *>::iterator iterArrayTel;

  map<int, GRootWriter *> *mRootWriter;
  map<int, GRootWriter *>::iterator iterRootWriter;
  string sHeaderTree;  //!< header tree name

  TTree *allTel;

  // primary details
  ROOT::Math::XYZVector vSCore; //!< core loc.vec. ground coors.(meters)
  ROOT::Math::XYZVector vSDcosGd; //!< core dir.cosines; grd coors.
  double fZnPrim;                   //!< primary zenith angle
  double fAzPrim;                   //!< primary azimuthal angle

  // variables for the root tree, copied from GRootWriter 
  unsigned int fEventNumber;
  unsigned int fPrimaryType;
  double        fPrimaryEnergy;
  //CD:2Mar2015  new parameters
  double fFirstIntHgt;
  double fFirstIntDpt;
  unsigned int iShowerID;

  //double        fXcore;
  //double         fYcore;
  //double         fXcos;
  //double         fYcos;
  double         fXsource;
  double        fYsource;
  double         fDelay;
  double        fPhotonToCameraTime;
  
  double fWobbleE;         //!< wobbleE in celestial coor. (radians)
  double fWobbleN;         //!< wobbleN in cel.coor.(toward Polaris)(radians)    
  double fWobbleR;   //!<  radius of wobble offset (radians)
  double fLatitude;  //!<  latitude of observatory (radians)

  double fWobbleTE; //!< wobble offset E, determined from fWobbleR
  double fWobbleTN; //!< wobble offset N, determined from fWobbleR

  double fFieldRot;  /*!< angle between arc from prim to zenith and arc from
                       primary and polaris */
  double fZnTel;  //!< zenith angle of telescope (radians)
  double fAzTel;  //!< azimuthal angle of telescop (radians)
  
  ROOT::Math::XYZVector *vTelDcosGrd;
  ROOT::Math::Rotation3D *rotGrdToTel; //!< rotation matrix: GrdCoor to TelCoor.
  double fEnergy;  //!< 

  bool primaryFlag;
  // photon details
  ROOT::Math::XYZVector vPhotonGrdLoc; // MAY WANT TO MAKE ALL VECTORS POINTERS
  ROOT::Math::XYZVector vPhotonDCosGd; // AND INSTANTIATE IN CONSTRUCTOR
  ROOT::Math::XYZVector vPhotonDCosTel; // AND INSTANTIATE IN CONSTRUCTOR

  ROOT::Math::XYZVector vTmp;  // just temporary for now


  double fAzPhot;
  double fZnPhot;
  double fPhotHgtEmiss;
  double fPhotGrdTime;
  double fPhotWaveLgt;
  int    iPhotType;
  int    iPhotTelHitNum;
  bool photonFlag;
  
  /*!
   */
  void fillAllTelTree();

  //static bool sortPair(const pair<int,unsigned> &i , const pair<int,unsigned> &j) {
  //bool test = (j.second < i.second);
  //return test;
  //};

 
 public:
  /*!
   */
  GSimulateOptics();

  /*!
   */
  GSimulateOptics(GReadPhotonBase *read,
		  map<int, GArrayTel *> *mATel,
		  map<int,GRootWriter *> *mRWriter,
		  const string &headerTree);

  /*!
   */
  ~GSimulateOptics();

  /*!
   */
  void getHeader() {
    sFileHeader = reader->getHeader();
  };

  /*!
   */
  void setWobble(const double &wobbleE, const double &wobbleN,
                 const double &wobbleR, const double &latitude) {
    fWobbleE = wobbleE;
    fWobbleN = wobbleN;
    fWobbleR = wobbleR;
    fLatitude = latitude;
  };

  /*!  returns true if complete simulations successfully
   */
  bool startSimulations(const int &numShowers,
                        const int &numPhotons);

  /*! \ brief make wobble offsets based on wobble radius
   */
  void makeWobbleOffset();

  /*!
   */
  void makeRotMatrixGrdToSky();

  /*!
   */
  void printDebugTelescope();

  /*!
   */
  void printDebugPhoton();

  /*!
   */
  void printDebugPrimary();
  
};

#endif

/*
VERSION3.1
2March2015
*/
/*! \brief  GReadPhotonGrISU class: concrete class for reading 
      cherenkov photons from GrISU formatted asci files

 */

#ifndef GREADPHOTONGRISU
#define GREADPHOTONGRISU

// forward declaratiosn
#include "Math/Vector3Dfwd.h"

class GReadPhotonGrISU: public GReadPhotonBase {

  istream *pInStream; //!< input file stream pointer (can be file or cin)
  string sInFileStr;  //!< name of input file

  string sInFileHeader;  //!< header string from input file

  double fObsHgt;      //!< observatory height (meters) from H record  
  double fGlobalEffic; //!< global efficiency from R record
  unsigned int iParticleType; //!< primary particle type from header record

  string sFileLine;  //!< current record, filled by getLine() method

  GrISURecType eRecType;  //<! enum for record most recently read

  // primary (the S is for shower) parameters (on S line)
  double fSEnergy;      //!< primary energy (TeV)
  double fSAz;          //!< primary azimuth (radians)
  double fSZn;          //!< primary zenith angle (radians)
  int iSSeed[3];        //!< up to three seeds from primary record
  // new parameters CD:2Mar2015
  double fFirstIntHgt;
  double fFirstIntDpt;
  unsigned int iShowerID;

  ROOT::Math::XYZVector vSCore; //!< core loc.vec. ground coors.(meters)
  ROOT::Math::XYZVector vSDcos; //!< core dir.cosines; grd coors.

  // photon on ground
  double fPHgtEmiss; //!< photon emission height (meters)
  double fPTime;     //!< photon grd. relative arrival time (nsec)
  double fPWaveLgt;  //!< photon wavelength (nmeters)
  int    iPType;     //!< particle type that produced photon
  int    iPTel;      //!< number of telescope sphere intersected by photon
  double fPAz;       //!< photon azimuth (radians)
  double fPZn;       //!< photon zenith angle (radians)
  ROOT::Math::XYZVector vPGrd;  //!< photon grd.loc relative to telescope
  ROOT::Math::XYZVector vPDcos; //!< photon dir.cosines grd.coors. system

  double fWobbleE;
  double fWobbleN;
  double fWobbleR;
  double fLatitude;

  /*!  \brief getLine.  Get the next input record. This record becomes
                     the current record.

      Place the string in the sFileLine variable. Set the eRecType
      enum variable to identify the record type for the current record.
      Only options are S line, P line, or EOF types.

      \return true new record read ok
      \return false no record available, set EOF record type

   */  
  bool getLine();

 protected:

 public:
  
  //! constructor
  GReadPhotonGrISU();

  //!  destructor
  ~GReadPhotonGrISU();

  /*! setInputfile 
         open input file and read preliminary records 
         up to first shower records.
      \param infile inputFile name, if infile = "",
              take input from cin.
      \return true file successfully opened
      \return false file can't be opened
   */
  bool setInputFile(const string &infile);

  //! get header string 
  string getHeader() {return sInFileHeader;};

  /*!  getPrimary
       get primary details from current record 
       then call getLine() to read next record
       \param pCore vector of grd.coor core location (meters)
       \param pDCos vector of primary direction cosines
       \param Az  primary azimuth (radians)
       \param Zn  primary zenith angle (radians)
       \param energy primary energy (TeV)
       \return true current record was a primary record
       \return false current record was not a primary record;
                     could be a photon record or EOF
   */
  bool getPrimary(ROOT::Math::XYZVector *pCore, 
                  ROOT::Math::XYZVector *pDCos,double *Az,
                  double *Zn, double *energy, unsigned int *particleType,
                  double *firstIntHgt, double *firstIntDpt,
                  unsigned int *showerid);


  /*!  getPhoton
       get cherenkov photon details from current record
       then call getLine() to read next record
       \param pGrd vector of p loc. with respect to telescope (meters)
       \param pDcos vector of photon direction cosines
       \param pAz photon azimuth (radians)
       \param pZn photon zenith angle (radians)
       \param pHgtEmiss photon emission height (meters)
       \param pTime photon relative time to ground (nsec)
       \param pWaveLgt photon wavelength (nm)
       \param pType particle type emitting photon
       \param pTel  telescope number of intersecting telescope  
       \return true current record was a photon record
       \return false current line was not a photon record; 
                     could be a primary record or EOF
   */
  bool getPhoton(ROOT::Math::XYZVector *pGrd,
                 ROOT::Math::XYZVector *pDcos,
                 double *pAz,double *pZn,
                 double *pHgtEmiss,double *pTime,
                 double *pWaveLgt, int *pType,
                 int *pTel);

  /*! getObsHeight: get observatory height
      \return obsHeight obervatory height about sea level (meters)
   */
  double getObsHeight() { return fObsHgt; };

  /*! getGlobalEffic: get global efficiency
      \return globalEffic global efficiency
   */
  double getGlobalEffic(){ return fGlobalEffic; };

};

#endif

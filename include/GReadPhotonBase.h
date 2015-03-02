/*
VERSION3.1
2March2015
*/
/*!  \brief GReadPhotonBase:  base class for reading Cherenkov photons
            
  Currently one concrete class for reading GrISU format files.
  Can read from stdin if no file is specified.
 */

#ifndef GREADPHOTONBASE
#define GREADPHOTONBASE

// forward declarations
#include "Math/Vector3Dfwd.h"

class GReadPhotonBase {

 protected:

 public:
  //! constructor
  GReadPhotonBase();

  //! virtual destructor
  virtual ~GReadPhotonBase();
  
  /*! setInputfile 
         open input file and read preliminary records 
         up to first shower records.
      \param infile inputFile name
      \return true file successfully opened
      \return false file can't be opened
   */
  virtual bool setInputFile(const string &infile) = 0;

  /*! get header string
    \return header as a string
    \return "" if no header found
						
   */
  virtual string getHeader() = 0;

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
  virtual bool getPrimary(ROOT::Math::XYZVector *pCore, 
                          ROOT::Math::XYZVector *pDCos,double *Az,
                          double *Zn, double *energy, 
                          unsigned int *particleType,
                          double *firstIntHgt, double *firstIntDpt,
                          unsigned int *showerid) = 0;

  /*!  getPhoton
       get cherenkov photon details from current record
       then call getLine() to read next record
       \param pGrd vector of photon loc. in ground coor. (meters)
       \param pDcos vector of photon direction cosines in ground coor.
       \param pAz photon azimuth (radians)
       \param pZn photon zenith angle (radians)
       \param pHgtEmiss photon emission height (meters)
       \param pTime photon relative time to ground (nsec)
       \param pWaveLgt photon wavelength (nm)
       \param pType particle type emitting photon
       \param pTel  telescope number of intersecting telescope  
       \return true current record was a photon record
       \return false current line was not a primary record; 
                     could be a primary record or EOF
   */
  virtual bool getPhoton(ROOT::Math::XYZVector *pGrd,
                 ROOT::Math::XYZVector *pDcos,
                 double *pAz,double *pZn,
                 double *pHgtEmiss,double *pTime,
                 double *pWaveLgt, int *pType,
                 int *pTel) = 0;

  /*! getObsHeight: get observatory height
      \return obsHeight obervatory height about sea level (meters)
   */
  virtual double getObsHeight() = 0;

  /*! getGlobalEffic: get global efficiency
      \return globalEffic global efficiency
   */
  virtual double getGlobalEffic() = 0;


};

#endif

/*
VERSION4.0
30May2016
*/
/*
ATelescope provides the base class for DC and SC telescope
concrete classes

     C. Duke
     Grinnell College
     August 30, 2010

*/
#ifndef ATELESCOPE
#define ATELESCOPE

//#include "TMath.h"
#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"


//class GGeometryBase;
//class GRayTracerBase;

/*! ATelescope base class for ACT telescopes
 */
class ATelescope/*: public TObject*/ {
 private:

 protected:

  ostream *oPrtStrm; //*< used by all telescope types for various printing modes
  int iPrtMode;

 public:
  int iRayPlotInt;
  
  //*! constructor
  ATelescope();

  //*! virtual destructor
  virtual ~ATelescope() {  };


  /*! inject photon into the telescope in telescope coordinate system.

    \param photonLocT 3D vector giving photon location, telescope coordinates
    \param photonDirT 3D vector giving photon dir.cosines in telescope coordinates
    \param photWaveLgt photon wavelength in nanometers

   */
  virtual void injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                            const ROOT::Math::XYZVector &photonDirT,
                            const double &photWaveLgt) = 0;

  /*! get photon location on telescope camera after ray tracing.

    \param photonLoc pointer to 3D vector giving camera location of photon
    \param photonDcos pointer to 3D vector giving dir.cosines of photon prior to 
    striking camera.
    \param photonTime transit time of photon through telescope.

   */
  virtual bool getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                                       ROOT::Math::XYZVector *photonDcos,
                                       double *photonTime) = 0;

  //virtual void setLogFile(const ofstream &logFile) = 0;

  virtual void printTelescope() = 0;

  virtual void drawTelescope(const int &option = 0) = 0;

  virtual void setPrintMode(ostream &oStr=cout,const int prtMode=0) = 0;

  virtual void setPhotonHistory(const string &rootFile,const string &treeName,
                                const int &option = 0) = 0;

  virtual void writePhotonHistory() = 0;

  virtual double getAvgTransitTime() = 0;
  
  virtual double getTelescopeRadius() = 0;

  virtual double getFocalLength() = 0;

  virtual double getIdealTransitTime() = 0;

  virtual double getPlateScaleFactor() = 0;

  //C.Duke  added two new methods so don't need to define enum in ATelescope.
  // correctly implemented for GSegSCTelescope, not for DC yet
  //virtual void setRayPlotMode(const enum RayPlotType &eRayPlot) = 0;
  virtual void setRayPlotModeInt(const int &iplt);

  virtual int getRayPlotModeInt();
    //ClassDef(ATelescope, 1)
};
#endif // ATELESCOPE
//Testing git

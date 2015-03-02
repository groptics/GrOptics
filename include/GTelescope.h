/*
VERSION3.1
2March2015
*/
/*! \brief GTelescope base class for ACT telescopes

GTelescope provides the base class for DC and SC telescope
concrete classes

     Version 0.0
     C. Duke
     Grinnell College
     August 30, 2010

*/
#ifndef GTELESCOPE
#define GTELESCOPE

//#include "TMath.h"
#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"


class GGeometryBase;
class GRayTracerBase;

class GTelescope {
 private:
  float _location;

 protected:

  // used by all telescope types for various printing modes
  ostream *oPrtStrm;
  int iPrtMode;

 public:

  GTelescope();

  virtual ~GTelescope() {  };


  virtual void injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                            const ROOT::Math::XYZVector &photonDirT,
                            const double &photWaveLgt) = 0;

  //virtual void getCameraPhotonLocation(double &x,double &y,double &z) = 0;
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

  virtual void setRayPlotMode(const enum RayPlotType &eRayPlot) = 0;
};

#endif

/*
VERSION3.1
2March2015
*/
/*!  GRayTracerBase class for raytracing though any telescope.
  for use with both ROOT geometry ray tracing and with no ROOT
  tracing.  Probably will be a pure virual class


     Version 0.0
     C. Duke
     Grinnell College
     September 13, 2010

*/

#ifndef GRAYTRACERBASE
#define GRAYTRACERBASE

class GRayTracerBase {

 protected:
  // used by all telescope types for various printing modes
  ostream *oPrtStrm;
  int iPrtMode;

 public:
  
  double itest;
  
  GRayTracerBase() {};
  
  virtual ~GRayTracerBase();

  virtual bool getPhotonLocCamera(ROOT::Math::XYZVector *pCam,
                                  ROOT::Math::XYZVector *pDcos,
                                  double *ptime) = 0;
 
  virtual void printRayTracer() = 0;

  virtual void setPrintMode(ostream &oStr=cout,const int prtMode=0) = 0;
};

#endif

/*
VERSION3.1
2March2015
*/
#ifndef GREADDCSTDBASE
#define GREADDCSTDBASE 


class GDCTelescopeFactory;

class GReadDCStdBase {

 protected:
  // used by all telescope types for various printing modes
  ostream *oPrtStrm;
  int iPrtMode;

  GDCTelescopeFactory *DCFac;

  virtual void setupDCFactory() = 0;

 public:

  GReadDCStdBase();

  virtual ~GReadDCStdBase() {};

  virtual void setDCTelescopeFactory(GDCTelescopeFactory *DCFactory) = 0;

  virtual void setPrintMode(ostream &oStr=cout,const int prtMode=0) = 0;
};

#endif



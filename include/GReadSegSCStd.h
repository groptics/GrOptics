/*
VERSION4.0
30May2016
*/
/*! \brief  GReadSegSCStd class: concrete class for reading 
      standard telescope configurations for use by telescope factory

 */

#ifndef GREADSEGSCSTD
#define GREADSEGSCSTD

class GPilot; 
class GSegSCTelescopeFactory;

class GReadSegSCStd {

  // set factory as a friend or the reverse?
 private: 

  // used by all telescope types for various printing modes
  //ostream *oPrtStrm;
  //int iPrtMode;

  GSegSCTelescopeFactory *SegSCFac;

  // for reading pilot file, can be appended files
  string spilotfile;   //!< pilot file
  GPilot *pi;        //!< pilot file reader
  vector<string> tokens; //!< vector of pilot record string tokens
  string flagline;       //!< record flag for pilot file reading

  //map<int, vector<double> > *mVReflWaveLgts;
  //map<int, vector<double> > *mVCoeffs;

  int iStdNum;  //!< active telescope number used in filling entries

  SegSCStdOptics *opt;  //!< working SCStdOptics from DCFac.

  void setupSegSCFactory();

  void getReflCoeff();

    void getTranCurve();

  void getPolyCoeffs();

  void readBasicRecord(const vector<string> &tokens,
                       const MirSeg &eMirPS,
                       SegSCStdOptics *opt);
 
void readChangeRecord(const vector<string> &tokens,
                       const MirSeg &eMirPS,
                       SegSCStdOptics *opt);
 public:

  GReadSegSCStd(const string &pilotfile,GSegSCTelescopeFactory *SegSCFactory );

  GReadSegSCStd(const string &pilotfile);

  ~GReadSegSCStd();

  void setSegSCTelescopeFactory(GSegSCTelescopeFactory *SegSCFactory);  

  //void setPrintMode(ostream &oStr=cout,const int prtMode=0);

};

#endif

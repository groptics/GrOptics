/*
VERSION3.1
2March2015
*/
/*! \brief  GReadSCStd class: concrete class for reading 
      standard telescope configurations for use by telescope factory

 */

#ifndef GREADSCSTD
#define GREADSCSTD

class GPilot;

class GReadSCStd {

  // set factory as a friend or the reverse?
 private: 

  // used by all telescope types for various printing modes
  //ostream *oPrtStrm;
  //int iPrtMode;

  GSCTelescopeFactory *SCFac;

  // for reading pilot file, can be appended files
  string spilotfile;   //!< pilot file
  GPilot *pi;        //!< pilot file reader
  vector<string> tokens; //!< vector of pilot record string tokens
  string flagline;       //!< record flag for pilot file reading

  //map<int, vector<double> > *mVReflWaveLgts;
  //map<int, vector<double> > *mVCoeffs;

  int iStdNum;  //!< active telescope number used in filling entries

  SCStdOptics *opt;  //!< working SCStdOptics from DCFac.

  void setupSCFactory();

  void getReflCoeff();

  void getPolyCoeffs();

 public:

  GReadSCStd(const string &pilotfile,GSCTelescopeFactory *SCFactory );

  GReadSCStd(const string &pilotfile);

  ~GReadSCStd();

  void setSCTelescopeFactory(GSCTelescopeFactory *SCFactory);  

  //void setPrintMode(ostream &oStr=cout,const int prtMode=0);

};

#endif

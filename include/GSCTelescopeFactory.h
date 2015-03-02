/*
VERSION3.1
2March2015
*/
/*! \brief GSCTelescopeFactory concrete class for creating ACT 
  Telescopes

Modified Factory design pattern where GTelescopeFactory provides 
the base class for DC and SC telescope concrete factories. 
GSCTelescopeFactory produces the SC telescopes.

*/
#ifndef GSCTELESCOPEFACTORY
#define GSCTELESCOPEFACTORY

// forward declarations
class TGraph;
class GPilot; 
class GSCTelescope;
class GTelescope;
class GReadSCStd;

/*!  /brief SCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
struct SCStdOptics {

  int iPrimReflID;
  int iSecReflID;
  ostream *oStr;
  int iPrtMode;

  TelType stdType; 
  int stdNum;
  double fFocLgt;
  double fAvgTransitTime;
  Double_t fRotationOffset;
  Double_t fPlateScaleFactor;

  // primary
  double fDp;
  double fDpinner;
  double fZp;
  Double_t fPrimXoffset;
  Double_t fPrimYoffset;
  Double_t fPrimZoffset;
  Double_t fPrimThetaOffset;
  Double_t fPrimPhiAngle;
  Double_t fPrimRoughSigma;
  Double_t fPrimRoughMax;

  // secondary
  double fZs;
  double fDs;
  double fDsinner;
  Double_t fSecondXoffset;
  Double_t fSecondYoffset;
  Double_t fSecondZoffset;
  Double_t fSecondThetaOffset;
  Double_t fSecondPhiAngle;
  Double_t fSecondRoughSigma;
  Double_t fSecondRoughMax;

  // focal plane
  double fZf;
  double fk1;
  double fk2;
  Double_t fFocalPlXoffset;
  Double_t fFocalPlYoffset;
  Double_t fFocalPlZoffset;
  Double_t fFocalPlThetaOffset;
  Double_t fFocalPlPhiAngle;

  // camera
  double fPixelSize;
  double fPixelPitch;
  double fMAPMTWidth;
  double fMAPMTLength;
  double fInputWindowThickness;
  double fMAPMTAngularSize;
  double fMAPMTOffset;
  double fMAPMTGap;
  double fMAPMTRefIndex;
  bool bCameraFlag;

  int iNParP;
  vector<double> fzp;
  int iNParS;
  vector<double> fzs;

  SCStdOptics(); 

  ~SCStdOptics();

  // copy constructor
  SCStdOptics(const SCStdOptics &sco);

  void setPrintOptions(ostream *outStr, const int &prtMode) {
    oStr = outStr;
    iPrtMode = prtMode;
  };

  void printSCStdOptics();

};

////////////////////////////////////////////////////////////////

/*! \brief  GSCTelescopeFactory concrete class for constructing 
     SC telescope, inherits from GTelescopeFactory to implement
     factory method
 */
class GSCTelescopeFactory : public GTelescopeFactory {
 private:

  friend class GReadSCStd; 

  GReadSCStd *readSC;  //!< SC base reader

  GPilot *pi;  //!< pilot reader pointer
  vector<string> tokens;  //!< string vector for GPilot use
  string sPilotEdit;  //!< pilot file string from edit record

  map<int,SCStdOptics*> mStdOptics;    //*< map of standard telescopes
  map<int,SCStdOptics*>::iterator itmStdOp;  //*< iterator for this map

  map<int, TGraph *> *mGRefl;
  map<int, TGraph *>::iterator itmGRefl;
  
  SCStdOptics *opt;  //*< working stdOptics for current telescope
  
  GSCTelescope *SCTel;  //*< pointer to working telescope
  int iNumSCTelMade;
  /*! \brief editWorkingTelescope makes edits based on 
           pilotfile entries to telescope currently 
           under construction

           \param SCTel pointer to current telescope
   */
  void editWorkingTelescope(GSCTelescope *SCTel1);
 
 public:

  /*!  GSCTelescopeFactory constructs a SC telescope from 
       standard telescopes obtained from reader. USE THIS CONSTRUCTOR

       \param dcReader SC telescope reader instance 
       \param editPilotFile pilotfile name containing telescope edit records
   */
  GSCTelescopeFactory(GReadSCStd &scReader,
		    const string &editPilotFile);
  /*! \brief makeTelescope constructs a SC telescope based on 
             instructions in the pilot file

             \param id telescope id within array
             \param std number of standard telescope
             \param xLoc x position within array (meters)
             \param yLoc y position within array (meters)
             \param zLoc z position within array (meters)
             \return GSCTelescope pointer to constructed telescope
  */

  ~GSCTelescopeFactory();

  GSCTelescope *makeTelescope(const int &id,
                               const int &std);
  
  /*! \brief printStdTelescope debug print based on prtMode

      \param iStd standard telescope id
      \param mode printmode (see SCStdOptics::printSCStdOptics
      in this file, default 0
      \param oStr output stream, default cout
   */
  void printStdTelescope(const int &iStd, const int &mode = 0,
                         ostream &oStr=cout);

  /*! \brief setPrintMode defines the print mode for 
          SCStdOptics::printSCStdOptics (documented in this file)

      \param oStr output stream, default cout
      \param prtMode  see SCStdOptics::printSCStdOptics earlier in this file
      \param prtMode  -1: print all standard telescopes
   */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
};

#endif

/*
VERSION3.1
2March2015
*/
/*! \brief GSegSCTelescopeFactory concrete class for creating ACT 
  Telescopes

Modified Factory design pattern where GTelescopeFactory provides 
the base class for DC and SC telescope concrete factories. 
GSegSCTelescopeFactory produces the SC telescopes.

*/
#ifndef GSEGSCTELESCOPEFACTORY
#define GSEGSCTELESCOPEFACTORY

// forward declarations
class TGraph;
class GPilot; 
class GSegSCTelescope;
class GTelescope;
class GReadSegSCStd;
class mirrorSegmentDetails;

// move following declaration to GDefinition.h
// if a structure/variable used in more than one file
// put the declaration in GDefinition.h
/*
  struct mirrorSegmentDetails {
  Double_t rmin;
  Double_t rmax;
  Double_t margin;
  Double_t delPhi;
  Int_t reflect;
  Double_t posErrorX;
  Double_t posErrorY;
  Double_t posErrorZ;
  Double_t rotErrorPhi;
  Double_t rotErrorTheta;
  Double_t rotErrorPsi;
  Double_t roughness;
  Int_t bRead; // if 0, set from BASIC; if 1, set from CHANGE
};
*/

/*!  /brief SegSCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
struct SegSCStdOptics {

  Int_t iPrimReflID;
  Int_t iSecReflID;
  ostream *oStr;
  Int_t iPrtMode;

  TelType stdType; 
  Int_t stdNum;
  Double_t fAvgTransitTime;
  Double_t fRotationOffset;
  Double_t fPlateScaleFactor;

  Double_t fF;     // Focal length
  Double_t fAlpha; // \alpha
  Double_t fQ;     // q
  Double_t fRpMax; // Primary radius max
  Double_t fRpMin; // Primary radius min
  Double_t fRsMax; // Secondary radius max
  Double_t fRsMin; // Secondary radius min

  Double_t fKappa1;// Focal plane sag constant
  Double_t fKappa2;// Focal plane sag constant
  Double_t fRf;  // focal surface radius

  Double_t fZp;  // primary position on z axis
  Double_t fZs;  // secondary position on z axis, foc.length units
  Double_t fZf;  // focal surf. position on z axis, foc.length units

  Int_t     fNp;   // Number of coefficients for the primary
  Int_t     fNs;   // Number of coefficients for the secondary
  Double_t* fP;   // Polynomial coefficients (p0, p1 ...)
  Double_t* fS;   // Polynomial coefficients (s0, s1 ...)

  // MAPMT Parameters
  bool bCameraFlag;
  Double_t fPixelSize;            // Width of a MAPMT pixel
  Double_t fMAPMTWidth;           // Housing width
  Double_t fMAPMTLength;          // between input window and anode pins
  Double_t fInputWindowThickness; // Thickness of the input window
  Double_t fMAPMTGap;
  Double_t fMAPMTRefIndex;
  Double_t fMAPMTOffset;
  Bool_t bSingleMAPMTmodule;

  Double_t fFocalSurfaceXOffset;     // xOffset, entered in mm
  Double_t fFocalSurfaceYOffset;     // yOffset, entered in mm
  Double_t fFocalSurfaceZOffset;     // zOffset, entered in mm 
  Double_t fFocalSurfacePhiOffset;   // Euler angle offset, degrees
  Double_t fFocalSurfaceThetaOffset; // Euler angle offset, degrees
  Double_t fFocalSurfacePsiOffset;   // Euler angle offset, degrees

  Int_t iNParP;
  vector<Double_t> fzp;
  Int_t iNParS;
  vector<Double_t> fzs;

  Int_t iNumP1Mirrors;
  Int_t iNumP2Mirrors;
  Int_t iNumS1Mirrors;
  Int_t iNumS2Mirrors;

  vector<mirrorSegmentDetails *> vSegP1;
  vector<mirrorSegmentDetails *> vSegP2;
  vector<mirrorSegmentDetails *> vSegS1;
  vector<mirrorSegmentDetails *> vSegS2;

  SegSCStdOptics(); 

  ~SegSCStdOptics();

  // copy constructor
  SegSCStdOptics(const SegSCStdOptics &sco);

  void setPrintOptions(ostream *outStr, const int &prtMode) {
    //oStr = outStr;
    //iPrtMode = prtMode;
  };

  void printSegSCStdOptics();
  void printSegVector (const vector<mirrorSegmentDetails *> &vec);
};

////////////////////////////////////////////////////////////////

/*! \brief  GSegSCTelescopeFactory concrete class for constructing 
     SC telescope, inherits from GTelescopeFactory to implement
     factory method
 */
class GSegSCTelescopeFactory : public GTelescopeFactory {
 private:

  friend class GReadSegSCStd; 

  GReadSegSCStd *readSegSC;  //!< SC base reader

  GPilot *pi;  //!< pilot reader pointer
  vector<string> tokens;  //!< string vector for GPilot use
  string sPilotEdit;  //!< pilot file string from edit record

  map<int,SegSCStdOptics*> mStdOptics;    //*< map of standard telescopes
  map<int,SegSCStdOptics*>::iterator itmStdOp;  //*< iterator for this map

  map<int, TGraph *> *mGRefl;
  map<int, TGraph *>::iterator itmGRefl;
  
  SegSCStdOptics *opt;  //*< working stdOptics for current telescope
  
  GSegSCTelescope *SCTel;  //*< pointer to working telescope
  int iNumSCTelMade;
  /*! \brief editWorkingTelescope makes edits based on 
           pilotfile entries to telescope currently 
           under construction

           \param SCTel pointer to current telescope
   */
  void editWorkingTelescope(GSegSCTelescope *SCTel1);

 public:

  /*!  GSegSCTelescopeFactory constructs a SC telescope from 
       standard telescopes obtained from reader. USE THIS CONSTRUCTOR

       \param dcReader SC telescope reader instance 
       \param editPilotFile pilotfile name containing telescope edit records
   */
  GSegSCTelescopeFactory(GReadSegSCStd &scReader,
		    const string &editPilotFile);
  /*! \brief makeTelescope constructs a SC telescope based on 
             instructions in the pilot file

             \param id telescope id within array
             \param std number of standard telescope
             \param xLoc x position within array (meters)
             \param yLoc y position within array (meters)
             \param zLoc z position within array (meters)
             \return GSegSCTelescope pointer to constructed telescope
  */

  ~GSegSCTelescopeFactory();

  GSegSCTelescope *makeTelescope(const int &id,
                               const int &std);
  
  /*! \brief printStdTelescope debug print based on prtMode

      \param iStd standard telescope id
      \param mode printmode (see SegSCStdOptics::printSegSCStdOptics
      in this file, default 0
      \param oStr output stream, default cout
   */
  void printStdTelescope(const int &iStd, const int &mode = 0,
                         ostream &oStr=cout);

  /*! \brief setPrintMode defines the print mode for 
          SegSCStdOptics::printSegSCStdOptics (documented in this file)

      \param oStr output stream, default cout
      \param prtMode  see SegSCStdOptics::printSegSCStdOptics earlier in this file
      \param prtMode  -1: print all standard telescopes
   */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
};

#endif

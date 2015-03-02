/*
VERSION3.1
2March2015
*/
/*! \brief GDCTelescopeFactory concrete class for creating ACT 
  Telescopes

Modified Factory design pattern where GTelescopeFactory provides 
the base class for DC and SC telescope concrete factories. 
GDCTelescopeFactory produces the DC telescopes.

*/
#ifndef GDCTELESCOPEFACTORY
#define GDCTELESCOPEFACTORY

// forward declarations
class GPilot;
class GReadDCStdBase; 
class GDCTelescope;
class GTelescope;
class GGeometryBase;
class GDCStdFacet;
class DCStdFacet;
class GOrderedGrid;

/*!  /brief DCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
struct DCStdOptics {

  /*! constructor 
   */
  DCStdOptics(); 

  /*! setPrintOptions debug printer
    only use for changing from *oLog or from iPrtMode in input file
    for standard telescope
    \param outStr output stream,
    \param prtMode 0  no vectors or classes printed in detail
    \param prtMode 1  add geoStructure print
    \param prtMode 2  add facet print
    \param prtMode 3  add reflection coefficient tables
  */

  ~DCStdOptics();

  void setPrintOptions(ostream *outStr, const int &prtMode) {
    oStr = outStr;
    iPrtMode = prtMode;
  }

  /*! printDCStdOptics debug printer
    \param oStr output stream, cout default
    \param prtMode 0  no vectors or classes printed in detail
    \param prtMode 1  add geoStructure print
    \param prtMode 2  add facet print
    \param prtMode 3  add reflection coefficient tables
    
   */
  void printDCStdOptics();

  ostream *oStr;
  int iPrtMode;

  TelType stdType;  
  double fAvgTransitTime;

  double defoc;           /*!< Defocusing >0 when camera is too far 
			   away from the dish <0 when too close.*/
  double foclength;       //!< Focal length (meters)
  double camRadius;       //!< Radius of camera
  double plateScaleFactor; //!< Plate scale factor (cm/deg).
  double radius;          //!< Radius of dish (meters) 
  double rotation_offset; /*!< distance between focal point and
                            rotation center (meters) */
  double xoff;            //!< x pointing offset in radians, coor.sys???
  double yoff;            //!< y pointing offset in radians

  double cameraRadius;    //!< radius of the camera (meters)
 
  int nbr_mir;          //!< Number of mirror elements

  vector<DCStdFacet *> facet; //!< mirror element array 

  int nbinsx;    //!< num of grid bins in x direction
  int nbinsy;    //!< num of grid bins in y direction
  int gridOption; //!< option for grid/file.
  string gridFile; //!< name of gridFile, default ""
  GOrderedGrid *grid; //!< class instance, ordered hash table

  int igeoStd;  /*!< telescope geoStructure stdandard number */

  RayTracerType rayTracerEnum;  //!< ray tracer identifier enum
};

////////////////////////////////////////////////////////////////

/*! \brief  GDCTelescopeFactory concrete class for constructing 
     DC telescope, inherits from GTelescopeFactory to implement
     factory method
 */
class GDCTelescopeFactory : public GTelescopeFactory {
 private:

  friend class GReadDCStdGrISU;

  //int iTelID;   //!<  telescope id within array 

  GReadDCStdBase *readDC;  //!< DC base reader

  GPilot *pi;  //!< pilot reader pointer
  vector<string> tokens;  //!< string vector for GPilot use
  string sPilotEdit;  //!< pilot file string from edit record

  int iNumTelMake;  //!< number of telescopes make by the factory

  map<int,DCStdOptics*> mStdOptics;    //*< map of standard telescopes
  map<int,DCStdOptics*>::iterator itmStdOp;  //*< iterator for this map

  DCStdOptics *opt;  //*< working stdOptics for current telescope

  map<int,GGeometryBase*> mDCGeo;

  //GDCTelescope *DCTel;  //*< pointer to working telescope

  // maps of reflection coeff. wavelgts and coeffs.
  map<int, vector<double> > *mVReflWaveLgts;
  map<int, vector<double> > *mVCoeffs;

  /*! \ brief makeRayTracer adds ray tracer algorithm to DCTel
    \param DCTel pointer to current telescope
   */
  void makeRayTracer(GDCTelescope *DCTel,DCStdOptics *opt1);

  /*! \brief editWorkingTelescope makes edits based on 
           pilotfile entries to telescope currently 
           under construction

           \param DCTel pointer to current telescope
   */
  void editWorkingTelescope(GDCTelescope *DCTel);

 public:

  /*!  GDCTelescopeFactory constructs a DC telescope from 
       standard telescopes obtained from reader. USE THIS CONSTRUCTOR

       \param dcReader DC telescope reader instance 
       \param editPilotFile pilotfile name containing telescope edit records
   */
  GDCTelescopeFactory(GReadDCStdBase &dcReader,
                         const string &editPilotFile);

  ~GDCTelescopeFactory();

  /*! \brief makeTelescope constructs a DC telescope based on 
             instructions in the pilot file

             \param id telescope id within array
             \param std number of standard telescope
             \param xLoc x position within array (meters)
             \param yLoc y position within array (meters)
             \param zLoc z position within array (meters)
             \return GDCTelescope pointer to constructed telescope
  */
  GDCTelescope *makeTelescope(const int &id,
                               const int &std);
  
  /*! \brief printStdTelescope debug print based on prtMode

      \param iStd standard telescope id
      \param mode printmode (see DCStdOptics::printDCStdOptics
      in this file, default 0
      \param oStr output stream, default cout
   */
  void printStdTelescope(const int &iStd, const int &mode = 0,
                         ostream &oStr=cout);

  /*! \brief setPrintMode defines the print mode for 
          DCStdOptics::printDCStdOptics (documented in this file)

      \param oStr output stream, default cout
      \param prtMode  see DCStdOptics::printDCStdOptics earlier in this file
      \param prtMode  -1: print all standard telescopes
   */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
};

#endif

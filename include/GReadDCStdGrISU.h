/*
VERSION3.1
2March2015
*/
/*! \brief  GReadDCStdGrISU class: concrete class for reading 
      standard telescope configurations for use by telescope factory

 */

#ifndef GREADDCSTDGRISU
#define GREADDCSTDGRISU

class GPilot;
class GGeometryBase;

class GReadDCStdGrISU : public GReadDCStdBase {

 private:

  // for reading pilot file, can be appended files
  string spilotfile;   //!< pilot file
  GPilot *pi;        //!< pilot file reader
  vector<string> tokens; //!< vector of pilot record string tokens
  string flagline;       //!< record flag for pilot file reading

  // reader gets a list of telescopes to read
  // and a list of Geometries to read 
  vector<int> vTelNum;     //!< list of telescopes names in TELSTD record
  vector<int> vGeoStd;     //!< list of GeoStandards from TELSTD record

  // reader also gets optics standard number and associated 
  // telnumbers and geostdnumbers from TELSTD record
  map<int,pair<int,int> > mOptics;   //!< map<std#,pair<GrISU Tel#,geom#>>
  map<int,pair<int,int> >::iterator itOpt; //!< map iterator

  map<int,RayTracerType> mRayTracerEnum;  //!< map<std#, raytracerEnum> 

  map<int,int> mPrintStdTel;
  map<int,double> mTransitTime;

  // map of telescope optics key = telescope number
  // an entry in this map for each entry in vTelNum vector
  map<int,DCStdOptics*> mTel;  //!<  map<GrISU tel#,DCStdOptics pointer>

  int iTelNum;  //!< active telescope number used in filling entries
  int iGeoNum; //!< active geometry number, used in filling telescope entries

  double dFocLgt;  //!< temporary storage, focal length
  double dCamRad;  //<! temporary storage, camera radius
  double dPlateScaleFactor; //<! temporary storage
  double dFocError; //!< temporary for now
  double dTelRadius; //!< temporary

  /*! \brief setupDCFactory: creates all standard array telescopes
               from GrISU telescope objects. Results loaded into the
               DCTelescopeFactory object *DCFac.

   */
  void setupDCFactory();

  /*! \brief makeStdTelescope constructs a standard telescope
             from GrISU configuration file and from GDCGeometry 
             object
   */
  void makeStdTelescope();

  /*! \brief  makeStdGeometry constructs GDCGeometry objects from
                  geometry configuration file

   */
  void makeStdGeometry();

  /*! \brief  getReflCoeff reads reflection table from GrISU
              configuration file

   */
  void getReflCoeff();

  /*! \brief makeStdOptics fills standard optics instance of
             DCStdOptics class

   */
  void makeStdOptics();

 public:

   /*! \brief GReadDCStdGrISU constructor
     \param pilotfile configuration file with GrISU configuration-type
            entries plus table showing std#,GrISUTel#,geom#,
            and, ray-tracing id
     \param DCFactory instance of GDCTelescopeFactory
   */
 GReadDCStdGrISU(const string pilotfile,GDCTelescopeFactory *DCFactory );

  /*! \brief GReadDCStdGrISU constructor
      \param pilotfile configuration file with GrISU configuration-type
            entries plus table showing std#,GrISUTel#,geom#,
            and, ray-tracing id

   */
  GReadDCStdGrISU(const string pilotfile);

  /*! \brief ~GReadDCStdGrISU destructor

   */
  ~GReadDCStdGrISU();

  /*! \brief setDCTelescopeFactory pass in instance of
          GDCTelescopeFactory if not passed in constructor
      \param DCFactory instance of GDCTelescopeFactory
   */
  void setDCTelescopeFactory(GDCTelescopeFactory *DCFactory);  

  /*! \brief setPrintMode sets the debug printmode (still needs work)
       \param oStr output stream for printer (default cout)
       \param prtMode printmode, if 0, no printing

   */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);

};

#endif


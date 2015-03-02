/*
VERSION3.1
2March2015
*/
/*! \brief  GDCGeometry concrete class for holding telescope
  structure parameters for VERITAS-like DC telescope structures.
  Other types of DC structures require individual concrete classes.
  When used, the base class must be downcasted to this class in order
  to extract the data.

*/

#ifndef GDCGEOMETRY
#define GDCGEOMETRY

class GDCGeometry : public GGeometryBase {

 public:
  // type is declared in base class

  double epsil;  //!< addition to top volume (I think)

  // focus box
  double focBoxDim[3];  //!< focus box dimensions
  double focBoxRot[3];  //!< focus rotation angles
  
  // edge boxes (4)
  double edgeX[4];      //!<  edge boxes x dimensions
  double edgeY[4];      //!<  edge boxes y dimensions
  double edgeZ[4];      //!<  edge boxes z dimensions
  double edgeOffset[4]; //!<  edge box offsets, normally 0.0
  double edgeRot1[4];   //!<  first rotation (degrees) edge boxes
  double edgeRot2[4];   //!<  second rotation (degrees) edge boxes
  double edgeRot3[4];   //!<  third rotation (degrees) edge boxes

  // shutter
  double shutterX;     //!< shutter x dimension  
  double shutterZ;     //!< shutter z dimension  
  double shutterRot1;  //!< shutter rotation 1  
  double shutterRot2;  //!< shutter rotation 2  
  double shutterRot3;  //!< shutter rotation 3  

  // quad arms
  double quadArmX;       //!< quadArm x dimension
  double quadArmY;       //!< quadArm y dimension 
  double quadArmOffset;  //!< quadArm offset

  double quadArmBottomX[4];  //!< quad arm bottom location X
  double quadArmBottomY[4];  //!< quad arm bottom location Y

  double crossBarX[2];
  double crossBarY[2];
  double crossBarDistBelowFocBox[2];

  double cameraRadius;  //!< camera radius

  /*! \brief GDCGeometry constructor
   */
  GDCGeometry();

  /*! \brief GDCGeometry copy constructor
   */
  GDCGeometry(const GDCGeometry &c);

  /*! \brief GDCGeometry deconstructor
   */
  ~GDCGeometry();

  /*! \brief getType method returns geotype enum
      \return GeoType enum for geometry type
   */
  GeoType getType();

  /*! \brief printGeometry prints geometry details
   */
  void printGeometry();

  /*! \brief setPrintMode sets stream and mode for printGeometry()

      \param oStr output stream
      \param prtMode print mode, to be documented later
   */
  void setPrintMode(ostream &oStr=cout,const int prtMode=0);
  
  /*! \brief printGeometry1 prints geometry details (TODO eliminate
         one of these printGeometry methods)

       \param oStr output stream
       \param prtMode print mode, to be documented later   
   */
  void printGeometry1(ostream &oStr=cout,const int &prtMode=0);

};

#endif

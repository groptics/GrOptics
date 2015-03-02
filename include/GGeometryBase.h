/*
VERSION3.1
2March2015
*/
/*!  GGeometryBase base class for holding telescope structure data
  for use in ROOT geometry and navigation classes. No facet information
  contained in this or subsequent concrete class.  This is a pure virtual
  class. Downcasting later makes available the different telescope structure 
  information; e.g. differing number of quad arms, etc.

*/

#ifndef GGEOMETRYBASE
#define GGEOMETRYBASE

class GGeometryBase {

 protected:

  // used by all telescope types for various printing modes
  ostream *oPrtStrm;
  int iPrtMode;

 public:
  GeoType type;

  /*!  \brief GGeometryBase constructor
   */
  GGeometryBase();
  
  /*!  \brief GGeometryBase virtual deconstructor
   */
  virtual ~GGeometryBase();

  /*!  \brief getType pure virtual method
       \return GeoType geoType Enum of Geometry type
   */
  virtual GeoType getType() = 0;

  /*!  \brief printGeometry pure virtual method to print 
              geometry parameters
   */
  virtual void printGeometry() = 0;

  /*!  \brief setPrintMode pure virtual method to set 
              print parameters 
       \param oStr print output stream
       \param prtMode determines printing, to be documented later
   */
  virtual void setPrintMode(ostream &oStr=cout,const int prtMode=0) = 0;

  /*!   \brief printGeometry1 pure virtual method to print
              geometry parameters
        \param oStr output stream
        \param prtMode determines printing, to be documented later
   */
  virtual void printGeometry1(ostream &oStr=cout,const int &prtMode=0) = 0;
};

#endif

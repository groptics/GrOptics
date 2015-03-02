/*
VERSION3.1
2March2015
*/
/*! \brief GTelescopeFactory base class for creating ACT 
  Telescopes

GTelescopeFactory provides the base class for DC and SC telescope
concrete factories, and for various readers

*/
#ifndef GTELESCOPEFACTORY
#define GTELESCOPEFACTORY

class GTelescope;  //!<  forward declaration

class GTelescopeFactory {
 private:

 protected:

  ostream *oPrtStrm; //<! ostream used by printing methods
  int iPrtMode; //<! print mode used by printing methods

 public:

  /*!  \brief GTelescopeFactory constructor

   */
  GTelescopeFactory();

  /*!  \brief ~GTelescopeFactory virtual destructor

   */
  virtual ~GTelescopeFactory() { };

  /*!  \brief makeTelescope pure virtual method to construct 
             a telescope based on TODO ADD DC/SC?

             \param id telescope id within array
             \param std number of standard telescope
             \param xLoc x position within array (meters)
             \param yLoc y position within array (meters)
             \param zLoc z position within array (meters)
             \return GDCTelescope pointer to constructed telescope

   */
  virtual GTelescope *makeTelescope(const int &id,
                                       const int &std) = 0;

  /*!  \brief printStdTelescope print pure virtual method

      \param iStd standard telescope id
      \param mode printmode (see above)
      in this file, default 0
      \param oStr output stream, default cout

   */
  virtual void printStdTelescope(const int &iStd, const int &mode = 0,
                                 ostream &oStr=cout) = 0; 

  /*!  \brief setPrintMode pure virtual method defines the print mode

      \param oStr output stream, default cout
      \param prtMode  print mode

   */
  virtual void setPrintMode(ostream &oStr=cout,const int prtMode=0) = 0;
};

#endif

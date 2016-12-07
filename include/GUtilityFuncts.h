/*
VERSION4.0
30May2016
*/
/*!  GUtilityFuncts
     utility functions for use in new optics code
         Version 0.0
           C. Duke
           Grinnell College
           July 13, 2010

  These utility functions fall under namespace GUtilityFuncts
*/

#ifndef GUTILITYFUNCTS
#define GUTILITYFUNCTS

// forward declarations
#include "Math/Rotation3D.h"
#include "Math/Vector3D.h"

extern ostream *oLog;
struct mirrorSegmentDetails;

namespace GUtilityFuncts {
  
  /*! \brief polyInside function of testing regular polygon pixel or
    facet shapes

    bool polyInside(const int &sides, const double &alpha, 
    const double &radius, 
    const double &x, const double &y) 
    
    poly_inside determines whether or not point, (x,y) with origin at the center
    of the polygon, lies inside a regular polygon described by:
    sides:  number of sides for the polygon
    radius: distance from the center of the polygon to one of its points
    alpha:  orientation of the polygon. When alpha is zero 
    a point of the polygon is along the y-axis. 
    For alpha > 0, alpha is the clockwise rotation of the 
    polygon about its center through angle
    alpha (in radians).(last tested 08/18/05 C. Duke)
    
    return: true if (x,y) inside polygon or false if not 
    
    if sides = 0, always return 0
    if sides = 1 or 2, determine if (x,y) lies with a circle of r= radius.
    if sides > 2, sides = number of sides of regular polygon
  */    
  
  bool polyInside(const int &sides, const double &alpha, 
                  const double &radius, 
                  const double &x, const double &y); 



  /*! \brief tokenizer function, converts to string to token vector.
       delimiters are blank and comma

     \param tokens vector for token strings
     \param str string for conversion to tokens
   */
  void tokenizer(const string& str, vector<string>& tokens);
 
  /*! \brief convert matlab string, [n:m] or [n m o p] to vector of 
    signed integers or if string is a single integer, 
    e.g. n, put n into the vector.

   */
  bool decodeMatlabString(const string &matlabStr,
                          vector<int> &numRead);

  /*! \brief convert matlab string, [n:m] or [n m o p] to vector of 
    unsigned integers
    or if string is a single integer, e.g. n, put n into the vector.

   */
  bool decodeMatlabString(const string &matlabStr,
                        vector<unsigned> &numRead);

  /*! \brief tokenize string with leading asterisk, omit the asterisk

   */
  bool getAsteriskTokens( const string &str, vector<string> & tokens);

  /*! \brief tokenize string, keeping intact token with matlab format

   */
   void tokenizeMatlab(const string &str, vector<string> &tokens);

 /*! \brief template function for printing a vector of any type
   example: printVector<string> vector;

  */
   template<class T>
     void printVector(const vector<T> &vec) 
     {
       for (unsigned i = 0;i< vec.size();i++) {
         *oLog << i << "   " << vec[i] << endl;
       }
     };

   /*! \brief template function for printing a root XYZVector or any
     3D vector that uses methods X(), Y(), Z() to access components.

     \param vec1 vector with components vec1.X(), vec1.Y(), and vec1.Z()

    */
   template<class Vector1> 
     void printGenVector(const Vector1 &vec1) 
     {
       *oLog << setw(12) << setprecision(8) << showpoint << vec1.X() << " " 
            << vec1.Y() << "   " << vec1.Z();
     };
   
   //string evalTelType(const TelType &teltype);

   /*! \brief template function for linear interpolation, using two vectors
     using double datatype, returns -1.0 if outside range 
    
    */
   double linearInterpolation2D(const vector<double> &xvec,
                                const vector<double> &yvec,const double xval);
   
   /*! \brief add a random misallignment vector to a unit vector.  The
    angle, max_angle, is the maximum conical angle between the unit vector 
    and the misaligned unit vector. 

    The misalignment directions have a constant density in the circular x/y space 
    perpendicular to the given unit vector.  Root functions are not used in this code.
    Not currently used in GrOptics. See GUtilityFuncts::addErrorRoot.
    \param v[3] given unit vector
    \param maxAngle maximum conical angle between the given unit vector and the 
           misaligned unit vector.
   */
   void addError(double v[3], const double &max_angle);

   /*! \brief multiply a vector by a matrix, v2 =  A*v1

     \param  a matrix of size 3 x 3
     \param  b 3D vector
     \param  c 3D vector result of matrix multiplication a*b.

    */
   void mtxMlt(double a[3][3], double b[3], double c[3]);

   /*! \brief Add a random misallignment vector to a unit vector. The
    angle, max_angle, is the maximum conical angle between the unit vector 
    and the misaligned unit vector. 

    The misalignment directions have a constant density in the circular x/y space 
    perpendicular to the given unit vector.
    \param vec pointer to the given unit vector
    \param maxAngle maximum conical angle between the given unit vector and the 
           misaligned unit vector.
   */
   void addErrorRoot(ROOT::Math::XYZVector *vec, 
                     const Double_t &maxAngle);

   /*! \brief determine (r,theta}, where r < rMax and 
        theta in radians, for a circular area with a
        Gaussian distribued probability density proportional
        to exp(-r^2/sigma^2). 

        Calculation based on probability transformation 
	to a uniformly distribued deviate.
	\param rMax  maximum value of r in (r,theta) plane
 	\param sigma Gaussian parameter in exp(-r^2/sigma^2)
	\param r     pointer to r value in (r,theta) plane
	\param theta pointer to theta value in (r,theta) plane
   */
   void getGaussianErrorPolarCoor(const double &rMax,
                                  const double &sigma,
                                  double *r,double *theta);

   /*!  \brief adds error to unit vector, Gaussian error with
     maximum Gaussian offset specified

     \param vec   pointer to XYZVector unit vector
     \param rMax  maximum Gaussian offset
     \param sigma Gaussian width in exp(-r^2/sigma^2)
   */
   void addErrorGaussianRoot(ROOT::Math::XYZVector *vec, 
                             const double &rMax,
                             const double &sigma);

   /*!  returns angle between the sky location to zenith arc
     and the sky location to polaris arc on the celestial
     sphere.
     
     \param zn zenith angle (radians) of sky location 
     \param az azimuthal angle (radians) of sky location
     \param latitude latitude (radians) of the observatory
     \return angle(radians) of Polaris with respect to zenith.
   */
   double fieldRot(const double &zn, const double &az,
                   const double &latitude);

   /*! \brief find azimuth and zenith angles from given x and y direction cosines

     \param xcos x-direction cosine of vector pointing to azimuth and zenith defined point
     \param ycos y-direction cosine of vector pointing to azimuth and zenith defined point
     \param[out] az  pointer to azimuthal angle determined from xcos and ycos
     \param[out] zn  pointer to zenith angle determined from xcos and ycos.
     \return 0
    */
   int XYcosToAzZn(const double &xcos, const double &ycos,
                                      double *az,double *zn);
   /*! \brief find unit vector pointing to given azimuth and zenith angles defined point on sphere.
     
     \param az  azimuthal angle of point on unit sphere
     \param zn  zenith angle of point on unit sphere
     \param xcos pointer to x direction cosine of vector to point defined by az and zn
     \param ycos pointer to y direction cosine of vector to point defined by az and zn.
     \return 0
    */
   int AzZnToXYcos(const double &az,const double &zn,
                   double *xcos, double *ycos);

   /*! \brief  Determines telescope azimuth and zenith angles from source location
     and telescope offsets/wobble angles using small angle approximations
     
     All angles in radians
     \param primAz       source azimuth
     \param primZn       source zenith angle
     \param wobbleN      telescope wobble North (toward Polaris)
     \param wobbleE      telescope wobble East (perpendicular to wobbleN)
     \param telOffsetX   telescope offset along X axis of telescope system
     \param tesOffsetY   telescope offset along Y axis of telescope system
     \param latitude     observatory latitude
     \param telAz        pointer to calculated telescope azimuth
     \param telZn        pointer to calculated telescope zenith angle
     \param sourceX      pointer to source location along X axis of telescope system
     \param sourceY      pointer to source location along Y axis of telescope system

    */
   void telescopeAzZn(const double &primAz,
                      const double &primZn,
                      const double &wobbleN,
                      const double &wobbleE,
                      const double &telOffsetX,
                      const double &telOffsetY,
                      const double &latitude,
                      double *telAz,double *telZn,
                      double *sourceX,double *sourceY);

   /*! \brief  Determines telescope azimuth and zenith angles from source location
     and telescope offsets/wobble angles with no small angle approximations
     
     All angles in radians
     \param primAz       source azimuth
     \param primZn       source zenith angle
     \param wobbleN      telescope wobble North (toward Polaris)
     \param wobbleE      telescope wobble East (perpendicular to wobbleN)
     \param telOffsetX   telescope offset along X axis of telescope system
     \param tesOffsetY   telescope offset along Y axis of telescope system
     \param latitude     observatory latitude
     \param telAz        pointer to calculated telescope azimuth
     \param telZn        pointer to calculated telescope zenith angle
     \param sourceX      pointer to source location along X axis of telescope system
     \param sourceY      pointer to source location along Y axis of telescope system
    */
   void telescopeAzZnNew(const double &primAz,
			 const double &primZn,
			 const double &wobbleN,
			 const double &wobbleE,
			 const double &telOffsetX,
			 const double &telOffsetY,
			 const double &latitude,
			 double *telAz,double *telZn,
			 double *sourceX,double *sourceY);
         
   /*! \brief Determines the rotation matrix from ground to telescope coordinates
     
     Usage: vecTelCoor = (*rotM)*vecGrdCoor where vecGrdCoor and vecTelCoor 
     are ROOT::Math::XYZVector vectors in ground and telescope coordinates.
     All angles in radians
     \param az    azimuth of point on unit sphere
     \param zn    zenith angle of point on unit sphere
     \param rotM  Rotation3D matrix to transform a vector from gr.coor. to telescope coordinates
    */
   void AzZnToRotMat(const double &az,
                     const double &zn,
                     ROOT::Math::Rotation3D *rotM);
  
   /*! \brief  Determine the azimuth and zenith of the offset point with coordinates 
      (offsetX, offsetY) with respect to tangent plane telescope axes.

      All angles are in radians.
      \param offsetX   tangent plane offset in X direction in telescope system
      \param offsetY   tangent plane offset in Y direction in telescope system
      \param az        azimuth of tangent plane origin
      \param zn        zenith angle of tangent plane origin
      \param azOffset  azimuth of offset point
      \param znOffset  zenith angle of offset point
      
    */
   void offsetXYToAzZn(const double &offsetX, const double &offsetY,
                       const double &az, const double &zn,
                       double *azOffset, double *znOffset);

   /*! \brief Determines source location on telescope plane given az/zn
     angles of the source and the telescope.

     All angles in radians.
     \param az_s azimuth of source
     \param zn_s zenith angle of source
     \param az_t azimuth of telescope
     \param zn_t zenith angle of telescope
     \param x    pointer to source x location in telescope coordinates
     \param y    pointer to source y location in telescope coordinates
    */
   void sourceOnTelescopePlane(const double &az_s, const double &zn_s,
			       const double &az_t, const double &zn_t,
			       double *x, double *y);
   
   /*! \brief Tests if a datatype double is within numeric_limits<double> 
     of zero. If true, sets value to 0.0 if zeroFlag == 1.

     Uses TMath::AreEqualAbs(*ax,0.0,numeric_limits<float>::epsilon())
     for zero test.  
     \param ax pointer to datatype double variable
     \param zeroFlag  set ax to 0.0 if zero test is true
     \return result of zero test.

    */
   bool testZeroFloat(double *ax, const int &zeroFlag);

   /*! \brief Tests if vector components are within numeric_limits<double>
     of 0.0.  If true, sets the component to 0.0.

     Uses TMath::AreEqualAbs(vec->X(),0.0,numeric_limits<double>::epsilon())
     for zero test.
     \param vec pointer to a ROOT::Math::XYZVector.
     \return true if any component is set to 0.0
    */
   bool zeroFloatVectorFix(ROOT::Math::XYZVector *vec);

   /*! \brief Determines if photon will reflect from a surface. 
     Uses reflection coefficient table plus random-number comparison
     
     Currently used with mirror elements of DC telescope model
     \param refCoeff   vector of reflection coefficients 
     \param refWaveLgt vector of wavelengths for coefficients 
     \param reflect    reflection degradation factor
     \param photonWaveLgt wavelength of incident photon
     \return true if photon reflects from surface.
    */
   bool photonReflect(const vector<double> &refCoeff,
                      const vector<double> &refWaveLgt,
                      const double &reflect,
                      const double &photonWaveLgt );
                         
   /*! \brief Determines direction of reflected photon. Surface roughness added 
     with GUtilityFuncts::addErrorRoot(vphotonReflDcos,roughness).

     \param vnormUnit direction cosine unit vector: normal to surface
     \param vphotonUnit direction cosine unit vector: photon direction 
     \param vphotonReflDcos direction cosine unit vector; reflected photon 
     \param roughness surface roughness
     \param option not used
    */
   void reflectDirection(const ROOT::Math::XYZVector &vnormUnit,
                         const ROOT::Math::XYZVector &vphotonUnit,
                         ROOT::Math::XYZVector *vphotonReflDcos,
                         const double &roughness, const int &option=0);

   /*! \brief criteria for sort function, returns true if j.second < i.second

     \param i  pair of integers
     \param j  pair of integers
     \return true if j.second < i.second
    */
   bool sortPair(const pair<int,int> i , const pair<int,int> j);

   /*! \brief prints to *oLog labeled table of characteristics of all mirror
     segments of DC telescope.

     \param vec pointer to vector of mirrorSegmentDetails structures
    */
   void printSegVector (const vector<mirrorSegmentDetails *> &vec);

   /*! \brief prints to *oLog a label followed by the three components of vec, 
     a ROOT::Math::XYZVector

     \param vec a ROOT::Math::XYZVector
     \param label the string initially printed that labels the vector
    */
   void printXYZVector(const ROOT::Math::XYZVector &vec,const string &label);

};
#endif

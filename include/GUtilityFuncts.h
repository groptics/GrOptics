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
   *   delimiters are blank and comma

     \param tokens vector for token strings
     \param str string for conversion to tokens
   */
  void tokenizer(const string& str, vector<string>& tokens);
 
  /*! \brief convert matlab string, [n:m] or [n m o p] to vector of 
   * signed integers
   *  or if string is a single integer, e.g. n, put n into the vector.

   */
  bool decodeMatlabString(const string &matlabStr,
                          vector<int> &numRead);

  /*! \brief convert matlab string, [n:m] or [n m o p] to vector of 
   * unsigned integers
   *  or if string is a single integer, e.g. n, put n into the vector.

   */
  bool decodeMatlabString(const string &matlabStr,
                        vector<unsigned> &numRead);

  /*! \brief tokenize string with leading asterisk, omit the asterisk
   *   
   *

   */
  bool getAsteriskTokens( const string &str, vector<string> & tokens);

  /*! \brief tokenize string, keeping intact token with matlab format
   *

   */
   void tokenizeMatlab(const string &str, vector<string> &tokens);


 /*! \brief template function for printing a vector of any type
  * example: printVector<string> vector;

  */
   template<class T>
     void printVector(const vector<T> &vec) 
     {
       for (unsigned i = 0;i< vec.size();i++) {
         *oLog << i << "   " << vec[i] << endl;
       }
     };
   
   
   template<class Vector1> 
     void printGenVector(const Vector1 &vec1) 
     {
       *oLog << setw(12) << setprecision(8) << showpoint << vec1.X() << " " 
            << vec1.Y() << "   " << vec1.Z();
     };
   
   //string evalTelType(const TelType &teltype);

   /*! \brief template function for linear interpolation, using two vectors
    *  using double datatype, returns -1.0 if outside range 
    
    */
   double linearInterpolation2D(const vector<double> &xvec,
                                const vector<double> &yvec,const double xval);
   
   /*! \brief add a random misallignment vector to a unit vector.  The
    angle, max_angle, is the maximum conical angle between the unit vector 
    and the misaligned unit vector. 
   */
   void addError(double v[3], const double &max_angle);

   /*! \brief multiply a vector by a matrix, v2 =  A*v1
    */
   void mtxMlt(double a[3][3], double b[3], double c[3]);


   void addErrorRoot(ROOT::Math::XYZVector *vec, 
                     const Double_t &maxAngle);

   /*! \brief determine (r,theta}, where r < rMax and 
        theta in radians, for a circular area with a
        Gaussian distribued probability density proportional
        to exp(-r^2/sigma^2). 

        Calculation based on probability
        transformation to a uniformly distribued deviate. 
    */
   void getGaussianErrorPolarCoor(const double &rMax,
                                  const double &sigma,
                                  double *r,double *theta);

   /*!  \brief adds error to unit vector, gaussian error with
     maximum Gaussian offset specified

   */
   void addErrorGaussianRoot(ROOT::Math::XYZVector *vec, 
                             const double &rMax,
                             const double &sigma);


   /*!  returns angle between sky location to zenith arc
     and the sky location to polaris arc on the celestial
     sphere.
     
     \param zn zenith angle (radians) of sky location 
     \param az azimuthal angle (radians) of sky location
     \param latitude latitude (radians) of the observatory
     \return angle(radian) of Polaris with respect to zenith.
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

   /*! \brief


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

   /*! \brief


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
         
   /*! \brief using this


    */
   void AzZnToRotMat(const double &az,
                     const double &zn,
                     ROOT::Math::Rotation3D *rotM);
  
   /*! \brief using this


    */
   void offsetXYToAzZn(const double &offsetX, const double &offsetY,
                       const double &az, const double &zn,
                       double *azOffset, double *znOffset);

   /*! \brief


    */
   void tangentPlaneOffset( const double &az  , const double &zn,
                            const double &az_t, const double &zn_t,
                            double       *xoff, double       *yoff,
                            int          *calc_status   ) ;
   
   void tangentPlaneOffsetNew( const double &az  , const double &zn,
                            const double &az_t, const double &zn_t,
                            double       *xoff, double       *yoff,
                            int          *calc_status   ) ;


   void sourceOnTelescopePlane(const double &az_s, const double &zn_s,
			       const double &az_t, const double &zn_t,
			       double *x, double *y);
   
   /*! \brief


    */
   bool testZeroFloat(double *ax, const int &zeroFlag);

   /*! \brief


    */
   bool zeroFloatVectorFix(ROOT::Math::XYZVector *vec);

   /*! \brief


    */
   bool photonReflect(const vector<double> &refCoeff,
                      const vector<double> &refWaveLgt,
                      const double &reflect,
                      const double &photonWaveLgt );
                         

   /*! \brief


    */
   void reflectDirection(const ROOT::Math::XYZVector &vnormUnit,
                         const ROOT::Math::XYZVector &vphotonUnit,
                         ROOT::Math::XYZVector *vphotonReflDcos,
                         const double &roughness, const int &option=0);

   /*! \brief


    */
   bool sortPair(const pair<int,int> i , const pair<int,int> j);

   /*! \brief


    */
   void printSegVector (const vector<mirrorSegmentDetails *> &vec);

   /*! \brief


    */
   void printXYZVector(const ROOT::Math::XYZVector &vec,const string &label);

};
#endif

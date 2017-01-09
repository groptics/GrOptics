/*
VERSION4.0
30May2
*/
/* GUtilityFuncts.cpp : set of utility functions
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <limits>

using namespace std;

#include "TMatrixD.h"
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/EulerAngles.h>
#include <Math/VectorUtil.h>
#include <TRandom.h>
#include <TRandom3.h>

using namespace ROOT::Math;

#define DEBUG(x) *oLog << #x << " = " << x << endl

#include "GDefinition.h"
#include "GUtilityFuncts.h"

bool GUtilityFuncts::polyInside(const int &sides, const double &alpha, 
                                  const double &radius, 
                                  const double &x, const double &y) {
  
  double PI = TMath::Pi();
  double theta,dist;
  double phi;

  static double delta[50],beta[50],gamma[50];
  static int ifirst=0;
  int iside;
  double x1,y1;

  /* initialize constant parameters */
  if (ifirst==0) {
    ifirst = 1;
    delta[0]=delta[1]=delta[2]=0.0;
    beta[0] = beta[1]= beta[2]=0.0;
    gamma[0]=gamma[1]=gamma[2]=0.0;

    for (iside=3;iside<50;iside++) {
      delta[iside] = (PI)*( 0.5 - (1/(double)iside) );
      gamma[iside] = 2*PI/(double)iside;
      beta[iside]  = gamma[iside]/2.;
    }
  }
  /* number of valid sides determined by size of arrays 
     I only went to 49.............. Maybe that's a circle....
  */
  if (sides > 49) {
    *oLog << " number of sides greater than 49! " << endl;
    *oLog << " exiting code" << endl;
    exit(0);
  }

  /* if at least three sides, can proceed */
  /*-------------------------------------------*/
  if (sides > 2) {

    dist =  sqrt((x*x) + (y*y));

    /* first see if the distance isn't too large */
    if (dist > radius) return 0;

    /* rotate x and y through alpha, dist doesn't change */
    x1 = (cos(alpha)*x)- (sin(alpha)*y);
    y1 = (sin(alpha)*x)+ (cos(alpha)*y);

    /* if you use atan to get theta, phi will be incorrect */
    theta = acos(y1 / dist);
  
    phi = fabs( (fmod((theta + beta[sides]), gamma[sides])) - beta[sides]); 

    if (dist<=((sin(delta[sides])*radius)/sin(delta[sides] + phi))) 
      return true;
    
  }
  /*----------------------------------------------*/
  /* if sides = 1 or 2, then shape is a circle rather than a polygon */

  else if (sides > 0 ) {
    dist =  sqrt(x * x + y * y);
    if (dist <= radius) return 1;
  }

  /* if sides = 0, then always return 0 */
  return false;
};
/********************  end of polyInside ***************/

void GUtilityFuncts::tokenizer(const string& str, 
                                     vector<string>& tokens) {
  tokens.clear();
  // tokenizer with multiple delimiters
  char delim[] = " ,";
  string::size_type lgt = 2;
  
  string::size_type idx = 0;
  string::size_type start;

  // find first non-delimiter character for start
  while ( (idx = str.find_first_not_of(delim,idx,lgt) ) !=string::npos ) {
    start = idx; // start of string

    idx = str.find_first_of(delim,idx,lgt);

    tokens.push_back(str.substr(start,idx - start));
  }
};
/********************  end of tokenizer ************/

bool GUtilityFuncts::decodeMatlabString(const string &matlabStr,
                                           vector<int> &numRead) {

  size_t d1,d2,d3;
  vector<string>tokens;
 
  // first look for matlab symbols
  d1 = matlabStr.find("[");
  d2 = matlabStr.find(":"); 
  d3 = matlabStr.find("]"); 
  
  if ( (d1!= string::npos) &&
       (d2!= string::npos) &&
       (d3!= string::npos) ) {
         
    // we have a colon delimited integers
    string loStr = matlabStr.substr(d1+1,d2-d1-1);
    string hiStr = matlabStr.substr(d2+1,d3-d2-1);

    tokens.clear();
    GUtilityFuncts::tokenizer(loStr,tokens);   
    if (tokens.size() != 1) {
      *oLog << "incorrect matlab format: " << loStr 
           << " in " << matlabStr << endl;
      exit(0);
    }
    tokens.clear();
    GUtilityFuncts::tokenizer(hiStr,tokens);   
    if (tokens.size() != 1) {
      *oLog << "incorrect matlab format: " << hiStr 
           << " in " << matlabStr << endl;
      exit(0);
   }
    tokens.clear();
     

    int lo = atoi(loStr.c_str());
    int hi = atoi(hiStr.c_str());
    
    for (int i = lo;i <= hi; i++) {
      numRead.push_back(i);
    }
    
  }
  else if ( (d1 != string::npos) &&
            (d3 != string::npos) && 
            (d2 == string::npos) ) {
    // we have a list of space delimited integers
    string list = matlabStr.substr(d1+1,d3-d1-1);
    *oLog << "list: " << list << endl;
    
    // tokenize and read into the vector
    GUtilityFuncts::tokenizer(list,tokens);
    for (unsigned i = 0;i< tokens.size(); i++) {
      numRead.push_back(atoi(tokens[i].c_str()) );
    }
    
    return true;
  }
  
  else if ( (d1 == string::npos) &&
            (d3 == string::npos) && 
            (d2 == string::npos) ) {
    // this should be a single integer, but check first
    tokens.clear();
    GUtilityFuncts::tokenizer(matlabStr,tokens);
    if (tokens.size() != 1) {
      *oLog << "incorrect matlab format: " 
           << matlabStr << endl;
      exit(0);
    }
    tokens.clear();
    

    int singleInt = atoi(matlabStr.c_str());
    numRead.push_back(singleInt);
  }
  else {
      *oLog << "incorrect matlab format: " 
           << matlabStr << endl;
      exit(0);
  }
  return true;
};
/************* end of decodeMatlabString **************************/

bool GUtilityFuncts::decodeMatlabString(const string &matlabStr,
                                           vector<unsigned> &numRead) {

  *oLog << " -- GUtilityFuncts::decodeMatlabString" << endl;
  *oLog << "       " << matlabStr << endl;

  size_t d1,d2,d3;
  vector<string>tokens;
 
  // first look for matlab symbols
  d1 = matlabStr.find("[");
  d2 = matlabStr.find(":"); 
  d3 = matlabStr.find("]"); 
  
  if ( (d1!= string::npos) &&
       (d2!= string::npos) &&
       (d3!= string::npos) ) {
         
    // we have a colon delimited integers
    string loStr = matlabStr.substr(d1+1,d2-d1-1);
    string hiStr = matlabStr.substr(d2+1,d3-d2-1);

    tokens.clear();
    GUtilityFuncts::tokenizer(loStr,tokens);   
    if (tokens.size() != 1) {
      *oLog << "incorrect matlab format: " << loStr 
           << " in " << matlabStr << endl;
      exit(0);
    }
    tokens.clear();
    GUtilityFuncts::tokenizer(hiStr,tokens);   
    if (tokens.size() != 1) {
      *oLog << "incorrect matlab format: " << hiStr 
           << " in " << matlabStr << endl;
      exit(0);
   }
    tokens.clear();
     

    int lo = atoi(loStr.c_str());
    int hi = atoi(hiStr.c_str());
    
    for (int i = lo;i <= hi; i++) {
      numRead.push_back(i);
    }
    
  }
  else if ( (d1 != string::npos) &&
            (d3 != string::npos) && 
            (d2 == string::npos) ) {
    // we have a list of space delimited integers
    string list = matlabStr.substr(d1+1,d3-d1-1);
    *oLog << "list: " << list << endl;
    
    // tokenize and read into the vector
    //    have to remove the [ and ] first

    GUtilityFuncts::tokenizer(list,tokens);
    for (unsigned i = 0;i< tokens.size(); i++) {
      numRead.push_back(atoi(tokens[i].c_str()) );
    }
    
    return true;
  }
  
  else if ( (d1 == string::npos) &&
            (d3 == string::npos) && 
            (d2 == string::npos) ) {
    // this should be a single integer, but check first
    tokens.clear();
    GUtilityFuncts::tokenizer(matlabStr,tokens);
    if (tokens.size() != 1) {
      *oLog << "incorrect matlab format: " 
           << matlabStr << endl;
      exit(0);
    }
    tokens.clear();
    

    int singleInt = atoi(matlabStr.c_str());
    numRead.push_back(singleInt);
  }
  else {
      *oLog << "incorrect matlab format: " 
           << matlabStr << endl;
      exit(0);
  }
  return true;
};
/************************end of decodeMatlabString********/

bool GUtilityFuncts::getAsteriskTokens(const string &str, 
                                          vector<string> & tokens) {
  // make sure that tokens is an empty vector
  tokens.clear();

  // do we have an asterisk

  if ( str.find("*") != string::npos ) {
    // tokenize the string if we have an asterisk

    // do we have a string containing a matlab formatted entry
    // look for left bracket 
    if ( str.find("[") != string::npos) {
      GUtilityFuncts::tokenizeMatlab(str,tokens);
    }
    else {
    GUtilityFuncts::tokenizer(str, tokens);
    }

    // erase the asterisk; if there is no space after the asterisk,
    // the remaining string is kept.

    if (tokens[0].size() == 1) {
      tokens.erase(tokens.begin());  // erase the entire vector element
    }
    else {
      tokens[0].erase(tokens[0].begin()); // erase the asterisk only
    }
    
   return true;   // successful tokenize asterisk record
  }                         
  return false;  // there was no asterisk
};
/*************** end of getAsteriskTokens ************/

void GUtilityFuncts::tokenizeMatlab(const string &str, vector<string> &tokens) {

  tokens.clear();

  string strc = str;
  
 // string can contain matlab entries, e.g. [1 3 5], so simple
  // tokenizer doesn't work.

  char chB[] = " ";
  char chL[] = "[";
  char chR[] = "]";
  string::size_type lgt = 1;

  string::size_type idx = 0;

  while ( (idx = strc.find_first_not_of(chB,idx,lgt) ) !=string::npos ) {
    // get position of next non-blank character, start of token
    string::size_type start = idx;
    string::size_type endt;
    // is this token a [, matlab bracket
    if ( strc.substr(idx,1) == chL ) {

     // get position of ending bracket
      if ( (idx = strc.find_first_of(chR,idx+1,lgt)) == string::npos)  {
        exit(0);
      }
      endt = idx;
      idx++;
    }
    else {
      // get position of next blank character
      idx = strc.find_first_of(chB,idx+1,lgt);
      endt = idx-1;
     }

    string tokenstr = strc.substr(start,endt-start+1);
    tokens.push_back(tokenstr);

  }    
};
/********************* end of tokenizeMatlab *****************/

double GUtilityFuncts::linearInterpolation2D(const vector<double> &xvec,
                                             const vector<double> &yvec,
                                             const double xval) {
  // Currently unused and untested
  bool debug = false;
  if (debug) {
    *oLog << " -- GUtilityFuncts::linearInterpolation2D" << endl;
    *oLog << "       xval: " << xval << endl;
  }
  // see if xval is within range of orderec vector xvec.
  int size = xvec.size();
  if ( (xval < xvec[0]) || 
       (xval > xvec[size - 1]) ) {
    if (debug) {
      *oLog << "       xval is outside xvec range, " << xvec[0] << "  to  " 
           << xvec[size-1] << endl;
      *oLog << "       returning -9999.9 for interpolated value" << endl;
    }
    return -9999.9;

  }

  // find "insertion point" for xval in xvec - to get position interpolation 
  //  indices
  // using a binary search.
  // insertion point is the index below which xval would be inserted 
  // to maintain an ordered collection

  // set initial limits
  int low = 0;
  int hi = size-1;
      
  while (low <= hi) {
    int mid = (low + hi) /2;
    if (xval > xvec[mid] )
      low = mid + 1;
    else if ( xval < xvec[mid])
      hi = mid -1;
    else {
      if (debug) {
        *oLog << "       found exact match at index " << mid 
             << "  at xvec value " << xvec[mid] << endl;
        *oLog << "       yvec value here is " << yvec[mid] << endl;  
      }
    
      return yvec[mid];
    }
  }

  // interpolate between indices, low and hi
  hi = low;
  low--;
  
  if (debug) {
    *oLog << "       interpolate between indices " << low << " and " << hi;
    *oLog << "  for xvec range " << xvec[low] << "  to  " << xvec[hi] << endl;
    *oLog << "                                           and yvec range " << yvec[low] 
         << "  to  " << yvec[hi] << endl;
  }
  // do interpolation
  double interVal = yvec[low] + (yvec[hi]-yvec[low])*( xval - xvec[low])/(xvec[hi]-xvec[low]);
  
  if (debug) {
    *oLog << "       interpolated value is " << interVal << endl;
  }

  return interVal;
};
/**************** end of linearInterpolation2D ****************/

void GUtilityFuncts::addError(double v[3], const double &max_angle) {
  /*
    This procedure adds a random misallignment vector to a unit vector. The
    angle, max_angle, is the maximum conical angle between the unit vector 
    and the misaligned unit vector. 
  */

  double maxAngle = max_angle;
 
  double theta;     /* actual angle of misalignment */
  double phi;       /* azimuthal angle */
  double mat[3][3]; /* rotation matrix */
  double temp[3];
  double sv;

  int i;         /* for loop index */

  if ( (maxAngle <= 0.0) ) return;

  maxAngle /= TMath::RadToDeg();

  theta = maxAngle * sqrt( TR3.Rndm() );  
  phi = (TR3.Rndm()) * (TMath::TwoPi() )  ;
  
  if (((v[0]*v[0]) + (v[1]*v[1])) == 0) {
    v[0] = sin(theta) * cos(phi);
    v[1] = sin(theta) * sin(phi);
    v[2] = cos (theta);
  }

  else { 
    sv = sqrt( (v[0]*v[0]) + (v[1]*v[1]) );
    mat[0][0] = (v[0] * v[2]) / sv;
    mat[1][0] = (v[1] * v[2]) / sv;
    mat[2][0] = -sv;
    mat[0][1] = -v[1] / sv;
    mat[1][1] = v[0] / sv;
    mat[2][1] = 0;
    mat[0][2] = v[0];
    mat[1][2] = v[1];
    mat[2][2] = v[2];
    v[0] = sin(theta) * cos(phi);
    v[1] = sin(theta) * sin(phi);
    v[2] = cos (theta);
    mtxMlt(mat, v, temp);
    for (i = 0; i < 3; i++) v[i] = temp[i]; 
  }
};
/*********** end of addError *******************/

void GUtilityFuncts::mtxMlt(double a[3][3], double b[3], double c[3])
/* 
       | a[0][0] a[0][1] a[0][2] |   | b[0] |   | c[0] |
       | a[1][0] a[1][1] a[1][2] | x | b[1] | = | c[1] |
       | a[2][0] a[2][1] a[2][2] |   | b[2] |   | c[2] |                
*/
     /*.....................................................................*/
{
  int i, j;
  for (i = 0; i < 3; i++)
    {
      c[i] = 0.0;
      for(j = 0; j < 3; j++)
        c[i] += a[i][j] * b[j]; 
    }
};
/************** end of mtxMlt********************/
void GUtilityFuncts::addErrorRoot(ROOT::Math::XYZVector *vec, 
                                  const Double_t &maxAngle) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::addErrorRoot " << endl;
    *oLog << "        maxAngle Deg " 
          << maxAngle*(TMath::RadToDeg()) << endl;
    *oLog << "        vector orig  ";
    GUtilityFuncts::printGenVector(*vec); *oLog << endl;
  }

  Double_t r = maxAngle*sqrt(TR3.Rndm());
  Double_t theta   = ( TMath::TwoPi())* TR3.Rndm() - (TMath::Pi());

  ROOT::Math::XYZVector n1(sin(r)*cos(theta),
                           sin(r)*sin(theta),
                           cos(r));
  if (debug) {
    *oLog << "        error r theta (Deg) " << r*(TMath::RadToDeg()) 
          << "  " << theta*(TMath::RadToDeg()) << endl;
    *oLog << "        n1  ";
    GUtilityFuncts::printGenVector(n1); *oLog << endl;
  }

  Double_t phiR = vec->Phi();
  Double_t thetaR = vec->Theta();

  ROOT::Math::RotationZ rz(-((TMath::PiOver2()) -phiR));
  ROOT::Math::RotationX rx(-thetaR);
  ROOT::Math::Rotation3D rotM = rz*rx;
  
  //n1.SetCoordinates(0.0,0.0,1.0);
  *vec = rotM*n1;

  if (debug) {
    *oLog << "        phiR thetaR (Deg) " << phiR*(TMath::RadToDeg())
          << "  " << thetaR*(TMath::RadToDeg()) << endl;
    *oLog << "        new vector ";
    GUtilityFuncts::printGenVector(*vec); *oLog << endl;
    *oLog << endl;
  }

};
/************** end of addErrorRoot ********************/

void GUtilityFuncts::getGaussianErrorPolarCoor(const double &rMax,
                                                  const double &sigma,
                                                  double *r,double *theta) {
  bool debug = false;
  if (debug) {
    *oLog << "  -- getGaussianErrorPolarCoor " << endl;
    *oLog << " rMax,sigma " << rMax << " " << sigma << endl;
  }
  double tmp1 = 1 - exp(-rMax*rMax/(sigma*sigma));
  double tmp2 = 1 - (TR3.Rndm())* tmp1;
  if (debug) {
    *oLog << "tmp1 tmp2 " << tmp1 << " " << tmp2 << endl;
    *oLog << " logtmp2 " << log(tmp2) << endl;
  }
  *r = sqrt(-log(tmp2)) * sigma;
  *theta = (TMath::TwoPi()* TR3.Rndm()) - (TMath::Pi());

  if (debug) {
    *oLog << "r , theta " << *r << " " << *theta << endl;
  }
  return;
};
/************** end of getGaussianErrorPolarCoor ********************/

void GUtilityFuncts::addErrorGaussianRoot(ROOT::Math::XYZVector *vec, 
                                          const double &rMax,
                                          const double &sigma) {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::addErrorGaussianRoot" << endl;
    *oLog << "       rMax sigme " << rMax << "  " << sigma << endl;
    *oLog << "       vec orig ";
    GUtilityFuncts::printGenVector(*vec); *oLog << endl;
  }

  //  get Gaussian offset in polar coordinates
  double theta;
  double r;

  GUtilityFuncts::getGaussianErrorPolarCoor(rMax,sigma,&r,&theta);

  ROOT::Math::XYZVector n1(sin(r)*cos(theta),
                           sin(r)*sin(theta),
                           cos(r));

  // transform back
  Double_t phiR = vec->Phi();
  Double_t thetaR = vec->Theta();

  ROOT::Math::RotationZ rz(-((TMath::PiOver2()) -phiR));
  ROOT::Math::RotationX rx(-thetaR);
  ROOT::Math::Rotation3D rotM = rz*rx;
  
  //n1.SetCoordinates(0.0,0.0,1.0);
  *vec = rotM*n1;

  if (debug) {
    *oLog << "        phiR thetaR (Deg) " << phiR*(TMath::RadToDeg())
          << "  " << thetaR*(TMath::RadToDeg()) << endl;
    *oLog << "        new vector ";
    GUtilityFuncts::printGenVector(*vec); *oLog << endl;
    *oLog << endl;
  }
};
/************** end of addErrorGaussianRoot ********************/

double GUtilityFuncts::fieldRot(const double &zn, const double &az,
               const double &latitude) {

  /* returns angle on the celestial sphere between the arc from the (zn,az)
     defined point to the zenith and the arc from the (zn,az) defined
     point to Polaris.
     
     all angles in radians
  */
  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::fieldRot " << endl;
    *oLog << "      zn  : " << (TMath::RadToDeg())*zn << endl;
    *oLog << "      az  : " << (TMath::RadToDeg())*az << endl;
    *oLog << "      lat : " << (TMath::RadToDeg())*latitude << endl;
  }

  double el = (TMath::PiOver2()) - zn;

  if (debug) {
    *oLog << "      el  : " << (TMath::RadToDeg())*el << endl;
  }

  double sinR= sin(az)*cos(latitude);
  double cosR= -cos(az)*sin(el)*cos(latitude)+cos(el)*sin(latitude);
  double norm=sqrt(sinR*sinR+cosR*cosR);
  double phiR=0;
  if( norm > 0.0 ){
    cosR/=norm;
    sinR/=norm;
    phiR=acos(cosR);
    if(sinR<0) phiR=-phiR; 
  }
  if (debug) {
    *oLog << "      phiR  : " << (TMath::RadToDeg())*phiR << endl;
  }
  if (debug) {
    *oLog << "        get phir from atan2(sinR,cosR) , atan2(0.0) = 0 " 
          << endl;
    *oLog << "         sinR cosR " << sinR << "  " << cosR << endl;
    phiR = ( (sinR!=0.0)||(cosR!=0.0) ? atan2(sinR,cosR) : 0.0);
    *oLog << "    phiR from atan2  " << phiR*(TMath::RadToDeg()) << endl;
  }

  return phiR;
};
/************** end of getGaussianErrorPolarCoor ********************/

int GUtilityFuncts::AzZnToXYcos(const double &az,const double &zn,
                                   double *xcos, double *ycos) {
  /*For a given azimuth az and zenith zn, computes the x and y cosine of 
    the vector to the az,zn point on the unit sphere. 
    Coor.sys: x(East), y(North), z(Up). assumes zn > 0.0.
    all angles in radians

    Tests for zero within double numerical limits 
  */

  double epsilon = numeric_limits<double>::epsilon();
 
  *xcos=sin(zn)*sin(az);
  Bool_t testZero = TMath::AreEqualAbs(*xcos,0.0,epsilon);
  if (testZero) *xcos = 0.0;

  *ycos=sin(zn)*cos(az);
  testZero = TMath::AreEqualAbs(*ycos,0.0,epsilon);
  if (testZero) *ycos = 0.0;

  return 0;
};
/************** end of AzZnToXYcos ********************/

int GUtilityFuncts::XYcosToAzZn(const double &xcos, const double &ycos,
                                   double *az,double *zn) {

  /* determine az and zn angles for direction cosines, assuming zcos > 0.0 (zn > 0.0)
     Coordinate system, x(East), y(North), z (Up).
     az values: 0.0 to 369 degrees. All values here are in radians.
     az and zn set to zero if within 1.0E-12 of 0.0.
     all angles in radians
   */
  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::XYcosToAzZnNew " << endl;
    *oLog << "     xcos / ycos " << xcos << " " << ycos << endl;
  }

  double epsilon = numeric_limits<double>::epsilon();

  double r = xcos*xcos+ycos*ycos;
  Bool_t testZero = TMath::AreEqualAbs(r,0.0,epsilon); 
  if (!testZero) {
    *zn=asin(sqrt(r));
    *az=atan2(xcos,ycos);

    testZero = TMath::AreEqualAbs(*az,0.0,epsilon); 
    if (testZero) {
      *az = 0.0;
    }
    else if (*az < 0.0) {
      *az=*az+(TMath::TwoPi());
    }

    testZero = TMath::AreEqualAbs(*zn,0.0,epsilon); 
    if (testZero) {
      *zn = 0.0;
    }
  }
  else {
    *zn = 0.0;
    *az = 0.0;
  }
  if (debug) {
    *oLog << "    az  zn  " << (*az)*(TMath::RadToDeg()) << "  " 
          << (*zn)*(TMath::RadToDeg()) << endl;
  }

  return 0;
};
/************** end of XYcosToAzZn ********************/

void GUtilityFuncts::offsetXYToAzZn(const double &offsetX, const double &offsetY,
                    const double &az, const double &zn,
                    double *azOffset, double *znOffset) {

  /*  function to determine the azimuth and zenith of the offset point with coordinates 
      (offsetX, offsetY) with respect to tangent plane axes (parallel to telescope axes,
      where x is to the East and y is down.)
   */

  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::offsetXYToAzZn " << endl;
    *oLog << "       offsetX/Y " << offsetX*(TMath::RadToDeg()) << "  " 
          << offsetY*(TMath::RadToDeg()) << endl;
    *oLog << "       az/zn     " << az*(TMath::RadToDeg()) << "  " 
          << zn*(TMath::RadToDeg()) << endl;
  }

  /* The vector from the origin of the unit sphere to the az,zn point is (0,0,1)
     with respect to tangent plane coordinates.
     The vector from the origin of the tangent plane to the offset point is
     (offsetX, offsetY, 0). Thus, the sum of these vectors is the vector
     from the center of the unit sphere to the offset point on the tangent plane.

     Construct this vector, rotate it to ground coordinates, and find its az,zn point on 
     the unit sphere.
   */
  // make the vector
  double zp = 1.0;
  ROOT::Math::XYZVector vTang(offsetX,offsetY,zp);
  ROOT::Math::XYZVector vRot;
  
  // get the rotation matrix
  ROOT::Math::Rotation3D *rotMat = new ROOT::Math::Rotation3D;
  GUtilityFuncts::AzZnToRotMat(az,zn,rotMat);

  if (0) {
    double xcos = 0.0; double ycos = 0.0; double zcos = 0.0;
    GUtilityFuncts::AzZnToXYcos(az,zn,&xcos,&ycos);
    zcos = (1 - xcos*xcos - ycos*ycos);
    if (zcos > 0.0) {
      zcos = sqrt(zcos);
    }
    else {
      zcos = 0.0;
    }
    *oLog << endl << "     check rotation matrix:  " << endl;
    *oLog << "        vector to az/zn point " << xcos << " " << ycos
          << "  " << zcos << endl;
    ROOT::Math::XYZVector vPoint(xcos,ycos,zcos);
     ROOT::Math::XYZVector vRot1 = (*rotMat)*vPoint;
    *oLog << "        rotate vector, expect (0,0,1): to get " ;
    GUtilityFuncts::printGenVector(vRot1); *oLog << endl;
    *oLog << "        start with (0,0,1) and apply inverse rotation" << endl;
    *oLog << "        expect above vector to az/zn point ";
    ROOT::Math::XYZVector vz(0.0,0.0,1.0);
    vRot1 = (rotMat->Inverse())*vz;
    GUtilityFuncts::printGenVector(vRot1); *oLog << endl << endl;
  }

  // rotate vector to offset point from center of sphere to ground coordinates
  vRot = (rotMat->Inverse())*vTang;
  vRot = vRot.Unit();
  if (debug) {
    *oLog << "   vTang ";
    GUtilityFuncts::printGenVector(vTang); *oLog << endl;   
  }
  // find az,zn for this vector in ground coordinates
  GUtilityFuncts::XYcosToAzZn(vRot.X(),vRot.Y(),azOffset,znOffset);
  if (debug) {
    *oLog << " az/zn Offset " << (*azOffset)*(TMath::RadToDeg()) << "  " 
          << (*znOffset)*(TMath::RadToDeg()) << endl;
    *oLog << " ------------- end of GUtilityFuncts::offsetXYToAzZn ------ " << endl;
  }
};
/************** end of offsetXYToAzZn  ********************/

void GUtilityFuncts::telescopeAzZnNew(const double &primAz,
                                   const double &primZn,
                                   const double &wobbleN,
                                   const double &wobbleE,
                                   const double &telOffsetX,
                                   const double &telOffsetY,
                                   const double &latitude,
                                   double *telAz,double *telZn,
                                   double *sourceX,double *sourceY) {

  /*  function to determine the az,zn of the telscope. The telescope is moved from 
      the primAz,primZn location through the vector (wobbleN,wobbleE) on the tangent 
      plane after adding the telescope offset vector.  The Source vector is the location
      of the source (primAz,primZn) location on the tangent plane centered at the new 
      location of the telescope.  The usual coordinate systems apply: x(East), Y(North),
      and z(Up) for the ground system.  The tangent plane coordinates are y(down from the 
      zenith), x (perpendicular to y in the direction of increasing azimuth).
   */
  bool debug = false;
  
  *sourceX = 0.0;
  *sourceY = 0.0;

  double epsilon = numeric_limits<double>::epsilon();

  if (debug) {
    *oLog << "------------------------------------------------" << endl;
    *oLog << "  -- GUtilityFuncts::telescopeAzZnNew"  << endl;
    *oLog << "         prim Az / Zn   " << primAz*(TMath::RadToDeg())
          << "  /  " << primZn*(TMath::RadToDeg()) << endl;
    *oLog << "         wobbleN / E    " << wobbleN*(TMath::RadToDeg())
          << "  /  " << wobbleE*(TMath::RadToDeg()) << endl;
    *oLog << "         latitude       " << latitude*(TMath::RadToDeg()) << endl;
    *oLog << "         telOffSetX / Y " << telOffsetX*(TMath::RadToDeg())
          << "  /  " << telOffsetY*(TMath::RadToDeg()) << endl;
    *oLog << "--------------------------------------------------" << endl;
  }
  bool haveWobbleFlag =  ! ((TMath::AreEqualAbs(wobbleN,0.0,epsilon)) && (TMath::AreEqualAbs(wobbleE,0.0,epsilon) ) );
  bool haveOffsetFlag =  ! ((TMath::AreEqualAbs(telOffsetX,0.0,epsilon)) && (TMath::AreEqualAbs(telOffsetY,0.0,epsilon) ) );

  *telAz = primAz;
  *telZn = primZn;

  if ( (haveWobbleFlag || haveOffsetFlag) ) {
    if (debug) {
      *oLog << " ready to offset/wobble " << endl;
      *oLog << "------------------------------" << endl;
    }
  

    /* find telescope az/zn after telescope offsets. The offsetx/y is the 
       arc length in radians along a great circle arc starting at the origin
       of the tangent plane (at the location of the primary on the unit 
       sphere.) If "a" is the arc length, then the equivalent distance, "b",
       along the tangent plane is "b = tan(a)". The offsets are along the
       telescope X and Y axes where a positive y offset increases the 
       telescope zenith angle and a positive x offset increases the telescope 
       azimuth and will slightly increase the zenith angle since the x 
       offset is along a great circle (all straight lines on the tangent
       plane are along a great circle on the unit sphere).
    */
    if (haveOffsetFlag) {
      // get offset distances on tangent plane (Tp)
      double telOffsetXTp = tan(telOffsetX);
      double telOffsetYTp = tan(telOffsetY);
      
      if (debug) {
	*oLog << "get offset distances on the tangent plane from tan(offsetx/y) " << endl;
	DEBUG(telOffsetX); DEBUG(telOffsetY);
	DEBUG(telOffsetXTp); DEBUG(telOffsetYTp);
	DEBUG(telOffsetX*(TMath::RadToDeg())); DEBUG(telOffsetY*(TMath::RadToDeg()));
	DEBUG(telOffsetXTp*(TMath::RadToDeg())); DEBUG(telOffsetYTp*(TMath::RadToDeg()));
      }
      // get new unit vector locating the telescope in telescope coordinates
      // by adding (0,0,1) to the x and y coordinates and finding the unit vector.
      // Use the unit vector to find the az/zn for the telescope.
      //        first the rotation matrix. 
      ROOT::Math::RotationZ rz1(primAz);
      ROOT::Math::RotationX rx1(primZn);
      ROOT::Math::Rotation3D rotM1 = rx1*rz1;
      
      //        then the direction cosines toward the unit sphere for the primary
      double xcos0,ycos0;
      GUtilityFuncts::AzZnToXYcos(primAz,primZn,&xcos0,&ycos0);
      
      if (debug) {
	*oLog << "---------" << endl;
	*oLog << "direction cosines for the primary" << endl;
	double zcos0 = sqrt(1 - xcos0*xcos0 - ycos0*ycos0);
	DEBUG(xcos0); DEBUG(ycos0); DEBUG(zcos0);
      }
      // 
      ROOT::Math::XYZVector vTel1(telOffsetXTp, telOffsetYTp,1.0);
      ROOT::Math::XYZVector vTel1Gr = rotM1.Inverse()*vTel1.Unit();
      double az1,zn1;
      GUtilityFuncts::XYcosToAzZn(vTel1Gr.X(),vTel1Gr.Y(),&az1,&zn1);
      
      if (debug) {
	*oLog << "----------" << endl;
	GUtilityFuncts::printGenVector(vTel1); *oLog << " vTel1; vector to offset point " << endl;
	*oLog << "            on tangent plane in tangent plane coor." << endl;
	GUtilityFuncts::printGenVector(vTel1Gr); *oLog << "  vTel1Gr Unit vector in ground coordinates" << endl;
	*oLog << "----------" << endl;
	*oLog << "new Az/Zn after teloffsets " << az1*(TMath::RadToDeg())
	      << "  " << zn1*(TMath::RadToDeg()) << endl;
      }
      // update telAz and telZn, no wobble offsets yet.
      *telAz = az1;
      *telZn = zn1;
    }
    
    if (haveWobbleFlag) {
      // get the rotation matrix from ground to telescope coordinates
      double wobbleNTp = tan(wobbleN);
      double wobbleETp = tan(wobbleE);

      if (debug) {
	*oLog << "----------------- starting wobble section " << endl;
	DEBUG(wobbleN); DEBUG(wobbleE);
	DEBUG(wobbleN*(TMath::RadToDeg())); DEBUG(wobbleE*(TMath::RadToDeg()));
	DEBUG(wobbleNTp); DEBUG(wobbleETp);
	DEBUG(wobbleNTp*(TMath::RadToDeg())); DEBUG(wobbleETp*(TMath::RadToDeg()));
      }
      
      ROOT::Math::RotationZ rz2(*telAz);
      ROOT::Math::RotationX rx2(*telZn);
      ROOT::Math::Rotation3D rotM2 = rx2*rz2;
      
      double wobbleX = 0.0, wobbleY = 0.0, fRot = 0.0;
      
      fRot = GUtilityFuncts::fieldRot(*telZn,*telAz,latitude);
      wobbleY = -wobbleNTp*cos(fRot) - wobbleETp*sin(fRot);
      wobbleX = -wobbleNTp*sin(fRot) + wobbleETp*cos(fRot);
      
      ROOT::Math::XYZVector vTel2(wobbleX, wobbleY,1.0);
      ROOT::Math::XYZVector vTel2Gr = rotM2.Inverse()*vTel2.Unit();
      double az2,zn2;
      double xcos2 = vTel2Gr.X();
      double ycos2 = vTel2Gr.Y();
      GUtilityFuncts::XYcosToAzZn(xcos2,ycos2,&az2,&zn2);
      
      *telAz = az2;
      *telZn = zn2;
      
      if (debug) {
	*oLog << "---------------" << endl;
	GUtilityFuncts::printGenVector(vTel2); *oLog << " vTel2 to wobble offset point on tangent plane" << endl;
	*oLog << "                        in tangent plane coordinates" << endl;
	*oLog << "                        includes prior telOffsets " << endl;
	GUtilityFuncts::printGenVector(vTel2Gr); *oLog << "  vTel2Gr Unit vector in ground coordinates" << endl;
	*oLog << "----------" << endl;
	*oLog << "new Az/Zn after teloffsets and wobble offsets " << az2*(TMath::RadToDeg())
	      << "  " << zn2*(TMath::RadToDeg()) << endl;
      }
    }
  }
  if (debug) {
    *oLog << " ======= final az/zn values ========= " << endl;
    DEBUG((*telAz)*(TMath::RadToDeg())); DEBUG((*telZn)*(TMath::RadToDeg()));
    
    // get source location relative to telescope.
    double xs,ys;
    GUtilityFuncts::sourceOnTelescopePlane(primAz, primZn,
					   *telAz, *telZn,
					   &xs, &ys);    
    *sourceX = xs;
    *sourceY = ys;
  }
};
/**************** end of telescopeAzZnNew ************************/

void GUtilityFuncts::telescopeAzZn(const double &primAz,
                                   const double &primZn,
                                   const double &wobbleN,
                                   const double &wobbleE,
                                   const double &telOffsetX,
                                   const double &telOffsetY,
                                   const double &latitude,
                                   double *telAz,double *telZn,
                                   double *sourceX,double *sourceY) {

  /*  function to determine the az,zn of the telscope. The telescope is moved from 
      the primAz,primZn location through the vector (wobbleN,wobbleE) on the tangent 
      plane after adding the telescope offset vector.  The Source vector is the location
      of the source (primAz,primZn) location on the tangent plane centered at the new 
      location of the telescope.  The usual coordinate systems apply: x(East), Y(North),
      and z(Up) for the ground system.  The tangent plane coordinates are y(down from the 
      zenith), x (perpendicular to y in the direction of increasing azimuth).
   */

  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::telescopeAzZn"  << endl;
    *oLog << "     prim Az / Zn   " << primAz*(TMath::RadToDeg())
          << "  /  " << primZn*(TMath::RadToDeg()) << endl;
    *oLog << "     wobbleN / E    " << wobbleN*(TMath::RadToDeg())
          << "  /  " << wobbleE*(TMath::RadToDeg()) << endl;
    *oLog << "     latitude       " << latitude*(TMath::RadToDeg()) << endl;
    *oLog << "     telOffSetX / Y " << telOffsetX*(TMath::RadToDeg())
          << "  /  " << telOffsetY*(TMath::RadToDeg()) << endl;
  }

  if ( (wobbleN==0.0)&&(wobbleE==0.0)&&
       (telOffsetX==0.0)&&(telOffsetY==0.0) ) {
    *telAz = primAz;
    *telZn = primZn;
  }
  else {
    // find total offset for telescope in the tangent plane.
    double wobbleX = 0.0;
    double wobbleY = 0.0;
    double fRot = 0.0;

    double xOffsetAll = 0.0;
    double yOffsetAll = 0.0;

    if ( (wobbleN !=0.0)|| (wobbleE != 0.0) ) {
      fRot = GUtilityFuncts::fieldRot(primZn,primAz,latitude);
      wobbleY = -wobbleN*cos(fRot) - wobbleE*sin(fRot);
      wobbleX = -wobbleN*sin(fRot) + wobbleE*cos(fRot);
    }

    xOffsetAll = wobbleX + telOffsetX;
    yOffsetAll = wobbleY + telOffsetY;
 
    // the Source vector is in telescope, not camera coordinates.
    *sourceX = -xOffsetAll;
    *sourceY = -yOffsetAll;

    GUtilityFuncts::offsetXYToAzZn(xOffsetAll, yOffsetAll,primAz,primZn,
                                   telAz,telZn);             
 
    if (debug) {
      *oLog << "  fRot " << fRot * (TMath::RadToDeg()) << endl;
      *oLog << "  xOffsetAll  yOffsetAll  primAz  primZn "  
            << xOffsetAll*(TMath::RadToDeg()) << "  " 
            << yOffsetAll*(TMath::RadToDeg()) << "  " 
            << primAz*(TMath::RadToDeg()) << "  " 
            << primZn*(TMath::RadToDeg()) << endl;
      *oLog << "       telAz  telZn    sourceX  sourceY "  
            << *telAz*(TMath::RadToDeg()) << "  "
            << *telZn*(TMath::RadToDeg()) << "  "
            << *sourceX*(TMath::RadToDeg()) << "  "
            << *sourceY*(TMath::RadToDeg()) << endl;
    }
  }
};
/***************** end of telescopeAzZn  ***************/

void GUtilityFuncts::AzZnToRotMat(const double &az,
                                  const double &zn,
                                  ROOT::Math::Rotation3D *rotM) {
  bool debug = false;
  ROOT::Math::RotationZ rz(az);
  ROOT::Math::RotationX rx(zn);
  *rotM = rx*rz;

  if (debug) {
    *oLog << "    in GUtilityFuncts::AzZnToRotMat" << endl;
    *oLog << "      az zn " << az*TMath::RadToDeg() << "   "
	  << zn*TMath::RadToDeg() << endl;
  }
};
/*****************  end of rotXYZcos ***************/
bool GUtilityFuncts::testZeroFloat(double *ax, const int &zeroFlag) {

  bool testEq = TMath::AreEqualAbs(*ax,0.0,numeric_limits<double>::epsilon());
  if ( (testEq ==true) && (zeroFlag ==1) ) {
      *ax = 0.0;
    }
  return testEq;
};
/*****************  end of testZeroFloat ***************/

bool GUtilityFuncts::zeroFloatVectorFix(ROOT::Math::XYZVector *vec) {
  
  bool test = false;
  
  if (TMath::AreEqualAbs(vec->X(),0.0,numeric_limits<double>::epsilon())) {
    vec->SetX(0.0);
    test = true;
  }
  if (TMath::AreEqualAbs(vec->Y(),0.0,numeric_limits<double>::epsilon())) {
    vec->SetY(0.0);
    test = true;
  }
  if (TMath::AreEqualAbs(vec->Z(),0.0,numeric_limits<double>::epsilon())) {
    vec->SetZ(0.0);
    test = true;
  }
  return test;

};
/*****************  end of zeroFloatVectorFix ***************/

bool GUtilityFuncts::photonReflect(const vector<double> &refCoeff,
                                   const vector<double> &refWaveLgt,
                                   const double &reflect,
                                   const double &photonWaveLgt ) {

  bool debug = false;

  bool reflectFlag = false;

  int numberRC = refCoeff.size();
  double minW = refWaveLgt[0];
  double maxW = refWaveLgt[numberRC - 1];
  double refl = 0.0;  // interpolated reflectivity from table

  if (debug) {
    *oLog << "  -- GUtilityFuncts::photonReflect" << endl;
    *oLog << "        reflect " << reflect << endl;
    *oLog << "        photonWaveLgt " << photonWaveLgt << endl;
    *oLog << "        refCoeff.size() " << refCoeff.size() << endl;
    *oLog << "        minW maxW  " << minW << "  " << maxW << endl;
  }

  // is the wavelength within range of the table
  if ( (photonWaveLgt >= minW) &&
       (photonWaveLgt <= maxW) ) {

    // loop on reflectivity curve points, no binary search for now since
    // points are not necessarily equidistant

    for (int n = 1;n<numberRC;n++) {
      double loWl = refWaveLgt[n-1];
      double hiWl = refWaveLgt[n];

     if ( (loWl <= photonWaveLgt) && (photonWaveLgt <=hiWl) ) {
        double loR = refCoeff[n-1];
        double hiR = refCoeff[n];

        double m = (hiR - loR)/(hiWl - loWl);   // slope within local limits
        double b = loR - (m*loWl);

        refl = ((m*photonWaveLgt) + b)*reflect;
	double trand = TR3.Rndm();
        if (trand < refl ) {
          reflectFlag = true;
	}
	if (debug) {
	  *oLog << " reflection test refl trand " << refl
		<< " " << trand << endl;
	}
        break;
      }
    }
  }
  if (debug) {
    *oLog << "       reflection completed, if refl = 0.0.  no reflection" 
          << endl;
    *oLog << "       refl = " << refl << endl;
    *oLog << "       reflectFlag " << reflectFlag << endl;
  }

  return reflectFlag;
};
/******************* end of photonReflect ***************/ 

void GUtilityFuncts::reflectDirection(const ROOT::Math::XYZVector &vnormUnit,
                      const ROOT::Math::XYZVector &vphotonUnit,
                      ROOT::Math::XYZVector *vphotonReflDcos,
                      const double &roughness, const int &option) {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::reflectDirection" << endl;
    *oLog << "       vnormUnit ";
    GUtilityFuncts::printGenVector(vnormUnit); *oLog << endl;
    *oLog << "       vphotonUnit ";
    GUtilityFuncts::printGenVector(vphotonUnit); *oLog << endl;
    *oLog << "       option " << option << endl;
  }

  double dot = vphotonUnit.Dot(vnormUnit);
  *vphotonReflDcos = vphotonUnit - 2*dot*vnormUnit;

  if (debug) {
    *oLog << "        dot " << dot << endl;
    *oLog << "        vphotonReflDcos no roughness ";
    GUtilityFuncts::printGenVector(*vphotonReflDcos); *oLog << endl;
  }
  GUtilityFuncts::addErrorRoot(vphotonReflDcos,roughness);
  
  if (debug) {
    *oLog << "        vphotonReflDcos with roughness ";
    GUtilityFuncts::printGenVector(*vphotonReflDcos); *oLog << endl;
  }
};
/******************* end of reflectDirection ***************/

bool GUtilityFuncts::sortPair(const pair<int,int> i , 
                              const pair<int,int> j)  {
  return (j.second < i.second);
};
/******************** end of sortPair *********************************/

void GUtilityFuncts::printSegVector (const vector<mirrorSegmentDetails *> &vec) {
  *oLog << "        num rmin    rmax marg dPhi refl rough errXYZ rotErrABC"
        << endl;
  int numElem = vec.size();
  for (int i = 0;i<numElem; i++) {
    mirrorSegmentDetails *t = vec[i];
    *oLog << "         " << i+1 << "  " 
          << t->rmin << "  " << t->rmax << "  " << t->margin << "  " 
          << t->delPhi << "  " << t->reflect << "   " << t->roughness
          << "      " << t->posErrorX
          << " " << t->posErrorY << " " << t->posErrorZ << "    " 
          << t->rotErrorPhi << " " << t->rotErrorTheta << " " << t->rotErrorPsi
          << endl;
  }
};
/***************************** end of printSegVector **********************/

void GUtilityFuncts::printXYZVector(const ROOT::Math::XYZVector &vec,const string &label) {
  *oLog << label << "   " << vec.X() << "  " << vec.Y() << " "
        << vec.Z() << endl;

}
/*************** end of printXYZVector *********************/
void GUtilityFuncts::sourceOnTelescopePlane(const double &az_s, const double &zn_s,
					    const double &az_t, const double &zn_t,
					    double *x, double *y) {
  *x = 0.0;
  *y = 0.0;
  
  bool debug = false;
  if (debug) {
    *oLog << "---------- in GUtilityFuncts::sourceOnTelescopePlane -------" << endl;
    *oLog << "    az_s  zn_x  " << az_s*TMath::RadToDeg() << "  " << zn_s*TMath::RadToDeg()
	  << endl;
    *oLog << "    az_t  zn_t  " << az_t*TMath::RadToDeg() << "  " << zn_t*TMath::RadToDeg()
	  << endl;
  }

  // find the vector in the telescope tangent plane from the plane center to the intersection of the
  // vector locating the source vector with the tangent plane.  Do this using vector algebra rather
  // than explicit spherical trig.

  // find the unit vectors locating the telescope and the source on the unit sphere.
  double xcos_t = 0.0;
  double ycos_t = 0.0;
  double zcos_t = 0.0;
  GUtilityFuncts::AzZnToXYcos(az_t,zn_t,&xcos_t, &ycos_t);
  zcos_t = sqrt(1 - xcos_t*xcos_t - ycos_t*ycos_t);
  ROOT::Math::XYZVector nUnit_t(xcos_t,ycos_t,zcos_t);
  
  // and now the rotation matrix to get to telescope coordinates
  ROOT::Math::Rotation3D rotM;
  GUtilityFuncts::AzZnToRotMat(az_t,zn_t,&rotM);
  
  double xcos_s = 0.0;
  double ycos_s = 0.0;
  double zcos_s = 0.0;
  GUtilityFuncts::AzZnToXYcos(az_s,zn_s,&xcos_s, &ycos_s);
  zcos_s = sqrt(1 - xcos_s*xcos_s - ycos_s*ycos_s);
  ROOT::Math::XYZVector nUnit_s(xcos_s,ycos_s,zcos_s);

  // some vector algebra to find the vector to the intersection point
  double dotP = nUnit_s.Dot(nUnit_t);
  ROOT::Math::XYZVector tTos_vecGC = (nUnit_s/dotP) - nUnit_t; 
  ROOT::Math::XYZVector tTos_vecTC = rotM*tTos_vecGC;

  // check for zero components, set to zero is within numeric limits
  GUtilityFuncts::zeroFloatVectorFix(&tTos_vecTC);
  // set return x and y components, z component should always be zero.
  *x = tTos_vecTC.X();
  *y = tTos_vecTC.Y();
  
  if (debug) {
    *oLog << "  unit vector locating source " << xcos_s << "  " << ycos_s << "  "
	  << zcos_s << endl;
    *oLog << "  unit vector locating telescope " << xcos_t << "  " << ycos_t << "  "
	  << zcos_t << endl;
    *oLog << "  dot product of tel and source unit vecs " << dotP << endl;
    *oLog << "  nUnit_s "; GUtilityFuncts::printGenVector(nUnit_s); *oLog << endl;
    *oLog << "  nUnit_t "; GUtilityFuncts::printGenVector(nUnit_t); *oLog << endl;
    *oLog << "  tTos_vecGC "; GUtilityFuncts::printGenVector(tTos_vecGC); *oLog << endl;
    *oLog << "  tTos_vecTC "; GUtilityFuncts::printGenVector(tTos_vecTC); *oLog << endl;
    *oLog << "---------- exiting GUtilityFuncts::sourceOnTelescopePlane -------" << endl;
  }
};
/**************************** end of sourceOnTelescopePlane ********************/

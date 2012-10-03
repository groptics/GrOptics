/*
VERSION2.3
3OCT2012
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
  double PI =  3.141592654;     /*   Value of pi */

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
}

/********************  end of polyInside ***************/

//void GUtilityFuncts::tokenizer(const string& str, vector<string>& tokens) {
  
//string buf; // Have a buffer string
//stringstream ss(str); // Insert the string into a stream
  
//while (ss >> buf)
//  tokens.push_back(buf);
//}
/********************  end of tokenizer ************/

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


}
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
//**************** end of linearInterpolation2D ****************

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
    *oLog << "        get phir from atan2(sinR,cosR) , zero = both are 0" 
          << endl;
    *oLog << "         sinR cosR " << sinR << "  " << cosR << endl;
    phiR = ( (sinR!=0.0)||(cosR!=0.0) ? atan2(sinR,cosR) : 0.0);
    *oLog << "      " << phiR*(TMath::RadToDeg()) << endl;
  }

  return phiR;
};
/************** end of getGaussianErrorPolarCoor ********************/

int GUtilityFuncts::AzZnToXYcos(const double &az,const double &zn,
                                   double *xcos, double *ycos) {
    //For a given azimuth f and elevation t, computes the x and y cosine of 
    //the vector directing the corresponding direction
  double t;
  t = (TMath::PiOver2()) - zn;

    *xcos=cos(t)*sin(az);
    *ycos=cos(t)*cos(az);
return 0;
};
/************** end of AzZnToXYcos ********************/

int GUtilityFuncts::XYcosToAzZn(const double &xcos, const double &ycos,
                                   double *az,double *zn) {
  /*
    NEED TO FIX THIS LATER, NEED ZCOS, SINCE SIGN IS IMPORTANT
    For a given azimuth az and zenith angle zn, computes the 
    x and y direction cosines pointing to the sky location
  */

  bool debug = false;  // tested for 0 to pi/2 azimuth
  if (debug) {
    *oLog << "  -- GUtilityFuncts::XYcosToAzZn " << endl;
    *oLog << "     xcos / ycos " << xcos << " " << ycos << endl;
  }

  double r = xcos*xcos+ycos*ycos;
  double t;

  if (r > 0.0) {
    t=acos(sqrt(r));
    *az=atan2(xcos,ycos);
 
    if(*az<0.0) {
      if (*az > -1.0E-10) {
        *az = 0.0;
      }
      else {
      *az=*az+(TMath::TwoPi());
      }
    }
  }
  else {
    t = (TMath::PiOver2());   // elevation is 90 degrees
    *az = 0.0;      // by definition azimuth is 0.0 degrees
  }  

  *zn = (TMath::PiOver2()) - t;
  if (debug) {
    *oLog << "      az / zn " << (TMath::RadToDeg())*(*az) 
          << "  " << (TMath::RadToDeg())*(*zn) << endl;
  }
  return 0;
};
/************** end of XYcosToAzZn ********************/

void GUtilityFuncts::wobbleToAzZn(const double &wobbleN, 
                                     const double &wobbleE,
                                     const double &latitude,
                                     const double &aZ, 
                                     const double &zN,
                                     double *aZNew, 
                                     double *zNNew) {

  // in the rotated coordinate system (from aZ,zN), the vector
  // to the offset point in the tangent plane is (


  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::wobbleToAzZn " << endl;
    *oLog << "     wobbleN / E " << wobbleN*(TMath::RadToDeg())
          << " " << wobbleE*(TMath::RadToDeg()) 
          << endl;
    *oLog << "     az zn       " << aZ*(TMath::RadToDeg())
          << " " << zN*(TMath::RadToDeg()) << endl;
  }

  double rot = GUtilityFuncts::fieldRot(zN,aZ,latitude);

  if (debug) {
    *oLog << "     field rotation " << rot*(TMath::RadToDeg())
          << endl;
  }

  double delZn = -wobbleN*cos(rot)+wobbleE*sin(rot);
  double delTheta = -wobbleN*sin(rot) - wobbleE*cos(rot);
  if (debug) {
    *oLog << "     delZn delTheta  " << delZn*(TMath::RadToDeg())
          << " " << delTheta*(TMath::RadToDeg())
          << endl;
  }

  *zNNew = zN + delZn;

  if (zN > 0.00001) {
    *aZNew = aZ + (delTheta/sin(zN));
  }
  else {
    *aZNew = aZ;
  }

  if (debug) {
    *oLog << "     zNNew aZNew " << (*zNNew)*(TMath::RadToDeg())
          << " " << (*aZNew)*(TMath::RadToDeg()) << endl; 
  }

};
/************** end of wobbleToAzZn ********************/

void GUtilityFuncts::telescopeAzZnRot(const double &primAz,
                                      const double &primZn,
                                      const double &wobbleN,
                                      const double &wobbleE,
                                      const double &telOffsetX,
                                      const double &telOffsetY,
                                      const double &latitude,
                                      double *telAz,double *telZn) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GUtilityFuncts::telescopeAzZnRot"  << endl;
    *oLog << "     prim Az / Zn   " << primAz*(TMath::RadToDeg())
          << "  /  " << primZn*(TMath::RadToDeg()) << endl;
    *oLog << "     wobbleN / E    " << wobbleN*(TMath::RadToDeg())
          << "  /  " << wobbleE*(TMath::RadToDeg()) << endl;
    *oLog << "     telOffSetX / Y " << telOffsetX*(TMath::RadToDeg())
          << "  /  " << telOffsetY*(TMath::RadToDeg()) << endl;
  }

  ROOT::Math::Rotation3D *rotPrim = new ROOT::Math::Rotation3D;
  GUtilityFuncts::AzZnToRotMat(primAz,primZn,rotPrim);

  if ( (wobbleN==0.0)&&(wobbleE==0.0)&&
       (telOffsetX==0.0)&&(telOffsetY==0.0) ) {
    *telAz = primAz;
    *telZn = primZn;
  }
  else {
    // find and add components in the tangent plane to the primary 
    // unit vector that locates the center of the tangent plane, all
    // The primary location vector is (0,0,1) in this frame.

    double wobbleX = 0.0;
    double wobbleY = 0.0;
    double fRot;

    double xOffsetAll = 0.0;
    double yOffsetAll = 0.0;
    if ( (wobbleN !=0.0)|| (wobbleE != 0.0) ) {
        fRot = GUtilityFuncts::fieldRot(primZn,primAz,latitude);
        wobbleY = -wobbleN*cos(fRot) + wobbleE*sin(fRot);
        wobbleX = - wobbleN*sin(fRot) - wobbleE*cos(fRot);
        xOffsetAll = wobbleX + telOffsetX;
        yOffsetAll = wobbleY + telOffsetY;
        
        if (debug) {
          *oLog << "     fRot           " << fRot*(TMath::RadToDeg()) << endl;
          *oLog << "     wobbleX / Y    " << wobbleX*(TMath::RadToDeg())
                << "  /  " << wobbleY*(TMath::RadToDeg()) << endl;
          *oLog << "     xOffsetAll / y " << xOffsetAll*(TMath::RadToDeg())
                << "  /  " << yOffsetAll*(TMath::RadToDeg()) << endl;
        }
        
        // vector to tangent point in prim. coor.system, 
        double xp,yp,zp;
        xp = xOffsetAll;
        yp = yOffsetAll;
        zp = 1.0;

        // rotate coor. system back to ground system
        ROOT::Math::XYZVector vTang(xp,yp,zp);
        ROOT::Math::XYZVector vGrd = (rotPrim->Inverse())*vTang;

        if (debug) {
          *oLog << "      vTang  ";
          GUtilityFuncts::printGenVector(vTang); *oLog << endl; 

          *oLog << "      vGrd   ";
          GUtilityFuncts::printGenVector(vGrd); *oLog << endl;
        }
        // normalize vGrd
        // This unit vector points at the telescope location on the sky.
        vGrd = vGrd.Unit();
       
        
        // get aZ and zN for telescope
        XYcosToAzZn(vGrd.X(), vGrd.Y(), telAz, telZn);

        if (debug) {
          *oLog << "     vGrdUnit   ";
          GUtilityFuncts::printGenVector(vGrd); *oLog << endl;
          *oLog << "     telAz / Zn " << *telAz << " " << *telZn << endl;
        }
        
        if (debug) {
          // print primary and tel Az,Zn
          *oLog << "    primary Az Zn  " << primAz*(TMath::RadToDeg())
                << "  " << primZn*(TMath::RadToDeg()) << endl;
          *oLog << "    teles   Az Zn  " << (*telAz)*(TMath::RadToDeg())
                << "  " << (*telZn)*(TMath::RadToDeg()) 
                << endl;
          
        }
    }
  }

  SafeDelete(rotPrim);
  
};


void GUtilityFuncts::XYZcosToRotMat(const double &xcos,
                               const double &ycos,
                               const double &zcos,
                               ROOT::Math::Rotation3D *rotM) {

  double az,zn;
  GUtilityFuncts::XYcosToAzZn(xcos,ycos,&az,&zn);

};
/*****************  end of rotXYZcos ***************/
void GUtilityFuncts::AzZnToRotMat(const double &az,
                                  const double &zn,
                                  ROOT::Math::Rotation3D *rotM) {
  
  ROOT::Math::RotationZ rz(az);
  ROOT::Math::RotationX rx(zn);
  *rotM = rx*rz;

};
/*****************  end of rotXYZcos ***************/
bool GUtilityFuncts::testZeroFloat(const double &ax) {
    return TMath::AreEqualAbs(ax,0.0,numeric_limits<float>::epsilon());

};
/*****************  end of testZeroFloat ***************/

bool GUtilityFuncts::zeroFloatVectorFix(ROOT::Math::XYZVector *vec) {
  
  bool test = false;
  if (GUtilityFuncts::testZeroFloat(vec->X()) ) {
    vec->SetX(0.0);
    test = true;
  }
  if (GUtilityFuncts::testZeroFloat(vec->Y()) ) {
    vec->SetY(0.0);
    test = true;
  }
  if (GUtilityFuncts::testZeroFloat(vec->Z()) ) {  
    vec->SetZ(0.0);
    test = true;
  }
  return true;

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

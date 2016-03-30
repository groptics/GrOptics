/*!  makePhotonsPtSrc.cpp
     version 2.0
     28March2016
     Charlie Duke
     Grinnell College

    makePhotonsOffset.cpp produces a grisudet or grOptics compatible .cph file 
   of cherenkov photons offset from a specified primary direction. 
   The photons are randomly placed on the telescope area with  
   radius increased by 10%.  Telescope locations and radii are from 
   a standard GrISU telescope configuration file. The code has 
   a ROOT dependency. All input information designated by this pilot file.

   All log information goes to cerr.

   Testing of this code is straightforward using a photonHistory.root debug
   file from grOptics. Or, when using grisudet, create a photonhistory file
   (see detector.pilot) and use the utility code to move the history data
   into a root file.

   usage:  makePhotonsOffset [name of pilotfile: default makePhotonsOffset.pilot

   Coor.System: ground X(East), Y(North), Z(up). When pointing to zenith,
   the telescope axes are parallel to these
   ground coordinate axes.  The camera axes are slightly different, camera
   X axis is parallel to the telescope x axis; the camera y axis is
   opposite to the telescope y axis (thus permitting x camera to the right and y
   camera up when facing the camera with telescope in stow position).


   usage:  makePhotonsPtSrc [name of pilotfile: default makePhotonsPtSrc.pilot]
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
//#include <deque>
//#include <list>
//#include <iterator>
//#include <algorithm>
//#include <bitset>
//#include <iomanip>

using namespace std;

#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TRandom3.h"

#include "VG_Pilot.h"
#include "GUtilityFuncts.h"

ostream *oLog;
TRandom3 TR3;

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

struct Pilot {
  string configFileName;
  string outFileName;  
  int waveLgt;         // int I think
  double phoTime;       // int or float
  int maxShowers;
  int maxPhotons;
  vector<ROOT::Math::XYZVector *> vTelLoc;
  vector<double> TelRadius;
 
  UInt_t seedr;

  double obser; // observatory height, default 1277.06

  double az; // az and zn of primary
  double zn;
  double wobbleN; // telescope offset from primary direction
  double wobbleE;
  double wobbleR;  // radius centered on wobbleN and wobbleE
  bool debug;      // when true, print intermediate results
};

void readPilot(struct Pilot *pilot,const string &pilotfile);
void readConfig(struct Pilot *pilot);
void printXYZVector(const ROOT::Math::XYZVector &vec,const string &label);

int main(int argc, char *argv[]) {


  
  oLog = &cerr;
  *oLog << "  -- start of makePhotonsPtSrc" << endl;

  string pilotfile("makePhotonsOffset.pilot");
  if (argc > 2) {
    *oLog << " TOO MANY COMMAND LINE ARGUMENTS  " << endl;
    *oLog << " usage: make_photons [name of pilot file]" << endl;
    *oLog << "     default pilot file name: makePhotons.pilot" << endl;
    exit(EXIT_FAILURE);
  }  
    else if ( argc == 2 ) {
      pilotfile = argv[1];
    }
  *oLog << "    pilotfile: " << pilotfile << endl;

  Pilot pilot;
  // read pilot and configuration files
  readPilot(&pilot,pilotfile);
  readConfig(&pilot);

  bool debug = pilot.debug;
  
  // write photon density to log file
  *oLog << "tel_num    telRadius*1.1     area       photon_density" << endl;

  for (unsigned tel=0;tel < pilot.TelRadius.size();tel++) {
    double radius = pilot.TelRadius[tel];
    double farea = (TMath::Pi())*radius*radius;
    double fdensity = (float)pilot.maxPhotons/farea;
    *oLog << "  " << tel+1 << "            " <<  radius
          << "          " << farea << "        " << fdensity << endl; 
  }

  // set up random number generator
  TR3.SetSeed(pilot.seedr);

  // open output file, use "C" commands for easier formatting
  FILE *outunit;   /*   Output Unit*/
  if (pilot.outFileName == "") {  
    outunit = stdout;
  }
  else {
    outunit =fopen(pilot.outFileName.c_str(),"w");
  }
  
  // set up unit vectors for primary (the position defined by WobbleX/Y).
  ////////////////////////////////////////////////////////////////////
  // get primary unit vector from its az/zn in ground coordinates,
  double xCosPrG = 0.0;  // primary direction cosines, ground coordinate system
  double yCosPrG = 0.0;  // pointing to the sky, not downward
  double zCosPrG = 0.0;
  double azPr = (pilot.az)*(TMath::DegToRad());  // primary az/zn
  double znPr = (pilot.zn)*(TMath::DegToRad());

  GUtilityFuncts::AzZnToXYcos(azPr, znPr,
                   &xCosPrG, &yCosPrG);

  zCosPrG = sqrt(1 - xCosPrG*xCosPrG - yCosPrG*yCosPrG);

  *oLog << " primary az zn: " << pilot.az << " " << pilot.zn << endl;
  *oLog << " primary direction cosines: " << xCosPrG << " " << yCosPrG << "  " 
        << zCosPrG << endl;

  //////////////////////////////////////////////////////////////////////
  // get photon unit vector when WobbleR = 0.0
  double xCosG = 0.0;  // photon direction cosines, ground coordinate system
  double yCosG = 0.0;  // pointing to the sky, not downward
  double zCosG = 0.0;
  double azPhoton   = 0.0;
  double znPhoton   = 0.0;
  double wobbleE = (pilot.wobbleE)*(TMath::DegToRad()); 
  double wobbleN = (pilot.wobbleN)*(TMath::DegToRad());
  double wobbleR = (pilot.wobbleR)*(TMath::DegToRad());
  double latitude= 90.0*(TMath::DegToRad());
 
  // changng the wobble sign just means that a positive wobbleN (E) in the pilot file results
  // in a larger zenith angle (larger azimuth) for the photon in comparison to the primary.
  GUtilityFuncts::wobbleToAzZn(-wobbleN,-wobbleE,latitude,azPr,znPr,&azPhoton,&znPhoton);

  *oLog << " photon with WobbleR = 0.0, Az Zn " << azPhoton*(TMath::RadToDeg()) << "  " 
        << znPhoton*(TMath::RadToDeg()) << endl;

  // get photon direction cosines (pointing to the sky)
  GUtilityFuncts::AzZnToXYcos(azPhoton, znPhoton,
                   &xCosG, &yCosG);
  zCosG = sqrt(1 - xCosG*xCosG - yCosG*yCosG);
  ROOT::Math::XYZVector photonGDir(xCosG,yCosG,zCosG);
  *oLog << " photon direction cosines (toward sky) WobbleR = 0.0:  " << xCosG << "  " << yCosG << " " << zCosG << endl;
  
  // get photon dir cosines in kascade coordinates.
  double xPhotonDirKas = -xCosG;
  double yphotonDirKas = yCosG;  // reverse direction, invert y axis.
  double zphotonDirKas = -zCosG;
  ROOT::Math::XYZVector photonDirKas(xPhotonDirKas,yphotonDirKas,zphotonDirKas);
  printXYZVector(photonDirKas," photon direction down WobbleR = 0.0 in kascade coordinates:");
  
////////////////////////////////////////////////

  // write header to output
  fprintf(outunit,"* HEADF  <--Start of header lines\n");
  fprintf(outunit,
          "   configuration file %s\n",pilot.configFileName.c_str());
  fprintf(outunit,"   nShowers nPhotons %d  %d\n",pilot.maxShowers,
          pilot.maxPhotons);
  fprintf(outunit,"   primary az zn %f %f\n",pilot.az,pilot.zn);
  fprintf(outunit,"   wobble E N %f %f %f\n",pilot.wobbleE,pilot.wobbleN,pilot.wobbleR);
  fprintf(outunit,"   photon az zn WobbleR = 0.0  %f %f\n",azPhoton*(TMath::RadToDeg()),
          znPhoton*(TMath::RadToDeg()));

  fprintf(outunit,"* DATAF  <-- Start of data\n");
  //////////////////////////////////////////////

  // print R and H lines
  fprintf(outunit,"R 1.0\n");
  fprintf(outunit,"H %f\n",pilot.obser);

  // get rotation matrix for rotating from telescope system to ground system
  //////////////////////////////////////////////////////////////////
  ROOT::Math::Rotation3D rotMat;
  GUtilityFuncts::AzZnToRotMat(azPr,znPr,&rotMat);
  ROOT::Math::Rotation3D rotMatInv = rotMat.Inverse();
  // check rotation matrix:
  //ROOT::Math::XYZVector telPt(xCosPrG,yCosPrG,zCosPrG);
  //printXYZVector(telPt,"  points to primary, ground coor");
  //ROOT::Math::XYZVector telPtTel = rotMat*telPt;
  //printXYZVector(telPtTel,"  points to primary, tel. coor");

  //ROOT::Math::XYZVector telch = rotMatInv*telPtTel;
  //printXYZVector(telch, " after rotMatInv ");
  
  // loop over showers, telescopes, photons
  //////////////////////////////////////////////////////////////////
  for (int j=0;j<pilot.maxShowers;j++) {
    // write S line
    fprintf(outunit,"S 0.0  0.0 0.0 %f %f 1277.06 -1 -1 -1 \n",-xCosPrG, yCosPrG);
    
    //for (unsigned tel = 0;tel < pilot.TelRadius.size(); tel++) {
    for (unsigned tel = 0;tel < 1; tel++) {
      
      for (int ip = 0;ip<pilot.maxPhotons;ip++) {
        
        // find photon hit randomly on telescope plane
        double r = pilot.TelRadius[tel]*sqrt(TR3.Rndm());
        double phi = TR3.Rndm()*(TMath::TwoPi());
        double x = r*cos(phi);
        double y = r*sin(phi);
        double z = 0.0;
        //x = 0.0;
        //y = 1.0;
        //z = 0.0;
        ROOT::Math::XYZVector telHit(x,y,z); // tel hit in telescope coordinates (on tel plane)
        //*oLog << "     telescope number: " << tel << endl;
        //printXYZVector(telHit,"photon hit on telescope");
        
        ROOT::Math::XYZVector telHitTG = rotMatInv*telHit;
        //printXYZVector(telHitTG," photon hit on telescope, telescope ground coor");
        // change origin to array center.
        ROOT::Math::XYZVector telHitG = telHitTG + *(pilot.vTelLoc[tel]);
        //printXYZVector(*(pilot.vTelLoc[tel]), " telescope location ");
        //printXYZVector(telHitG,"      telhit from array Center"); 

        ROOT::Math::XYZVector telHitGKas(telHitG.X(),-telHitG.Y(),telHitG.Z() );
        ////////////////////////////////////////////////
        // photon hit on ground, given the photon hit on the telescope
        ROOT::Math::XYZVector photonHitG;
        photonHitG = telHitG - (telHitG.Z()/zCosG)*photonGDir;
        //printXYZVector(photonHitG," photon hit on the ground wrt array center");

        double xHitKas = telHitGKas.X();
        double yHitKas = telHitGKas.Y();
        double emissionHt = 10000.0;
        if (debug) *oLog << "  photon hit in kas coor. " << xHitKas << " " << yHitKas << endl;

        ///////////// get photon direction /////////////
        double wobbleE1 = wobbleE;  // tmp variables to keep wobbleE/N constant
        double wobbleN1 = wobbleN;
        
        // get wobbleE/N by adding wobbleRand
        if (pilot.wobbleR > 0.01) {
          if (debug) *oLog << " ======== wobbleE/N before " << wobbleE << "  " << wobbleN << endl;
          double wobRrand = wobbleR*sqrt(TR3.Rndm());
          double wobTrand = (TMath::TwoPi())*TR3.Rndm();
          double delE     = wobRrand*cos(wobTrand);
          double delN     = wobRrand*sin(wobTrand);
          wobbleE1     = wobbleE + delE;
          wobbleN1     = wobbleN + delN;
          if (debug) *oLog << " ======== delE/N wobbleE/N after " << delE << "  " << delN
               << "  " << wobbleE1 << "  " << wobbleN1 << endl;
        }
        else {
          if (debug) *oLog << " ===============  wobbleE/N wobbleR=0 " << wobbleE1
                           << "  " << wobbleN1 << endl;
        }
        // changng the wobble sign just means that a positive wobbleN (E) in the pilot file results
        // in a larger zenith angle (larger azimuth) for the photon in comparison to the primary.
        double azPhoton1 = 0.0;
        double znPhoton1 = 0.0;
        GUtilityFuncts::wobbleToAzZn(-wobbleN1,-wobbleE1,latitude,azPr,znPr,&azPhoton1,&znPhoton1);

        if (debug) *oLog << " photon1 Az Zn: " << azPhoton1*(TMath::RadToDeg()) << "  " 
              << znPhoton1*(TMath::RadToDeg()) << endl;

        // get photon direction cosines (pointing to the sky)
        GUtilityFuncts::AzZnToXYcos(azPhoton1, znPhoton1,
                                    &xCosG, &yCosG);
        zCosG = sqrt(1 - xCosG*xCosG - yCosG*yCosG);
        ROOT::Math::XYZVector photonGDir(xCosG,yCosG,zCosG);
        if (debug) *oLog << " photon dir.cosines (toward sky):  " << xCosG
                         << "  " << yCosG << " " << zCosG << endl;
  
        // get photon dir cosines in kascade coordinates.
        double xPhotonDirKas = -xCosG;
        double yphotonDirKas = yCosG;  // reverse direction, invert y axis.
        double zphotonDirKas = -zCosG;
        ROOT::Math::XYZVector photonDirKas(xPhotonDirKas,yphotonDirKas,zphotonDirKas);
        if (debug) printXYZVector(photonDirKas," photon direction down in kascade coordinates");
             
        ///////////////////////////////////////////////
        
        fprintf(outunit,"P %f %f %f %f %f 6 %d 2 %d\n",xHitKas,
                yHitKas, xPhotonDirKas,yphotonDirKas,
                emissionHt,pilot.waveLgt,tel+1);
        
        
      }  // end of photon loop
      
    }  // end of telescope loop
    
  }  // end of shower loop
  
  
};
/*************** end of main *********************/

void readPilot(struct Pilot *pilot,const string &pilotfile) {

  string pilotName = pilotfile;

  pilot->configFileName = "";
  pilot->outFileName    = "";
  pilot->maxShowers     = 1;
  pilot->maxPhotons     = 10000;
  pilot->seedr          = 0;
  pilot->waveLgt        = 500;
  pilot->obser          =  1277.06;
  pilot->az             = 0.0;
  pilot->zn             = 0.0;
  pilot->wobbleN        = 0.0;
  pilot->wobbleE        = 0.0;
  pilot->wobbleR        = 0.0;
  pilot->debug          = 0;

  vector<string> tokens;
  string flag;
  VG_Pilot *pi = new VG_Pilot(pilotName);

  flag = "OUTFL";

  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->outFileName = tokens[0];
  }

  flag = "CONFG";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->configFileName = tokens[0];
  }

  flag = "NUMBR";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->maxShowers = atoi(tokens[0].c_str());
    pilot->maxPhotons = atoi(tokens[1].c_str());
  }

  flag = "SEEDR";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->seedr = (UInt_t) atoi(tokens[0].c_str());
  }
  flag = "OBSER";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->obser = atof(tokens[0].c_str());
  }
  flag = "AZZN";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->az = atof(tokens[0].c_str());
    pilot->zn = atof(tokens[1].c_str());
  }
  flag = "OFFSET";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->wobbleE = atof(tokens[0].c_str());
    pilot->wobbleN = atof(tokens[1].c_str());
    if (tokens.size() < 3) {
      *oLog << " **********************  NEED THREE ENTRIES IN THE OFFSET RECORD"
            << endl;
      *oLog << " ENDING PROGRAM" << endl;
      exit(0);
    }
    pilot->wobbleR = atof(tokens[2].c_str());
  }
  flag = "DEBUG";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    pilot->debug = atoi(tokens[0].c_str());
  }

  *oLog << "     from pilot file: including defaults" << endl;
  *oLog << "          outFileName    " << pilot->outFileName << endl;  
  *oLog << "          configFileName " << pilot->configFileName << endl;  
  *oLog << "          maxShowers     " << pilot->maxShowers << endl;  
  *oLog << "          maxPhotons     " << pilot->maxPhotons << endl;  
  *oLog << "          seed           " << pilot->seedr << endl;
  *oLog << "          obser.height   " << pilot->obser << endl;
  *oLog << "          star azimuth   " << pilot->az << endl;
  *oLog << "          star zenith    " << pilot->zn << endl;
  *oLog << "          wobbleE/N/R    " << pilot->wobbleE << "  "
        << pilot->wobbleN << "  " << pilot->wobbleR << endl;
  *oLog << "          debug          " << pilot->debug << endl;
  delete pi;
};
/*************** end of readPilot *********************/

void readConfig(struct Pilot *pilot) {

  string configFile = pilot->configFileName;

  vector<string> tokens;
  string flag;
  VG_Pilot *pi = new VG_Pilot(configFile);

  flag = "TLLOC";

  // get number of telescopes and make space available
  // in case telescope numbering is out of order
  int numTel = 0;
  numTel = pi->set_flag(flag);
  
  // initialze vectors 
  for (int i = 0;i<numTel;i++) {
    pilot->vTelLoc.push_back(new ROOT::Math::XYZVector );
    pilot->vTelLoc[i]->SetXYZ(0.0,0.0,0.0);
    pilot->TelRadius.push_back(0.0);
  }

  while (pi->get_line_vector(tokens) >= 0) {
    int iTel = atoi(tokens[0].c_str());
    double x = atof(tokens[1].c_str());
    double y = atof(tokens[2].c_str());
    double z = atof(tokens[3].c_str());
   (pilot->vTelLoc[iTel-1])->SetXYZ(x,y,z);
  
  }

  // get mirror radii
  flag = "MIROR";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    int iTel = atoi(tokens[0].c_str());
    double radius = atof(tokens[1].c_str());
    pilot->TelRadius[iTel-1] = radius*1.1;  // increase radius by 10%
                         
  }

  /*
   *oLog << "      configuration file tel.loc / radii" << endl;
   for (int iTel = 0;iTel<numTel;iTel++) {
   double x = pilot->vTelLoc[iTel]->X();
   double y = pilot->vTelLoc[iTel]->Y();
   double z = pilot->vTelLoc[iTel]->Z();
   *oLog << "   " << iTel << " " << x 
   << "  " << y << "  " << z << "  /  " 
   << pilot->TelRadius[iTel] << endl;
   }
  */
  
};
/*************** end of readConfig *********************/

void printXYZVector(const ROOT::Math::XYZVector &vec,const string &label) {
  *oLog << label << "   " << vec.X() << "  " << vec.Y() << " "
        << vec.Z() << endl;

}
/*************** end of printXYZVector *********************/

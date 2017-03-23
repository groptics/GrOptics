// version: 6.1.0    date: 2014-05-06
/*!  makePhotonsPtSrc.cpp
     version 1.0
     Charlie Duke
     Grinnell College

   makePhotonsPtSrc produces .cph file of cherenkov photons as input
   to grisudet for a designated point source of photons. Telescope
   locations and radii are from a standard GrISU telescope configuration
   file. The code has a ROOT dependency.

   All log information goes to cerr.

   Testing of this code is straightforward using a photon_history debug
   file (see the DEBUG record in detector.pilot). This history file may
   be converted to a root tree in GrISU/Utilities/Diagnostics.

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
#include <deque>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>

using namespace std;

#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TRandom3.h"

#include "VG_Pilot.h"

ostream *oLog;

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

struct Pilot {
  string configFileName;
  string outFileName;  // "" means stdout
  int waveLgt;         // int I think
  double phoTime;       // int or float
  int maxShowers;
  int maxPhotons;
  vector<ROOT::Math::XYZVector *> vTelLoc;
  vector<double> TelRadius;
 
  ROOT::Math::XYZVector srcPt; 
  UInt_t seedr;

  double obser; // observatory height, default 1277.06
};

void readPilot(struct Pilot *pilot,const string &pilotfile);
void readConfig(struct Pilot *pilot);
void printXYZVector(const ROOT::Math::XYZVector &vec,const string &label);

int main(int argc, char *argv[]) {

  oLog = &cerr;
  *oLog << "  -- start of makePhotonsPtSrc" << endl;

  string pilotfile("makePhotonsPtSrc.pilot");
  if (argc > 2) {
    cout << " TOO MANY COMMAND LINE ARGUMENTS  " << endl;
    cout << " usage: make_photons [name of pilot file]" << endl;
    cout << "     default pilot file name: makePhotons.pilot" << endl;
    exit(EXIT_FAILURE);
  }  
    else if ( argc == 2 ) {
      pilotfile = argv[1];
    }
  cout << "    pilotfile: " << pilotfile << endl;

  Pilot pilot;
  // read pilot and configuration files
  readPilot(&pilot,pilotfile);
  readConfig(&pilot);

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
  TRandom3 TR3;
  TR3.SetSeed(pilot.seedr);

  // open output file, use "C" commands for easier formatting
  FILE *outunit;   /*   Output Unit*/
  outunit =fopen(pilot.outFileName.c_str(),"w");

  // write header to output
  fprintf(outunit,"* HEADF  <--Start of header lines\n");
  fprintf(outunit,
          "   configuration file %s\n",pilot.configFileName.c_str());
  fprintf(outunit,"   source point %f %f %f\n",pilot.srcPt.X(),
          pilot.srcPt.Y(),pilot.srcPt.Z());
  fprintf(outunit,"   nShowers nPhotons %d  %d\n",pilot.maxShowers,
          pilot.maxPhotons);
  fprintf(outunit,"* DATAF  <-- Start of data\n");

  // print R and H lines
  fprintf(outunit,"R 1.0\n");
  fprintf(outunit,"H %f\n",pilot.obser);
  
  // loop over showers, telescopes, photons
    for (int j=0;j<pilot.maxShowers;j++) {
      // write S line
      fprintf(outunit,"S 0.0 %f %f 0.000 0.0000 1277.06 -1 -1 -1 \n",
              pilot.srcPt.X(),-pilot.srcPt.Y());

      for (unsigned tel = 0;tel < pilot.TelRadius.size(); tel++) {

        for (int ip = 0;ip<pilot.maxPhotons;ip++) {

          // find photon hit randomly on telescope plane
          double r = pilot.TelRadius[tel]*sqrt(TR3.Rndm());
          double phi = TR3.Rndm()*(TMath::TwoPi());
          double x = r*cos(phi) + pilot.vTelLoc[tel]->X();
          double y = r*sin(phi) + pilot.vTelLoc[tel]->Y();
          double z = pilot.vTelLoc[tel]->Z();
          ROOT::Math::XYZVector telHit(x,y,z);  // hit location on telescope

          // find ground hit vector, grdHit
          // vector from srcPt to telescope
          double xSrcToTel = x - pilot.srcPt.X();
          double ySrcToTel = y - pilot.srcPt.Y();
          double zSrcToTel = pilot.vTelLoc[tel]->Z() - pilot.srcPt.Z();
          ROOT::Math::XYZVector srcToTel(xSrcToTel,ySrcToTel,zSrcToTel);

          // unit vector from srcPt to telHit location
          ROOT::Math::XYZVector nSrcToTel = srcToTel.Unit();

          // find ground location
          ROOT::Math::XYZVector grdHit = telHit - 
            (telHit.Z()/nSrcToTel.Z())*nSrcToTel;
          
          //*oLog << endl;
          //printXYZVector(pilot.srcPt,"source location");
          //printXYZVector(*(pilot.vTelLoc[tel]),"telescope location");
          //printXYZVector(telHit,"photon hit on telescope");
          //printXYZVector(nSrcToTel,"photon to Tel unit vector");
          //printXYZVector(grdHit,"ground photon hit");
          
          fprintf(outunit,"P %f %f %f %f %f 6 %d 2 %d\n",grdHit.X(),
                  -grdHit.Y(), nSrcToTel.X(),-nSrcToTel.Y(),
                  pilot.srcPt.Z(),pilot.waveLgt,tel+1);
                  
       
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

  flag = "SRCPT";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    float x = atof(tokens[0].c_str());
    float y = atof(tokens[1].c_str());
    float z = atof(tokens[2].c_str());
    pilot->srcPt.SetXYZ(x,y,z);
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

  *oLog << "     from pilot file: including defaults" << endl;
  *oLog << "          outFileName    " << pilot->outFileName << endl;  
  *oLog << "          configFileName " << pilot->configFileName << endl;  
  *oLog << "          maxShowers     " << pilot->maxShowers << endl;  
  *oLog << "          maxPhotons     " << pilot->maxPhotons << endl;  
  *oLog << "          srcPoint:x/y/z " << pilot->srcPt.X() << " " 
        << pilot->srcPt.Y() << " " << pilot->srcPt.Z() << endl;
  *oLog << "          seed           " << pilot->seedr << endl;
  *oLog << "          obser.height   " << pilot->obser << endl;

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

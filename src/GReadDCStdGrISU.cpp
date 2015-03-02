/*
VERSION3.1
2March2015
*/
/*!  GReadDCStdGrISU.cpp
     Charlie Duke
     Grinnell College
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

#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"
#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "GGeometryBase.h"
#include "GDCGeometry.h"

#include "GTelescope.h"
#include "GDCTelescope.h"

#include "GTelescopeFactory.h"
#include "GDCTelescopeFactory.h"

#include "GOrderedGrid.h"
#include "GReadDCStdBase.h"
#include "GReadDCStdGrISU.h"


#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

GReadDCStdGrISU::GReadDCStdGrISU(const string pilotfile,GDCTelescopeFactory *DCFactory  ) {

  spilotfile = "";
  pi = 0;
  flagline = "";
  iTelNum = 0;
  iGeoNum = 0;
  dFocLgt    = 0.0;
  dCamRad    = 0.0;
  dPlateScaleFactor = 0.0;
  dFocError  = 0.0;
  dTelRadius = 0.0;
  //mVReflWaveLgts = 0;
  //mVCoeffs       = 0;

  DCFac = DCFactory;
  spilotfile = pilotfile;

  setupDCFactory();

};
/****************** end of GReadDCStdGrISU **********/

GReadDCStdGrISU::GReadDCStdGrISU(const string pilotfile) {

  DCFac = 0;
  spilotfile = pilotfile;

};
/****************** end of GReadDCStdGrISU **********/

GReadDCStdGrISU::~GReadDCStdGrISU() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GReadDCStdGrISU::~GReadDCStdGrISU " << endl;
  }
  // factory or telescopes own everything except pi and mtel.
  if (pi) SafeDelete(pi);
 
  map<int,DCStdOptics*>::iterator mTelIter;
  for (mTelIter = mTel.begin();
       mTelIter != mTel.end();
       mTelIter++) {
    
    if (mTelIter->second ) SafeDelete(mTelIter->second);
  }
 
};
/****************** end of ~GReadDCStdGrISU **********/

void GReadDCStdGrISU::
setDCTelescopeFactory(GDCTelescopeFactory *DCFactory) {

  if (DCFac != 0 ) {
    cerr << "GReadDCStdGrISU::setDCTelescopeFactory " << endl;
    cerr << "   ERROR: DCFac pointer to DCFactory previously set" << endl;
    cerr << "   stopping code " << endl;
    exit(0);
  }

  DCFac = DCFactory;
  setupDCFactory();
  
};
/****************** end of setDCTelescopeFactory **********/

void GReadDCStdGrISU::setupDCFactory() {

  int debug = 0;

  if (debug) {
    *oLog << "  -- GReadDCStdGrISU::setupDCFactory" << endl;
  }

  if (DCFac == 0) {
    cerr << " GReadDCStdGrISU doesn't have a DCFactory object" << endl;
    cerr << "exiting code" << endl;
    exit(0);
  }

  pi = new GPilot(spilotfile);

  // check for additional pilot files to append to this file in the pilotreader
  string flag = "PILOTF";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    string newPilotFile = tokens[0];
    pi->addPilotFile(newPilotFile);
  }

  // get pairs of std number and associated telescope number
  // flags are TELSTD followed by telescope number and geometry number

  flag = "TELSTD";
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0) {
    int iStdOptNum = atoi(tokens[0].c_str());
    int iStdOptTel = atoi(tokens[1].c_str());
    int iStdGeoNum = atoi(tokens[2].c_str());
    int iraytracer = atoi(tokens[3].c_str());
    double fAvgTransitTime = atof(tokens[4].c_str());
    mTransitTime[iStdOptNum] = fAvgTransitTime;
    if (tokens.size() >= 6) {
      int iprintmode = atoi(tokens[5].c_str());
      mPrintStdTel[iStdOptNum] = iprintmode;
    }

    // store telid and geoid for each optics standard.
    iStdGeoNum = atoi(tokens[2].c_str());
    mOptics[iStdOptNum] = make_pair(iStdOptTel,iStdGeoNum);
    RayTracerType rte = (RayTracerType)iraytracer;
    mRayTracerEnum[iStdOptNum] = rte;

    if (debug==1){
      DEBUGS(iStdOptNum); DEBUGS(iStdOptTel); 
      DEBUGS(iStdGeoNum); DEBUGS(iraytracer); *oLog << endl;
    }
   
  }

  if (debug==1) {
    DEBUGS(mRayTracerEnum.size());
    DEBUGS(mOptics.size());
  }

  // get list of geometry standards to make
  // get list of telescope numbers used in optstd

  vector<int>::iterator itV;
  for (itOpt=mOptics.begin();itOpt!=mOptics.end();++itOpt) {
    int telnum = (*itOpt).second.first;
    int geonum = (*itOpt).second.second;
    
    itV = find(vTelNum.begin(),vTelNum.end(),telnum);
    if ( itV==vTelNum.end() ) vTelNum.push_back(telnum);

    itV = find(vGeoStd.begin(),vGeoStd.end(),geonum);
    if ( itV==vGeoStd.end() ) vGeoStd.push_back(geonum);
  }

  // make all standard geometries, store in mDCGeo map
  // as GDCGeometry object
  for (unsigned i = 0;i<vGeoStd.size(); i++) {
    iGeoNum = vGeoStd[i];
 
   if (debug==1) {
      *oLog << "       ready to call makeStdGeometry for this iGeoNum" << endl;
      DEBUGS(iGeoNum);
    }
    makeStdGeometry();
  }

  if (debug==1) {
    // see what's made
    *oLog << endl << "       list of  map mGeo keys and types" << endl;
    map<int,GGeometryBase*>::iterator itGeo;
    for (itGeo=DCFac->mDCGeo.begin();itGeo!=DCFac->mDCGeo.end();++itGeo) {
      *oLog << "         " << (*itGeo).first << " "
	    << (*itGeo).second->type << endl;
      (*itGeo).second->printGeometry();
      
    }
    *oLog << "     finished debug print of non-NOSTRUCT geometry classes ";
    *oLog << "from setupDCFactory" << endl<< endl;
  }


  // now make all standard telescopes
  for (unsigned i = 0;i<vTelNum.size(); i++) {
    iTelNum = vTelNum[i];
    makeStdTelescope();
  }

  delete pi;
  pi = 0;
  // get the reflection coefficients, put into mReflC map in DCTelescopeFactory
  getReflCoeff();

  // make standard optics objects for GDCTelescopeFactory
  makeStdOptics();

  if (debug) {
    *oLog << "  -- end of setupDCFactory" << endl;
  }

};
/******************** end of setupDCFactory ****************/

void GReadDCStdGrISU::makeStdTelescope() {

  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GReadDCStdGrISU::makeStdTelescopes" << endl;
    *oLog << "       starting to make telescope  "  
         << iTelNum << endl;
  }
  // make map entry for this telescope standard
  mTel[iTelNum] = new DCStdOptics;
  // just a shortcut to the DCStdOptics pointer in the map.
  DCStdOptics * telp = mTel[iTelNum];
  //telp->geoStruct = 0;  // set pointer to zero. not used here.

  // now fill the *telp class with entries for this standard telescope.

  bool readtelflag = false;  // check for telescope in config file

   string flag = "CAMRAD"; 
   pi->set_flag(flag);
   while (pi->get_line_vector(tokens) >=0) {
     if ( atoi(tokens[0].c_str()) == iTelNum ) {
       telp->camRadius = atof(tokens[1].c_str());
     }
   }

   flag = "PLATESCALE"; 
   pi->set_flag(flag);
   while (pi->get_line_vector(tokens) >=0) {
     if ( atoi(tokens[0].c_str()) == iTelNum ) {
       telp->plateScaleFactor = atof(tokens[1].c_str());
    }
   }

  flag = "GRIDF"; 
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    if ( atoi(tokens[0].c_str()) == iTelNum ) {
      telp->gridOption = atoi(tokens[1].c_str());
      telp->nbinsx = atoi(tokens[2].c_str());
      telp->nbinsy = atoi(tokens[3].c_str());
      if (tokens.size()==5) {
        telp->gridFile = tokens[4];
      }
      //DEBUG(iTelNum); DEBUG(telp->gridOption);DEBUG(telp->nbinsx);
      //DEBUG(telp->nbinsy);DEBUG(telp->gridFile);
    }
  }

  flag = "TLLOC";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    if ( atoi(tokens[0].c_str()) == iTelNum ) {
      // get rotationOffset, pointing offsets
      telp->rotation_offset = atof(tokens[4].c_str());
      // read pointing offsets and convert to radians (from degrees)
      telp->xoff   = (atof(tokens[5].c_str())) / (TMath::RadToDeg());
      telp->yoff   = (atof(tokens[6].c_str())) / (TMath::RadToDeg());

      readtelflag = true;
    }
  }

  if (!readtelflag) {
    cerr << "TELESCOPE NUMBER " << iTelNum << " not present in config.file" << endl;
    cerr << "check config file and TELSTD record " << endl;
    cerr << "exiting code" << endl;
    exit(0);
  }

  flag = "MIROR";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    if ( atoi(tokens[0].c_str()) == iTelNum ) {
      telp->radius = atof(tokens[1].c_str());
      telp->foclength = atof(tokens[2].c_str());
      telp->defoc = atof(tokens[3].c_str());
      telp->nbr_mir =  atoi(tokens[5].c_str());
    }
  }
  
  // get facet details, read into vector
  // create enough space in the vectors for ordered insertion

  (telp->facet).resize(telp->nbr_mir);

  unsigned numf = 0;

  flag = "MIREL";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >= 0) {
    if ( atoi(tokens[0].c_str()) == iTelNum ) {
     DCStdFacet *f = new DCStdFacet;
      numf = atoi(tokens[1].c_str());
      f->facNum = ((int)numf) -1;  // for numbering to start at 0
      f->type   = atoi(tokens[2].c_str());

      if (f->type == 1) {
        f->sides = 1;
      }
      else if (f->type==2) {
        f->sides = 6;
      }
      else if (f->type==3) {
        f->sides = 4;
      }
      else {
        *oLog << "incorrect facet type (not 1,2,or 3)" << endl;
        *oLog << "stopping code " << endl;
        exit(0);
      }

      f->ftprot    = atof(tokens[3].c_str())/(TMath::RadToDeg());
      f->radius = atof(tokens[4].c_str());
      f->curv   = atof(tokens[5].c_str());
      f->xm     = atof(tokens[6].c_str());
      f->ym     = atof(tokens[7].c_str());
      f->zm     = 0.0;
      f->mis_align = atof(tokens[8].c_str())*(TMath::DegToRad());
      f->ffprot    = 0.0;
      f->roughness   = atof(tokens[9].c_str())*(TMath::DegToRad());
      f->reflect     = atof(tokens[10].c_str());
      f->rflctid     = atoi(tokens[11].c_str());

      (telp->facet)[numf - 1] = f;
    }
  }
  //telp->foclength
  //telp->plateScaleFactor

   double tmpPlate = telp->plateScaleFactor;
   if (tmpPlate == 0.0) {
     tmpPlate = tan( 1.0*(TMath::DegToRad())) * (telp->foclength)*100.0;
     telp->plateScaleFactor = tmpPlate;
   }
 

  if ( (int)numf != telp->nbr_mir) {
    *oLog << "incorrect number of facets in config.file" << endl;
    *oLog << numf << " " << telp->nbr_mir << endl;
    exit(0);
  }

};
/******************* end of makeStdTelescope ***************/

void GReadDCStdGrISU::makeStdGeometry() {
  int debug = 0;
  if (debug>0) {
    *oLog << " -- GReadDCStdGrISU::makeStdGeometry" << endl;
    *oLog << "       making standard geometry number: " << iGeoNum << endl;
  }

  GDCGeometry *gp = new GDCGeometry();
  DCFac->mDCGeo[iGeoNum] = gp;
  //mDCGeo[iGeoNum] = gp;
  // in case there are no GEOST records for this standard, set gp->type to NOSTRUCT
  // as a default
  gp->type = NOSTRUCT;

  int geostd = -1;  // used to do a read test later

  string flag = "GEOST";
  pi->set_flag(flag);

  while (pi->get_line_vector(tokens) >=0)  {
    if ( atoi(tokens[0].c_str()) == iGeoNum) {
      if (tokens[1]=="TOPVOL") {
       gp->epsil = atof(tokens[2].c_str());
      }      
      else if (tokens[1]=="FOCBOX") {
        gp->focBoxDim[0] = atof(tokens[2].c_str()); 
        gp->focBoxDim[1] = atof(tokens[3].c_str());
        gp->focBoxDim[2] = atof(tokens[4].c_str());

        gp->focBoxRot[0] = atof(tokens[5].c_str());
        gp->focBoxRot[1] = atof(tokens[6].c_str());
        gp->focBoxRot[2] = atof(tokens[7].c_str());
 
      }
      else if (tokens[1]=="EDGEBOX") {
        int numedge = atoi(tokens[2].c_str()) -1;

        gp->edgeX[numedge] = atof(tokens[3].c_str());
        gp->edgeY[numedge] = atof(tokens[4].c_str());
        gp->edgeZ[numedge] = atof(tokens[5].c_str());
        gp->edgeOffset[numedge] = atof(tokens[6].c_str());
        gp->edgeRot1[numedge] = atof(tokens[7].c_str());
        gp->edgeRot2[numedge] = atof(tokens[8].c_str());
        gp->edgeRot3[numedge] = atof(tokens[9].c_str());

     }
      else if (tokens[1]=="SHUTTER") {

        gp->shutterX = atof(tokens[2].c_str());
        gp->shutterZ = atof(tokens[3].c_str());
        gp->shutterRot1 = atof(tokens[4].c_str()); 
        gp->shutterRot2 = atof(tokens[5].c_str()); 
        gp->shutterRot3 = atof(tokens[6].c_str()); 
      }
      else if (tokens[1]=="QUADARMSIZE") {
        gp->quadArmX = atof(tokens[2].c_str());
        gp->quadArmY = atof(tokens[3].c_str());
        gp->quadArmOffset = atof(tokens[4].c_str());
      }
      else if (tokens[1]=="QUADARM") {
        int numArm = atoi(tokens[2].c_str()) - 1;
        gp->quadArmBottomX[numArm] = atof(tokens[3].c_str());
        gp->quadArmBottomY[numArm] = atof(tokens[4].c_str());
     }
      else if (tokens[1]=="CROSSBAR") {
	int numBar = atoi(tokens[2].c_str()) - 1;
	gp->crossBarX[numBar] = atof(tokens[3].c_str());
	gp->crossBarY[numBar] = atof(tokens[4].c_str());
	gp->crossBarDistBelowFocBox[numBar] = atof(tokens[5].c_str());
      }
      else if (tokens[1]=="STANDARD") {
        geostd = atoi(tokens[2].c_str() );
        GeoType gt = (GeoType)geostd;
        gp->type = gt;
      }      
      else if (tokens[1]=="CAMERARADIUS") {
        gp->cameraRadius = atof(tokens[2].c_str());
      }
    }

  }

  if (geostd < 0) {
    *oLog << "    -- GReadDCStdGrISU::makeStdGeometry " << endl;
    *oLog << "         DID NOT FIND GEOSTD: " << iGeoNum << endl;
    *oLog << "         setting GDCGeometry::type = NOSTRUCT " << endl;
  }
};
/**************** end of makeStdGeometry *************/

void GReadDCStdGrISU::getReflCoeff() {
  int debug = 0;
  if (debug>0) {
    *oLog << " -- GReadDCStdGrISU::getReflCoeff" << endl;
    *oLog << "       opening config file: " << spilotfile << endl;
  }

  // create reflection curve map directly in DCFac.

  DCFac->mVReflWaveLgts = new map<int, vector<double> >;
  DCFac->mVCoeffs  = new map<int, vector<double> >;

  map<int, vector<double> >  *mVR = DCFac->mVReflWaveLgts; 
    map<int, vector<double> > *mVC = DCFac->mVCoeffs; 

  ifstream inFile(spilotfile.c_str(),ios::in);
  if (!inFile) {
    *oLog << " GReadDCStdGrISU::getReflCoeff: cannot open file: " 
	  << spilotfile << endl;
    exit(0);
  }

  string pilotline = "";
  while( getline(inFile,pilotline,'\n')) {

    vector<string> tokens1;
    GUtilityFuncts::tokenizer(pilotline,tokens1);
    if ( (tokens1.size() > 0) &&
         (tokens1[0] == "*")  &&
         (tokens1[1] == "RFCRV") )  {

      int index = atoi(tokens1[2].c_str());
      int number = atoi(tokens1[3].c_str());
 
      map<int, vector<double> >::iterator ivtm;

      ivtm = mVC->find(index);
      if (ivtm!=mVC->end()) {
        cerr << "trying to load reflection curve " << index << endl;
        cerr << "but curve already in the map: check pilot files " << endl;
        exit(0);
      }

      // get ready to read and load reflection coefficients

      (*mVR)[index].resize(number); 
      (*mVC)[index].resize(number);

      // read and store the reflection coefficients
      for (int i = 0;i<number;i++) {
        double wavel     = 0.0;
        double reflcoeff = 0.0;
        inFile >> wavel >> reflcoeff;
        
        ((*mVR)[index])[i] = wavel;
        ((*mVC)[index])[i] = reflcoeff;
      }

      getline(inFile,pilotline);// finish the end of the line
     
      if (debug > 0) {
        *oLog << "      print reflectance vectors for index: " << index 
             << "  number of points " << number << endl;
        for (int i=0;i<number;i++) {
          *oLog << "          " << setw(4) << i << "   " 
               << ((*mVR)[index])[i] << "     "
               <<  ((*mVC)[index])[i] << endl;
        }
        *oLog << endl;
      }        
    }
  }
  if (mVC==0) {
    *oLog << "no reflection coeffient table found in config file" << endl;
    *oLog << "    STOPPING CODE" << endl;
    exit(0);
  }

};
/************* end of getReflCoeff *************/

void GReadDCStdGrISU::makeStdOptics() {

  bool debug = false;
  if (debug) {
    *oLog << "   -- GReadDCStdGrISU::makeStdOptics" << endl;  

    // loop over optics standard, creating each one in DCTelescopeFactory
    *oLog << "        optics_std   tel_std   geo_std" << endl;
    for (itOpt=mOptics.begin();itOpt!=mOptics.end();++itOpt) {
      *oLog << "            " << (*itOpt).first << "           " 
           << (*itOpt).second.first << "         "
           << (*itOpt).second.second << endl;
    }
  }

  for (itOpt=mOptics.begin();itOpt!=mOptics.end();++itOpt) {
    int stdOpt = (*itOpt).first;
    int stdTel = (*itOpt).second.first;
    int stdGeo = (*itOpt).second.second;

    if (debug) {
      *oLog << endl;
      *oLog << "        starting to fill telescope entries for stdOptics: " 
           << stdOpt << endl;
      *oLog << "             stdTel: " << stdTel << "      stdGeo: " 
	    << stdGeo << endl << endl;
    }
    // create new standard optics entry in DCTelescopeFactory
    DCFac->mStdOptics[stdOpt] = new DCStdOptics;
    
    // fill standard optics entry with telescope parameters
    map<int,DCStdOptics*>::iterator itTel;
    itTel = mTel.find(stdTel);
    if (itTel == mTel.end() ) {
      *oLog << "       could not find mTel entry with stdTel = " 
           << stdTel << endl;
      *oLog << "exiting code" << endl;
      exit(0);
    }
 
    // shortcut notation
    DCStdOptics *dco = DCFac->mStdOptics[stdOpt];

    dco->stdType = DC;
    dco->defoc = (*itTel).second->defoc;
    dco->foclength = (*itTel).second->foclength;
    dco->camRadius = (*itTel).second->camRadius;
    dco->plateScaleFactor = (*itTel).second->plateScaleFactor;
    dco->radius = (*itTel).second->radius;
    dco->rotation_offset = (*itTel).second->rotation_offset;
    dco->xoff = (*itTel).second->xoff;
    dco->yoff = (*itTel).second->yoff;
    dco->nbr_mir = (*itTel).second->nbr_mir;
    dco->nbinsx = (*itTel).second->nbinsx;
    dco->nbinsy = (*itTel).second->nbinsy;
    dco->gridOption = (*itTel).second->gridOption;
    dco->gridFile = (*itTel).second->gridFile;
    dco->facet.resize(dco->nbr_mir);

    dco->rayTracerEnum = mRayTracerEnum[stdOpt];
    dco->fAvgTransitTime = mTransitTime[stdOpt];
    dco->iPrtMode = mPrintStdTel[stdOpt];

    for (int i = 0;i<dco->nbr_mir;++i) {
      
      DCStdFacet *ftmp = new DCStdFacet();// for insertion into fac. std optics

      DCStdFacet *of = (*itTel).second->facet[i];  // shortcut notation      
      ftmp->facNum = of->facNum;
      ftmp->type = of->type;
      ftmp->sides = of->sides;
      ftmp->radius = of->radius;

      ftmp->curv = of->curv;
      ftmp->xm = of->xm;

      ftmp->ym = of->ym;
      ftmp->zm = of->zm;
      double fl = dco->foclength;
      double xm = ftmp->xm;
      double ym = ftmp->ym;

      double zm = (fl*fl)- ((xm*xm)+(ym*ym));
      zm = -sqrt(zm); 
      ftmp->zm = zm;
      ftmp->vFacLoc.SetCoordinates(xm,ym,zm);
      
      ftmp->mis_align = of->mis_align;
      ftmp->ftprot = of->ftprot;
      ftmp->ffprot = of->ffprot;
      ftmp->roughness = of->roughness;
      ftmp->reflect = of->reflect;
      ftmp->rflctid = of->rflctid;
      dco->facet[i] = ftmp;
    }

    // add ordered grid
    // need vectors of mirror locations and radii

    vector<double> xm;
    vector<double> ym;
    vector<double> rm;
    int nbinsx = dco->nbinsx;
    int nbinsy = dco->nbinsy;
    string fileNameOrdGrid = dco->gridFile;

    for (unsigned i = 0;i<dco->facet.size();i++) {
      xm.push_back(dco->facet[i]->xm);
      ym.push_back(dco->facet[i]->ym);
      rm.push_back(dco->facet[i]->radius);
    } 

    if (debug) {
      *oLog << "         instantiating grid:  ";
      *oLog << "     nbinsx nbinsy  " << nbinsx << " " 
            << nbinsy << endl;
    }

    int optiongrid = dco->gridOption;
    dco->grid = new GOrderedGrid(xm,ym,rm,nbinsx,nbinsy,
                                 optiongrid, fileNameOrdGrid);

    // ok telescope stuff transfered. Still have to do geometry and ray tracer.
    if (debug) {
      *oLog << "         starting to move Geometry base pointer ";
      *oLog << "entries to stdOptics: " << stdOpt << endl;
      *oLog << "             stdTel: " << stdTel << "      stdGeo: " 
           << stdGeo << endl;
    }
    //dco->geoStruct = mDCGeo[stdGeo];
    dco->igeoStd   = stdGeo;
    if (debug) {
      if ( (DCFac->mDCGeo[stdGeo]->type == DCVERITAS) ) {
	  GDCGeometry *geotest = 
	    dynamic_cast<GDCGeometry*>(DCFac->mDCGeo[stdGeo]);
	  geotest->printGeometry();
      }     
    }
    // print results (see iPrtMode in each DCStrOptics instance
    dco->printDCStdOptics();
  }
};
/************** end of makeStdOptics *****************/

void GReadDCStdGrISU::setPrintMode(ostream &oStr,
                                          const int prtMode) {
  int debug = 0;
  if (debug > 0) {
    *oLog << " -- GReadDCStdGrISU::setPrintMode" << endl;
  }
  oPrtStrm = &oStr;
  iPrtMode = prtMode;
}; 
/******************* end of setPrintMode */

/*  GSCTelescope.cpp

    Akira Okumura 

        and

    Charlie Duke
    Grinnell College
    May 2011
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

using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TRandom3.h"

#include "AOpticsManager.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
#include "AGeoAsphericDisk.h"
#include "AGlassCatalog.h"
#include "TView.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVector3.h"
#include "TGeoPgon.h"
#include "TVersionCheck.h"
#include "TGraph.h"

#include "GUtilityFuncts.h"
#include "GDefinition.h"
#include "GTelescope.h"
#include "GSCTelescope.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGW(x) *oLog << "         " << #x << " = " << x << endl

// define useful units
static const double cm = AOpticsManager::cm();
static const double mm = AOpticsManager::mm();
static const double um = AOpticsManager::um();
static const double nm = AOpticsManager::nm();
static const double  m = AOpticsManager::m();

GSCTelescope::GSCTelescope() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSCTelescope::GSCTelescope() " << endl;
  }

  // initialize parameters
  manager = 0;
  ray = 0;

  fTX = 15.0;
  fTY =  15.0;
  fTZ = 15.0; 

  iTelID          = 1;
  iStdID          = 1;
  fAvgTransitTime = 0.0;
  fphotWaveLgt = 0.0;
  fphotonToTopVolTime = 0.0;

  historyFileName = "";
  historyTreeName = "";
  iHistoryOption   = 0;
  bPhotonHistoryFlag = false;
  fNPoints = 0;
  fStatusLast = 0;
  fRotationOffset = 0.0;

  hisF = 0;
  hisT = 0;
  for (int i = 0;i<3;i++) {
    fInjectLoc[i] = 0.0;
    fInjectDir[i] = 0.0;
  }

  fInjectTime   = 0.0;
  fInjectLambda = 0.0;

  oPrtStrm = oLog;
  iPrtMode = 0;

  gPrimRefl  = 0;
  gSeconRefl = 0;
  fSeconMaxLmda = 0.0;
  fSeconMinLmda = 0.0;
  fPrimMaxLmda  = 0.0;
  fPrimMinLmda  = 0.0;

  eTelType     = SC;
  //initialization();
  fFocLgt = 0.0;
  // primary
  fDp = 0.0;
  fDpinner = 0.0;
  fZp = 0.0;
  
  // secondary
  fZs = 0.0;
  fDs = 0.0;
  fDsinner = 0.0;
  
  // focal plane
  fZf = 0.0;
  fk1 = 0.0;
  fk2 = 0.0;
  
  // camera
  fPixelSize = 0.0;
  fPixelPitch = 0.0;
  fMAPMTWidth = 0.0;
  fMAPMTLength = 0.0;
  fInputWindowThickness = 0.0;
  fMAPMTAngularSize = 0.0;
  
  iNParP = 0.0;
  iNParS = 0.0;
  
};
/********************** end of GSCTelescope *****************/
GSCTelescope::~GSCTelescope() {

  bool debug=false;
  if (debug) {
    *oLog << "  -- GSCTelescope::~GSCTelescope " << endl;
  }
  // close history file if not already closed
  if (hisF != 0) {
    hisF->cd();  
    hisT->Write();
    // deleting the file will close the file;
    // the file owns the tree and will delete it
    SafeDelete(hisF);
 }
  /*
  if (manager != 0) {
    gGeoManager = manager;
    delete manager;
    gGeoManager = 0;
    manager = 0;
  }
  */
 //if (ray) SafeDelete(ray);

};
/********************** end of ~GSCTelescope *****************/

void GSCTelescope::buildTelescope(bool os8)
{

  //TGeoManager *oldgeom = gGeoManager;
  gGeoManager = 0;
  manager = new AOpticsManager("manager","My 3D Project");
  
  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  vacuum->SetTransparency(80);   
  TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);

  /////////////////////////// top volume ///////////////////////////////////
  TGeoVolume *top = manager->MakeBox("top",Air,fTX*m,fTY*m,fTZ*m);
  manager->SetTopVolume(top);
  manager->SetTopVisible(1); 
 
   ////////////// telescope parameters ////////////////////
   // Note: The center of the primary mirror is located at (X, Y, Z) = (0, 0, 0)
   const double kFp      = fFocLgt*m;   // Fp
   const double kDp      = fDp*2*m; // Dp
   const double kDpinner = fDpinner*2*m; // Dpinner
   const double kZp      = fZp*m;        // Primary mirror's position
   const double kZs      = fZs*kFp;  // Secondary mirror's position
   const double kDs      = fDs*2*m; // Ds
   const double kDsinner = fDsinner*2*m; // Dsinner
   const double kZf      = fZf*kFp;  // Focal plane position
   const double kk1      = fk1*mm; // kappa1
   const double kk2      = fk2*mm; // kappa2

  // The actual dead area should be considered carefully
  // See the spec sheet of H8500
  // const double kPixelSize = 5.8*mm; // width of a MAPMT pixel
  // We do not care the dead area for now
  double kPixelSize   = fPixelSize*mm; // width of an MAPMT pixel
  double kPixelPitch  = fPixelPitch*mm; // pitch of MAPMT channels
  double kMAPMTWidth  = fMAPMTWidth*mm;
  double kMAPMTLength = fMAPMTLength*mm; // between input window and anode pins
  double kInputWindowThickness = fInputWindowThickness*mm;

  double kMAPMTAngularSize = fMAPMTAngularSize/60.; // (deg)

   // 14th order polynomials
   // fitted with gnuplot by Akira
   const int kNPar = iNParP;
   //const int kNPar = iNParP;
   
   const double kzp[] = {TMath::Power(kFp,  1)*fzp[0],
   //double kzp[kNPar] = {TMath::Power(kFp,  1)*fzp[0],
   TMath::Power(kFp, -1)*fzp[1],
   TMath::Power(kFp, -3)*fzp[2],
   TMath::Power(kFp, -5)*fzp[3],
   TMath::Power(kFp, -7)*fzp[4],
   TMath::Power(kFp, -9)*fzp[5],
   TMath::Power(kFp,-11)*fzp[6],
   TMath::Power(kFp,-13)*fzp[7]};
   //const double kzs[kNPar] = {TMath::Power(kFp,  1)*6.45608e-08,
   const double kzs[] = {TMath::Power(kFp,  1)*6.45608e-08,
			      TMath::Power(kFp, -1)*-0.416688,
			      TMath::Power(kFp, -3)*-0.144035,
			      TMath::Power(kFp, -5)*0.647955,
			      TMath::Power(kFp, -7)*-2.96087,
			      TMath::Power(kFp, -9)*9.39256,
			      TMath::Power(kFp,-11)*-18.0811,
			      TMath::Power(kFp,-13)*15.3711};
 
 //////////////// Make the primary mirror ////////////////////////
   AGeoAsphericDisk* primaryV = 0;
   primaryV = new AGeoAsphericDisk("primaryV", kZp + kzp[0] - 1*um, 
				   0, kZp + kzp[
   0] , 0, kDp/2., kDpinner/2.); // curvatures are set to infinity (0 = 1./inf)
   primaryV->SetPolynomials(kNPar - 1, &kzp[1], kNPar - 1, &kzp[1]);
   

   // need to set reflection graph here, after each mirror.
   // have default 1.0
   AMirror* primaryMirror = new AMirror("primaryMirror", primaryV);
   if (gPrimRefl!=0) {
     primaryMirror->SetReflectivity(gPrimRefl);
   }

   top->AddNode(primaryMirror, 1);

   ////////////////// make secondary mirror /////////////////////
   AGeoAsphericDisk* secondaryV = 0;
   secondaryV = new AGeoAsphericDisk("secondaryV", 
				     kZs + kzs[0], 0, 
				     kZs + kzs[0]  + 1*um, 
				     0, kDs/2., kDsinner/2.); 
   // curvatures are set to infinity (0 = 1./inf)
   secondaryV->SetPolynomials(kNPar - 1, &kzs[1], kNPar - 1, &kzs[1]);
   
   AMirror* secondaryMirror = new AMirror("secondaryMirror", secondaryV); 
   if (gSeconRefl!=0) {
     secondaryMirror->SetReflectivity(gSeconRefl);
   }
   //secondaryMirror->SetLineColor(kBlue);
   //secondaryMirror->SetFillColor(kBlue);
  top->AddNode(secondaryMirror, 1);

  /////////////////////////////////////////////////////
  ////// make Dummy transparent material //////////////
  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(0); // needed in OpenGL view
  TGeoMedium* med = new TGeoMedium("med", 1, mat);

  //////////////////////////////////////////////////////////////////
  ///// Make a dummy obscuration before the secondary mirror ////////
  AGeoAsphericDisk* secondaryObsV = 0;

  secondaryObsV = new AGeoAsphericDisk("secondaryObsV", kZs + kzs[0] + 10*um, 
				       0, kZs + kzs[0] + 1*um + 10*um, 
				       0, kDs/2., 0);

  secondaryObsV->SetPolynomials(kNPar - 1, &kzs[1], kNPar - 1, &kzs[1]);

  AObscuration* secondaryObs = new AObscuration("secondaryObs", secondaryObsV, med);

  top->AddNode(secondaryObs, 1);

  //////////// Make a single MAPMT photocathode /////////////////////
  //a very thin box
  TGeoBBox* mapmtCathodeV = new TGeoBBox("mapmtCathodeV", kPixelSize/2., 
					 kPixelSize/2., 100*um); 

  AFocalSurface* mapmtCathode = new AFocalSurface("mapmtCathode", 
						  mapmtCathodeV);

  //////////////// Make a single MAPMT //////////////////////////
  TGeoBBox* mapmtV = new TGeoBBox("mapmtV", kMAPMTWidth/2., kMAPMTWidth/2.,
                                  kMAPMTLength/2.);
  AOpticalComponent* mapmt = new AOpticalComponent("mapmt", mapmtV);
  TGeoBBox* mapmtInputWindowV = new TGeoBBox("mapmtInputWindowV",
                                             kMAPMTWidth/2., kMAPMTWidth/2.,
                                             kInputWindowThickness/2.);

  ALens* mapmtInputWindow = new ALens("mapmtInputWindow", 
				      mapmtInputWindowV, med);

  // The refractive index of the input window must be checked
  ARefractiveIndex* bk7 = AGlassCatalog::GetRefractiveIndex("N-BK7");
  mapmtInputWindow->SetRefractiveIndex(bk7);
  mapmt->AddNode(mapmtInputWindow, 1, new TGeoTranslation(0, 
							  0, 
			kMAPMTLength/2. - kInputWindowThickness/2.));

  int n = 1;
  for(int i = 0; i < 8; i++){
    double dx = (i - 3.5)*kPixelPitch;
    for(int j = 0; j < 8; j++){
      double dy = (j - 3.5)*kPixelPitch;
      mapmt->AddNode(mapmtCathode, n, new TGeoTranslation(dx, dy, kMAPMTLength/2. - kInputWindowThickness - 100*um));
      n++;
    } // j
  } // i

  TGeoBBox* mapmtBackObsV = new TGeoBBox("mapmtBackObsV",
                                         kMAPMTWidth/2., kMAPMTWidth/2.,
                                         15*mm);
  AObscuration* mapmtBackObs = new AObscuration("mapmtBackObs", mapmtBackObsV);
  mapmt->AddNode(mapmtBackObs, 1, new TGeoTranslation(0, 
						      0, 
				 -kMAPMTLength/2. + 15*mm));

  /////////////////// Make the focal plane ///////////////////////
  n = 1;
  for(int i = -7; i <= 7; i++){
    double dx = i*kMAPMTWidth;
    for(int j = -7; j <= 7; j++){
      if((TMath::Abs(i) + TMath::Abs(j) >= 11) || (TMath::Abs(i)*TMath::Abs(j) == 21)){
        continue;
      } // if
      double dy = j*kMAPMTWidth;
      double r2 = (i/7.)*(i/7.) + (j/7.)*(j/7.);
      double dz = kk1*r2 + kk2*r2*r2;
      top->AddNode(mapmt, n, new TGeoTranslation(dx, dy, kZf - kMAPMTLength/2. + dz));
      n++;
    } // y
  } // x

  //////////////// close Geometry and reset gGeoManager ////////////////////    
  manager->CloseGeometry();
  
   return;
};
/*************************************************************************************/

void GSCTelescope::injectPhoton(const ROOT::Math::XYZVector &photonLocT,
                                const ROOT::Math::XYZVector &photonDirT,
				const double &photWaveLgt) {
  gGeoManager = manager;

  bool debug = false;
  if (debug) {
    *oLog << " -- GSCTelescope::injectPhoton " << endl;
    *oLog << "       photonLocT  ";
    GUtilityFuncts::printGenVector(photonLocT); *oLog << endl;
    *oLog << "       photonDirT  ";
    GUtilityFuncts::printGenVector(photonDirT); *oLog << endl;
    *oLog << "      bPhotonHistoryFlag " << bPhotonHistoryFlag << endl;
  } 
  // for debugging and testing only
  /*
  fphotonInjectLoc[0] = 3.0;
  fphotonInjectLoc[1] = 0.0;
  fphotonInjectLoc[2] = 0.0;

  fphotonInjectDir[0] = 0.0;
  fphotonInjectDir[1] = 0.0;
  fphotonInjectDir[2] = -1.0;

  *oLog << "        TESTING WITH THESE VALUES" << endl;
  *oLog << "             specified location and direction " << endl;
  for (int i = 0;i<3;i++) {
    *oLog << i << "  " << fphotonInjectLoc[i] << "  " 
  	  << fphotonInjectDir[i] << endl;
  }

  */
  photonLocT.GetCoordinates(fInitialInjectLoc);
  photonLocT.GetCoordinates(fphotonInjectLoc);
  photonDirT.GetCoordinates(fphotonInjectDir); 
  fphotWaveLgt = photWaveLgt;

  // move base of injection location vector to center of primary mirror
  // (add rotation offset)
  fphotonInjectLoc[2] = fphotonInjectLoc[2] + fRotationOffset;

  movePositionToTopOfTopVol();

  // initialize photonHistory parameters if necessary
  if (bPhotonHistoryFlag) {
    initializePhotonHistoryParms();    
  }
  
  // Assuming that three arguments are given in units of (m), (m), (nm)
  double t = 0;
  double x  = fphotonInjectLoc[0];
  double y  = fphotonInjectLoc[1];
  double z  = fphotonInjectLoc[2];
  double dx = fphotonInjectDir[0];
  double dy = fphotonInjectDir[1];
  double dz = fphotonInjectDir[2];

  //if (ray!=0) SafeDelete(ray);
  ray = new ARay(0, photWaveLgt, x*m, y*m, z*m, t, dx, dy, dz);

  gGeoManager = manager;
  //manager->AddTrack(ray);
  manager->TraceNonSequential(*ray);
 
  // Here you can get the traced result
  double new_xyzt[4];
  ray->GetLastPoint(new_xyzt);
  fNPoints = ray->GetNpoints();
  double new_x = new_xyzt[0];
  double new_y = new_xyzt[1];
  double new_z = new_xyzt[2];
  double new_t = new_xyzt[3];
  if (debug) {
    *oLog << " from GetLastPoint " << new_x << " " << new_y 
	  << " " << new_z << " " << new_t << endl;
  }
  //if (manager->GetNtracks() == 27) {
  //manager->DrawTracks();
  //}
};
/********************** end of injectPhoton *****************/
void GSCTelescope::movePositionToTopOfTopVol() {

  gGeoManager = manager;

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSCTelescope::movePositionToTopOfTopVol " << endl;
    *oLog << "        position prior to move to top ";
    *oLog << fphotonInjectLoc[0] << "  " << fphotonInjectLoc[1] << "  " 
	  << fphotonInjectLoc[2] << endl;
  }

  double rfx = fphotonInjectLoc[0];
  double rfy = fphotonInjectLoc[1];
  double rfz = fphotonInjectLoc[2];

  double Z = fTZ; // just inside top volume

  double dl = fphotonInjectDir[0];
  double dm = fphotonInjectDir[1];
  double dn = fphotonInjectDir[2];

  double Rx = rfx - (rfz - Z)*dl/dn;
  double Ry = rfy - (rfz - Z)*dm/dn;
  double Rz = Z;

  fphotonInjectLoc[0] = Rx;
  fphotonInjectLoc[1] = Ry;
  fphotonInjectLoc[2] = Rz;

  // distance traveled from top to inject location
  double dist = (rfz - Z)/dn;
  fphotonToTopVolTime = - dist/(TMath::C());
  if (debug) {
    *oLog << "        TopVolPos in focal point coor.  ";
    for (int i = 0;i<3;i++) {
      *oLog << fphotonInjectLoc[i] << " ";
    }
    *oLog << endl;
    *oLog << "        distance from top to inject loc " << dist << endl;
    *oLog << "        fphotonToTopVolTime " << fphotonToTopVolTime << endl;
    *oLog << endl;
  }
  
};
//****************************************************
bool GSCTelescope::getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
                                           ROOT::Math::XYZVector *photonDcos,
                                           double *photonTime) {
  
  gGeoManager = manager;

  bool debug = false;
  if (debug) {
    *oLog << "  -- GSCTelescope::getCameraPhotonLocation " << endl;
  }
  // just return zeros
  photonLoc->SetCoordinates(0.0,0.0,0.0);
  photonDcos->SetCoordinates(0.0,0.0,-1.0);
  *photonTime = 10.0;

  double x[4],dir[3];
  for (int i = 0;i<3;i++) {
    x[i] = 0.0;
    dir[i] = 0.0;
  }
  x[3] = 0.0;

  ray->GetLastPoint(x);
  ray->GetDirection(dir);
  fStatusLast = -1;

  // convert distances to mm from cm.
  for (int i = 0;i<3;i++) {
    x[i] = x[i] * 10.0;
  }
  photonLoc->SetCoordinates(x[0],x[1],x[2]);
  photonDcos->SetCoordinates(dir[0],dir[1],dir[2]);
  *photonTime = (x[3] + fphotonToTopVolTime)*1.0e09;
 
  enum {kRun, kStop, kExit, kFocus, kSuspend, kAbsorb};
  
  if (ray->IsExited()) fStatusLast = kExit;
  else if (ray->IsFocused()) {
    fStatusLast = kFocus;
  }
  else if (ray->IsStopped()) fStatusLast = kStop;
  else if (ray->IsSuspended()) fStatusLast = kSuspend;
  else if (ray->IsAbsorbed() ) fStatusLast = kAbsorb;
  
  if (debug) {
    *oLog << "       GetLastPOint " << endl;
    for (int i = 0; i < 3; i++) {
      *oLog << "            " << i << "  loc " << x[i] << endl;
      *oLog << "            " << i << "  dir " << dir[i] << endl;
    }
    *oLog <<   "            time " << x[3] << endl;
    *oLog <<   "            fStatus " << fStatusLast << endl;
    *oLog <<   "            enum for fStatus {kRun, kStop, kExit, kFocus, kSuspend, kAbsorb} " << endl;
    *oLog <<   "            *photonTime " << *photonTime << endl;
  }

  if (bPhotonHistoryFlag) {
    for (int i = 0;i<3;i++) {
      fLocLast[i] = x[i];
      fDirLast[i] = dir[i];
    }
    
    fTimeLast = (x[3] + fphotonToTopVolTime)*1.0e09;
    fillPhotonHistory();
  }

  if (ray->IsFocused() ) {
    manager->AddTrack(ray);
  *oLog << "GETNTRACKS  " << manager->GetNtracks() << endl;
  }
  if ( manager->GetNtracks() == 20) {
    manager->DrawTracks("gl");
  }

  return ray->IsFocused();
};
/********************** end of getCameraPhotonLocation *****************/

//void GSCTelescope::setLogFile(const ofstream &logFile) {
//bool debug = false;
//if (debug) {
//  *oLog << " -- GSCTelescope::setLogFile " << logFile << endl;
//}
//
//};
/********************** end of setLogFile *****************/

void GSCTelescope::printTelescope() {

  bool debug = false;
  if (debug) {
    *oPrtStrm << " -- GSCTelescope::printTelescope" << endl;
  }
  if (iPrtMode == 0)  return;

  *oPrtStrm << "       GSCTelescope::printTelescope" << endl;
  *oPrtStrm << "        Telescope Number: " << iTelID << endl;
  *oPrtStrm << "        constructed from standard SC telescope: " 
	    << iStdID << endl << endl;;
  *oPrtStrm << "        AFTER ALL EDITS " << endl;
  *oPrtStrm << "        fAvgTransitTime " << fAvgTransitTime << endl;
  *oPrtStrm << "        fRotationOffset " << fRotationOffset << endl;
  *oPrtStrm << "        Primary   fDp fDpinner fZp " << fDp << " " 
	    << fDpinner << " " << fZp << endl;
  *oPrtStrm << "        Secondary fDs fDsinner fZs " << fDs << " "
	    << fDsinner << " " << fZs << endl;
  *oPrtStrm << "        Foc.Plane fk1 fk2 fZf " << fk1 << " " 
	    << fk2 << " " << fZf << endl;
  *oPrtStrm << "        iNParP " << iNParP << endl;
  if (iNParP > 0) {
    *oPrtStrm << "          primPolyCoeff " << endl;
    for (int i = 0;i<iNParP;i++) {
      *oPrtStrm << "               " << i << "    "  << fzp[i] << endl;
    }
  }
  *oPrtStrm << "        iNParS " << iNParS << endl;
  if (iNParS > 0) {
    *oPrtStrm << "          secondaryPolyCoeff " << endl;
    for (int i = 0;i<iNParS;i++) {
      *oPrtStrm << "               " << i << "    " << fzs[i] << endl;
    }
  }
       
  *oPrtStrm << "        fPixelSize   " << fPixelSize << endl;
  *oPrtStrm << "        fPixelPitch  " << fPixelPitch  << endl;
  *oPrtStrm << "        fMAPMTWidth   " << fMAPMTWidth << endl;
  *oPrtStrm << "        fMAPMTLength   " << fMAPMTLength << endl;
  *oPrtStrm << "        fInputWindowThickness   " << fInputWindowThickness << endl;
  *oPrtStrm << "        fMAPMTAngularSize   " << fMAPMTAngularSize << endl;
  if (gPrimRefl) {
    *oPrtStrm << "      Primary reflection coefficients " << endl;
    Int_t numP = gPrimRefl->GetN();
    for (Int_t i =0;i<numP;i++) {
      double x=0.0; 
      double y=0.0;
      gPrimRefl->GetPoint(i,x,y);
      *oPrtStrm << "        " << i << "     " << x 
		<< "    " << y << endl;
    }
  }
  else {
    *oPrtStrm << "      gPrimRefl TGraph pointer is zero " << endl;      
  }
  if (gSeconRefl) {
    *oPrtStrm << "      Secondary reflection coefficients " << endl;
    Int_t numP = gSeconRefl->GetN();
    for (Int_t i =0;i<numP;i++) {
      double x=0.0; 
      double y=0.0;
      gSeconRefl->GetPoint(i,x,y);
      *oPrtStrm << "        " << i << "     " << x 
		<< "    " << y << endl;
    }
  }
  else {
    *oPrtStrm << "      gSeconRefl TGraph pointer is zero " << endl;      
  }
  
  *oPrtStrm << endl;
};
/********************** end of printTelescope *****************/

void GSCTelescope::drawTelescope() {

  bool debug = false;
  if (debug) {
    *oLog << " ready to draw the SC telescope" << endl; 
  }
  //TGeoManager *managerOld = gGeoManager;
  //gGeoManager = manager;
  gGeoManager = manager;
  gGeoManager->GetMasterVolume()->Draw("ogl");
  //TGeoVolume *mv = manager->GetMasterVolume();
  //mv->Draw("ogl");
  //gGeoManager = managerOld;
};
/********************** end of drawTelescope *****************/

void GSCTelescope::setPrintMode(ostream &oStr,const int prtMode) {
  bool debug = false;
  if (debug) {
    *oLog << " -- GSCTelescope::setPrintMode" << endl;
  }

  oPrtStrm = &oStr;
  iPrtMode = prtMode;
};
/********************** end of setPrintMode *****************/

void GSCTelescope::setPhotonHistory(const string &rootFile,
                                       const string &treeName,
                                       const int &option) {
  bool debug = false;

  bPhotonHistoryFlag = true;

  historyFileName = rootFile;
  historyTreeName = treeName;
  iHistoryOption = option;
  
  if (debug) {
    *oLog << "  -- setPhotonHistory " << endl;
    *oLog << "        historyFileName / rootTree " <<  historyFileName 
          << "  " << historyTreeName << " /option " << option 
          << endl;
  }
  // add telID to historyFileName
  //   make string to insert
  stringstream osTmp;
  osTmp << "_Tel" << iTelID;
  string strInsert = osTmp.str();

  // find insertion point, append if . isn't found
  size_t idx = historyFileName.rfind('.');
  
  if (idx != string::npos) {
    historyFileName.insert(idx,strInsert);
    if (debug) *oLog << historyFileName << endl;
  }
  else {
    historyFileName = historyFileName + strInsert;
  }

  if (debug) {
    *oLog << "     opening photon history file / tree:  " 
          << historyFileName << " / " << historyTreeName << endl;
  }
  hisF = new TFile(historyFileName.c_str(),"RECREATE");
  hisT = new TTree(historyTreeName.c_str(),historyTreeName.c_str());

  makePhotonHistoryBranches();
  initializePhotonHistoryParms();
};
/********************** end of setPhotonHistory *****************/

void GSCTelescope::makePhotonHistoryBranches() {

  hisF->cd();
  bool debug = false;
  if (debug) {
    *oLog << "  -- makePhotonHistoryBranches " << endl;
  }

  hisT->Branch("status",&fStatusLast,"status/I");
  hisT->Branch("nPoints",&fNPoints,"nPoints/I");
  hisT->Branch("injectX",&fInitialInjectLoc[0],"injectX/D");
  hisT->Branch("injectY",&fInitialInjectLoc[1],"injectY/D");
  hisT->Branch("injectZ",&fInitialInjectLoc[2],"injectZ/D");

  hisT->Branch("xLast",&fLocLast[0],"xLast/D");
  hisT->Branch("yLast",&fLocLast[1],"xLast/D");
  hisT->Branch("zLast",&fLocLast[2],"xLast/D");
  hisT->Branch("xLastDir",&fDirLast[0],"xLastDir/D");
  hisT->Branch("yLastDir",&fDirLast[1],"xLastDir/D");
  hisT->Branch("zLastDir",&fDirLast[2],"xLastDir/D");
  hisT->Branch("timeLast",&fTimeLast,"timeLast/D");
};
/************************* end of makePhotonHistoryBranches *****/

void GSCTelescope::fillPhotonHistory() {

  // careful have to move to correct hisF ?, do a changedirectory?

  bool debug = false;
  if (debug) {
    *oLog << "  -- fillPhotonHistory " << endl;
  }
  hisF->cd();
  hisT->Fill();
};
/************************* end of fillPhotonHistory *****/

void GSCTelescope::initializePhotonHistoryParms() {

};
/************************* end of initializePhotonHistoryParms *****/

void GSCTelescope::writePhotonHistory() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- in GSCTelescope::writePhotonHistory " << endl;
  }
  if (hisF != 0) {
    hisF->cd();  
    hisT->Write();
    // deleting the file will close the file;
    // the file owns the tree and will delete it
    SafeDelete(hisF);
 }
};
/************************* end of writePhotonHistory *****/
void GSCTelescope::setReflCv(const int &primID, const int &secondID,
			     map<int, TGraph *> *mGRefl) {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GSCTelescope::setReflCv" << endl;
  }
  // extract the graphs from the maps 
  //USE COPY CONSTRUCTOR TO PRODUCE INDEPENDENT GRAPHS
  TGraph *tmp1 = (*mGRefl)[primID];
  gPrimRefl = new TGraph(*tmp1);
  gPrimRefl->SetTitle("primary reflectance");

  TGraph *tmp2 = (*mGRefl)[secondID];
  gSeconRefl = new TGraph(*tmp2);
  gSeconRefl->SetTitle("secondary reflectance");

  // get the number of points
  int primPts =  gPrimRefl->GetN();
  int seconPts = gSeconRefl->GetN();

  // get min and max wavelength (point 0 and n-1)
  double x = 0.0;
  double y = 0.0;
  gPrimRefl->GetPoint(0,x,y);  
  fPrimMinLmda = x;

  gPrimRefl->GetPoint(primPts-1,x,y);  
  fPrimMaxLmda = x;
  gSeconRefl->GetPoint(0,x,y);  
  fSeconMinLmda = x;
  gSeconRefl->GetPoint(seconPts-1,x,y);  
  fSeconMaxLmda = x;
  
  if (debug) {
    *oLog << "       primary id/size/min/maxlmda   " << primID << "  " 
  	  << primPts << " " << fPrimMinLmda << " " 
  	  << fPrimMaxLmda << endl;
    
    *oLog << "       secondary id/size/ " << secondID << "  " 
  	  << seconPts << " " << fSeconMinLmda << " " 
  	  << fSeconMaxLmda << endl;

    // plot the reflectance graphs
    TCanvas *cref = new TCanvas("cref","reflCoeff");
   
    cref->Divide(1,2);
    cref->cd(1);
    gPrimRefl->Draw("AC*");

    cref->cd(2);
    gSeconRefl->Draw("AC*");
  }

};
/************************* end of writePhotonHistory *****/

void GSCTelescope::initialization() {
  bool debug = false;

  if (debug) {
    *oLog << "  -- GSCTelescope::initialization " << endl;
  }
  //for (int i = 0;i<3;i++) {
  //fLastLoc[i] = 0.0;
  //}

};

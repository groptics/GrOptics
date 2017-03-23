/*
VERSION4.0
30May2016
*/
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
#include "TPolyLine3D.h"

#include "AOpticsManager.h"
#include "ABorderSurfaceCondition.h"

#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
//#include "TGeoRotation.h"
#include "TGeoMatrix.h"
//#include "TGeoCombiTrans.h"
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
#include "ADefinition.h"
#include "ATelescope.h"
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
  /*
  for (int i = 0;i<3;i++) {
    fInjectLoc[i] = 0.0;
    fInjectDir[i] = 0.0;
  }

  fInjectLambda = 0.0;
  */
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
  fPlateScaleFactor = 0.0;
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
  fMAPMTOffset = 0.0;
  fMAPMTGap = 0.0;
  fMAPMTRefIndex = 0.0;
  bCameraFlag = true;

  bRayPlotModeFlag = false;

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
  if (manager != 0) {
    gGeoManager = manager;
    SafeDelete(manager);
  }

  SafeDelete(ray);

};
/********************** end of ~GSCTelescope *****************/

void GSCTelescope::buildTelescope()
{
  
  gGeoManager = 0;
  manager = new AOpticsManager("manager","My 3D Project");
  
  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  vacuum->SetTransparency(80);   
  TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);

  /////////////////////////// top volume ///////////////////////////////////
  // fTX/Y/Z set to 15 meters in constructor
  TGeoVolume *top = manager->MakeBox("top",Air,fTX*m,fTY*m,fTZ*m);
  manager->SetTopVolume(top);
  manager->SetTopVisible(0); 
  //manager->SetTopVisible(1); 
 
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

  // used variable.
  //double kMAPMTAngularSize = fMAPMTAngularSize/60.; // (deg)
  double kMAPMTOffset = fMAPMTOffset*mm;
  const double kgap     = fMAPMTGap*mm; // gap between input window and photocathode
  double kMAPMTRefIndex = fMAPMTRefIndex;
  bool kCameraFlag = bCameraFlag;

  *oLog << " CameraFlag " << bCameraFlag << endl;
  // 14th order polynomials
  // fitted with gnuplot by Akira
  const int kNPar = iNParP;

  bool printcoeff = true;
 
  double kzp[iNParP];
  kzp[0] = TMath::Power(kFp,  1)*fzp[0];
  if (printcoeff) {
    *oLog << " primary ip fpPower  fzp[ip] kzp[ip] " << "0" 
          << "   " << "1" << "   " << fzp[0] << "   " 
          << kzp[0] << endl;
  }
   int fpPower = -1;
  for (int ip = 1;ip<iNParP;ip++) {
    kzp[ip] = TMath::Power(kFp, fpPower)*fzp[ip];
    if (printcoeff) {
      *oLog << " primary ip fpPower  fzp[ip] kzp[ip] " << ip 
            << "   " << fpPower << "   " << fzp[ip] << "   " 
            << kzp[ip] << endl;
    }
    fpPower = fpPower -2;
  }

  //const double kzp[] = {TMath::Power(kFp,  1)*fzp[0],
  //double kzp[] = {TMath::Power(kFp,  1)*fzp[0],
  //			//double kzp[kNPar] = {TMath::Power(kFp,  1)*fzp[0],
  //			TMath::Power(kFp, -1)*fzp[1],
  //			TMath::Power(kFp, -3)*fzp[2],
  //			TMath::Power(kFp, -5)*fzp[3],
  //			TMath::Power(kFp, -7)*fzp[4],
  //			TMath::Power(kFp, -9)*fzp[5],
  //			TMath::Power(kFp,-11)*fzp[6],
  //			TMath::Power(kFp,-13)*fzp[7]};
  //const double kzs[kNPar] = {TMath::Power(kFp,  1)*6.45608e-08,


  double kzs[iNParS];
  kzs[0] = TMath::Power(kFp,  1)*fzs[0];
  if (printcoeff) {
    *oLog << " secondary ip fsPower  fzs[ip] kzs[ip] " << "0" 
          << "   " << "1" << "   " << fzs[0] << "   " 
          << kzs[0] << endl;
  }
  int fsPower = -1;
  for (int ip = 1;ip<iNParS;ip++) {
    kzs[ip] = TMath::Power(kFp, fsPower)*fzs[ip];
    
    if (printcoeff) {
    *oLog << " secondary ip fsPower  fzs[ip] kzs[ip] " << ip 
          << "   " << fsPower << "   " << fzs[ip] << "   " 
          << kzs[ip] << endl;
    }
    fsPower = fsPower -2;
  }

  //  const double kzs[] = {TMath::Power(kFp,  1)*fzs[0],
  //			TMath::Power(kFp, -1)*fzs[1],
  //			TMath::Power(kFp, -3)*fzs[2],
  //			TMath::Power(kFp, -5)*fzs[3],
  //			TMath::Power(kFp, -7)*fzs[4],
  //			TMath::Power(kFp, -9)*fzs[5],
  //			TMath::Power(kFp,-11)*fzs[6],
  //			TMath::Power(kFp,-13)*fzs[7]};

  //  const double kzs[] = {TMath::Power(kFp,  1)*6.45608e-08,
  //			TMath::Power(kFp, -1)*-0.416688,
  //			TMath::Power(kFp, -3)*-0.144035,
  //			TMath::Power(kFp, -5)*0.647955,
  //			TMath::Power(kFp, -7)*-2.96087,
  //			TMath::Power(kFp, -9)*9.39256,
  //			TMath::Power(kFp,-11)*-18.0811,
  //			TMath::Power(kFp,-13)*15.3711};
  
  //////////////// Make the primary mirror ////////////////////////
  AGeoAsphericDisk* primaryV = 0;
  
  // kZp is normally 0.0 from stdSCTelescopes.cfg configuration file
  // include kZp in TGeoCombiTrans below
  primaryV = new AGeoAsphericDisk("primaryV", 
				  kzp[0] - 1*um, 0, 
				  kzp[0] , 0, 
				  kDp/2., kDpinner/2.); 
  // curvatures are set to infinity (0 = 1./inf)
  
  primaryV->SetPolynomials(kNPar - 1, &kzp[1], kNPar - 1, &kzp[1]);
  
  
  // need to set reflection graph here, after each mirror.
  // have default 1.0
  AMirror* primaryMirror = new AMirror("primaryMirror", primaryV);
  if (gPrimRefl!=0) {
    primaryMirror->SetReflectivity(gPrimRefl);
  }

  /*  
  TGeoRotation *rPrim = new TGeoRotation("rPrim", 
					 fPrimPhiAngle,fPrimThetaOffset,0.0);
  TGeoCombiTrans *cPrim = new TGeoCombiTrans("cPrim",
					     fPrimXoffset*mm,
					     fPrimYoffset*mm,
					     kZp + fPrimZoffset*mm,rPrim);
  top->AddNode(primaryMirror, 1,cPrim);
  */

  // valgrind shows a memory leak here, but not in other similar addNode calls.
  top->AddNode(primaryMirror, 1, new TGeoCombiTrans("cPrim",
						    fPrimXoffset*mm,
						    fPrimYoffset*mm,
						    kZp + fPrimZoffset*mm,
						    new TGeoRotation("rPrim", 
								     fPrimPhiAngle,fPrimThetaOffset,0.0) ));


  cout << " adding primary roughness sigma = " << fPrimRoughSigma << endl;
  //ABorderSurfaceCondition * condition = new ABorderSurfaceCondition(top,primaryMirror);
  //condition->SetGaussianRoughness((fPrimRoughSigma/60.)*TMath::DegToRad());
  ////////////////// make secondary mirror /////////////////////
  // make a new volume to place secondary mirror and obscuration later
  // may want to reduce the depth of the box; check distance to focal surface 
  TGeoVolume *secVol = gGeoManager->MakeBox("secVol",Air,kDs/2.,kDs/2.,150.*cm);
  
  AGeoAsphericDisk* secondaryV = 0;

  secondaryV = new AGeoAsphericDisk("secondaryV", 
				    kzs[0], 0, 
				    kzs[0]  + 1*um, 0, 
				    kDs/2., kDsinner/2.); 
  // curvatures are set to infinity (0 = 1./inf)
  secondaryV->SetPolynomials(kNPar - 1, &kzs[1], kNPar - 1, &kzs[1]);
  
  AMirror* secondaryMirror = new AMirror("secondaryMirror", secondaryV); 
  if (gSeconRefl!=0) {
    secondaryMirror->SetReflectivity(gSeconRefl);
  }
  //secondaryMirror->SetLineColor(kBlue);
  //secondaryMirror->SetFillColor(kBlue);
  secVol->AddNode(secondaryMirror, 1);

  cout << " adding secondary roughness sigma = " << fSecondRoughSigma << endl;    
  //ABorderSurfaceCondition * conditionSec = new ABorderSurfaceCondition(secVol,secondaryMirror);
  //conditionSec->SetGaussianRoughness((fSecondRoughSigma/60.)*TMath::DegToRad());
   /////////////////////////////////////////////////////
  ////// make Dummy transparent material //////////////

  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(0); // needed in OpenGL view
  TGeoMedium* med = new TGeoMedium("med", 1, mat);

  //////////////////////////////////////////////////////////////////
  ///// Make a dummy obscuration before the secondary mirror ////////
  AGeoAsphericDisk* secondaryObsV = 0;

  secondaryObsV = new AGeoAsphericDisk("secondaryObsV", kzs[0] + 10*um, 
				       0, kzs[0] + 1*um + 10*um, 
				       0, kDs/2., 0);

  secondaryObsV->SetPolynomials(kNPar - 1, &kzs[1], kNPar - 1, &kzs[1]);

  AObscuration* secondaryObs = new AObscuration("secondaryObs", secondaryObsV, med);

  
  secVol->AddNode(secondaryObs, 1);
 
  /*
  TGeoRotation *rSec = new TGeoRotation("rSec", 
  					 fSecondPhiAngle,
  					 fSecondThetaOffset,0.0);
  TGeoCombiTrans *cSec = new TGeoCombiTrans("cSec",
  					     fSecondXoffset*mm,
  					     fSecondYoffset*mm,
  					     kZs + fSecondZoffset*mm,rSec);
   top->AddNode(secVol, 1,cSec);
  */
  top->AddNode(secVol,1, new TGeoCombiTrans("cSec",
					    fSecondXoffset*mm,
					    fSecondYoffset*mm,
					    kZs + fSecondZoffset*mm,
					    new TGeoRotation("rSec", 
							     fSecondPhiAngle,
							     fSecondThetaOffset,0.0) ));

  ////////////////// make the camera or the focal surface with no camera ///////////////////
  // make a new volume for the camera or focal plane
  // size adequately covers os8 camera/focal surface 
  TGeoVolume *focVol = gGeoManager->MakeBox("focVol",Air,50.*cm,50.*cm,10.*cm);
 
  if (kCameraFlag ) {
    //////////// Make a single MAPMT photocathode /////////////////////
    //a very thin box
    *oLog << "      ----  building a single MAPMT " << endl;

    TGeoBBox* mapmtCathodeV = new TGeoBBox("mapmtCathodeV", kPixelSize/2., 
					   kPixelSize/2., 100*um); 

    *oLog << "       1/2 size of mapmtCathodeV " << setw(12) << setprecision(6) 
	  << kPixelSize/2. << "  " 
	  << kPixelSize/2. << "  " << 100*um << endl;

    AFocalSurface* mapmtCathode = new AFocalSurface("mapmtCathode", 
						    mapmtCathodeV);

    //////////////// Make a single MAPMT //////////////////////////
    TGeoBBox* mapmtV = new TGeoBBox("mapmtV", kMAPMTWidth/2., kMAPMTWidth/2.,
				    kMAPMTLength/2.);
    AOpticalComponent* mapmt = new AOpticalComponent("mapmt", mapmtV);
    TGeoBBox* mapmtInputWindowV = new TGeoBBox("mapmtInputWindowV",
					       kMAPMTWidth/2., kMAPMTWidth/2.,
					       kInputWindowThickness/2.);
 
    *oLog << "       1/2 size of mapmt " << setw(12) << setprecision(6)
	  << kMAPMTWidth/2. << "  " << kMAPMTWidth/2. << "  " << kMAPMTLength/2.
	  << endl;
    *oLog << "       1/2 size of mapmtInputWindowV " << setw(12) << setprecision(6) 
	  << kMAPMTWidth/2. << "  " 
	  << kMAPMTWidth/2. << "  " << kInputWindowThickness/2. << endl;
    *oLog << "       gap between bottom of window and cathode " << kgap << endl;
 
    ALens* mapmtInputWindow = new ALens("mapmtInputWindow", 
					mapmtInputWindowV, med);
    
    // The refractive index of the input window must be checked
    *oLog << "       kMAPMTRefIndex  " << kMAPMTRefIndex << endl;
   if (kMAPMTRefIndex) {
      mapmtInputWindow->SetConstantRefractiveIndex(kMAPMTRefIndex);
    }
    else {
      ARefractiveIndex* bk7 = AGlassCatalog::GetRefractiveIndex("N-BK7");
      mapmtInputWindow->SetRefractiveIndex(bk7);
    }

    mapmt->AddNode(mapmtInputWindow, 1, new TGeoTranslation(0, 
							    0, 
		     kMAPMTLength/2. - kInputWindowThickness/2.));
 
   *oLog << "       placement of input window  " 
	  << setw(12) << setprecision(6)
 	  << kMAPMTLength/2. - kInputWindowThickness/2. << endl;

    int n = 1;
    for(int i = 0; i < 8; i++){
      double dx = (i - 3.5)*kPixelPitch;
      for(int j = 0; j < 8; j++){
	double dy = (j - 3.5)*kPixelPitch;
	mapmt->AddNode(mapmtCathode, n, new TGeoTranslation(dx, dy, kMAPMTLength/2. 
							    - kInputWindowThickness - kgap - (100.)*um));
	n++;
      } // j
    } // i
    *oLog << "       placement of photocathode " 
	  << setw(12) << setprecision(6) 
	  << kMAPMTLength/2. - kInputWindowThickness - kgap - (100.)*um << endl;
 
    
    TGeoBBox* mapmtBackObsV = new TGeoBBox("mapmtBackObsV",
					   kMAPMTWidth/2., kMAPMTWidth/2.,
					   0.5*mm);

    *oLog << "       1/2 size of backObscuration " << setw(12) << setprecision(6) 
	  << kMAPMTWidth/2.
	  << "   " << kMAPMTWidth/2. << "  " << 0.5*mm << endl;

    AObscuration* mapmtBackObs = new AObscuration("mapmtBackObs", mapmtBackObsV);
    mapmt->AddNode(mapmtBackObs, 1, new TGeoTranslation(0, 
							0, 
							-kMAPMTLength/2. + 0.5*mm));

    *oLog << "       placement of backObscuration" << setw(12) << setprecision(6) 
	  << -kMAPMTLength/2. + 0.5*mm << endl;
    /////////////////// Make the focal plane ///////////////////////
    n = 1;
    ofstream ofile("mapmtLocations.txt");
    bool prtMAPMT = true;
    if (prtMAPMT) {
      ofile << " -- GSCTelescope::buildTelescope, prtMAPMT " << endl;
      ofile << "               all distances in centimeters" << endl;
      ofile << "        kMAPMTWidth / kMAPMTLength " << kMAPMTWidth 
	    << "  /  " << kMAPMTLength << endl;
      ofile << "        i   j     dx        dy      dz    r" << endl;
    }
    for(int i = -7; i <= 7; i++){
      double dx = i*kMAPMTWidth;
      for(int j = -7; j <= 7; j++){
	if((TMath::Abs(i) + TMath::Abs(j) >= 11) || (TMath::Abs(i)*TMath::Abs(j) == 21)){
	  continue;
	} // if
	double dy = j*kMAPMTWidth;
	double r2 = (i/7.)*(i/7.) + (j/7.)*(j/7.);
	double dz = kk1*r2 + kk2*r2*r2;

	if (prtMAPMT) {
	  ofile << "      " << i << "   " << j  << "      " << dx
		<< "     " << dy << "     " << dz << "     "
		<< sqrt(dx*dx + dy*dy) << endl;
	}
	
	//top->AddNode(mapmt, n, new TGeoTranslation(dx, dy, kZf - kMAPMTLength/2. + dz));
	
	// add kZf later
	//top->AddNode(mapmt, n, new TGeoTranslation(dx, dy, kZf + dz 
	//					   - kMAPMTLength/2.
	//					   + kInputWindowThickness + kgap 
	//					   + kMAPMTOffset) );
	focVol->AddNode(mapmt, n, new TGeoTranslation(dx, dy, dz 
						   - kMAPMTLength/2.
						   + kInputWindowThickness + kgap 
						   + kMAPMTOffset) );
	n++;
      } // y
    } // x

    TGeoRotation *rFoc = new TGeoRotation("rFoc", 
					 fFocalPlPhiAngle,fFocalPlThetaOffset,0.0);
    TGeoCombiTrans *cFoc = new TGeoCombiTrans("cFoc",
					     fFocalPlXoffset*mm,
					     fFocalPlYoffset*mm,
					     kZf + fFocalPlZoffset*mm,rFoc);
    top->AddNode(focVol, 1,cFoc);
 

    //*oLog << "     --  placing modules " << endl;
    //*oLog << "       kZf   " << setw(12) << setprecision(6) << kZf << endl;
    //*oLog << "       kMAPMTOffset "  << kMAPMTOffset << endl;
    
    //*oLog << "       location of module with dz = 0  " << setw(15) << setprecision(8)
    //	  << kZf  - kMAPMTLength/2. + kInputWindowThickness + kgap + kMAPMTOffset << endl;
  }
  else {
    //*oLog << " ready to make focal plane with no camera " << endl;
    //*oLog << " kZf " << kZf << endl;
    //*oLog << " kZf - 1*mm " << kZf-1*mm << endl;
    //*oLog << " 38.*cm " << 38.*cm << endl;
    //*oLog << " kMAPMTWidth " << kMAPMTWidth << endl;
    //*oLog << " kk1  kk2 " << kk1 << "  " << kk2 << endl;
    Double_t tmp = 1/(7.*kMAPMTWidth);
    Double_t tmp2 = tmp*tmp;
    Double_t beta[2];
    beta[0] = tmp2*kk1;
    beta[1] = tmp2*tmp2*kk2;
    //*oLog << "  beta[0]  beta[1] " << beta[0] << "  " << beta[1] << endl;
    //double planeRad = 38.*cm;
    double planeRad = 38.*cm;
    //*oLog << " planeRad " << planeRad << endl;

    //AGeoAsphericDisk* focalV = new AGeoAsphericDisk("focalV", kZf + kzs[0] - 1*mm, 0, 
    //						    kZf + kzs[0], 0, planeRad, 0.);
    AGeoAsphericDisk* focalV = new AGeoAsphericDisk("focalV", kzs[0] - 1*mm, 0, 
						    kzs[0], 0, planeRad, 0.);
    
     // curvatures are set to infinity (0 = 1./inf)
    focalV->SetPolynomials(2, beta, 2, beta);
    AFocalSurface* focalPlane = new AFocalSurface("focalPlane", focalV);
    focVol->AddNode(focalPlane, 1);

    // make an obscuration disk before the focal plane
    //AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV",
    //                  kZf + kzs[0] - 2*mm, 0, kZf + kzs[0] - 1*mm, 0, planeRad, 0.);
    AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV",
                      kzs[0] - 2*mm, 0, kzs[0] - 1*mm, 0, planeRad, 0.);
 
    focalObsV->SetPolynomials(2, beta, 2, beta);

    AObscuration* focalObs = new AObscuration("focalObs", focalObsV);

    focVol->AddNode(focalObs, 1);
    /*
    TGeoRotation *rFocS = new TGeoRotation("rFocS", 
					   fFocalPlPhiAngle,
					   fFocalPlThetaOffset,0.0);
    TGeoCombiTrans *cFocS = new TGeoCombiTrans("cFocS",
					       fFocalPlXoffset*mm,
					       fFocalPlYoffset*mm,
					       kZf + fFocalPlZoffset*mm,rFocS);
    top->AddNode(focVol, 1,cFocS);
    */

    top->AddNode(focVol, 1,new TGeoCombiTrans("cFocS",
					      fFocalPlXoffset*mm,
					      fFocalPlYoffset*mm,
					      kZf + fFocalPlZoffset*mm,
					      new TGeoRotation("rFocS", 
							       fFocalPlPhiAngle,
							       fFocalPlThetaOffset,0.0)));

    //*oLog << " focal plane positions: " << kZf + kzs[0] - 1*mm << "  " << kZf + kzs[0] << endl;
    //*oLog << " obs plane positions:   " << kZf + kzs[0] - 2*mm << "  " << kZf + kzs[0]-1*mm << endl;

    
  }

  //////////////// close Geometry and reset gGeoManager ////////////////////    
  manager->CloseGeometry();
  //drawTelescope();
  return;
};
/*************************************************************************************/
void GSCTelescope::testFocalPlane() {

  // inject photons 

};

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

  // convert to cm as required for robast
  fphotWaveLgt = photWaveLgt*nm;

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

  SafeDelete(ray);
  ray = new ARay(0, fphotWaveLgt, x*m, y*m, z*m, t, dx, dy, dz);

  gGeoManager = manager;

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

  // the following will draw polylines that end up on the camera/focal plane.
  // you have to instantiate app and do a app.run in grOptics.cpp by 
  // uncommenting these lines in grOptics (shortly after the start of main
  // shortly before the end of main.
  // This also works with the testtel option; however, you have to set
  // nPhotons to a small number, e.g. 10 or less in GArrayTel.cpp (in 
  // the testTelescope method. The TPolyLine3D::Print("all") will print
  // the start of the line, the intermediate points, and the end of the 
  // line.  A good way to test the code.

  //*oLog << "bRayPlotModeFlag " << bRayPlotModeFlag << endl;

  if (bRayPlotModeFlag) {
    static int idraw = 1;
    if (idraw) {
      drawTelescope();
      idraw = 0;
    }
    // do we draw the ray.
    if ( ( (eRayPlotType == FOCUSONLY) && (fStatusLast == 3) ) ||
         (eRayPlotType == ALLSURFACES) ) {
      

      /*  following was in SCTelescope, replace with SegSC code, 1Jan2017 C.Duke
      *oLog << " ready to draw polyline" << endl;
      TPolyLine3D *pol = ray->MakePolyLine3D();
      pol->Print("all");
      cout << " fStatusLast " << fStatusLast << endl;
      pol->SetLineColor(2);
      pol->Draw();
      */

      *oLog << " ready to draw polyline" << endl;
      *oLog << "      ray->GetNpoints(): " << ray->GetNpoints() << endl;
      TPolyLine3D *pol = ray->MakePolyLine3D();
      for (Int_t ii = 0 ;ii < pol->GetN(); ii++ ) {
	*oLog << ii << "  " << (pol->GetP())[ii] << endl;
      }
      pol->Print("all");
      *oLog << " fStatusLast " << fStatusLast << "  ";
      if (fStatusLast == 0) {
	*oLog << "kRun" << endl;
      }
      else if (fStatusLast == 1) {
	*oLog << "kStop" << endl;
      }
      else if (fStatusLast == 2) {
	*oLog << "kExit" << endl;
      }
      else if (fStatusLast == 3) {
	*oLog << "kFocus" << endl;
      }
      else if (fStatusLast == 4) {
	*oLog << "kSuspend" << endl;
      }
      else if (fStatusLast == 5) {
	*oLog << "kAbsorb" << endl;
      }
      else if (fStatusLast == 1) {
	*oLog << "kStop" << endl;
      }
      else {
	*oLog << " can't interpret fStatusLast" << endl;
      }
      *oLog << endl;
      
      pol->SetLineColor(2);
      pol->Draw();
      gPad->Update();
    }
  }

  return ray->IsFocused();

};
/********************** end of getCameraPhotonLocation *****************/

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
  *oPrtStrm << "        fFocLgt " << fFocLgt << endl;
  *oPrtStrm << "        fPlateScaleFactor " << fPlateScaleFactor << endl;
  *oPrtStrm << "        fAvgTransitTime " << fAvgTransitTime << endl;
  *oPrtStrm << "        fRotationOffset " << fRotationOffset << endl;
  *oPrtStrm << "        Primary   fDp fDpinner fZp " << fDp << " " 
	    << fDpinner << " " << fZp << endl;
  *oPrtStrm << "        Primary translation offsetx/y/z " << fPrimXoffset
	    << "  " << fPrimYoffset << "  " << fPrimZoffset << endl;
  *oPrtStrm << "        Primary rotation offset theta phi " 
	    << fPrimThetaOffset << "  " << fPrimPhiAngle << endl;
  *oPrtStrm << "        Primary Roughness Sigma " << fPrimRoughSigma << endl;
  //*oPrtStrm << "        Primary Roughness Max " << fPrimRoughMax << endl;
  *oPrtStrm << "        Secondary fDs fDsinner fZs " << fDs << " "
	    << fDsinner << " " << fZs << endl;
  *oPrtStrm << "        Secondary translation offsetx/y/z " << fSecondXoffset
	    << "  " << fSecondYoffset << "  " << fSecondZoffset << endl;
  *oPrtStrm << "        Secondary rotation offset theta phi " 
	    << fSecondThetaOffset << "  " << fSecondPhiAngle << endl;
  *oPrtStrm << "        Secondary Roughness Sigma " << fSecondRoughSigma 
	    << endl;
  //*oPrtStrm << "        Secondary Roughness Max " << fSecondRoughMax << endl;
  *oPrtStrm << "        Focal Surface fk1 fk2 fZf " << fk1 << " " 
	    << fk2 << " " << fZf << endl;
  *oPrtStrm << "        Focal Surfac2 translation offsetx/y/z " << fFocalPlXoffset
	    << "  " << fFocalPlYoffset << "  " << fFocalPlZoffset << endl;
  *oPrtStrm << "        Focal Surface rotation offset theta phi " 
	    << fFocalPlThetaOffset << "  " << fFocalPlPhiAngle << endl;
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
  *oPrtStrm << "        fMAPMTOffset        " << fMAPMTOffset << endl;
  *oPrtStrm << "        fMAPMTGap        " << fMAPMTGap << endl;
  *oPrtStrm << "        fMAPMTRefIndex        " << fMAPMTRefIndex << endl;
  if (gPrimRefl) {
    *oPrtStrm << "      Primary reflection coefficients " << endl;
    Int_t numP = gPrimRefl->GetN();
    for (Int_t i =0;i<numP;i++) {
      double x=0.0; 
      double y=0.0;
      gPrimRefl->GetPoint(i,x,y);
      *oPrtStrm << "        " << i << "     " << x/nm 
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
      *oPrtStrm << "        " << i << "     " << x/nm 
		<< "    " << y << endl;
    }
  }
  else {
    *oPrtStrm << "      gSeconRefl TGraph pointer is zero " << endl;      
  }
  
  *oPrtStrm << endl;
};
/********************** end of printTelescope *****************/

void GSCTelescope::drawTelescope(const int &option) {
  (void) option; // unused parameter
  
  bool debug = false;
  if (debug) {
    *oLog << " ready to draw the SC telescope" << endl; 
  }
  gGeoManager = manager;
  gGeoManager->GetMasterVolume()->Draw("ogl");
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
  hisT->Branch("yLast",&fLocLast[1],"yLast/D");
  hisT->Branch("zLast",&fLocLast[2],"zLast/D");
  hisT->Branch("xLastDir",&fDirLast[0],"xLastDir/D");
  hisT->Branch("yLastDir",&fDirLast[1],"yLastDir/D");
  hisT->Branch("zLastDir",&fDirLast[2],"zLastDir/D");
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

/*
VERSION3.1
2March2015
*/
/*!  GTelArray.cpp
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
#include <deque>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <limits>

using namespace std;

#include "TFile.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Rotation3D.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TAxis.h"
#include "TPad.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"
#include "TLatex.h"
#include "TTree.h"
#include "TBranch.h"

#include "GDefinition.h"

#include "GUtilityFuncts.h"
#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GSCTelescope.h"
#include "GTelescopeFactory.h"
#include "GDCTelescopeFactory.h"
#include "GReadDCStdBase.h"
#include "GReadDCStdGrISU.h"
#include "GReadPhotonBase.h"
#include "GReadPhotonGrISU.h"
#include "GArrayTel.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

GArrayTel::GArrayTel() {
  int debug = false;
  if (debug) {
    *oLog << "  --GArrayTel default constructor" << endl;
  }

  tel = 0;
  fDelay = 0;
  
};
/**************end of GArrayTel ***************************/


GArrayTel::GArrayTel(const ROOT::Math::XYZVector telLocGrd,
                     const double &telOffsetX,
                     const double &telOffsetY,
                     const TelType &teltype,
		     const int &telid,
                     const int &telstd,
		     const int &printMode,
		     GTelescope *ctel )
  
  : telLocGrdGC(telLocGrd),fpointingOffsetX(telOffsetX),
    fpointingOffsetY(telOffsetY), telType(teltype),
    telID(telid), telStd(telstd), iPrintMode(printMode),
    tel(ctel) {

  fDelay = 0;
  fWobbleN = 0.0;
  fWobbleE = 0.0;
  fLatitude = 0.0;
  fAzPrim = 0.0;
  fZnPrim = 0.0;
  fEnergy = 0.0;
  fAzTel = 0.0;
  fZnTel = 0.0;
  fAzPhotGr = 0.0;
  fZnPhotGr = 0.0;
  fPhotHgtEmiss = 0.0;
  fPhotGrdTime = 0.0;
  fPhotWaveLgt = 0.0;
  iPhotType = 0;
  iPhotTelHitNum = 0;
  fTotalPhotonTime = 0.0;

  dStoreLoc = 0;
  dStorePix = 0;

  // initialize rotation matrix
  // will use later in setting up rotation matrix for telCoors.
  //tel->getPointingOffset(&fpointingOffsetX,&fpointingOffsetX);
  tel->setPrintMode(*oLog,iPrintMode);
  
  tel->printTelescope();
};
/**************end of GArrayTel ***************************/

GArrayTel::~GArrayTel() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GArrayTel::~GArrayTel: telID = " 
	  << telID << endl;
  } 
  if (tel) SafeDelete(tel);
  if (dStoreLoc) SafeDelete(dStoreLoc);
  if (dStorePix) SafeDelete(dStorePix);

};
/**************end of GArrayTel ***************************/

void GArrayTel::setPrimary(const ROOT::Math::XYZVector &vSCorec,
                           const ROOT::Math::XYZVector &vSDcosGdc,
                           const double &AzPrimc,const double &ZnPrimc,
                           const double &energyc,const double &WobbleTNc,
                           const double &WobbleTEc,
                           const double &Latitude) {
  bool debug = false;
  
  vSCoreGC = vSCorec;
  vSDcosGC = vSDcosGdc;
  fAzPrim = AzPrimc;
  fZnPrim = ZnPrimc;
  fEnergy = energyc;
  fWobbleN = WobbleTNc;
  fWobbleE = WobbleTEc;
  fLatitude = Latitude;
  
  if (debug) {
    *oLog << "  -- GArrayTel::setPrimary: telID = " << telID << endl;
    *oLog << "        core loc.Grd.Coor:  ";
    GUtilityFuncts::printGenVector(vSCoreGC); *oLog << endl;
    *oLog << "        core Az Zn. GrdCoor:  ";
    GUtilityFuncts::printGenVector(vSDcosGC); *oLog << endl;
    *oLog << "        wobbleN wobbleE  " 
          << fWobbleN*(TMath::RadToDeg()) << "  " 
          << fWobbleE*(TMath::RadToDeg()) << endl;
    *oLog << "        primary energy  " << fEnergy << endl;
    *oLog << "        fZnPrim   " << fZnPrim << endl;
    *oLog << "        telLocGrdGC.Z() " << telLocGrdGC.Z() << endl;

    *oLog << "              parameters for telescopeAzZnRot " << endl;
    *oLog << "        fAzPrim, fZnPrim   " << fAzPrim << " " << fZnPrim << endl;
    *oLog << "        fWobbleN, fWobbleE " << fWobbleN << " " << fWobbleE << endl;
    *oLog << "        offsetX, offsetY   " << fpointingOffsetX << " "
	  << fpointingOffsetY << endl;
    *oLog << "        fLatitude " << fLatitude << endl;
  }
  
  // get telescope az and zn
  GUtilityFuncts::telescopeAzZnRot(fAzPrim,
                                   fZnPrim,
                                   fWobbleN,
                                   fWobbleE,
                                   fpointingOffsetX,
                                   fpointingOffsetY,
                                   fLatitude,
                                   &fAzTel,&fZnTel);
  // get telescope dir.cosines
  double xcos = 0.0;
  double ycos = 0.0;
  double zcos = 0.0;

  GUtilityFuncts::AzZnToXYcos(fAzTel,fZnTel,
                   &xcos,&ycos);
  zcos = sqrt(1-xcos*xcos - ycos*ycos);
  vTelDcosGC.SetCoordinates(xcos,ycos,zcos); 
      
  if (debug) {
    *oLog << "        telescope Az Zn " 
          << fAzTel*(TMath::RadToDeg()) << "  " 
          << "  " << fZnTel*(TMath::RadToDeg()) << endl;
    *oLog << "        telescope DCos ";
    GUtilityFuncts::printGenVector(vTelDcosGC); *oLog << endl;
  }

  // get ground to telescope rotation matrix
  GUtilityFuncts::AzZnToRotMat(fAzTel,fZnTel,&rotGrdToTel);

  // get telescope grd.location in telescope coordinates
  telLocGrdTC = rotGrdToTel*telLocGrdGC;
  
  // get delay time for this telescope
  fDelay = telLocGrdTC.Z();  // convert to ns later

  // test rotation matrix, rotate telescope
  if (debug) {
    ROOT::Math::XYZVector vTDCosT = rotGrdToTel*vTelDcosGC;
    GUtilityFuncts::zeroFloatVectorFix(&vTDCosT);
    *oLog << "        telescope dCos TelSys ";
    GUtilityFuncts::printGenVector(vTDCosT); *oLog << endl;
    *oLog << "           should be 0 0 1 " << endl;
    *oLog << "        telescope location grd.coor  ";
    GUtilityFuncts::printGenVector(telLocGrdGC); *oLog << endl;
    *oLog << "        telescope location tel.coor  ";
    GUtilityFuncts::printGenVector(telLocGrdTC); *oLog << endl;
    *oLog << "        fDelay(meters): " << fDelay << endl;
  }

};
/**************end of setPrimary ***************************/
  
void GArrayTel::printArrayTel() {
  *oLog << "     GArrayTel::printArrayTel " << endl;
  *oLog << "        telID  telStd / TelType / location / telOffsetX/Y "
	<< endl;
  *oLog << "          " << telID << "       " << telStd << " / " << getTelType(telType)
        << " / " << telLocGrdGC.X() << " " << telLocGrdGC.Y() 
        << " " << telLocGrdGC.Z() << " /  " 
        << fpointingOffsetX << "  " << fpointingOffsetY << endl;

};
/**************end of printGArrayTel ***************************/

void GArrayTel::setPhoton(const ROOT::Math::XYZVector &pGrd,
                          const ROOT::Math::XYZVector &pDcos,
                          const double &pAz,const double &pZn,
                          const double &pHgtEmiss,const double &pTime,
                          const double &pWaveLgt, const int &pType,
                          const int &pTel) {
  
  
  bool debug = false;
  if (telID==1) debug = false;
  bool debugShort = false;
  if ( (debug) || (debugShort) ){
    *oLog << "  -- GArrayTel::setPhoton: telID: " << telID << endl;
  }

  vPhotonGrdLocGC  = pGrd;
  vPhotonDcosGC  = pDcos;
  fAzPhotGr      = pAz;
  fZnPhotGr      = pZn;
  fPhotHgtEmiss  = pHgtEmiss;
  fPhotGrdTime   = pTime;
  fPhotWaveLgt   = pWaveLgt;
  iPhotType      = pType;
  iPhotTelHitNum = pTel;

  // rotate photon Dcos and grd.hit to telescope coor.
  vPhotonGrdLocTC  = rotGrdToTel*vPhotonGrdLocGC;
  vPhotonDcosTC = rotGrdToTel*vPhotonDcosGC; 
 
  vPhotonRelLocTC = vPhotonGrdLocTC - telLocGrdTC;

  if (debug) {
    *oLog << "         vPhotonGrdLoc ";       
    GUtilityFuncts::printGenVector(vPhotonGrdLocGC); *oLog << endl;    
    *oLog << "         vPhotonDcosGC ";       
    GUtilityFuncts::printGenVector(vPhotonDcosGC); *oLog << endl; 
    *oLog << "         az / zn:       " << fAzPhotGr*(TMath::RadToDeg())
          << "  /  "<< fZnPhotGr*(TMath::RadToDeg()) << endl; 
    *oLog << "         hgtEmiss  arrivTime: " << fPhotHgtEmiss
          << "  " << fPhotGrdTime << endl;
    *oLog << "         wavelgt  partType: " << fPhotWaveLgt 
          << "  " << iPhotType << endl;
    *oLog << "         Telescope hit number: " << iPhotTelHitNum << endl
          << endl;     
    *oLog << "         vPhotonGrdLocTC ";       
    GUtilityFuncts::printGenVector(vPhotonGrdLocTC); *oLog << endl; 
    *oLog << "         vPhotonDcosTC ";       
    GUtilityFuncts::printGenVector(vPhotonDcosTC); *oLog << endl;
    *oLog << "         projection of photon on tangent plane x/y(deg) " 
          << vPhotonDcosTC.X()*(TMath::RadToDeg()) << "  " 
          << vPhotonDcosTC.Y()*(TMath::RadToDeg()) << endl;
    *oLog << "         relative photon loc.TelCoor  ";       
    GUtilityFuncts::printGenVector(vPhotonRelLocTC); *oLog << endl;
    *oLog << "         telescope loc. TelCoor  ";       
     GUtilityFuncts::printGenVector(telLocGrdTC); *oLog << endl;

     *oLog << "    ----  end of debug GArrayTel photon" << endl;
    
  }
  // inject photon into telescope

  if (0) {
    *oLog << "         vPhotonGrdLoc ";       
    GUtilityFuncts::printGenVector(vPhotonGrdLocGC); *oLog << endl;    
     *oLog << "         vPhotonGrdLocTC ";       
    GUtilityFuncts::printGenVector(vPhotonGrdLocTC); *oLog << endl; 
     
    *oLog << "         telLocGrdTC  ";       
     GUtilityFuncts::printGenVector(telLocGrdTC); *oLog << endl;

    *oLog << "         vPhotonRelLocTC ";
    GUtilityFuncts::printGenVector(vPhotonRelLocTC); *oLog << endl << endl;;
  }    

  tel->injectPhoton(vPhotonRelLocTC,vPhotonDcosTC,fPhotWaveLgt);

  double netTelescopeTime;
  bOnCamera = tel->getCameraPhotonLocation(&vPhotonCameraLoc,
					   &vPhotonCameraDcos,
					   &netTelescopeTime);

  fTotalPhotonTime = netTelescopeTime + fPhotGrdTime;

  // get the ray-tracing results from tel and store in deque.
  // could have a get photon for this class as well.
  if (debugShort) {
    *oLog << endl;
    *oLog << "          bOnCamera  " << bOnCamera << endl;
    *oLog << "       GArrayTel:  from GTelescope::getCameraPhotonLocation ";
    GUtilityFuncts::printGenVector(vPhotonCameraLoc);  *oLog << endl;
    *oLog << "          vPhotonCameraDcos ";
    GUtilityFuncts::printGenVector(vPhotonCameraDcos); *oLog << endl;
    *oLog << "          netTelescopeTime " << netTelescopeTime << endl;
    *oLog << "          fPhotGrdTime     " << fPhotGrdTime << endl;
    *oLog << "          fTotalPhotonTime " << fTotalPhotonTime << endl;
  }

};  
/************** end of setPhoton ***************************/

bool GArrayTel::getCameraPhotonLocation(ROOT::Math::XYZVector *photonLoc,
					ROOT::Math::XYZVector *photonDcos,
					double *photonTime) {

  bool debug = false;
  if (debug && bOnCamera) {
    *oLog << "  -- GArrayTel::getCameraPhotonLocation " << endl;
  }
  *photonLoc = vPhotonCameraLoc;
  *photonDcos = vPhotonCameraDcos;
  *photonTime = fTotalPhotonTime;
  if (debug && bOnCamera) {
    *oLog << "     telID fTotalPhotonTime " << telID << " "
	  << fTotalPhotonTime << endl;  
  }
  return bOnCamera;
};
/************** end of getCameraPhotonLocation ***************************/

void GArrayTel::makeTelescopeTest(const string& testfile) {
  
  string baseName = testfile;

  if (baseName=="") {
    *oLog << "  NO BASENAME GIVEN IN OPTICSSIMULATION TESTTEL RECORD" << endl;
    *oLog << "  setting default base name to testTel" << endl;
    baseName = "testTel";
  }

  string testTelRootFile = baseName + ".root";
  TFile *fOut = new TFile(testTelRootFile.c_str(),"recreate");
  TTree *ttree = new TTree("rmsSpot","rmsSpot");

  Double_t rmsXT = 0.0;
  Double_t rmsYT = 0.0;
  Double_t degT  = 0.0;
  ttree->Branch("rmsX",&rmsXT,"rmsX/D");
  ttree->Branch("rmsY",&rmsYT,"rmsY/D");
  ttree->Branch("deg",&degT,"deg/D");

  bool debug = false;
  if (1) {
    *oLog << "  -- GArrayTel::makeTelescopeTest " << endl;
    *oLog << "      baseName " << baseName << endl;
  }

  /*
  
	  make array of histograms for x/y on camera and transit time from ground,
                     correct for time offset and for location from plate scale.

	  Could put histograms in a root file for plotting later???? and/or make plots.
          make plots first.
   */

  // make vector of directions, must come from method parameter above(make vector from min,max,delta)
  // numDegBins must not exceed 8
  int numDegBins;
  *oLog << "   telType  " << telType << endl;

  double rangeGraphTimeX = 5.0;
  double rangeGraphTimeY = 0.8;
  double rangeGraphPsfX = 5.0;
  double rangeGraphPsfY = 0.8;
  double spotHistXY = 8.0;
  double timeHistXY = 6.0;

  double deltaDeg = 0.0;
  //string baseFileName;
  string imageFileType = "png";

  // use these values for SC telescope (telType ==1)
  if (telType ==1) {
    numDegBins = 9;
    deltaDeg = 0.5;
 
    rangeGraphTimeX = 5.0;
    rangeGraphTimeY = 0.8;
    rangeGraphPsfX = 5.0;
    rangeGraphPsfY = 6.0;
    spotHistXY = 8.0;
    timeHistXY = 4.0;    
    // fill in the rest later if have time
  }

  // these values are for the veritas DC telescope
  if (telType ==0) {
    numDegBins = 5;
    deltaDeg = 0.5;

    rangeGraphTimeX = 3.0;
    rangeGraphTimeY = 10.0;
    rangeGraphPsfX = 2.5;
    rangeGraphPsfY = 10.0;
    spotHistXY = 20.0;
    timeHistXY = 10.0;
     // fill in the rest later if have time
  }
  // set number of photons
  int nPhotons = 25000;

  vector<double> vDeg;
  double setdeg = 0.0;
  for (int i = 0; i < numDegBins; i++) {
    vDeg.push_back(setdeg);
    setdeg+=deltaDeg;
  }

  // create hist point arrays
  const int kN = vDeg.size();

  //TH1D* histT[kN];   // time histograms
  TGraph* graRMS = new TGraph;
  graRMS->SetTitle(";Field angle (deg);2 #times max{RMS_{sagital}, RMS_{tangential}} (arcmin)");

  TGraph* graT = new TGraph;
  graT->SetTitle(";Field angle (deg);Photon propagation time spread (RMS) (ns)");
  double spot = spotHistXY;
  for(Int_t n = 0; n < kN; n++){
    
    double deg = vDeg[n];

    vHist.push_back(new TH2D(Form("hist%d", n), Form("#it{#theta} = %4.1f (deg);X (arcmin);Y (arcmin)", deg), 1000, -spot, spot, 1000, -spot, spot) );
 
    vHistT.push_back(new TH1D(Form("histT%d",n), Form("#it{#theta} = %4.1f (deg);Propagation delay (ns);Entries", deg), 120, -timeHistXY, timeHistXY) );
  }
  // get estimated transit time
  double idealTransitTime = tel->getIdealTransitTime();
  *oLog << "    estimated ideal transit time, normal incidence  " << idealTransitTime  << endl;

  // get focal length of telescope
  double fLgt = tel->getFocalLength();
  *oLog << " focal length " << fLgt << endl;
  // get array of photon locations on circle on ground plane
  double fRadius = tel->getTelescopeRadius();
  fRadius = fRadius*1.10;  // increase by 10 percent
  if (1) {
    *oLog << "         telescope radius " << fRadius << endl;
    *oLog << "         number of photons " << nPhotons << endl;
  }
  
  // vectors for incoming photons
  vector<double> vX;
  vector<double> vY;
  vector<double> vZ;
  vX.reserve(nPhotons);
  vY.reserve(nPhotons);
  vZ.reserve(nPhotons);
  
  for (int j=0;j<nPhotons;j++) {
    double r = fRadius*sqrt(TR3.Rndm());
    double phi =  TR3.Rndm()*(TMath::TwoPi());
    double x  = r*cos(phi);
    double y  = r*sin(phi);
    double z  = 0.0;
    vX.push_back(x);
    vY.push_back(y);
    vZ.push_back(z);
  }
  if (debug) {
     *oLog <<      " vX.size() " << vX.size() << endl;
     for (int j=0;j<nPhotons;j++) {
       *oLog <<      "  vX  vY  vZ  " << vX.at(j) << " " << vY.at(j) << " " << vZ.at(j) << endl;
    }
  }

  // set direction and wavelength of incoming photons
  double lambda = 400.0;

  for (unsigned ideg=0;ideg<vDeg.size();ideg++) {
    // set directions
    double theta = vDeg[ideg]*(TMath::DegToRad());
    if (1) {
      *oLog << endl << "           start of ideg loop for deg = " << vDeg[ideg] << endl;
    }
    double dl = sin(theta);
    double dm = 0.0;
    double dn = -cos(theta);
    ROOT::Math::XYZVector photonDir(dl,dm,dn);
    
    for (int iph = 0;iph<nPhotons;iph++) {
      
     // set locations of photon (put inside loop
      ROOT::Math::XYZVector photonLoc(vX.at(iph),vY.at(iph),vZ.at(iph));
      
      if (0) {
	*oLog << "         photonLoc  ";       
      GUtilityFuncts::printGenVector(photonLoc); *oLog << endl;
      *oLog << "         photonDir  ";       
      GUtilityFuncts::printGenVector(photonDir); *oLog << endl;
      } 
    
      tel->injectPhoton(photonLoc,photonDir,lambda);
      
      ROOT::Math::XYZVector cameraLoc(0.0,0.0,0.0);
      ROOT::Math::XYZVector cameraDir(0.0,0.0,0.0);
      double cameraTime = 0.0;
      bool onCamera = false;
      
      onCamera = tel->getCameraPhotonLocation(&cameraLoc,&cameraDir,&cameraTime);

      if (onCamera) {
	// convert to minutes of arc, first subtract focus location since we just want spread
	double xMin = (cameraLoc.X() / (fLgt*1000.))*TMath::RadToDeg();
	double yMin = (cameraLoc.Y() / (fLgt*1000.))*TMath::RadToDeg();
	xMin = (xMin - theta*TMath::RadToDeg()) * 60.0;
	yMin = yMin*60.0;
	if (0) {
	  *oLog << "   onCamera location in minutes of arc " << xMin << "  " << yMin << endl;
	  *oLog << "         cameraDir  ";       
	  GUtilityFuncts::printGenVector(cameraDir); *oLog << endl;
	  *oLog << "         cameraTime " << cameraTime << endl;
	}
	vHist[ideg]->Fill(xMin,yMin);
	vHistT[ideg]->Fill(cameraTime - idealTransitTime);
	//*oLog << "hist fill " << ideg << " " << xMin << " " << yMin << endl;
	//histT[ideg]->Fill(cameraTime - idealTransitTime);
      }  // onCamera
    }  // iph
    Double_t rmsx = vHist[ideg]->GetRMS(1);
    Double_t rmsy = vHist[ideg]->GetRMS(2);

    rmsXT = rmsx;
    rmsYT = rmsy;
    degT = vDeg[ideg];
    ttree->Fill();
    if (1) {
      *oLog << " rmsx  rmsy " << rmsx << "  " << rmsy << endl;
      *oLog << " graRMS->GetN() " << graRMS->GetN() << endl;
    }
      graRMS->SetPoint(graRMS->GetN(),vDeg[ideg], (rmsx > rmsy ? rmsx: rmsy)*2);
      graT->SetPoint(graT->GetN(),vDeg[ideg],vHistT[ideg]->GetRMS());


  } //ideg
  //ttree->Write();

  /*
  string imageFilename;


  TCanvas* canSpot = new TCanvas("canSpot", "canSpot", 1200, 600);
  canSpot->Divide(4, 2);
  TCanvas* canTime = new TCanvas("canTime", "canTime", 1200, 600);
  canTime->Divide(4, 2);

  for(Int_t i = 0; i < kN; i++){
    canSpot->cd(i + 1);   
    vHist[i]->DrawCopy("colz");
    canTime->cd(i + 1);
    vHistT[i]->DrawCopy();
  }
  imageFilename = baseName + "_Spot." + imageFileType;
  canSpot->SaveAs(imageFilename.c_str());
  
  imageFilename = baseName + "_Time." + imageFileType;
  canTime->SaveAs(imageFilename.c_str());
  
  TCanvas* canGrPSF = new TCanvas("canGrPSF", "canGrPSF", 600, 300);
  //canFig5->Divide(1,2);
  
  //canFig5->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graRMS->SetMarkerStyle(25);
  graRMS->GetXaxis()->SetLimits(0, rangeGraphPsfX);
  graRMS->GetYaxis()->SetRangeUser(0, rangeGraphPsfY);
  graRMS->Draw("apl");

  imageFilename = baseName +"_grPSF." + imageFileType;
  canGrPSF->SaveAs(imageFilename.c_str());

  TCanvas* canGrTime = new TCanvas("canGrTime", "canGrTime", 600, 300);

  //canFig5->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  graT->Draw("apl");
  graT->SetMarkerStyle(25);
  graT->GetXaxis()->SetLimits(0, rangeGraphTimeX);
  graT->GetYaxis()->SetRangeUser(0, rangeGraphTimeY);

  imageFilename = baseName +"_GrTime." + imageFileType;
  canGrTime->SaveAs(imageFilename.c_str());
  */
  fOut->Write();
  fOut->Close();
};
/************** end of makeTelescopeTest ***************************/

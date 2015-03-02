/*
VERSION3.1
2March2015
*/
/*! \class GRootWriter
    \brief write for pe data to root files

    Revision $Id: GRootWriter.cpp,v 1.1.2.1.10.3 2010/08/27 07:04:04 prokoph Exp $

    \author Gernot Maier
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
#include <ctime>

using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRint.h"
#include "TROOT.h"
#include "Math/Vector3D.h"

#include "GDefinition.h"
#include "GRootWriter.h"

ClassImp(GRootWriter);

GRootWriter::GRootWriter( TFile *tfile,const unsigned int &iTelID,
			  const string &treeBaseName,
			  const bool &storePhotonDcos,
			  const unsigned int &iNInitEvents,
                          const bool &debugBranchesFlag) 
  :fFile(tfile),fTelID(iTelID), treeBaseName(treeBaseName),
   bStoreDcos(storePhotonDcos),bDebugBranchesFlag(debugBranchesFlag),
   iNInitReserve(iNInitEvents)
{

  
  strcpy(version,VERSN.c_str());

  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootWriter::GRootWriter " << endl;
    *oLog << "     version: " << version << endl;
  }

  fPE_photonX = 0;
  fPE_photonY = 0;
  fPE_DcosX = 0;
  fPE_DcosY = 0;
  fPE_time = 0;
  fPE_wl = 0;
  fSrcRelTelX = 0.0;
  fSrcRelTelY = 0.0;
  fSrcRelToCameraX = 0.0;
  fSrcRelToCameraY = 0.0;  

  fXcoreTC = 0.0;
  fYcoreTC = 0.0;
  fXcosTC  = 0.0;
  fYcosTC  = 0.0;
  fXcoreSC = 0.0;
  fYcoreSC = 0.0;
  fXcosSC  = 0.0;
  fYcosSC  = 0.0;
  fXTelTC  = 0.0;
  fYTelTC  = 0.0;
  fZTelTC  = 0.0;
  numPhotonX = 0.0;
  bReserveFlag = false;
    
  // data vectors
  if (iNInitReserve > 1) bReserveFlag = true;

  fPE_photonX = new std::vector< float >();
  fPE_photonY = new std::vector< float >();
  fPE_time = new std::vector< float >();
  fPE_wl = new std::vector< float >();

  if(bReserveFlag) {
    fPE_photonX->reserve( iNInitReserve );
    fPE_photonY->reserve( iNInitReserve );
    fPE_time->reserve( iNInitReserve );
    fPE_wl->reserve( iNInitReserve );

  }

  if (bStoreDcos) {
    fPE_DcosX = new std::vector< float >(); 
    fPE_DcosY = new std::vector< float >(); 
    if(bReserveFlag) {
      fPE_DcosX->reserve( iNInitReserve );
      fPE_DcosY->reserve( iNInitReserve );
    }
  }

 // define trees with pe data
  char hname[400];
  char htitle[400];
 
  sprintf(hname,"%s%d",treeBaseName.c_str(),fTelID);
  sprintf( htitle, "photon data for telescope %d", fTelID );

  if (debug) {  
    *oLog << "      making  TTree: " << hname << " " << htitle << endl;
  }
  fFile->cd();
  fTree = new TTree( hname, htitle );
  
  fTree->Branch( "eventNumber", &fEventNumber, "eventNumber/i" );
  fTree->Branch( "primaryType", &fPrimaryType, "primaryType/i" );
  fTree->Branch( "primaryEnergy", &fPrimaryEnergy, "primaryEnergy/F"  );

  //CD:2Mar2015  3 new branches
  fTree->Branch( "FirstIntHgt", &fFirstIntHgt,"FirstIntHgt/F" );
  fTree->Branch( "FirstIntDpt", &fFirstIntDpt,"FirstIntDpt/F" );
  fTree->Branch( "ShowerID",&iShowerID,"ShowerID/i" );

  fTree->Branch( "Xcore", &fXcore, "Xcore/F" );
  fTree->Branch( "Ycore", &fYcore, "Ycore/F" );
  fTree->Branch( "Xcos", &fXcos, "Xcos/F" );
  fTree->Branch( "Ycos", &fYcos, "Ycos/F" );
  fTree->Branch( "Xsource", &fSrcRelToCameraX, "Xsource/F" );
  fTree->Branch( "Ysource", &fSrcRelToCameraY, "Ysource/F" );
  fTree->Branch( "AzPrim", &fAzPrim, "AzPrim/F" );
  fTree->Branch( "ZnPrim", &fZnPrim, "ZnPrim/F" );
  fTree->Branch( "AzTel", &fAzTel, "AzTel/F" );
  fTree->Branch( "ZnTel", &fZnTel, "ZnTel/F" );

  fTree->Branch( "delay", &fDelay, "delay/F" );
  
  fTree->Branch( "photonX", &fPE_photonX );
  fTree->Branch( "photonY", &fPE_photonY );
  fTree->Branch( "time", &fPE_time );
  fTree->Branch( "wavelength", &fPE_wl );

  if (bStoreDcos) {
    fTree->Branch( "photonDcosX", &fPE_DcosX );    
    fTree->Branch( "photonDcosY", &fPE_DcosY );    
  }

  if (bDebugBranchesFlag) {
    *oLog << "making debug branches" << endl;
    fTree->Branch( "TelID", &fTelID, "TelID/i" );
    fTree->Branch( "XcoreTC", &fXcoreTC, "XcoreTC/F" );
    fTree->Branch( "YcoreTC", &fYcoreTC, "YcoreTC/F" );
    fTree->Branch( "XcosTC", &fXcosTC, "XcosTC/F" );
    fTree->Branch( "YcosTC", &fYcosTC, "YcosTC/F" );
    fTree->Branch( "XcoreSC", &fXcoreSC, "XcoreSC/F" );
    fTree->Branch( "YcoreSC", &fYcoreSC, "YcoreSC/F" );
    fTree->Branch( "XcosSC", &fXcosSC, "XcosSC/F" );
    fTree->Branch( "YcosSC", &fYcosSC, "YcosSC/F" );
    fTree->Branch( "XTelTC", &fXTelTC, "XTelTC/F" );
    fTree->Branch( "YTelTC", &fYTelTC, "YTelTC/F" );
    fTree->Branch( "ZTelTC", &fZTelTC, "ZTelTC/F" );
    
  }  
};
//******************************** end of GRootWriter ********************

GRootWriter::~GRootWriter() {
  // worry about TTree later (fTree)
  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootWriter::~GRootWriter" << endl;
  }
  SafeDelete(fPE_photonX);
  SafeDelete(fPE_photonY);
  SafeDelete(fPE_DcosX);
  SafeDelete(fPE_DcosY);
  SafeDelete(fPE_time);
  SafeDelete(fPE_wl);

}
//******************************** end of ~GRootWriter ********************

int GRootWriter::addEvent(const unsigned int &eventNumber, const unsigned int &primaryType, 
			  const double &primaryEnergy,const ROOT::Math::XYZVector &vSCore,
			  const ROOT::Math::XYZVector &vSDCore, const double &xSource,
			  const double &ySource,const double &delayTime, 
			  const double &transitTime, const double &azTel,const double &znTel,
                          const double &azPrim, const double &znPrim,
                          const double &srcX,const double &srcY,
                          const double &firstIntHgt, const double &firstIntDpt,
                          const unsigned int &showerID,

                          const ROOT::Math::XYZVector &vSCoreTC,
                          const ROOT::Math::XYZVector &vSDcosTC,
                          const ROOT::Math::XYZVector &vSCoreSC,
                          const ROOT::Math::XYZVector &vSDcosSC,
                          
                          const ROOT::Math::XYZVector &vTelLocTC
                          ) {


  bool debug1 = false;
  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootWriter::addEvent; telnumber  " << fTelID << endl;
  }

  if( !fTree ) {
    if (debug) {
      *oLog << "         NO TREE, fTree is 0: returning " << endl;
    }
    return 0;
  }
  
  fEventNumber = eventNumber;
  fPrimaryType = primaryType;
  fPrimaryEnergy = (float)primaryEnergy;
  fXcore = (float)vSCore.X();
  fYcore = (float)vSCore.Y();
  fXcos  = (float)vSDCore.X();
  fYcos  = (float)vSDCore.Y();
  fXsource = (float)xSource*(TMath::RadToDeg());
  fYsource = (float)ySource*(TMath::RadToDeg());
  fDelay = (float)delayTime;
  fTransit = (float)transitTime;
  fAzPrim = (float)azPrim*(TMath::RadToDeg());
  fZnPrim = (float)znPrim*(TMath::RadToDeg());
  fAzTel = (float)azTel*(TMath::RadToDeg());
  fZnTel = (float)znTel*(TMath::RadToDeg());
  fSrcRelToCameraX = (float)srcX*(TMath::RadToDeg());
  fSrcRelToCameraY = (float)srcY*(TMath::RadToDeg());

  fFirstIntHgt = (float)firstIntHgt;
  fFirstIntDpt = (float)firstIntDpt;
  iShowerID    = showerID;
    
  // debug branches
  fXcoreTC = (float)vSCoreTC.X();
  fYcoreTC = (float)vSCoreTC.Y();
  fXcosTC = (float)vSDcosTC.X();
  fYcosTC = (float)vSDcosTC.Y();
  fXcoreSC = (float)vSCoreSC.X();
  fYcoreSC = (float)vSCoreSC.Y();
  fXcosSC = (float)vSDcosSC.X();
  fYcosSC = (float)vSDcosSC.Y();
  fXTelTC = (float)vTelLocTC.X();
  fYTelTC = (float)vTelLocTC.Y();
  fZTelTC = (float)vTelLocTC.Z();

  fFile->cd();

  numPhotonX = fPE_photonX->size();

  int r = fTree->Fill();
  if (debug ) {
    *oLog << "       r after tree fill " << r << endl;
    *oLog << "       photon vector sizes " << endl;
    *oLog << "       fPE_photonX->size():  " << fPE_photonX->size() << endl;
    *oLog << "       fPE_photonY->size():  " << fPE_photonY->size() << endl;
    *oLog << "       fPE_time->size():  " << fPE_time->size() << endl;
    *oLog << "       fPE_wl->size():  " << fPE_wl->size() << endl;
  }

  if (debug1) {
    
    *oLog << "telid photonX.size " << fTelID << " " << fPE_photonX->size() << endl;
  }

  fPE_photonX->clear();
  fPE_photonY->clear();
  fPE_time->clear();
  fPE_wl->clear();

  if (bStoreDcos) {
    fPE_DcosX->clear();
    fPE_DcosY->clear();
  }

  // set vector capacities if necessary
  if(bReserveFlag) {
    if ( (unsigned)fPE_photonX->capacity() > iNInitReserve ) {
      vector<float>().swap(*fPE_photonX);
      fPE_photonX->reserve( iNInitReserve );

      vector<float>().swap(*fPE_photonY);
      fPE_photonY->reserve( iNInitReserve );

      vector<float>().swap(*fPE_time);
      fPE_time->reserve( iNInitReserve );

      vector<float>().swap(*fPE_wl);
      fPE_wl->reserve( iNInitReserve );
    
      if (bStoreDcos) {
       vector<float>().swap(*fPE_DcosX);
       fPE_DcosX->reserve( iNInitReserve );

       vector<float>().swap(*fPE_DcosY);
       fPE_DcosY->reserve( iNInitReserve );
        
      }
    }
  }
    
    return r;
};
//******************************** end of add_event **********************************

//bool GRootWriter::sortPair(const pair<int,int> i , const pair<int,int> j) {
//bool test = (i.second < j.second);
//return test;
//};

//******************************** end of add_event **********************************
void GRootWriter::addPhoton(const ROOT::Math::XYZVector &PhotonCameraLoc,
			    const ROOT::Math::XYZVector &PhotonCameraDcos,
			    const double &iPE_time, 
			    const double &iPE_wl) {
 
  //return;
  float phox = (float) PhotonCameraLoc.X();
  // camera y-axis is opposite to telescope y-axis.
  float phoy = -(float)PhotonCameraLoc.Y();
  fPE_photonX->push_back( phox );
  fPE_photonY->push_back( phoy );
  float time = (float)iPE_time;
  fPE_time->push_back( time );
  float wl = (float) iPE_wl;
  fPE_wl->push_back( wl );

  if (bStoreDcos) {
    float phoDcosX = (float)PhotonCameraDcos.X();
    // camera y-axis is opposite to telescope y-axis
    float phoDcosY = (float)PhotonCameraDcos.Y();
    fPE_DcosX->push_back( phoDcosX);
    fPE_DcosY->push_back( phoDcosY);
  }
};
//******************************** end of add_pe **********************************

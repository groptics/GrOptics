/*
VERSION2.2
10MAY2012
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
			  const unsigned int &iNInitEvents) 
  :fFile(tfile),fTelID(iTelID), treeBaseName(treeBaseName),
   bStoreDcos(storePhotonDcos) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootWriter::GRootWriter " << endl;
  }
  fPE_photonX = 0;
  fPE_photonY = 0;
  fPE_DcosX = 0;
  fPE_DcosY = 0;
  fPE_time = 0;
  fPE_wl = 0;
  
  //fFile = tfile;
  //fTelID = iTelID;
  //bStoreDcos = storePhotonDcos;
  
  // data vectors
  fPE_photonX = new std::vector< float >();
  if( iNInitEvents < fPE_photonX->max_size() ) fPE_photonX->reserve( iNInitEvents );
  else                                         fPE_photonX->reserve( fPE_photonX->max_size() - 1 );
  fPE_photonY = new std::vector< float >();
  
  if( iNInitEvents < fPE_photonY->max_size() ) fPE_photonY->reserve( iNInitEvents );
  else                                         fPE_photonY->reserve( fPE_photonY->max_size() - 1 );
  
  fPE_time = new std::vector< float >();
  if( iNInitEvents < fPE_time->max_size() ) fPE_time->reserve( iNInitEvents );
  else                                      fPE_time->reserve( fPE_time->max_size() - 1 );
  
  fPE_wl = new std::vector< float >();
  if( iNInitEvents < fPE_wl->max_size() ) fPE_wl->reserve( iNInitEvents );
  else                                    fPE_wl->reserve( fPE_wl->max_size() - 1 );
  
  if (bStoreDcos) {
    fPE_DcosX = new std::vector< float >(); 
    if( iNInitEvents < fPE_DcosX->max_size() ) fPE_DcosX->reserve( iNInitEvents );
    else                                      fPE_DcosX->reserve( fPE_DcosX->max_size() - 1 );
    fPE_DcosY = new std::vector< float >(); 
    if( iNInitEvents < fPE_DcosY->max_size() ) fPE_DcosY->reserve( iNInitEvents );
    else                                      fPE_DcosY->reserve( fPE_DcosY->max_size() - 1 );
  } 

 // define trees with pe data
  char hname[400];
  char htitle[400];
 
  sprintf(hname,"%s%d",treeBaseName.c_str(),fTelID);
  //sprintf( hname, "tPH_T%d", fTelID );
  sprintf( htitle, "photon data for telescope %d", fTelID );

  if (debug) {  
    *oLog << "      making  TTree: " << hname << " " << htitle << endl;
  }
  fFile->cd();
  fTree = new TTree( hname, htitle );
  
  fTree->Branch( "eventNumber", &fEventNumber, "eventNumber/i" );
  fTree->Branch( "primaryType", &fPrimaryType, "primaryType/i" );
  fTree->Branch( "primaryEnergy", &fPrimaryEnergy, "primaryEnergy/F"  );
  fTree->Branch( "Xcore", &fXcore, "Xcore/F" );
  fTree->Branch( "Ycore", &fYcore, "Ycore/F" );
  fTree->Branch( "Xcos", &fXcos, "Xcos/F" );
  fTree->Branch( "Ycos", &fYcos, "Ycos/F" );
  fTree->Branch( "Xsource", &fXsource, "Xsource/F" );
  fTree->Branch( "Ysource", &fYsource, "Ysource/F" );
  fTree->Branch( "delay", &fDelay, "delay/F" );
  //fTree->Branch( "transit", &fTransit, "transit/F" );
  
  fTree->Branch( "photonX", &fPE_photonX );
  fTree->Branch( "photonY", &fPE_photonY );
  fTree->Branch( "time", &fPE_time );
  fTree->Branch( "wavelength", &fPE_wl );

  if (bStoreDcos) {
    fTree->Branch( "photonDcosX", &fPE_DcosX );    
    fTree->Branch( "photonDcosY", &fPE_DcosY );    
  }
  
};
//******************************** end of GRootWriter ********************

GRootWriter::~GRootWriter() {
  // worry about TTree later (fTree)
  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootWriter::~GRootWriter" << endl;
  }
  delete fPE_photonX;
  delete fPE_photonY;
  delete fPE_DcosX;
  delete fPE_DcosY;
  delete fPE_time;
  delete fPE_wl;

}
//******************************** end of ~GRootWriter ********************

int GRootWriter::addEvent(const unsigned int &eventNumber, const unsigned int &primaryType, 
			  const double &primaryEnergy,const ROOT::Math::XYZVector &vSCore,
			  const ROOT::Math::XYZVector &vSDCore, const double &xSource,
			  const double &ySource,const double &delayTime, 
			  const double &transitTime) {


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
  fXsource = (float)xSource;
  fYsource = (float)ySource;
  fDelay = (float)delayTime;
  fTransit = (float)transitTime;
  
  fFile->cd();
  int r = fTree->Fill();
  if (debug ) {
    *oLog << "       r after tree fill " << r << endl;
    *oLog << "       photon vector sizes " << endl;
    *oLog << "       fPE_photonX->size():  " << fPE_photonX->size() << endl;
    *oLog << "       fPE_photonY->size():  " << fPE_photonY->size() << endl;
    *oLog << "       fPE_time->size():  " << fPE_time->size() << endl;
    *oLog << "       fPE_wl->size():  " << fPE_wl->size() << endl;
  }

  fPE_photonX->clear();
  fPE_photonY->clear();
  fPE_time->clear();
  fPE_wl->clear();

  if (bStoreDcos) {
    fPE_DcosX->clear();
    fPE_DcosY->clear();
  }
  
  return r;
};
//******************************** end of add_event **********************************
void GRootWriter::addPhoton(const ROOT::Math::XYZVector &PhotonCameraLoc,
			    const ROOT::Math::XYZVector &PhotonCameraDcos,
			    const double &iPE_time, 
			    const double &iPE_wl) {
 
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

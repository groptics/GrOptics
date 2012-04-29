/*
VERSION2.1
1MARCH2012
*/
//! GROOTWRITER writer for photon data to root files
// adopted from VPEWriter in EVNDISP package
// GRootWriter written by Gernot Meier
// 

#ifndef GROOTWRITER
#define GROOTWRITER

// forward declarations (use include files for CINT
class TTree;
class TFile;

#include "TTree.h"
#include "TFile.h"

#include "Math/Vector3Dfwd.h"

class GRootWriter
{
 private:
   TFile *fFile;

   unsigned int fTelID;

   TTree *fTree;
   string treeBaseName;

   unsigned int fEventNumber;
   unsigned int fPrimaryType;
   float        fPrimaryEnergy;
   float        fXcore;
   float        fYcore;
   float        fXcos;
   float        fYcos;
   float        fXsource;
   float        fYsource;
   float        fDelay;
   float        fTransit;

   bool bStoreDcos;

   std::vector< float > *fPE_photonX;
   std::vector< float > *fPE_photonY;
   std::vector< float > *fPE_DcosX;
   std::vector< float > *fPE_DcosY;
   std::vector< float > *fPE_time;
   std::vector< float > *fPE_wl;

 public:

   GRootWriter( TFile *tfile, const unsigned int &iTelID, 
		const string &treeBaseName, const bool &storePhotonDcos = false, 
		const unsigned int &iNInitEvents = 100000 );

   ~GRootWriter();

  void addPhoton(const ROOT::Math::XYZVector &PhotonCameraLoc,
		  const ROOT::Math::XYZVector &PhotonCameraDcos,
		  const double &iPE_time, 
		  const double &iPE_wl = 0. );


   int addEvent(const unsigned int &eventNumber, const unsigned int &primaryType, 
		const double &primaryEnergy,const ROOT::Math::XYZVector &vSCore,
		const ROOT::Math::XYZVector &vSDCore, const double &xSource,
		const double &ySource,const double &delayTime, const double &transitTime);

   TTree*  getDataTree() { return fTree; };

   void write() {
     fFile->cd();
     fTree->Write(); 
   };

   unsigned int getTelID () { 
     return fTelID;
   };

   void cdToWriteRootFile() {
     fFile->cd();
   }
  ClassDef(GRootWriter,1);
};

#endif


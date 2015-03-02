/*
VERSION3.1
2March2015
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

   char version[20];
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
 
   float        fAzPrim;
   float        fZnPrim;
   float        fAzTel;
   float        fZnTel;
   float        fSrcRelTelX;
   float        fSrcRelTelY;
   float        fSrcRelToCameraX;
   float        fSrcRelToCameraY;

  //CD:2Mar2015  new parameters
  float fFirstIntHgt;
  float fFirstIntDpt;
  unsigned int iShowerID;

   // debug branch variables
   float        fXcoreTC;
   float        fYcoreTC;
   float        fXcosTC;
   float        fYcosTC;
   float        fXcoreSC;
   float        fYcoreSC;
   float        fXcosSC;
   float        fYcosSC;
   float        fXTelTC;
   float        fYTelTC;
   float        fZTelTC;


   bool bStoreDcos;
   bool bDebugBranchesFlag;

   std::vector< float > *fPE_photonX;
   std::vector< float > *fPE_photonY;
   std::vector< float > *fPE_DcosX;
   std::vector< float > *fPE_DcosY;
   std::vector< float > *fPE_time;
   std::vector< float > *fPE_wl;
   int iNInitReserve;
   bool bReserveFlag;

   unsigned numPhotonX;

 public:

   GRootWriter( TFile *tfile, const unsigned int &iTelID, 
		const string &treeBaseName, const bool &storePhotonDcos = false, 
		const unsigned int &iNInitEvents = 100000, 
                const bool &debugBranchesFlag = false);

   ~GRootWriter();

   void addPhoton(const ROOT::Math::XYZVector &PhotonCameraLoc,
		  const ROOT::Math::XYZVector &PhotonCameraDcos,
		  const double &iPE_time, 
		  const double &iPE_wl = 0. );


   int addEvent(const unsigned int &eventNumber, const unsigned int &primaryType, 
		const double &primaryEnergy,const ROOT::Math::XYZVector &vSCore,
		const ROOT::Math::XYZVector &vSDCore, const double &xSource,
		const double &ySource,const double &delayTime, const double &transitTime,
                const double &azTel,const double &znTel,
                const double &azPrim, const double &znPrim,
                const double &srcX,const double &srcY,
                const double &firstIntHgt, const double &firstIntDpt,
                const unsigned int &showerID,
                
                const ROOT::Math::XYZVector &vSCoreTC,
                const ROOT::Math::XYZVector &vSDcosTC,
                const ROOT::Math::XYZVector &vSCoreSC,
                const ROOT::Math::XYZVector &vSDcosSC,
                const ROOT::Math::XYZVector &vTelLocTC
                );

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
   
   unsigned getPhotonXSize() {
     return numPhotonX;
   }
   
   ClassDef(GRootWriter,1);
};

#endif


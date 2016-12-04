// compare output variables from two photonLocation trees created by GrOptics.  This is useful in
// testing code changes that should not change the output variables.
//   code somewhat based on treefriend.C in root tutorials directory.
//       Charlie Duke, Grinnell College
//       December 2, 2016

// The default names of the input files are:
//    photonLocationOrig.root and photonLocation.root.
// The default name of the tree in each of these files is T1.

// execute as a root script, i.e. "root comparePhotonLocationTrees.C" 

// The .cfg and .pilot files used to produce the photonLocationOrig.root file are in the
// GrOptics/Test/Config directory.  Thus, these files may be copied to the GrOptics/Config
// directory for use in creating a photonLocation.root file for comparison with the
// photonLocationOrig.root output file (both must use identical input files).

#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TMath.h"

#include "TCanvas.h"
#include "TStyle.h"

#include <vector>
#include <string>

// make these as global variables for now
std::vector<float> *photonX = 0;
std::vector<float> *photonY = 0;
std::vector<float> *photonDcosX = 0;
std::vector<float> *photonDcosY = 0;

void CompareTrees() {
  
  float episilon = 0.0001; // for equality tests
  std::string fileOrig = "photonLocationOrig.root";
  std::string fileTest = "photonLocation.root";
  std::string treeName = "T1";

  cout << "   original filename " << fileOrig << endl;
  cout << "   test filename     " << fileTest << endl;
  cout << "   treename          " << treeName << endl<<endl;
  
  // open the files and trees
  TFile *f = new TFile(fileOrig.c_str());
  TTree *T  = (TTree*)f->Get(treeName.c_str());
  TFile *ff = new TFile(fileTest.c_str());
  TTree *TF = (TTree*)ff->Get(treeName.c_str());
   
  std::vector<float> *fphotonX = 0;
  std::vector<float> *fphotonY = 0;
  std::vector<float> *fphotonDcosX = 0;
  std::vector<float> *fphotonDcosY = 0;
 
   T->SetBranchAddress("photonX",&photonX);
   T->SetBranchAddress("photonY",&photonY);
   T->SetBranchAddress("photonDcosX",&photonDcosX);
   T->SetBranchAddress("photonDcosY",&photonDcosY);
   
   TF->SetBranchAddress("photonX",&fphotonX);
   TF->SetBranchAddress("photonY",&fphotonY);
   TF->SetBranchAddress("photonDcosX",&fphotonDcosX);
   TF->SetBranchAddress("photonDcosY",&fphotonDcosY);
   
   T->AddFriend(TF);
   
   Long64_t nentries = T->GetEntries();
   cout << "T->GetEntries " << nentries << endl;
 
   //for (Long64_t i=0;i<nentries;i++) {
   T->GetEntry(0);
   cout << "Each vector should have the same size" << endl;
   cout << "------- entry 0 in each tree: orig tree first -------" << endl;
   cout << "  photonX->size() " << photonX->size()
	<< "   " << fphotonX->size() << endl;
   cout << "  photonY->size() " << photonY->size()
	<< "   " << fphotonY->size() << endl;
   cout << "  photonDcosX->size() " << photonDcosX->size()
	<< "   " << fphotonDcosX->size() << endl;
   cout << "  photonDcosY->size() " << photonDcosY->size()
	<< "   " << fphotonDcosY->size() << endl;

   if (photonX->size() != fphotonX->size()) {
     cout << "vector sizes are not the same " << endl;
     cout << "    exiting code " << endl;
     return;
   }
   //plot vector in T vs vector in TF 
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(2,2);
   
   c1->cd(1);
   T->Draw("T1.photonX:photonX","","", 1, 0);
   c1->cd(2);
   T->Draw("T1.photonY:photonY","","", 1, 0);
   c1->cd(3);
   T->Draw("T1.photonDcosX:photonDcosX","","", 1, 0);
   c1->cd(4);
   T->Draw("T1.photonDcosY:photonDcosY","","", 1, 0);

   // make comparison directly with epison set in code 
   Int_t numT = 30;
   cout << endl << "photonX   fphotonX       photonY  fphotonY " << endl;
   for (int i = 0;i<numT;i++) {
     cout << (*photonX)[i] << "    " << (*fphotonX)[i] << "       "
	  << (*photonY)[i] << "    " << (*fphotonY)[i] << endl;
   }
   cout << endl << "photonDcosX  fphotonDcosX   photonDcosY  fphotonDcosY " << endl;
   for (int i = 0;i<numT;i++) {
     cout << (*photonDcosX)[i] << "    " << (*fphotonDcosX)[i] << "       "
	  << (*photonDcosY)[i] << "    " << (*fphotonDcosY)[i] << endl;
   }
   cout << endl;
   cout << "detailed check for equality using episilon = " << episilon << endl;
   cout << "    number of photon locations and directions to test " << photonX->size() << endl;
   int icountNotEqualX = 0, icountNotEqualY = 0;
   int icountNotEqualCosX = 0, icountNotEqualCosY = 0;
 
   for (i=0;i<photonX->size();i++) {
     if (!TMath::AreEqualAbs( (*photonX)[i], (*fphotonX)[i],episilon) ) {
       icountNotEqualX += 1;
     }
     if (!TMath::AreEqualAbs( (*photonY)[i], (*fphotonY)[i],episilon) ) {
       icountNotEqualY += 1;
     }
     if (!TMath::AreEqualAbs( (*photonDcosX)[i], (*fphotonDcosX)[i],episilon) ) {
       icountNotEqualCosX += 1;
     }
     if (!TMath::AreEqualAbs( (*photonDcosY)[i], (*fphotonDcosY)[i],episilon) ) {
       icountNotEqualCosY += 1;
     }
   }
   
   cout << "icountNotEqualX    " << icountNotEqualX << endl;
   cout << "icountNotEqualY    " << icountNotEqualY << endl;
   cout << "icountNotEqualCosX " << icountNotEqualCosX << endl;
   cout << "icountNotEqualCosY " << icountNotEqualCosY << endl;
      
   delete f;
   delete ff;
}

void comparePhotonLocationTrees() {
   CompareTrees();
   cout << " ------------- execution complete ---------------" << endl;
}

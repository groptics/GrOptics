/* 
 VERSION4.0 
 30May2016

   root script to produce text output of the GrOptics output rootfile in
   a single file.  Do this for two different versions and then run 
   unix "diff" or "meld", etc. on these two files. If they are the same,
   then the versions produce the same output.

   Procedure: save outside of the git repository, the files in the
   Config directory.  Move to the git commit of the old version; copy
   the Config files into the Config directory, and run grOptics. Execute
   this file to produce the text output.  Move this file to a non-git
   directory.

   Do the same for the new version. Compare the two text files.

   BE SURE TO SET THE FILE NAME FOR THE OUTPUT TEXT FILE BELOW.
   THIS SCRIPT ASSUMES THERE ARE AT LEAST 4 TELESCOPES!

   Charlie Duke
   Grinnell College
   December 7, 2012
*/


{
  gROOT->Reset();
  
  cout << "  -- compareVersions " << endl;
  
  // SET FILENAMES, textFile is base name
  string rootFile = "photonLocation.root";
  string textFile = "grOpticsRootVersion2.5Beta";
  
  TFile *fr = new TFile(rootFile.c_str());
  TTree *ta[4];
  
  ta[0] =  (TTree*)fr->Get("T1");
  ta[1] =  (TTree*)fr->Get("T2");
  ta[2] =  (TTree*)fr->Get("T3");
  ta[3] =  (TTree*)fr->Get("T4");
  

  for (int i = 0;i<4;i++) {
    string scanString = "eventNumber:primaryType:primaryEnergy:Xcore:Ycore:Xcos:";
    string scanString1 = "Ycos:Xsource:Ysource:AzPrim:ZnPrim:AzTel:ZnTel:delay";
    ostringstream fi;
    fi << "testFile" << i <<".txt";
    string filename = fi.str();
    cout << filename << endl;
    scanString = scanString + scanString1;
    //cout << scanString << endl;
    string cmd = "T1->Scan(" + scanString + ");";
    //cout << cmd << endl;
    //ta[i]->Scan(scanString.c_str()); s 
    ( (TTreePlayer *)( ta[i]->GetPlayer()))->SetScanRedirect(1);
    ( (TTreePlayer *)( ta[i]->GetPlayer()))->SetScanFileName(filename.c_str());
    ta[i]->Scan("*","","",100000);
    //ta[i]->Show(0);
  }
}



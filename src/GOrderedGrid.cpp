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

#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Rotation3D.h"

#include "GUtilityFuncts.h"
#include "GDefinition.h"
#include "GOrderedGrid.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl;

bool gridBinSort(const GridFacet &s1, const GridFacet &s2) {
  return (s1.dist < s2.dist);
}


GOrderedGrid::GOrderedGrid( const vector<double> &xe,const vector<double> &ye,
                            const vector<double> &re,
                            const int &nbinsx1, const int &nbinsy1,
                            const int &option,
                            const string &filename) :vmx(xe),vmy(ye),vmr(re),
                                                     nbinsx(nbinsx1),
                                                     nbinsy(nbinsy1),
                                                     iGridOption(option),
                                                     sFileName(filename)
 {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GOrderedGrid::GOrderedGrid " << endl;
    *oLog << "       nbinsx   " << nbinsx << endl;
    *oLog << "       nbinsy   " << nbinsy << endl;
    *oLog << "       iOption  " << iGridOption << endl;
    *oLog << "       gridfilename " << sFileName << endl;
  }

  if ( (nbinsx == 0) || (nbinsy == 0) ) {
    *oLog << "     returning from GOrderedGrid constructor: " << endl; 
    *oLog << "          nbinsx == 0 or nbinsy == 0. " << endl;
    *oLog << "          full facet loop search forced " << endl;
  }
  else {

    if (iGridOption == 0) {
      return;
    }
    else if (iGridOption==1) {
      initialize();
      makeGridParameters();
      makeGrid();
    }
    else if (iGridOption==2) {
      initialize();
      makeGridParameters();
      makeGrid();
      printGrid();
    }
    else if (iGridOption==3) {
      if (!readGrid()) {
        initialize();
        makeGridParameters();
        makeGrid();
        printGrid();
      }
    }
    else {
      *oLog << "unknown iGridOption in GOrderGrid constructor: "
            << iGridOption << endl;
      exit(0);
    }
  }
};
/********************* end of GOrderedGrid ***********/

GOrderedGrid::~GOrderedGrid() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GOrderedGrid::~GOrderedGrid " << endl; 
  }

  for (unsigned i = 0;i<vGrid.size();i++) {
    SafeDelete(vGrid[i]);

  }

};
/********************* end of ~GOrderedGrid ***********/

//GOrderedGrid::GOrderedGrid(const GOrderedGrid &orderedGrid) {
//
//};
void GOrderedGrid::initialize() {
  fXmin =  99999999.9;
  fXmax = -99999999.9;
  fYmin =  99999999.9;
  fYmax = -99999999.9;

  fDelX = 0.0;
  fDelY = 0.0;

};
/********************* end of initialize ***********************/

bool GOrderedGrid::makeGridParameters() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- GOrderedGrid::makeGridParameters " << endl;
  }

  // find fXmin,fXmax,fYmin,fYmax
  for (unsigned i = 0;i<vmx.size();i++) {
    double diff = 0.0;
    diff = vmx[i] - vmr[i];
    if ( diff < fXmin ) fXmin = diff;
    diff = vmx[i] + vmr[i];
    if (diff > fXmax) fXmax = diff;
    diff = vmy[i] - vmr[i];
    if ( diff < fYmin ) fYmin = diff;
    diff = vmy[i] + vmr[i];
    if (diff > fYmax) fYmax = diff;
  }
  fDelX = (fXmax-fXmin)/nbinsx;
  fDelY = (fYmax-fYmin)/nbinsy;

  if (debug) {
    DEBUGS(fXmin);DEBUGS(fXmax);
    DEBUGS(fYmin);DEBUGS(fYmax);
    DEBUGS(fDelX);DEBUGS(fDelY);
  }

  return true;
};
/******************** end of makeGridParameters **********************/

bool GOrderedGrid::makeGrid() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GOrderedGrid::makeGrid " << endl;
  }

  int numBins = nbinsx*nbinsy;

  // prepare the grid since we know the number of bins
  for (int i = 0;i<numBins;i++) {
    vGrid.push_back(new list<GridFacet> );
  }

  if (debug) {
    *oLog << "vGrid size " << vGrid.size() << endl;
  }  

  double grid_dia = sqrt(fDelX*fDelX + fDelY*fDelY) / 2.0;

  // loop over grid bins
  for (int gb = 0;gb < numBins; ++gb) {
    // get bin center location
    int ybin = gb/nbinsx;
    int xbin = gb % nbinsx;

    double xc = (fXmin + xbin*fDelX) + (fDelX/2);
    double yc = (fYmin + ybin*fDelY) + (fDelY/2);

    // loop over elements
    for (unsigned elem=0;elem<vmx.size();++elem) {
      
      // get x,y,radius for this element
      double xe = vmx[elem];
      double ye = vmy[elem];
      double ra = vmr[elem];

      // get distance from bin center to element location
      double dst = sqrt( ((xe-xc)*(xe-xc)) + ((ye-yc)*(ye-yc)) );

      // test distance
      if (dst < (ra + grid_dia) ) {
        if (debug) {
          *oLog << "      insert following into grid " << endl;
	  *oLog << " gb elem dst  " << gb << " " << elem 
		<< " " << dst << endl;
        }
      
	GridFacet grdf(elem,dst);
        vGrid[gb]->push_back(grdf);
      }
    }
  }

  // sort the lists
  for (int gb = 0;gb < numBins; ++gb) {
    vGrid[gb]->sort(gridBinSort);
  }

  return true;
};
/********************* end of makeGrid ************************/

bool GOrderedGrid::getGridBinList(const double &x, const double &y,
                      list<GridFacet> *&gridList) {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GOrderedGrid::getGridBinList" << endl;
  }
  if ( (nbinsx==0) || (nbinsy==0) ) return false;

  bool gridFound = true;

  int xbin = 0;
  int ybin = 0;

  xbin = (int)( floor( (x - fXmin)/fDelX ));
  ybin = (int)( floor( (y - fYmin)/fDelY ));

  int gridkey = (xbin + (nbinsx*ybin) );

  if (debug) {
    DEBUGS(x);DEBUGS(xbin);
    DEBUGS(y);DEBUGS(ybin);
    DEBUGS(gridkey);
  }
 
  if ( (gridkey < 0) || (gridkey > (nbinsx*nbinsy) - 1) ) {
    gridFound = false;
  }
  else {
    gridList = vGrid[gridkey];
    if (debug) {
      *oLog << "       gridList->size() " << gridList->size() << endl;
    }
  }

  if (debug) {
    *oLog << "        gridFound " << gridFound << endl;
  }
  //exit(0);
  return gridFound;
};
/********************* end of getGridBinList ************************/

bool GOrderedGrid::readGrid() {

  bool debug = false;

  if (debug) {
    *oLog << "  -- GOrderedGrid::readGrid" << endl;
  }

  ifstream inFile(sFileName.c_str(),ios::in);
  if (! inFile) {
    return false;   // no such file
  }
  string fileRec = "";
  vector<string> tokens;

  // skip the first three documentation records
  getline(inFile,fileRec,'\n');
  getline(inFile,fileRec,'\n');
  getline(inFile,fileRec,'\n');

  // read grid parameters
  getline(inFile,fileRec,'\n');
  GUtilityFuncts::tokenizer(fileRec,tokens);

  nbinsx = atoi(tokens[0].c_str());
  fDelX = atof(tokens[1].c_str());
  fXmin = atof(tokens[2].c_str());
  nbinsy = atoi(tokens[3].c_str());
  fDelY = atof(tokens[4].c_str());
  fYmin = atof(tokens[5].c_str());

  if (debug) {
    DEBUGS(nbinsx);DEBUGS(fDelX);DEBUGS(fXmin);
    DEBUGS(nbinsy);DEBUGS(fDelY);DEBUGS(fYmin);
  }

  int numBins = nbinsx*nbinsy;
    
  // start with a new grid if it's already filled
  for (unsigned i = 0;i<vGrid.size();i++) {
    SafeDelete(vGrid[i]);
    *oLog << " vGrid " << vGrid[i] << endl;
  }
  vGrid.clear();

  // prepare the grid since we know the number of bins
  for (int i = 0;i<numBins;i++) {
    vGrid.push_back(new list<GridFacet> );
  } 
  
  for (int i = 0;i<nbinsx*nbinsy;i++) {

    string elements  = "";
    string distances = "";
    vector<string> tokElem;
    vector<string> tokDist;

    // read element list
    getline(inFile,elements,'\n');
    GUtilityFuncts::tokenizer(elements,tokElem);

    int gridkey = atoi(tokElem[0].c_str());
    int numElem = atoi(tokElem[3].c_str());

    // read distance list
    getline(inFile,distances,'\n');
    GUtilityFuncts::tokenizer(distances,tokDist);

    // loop over elements
    for (int el=0;el<numElem;++el) {

      int elem = atoi(tokElem[4+el].c_str());
      double dist = atof(tokDist[4+el].c_str());
 
      GridFacet grdf(elem-1,dist);
      vGrid[gridkey]->push_back(grdf);
    }
  }

  return true;
};
/********************* end of readGrid ************************/

bool GOrderedGrid::printGrid() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GOrderedGrid::printGrid " << endl;
    *oLog << "      opening file " << sFileName << endl;
  }

  if (sFileName=="") {
    return false;
  }

  FILE *fg = 0;
  if ( (fg= fopen(sFileName.c_str(),"w")) == NULL) {
    cerr << "    failed to open facet grid file "
         << sFileName << endl;
    exit(0);
  }
  
  fprintf(fg,"first data record: nbinsx delx xmin nbinsy dely ymin\n");
  fprintf(fg,"listing 1: gridKey,xbin,ybin,numElements,list of elements\n");
  fprintf(fg,"  followed by same format, list of distances\n");

  fprintf(fg,"%d %lf %lf %d %lf %lf\n",nbinsx,fDelX,fXmin,
          nbinsy,fDelY,fYmin);

  int numBins = nbinsx*nbinsy;
  
  for (int i = 0;i<numBins;++i) {
    int yb = i / nbinsx;
    int xb = i % nbinsx;

    // get list for this element
    list<GridFacet> *lst = vGrid[i];
    list<GridFacet>::iterator iterlst;
    
    unsigned numElem = vGrid[i]->size();
    fprintf(fg,"%4d %3d %3d %2d ",i, xb,yb,numElem);
 
    if (numElem == 0) {
      fprintf(fg,"\n");
    }
    else {
      
      for (iterlst = lst->begin();
           iterlst!= lst->end(); ++iterlst) {

        int facetN = iterlst->facetNum;
        fprintf(fg,"%3d ",facetN + 1);
        
      }
      fprintf(fg,"\n");
    
    }
    fprintf(fg,"%4d %3d %3d %2d ",i, xb,yb,numElem);
    if (numElem == 0) {
      fprintf(fg,"\n");
    }
    else {
      
      for (iterlst = lst->begin();
           iterlst!= lst->end(); ++iterlst) {

        double distg = iterlst->dist;
        fprintf(fg,"%5.2f ",distg);
        
      }
      fprintf(fg,"\n");
     
    }

   }

  fclose(fg);
  return true;
};

/********************* end of printGrid ************************/


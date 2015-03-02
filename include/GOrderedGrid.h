/*
VERSION3.1
2March2015
*/
/*!  /brief GOrderedGrid class for created an ordered 2D
     hash table, e.g. used for quickly finding facet number
     given hit location on the mirror.
 */

#ifndef GORDEREDGRID
#define GORDEREDGRID

#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Rotation3Dfwd.h"

// forward declarations

struct GridFacet {
  int facetNum;
  double dist;

  GridFacet(const int &facetNum1=0, const double &dist1=0.0) {
    facetNum = facetNum1;
    dist = dist1;
  };
};


class GOrderedGrid {
  
  vector<double> vmx;
  vector<double> vmy;
  vector<double> vmr;
  int nbinsx;
  int nbinsy;  
  int iGridOption;

  string sFileName;

  double fXmin;
  double fXmax;
  double fYmin;
  double fYmax;

  double fDelX;
  double fDelY;

  vector< list<GridFacet> * > vGrid;

  void initialize();

  bool makeGridParameters();

  bool makeGrid();

  bool readGrid();

  bool printGrid();

 public:
  GOrderedGrid(const vector<double> &xe,const vector<double> &ye,
               const vector<double> &re,
               const int &nbinsx, const int &nbinsy, const int &option = 1,
               const string &filename = "");

  ~GOrderedGrid();

  //GOrderedGrid(const GOrderedGrid &orderedGrid);

  bool getGridBinList(const double &x, const double &y,
                 list<GridFacet> *&gridList);

};

#endif

/*
VERSION3.1
2March2015
*/
/*!  /brief GArrayTel class contains all telescope details 
            including a pointer to a GTelescope instance
 */

#ifndef GSEGMENTEDMIRROR
#define GSEGMENTEDMIRROR

// forward declarations
class AMirror;
class AGeoAsphericDisk;
class TGeoCombiTrans;

class SegmentedMirror
{
  // The base class for segmented mirrors. This class holds some Parameters
  // to define the shape and alignment errors of a mirror. But the surfaces of
  // the primary and secondary mirrors are given in NewSCT class.

private:

protected:
  Double_t fRmin; // The innner radius of a sector
  Double_t fRmax; // The outer radius of a sector
  Double_t fPhimin; // The minimum azimuthal angle
  Double_t fPhimax; // The maximum azimuthal angle

  Double_t fMargin; // Margin around a mirror
  Double_t fPositionErrorX; // Alignment error along X axis
  Double_t fPositionErrorY; // Alignment error along Y axis
  Double_t fPositionErrorZ; // Alignment error along Z axis
  Double_t fRotationErrorXY; // Rotation error in X-Y plane
  Double_t fRotationErrorSagittal; // Rotation error in sagittal plane
  Double_t fRotationErrorTangential; // Rotation error in tangential plane
  Double_t fRoughness; // Mirror roughness (deg)

public:
  SegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  virtual ~SegmentedMirror() {}

  virtual AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary) = 0;
  TGeoCombiTrans* BuildMirrorCombiTrans(AGeoAsphericDisk* disk, Bool_t isPrimary);

  Double_t GetRoughness() const {return fRoughness;}

  void SetMargin(Double_t margin);
  void SetPositionErrors(Double_t x, Double_t y, Double_t z);
  void SetRotationErrors(Double_t xy, Double_t sag, Double_t tan);
  void SetRougness(Double_t roughness);
};

//______________________________________________________________________________
class SectorSegmentedMirror : public SegmentedMirror
{
private:

public:
  SectorSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~SectorSegmentedMirror(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};

//______________________________________________________________________________
class TetragonSegmentedMirror : public SegmentedMirror
{
private:

public:
  TetragonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~TetragonSegmentedMirror(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};
//______________________________________________________________________________
class PentagonSegmentedMirror : public SegmentedMirror
{
private:

public:
  PentagonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~PentagonSegmentedMirror(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};







#endif

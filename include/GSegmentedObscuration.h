/*
VERSION3.1
2March2015
*/
/*!  /brief GArrayTel class contains all telescope details 
            including a pointer to a GTelescope instance
 */

#ifndef GSEGMENTEDOBSCURATION
#define GSEGMENTEDOBSCURATION

// forward declarations
class AObscuration;
class AGeoAsphericDisk;
class TGeoCombiTrans;

class SegmentedObscuration
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

public:
  SegmentedObscuration(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  virtual ~SegmentedObscuration() {}

  virtual AObscuration* BuildObscuration(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary) = 0;
  TGeoCombiTrans* BuildObscurationCombiTrans(AGeoAsphericDisk* disk, Bool_t isPrimary);

  void SetMargin(Double_t margin);
  void SetPositionErrors(Double_t x, Double_t y, Double_t z);
  void SetRotationErrors(Double_t xy, Double_t sag, Double_t tan);
};

//______________________________________________________________________________
class SectorSegmentedObscuration : public SegmentedObscuration
{
private:

public:
  SectorSegmentedObscuration(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~SectorSegmentedObscuration(){}

  AObscuration* BuildObscuration(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};

//______________________________________________________________________________
class TetragonSegmentedObscuration : public SegmentedObscuration
{
private:

public:
  TetragonSegmentedObscuration(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~TetragonSegmentedObscuration(){}

  AObscuration* BuildObscuration(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};
//______________________________________________________________________________
class PentagonSegmentedObscuration : public SegmentedObscuration
{
private:

public:
  PentagonSegmentedObscuration(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax);
  ~PentagonSegmentedObscuration(){}

  AObscuration* BuildObscuration(const char* name, AGeoAsphericDisk* disk, Bool_t isPrimary);
};







#endif

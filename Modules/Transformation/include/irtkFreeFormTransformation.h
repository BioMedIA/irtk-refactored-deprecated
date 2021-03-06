/* The Image Registration Toolkit (IRTK)
 *
 * Copyright 2008-2015 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. */

#ifndef _IRTKFREEFORMTRANSFORMATION_H
#define _IRTKFREEFORMTRANSFORMATION_H

#include <irtkTransformationJacobian.h>
#include <irtkImageFunction.h>


/**
 * Base class for free-form transformations.
 */

class irtkFreeFormTransformation : public irtkTransformation
{
  irtkAbstractTransformationMacro(irtkFreeFormTransformation);

public:

  /// Default attributes of free-form transformation lattice
  ///
  /// \param[in] attr (Target) image attributes.
  /// \param[in] dx   Control point spacing in x.
  /// \param[in] dy   Control point spacing in y.
  /// \param[in] dz   Control point spacing in z.
  /// \param[in] dt   Control point spacing in t.
  ///
  /// \note If a control point number/spacing is negative, it is set to the specified number/size.
  static irtkImageAttributes DefaultAttributes(const irtkImageAttributes &attr,
                                               double dx = -1.0,
                                               double dy = -1.0,
                                               double dz = -1.0,
                                               double dt = -1.0);

  /// Default attributes of free-form transformation lattice
  static irtkImageAttributes DefaultAttributes(double, double, double, double,
                                               double, double, double, double,
                                               double, double, double, double,
                                               const double *, const double *,
                                               const double *);

  // ---------------------------------------------------------------------------
  // Types

  /// Type of vector storing the coefficients at a control point
  typedef irtkVector3D<DOFValue>                       CPValue;

  /// Type of image representation of free-form transformation
  typedef irtkGenericImage<CPValue>                    CPImage;

  /// Type of vector storing status of control point coefficients
  typedef irtkVector3D<DOFStatus>                      CPStatus;

  /// Base type used for interpolation of control point data
  /// \note Subclasses may use a particular specialized interpolation type.
  typedef irtkGenericInterpolateImageFunction<CPImage> CPInterpolator;

  /// Base type used for extrapolation of control point data
  typedef irtkGenericExtrapolateImageFunction<CPImage> CPExtrapolator;

  /// Alias for CPValue type which is preferred when referring to any FFD vector,
  /// not only those at control point locations
  typedef CPValue                                      Vector;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Mode used for extrapolation of control point image
  irtkReadOnlyAttributeMacro(irtkExtrapolationMode, ExtrapolationMode);

  /// Speedup factor for gradient computation
  irtkPublicAttributeMacro(double, SpeedupFactor);

public:

  /// Set extrapolation mode
  virtual void ExtrapolationMode(irtkExtrapolationMode);

  /// Get access to control point coefficients extrapolator
  /// Intended for use by non-member auxiliary functions of subclass implementations.
  const CPExtrapolator *Extrapolator() const;

protected:

  /// Finite discrete representation of free-form transformation coefficients
  ///
  /// Note that this image instance stores all the attributes of the control
  /// point lattice. The actual control point data is stored in the parameters
  /// of the transformation, the memory of which is managed by the base class.
  /// The image only wraps this memory and interprets it as an image. This allows
  /// the use of image functions such as interpolators and efficient image filters.
  CPImage _CPImage;

  /// Infinite discrete representation of free-form transformation coefficients
  CPExtrapolator *_CPValue;

  /// Infinite continuous function of free-form transformation
  ///
  /// Note that the actual interpolator object is (at the moment) attribute of
  /// the specific FFD implementation. The base class pointer here is only used
  /// to update the interpolator when necessary, e.g., when the extrapolation
  /// mode has changed. It is set by the base class constructor to the
  /// respective interpolate image function attribute of the subclass.
  CPInterpolator *_CPFunc;
  CPInterpolator *_CPFunc2D; // TODO: Remove when 2D FFDs are separate classes

  /// Status of control points
  ///
  /// Note that the actual status values are stored by the base class member
  /// irtkTransformation::_Status in a contiguous memory block.
  CPStatus ****_CPStatus;

  const irtkImageAttributes &_attr; ///< Control point lattice attributes

  const int    &_x;  ///< Read-only reference to _x  attribute of _CPImage
  const int    &_y;  ///< Read-only reference to _y  attribute of _CPImage
  const int    &_z;  ///< Read-only reference to _z  attribute of _CPImage
  const int    &_t;  ///< Read-only reference to _t  attribute of _CPImage

  const double &_dx; ///< Read-only reference to _dx attribute of _CPImage
  const double &_dy; ///< Read-only reference to _dy attribute of _CPImage
  const double &_dz; ///< Read-only reference to _dz attribute of _CPImage
  const double &_dt; ///< Read-only reference to _dt attribute of _CPImage

  const irtkMatrix &_matW2L; ///< Read-only reference to _matW2I of _CPImage
  const irtkMatrix &_matL2W; ///< Read-only reference to _matI2W of _CPImage

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Initialize interpolator of control points
  void InitializeInterpolator();

  /// Initialize status of control points
  void InitializeStatus();

  /// Initialize control points
  void InitializeCPs(const irtkImageAttributes &, bool = true);

  /// Copy control points from other transformation
  void InitializeCPs(const irtkFreeFormTransformation &, bool = true);

  /// Default constructor
  irtkFreeFormTransformation(CPInterpolator &, CPInterpolator * = NULL);

  /// Copy constructor
  irtkFreeFormTransformation(const irtkFreeFormTransformation &,
                             CPInterpolator &, CPInterpolator * = NULL);

public:

  /// Destructor
  virtual ~irtkFreeFormTransformation();

  /// Initialize free-form transformation
  virtual void Initialize(const irtkImageAttributes &) = 0;

  /// Initialize free-form transformation
  void Initialize(const irtkImageAttributes &, double, double, double = -1.0, double = -1.0);

  /// Initialize transformation from existing vector field
  void Initialize(const CPImage &, bool = false);

  /// Initialize transformation from existing 3D+t vector field
  void Initialize(const irtkGenericImage<double> &, bool = false);

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other Parameter overloads
  using irtkTransformation::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  using irtkTransformation::EvaluateRMSError;
  using irtkTransformation::Approximate;
  using irtkTransformation::ApproximateAsNew;
  using irtkTransformation::ApproximateGradient;

  /// Get image domain on which this free-form transformation should be defined
  /// in order to reduce the error when the given transformation is approximated
  /// and the resulting free-form transformation is applied to resample an image
  /// that is defined on the specified lattice.
  virtual irtkImageAttributes ApproximationDomain(const irtkImageAttributes &,
                                                  const irtkTransformation *);

  /// Evaluates RMS error of transformation at control points compared to another.
  double EvaluateRMSError(const irtkTransformation *) const;

  /// Approximate transformation: This function takes a transformation and finds
  /// a FFD which approximates this transformation at the control point locations.
  /// Returns the approximation error of the resulting FFD at the control point locations.
  virtual double Approximate(const irtkTransformation *, int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const irtkImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(const double *, const double *, const double *, const double *,
                             double *,       double *,       double *,       int,
                             int = 1, double = .0);

  /// Approximate another transformation and return approximation error
  virtual double ApproximateAsNew(const irtkTransformation *, int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const irtkImageAttributes &, double *, double *, double *,
                                  int = 1, double = .0);

  /// Interpolates displacements: This function takes a set of displacements defined
  /// at the control points and finds a FFD which interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double *) = 0;

  // ---------------------------------------------------------------------------
  // Lattice

  /// Returns attributes of control point grid
  const irtkImageAttributes &Attributes() const;

  /// Actual number of DoFs, i.e.,
  /// 2 * NumberOfDOFs for 2D FFD or 3 * NumberOfDOFs for 3D(+t) FFD
  int ActualNumberOfDOFs() const;
  
  /// Number of control points
  int NumberOfCPs() const;
  
  /// Number of active control points
  int NumberOfActiveCPs() const;
  
  /// Number of non-active control points
  int NumberOfPassiveCPs() const;

  /// Returns the number of control points in x
  int X() const;
  
  /// Returns the number of control points in y
  int Y() const;
  
  /// Returns the number of control points in z
  int Z() const;
  
  /// Returns the number of control points in t
  int T() const;

  /// Returns the number of control points in x
  int GetX() const;

  /// Returns the number of control points in y
  int GetY() const;

  /// Returns the number of control points in z
  int GetZ() const;

  /// Returns the number of control points in t
  int GetT() const;

  /// Returns the of control point spacing in x
  double GetXSpacing() const;
  
  /// Returns the of control point spacing in y
  double GetYSpacing() const;
  
  /// Returns the of control point spacing in z
  double GetZSpacing() const;
  
  /// Returns the of control point spacing in t
  double GetTSpacing() const;

  /// Gets the control point spacing (in mm)
  void GetSpacing(double &, double &, double &) const;

  /// Gets the control point spacing (in mm)
  void GetSpacing(double &, double &, double &, double &) const;

  /// Puts the orientation of the free-form deformation lattice
  void PutOrientation(double *, double *, double *);
  
  /// Gets the orientation of the free-form deformation lattice
  void GetOrientation(double *, double *, double *) const;

  /// Get indices of transformation parameters (DoFs)
  void IndexToDOFs(int, int &, int &) const;

  /// Get indices of transformation parameters (DoFs)
  void IndexToDOFs(int, int &, int &, int &) const;

  /// Get index of control point corresponding to transformation parameter (DoFs)
  int DOFToIndex(int) const;

  /// Get index of dimension corresponding to transformation parameter (DoFs)
  int DOFToDimension(int) const;

  /// Get control point index from lattice coordinates
  int LatticeToIndex(int, int, int = 0, int = 0) const;

  /// Get control point lattice coordinates from index
  void IndexToLattice(int, int &, int &) const;

  /// Get control point lattice coordinates from index
  void IndexToLattice(int, int &, int &, int &) const;

  /// Get control point lattice coordinates from index
  void IndexToLattice(int, int &, int &, int &, int &) const;

  /// Get world coordinates (in mm) of control point
  void IndexToWorld(int, double &, double &) const;

  /// Get world coordinates (in mm) of control point
  void IndexToWorld(int, double &, double &, double &) const;

  /// Get world coordinates (in mm) of control point
  void IndexToWorld(int, irtkPoint &) const;

  /// Get world coordinates (in mm) of control point
  irtkPoint IndexToWorld(int) const;

  /// Transforms world coordinates (in mm) to lattice coordinates
  virtual void WorldToLattice(double &, double &) const;

  /// Transforms world coordinates (in mm) to lattice coordinates
  virtual void WorldToLattice(double &, double &, double &) const;

  /// Transforms world coordinates (in mm) to lattice coordinates
  virtual void WorldToLattice(irtkPoint &) const;

  /// Transforms lattice coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &, double &) const;

  /// Transforms lattice coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &, double &, double &) const;

  /// Transforms lattice coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(irtkPoint &) const;

  /// Gets the location of the given control point (in mm)
  void ControlPointLocation(int, double &, double &) const;

  /// Gets the location of the given control point (in mm)
  void ControlPointLocation(int, double &, double &, double &) const;

  /// Returns the location of the given control point (in mm)
  irtkPoint ControlPointLocation(int) const;

  /// Transforms time (in ms) to temporal lattice coordinate
  virtual double TimeToLattice(double) const;

  /// Transforms temporal lattice coordinate to time (in ms)
  virtual double LatticeToTime(double) const;

  /// Returns the number of control points in x after subdivision
  virtual int GetXAfterSubdivision() const;

  /// Returns the number of control points in y after subdivision
  virtual int GetYAfterSubdivision() const;

  /// Returns the number of control points in z after subdivision
  virtual int GetZAfterSubdivision() const;

  /// Returns the number of control points in t after subdivision
  virtual int GetTAfterSubdivision() const;

  /// Returns the control point spacing in x after the subdivision
  virtual double GetXSpacingAfterSubdivision() const;

  /// Returns the control point spacing in y after the subdivision
  virtual double GetYSpacingAfterSubdivision() const;

  /// Returns the control point spacing in z after the subdivision
  virtual double GetZSpacingAfterSubdivision() const;

  /// Returns the control point spacing in t after the subdivision
  virtual double GetTSpacingAfterSubdivision() const;

  /// Subdivide lattice of free-form transformation
  virtual void Subdivide(bool = true, bool = true, bool = true, bool = true);

  /// Subdivide lattice in first two dimensions
  void Subdivide2D();

  /// Subdivide lattice in first three dimensions
  void Subdivide3D();

  /// Subdivide lattice in all four dimensions
  void Subdivide4D();

  /// Crop/pad lattice to discard passive control points at the boundary,
  /// keeping only a layer of passive control points of given width.
  /// The DoF values of passive control points are optionally reset to zero.
  bool CropPadPassiveCPs(int = 0, bool = false);

  /// Crop/pad lattice to discard passive control points at the boundary,
  /// keeping only a layer of passive control points of given width.
  /// The DoF values of passive control points are optionally reset to zero.
  virtual bool CropPadPassiveCPs(int, int, int = 0, int = 0, bool = false);

  // ---------------------------------------------------------------------------
  // Bounding box

  /// Size of support region of the used kernel
  virtual int KernelSize() const = 0;

  /// Radius of support region of the used kernel
  int KernelRadius() const;

  /// Puts the spatial bounding box for the free-form deformation (in mm)
  void PutBoundingBox(double, double, double,
                      double, double, double);
  
  /// Puts the temporal bounding box for the free-form deformation (in mm)
  void PutBoundingBox(const irtkPoint &, const irtkPoint &);

  /// Puts the temporal bounding box of the free-form deformation (in ms)
  void PutBoundingBox(double, double);

  /// Puts the spatio-temporal bounding box for the free-form deformation (in mm)
  void PutBoundingBox(double, double, double, double,
                      double, double, double, double);

  /// Gets the temporal bounding box of the free-form deformation (in ms)
  void BoundingBox(double &, double &) const;

  /// Gets the spatial bounding box of the free-form deformation (in mm)
  void BoundingBox(double &, double &, double &,
                   double &, double &, double &) const;

  /// Gets the spatial bounding box of the free-form deformation (in mm)
  void BoundingBox(irtkPoint &, irtkPoint &) const;

  /// Gets the spatio-temporal bounding box of the free-form deformation (in mm and ms)
  void BoundingBox(double &, double &, double &, double &,
                   double &, double &, double &, double &) const;

  /// Gets the spatio-temporal bounding box of the free-form deformation (in mm and ms)
  void BoundingBox(irtkPoint &, double &, irtkPoint &, double &) const;

  /// Gets the temporal bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double = 1) const;

  /// Gets the spatial bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double &,
                                double &, double &, double &, double = 1) const = 0;

  /// Gets the spatio-temporal bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  virtual void BoundingBox(int, double &, double &, double &, double &,
                                double &, double &, double &, double &, double = 1) const;

  /// Gets the spatial bounding box for a control point. The last parameter
  /// specifies what fraction of the bounding box to return. The default
  /// is 1 which equals 100% of the bounding box.
  void BoundingBox(int, irtkPoint &, irtkPoint &, double = 1) const;

  /// Gets the spatial bounding box for a control point in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  bool BoundingBox(const irtkImage *, int, int &, int &, int &,
                                           int &, int &, int &, double = 1) const;

  /// Gets the spatio-temporal bounding box for a control point in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  bool BoundingBox(const irtkImage *, int, int &, int &, int &, int &,
                                           int &, int &, int &, int &, double = 1) const;

  /// Gets the spatial bounding box for a transformation parameter in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  bool DOFBoundingBox(const irtkImage *, int, int &, int &, int &,
                                              int &, int &, int &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)
  using irtkTransformation::Put;
  using irtkTransformation::Get;
  using irtkTransformation::PutStatus;
  using irtkTransformation::GetStatus;

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const irtkTransformation *);

  /// Puts values of the parameters at a control point
  void Put(int, const Vector &);

  /// Puts values of the parameters at a control point
  void Put(int, int, int, double, double, double);
  
  /// Puts values of the parameters at a control point
  void Put(int, int, int, int, double, double, double);

  /// Gets values of the parameters at a control point
  void Get(int, Vector &) const;

  /// Gets values of the parameters at a control point
  void Get(int, int, int, double &, double &, double &) const;
  
  /// Gets values of the parameters at a control point
  void Get(int, int, int, int, double &, double &, double &) const;

  /// Puts status of the parameters at a control point
  void PutStatus(int, const CPStatus &);

  /// Puts status of the parameters at a control point
  void PutStatus(int, int, int, DOFStatus, DOFStatus, DOFStatus);

  /// Puts status of the parameters at a control point
  void PutStatus(int, int, int, int, DOFStatus, DOFStatus, DOFStatus);

  /// Gets status of the parameters at a control point
  void GetStatus(int, CPStatus &) const;

  /// Gets status of the parameters at a control point
  void GetStatus(int, int, int, DOFStatus &, DOFStatus &, DOFStatus &) const;

  /// Gets status of the parameters at a control point
  void GetStatus(int, int, int, int, DOFStatus &, DOFStatus &, DOFStatus &) const;

  /// Whether the control point at given lattice index is active
  virtual bool IsActive(int) const;

  /// Whether the control point at given lattice coordinates is active
  virtual bool IsActive(int, int, int = 0, int = 0) const;

  // ---------------------------------------------------------------------------
  // Point transformation
  using irtkTransformation::Transform;
  using irtkTransformation::Inverse;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using irtkTransformation::Jacobian;
  using irtkTransformation::GlobalJacobian;
  using irtkTransformation::JacobianDOFs;
  using irtkTransformation::ParametricGradient;

  /// Convert 1st order derivatives computed w.r.t 2D lattice coordinates to
  /// derivatives w.r.t world coordinates
  void JacobianToWorld(double &, double &) const;

  /// Convert 1st order derivatives computed w.r.t 3D lattice coordinates to
  /// derivatives w.r.t world coordinates
  void JacobianToWorld(double &, double &, double &) const;

  /// Convert 1st order derivatives computed w.r.t lattice coordinates to
  /// derivatives w.r.t world coordinates
  void JacobianToWorld(irtkMatrix &) const;

  /// Convert 2nd order derivatives computed w.r.t 2D lattice coordinates to
  /// derivatives w.r.t world coordinates
  void HessianToWorld(double &, double &, double &) const;

  /// Convert 2nd order derivatives computed w.r.t 3D lattice coordinates to
  /// derivatives w.r.t world coordinates
  void HessianToWorld(double &, double &, double &, double &, double &, double &) const;

  /// Convert 2nd order derivatives of single transformed coordinate computed
  /// w.r.t lattice coordinates to derivatives w.r.t world coordinates
  void HessianToWorld(irtkMatrix &) const;

  /// Convert 2nd order derivatives computed w.r.t lattice coordinates to
  /// derivatives w.r.t world coordinates
  void HessianToWorld(irtkMatrix [3]) const;

  /// Calculates the Jacobian of the transformation w.r.t either control point displacements or velocities
  virtual void FFDJacobianWorld(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the global transformation w.r.t world coordinates
  virtual void GlobalHessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters of a control point
  virtual void JacobianDOFs(irtkMatrix &, int, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(irtkTransformationJacobian &, double, double, double, double = 0, double = -1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  ///
  /// If the transformation itself is non-parametric, the gradient will be passed through
  /// unchanged. The default implementation uses the full Jacobian matrix computed for each
  /// DoF separately (i.e., calls JacobianDOFs for each DoF).
  ///
  /// For 4D transformations, the temporal coordinate t used for the computation of the Jacobian
  /// of the transformation w.r.t the transformation parameters for all spatial (x, y, z) voxel
  /// coordinates in world units, is assumed to correspond to the temporal origin of the given
  /// gradient image. For 4D transformations parameterized by velocities, a second time for the
  /// upper integration bound can be provided as last argument to this method. This last
  /// argument is ignored by transformations parameterized by displacements.
  ///
  /// \sa irtkImageSimilarityMetric::EvaluateGradient
  virtual void ParametricGradient(const irtkGenericImage<double> *, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  double = -1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkPointSet &, const irtkVector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Calculates the bending of the transformation given the 2nd order derivatives
  static double Bending3D(const irtkMatrix [3]);

  /// Calculates the bending of the transformation
  virtual double BendingEnergy(double, double, double, double = 0, double = -1, bool = true) const;

  /// Approximates the bending energy on the control point lattice
  virtual double BendingEnergy(bool = false, bool = true) const;

  /// Approximates the bending energy on the specified discrete domain
  virtual double BendingEnergy(const irtkImageAttributes &attr, double = -1, bool = true) const;

  /// Approximates the gradient of the bending energy on the control point
  /// lattice w.r.t the transformation parameters and adds it with the given weight
  virtual void BendingEnergyGradient(double *, double = 1.0, bool = false, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

protected:

  /// Writes the control point and status information to a file stream
  irtkCofstream &WriteCPs(irtkCofstream &) const;

  // ---------------------------------------------------------------------------
  // Backwards compatibility

public:

  /// \deprecated Use BendingEnergy() instead.
  virtual double Bending(double, double, double) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline const irtkFreeFormTransformation::CPExtrapolator *irtkFreeFormTransformation::Extrapolator() const
{
  return _CPValue;
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::KernelRadius() const
{
  return this->KernelSize() / 2;
}

// -----------------------------------------------------------------------------
inline const irtkImageAttributes &irtkFreeFormTransformation::Attributes() const
{
  return _CPImage.Attributes();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::NumberOfCPs() const
{
  return _x * _y * _z * _t;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::NumberOfActiveCPs() const
{
  int nactive = 0;
  for (int cp = 0; cp < NumberOfCPs(); ++cp) {
    if (this->IsActive(cp)) ++nactive;
  }
  return nactive;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::NumberOfPassiveCPs() const
{
  return NumberOfCPs() - NumberOfActiveCPs();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::ActualNumberOfDOFs() const
{
  if (_z == 1) return 2 * NumberOfCPs();
  else         return 3 * NumberOfCPs();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::X() const
{
  return _x;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::Y() const
{
  return _y;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::Z() const
{
  return _z;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::T() const
{
  return _t;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetX() const
{
  return X();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetY() const
{
  return Y();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetZ() const
{
  return Z();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetT() const
{
  return T();
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetXSpacing() const
{
  return _dx;
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetYSpacing() const
{
  return _dy;
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetZSpacing() const
{
  return _dz;
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetTSpacing() const
{
  return _dt;
}

// ---------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GetSpacing(double &dx, double &dy, double &dz) const
{
  dx = _dx;
  dy = _dy;
  dz = _dz;
}

// ---------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GetSpacing(double &dx, double &dy, double &dz, double &dt) const
{
  dx = _dx;
  dy = _dy;
  dz = _dz;
  dt = _dt;
}

// ---------------------------------------------------------------------------
inline void irtkFreeFormTransformation::PutOrientation(double *xaxis, double *yaxis, double *zaxis)
{
  _CPImage.PutOrientation(xaxis, yaxis, zaxis);
}

// ---------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GetOrientation(double *xaxis, double *yaxis, double *zaxis) const
{
  _CPImage.GetOrientation(xaxis, yaxis, zaxis);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToDOFs(int cp, int &x, int &y) const
{
  x = 3 * cp;
  y = x + 1;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToDOFs(int cp, int &x, int &y, int &z) const
{
  x = 3 * cp;
  y = x + 1;
  z = x + 2;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::DOFToIndex(int dof) const
{
  return dof / 3;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::DOFToDimension(int dof) const
{
  return dof % 3;
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::LatticeToIndex(int i, int j, int k, int l) const
{
  return _CPImage.VoxelToIndex(i, j, k, l);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToLattice(int index, int &i, int &j) const
{
  _CPImage.IndexToVoxel(index, i, j);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToLattice(int index, int &i, int &j, int &k) const
{
  _CPImage.IndexToVoxel(index, i, j, k);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToLattice(int index, int &i, int &j, int &k, int &l) const
{
  _CPImage.IndexToVoxel(index, i, j, k, l);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToWorld(int index, double &x, double &y) const
{
  _CPImage.IndexToWorld(index, x, y);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToWorld(int index, double &x, double &y, double &z) const
{
  _CPImage.IndexToWorld(index, x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::IndexToWorld(int index, irtkPoint &p) const
{
  _CPImage.IndexToWorld(index, p);
}

// -----------------------------------------------------------------------------
inline irtkPoint irtkFreeFormTransformation::IndexToWorld(int index) const
{
  return _CPImage.IndexToWorld(index);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::WorldToLattice(double &x, double &y) const
{
  _CPImage.WorldToImage(x, y);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::WorldToLattice(double &x, double &y, double &z) const
{
  _CPImage.WorldToImage(x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::WorldToLattice(irtkPoint &p) const
{
  this->WorldToLattice(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::LatticeToWorld(double &x, double &y) const
{
  _CPImage.ImageToWorld(x, y);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::LatticeToWorld(double &x, double &y, double &z) const
{
  _CPImage.ImageToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::LatticeToWorld(irtkPoint &p) const
{
  this->LatticeToWorld(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline irtkPoint irtkFreeFormTransformation::ControlPointLocation(int cp) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);
  irtkPoint p(i, j, k);
  this->LatticeToWorld(p);
  return p;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::ControlPointLocation(int cp, double &x, double &y) const
{
  int i, j;
  IndexToLattice(cp, i, j);
  x = i, y = j;
  this->LatticeToWorld(x, y);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::ControlPointLocation(int cp, double &x, double &y, double &z) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);
  x = i, y = j, z = k;
  this->LatticeToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::TimeToLattice(double t) const
{
  return _CPImage.TimeToImage(t);
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::LatticeToTime(double t) const
{
  return _CPImage.ImageToTime(t);
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetXAfterSubdivision() const
{
  return this->GetX();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetYAfterSubdivision() const
{
  return this->GetY();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetZAfterSubdivision() const
{
  return this->GetZ();
}

// -----------------------------------------------------------------------------
inline int irtkFreeFormTransformation::GetTAfterSubdivision() const
{
  return this->GetT();
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetXSpacingAfterSubdivision() const
{
  return this->GetXSpacing();
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetYSpacingAfterSubdivision() const
{
  return this->GetYSpacing();
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetZSpacingAfterSubdivision() const
{
  return this->GetZSpacing();
}

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::GetTSpacingAfterSubdivision() const
{
  return this->GetTSpacing();
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Subdivide(bool subdivide_x, bool subdivide_y, bool subdivide_z, bool subdivide_t)
{
  if (subdivide_x || subdivide_y || subdivide_z || subdivide_t) {
    cerr << this->NameOfClass() << "::Subdivide: Not implemented (or not possible?)" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Subdivide2D()
{
  this->Subdivide(true, true, false, false);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Subdivide3D()
{
  this->Subdivide(true, true, true, false);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Subdivide4D()
{
  this->Subdivide(true, true, true, true);
}

// =============================================================================
// Bounding box
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(double &t1, double &t2) const
{
  t1 = this->LatticeToTime(.0);
  t2 = this->LatticeToTime(_t - 1);
  if (t1 > t2) swap(t1, t2);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(double &x1, double &y1, double &z1,
                                                    double &x2, double &y2, double &z2) const
{
  x1 = y1 = z1 = .0;
  this->LatticeToWorld(x1, y1, z1);
  x2 = _x - 1, y2 = _y - 1, z2 = _z - 1;
  this->LatticeToWorld(x2, y2, z2);
  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(irtkPoint &p1, irtkPoint &p2) const
{
  BoundingBox(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(double &x1, double &y1, double &z1, double &t1,
                                                    double &x2, double &y2, double &z2, double &t2) const
{
  BoundingBox(x1, y1, z1, x2, y2, z2);
  BoundingBox(t1, t2);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(irtkPoint &p1, double &t1, irtkPoint &p2, double &t2) const
{
  BoundingBox(p1._x, p1._y, p1._z, t1, p2._x, p2._y, p2._z, t2);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(int, double &t1, double &t2, double) const
{
  t1 = this->LatticeToTime(0);
  t2 = this->LatticeToTime(_t - 1);
  if (t1 > t2) swap(t1, t2);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(int, double &x1, double &y1, double &z1,
                                                         double &x2, double &y2, double &z2, double) const
{
  x1 = 0,      y1 = 0,      z1 = 0;
  x2 = _x - 1, y2 = _y - 1, z2 = _z - 1;
  this->LatticeToWorld(x1, y1, z1);
  this->LatticeToWorld(x2, y2, z2);
  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(int cp, double &x1, double &y1, double &z1, double &t1,
                                                            double &x2, double &y2, double &z2, double &t2,
                                                            double fraction) const
{
  this->BoundingBox(cp, x1, y1, z1, x2, y2, z2, fraction);
  this->BoundingBox(cp, t1, t2,                 fraction);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::BoundingBox(int cp, irtkPoint &p1, irtkPoint &p2, double fraction) const
{
  BoundingBox(cp, p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, fraction);
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::BoundingBox(const irtkImage *image, int cp,
                                                    int &i1, int &j1, int &k1,
                                                    int &i2, int &j2, int &k2,
                                                    double fraction) const
{
  // Calculate bounding box in world coordinates parallel to world axes
  double x[2], y[2], z[2];
  this->BoundingBox(cp, x[0], y[0], z[0], x[1], y[1], z[1], fraction);

  // Map bounding box into image space and calculate minimal axes-alinged
  // bounding box parallel to image axes which need not coincide with world axes
  irtkPoint p;
  double x1, y1, z1, x2, y2, z2;
  x1 = y1 = z1 = + numeric_limits<double>::infinity();
  x2 = y2 = z2 = - numeric_limits<double>::infinity();

  for (int c = 0; c <= 1; ++c)
  for (int b = 0; b <= 1; ++b)
  for (int a = 0; a <= 1; ++a) {
    p = irtkPoint(x[a], y[b], z[c]);
    image->WorldToImage(p);
    if (p._x < x1) x1 = p._x;
    if (p._x > x2) x2 = p._x;
    if (p._y < y1) y1 = p._y;
    if (p._y > y2) y2 = p._y;
    if (p._z < z1) z1 = p._z;
    if (p._z > z2) z2 = p._z;
  }

  // Round to nearest voxel in image domain
  i1 = static_cast<int>(round(x1));
  i2 = static_cast<int>(round(x2));
  j1 = static_cast<int>(round(y1));
  j2 = static_cast<int>(round(y2));
  k1 = static_cast<int>(round(z1));
  k2 = static_cast<int>(round(z2));

  // When both indices are outside in opposite directions,
  // use the full range [0, N[. If they are both outside in
  // the same direction, the condition i1 <= i2 is false which
  // indicates that the bounding box is empty in this case
  i1 = (i1 < 0 ?  0 : (i1 >= image->X() ? image->X()     : i1));
  i2 = (i2 < 0 ? -1 : (i2 >= image->X() ? image->X() - 1 : i2));
  j1 = (j1 < 0 ?  0 : (j1 >= image->Y() ? image->Y()     : j1));
  j2 = (j2 < 0 ? -1 : (j2 >= image->Y() ? image->Y() - 1 : j2));
  k1 = (k1 < 0 ?  0 : (k1 >= image->Z() ? image->Z()     : k1));
  k2 = (k2 < 0 ? -1 : (k2 >= image->Z() ? image->Z() - 1 : k2));
  return i1 <= i2 && j1 <= j2 && k1 <= k2;
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::BoundingBox(const irtkImage *image, int cp,
                                                    int &i1, int &j1, int &k1, int &l1,
                                                    int &i2, int &j2, int &k2, int &l2,
                                                    double fraction) const
{
  // Calculate spatial bounding box in image coordinates
  bool bbvalid = BoundingBox(image, cp, i1, j1, k1, i2, j2, k2, fraction);

  // Calculate temporal bounding box
  double t1, t2;
  this->BoundingBox(cp, t1, t2, fraction);

  // Convert to image coordinates
  t1 = image->TimeToImage(t1);
  t2 = image->TimeToImage(t2);
  if (t2 < t1) swap(t1, t2);

  // Round to nearest voxel in image domain
  l1 = static_cast<int>(round(t1));
  l2 = static_cast<int>(round(t2));

  // When both indices are outside in opposite directions,
  // use the full range [0, N[. If they are both outside in
  // the same direction, the condition l1 <= l2 is false which
  // indicates that the bounding box is empty in this case
  l1 = (l1 < 0 ?  0 : (l1 >= image->T() ? image->T()     : l1));
  l2 = (l2 < 0 ? -1 : (l2 >= image->T() ? image->T() - 1 : l2));
  return bbvalid && l1 <= l2;
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::DOFBoundingBox(const irtkImage *image, int dof,
                                                       int &i1, int &j1, int &k1,
                                                       int &i2, int &j2, int &k2,
                                                       double fraction) const
{
  return BoundingBox(image, this->DOFToIndex(dof), i1, j1, k1, i2, j2, k2, fraction);
}

// =============================================================================
// Transformation parameters (DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::DOFGradientNorm(const double *gradient) const
{
  double norm, max = .0;
  int x, y, z;
  const int ncps = this->NumberOfCPs();
  for (int cp = 0; cp < ncps; ++cp) {
    this->IndexToDOFs(cp, x, y, z);
    norm = sqrt(gradient[x] * gradient[x] + gradient[y] * gradient[y] + gradient[z] * gradient[z]);
    if (norm > max) max = norm;
  }
  return max;
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::CopyFrom(const irtkTransformation *other)
{
  const irtkFreeFormTransformation *ffd;
  ffd = dynamic_cast<const irtkFreeFormTransformation *>(other);
  if (ffd && typeid(*ffd) == typeid(*this) && ffd->Attributes() == this->Attributes()) {
    return irtkTransformation::CopyFrom(other);
  } else {
    return false;
  }
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Put(int cp, const Vector &x)
{
  _CPImage(cp) = x;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Put(int i, int j, int k, int l,
                                            double x, double y, double z)
{
  _CPImage(i, j, k, l) = Vector(x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Put(int i, int j, int k,
                                            double x, double y, double z)
{
  Put(i, j, k, 0, x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Get(int cp, Vector &x) const
{
  x = _CPImage(cp);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Get(int i, int j, int k, int l,
                                            double &x, double &y, double &z) const
{
  const Vector &param = _CPImage(i, j, k, l);
  x = param._x, y = param._y, z = param._z;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Get(int i, int j, int k,
                                            double &x, double &y, double &z) const
{
  Get(i, j, k, 0, x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::PutStatus(int cp, const CPStatus &status)
{
  _CPStatus[0][0][0][cp] = status;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::PutStatus(int i, int j, int k, int l,
                                                  DOFStatus sx, DOFStatus sy, DOFStatus sz)
{
  CPStatus &status = _CPStatus[l][k][j][i];
  status._x = sx;
  status._y = sy;
  status._z = sz;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::PutStatus(int i, int j, int k,
                                                  DOFStatus sx, DOFStatus sy, DOFStatus sz)
{
  PutStatus(i, j, k, 0, sx, sy, sz);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GetStatus(int cp, CPStatus &status) const
{
  const CPStatus &s = _CPStatus[0][0][0][cp];
  status._x = s._x;
  status._y = s._y;
  status._z = s._z;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GetStatus(int i, int j, int k, int l,
                                                  DOFStatus &sx, DOFStatus &sy, DOFStatus &sz) const
{
  const CPStatus &status = _CPStatus[l][k][j][i];
  sx = status._x;
  sy = status._y;
  sz = status._z;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GetStatus(int i, int j, int k,
                                                  DOFStatus &sx, DOFStatus &sy, DOFStatus &sz) const
{
  GetStatus(i, j, k, 0, sx, sy, sz);
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::IsActive(int cp) const
{
  const CPStatus &status = _CPStatus[0][0][0][cp];
  return status._x == Active || status._y == Active || status._z == Active;
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::IsActive(int i, int j, int k, int l) const
{
  const CPStatus &status = _CPStatus[l][k][j][i];
  return status._x == Active || status._y == Active || status._z == Active;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GlobalTransform(double &, double &, double &, double, double) const
{
  // T_global(x) = x
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GlobalInverse(double &, double &, double &, double, double) const
{
  // T_global(x) = x
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Transform(double &x, double &y, double &z, double t, double t0) const
{
  this->LocalTransform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkFreeFormTransformation::Inverse(double &x, double &y, double &z, double t, double t0) const
{
  return this->LocalInverse(x, y, z, t, t0);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::JacobianToWorld(double &du, double &dv) const
{
  double dx = du * _matW2L(0, 0) + dv * _matW2L(1, 0);
  double dy = du * _matW2L(0, 1) + dv * _matW2L(1, 1);
  du = dx, dv = dy;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::JacobianToWorld(double &du, double &dv, double &dw) const
{
  double dx = du * _matW2L(0, 0) + dv * _matW2L(1, 0) + dw * _matW2L(2, 0);
  double dy = du * _matW2L(0, 1) + dv * _matW2L(1, 1) + dw * _matW2L(2, 1);
  double dz = du * _matW2L(0, 2) + dv * _matW2L(1, 2) + dw * _matW2L(2, 2);
  du = dx, dv = dy, dw = dz;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::JacobianToWorld(irtkMatrix &jac) const
{
  if (jac.Rows() == 2) {
    JacobianToWorld(jac(0, 0), jac(0, 1));
    JacobianToWorld(jac(1, 0), jac(1, 1));
  } else {
    JacobianToWorld(jac(0, 0), jac(0, 1), jac(0, 2));
    JacobianToWorld(jac(1, 0), jac(1, 1), jac(1, 2));
    JacobianToWorld(jac(2, 0), jac(2, 1), jac(2, 2));
  }
}

// -----------------------------------------------------------------------------
inline void
irtkFreeFormTransformation::HessianToWorld(double &duu, double &duv, double &dvv) const
{
  // The derivatives of the world to lattice coordinate transformation
  // w.r.t the world coordinates which are needed for the chain rule below
  const double &dudx = _matW2L(0, 0);
  const double &dudy = _matW2L(0, 1);
  const double &dvdx = _matW2L(1, 0);
  const double &dvdy = _matW2L(1, 1);
  // Expression computed here is transpose(R) * Hessian * R = transpose(Hessian * R) * R
  // where R is the 2x2 world to lattice reorientation and scaling matrix
  double du, dv, dxx, dxy, dyy;
  du  = duu * dudx + duv * dvdx;
  dv  = duv * dudx + dvv * dvdx;
  dxx = du  * dudx + dv  * dvdx;
  dxy = du  * dudy + dv  * dvdy;
  du  = duu * dudy + duv * dvdy;
  dv  = duv * dudy + dvv * dvdy;
  dyy = du  * dudy + dv  * dvdy;
  // Return computed derivatives
  duu = dxx, duv = dxy, dvv = dyy;
}

// -----------------------------------------------------------------------------
inline void
irtkFreeFormTransformation::HessianToWorld(double &duu, double &duv, double &duw,
                                                        double &dvv, double &dvw,
                                                                     double &dww) const
{
  // The derivatives of the world to lattice coordinate transformation
  // w.r.t the world coordinates which are needed for the chain rule below
  const double &dudx = _matW2L(0, 0);
  const double &dudy = _matW2L(0, 1);
  const double &dudz = _matW2L(0, 2);
  const double &dvdx = _matW2L(1, 0);
  const double &dvdy = _matW2L(1, 1);
  const double &dvdz = _matW2L(1, 2);
  const double &dwdx = _matW2L(2, 0);
  const double &dwdy = _matW2L(2, 1);
  const double &dwdz = _matW2L(2, 2);
  // Expression computed here is transpose(R) * Hessian * R = transpose(Hessian * R) * R
  // where R is the 3x3 world to lattice reorientation and scaling matrix
  double du, dv, dw, dxx, dxy, dxz, dyy, dyz, dzz;
  du  = duu * dudx + duv * dvdx + duw * dwdx;
  dv  = duv * dudx + dvv * dvdx + dvw * dwdx;
  dw  = duw * dudx + dvw * dvdx + dww * dwdx;
  dxx = du  * dudx + dv  * dvdx + dw  * dwdx;
  dxy = du  * dudy + dv  * dvdy + dw  * dwdy;
  dxz = du  * dudz + dv  * dvdz + dw  * dwdz;
  du  = duu * dudy + duv * dvdy + duw * dwdy;
  dv  = duv * dudy + dvv * dvdy + dvw * dwdy;
  dw  = duw * dudy + dvw * dvdy + dww * dwdy;
  dyy = du  * dudy + dv  * dvdy + dw  * dwdy;
  dyz = du  * dudz + dv  * dvdz + dw  * dwdz;
  du  = duu * dudz + duv * dvdz + duw * dwdz;
  dv  = duv * dudz + dvv * dvdz + dvw * dwdz;
  dw  = duw * dudz + dvw * dvdz + dww * dwdz;
  dzz = du  * dudz + dv  * dvdz + dw  * dwdz;
  // Return computed derivatives
  duu = dxx, duv = dxy, duw = dxz, dvv = dyy, dvw = dyz, dww = dzz;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::HessianToWorld(irtkMatrix &hessian) const
{
  if (hessian.Rows() == 2) {
    HessianToWorld(hessian(0, 0), hessian(0, 1), hessian(1, 1));
    hessian(1, 0) = hessian(0, 1);
  } else {
    HessianToWorld(hessian(0, 0), hessian(0, 1), hessian(0, 2),
                   hessian(1, 1), hessian(1, 2),
                   hessian(2, 2));
    hessian(1, 0) = hessian(0, 1);
    hessian(2, 0) = hessian(0, 2);
    hessian(2, 1) = hessian(1, 2);
  }
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::HessianToWorld(irtkMatrix hessian[3]) const
{
  HessianToWorld(hessian[0]);
  HessianToWorld(hessian[1]);
  HessianToWorld(hessian[2]);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::FFDJacobianWorld(irtkMatrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::FFDJacobianWorld: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GlobalJacobian(irtkMatrix &jac, double, double, double, double, double) const
{
  // T_global(x) = x
  jac.Initialize(3, 3);
  jac(0, 0) = 1;
  jac(1, 1) = 1;
  jac(2, 2) = 1;
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Jacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->LocalJacobian(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::GlobalHessian(irtkMatrix hessian[3], double, double, double, double, double) const
{
  // T_global(x) = x
  hessian[0].Initialize(3, 3);
  hessian[1].Initialize(3, 3);
  hessian[2].Initialize(3, 3);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::Hessian(irtkMatrix hessian[3], double x, double y, double z, double t, double t0) const
{
  this->LocalHessian(hessian, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::JacobianDOFs(irtkMatrix &jac, int cp, double x, double y, double z, double t, double t0) const
{
  int xdof, ydof, zdof;
  double dofgrad[3];

  // initialize 3x3 Jacobian matrix
  jac.Initialize(3, 3);
  // compute indices of control point DoFs
  this->IndexToDOFs(cp, xdof, ydof, zdof);
  // derivatives w.r.t. x component of control point
  if (this->GetStatus(xdof) == Active) {
    this->JacobianDOFs(dofgrad, xdof, x, y, z, t, t0);
    jac(0, 0) = dofgrad[0];
    jac(1, 0) = dofgrad[1];
    jac(2, 0) = dofgrad[2];
  }
  // derivatives w.r.t. y component of control point
  if (this->GetStatus(ydof) == Active) {
    this->JacobianDOFs(dofgrad, ydof, x, y, z, t, t0);
    jac(0, 1) = dofgrad[0];
    jac(1, 1) = dofgrad[1];
    jac(2, 1) = dofgrad[2];
  }
  // derivatives w.r.t. z component of control point
  if (this->GetStatus(zdof) == Active) {
    this->JacobianDOFs(dofgrad, zdof, x, y, z, t, t0);
    jac(0, 2) = dofgrad[0];
    jac(1, 2) = dofgrad[1];
    jac(2, 2) = dofgrad[2];
  }
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation::JacobianDOFs(irtkTransformationJacobian &jac, double x, double y, double z, double t, double t0) const
{
  int xdof, ydof, zdof;
  double dofgrad[3];

  for (int cp = 0; cp < NumberOfCPs(); ++cp) {
    // compute indices of control point DoFs
    this->IndexToDOFs(cp, xdof, ydof, zdof);
    // derivatives w.r.t. x component of control point
    if (this->GetStatus(xdof) == Active) {
      this->JacobianDOFs(dofgrad, xdof, x, y, z, t, t0);
      if (dofgrad[0] != .0 && dofgrad[1] != .0 && dofgrad[2] != .0) {
        irtkTransformationJacobian::ColumnType &xdofgrad = jac(xdof);
        xdofgrad._x = dofgrad[0];
        xdofgrad._y = dofgrad[1];
        xdofgrad._z = dofgrad[2];
      }
    }
    // derivatives w.r.t. y component of control point
    if (this->GetStatus(ydof) == Active) {
      this->JacobianDOFs(dofgrad, ydof, x, y, z, t, t0);
      if (dofgrad[0] != .0 && dofgrad[1] != .0 && dofgrad[2] != .0) {
        irtkTransformationJacobian::ColumnType &ydofgrad = jac(ydof);
        ydofgrad._x = dofgrad[0];
        ydofgrad._y = dofgrad[1];
        ydofgrad._z = dofgrad[2];
      }
    }
    // derivatives w.r.t. z component of control point
    if (this->GetStatus(zdof) == Active) {
      this->JacobianDOFs(dofgrad, zdof, x, y, z, t, t0);
      if (dofgrad[0] != .0 && dofgrad[1] != .0 && dofgrad[2] != .0) {
        irtkTransformationJacobian::ColumnType &zdofgrad = jac(zdof);
        zdofgrad._x = dofgrad[0];
        zdofgrad._y = dofgrad[1];
        zdofgrad._z = dofgrad[2];
      }
    }
  }
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkFreeFormTransformation::Bending3D(const irtkMatrix hessian[3])
{
  const irtkMatrix &hx = hessian[0];
  const irtkMatrix &hy = hessian[1];
  const irtkMatrix &hz = hessian[2];

  const double &x_ii = hx(0, 0);
  const double &x_ij = hx(0, 1);
  const double &x_ik = hx(0, 2);
  const double &x_jj = hx(1, 1);
  const double &x_jk = hx(1, 2);
  const double &x_kk = hx(2, 2);

  const double &y_ii = hy(0, 0);
  const double &y_ij = hy(0, 1);
  const double &y_ik = hy(0, 2);
  const double &y_jj = hy(1, 1);
  const double &y_jk = hy(1, 2);
  const double &y_kk = hy(2, 2);

  const double &z_ii = hz(0, 0);
  const double &z_ij = hz(0, 1);
  const double &z_ik = hz(0, 2);
  const double &z_jj = hz(1, 1);
  const double &z_jk = hz(1, 2);
  const double &z_kk = hz(2, 2);

  return         (  x_ii * x_ii + x_jj * x_jj + x_kk * x_kk
                  + y_ii * y_ii + y_jj * y_jj + y_kk * y_kk
                  + z_ii * z_ii + z_jj * z_jj + z_kk * z_kk)
         + 2.0 * (  x_ij * x_ij + x_ik * x_ik + x_jk * x_jk
                  + y_ij * y_ij + y_ik * y_ik + y_jk * y_jk
                  + z_ij * z_ij + z_ik * z_ik + z_jk * z_jk);
}


#endif

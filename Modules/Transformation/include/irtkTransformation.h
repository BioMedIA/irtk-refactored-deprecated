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

#ifndef _IRTKTRANSFORMATION_H
#define _IRTKTRANSFORMATION_H

#include <irtkBSpline.h>
#include <irtkGeometry.h>
#include <irtkImage.h>
#include <irtkIndent.h>
#include <irtkTransformationJacobian.h>


////////////////////////////////////////////////////////////////////////////////
// Auxiliary types, constants, and macros
////////////////////////////////////////////////////////////////////////////////

/// Enumeration of transformation types
///
/// Each transformation class has its own enumeration value which is written
/// to the transformation file right after the magic number. Different versions
/// are distinguished by different type IDs.
///
/// Enumeration values 50-60 were assigned to the refactored FFD types
/// which store the control point data in an instance of irtkGenericImage.
///
/// Enumeration values 70-80 were assigned to the refactored FFD types
/// which store also the method used to extrapolate control point coefficients
/// outside the finite discrete lattice on which the FFD is defined.
///
/// \attention Do not change the enumeration value of existing entries as then
///            already saved transformation files are not identified correctly.
///            To yet allow a better ordering of the entries, the enumeration
///            values are thus assigned explicitly, also to remind of this.
enum irtkTransformationType
{
  IRTKTRANSFORMATION_MAGIC                           = 815007,
  IRTKTRANSFORMATION_UNKNOWN                         =      0,
  // linear transformations
  IRTKTRANSFORMATION_HOMOGENEOUS                     =      1,
  IRTKTRANSFORMATION_RIGID                           =      2,
  IRTKTRANSFORMATION_SIMILARITY                      =     22,
  IRTKTRANSFORMATION_AFFINE                          =      3,
  IRTKTRANSFORMATION_QUATERNION                      =     11,
  IRTKTRANSFORMATION_HOMO_TEMPORAL                   =     30,
  IRTKTRANSFORMATION_RIGID_TEMPORAL                  =     31,
  IRTKTRANSFORMATION_AFFINE_TEMPORAL                 =     32,
  // linear FFD
  IRTKTRANSFORMATION_LINEAR_FFD_2D_v1                =     70,
  IRTKTRANSFORMATION_LINEAR_FFD_2D = IRTKTRANSFORMATION_LINEAR_FFD_2D_v1,
  IRTKTRANSFORMATION_LINEAR_FFD_3D_v1                =      5,
  IRTKTRANSFORMATION_LINEAR_FFD_3D_v2                =     13,
  IRTKTRANSFORMATION_LINEAR_FFD_3D_v3                =     51,
  IRTKTRANSFORMATION_LINEAR_FFD_3D_v4                =     71,
  IRTKTRANSFORMATION_LINEAR_FFD_3D = IRTKTRANSFORMATION_LINEAR_FFD_3D_v4,
  IRTKTRANSFORMATION_LINEAR_FFD_4D_v1                =     17,
  IRTKTRANSFORMATION_LINEAR_FFD_4D_v2                =     52,
  IRTKTRANSFORMATION_LINEAR_FFD_4D_v3                =     72,
  IRTKTRANSFORMATION_LINEAR_FFD_4D = IRTKTRANSFORMATION_LINEAR_FFD_4D_v3,
  IRTKTRANSFORMATION_LINEAR_FFD_SV_v1                =     73,
  IRTKTRANSFORMATION_LINEAR_FFD_SV = IRTKTRANSFORMATION_LINEAR_FFD_SV_v1,
  IRTKTRANSFORMATION_LINEAR_FFD_TD_v1                =     18,
  IRTKTRANSFORMATION_LINEAR_FFD_TD_v2                =     54,
  IRTKTRANSFORMATION_LINEAR_FFD_TD_v3                =     74,
  IRTKTRANSFORMATION_LINEAR_FFD_TD = IRTKTRANSFORMATION_LINEAR_FFD_TD_v3,
  // B-spline FFD
  IRTKTRANSFORMATION_BSPLINE_FFD_2D_v1               =     75,
  IRTKTRANSFORMATION_BSPLINE_FFD_3D_v1               =      4,
  IRTKTRANSFORMATION_BSPLINE_FFD_3D_v2               =     12,
  IRTKTRANSFORMATION_BSPLINE_FFD_3D_v3               =     56,
  IRTKTRANSFORMATION_BSPLINE_FFD_3D_v4               =     76,
  IRTKTRANSFORMATION_BSPLINE_FFD_3D = IRTKTRANSFORMATION_BSPLINE_FFD_3D_v4,
  IRTKTRANSFORMATION_BSPLINE_FFD_4D_v1               =     14,
  IRTKTRANSFORMATION_BSPLINE_FFD_4D_v2               =     57,
  IRTKTRANSFORMATION_BSPLINE_FFD_4D_v3               =     77,
  IRTKTRANSFORMATION_BSPLINE_FFD_4D = IRTKTRANSFORMATION_BSPLINE_FFD_4D_v3,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v1               =     16,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v2               =     23,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v3               =     24,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v4               =     25,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v5               =     27,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v6               =     58,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v7               =     65,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV_v8               =     78,
  IRTKTRANSFORMATION_BSPLINE_FFD_SV = IRTKTRANSFORMATION_BSPLINE_FFD_SV_v8,
  IRTKTRANSFORMATION_BSPLINE_FFD_TD_v1               =     15,
  IRTKTRANSFORMATION_BSPLINE_FFD_TD_v2               =     21,
  IRTKTRANSFORMATION_BSPLINE_FFD_TD_v3               =     59,
  IRTKTRANSFORMATION_BSPLINE_FFD_TD_v4               =     79,
  IRTKTRANSFORMATION_BSPLINE_FFD_TD = IRTKTRANSFORMATION_BSPLINE_FFD_TD_v4,
  IRTKTRANSFORMATION_EIGEN_FFD_3D_v1                 =      6,
  IRTKTRANSFORMATION_EIGEN_FFD_3D_v2                 =     60,
  IRTKTRANSFORMATION_EIGEN_FFD_3D_v3                 =     80,
  IRTKTRANSFORMATION_EIGEN_FFD_3D = IRTKTRANSFORMATION_EIGEN_FFD_3D_v3,
  IRTKTRANSFORMATION_PERIODIC_v1                     =     20, // obsolete
  IRTKTRANSFORMATION_PERIODIC = IRTKTRANSFORMATION_PERIODIC_v1,
  // "decorating" transformations
  IRTKTRANSFORMATION_BSPLINE_FFD_STATISTICAL         =     61,
  // composite transformations
  IRTKTRANSFORMATION_MFFD                            =      7,
  IRTKTRANSFORMATION_FLUID_v1                        =      8,
  IRTKTRANSFORMATION_FLUID_v2                        =     81,
  IRTKTRANSFORMATION_FLUID = IRTKTRANSFORMATION_FLUID_v2,
  IRTKTRANSFORMATION_MFFD_SV                         =     26,
  // others
  IRTKTRANSFORMATION_LATTICE_FFD                     =      9,
  IRTKTRANSFORMATION_MULTI_FRAME_LATTICE_FFD         =     10
};

/// Size of lookup table of 1D B-spline function values
const int FFDLOOKUPTABLESIZE = irtkBSplineFunction::LookupTableSize;

////////////////////////////////////////////////////////////////////////////////
// Abstract transformation class
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class for general transformations.
 *
 * This is the abstract base class which defines a common interface for all
 * transformations. Each derived class has to implement at least the abstract
 * methods and some of the virtual ones. Most other methods call these virtual
 * methods and should not be required to be overwritten in subclasses.
 *
 * The second time argument to the interface methods corresponds to the time
 * of the untransformed source image. It is only considered by some 3D+t
 * transformations, in particular those which parameterize the transformation
 * using a non-stationary velocity field.
 */

class irtkTransformation : public irtkObservable
{
  irtkAbstractMacro(irtkTransformation);

public:

  /// Type of transformation parameter value
  typedef double   DOFValue;

protected:

  /// Number of transformation parameters
  int _NumberOfDOFs;

  /// Value of each transformation parameter
  DOFValue *_Param;

  /// Status of each transformation parameter (Active or Passive)
  DOFStatus *_Status;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkTransformation(int = 0);

  /// Copy constructor
  irtkTransformation(const irtkTransformation &);

  /// Copy constructor
  irtkTransformation(const irtkTransformation &, int);

  /// Initialize transformation parameters
  void InitializeDOFs(int);

  /// Copy transformation parameters (DoFs) and their status
  void InitializeDOFs(const irtkTransformation &, int = -1);

public:

  /// Static constructor. This function returns a pointer to a concrete
  /// new transformation of the specified type (e.g., IRTKTRANSFORMATION_RIGID).
  static irtkTransformation *New(irtkTransformationType);
  static irtkTransformation *New(unsigned int type) {
    return New(static_cast<irtkTransformationType>(type));
  }

  /// Static constructor. This function returns a pointer to a concrete
  /// transformation by copying the transformation passed to it.
  static irtkTransformation *New(const irtkTransformation *);

  /// Static constructor. This function returns a pointer to a concrete
  /// transformation by reading the transformation parameters from a file
  /// and creating the appropriate transformation.
  static irtkTransformation *New(const char *);

  /// Default destructor.
  virtual ~irtkTransformation();

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const irtkTransformation *);

  /// Get number of transformation parameters
  virtual int NumberOfDOFs() const;

  /// Get number of active transformation parameters
  int NumberOfActiveDOFs() const;

  /// Get number of passive transformation parameters
  int NumberOfPassiveDOFs() const;

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

  /// Put value of transformation parameter
  virtual void Put(int, double);

  /// Put values of transformation parameters
  virtual void Put(const DOFValue *);

  /// Add change to transformation parameters
  virtual void Add(const DOFValue *);

  /// Update transformation parameters given parametric gradient
  virtual double Update(const DOFValue *);

  /// Get value of transformation parameter
  virtual double Get(int) const;

  /// Get values of transformation parameters
  virtual void Get(DOFValue *) const;

  /// Put status of transformation parameter
  virtual void PutStatus(int, DOFStatus);

  /// Get status of transformation parameter
  virtual DOFStatus GetStatus(int) const;

  /// Checks whether transformation depends on the same vector of parameters
  virtual bool HasSameDOFsAs(const irtkTransformation *) const;

  /// Checks whether the transformation is an identity mapping
  virtual bool IsIdentity() const;

  /// Gets the spatial bounding box for a transformation parameter in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  virtual bool DOFBoundingBox(const irtkImage *, int, int &, int &, int &,
                                                      int &, int &, int &, double = 1) const;

  /// Reset transformation
  virtual void Reset();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Evaluates RMS error of transformation compared to another
  double EvaluateRMSError(const irtkImageAttributes &, const irtkTransformation *) const;

  /// Evaluates RMS error of transformation compared to given displacement field
  ///
  /// This overloaded version of EvaluateRMSError is recommended for displacement
  /// fields defined on a regular lattice. It does not require memory for the
  /// explicit storage of the locations of each displacement vector. Moreover,
  /// if this transformation requires the caching of the displacements, the other
  /// overloads are either very slow or cannot be used.
  double EvaluateRMSError(const irtkImageAttributes &, double *, double *) const;

  /// Evaluates RMS error of transformation compared to given displacement field
  ///
  /// This overloaded version of EvaluateRMSError is recommended for displacement
  /// fields defined on a regular lattice. It does not require memory for the
  /// explicit storage of the locations of each displacement vector. Moreover,
  /// if this transformation requires the caching of the displacements, the other
  /// overloads are either very slow or cannot be used.
  double EvaluateRMSError(const irtkImageAttributes &, double *, double *, double *) const;

  /// Evaluates RMS error of transformation compared to displacement field
  double EvaluateRMSError(const double *, const double *, const double *, double,
                          double       *, double       *, double       *, int no) const;

  /// Evaluates RMS error of transformation compared to given displacement field
  double EvaluateRMSError(const double *,  const double *, const double *, const double *,
                          double       *, double        *, double       *, int no) const;

  /// Approximate another transformation and return approximation error
  virtual double Approximate(const irtkImageAttributes &, const irtkTransformation *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double Approximate(irtkGenericImage<double> &, int = 1, double = .0);

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
  virtual double ApproximateAsNew(const irtkImageAttributes &, const irtkTransformation *,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(irtkGenericImage<double> &, int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const irtkImageAttributes &, double *, double *, double *,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const double *, const double *, const double *,
                                  double *,       double *,       double *, int,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  virtual double ApproximateAsNew(const double *, const double *, const double *, const double *,
                                  double *,       double *,       double *,       int,
                                  int = 1, double = .0);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors and finds a gradient to minimize the L2 norm of the error.
  virtual void ApproximateGradient(const irtkImageAttributes &,
                                   const double *, const double *, const double *,
                                   double *, double = 1.0) const;

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors and finds a gradient to minimize the L2 norm of the error.
  virtual void ApproximateGradient(const double *, const double *, const double *,
                                   const double *, const double *, const double *, int,
                                   double *, double = 1.0) const;

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors and finds a gradient to minimize the L2 norm of the error.
  virtual void ApproximateGradient(const double *, const double *, const double *, const double *,
                                   const double *, const double *, const double *, int,
                                   double *, double = 1.0) const;

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors. It finds a gradient w.r.t. the transformation parameters
  /// which minimizes the L2 norm of the approximation error and adds it to the
  /// input gradient with the given weight.
  virtual void ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                                       const double *, const double *, const double *, int,
                                       double *, double = 1.0) const;

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other overloads
  using irtkObservable::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = -1) const = 0;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = -1) const = 0;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = -1) const = 0;

  /// Transforms a single point
  virtual void Transform(irtkPoint &, double = 0, double = -1) const;

  /// Transforms a set of points
  virtual void Transform(irtkPointSet &, double = 0, double = -1) const;

  /// Transforms a set of points
  virtual void Transform(int, double *, double *, double *, double = 0, double = -1) const;

  /// Transforms a set of points
  virtual void Transform(int, double *, double *, double *, const double *, double = -1) const;

  /// Transforms world coordinates of image voxels
  virtual void Transform(irtkWorldCoordsImage &, double = -1) const;

  /// Calculates the displacement of a single point using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point
  virtual void Displacement(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement at specified lattice points
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each point. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the lattice points.
  virtual void Displacement(const irtkImageAttributes &, double *, double *, double *) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(irtkGenericImage<double> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(irtkGenericImage<float> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(irtkGenericImage<double> &, double, double, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(irtkGenericImage<float> &, double, double, const irtkWorldCoordsImage * = NULL) const;

  /// Whether this transformation implements a more efficient update of a given
  /// displacement field given the desired change of a transformation parameter
  virtual bool CanModifyDisplacement(int = -1) const;

  /// Updates the displacement vectors for a whole image domain
  ///
  /// \param[in]     dof Transformation parameter.
  /// \param[in]     dv  Change of transformation parameter value.
  /// \param[in,out] dx  Displacement field to be updated.
  /// \param[in]     t   Time point of start point.
  /// \param[in]     t0  Time point of end point.
  /// \param[in]     i2w Pre-computed world coordinates.
  virtual void DisplacementAfterDOFChange(int dof, double dv,
                                          irtkGenericImage<double> &dx,
                                          double t, double t0 = -1,
                                          const irtkWorldCoordsImage *i2w = NULL) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(irtkPoint &, double = 0, double = -1) const;

  /// Transforms a set of points using the inverse of the transformation
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int Inverse(irtkPointSet &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the global transformation only
  virtual void GlobalInverseDisplacement(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the local transformation only
  virtual bool LocalInverseDisplacement(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the transformation
  virtual bool InverseDisplacement(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the inverse displacement at specified lattice points
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each point. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the lattice points.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(const irtkImageAttributes &, double *, double *, double *) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(irtkGenericImage<double> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(irtkGenericImage<float> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(irtkGenericImage<double> &, double, double, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(irtkGenericImage<float> &, double, double, const irtkWorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the determinant of the Jacobian of the global transformation w.r.t world coordinates
  virtual double GlobalJacobian(double, double, double, double = 0, double = -1) const;

  /// Calculates the determinant of the Jacobian of the local transformation w.r.t world coordinates
  virtual double LocalJacobian(double, double, double, double = 0, double = -1) const;

  /// Calculates the determinant of the Jacobian of the transformation w.r.t world coordinates
  virtual double Jacobian(double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the global transformation w.r.t world coordinates
  virtual void GlobalHessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t a transformation parameter
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = -1) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(irtkMatrix &, int, double, double, double, double = 0, double = -1) const;

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

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const irtkGenericImage<double> *, double *,
                          const irtkWorldCoordsImage *,
                          double = -1, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const irtkGenericImage<double> *, double *,
                          double = -1, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkGenericImage<double> **, int, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  const double * = NULL, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const irtkGenericImage<double> **, int, double *,
                          const irtkWorldCoordsImage *,
                          const double * = NULL, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  void ParametricGradient(const irtkGenericImage<double> **, int, double *,
                          const double * = NULL, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkPointSet &, const irtkVector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints information about the transformation
  virtual void Print(irtkIndent = 0) const = 0;

  /// Reads a transformation from a file
  virtual void Read(const char *);

  /// Writes a transformation to a file
  virtual void Write(const char *) const;

  /// Reads a transformation from a file stream
  virtual irtkCifstream &Read(irtkCifstream &);

  /// Writes a transformation to a file stream
  virtual irtkCofstream &Write(irtkCofstream &) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(irtkTransformationType) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual irtkCifstream &ReadDOFs(irtkCifstream &, irtkTransformationType);

  /// Writes transformation parameters to a file stream
  virtual irtkCofstream &WriteDOFs(irtkCofstream &) const;

public:

  // ---------------------------------------------------------------------------
  // Others

  /// Returns type ID corresponding to transformations of the named class
  static irtkTransformationType TypeOfClass(const char *);

  /// Returns type ID of the instantiated transformation class
  virtual irtkTransformationType TypeOfClass() const;

  /// Verifies that the transformation is well constructed
  /// according to class-specific rules
  virtual void Verify();

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for transformation implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define irtkAbstractTransformationMacro(name)                                  \
  irtkAbstractMacro(name)

// -----------------------------------------------------------------------------
#define irtkTransformationMacro(name)                                          \
  irtkObjectMacro(name)

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkTransformation::NumberOfDOFs() const
{
  return _NumberOfDOFs;
}

// -----------------------------------------------------------------------------
inline int irtkTransformation::NumberOfActiveDOFs() const
{
  int nactive = 0;
  for (int dof = 0; dof < this->NumberOfDOFs(); ++dof) {
    if (this->GetStatus(dof) == Active) ++nactive;
  }
  return nactive;
}

// -----------------------------------------------------------------------------
inline int irtkTransformation::NumberOfPassiveDOFs() const
{
  return this->NumberOfDOFs() - this->NumberOfActiveDOFs();
}

// -----------------------------------------------------------------------------
inline double irtkTransformation::DOFGradientNorm(const double *gradient) const
{
  double norm, max = .0;
  for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
    norm = fabs(gradient[dof]);
    if (norm > max) max = norm;
  }
  return max;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Put(int idx, double x)
{
  if (_Param[idx] != static_cast<DOFValue>(x)) {
    _Param[idx] = static_cast<DOFValue>(x);
    this->Changed(true);
  }
}

// -----------------------------------------------------------------------------
inline double irtkTransformation::Get(int idx) const
{
  return static_cast<double>(_Param[idx]);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Put(const DOFValue *x)
{
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    if (_Param[idx] != x[idx]) {
      _Param[idx] = x[idx];
      this->Changed(true);
    }
  }
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Add(const DOFValue *dx)
{
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    if (dx[idx] != .0) {
      _Param[idx] += dx[idx];
      this->Changed(true);
    }
  }
}

// -----------------------------------------------------------------------------
inline double irtkTransformation::Update(const DOFValue *dx)
{
  double delta, max_delta = .0;
  for (int idx = 0; idx < _NumberOfDOFs; ++idx) {
    _Param[idx] += dx[idx];
    delta = fabs(dx[idx]);
    if (delta > max_delta) max_delta = delta;
  }
  if (max_delta > .0) this->Changed(true);
  return max_delta;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Get(DOFValue *x) const
{
  memcpy(x, _Param, _NumberOfDOFs * sizeof(DOFValue));
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::PutStatus(int idx, DOFStatus s)
{
  _Status[idx] = s;
}

// -----------------------------------------------------------------------------
inline DOFStatus irtkTransformation::GetStatus(int idx) const
{
  return _Status[idx];
}

// -----------------------------------------------------------------------------
inline bool irtkTransformation::HasSameDOFsAs(const irtkTransformation *t) const
{
  return (_Param == t->_Param);
}

// -----------------------------------------------------------------------------
inline bool HaveSameDOFs(const irtkTransformation *t1, const irtkTransformation *t2)
{
  // If the virtual HasSameDOFsAs function was not overriden by the specialized
  // type of one transformation, it is assumed that it uses the _Param memory to
  // store its transformation parameters. However, the other transformation
  // might be a specialized type which does not make use of it. In this case,
  // even if the two _Param pointers do not reference the same memory, the two
  // transformations may still use the same parameters because the other
  // transformation directly or indirectly wraps this transformation. Therefore,
  // the check with the transformations exchanged and the boolean OR (not AND).
  return t1->HasSameDOFsAs(t2) || t2->HasSameDOFsAs(t1);
}

// -----------------------------------------------------------------------------
inline bool irtkTransformation::IsIdentity() const
{
  for (int i = 0; i < this->NumberOfDOFs(); ++i) {
    if (this->Get(i) != .0) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
inline bool irtkTransformation
::DOFBoundingBox(const irtkImage *image, int, int &i1, int &j1, int &k1,
                                              int &i2, int &j2, int &k2, double fraction) const
{
  i1 = j1 = k1 = 0;
  i2 = image->X() - 1, j2 = image->Y() -1, k2 = image->Z() - 1;
  return !image->IsEmpty();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkTransformation::RequiresCachingOfDisplacements() const
{
  return false;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Transform(irtkPoint &p, double t, double t0) const
{
  this->Transform(p._x, p._y, p._z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Transform(irtkPointSet &pset, double t, double t0) const
{
  for (int i = 0; i < pset.Size(); i++) this->Transform(pset(i), t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkTransformation::Inverse(irtkPoint &p, double t, double t0) const
{
  return this->Inverse(p._x, p._y, p._z, t, t0);
}

// -----------------------------------------------------------------------------
inline int irtkTransformation::Inverse(irtkPointSet &pset, double t, double t0) const
{
  int n = 0;
  for (int i = 0; i < pset.Size(); ++i) {
    if (!this->Inverse(pset(i), t, t0)) ++n;
  }
  return n;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::GlobalDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->GlobalTransform(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::LocalDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->LocalTransform(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Displacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->Transform(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::GlobalInverseDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->GlobalInverse(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline bool irtkTransformation::LocalInverseDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  bool ok = this->LocalInverse(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;

  return ok;
}

// -----------------------------------------------------------------------------
inline bool irtkTransformation::InverseDisplacement(double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  bool ok = this->Inverse(x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;

  return ok;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkTransformation::GlobalJacobian(irtkMatrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::GlobalJacobian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::LocalJacobian(irtkMatrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalJacobian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Jacobian(irtkMatrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::Jacobian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::DeriveJacobianWrtDOF(irtkMatrix &, int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::DeriveJacobianWrtDOF: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline double irtkTransformation::GlobalJacobian(double x, double y, double z, double t, double t0) const
{
  irtkMatrix jac(3, 3);
  this->GlobalJacobian(jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline double irtkTransformation::LocalJacobian(double x, double y, double z, double t, double t0) const
{
  irtkMatrix jac(3, 3);
  this->LocalJacobian(jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline double irtkTransformation::Jacobian(double x, double y, double z, double t, double t0) const
{
  irtkMatrix jac(3, 3);
  this->Jacobian(jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::GlobalHessian(irtkMatrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::GlobalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::LocalHessian(irtkMatrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::Hessian(irtkMatrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::Hessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation::JacobianDOFs(double [3], int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::JacobianDOFs: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation
::ParametricGradient(const irtkGenericImage<double> *in, double *out,
                     const irtkWorldCoordsImage *i2w, double t0, double w) const
{
  this->ParametricGradient(in, out, i2w, NULL, t0, w);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation
::ParametricGradient(const irtkGenericImage<double> *in, double *out, double t0, double w) const
{
  this->ParametricGradient(in, out, NULL, NULL, t0, w);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation
::ParametricGradient(const irtkGenericImage<double> **in, int n, double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     const double *t0, double w) const
{
  for (int i = 0; i < n; ++i) {
    this->ParametricGradient(in[i], out, i2w, wc, t0 ? t0[i] : 1.0, w);
  }
}

// -----------------------------------------------------------------------------
inline void irtkTransformation
::ParametricGradient(const irtkGenericImage<double> **in, int n, double *out,
                     const irtkWorldCoordsImage *i2w, const double *t0, double w) const
{
  this->ParametricGradient(in, n, out, i2w, NULL, t0, w);
}

// -----------------------------------------------------------------------------
inline void irtkTransformation
::ParametricGradient(const irtkGenericImage<double> **in, int n, double *out, const double *t0, double w) const
{
  this->ParametricGradient(in, n, out, NULL, NULL, t0, w);
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkTransformationType irtkTransformation::TypeOfClass() const
{
  return irtkTransformation::TypeOfClass(this->NameOfClass());
}

////////////////////////////////////////////////////////////////////////////////
// Auxiliary functions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Check whether a named file is an IRTK transformation file
///
/// \param[in] name File name.
///
/// \returns Whether the named file exists and stores an IRTK transformation.
bool IsTransformation(const char *name);

////////////////////////////////////////////////////////////////////////////////
// Transformations
////////////////////////////////////////////////////////////////////////////////

// Homogeneous transformations
#include <irtkHomogeneousTransformation.h>
#include <irtkRigidTransformation.h>
#include <irtkSimilarityTransformation.h>
#include <irtkAffineTransformation.h>

// Free-form transformations
#include <irtkFreeFormTransformation.h>
#include <irtkFreeFormTransformation3D.h>
#include <irtkFreeFormTransformation4D.h>

#include <irtkBSplineFreeFormTransformation3D.h>
#include <irtkBSplineFreeFormTransformation4D.h>
#include <irtkBSplineFreeFormTransformationSV.h>
#include <irtkBSplineFreeFormTransformationTD.h>
#include <irtkBSplineFreeFormTransformationStatistical.h>
#include <irtkLinearFreeFormTransformation3D.h>
#include <irtkLinearFreeFormTransformation4D.h>
#include <irtkLinearFreeFormTransformationTD.h>
#include <irtkEigenFreeFormTransformation.h>

// Composite transformations
#include <irtkMultiLevelTransformation.h>
#include <irtkMultiLevelFreeFormTransformation.h>
#include <irtkMultiLevelStationaryVelocityTransformation.h>
#include <irtkFluidFreeFormTransformation.h>

// Decorators (i.e., wrappers)
#include <irtkInverseAffineTransformation.h>
#include <irtkPartialAffineTransformation.h>
#include <irtkPartialBSplineFreeFormTransformationSV.h>
#include <irtkPartialMultiLevelStationaryVelocityTransformation.h>

// Image transformation filters
#include <irtkImageTransformation.h>
#include <irtkImageTransformation2.h>
#include <irtkImageHomogeneousTransformation.h>
//#include <irtkMultipleImageTransformation.h>

// Typedefs for backwards compatibility
typedef irtkBSplineFreeFormTransformation3D irtkBSplineFreeFormTransformation;
typedef irtkLinearFreeFormTransformation3D  irtkLinearFreeFormTransformation;


#endif

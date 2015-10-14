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

#ifndef _IRTKMULTILEVELTRANSFORMATION_H
#define _IRTKMULTILEVELTRANSFORMATION_H


const int MAX_TRANS = 200;


/**
 * Base class for multi-level transformations.
 *
 * This is the abstract base class which defines a common interface for all
 * multi-level transformations. Each derived class has to implement at least
 * the abstract methods and some of the virtual ones. Most other methods call
 * these virtual methods and should not be required to be overridden in
 * subclasses.
 */

class irtkMultiLevelTransformation : public irtkTransformation
{
  irtkAbstractTransformationMacro(irtkMultiLevelTransformation);

public:

  /// Type of local transformation status
  typedef DOFStatus FFDStatus;

protected:

  // ---------------------------------------------------------------------------
  // Data members

  /// Global transformation
  irtkAffineTransformation _GlobalTransformation;

  /// Local transformations
  irtkFreeFormTransformation *_LocalTransformation[MAX_TRANS];

  /// Status of local transformations
  FFDStatus _LocalTransformationStatus[MAX_TRANS];

  /// Number of local transformations
  int _NumberOfLevels;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkMultiLevelTransformation();

  /// Construct multi-level transformation given a rigid transformation
  irtkMultiLevelTransformation(const irtkRigidTransformation &);

  /// Construct multi-level transformation given an affine transformation
  irtkMultiLevelTransformation(const irtkAffineTransformation &);

  /// Copy constructor
  irtkMultiLevelTransformation(const irtkMultiLevelTransformation &);

public:

  /// Destructor
  virtual ~irtkMultiLevelTransformation();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const irtkImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const double *, const double *, const double *, const double *,
                             double *,       double *,       double *,       int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function also modifies the global transformation.
  ///       Use Reset and Approximate instead if this is not desired.
  virtual double ApproximateAsNew(const irtkImageAttributes &, double *, double *, double *,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function also modifies the global transformation.
  ///       Use Reset and Approximate instead if this is not desired.
  virtual double ApproximateAsNew(const double *, const double *, const double *,
                                  double *,       double *,       double *, int,
                                  int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function also modifies the global transformation.
  ///       Use Reset and Approximate instead if this is not desired.
  virtual double ApproximateAsNew(const double *, const double *, const double *, const double *,
                                  double *,       double *,       double *,       int,
                                  int = 1, double = .0);

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
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const irtkTransformation *);

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

  /// Get number of transformation parameters
  virtual int NumberOfDOFs() const;

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

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity() const;

  /// Reset transformation (does not remove local transformations)
  virtual void Reset();

  /// Reset transformation and remove all local transformations
  virtual void Clear();

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other overloads
  using irtkTransformation::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Levels

  /// Returns the number of levels
  virtual int NumberOfLevels() const;

  /// Returns the number of active levels
  virtual int NumberOfActiveLevels() const;

  /// Returns the number of passive levels
  virtual int NumberOfPassiveLevels() const;

  /// Returns the total number of control points
  virtual int NumberOfCPs(bool = false) const;

  /// Returns the total number of active control points on all active levels
  virtual int NumberOfActiveCPs() const;

  /// Gets global transformation
  virtual irtkAffineTransformation *GetGlobalTransformation();

  /// Get global transformation
  virtual const irtkAffineTransformation *GetGlobalTransformation() const;

  /// Gets local transformation
  virtual irtkFreeFormTransformation *GetLocalTransformation(int);

  /// Gets local transformation
  virtual const irtkFreeFormTransformation *GetLocalTransformation(int) const;

  /// Put local transformation and return pointer to previous one (needs to be deleted if not used)
  virtual irtkFreeFormTransformation *PutLocalTransformation(irtkFreeFormTransformation *, int);

  /// Push local transformation on stack (append transformation)
  virtual void PushLocalTransformation(irtkFreeFormTransformation *);

  /// Insert local transformation
  virtual void InsertLocalTransformation(irtkFreeFormTransformation *, int = 0);

  /// Pop local transformation from stack (remove last transformation)
  virtual irtkFreeFormTransformation *PopLocalTransformation();

  /// Remove local transformation and return the pointer (need to be deleted if not used)
  virtual irtkFreeFormTransformation *RemoveLocalTransformation(int = 0);

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Put status of local transformation
  virtual void LocalTransformationStatus(int, FFDStatus);

  /// Get status of local transformation
  virtual FFDStatus LocalTransformationStatus(int) const;

  /// Get whether local transformation is active
  virtual bool LocalTransformationIsActive(int) const;

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

protected:

  /// Checks whether a given transformation is supported as local transformation
  virtual void CheckTransformation(irtkFreeFormTransformation *) const;

  /// Helper function for MergeGlobalIntoLocalDisplacement
  void InterpolateGlobalDisplacement(irtkFreeFormTransformation *);

public:

  // ---------------------------------------------------------------------------
  // Point transformation

  // Do not hide methods of base class
  using irtkTransformation::Transform;
  using irtkTransformation::LocalDisplacement;
  using irtkTransformation::Displacement;
  using irtkTransformation::Inverse;
  using irtkTransformation::LocalInverseDisplacement;
  using irtkTransformation::InverseDisplacement;

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = -1) const = 0;

  /// Transforms a single point
  virtual void Transform(int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(int, int, irtkPoint &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(int, irtkPoint &, double = 0, double = -1) const;

  /// Transforms a set of points
  virtual void Transform(int, int, irtkPointSet &, double = 0, double = -1) const;

  /// Transforms a set of points
  virtual void Transform(int, irtkPointSet &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the local transformation component only
  virtual void LocalDisplacement(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the local transformation component only
  virtual void LocalDisplacement(int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point
  virtual void Displacement(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point
  virtual void Displacement(int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(int, int, irtkGenericImage<double> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(int, int, irtkGenericImage<float> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(int, irtkGenericImage<double> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  virtual void Displacement(int, irtkGenericImage<float> &, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the local transformation only
  virtual bool LocalInverseDisplacement(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the local transformation only
  virtual bool LocalInverseDisplacement(int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the transformation
  virtual bool InverseDisplacement(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the displacement of a single point using the inverse of the transformation
  virtual bool InverseDisplacement(int, double &, double &, double &, double = 0, double = -1) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(int, int, irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(int, int, irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(int, irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(int, irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points at which transformation is non-invertible.
  virtual int InverseDisplacement(irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  // Do not hide methods of base class
  using irtkTransformation::GlobalJacobian;
  using irtkTransformation::LocalJacobian;
  using irtkTransformation::Jacobian;

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(int, irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(int, int, irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(int, irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the determinant of the Jacobian of the local transformation w.r.t world coordinates
  virtual double LocalJacobian(int, double, double, double, double = 0, double = -1) const;

  /// Calculates the determinant of the Jacobian of the transformation w.r.t world coordinates
  virtual double Jacobian(int, int, double, double, double, double = 0, double = -1) const;

  /// Calculates the determinant of the Jacobian of the transformation w.r.t world coordinates
  virtual double Jacobian(int, double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the global transformation w.r.t world coordinates
  virtual void GlobalHessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(int, irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(int, int, irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(int, irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(irtkMatrix [3], double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = -1) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(irtkMatrix &, int, double, double, double, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Calculates the bending energy of the transformation
  virtual double BendingEnergy(int, int, double, double, double, double = 0, double = -1, bool = true) const;

  /// Calculates the bending energy of the transformation
  virtual double BendingEnergy(int, double, double, double, double = 0, double = -1, bool = true) const;

  /// Calculates the bending energy of the transformation
  virtual double BendingEnergy(double, double, double, double = 0, double = -1, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual irtkCifstream &ReadDOFs(irtkCifstream &, irtkTransformationType);

  /// Writes transformation parameters to a file stream
  virtual irtkCofstream &WriteDOFs(irtkCofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline int DOFIndexToLocalTransformation(const irtkMultiLevelTransformation *mffd, int  idx,
                                         const irtkFreeFormTransformation   *&ffd, int &dof)
{
  irtkAssert(idx >= 0, "DoF index must be positive");
  dof = idx;
  for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
    if (!mffd->LocalTransformationIsActive(l)) continue;
    ffd = mffd->GetLocalTransformation(l);
    if (dof < ffd->NumberOfDOFs()) return l;
    dof -= ffd->NumberOfDOFs();
  }
  ffd = NULL;
  irtkAssert(false, "DoF index out-of-bounds");
  return -1;
}

// -----------------------------------------------------------------------------
inline int DOFIndexToLocalTransformation(irtkMultiLevelTransformation *mffd, int  idx,
                                         irtkFreeFormTransformation   *&ffd, int &dof)
{
  irtkAssert(idx >= 0, "DoF index must be positive");
  dof = idx;
  for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
    if (!mffd->LocalTransformationIsActive(l)) continue;
    ffd = mffd->GetLocalTransformation(l);
    if (dof < ffd->NumberOfDOFs()) return l;
    dof -= ffd->NumberOfDOFs();
  }
  ffd = NULL;
  irtkAssert(false, "DoF index out-of-bounds");
  return -1;
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::NumberOfLevels() const
{
  return _NumberOfLevels;
}

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::NumberOfActiveLevels() const
{
  int n = 0;
  for (int i = 0; i < _NumberOfLevels; ++i) {
    if (_LocalTransformationStatus[i] == Active) ++n;
  }
  return n;
}

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::NumberOfPassiveLevels() const
{
  return _NumberOfLevels - this->NumberOfActiveLevels();
}

// -----------------------------------------------------------------------------
inline irtkAffineTransformation *irtkMultiLevelTransformation::GetGlobalTransformation()
{
  return &_GlobalTransformation;
}

// -----------------------------------------------------------------------------
inline const irtkAffineTransformation *irtkMultiLevelTransformation::GetGlobalTransformation() const
{
  return &_GlobalTransformation;
}

// -----------------------------------------------------------------------------
inline irtkFreeFormTransformation *irtkMultiLevelTransformation::GetLocalTransformation(int i)
{
  if (i < 0) i += _NumberOfLevels;
  return (0 <= i && i < _NumberOfLevels) ? _LocalTransformation[i] : NULL;
}

// -----------------------------------------------------------------------------
inline const irtkFreeFormTransformation *irtkMultiLevelTransformation::GetLocalTransformation(int i) const
{
  if (i < 0) i += _NumberOfLevels;
  return (0 <= i && i < _NumberOfLevels) ? _LocalTransformation[i] : NULL;
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalTransformationStatus(int i, FFDStatus status)
{
  if (i < 0) i += _NumberOfLevels;
  _LocalTransformationStatus[i] = status;
}

// -----------------------------------------------------------------------------
inline irtkMultiLevelTransformation::FFDStatus
irtkMultiLevelTransformation::LocalTransformationStatus(int i) const
{
  if (i < 0) i += _NumberOfLevels;
  return _LocalTransformationStatus[i];
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::LocalTransformationIsActive(int i) const
{
  if (i < 0) i += _NumberOfLevels;
  return (_LocalTransformationStatus[i] == Active);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::RequiresCachingOfDisplacements() const
{
  for (int i = 0; i < this->NumberOfLevels(); ++i) {
    if (this->GetLocalTransformation(i)->RequiresCachingOfDisplacements()) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::GlobalTransform(double &x, double &y, double &z, double t, double t0) const
{
  this->GetGlobalTransformation()->Transform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
// Note: The irtkTransformation interface requires that every transformation is
//       the sum of global and local transformations, i.e., T = T_global + T_local.
//       This allows a common interpretation of what LocalTransform does.
inline void irtkMultiLevelTransformation::LocalTransform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  if (m < 0) {
    // Compute y = x + (T(x) - T_global(x))
    double x1 = x, y1 = y, z1 = z;
    double x2 = x, y2 = y, z2 = z;
    this->GetGlobalTransformation()->Transform(x1, y1, z1);
    this->Transform(-1, n, x2, y2, z2, t, t0);
    x += (x2 - x1);
    y += (y2 - y1);
    z += (z2 - z1);
  } else {
    this->Transform(m, n, x, y, z, t, t0);
  }
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalTransform(int n, double &x, double &y, double &z, double t, double t0) const
{
  this->LocalTransform(0, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  this->LocalTransform(-1, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Transform(int n, double &x, double &y, double &z, double t, double t0) const
{
  this->Transform(-1, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Transform(double &x, double &y, double &z, double t, double t0) const
{
  this->Transform(-1, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Transform(int l, int n, irtkPoint &p, double t, double t0) const
{
  this->Transform(l, n, p._x, p._y, p._z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Transform(int n, irtkPoint &p, double t, double t0) const
{
  this->Transform(0, n, p, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Transform(int l, int n, irtkPointSet &pset, double t, double t0) const
{
  for (int i = 0; i < pset.Size(); i++) this->Transform(l, n, pset(i), t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Transform(int n, irtkPointSet &pset, double t, double t0) const
{
  this->Transform(0, n, pset, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalDisplacement(int l, int n, double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->LocalTransform(l, n, x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalDisplacement(int n, double &x, double &y, double &z, double t, double t0) const
{
  this->LocalDisplacement(0, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int l, int n, double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  this->Transform(l, n, x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int n, double &x, double &y, double &z, double t, double t0) const
{
  this->Displacement(-1, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int l, int n, irtkGenericImage<double> &disp, double t0, const irtkWorldCoordsImage *i2w) const
{
  disp = .0;
  this->Displacement(l, n, disp, disp.GetTOrigin(), t0, i2w);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int l, int n, irtkGenericImage<float> &disp, double t0, const irtkWorldCoordsImage *i2w) const
{
  disp = .0f;
  this->Displacement(l, n, disp, disp.GetTOrigin(), t0, i2w);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int n, irtkGenericImage<double> &disp, double t0, const irtkWorldCoordsImage *i2w) const
{
  this->Displacement(-1, n, disp, t0, i2w);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int n, irtkGenericImage<float> &disp, double t0, const irtkWorldCoordsImage *i2w) const
{
  this->Displacement(-1, n, disp, t0, i2w);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int n, irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  this->Displacement(-1, n, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(int n, irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  this->Displacement(-1, n, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  this->Displacement(-1, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Displacement(irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  this->Displacement(-1, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::GlobalInverse(double &x, double &y, double &z, double t, double t0) const
{
  this->GetGlobalTransformation()->Inverse(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
// Note: The irtkTransformation interface requires that every transformation is
//       the sum of global and local transformations, i.e., T = T_global + T_local.
//       This allows a common interpretation of what LocalTransform does.
inline bool irtkMultiLevelTransformation::LocalInverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  if (m < 0) {
    // Compute y = x + (inv T(x) - inv T_global(x))
    double x1 = x, y1 = y, z1 = z;
    double x2 = x, y2 = y, z2 = z;
    this->GetGlobalTransformation()->Inverse(x1, y1, z1);
    bool ok = this->Inverse(-1, n, x2, y2, z2, t, t0);
    x += (x2 - x1);
    y += (y2 - y1);
    z += (z2 - z1);
    return ok;
  } else {
    return this->Inverse(m, n, x, y, z, t, t0);
  }
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::LocalInverse(int n, double &x, double &y, double &z, double t, double t0) const
{
  return this->LocalInverse(0, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  return this->LocalInverse(-1, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::Inverse(int n, double &x, double &y, double &z, double t, double t0) const
{
  return this->Inverse(-1, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::Inverse(double &x, double &y, double &z, double t, double t0) const
{
  return this->Inverse(-1, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::LocalInverseDisplacement(int l, int n, double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  bool ok = this->LocalInverse(l, n, x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;

  return ok;
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::LocalInverseDisplacement(int n, double &x, double &y, double &z, double t, double t0) const
{
  return this->LocalInverseDisplacement(0, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::InverseDisplacement(int l, int n, double &x, double &y, double &z, double t, double t0) const
{
  const double u = x;
  const double v = y;
  const double w = z;

  bool ok = this->Inverse(l, n, x, y, z, t, t0);

  x -= u;
  y -= v;
  z -= w;

  return ok;
}

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelTransformation::InverseDisplacement(int n, double &x, double &y, double &z, double t, double t0) const
{
  return this->InverseDisplacement(-1, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::InverseDisplacement(int n, irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  return this->InverseDisplacement(-1, n, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::InverseDisplacement(int n, irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  return this->InverseDisplacement(-1, n, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::InverseDisplacement(irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  return this->InverseDisplacement(-1, disp, t, t0, wc);
}

// -----------------------------------------------------------------------------
inline int irtkMultiLevelTransformation::InverseDisplacement(irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  return this->InverseDisplacement(-1, disp, t, t0, wc);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Jacobian(int, int, irtkMatrix &, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::Jacobian(int, int, ...): Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Jacobian(int n, irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->Jacobian(-1, n, jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Jacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->Jacobian(-1, -1, jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::GlobalJacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->GetGlobalTransformation()->Jacobian(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalJacobian(int n, irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->Jacobian(0, n, jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->Jacobian(0, -1, jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline double irtkMultiLevelTransformation::Jacobian(int m, int n, double x, double y, double z, double t, double t0) const
{
  irtkMatrix jac(3, 3);
  this->Jacobian(m, n, jac, x, y, z, t, t0);
  return jac.Det3x3();
}

// -----------------------------------------------------------------------------
inline double irtkMultiLevelTransformation::Jacobian(int n, double x, double y, double z, double t, double t0) const
{
  return this->Jacobian(-1, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline double irtkMultiLevelTransformation::LocalJacobian(int n, double x, double y, double z, double t, double t0) const
{
  return this->Jacobian(0, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Hessian(int, int, irtkMatrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::Hessian(int, int,...): Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Hessian(int n, irtkMatrix hessian[3], double x, double y, double z, double t, double t0) const
{
  this->Hessian(-1, n, hessian, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::Hessian(irtkMatrix hessian[3], double x, double y, double z, double t, double t0) const
{
  this->Hessian(-1, -1, hessian, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::GlobalHessian(irtkMatrix hessian[3], double x, double y, double z, double t, double t0) const
{
  this->GetGlobalTransformation()->Hessian(hessian, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalHessian(int n, irtkMatrix hessian[3], double x, double y, double z, double t, double t0) const
{
  this->Hessian(0, n, hessian, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::LocalHessian(irtkMatrix hessian[3], double x, double y, double z, double t, double t0) const
{
  this->Hessian(0, -1, hessian, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkMultiLevelTransformation::DeriveJacobianWrtDOF(irtkMatrix &jac, int dof, double x, double y, double z, double t, double t0) const
{
  cerr << this->NameOfClass() << "::DeriveJacobianWrtDOF: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkMultiLevelTransformation::BendingEnergy(int m, int n, double x, double y, double z, double t, double t0, bool wrt_world) const
{
  double bending = .0;
  const int l1 = (m < 0 ? 0 : m);
  const int l2 = (n < 0 || n > this->NumberOfLevels() ? this->NumberOfLevels() : n);
  for (int l = l1; l < l2; ++l) {
    bending += this->GetLocalTransformation(l)->BendingEnergy(x, y, z, t, t0, wrt_world);
  }
  return bending;
}

// -----------------------------------------------------------------------------
inline double irtkMultiLevelTransformation::BendingEnergy(int n, double x, double y, double z, double t, double t0, bool wrt_world) const
{
  return this->BendingEnergy(0, n, x, y, z, t, t0, wrt_world);
}

// -----------------------------------------------------------------------------
inline double irtkMultiLevelTransformation::BendingEnergy(double x, double y, double z, double t, double t0, bool wrt_world) const
{
  return this->BendingEnergy(0, -1, x, y, z, t, t0, wrt_world);
}


#endif

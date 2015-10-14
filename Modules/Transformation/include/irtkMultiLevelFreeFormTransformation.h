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

#ifndef _IRTKMULTILEVELFREEFORMTRANSFORMATION_H
#define _IRTKMULTILEVELFREEFORMTRANSFORMATION_H


/**
 * Class for multi-level FFD where global and local transformations are summed up.
 *
 * T_mffd(x) = T_global(x) + sum_i T_local^i(x)
 */

class irtkMultiLevelFreeFormTransformation : public irtkMultiLevelTransformation
{
  irtkTransformationMacro(irtkMultiLevelFreeFormTransformation);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  irtkMultiLevelFreeFormTransformation();

  /// Construct multi-level transformation given a rigid transformation
  irtkMultiLevelFreeFormTransformation(const irtkRigidTransformation &);

  /// Construct multi-level transformation given an affine transformation
  irtkMultiLevelFreeFormTransformation(const irtkAffineTransformation &);

  /// Copy constructor
  irtkMultiLevelFreeFormTransformation(const irtkMultiLevelFreeFormTransformation &);

  /// Destructor
  virtual ~irtkMultiLevelFreeFormTransformation();

  // ---------------------------------------------------------------------------
  // Levels

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

  // ---------------------------------------------------------------------------
  // Bounding box

  /// Gets the spatial bounding box for a transformation parameter in image coordinates.
  /// The last parameter specifies what fraction of the bounding box to return.
  /// The default is 1 which equals 100% of the bounding box.
  bool DOFBoundingBox(const irtkImage *, int, int &, int &, int &,
                                              int &, int &, int &, double = 1) const;

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and passive levels.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const irtkImageAttributes &, double *, double *, double *,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and passive levels.
  ///       Use ApproximateAsNew to also approximate a new global transformation.
  virtual double Approximate(const double *, const double *, const double *,
                             double *,       double *,       double *, int,
                             int = 1, double = .0);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds a transformation which approximates these
  /// displacements. After approximation, the displacements are replaced by the
  /// residual displacement errors at the points.
  ///
  /// \note This function does not change the global transformation and passive levels.
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
  // Point transformation

  // Do not hide base class methods
  using irtkMultiLevelTransformation::LocalTransform;
  using irtkMultiLevelTransformation::Transform;
  using irtkMultiLevelTransformation::Displacement;
  using irtkMultiLevelTransformation::InverseDisplacement;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(int, int, double &, double &, double &, double = 0, double = -1) const;

  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = -1) const;

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

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points for which transformation is non-invertible.
  virtual int InverseDisplacement(int, int, irtkGenericImage<double> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Number of points for which transformation is non-invertible.
  virtual int InverseDisplacement(int, int, irtkGenericImage<float> &, double, double = -1, const irtkWorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives

  // Do not hide base class methods
  using irtkMultiLevelTransformation::LocalJacobian;
  using irtkMultiLevelTransformation::Jacobian;
  using irtkMultiLevelTransformation::LocalHessian;
  using irtkMultiLevelTransformation::Hessian;
  using irtkMultiLevelTransformation::DeriveJacobianWrtDOF;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(int, int, irtkMatrix &, double, double, double, double = 0, double = -1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(int, int, irtkMatrix [3], double, double, double, double = 0, double = -1) const;

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
                                  const irtkWorldCoordsImage * = NULL,
                                  const irtkWorldCoordsImage * = NULL,
                                  double = -1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkPointSet &, const irtkVector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

  // ---------------------------------------------------------------------------
  // Backwards compatibility

  /// Bending of multi-level free-form deformation
  double Bending(double, double, double) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Bounding box
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkMultiLevelFreeFormTransformation
::DOFBoundingBox(const irtkImage *image, int dof, int &i1, int &j1, int &k1,
                                                  int &i2, int &j2, int &k2, double fraction) const
{
  const irtkFreeFormTransformation *ffd;
  DOFIndexToLocalTransformation(this, dof, ffd, dof);
  return ffd->DOFBoundingBox(image, dof, i1, j1, k1, i2, j2, k2, fraction);
}


#endif

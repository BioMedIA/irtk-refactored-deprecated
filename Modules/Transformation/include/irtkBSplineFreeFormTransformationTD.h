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

#ifndef _IRTKBSPLINEFREEFORMTRANSFORMATIONTD_H
#define _IRTKBSPLINEFREEFORMTRANSFORMATIONTD_H


/**
 * Temporal diffeomorphic free-form transformation.
 *
 * This class implements a free-form transformation which is represented
 * by a time-varying velocity field (3D+t). The 3D displacement field at a
 * specific time is obtained by integrating the velocity field starting at the
 * reference time point. The integration steps are adjusted if necessary in
 * order to ensure that the resulting spatial transformation is diffeomorphic.
 *
 * For more details about the implementation see De Craene et al. (2012).
 * Temporal diffeomorphic free-form deformation: application to motion and
 * strain estimation from 3D echocardiography.
 * Medical image analysis, 16(2), 427, 2012. doi:10.1016/j.media.2011.10.006
 */

class irtkBSplineFreeFormTransformationTD : public irtkBSplineFreeFormTransformation4D
{
  irtkTransformationMacro(irtkBSplineFreeFormTransformationTD);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Numerical integration method
  irtkPublicAttributeMacro(FFDIntegrationMethod, IntegrationMethod);

  /// Minimum length of temporal steps (in ms)
  irtkPublicAttributeMacro(double, MinTimeStep);

  /// Maximum length of temporal steps (in ms)
  irtkPublicAttributeMacro(double, MaxTimeStep);

  /// Local integration error tolerance (in mm)
  irtkPublicAttributeMacro(double, Tolerance);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  irtkBSplineFreeFormTransformationTD();

  /// Construct free-form transformation for given image domain and lattice spacing
  explicit irtkBSplineFreeFormTransformationTD(const irtkImageAttributes &,
                                               double = -1, double = -1, double = -1, double = -1);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit irtkBSplineFreeFormTransformationTD(const irtkBaseImage &,
                                               double, double, double, double);

  /// Copy constructor
  irtkBSplineFreeFormTransformationTD(const irtkBSplineFreeFormTransformationTD &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformationTD();

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  using irtkBSplineFreeFormTransformation4D::ApproximateAsNew;

  /// Approximate displacements: This function takes a set of 3D displacement fields
  /// and corresponding time point and temporal interval. Given these inputs,
  /// it finds a time-varying velocity field which approximates these displacements.
  virtual void ApproximateDOFs(const irtkGenericImage<double> * const *,
                               const double *, const double *, int,
                               bool = false, int = 3, int = 8);

  /// Approximates displacements: This function takes a set of points and a set
  /// of displacements and finds a !new! FFD which approximates these displacements.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors. It finds a gradient w.r.t. the transformation parameters
  /// which minimizes the L2 norm of the approximation error and adds it to the
  /// input gradient with the given weight.
  virtual void ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                                       const double *, const double *, const double *, int,
                                       double *, double = 1.0) const;

  /// Approximate displacements: This function takes a set of 3D displacement fields
  /// and corresponding time point and temporal interval. Given these inputs,
  /// it finds a time-varying velocity field which approximates these displacements.
  /// The displacements are replaced by the residual displacements of the newly
  /// approximated transformation. Returns the approximation error of the resulting FFD.
  virtual double ApproximateAsNew(irtkGenericImage<double> **,
                                  const double *, const double *, int,
                                  bool = false, int = 3, int = 8);

  /// Interpolates displacements: This function takes a set of displacements defined at
  /// the control points and finds a time-varying velocity field such that the
  /// resulting transformation interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double *);

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation
  using irtkBSplineFreeFormTransformation4D::Displacement;

  /// Transforms a single point
  virtual void LocalTransform(double &, double &, double &, double, double) const;

  /// Transforms a single point using the inverse transformation
  virtual bool LocalInverse(double &, double &, double &, double, double) const;

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

  // ---------------------------------------------------------------------------
  // Derivatives
  using irtkBSplineFreeFormTransformation4D::JacobianDOFs;
  using irtkBSplineFreeFormTransformation4D::ParametricGradient;

  /// Calculates the Jacobian of the (local) transformation w.r.t world coordinates
  /// and transforms the given point at the same time
  virtual void TransformAndJacobian(irtkMatrix &, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  /// of the specified control point and transforms the given point at the same time
  virtual void TransformAndJacobianDOFs(irtkMatrix &, int, int, int, int, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  /// of the specified control point and transforms the given point at the same time
  virtual void TransformAndJacobianDOFs(irtkMatrix &, int, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t all transformation parameters
  /// and transforms the given point at the same time
  virtual void TransformAndJacobianDOFs(irtkTransformationJacobian &, double &, double &, double &, double, double) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double, double) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(irtkMatrix [3], double, double, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(irtkMatrix &, int, int, int, int, double, double, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, int, int, int, double, double, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(irtkTransformationJacobian &, double, double, double, double, double) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const irtkGenericImage<double> *, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  double = 1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkPointSet &, const irtkVector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

  using irtkBSplineFreeFormTransformation4D::Read;

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

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

  /// Invert the transformation
  virtual void Invert();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////
  
// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkBSplineFreeFormTransformationTD
::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  this->LocalTransform(x, y, z, t, t0);
  return true;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationTD
::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  this->TransformAndJacobian(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationTD
::LocalHessian(irtkMatrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationTD
::JacobianDOFs(irtkMatrix &jac, int    i, int    j, int    k, int    l,
                                double x, double y, double z, double t, double t0) const
{
  this->TransformAndJacobianDOFs(jac, i, j, k, l, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationTD
::JacobianDOFs(irtkTransformationJacobian &jac, double x, double y, double z, double t, double t0) const
{
  jac.Clear();
  this->TransformAndJacobianDOFs(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationTD
::JacobianDOFs(double [3], int, int, int, int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::JacobianDOFs: Jacobian is full symmetric 3x3 matrix, not only a diagonal matrix" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationTD
::ParametricGradient(const irtkGenericImage<double> *in, double *out,
                     const irtkGenericImage<double> *i2w, const irtkWorldCoordsImage *wc,
                     double t0, double w) const
{
  // Use general implementation provided by FFD base class, not the one
  // of irtkFreeFormTransformation4D which is only valid for FFDs that are
  // parameterized by displacements instead of velocities.
  irtkFreeFormTransformation::ParametricGradient(in, out, i2w, wc, t0, w);
}

// -----------------------------------------------------------------------------
/*
void irtkBSplineFreeFormTransformationTD
::ParametricGradient(const irtkGenericImage<double> **in, int n, double *out,
                     const irtkGenericImage<double> *i2w, const irtkWorldCoordsImage *wc,
                     double *t0, double w) const
{
  // TODO: Can be more efficient by re-using already computed trajectories
  //       especially when different source images are frames of a temporal
  //       sequence. See irtkImageTDFFDRegistration.
}
*/


#endif

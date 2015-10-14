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

#ifndef _IRTKPARTIALBSPLINEFREEFORMTRANSFORMATIONSV_H
#define _IRTKPARTIALBSPLINEFREEFORMTRANSFORMATIONSV_H


/**
 * Decorator for SV FFD transformation.
 *
 * This decorator wraps a SV FFD transformation but defines it's own attributes
 * regarding integration and whether or not to invert the transformation by
 * default. As the SV FFD is parameterized by stationary velocities, inverting
 * the transformation is simply achieved by a negative upper integration limit.
 * Two instances of this decorator are in particular used by
 * irtkSymmetricBSplineFreeFormTransformation for the two half transformations
 * to deform both input images.
 */

class irtkPartialBSplineFreeFormTransformationSV : public irtkTransformation
{
  irtkTransformationMacro(irtkPartialBSplineFreeFormTransformationSV);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Pointer to decorated transformation
  irtkPublicAggregateMacro(irtkBSplineFreeFormTransformationSV, Transformation);

  /// Fraction of decorated transformation, negative value corresponds to inverse
  irtkPublicAttributeMacro(double, Fraction);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkPartialBSplineFreeFormTransformationSV(irtkBSplineFreeFormTransformationSV * = NULL, double = 1.0);

  /// Copy constructor
  irtkPartialBSplineFreeFormTransformationSV(const irtkPartialBSplineFreeFormTransformationSV &);

  /// Destructor
  virtual ~irtkPartialBSplineFreeFormTransformationSV();

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const irtkTransformation *);

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

  /// Checks whether the transformation is an identity mapping
  virtual bool IsIdentity() const;

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other overloads
  using irtkTransformation::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation

private:

  /// Partial upper integration limit used given the temporal origin of both target
  /// and source images. If both images have the same temporal origin, the sign
  /// of the fraction determines the direction of integration. Otherwise, if the
  /// temporal origin of the two images differs, the signed difference between
  /// these determines the direction of integration.
  double UpperIntegrationLimit(double t, double t0) const;

public:

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(double &, double &, double &, double = 0, double = 1) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(irtkGenericImage<double> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(irtkGenericImage<float> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(irtkGenericImage<double> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(irtkGenericImage<float> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using irtkTransformation::ParametricGradient;
  using irtkTransformation::GlobalJacobian;
  using irtkTransformation::LocalJacobian;
  using irtkTransformation::Jacobian;

  /// Calculates the Jacobian of the global transformation w.r.t world coordinates
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0, double = 1) const;

  /// Calculates the Jacobian of the local transformation w.r.t world coordinates
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0, double = 1) const;

  /// Calculates the Jacobian of the transformation w.r.t world coordinates
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0, double = 1) const;

  /// Calculates the Hessian for each component of the global transformation w.r.t world coordinates
  virtual void GlobalHessian(irtkMatrix [3], double, double, double, double = 0, double = 1) const;

  /// Calculates the Hessian for each component of the local transformation w.r.t world coordinates
  virtual void LocalHessian(irtkMatrix [3], double, double, double, double = 0, double = 1) const;

  /// Calculates the Hessian for each component of the transformation w.r.t world coordinates
  virtual void Hessian(irtkMatrix [3], double, double, double, double = 0, double = 1) const;

  /// Calculates the Jacobian of the transformation w.r.t the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = 1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkGenericImage<double> *, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  double = 1, double = 1) const;


  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkGenericImage<double> **, int, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  const double * = NULL, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O
  using irtkTransformation::Read;
  using irtkTransformation::Write;

  /// Prints information about the transformation
  virtual void Print(irtkIndent = 0) const;

  /// Reads transformation from a file stream
  virtual irtkCifstream &Read(irtkCifstream &);

  /// Writes transformation to a file stream
  virtual irtkCofstream &Write(irtkCofstream &) const;

};


#endif

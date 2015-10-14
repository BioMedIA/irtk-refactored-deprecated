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

#ifndef _IRTKBSPLINEFREEFORMTRANSFORMATIONSTATISTICAL_H
#define _IRTKBSPLINEFREEFORMTRANSFORMATIONSTATISTICAL_H

#include <irtkImageFunction.h>


/**
 * Class for free-form transformations based on tensor product B-splines.
 *
 * This class implements 3D statistical free-form transformation using B-splines.
 *
 */

class irtkBSplineFreeFormTransformationStatistical : public irtkBSplineFreeFormTransformation3D
{
  irtkTransformationMacro(irtkBSplineFreeFormTransformationStatistical);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Basis vectors (colums of the irtkMatrix object)
  irtkReadOnlyAttributeMacro(irtkMatrix, BasisVectors);

  /// Mean vector
  irtkReadOnlyAttributeMacro(irtkVector, MeanVector);

  /// Name of file from which statistical deformation model was read
  irtkAttributeMacro(string, ModelFile);

public:

  /// Default constructor
  irtkBSplineFreeFormTransformationStatistical();

  /// Constructor based on a basis vectors matrix and a mean vector
  irtkBSplineFreeFormTransformationStatistical(const irtkImageAttributes &,
                                               irtkVector3D<DOFStatus> ****,
                                               const irtkMatrix &,
                                               const irtkVector &);

  /// Copy Constructor
  irtkBSplineFreeFormTransformationStatistical(const irtkBSplineFreeFormTransformationStatistical &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformationStatistical();

  // Import other Initialize overloads
  using irtkBSplineFreeFormTransformation3D::Initialize;

  /// Initialize free-form transformation
  virtual void Initialize(const irtkImageAttributes &);

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

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

  /// Interpolates displacements: This function takes a set of displacements defined
  /// at the control points and finds a FFD which interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double * = NULL);

  // ---------------------------------------------------------------------------
  // Lattice

  using irtkBSplineFreeFormTransformation3D::CropPadPassiveCPs;

  /// Crop/pad lattice to discard passive control points at the boundary,
  /// keeping only a layer of passive control points of given width.
  /// The DoF values of passive control points are optionally reset to zero.
  virtual bool CropPadPassiveCPs(int, int, int = 0, int = 0, bool = false);

  // ---------------------------------------------------------------------------
  // Transformation parameters (DOFs)

  /// Get norm of the gradient vector
  virtual double DOFGradientNorm(const double *) const;

  /// Puts a transformation parameter
  virtual void Put(int, DOFValue);

  /// Puts transformation parameters
  virtual void Put(const DOFValue *);

  /// Add change to transformation parameters
  virtual void Add(const DOFValue *);

  // ---------------------------------------------------------------------------
  // Update

  /// Update control point displacements after change of parameters
  void UpdateCPs();

  /// Update parameters after change of control point displacements
  void UpdateDOFs();

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  using irtkBSplineFreeFormTransformation3D::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Derivatives

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const irtkGenericImage<double> *, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  double = 1, double = 1) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Calculates the gradient of the bending energy w.r.t the transformation parameters
  virtual void BendingEnergyGradient(double *, double = 1, bool = false, bool = true) const;

  // ---------------------------------------------------------------------------
  // I/O

  using irtkBSplineFreeFormTransformation3D::Write;

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

  /// Writes a transformation to a file stream
  virtual irtkCofstream &Write(irtkCofstream &) const;

  /// Reads statistical deformation model from a file
  virtual void ReadSDM(const char *);

  /// Writes statistical deformation model to a file
  virtual void WriteSDM(const char *);

  // ---------------------------------------------------------------------------
  // Other

  /// Verifies that the transformation is well constructed
  /// according to class-specific rules
  virtual void Verify();

protected:

  /// Reads transformation parameters from a file stream
  virtual irtkCifstream &ReadDOFs(irtkCifstream &, irtkTransformationType);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkBSplineFreeFormTransformationStatistical::DOFGradientNorm(const double *gradient) const
{
  return irtkTransformation::DOFGradientNorm(gradient);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationStatistical::Put(int i, DOFValue x)
{
  irtkTransformation::Put(i, x);
  this->UpdateCPs();
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationStatistical::Put(const DOFValue *x)
{
  irtkTransformation::Put(x);
  this->UpdateCPs();
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformationStatistical::Add(const DOFValue *dx)
{
  irtkTransformation::Add(dx);
  this->UpdateCPs();
}

#endif

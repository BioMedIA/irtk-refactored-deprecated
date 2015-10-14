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

#ifndef _IRTKDEFORMABLESURFACEMODEL_H
#define _IRTKDEFORMABLESURFACEMODEL_H

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>

#include <irtkObjectiveFunction.h>
#include <irtkRegisteredPointSet.h>
#include <irtkRegisteredImage.h>
#include <irtkExternalForce.h>
#include <irtkPointSetConstraint.h>
#include <irtkTransformationConstraint.h>
#include <irtkEdgeTable.h>


/**
 * Energy function describing a deformable surface model
 */
class irtkDeformableSurfaceModel : public irtkObjectiveFunction
{
  irtkObjectMacro(irtkDeformableSurfaceModel);

private:

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input surface mesh
  irtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, Input);

  /// Intensity image
  irtkPublicAggregateMacro(irtkRegisteredImage, Image);

  /// Transformation to deform the surface with (optional)
  irtkPublicAggregateMacro(irtkTransformation, Transformation);

  /// Deformed surface mesh
  irtkAttributeMacro(irtkRegisteredPointSet, Surface);

  /// Number of energy terms
  irtkReadOnlyAttributeMacro(int, NumberOfTerms);

  /// Number of gradient smoothing iterations
  irtkPublicAttributeMacro(int, GradientSmoothing);

  /// Minimum (average) output mesh edge length
  irtkPublicAttributeMacro(double, MinEdgeLength);

  /// Maximum (average) output mesh edge length
  irtkPublicAttributeMacro(double, MaxEdgeLength);

  /// Minimum edge end point normal angle of feature edges
  irtkPublicAttributeMacro(double, MinFeatureAngle);

  /// Maximum edge end point normal angle of feature edges
  irtkPublicAttributeMacro(double, MaxFeatureAngle);

  /// Remesh deformed surface every n-th iteration
  irtkPublicAttributeMacro(int, RemeshInterval);

  /// Number of iterations since last performed remeshing
  irtkAttributeMacro(int, RemeshCounter);

  /// Energy terms corresponding to external forces
  vector<irtkExternalForce *> _ExternalForce;

  /// Energy terms corresponding to internal forces
  vector<irtkPointSetConstraint *> _InternalForce;

  /// Energy terms which regularize the surface transformation
  vector<irtkTransformationConstraint *> _Constraint;

  /// Cached values of individual energy terms
  vector<double> _Value;

public:

  /// Output surface mesh
  vtkSmartPointer<vtkPointSet> Output() const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy constructor
  /// \note Intentionally not implemented!
  irtkDeformableSurfaceModel(const irtkDeformableSurfaceModel &);

  /// Assignment operator
  /// \note Intentionally not implemented!
  irtkDeformableSurfaceModel &operator =(const irtkDeformableSurfaceModel &);

public:

  /// Constructor
  irtkDeformableSurfaceModel();

  /// Destructor
  virtual ~irtkDeformableSurfaceModel();

  // ---------------------------------------------------------------------------
  // Energy terms

  /// Initialize energy terms once input and parameters have been set
  virtual void Initialize();

  /// Delete previously added energy terms
  void Clear();

  /// Whether energy function has no terms
  bool Empty() const;

  /// Number of energy terms
  int NumberOfForces() const;

  /// Number of internal force terms
  int NumberOfInternalForces() const;

  /// Number of external force terms
  int NumberOfExternalForces() const;

  /// Add external force term and take over ownership of the object
  void Add(irtkExternalForce *);

  /// Remove external force term and revoke ownership of the object
  void Sub(irtkExternalForce *);

  /// Add internal force term and take over ownership of the object
  void Add(irtkPointSetConstraint *);

  /// Remove internal force term and revoke ownership of the object
  void Sub(irtkPointSetConstraint *);

  /// Add transformation regularization term and take over ownership of the object
  void Add(irtkTransformationConstraint *);

  /// Remove transformation regularization term and revoke ownership of the object
  void Sub(irtkTransformationConstraint *);

  /// Get the n-th energy term
  irtkEnergyTerm *Term(int);

  /// Get the n-th energy term
  const irtkEnergyTerm *Term(int) const;

  /// Get the n-th external force term
  irtkExternalForce *ExternalForce(int);

  /// Get the n-th external force term
  const irtkExternalForce *ExternalForce(int) const;

  /// Get the n-th internal force term
  irtkPointSetConstraint *InternalForce(int);

  /// Get the n-th internal force term
  const irtkPointSetConstraint *InternalForce(int) const;

  // ---------------------------------------------------------------------------
  // Settings

  // Import other overloads
  using irtkObjectiveFunction::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Function parameters

  /// Get number of deformable surface parameters
  virtual int NumberOfDOFs() const;

  /// Get number of deformable surface points
  int NumberOfPoints() const;

  /// Set deformable surface parameters
  ///
  /// This function can be used to restore the deformable surface parameters
  /// after a failed update which did not result in the desired improvement.
  ///
  /// \param[in] x Value of deformable surface parameters (DoFs).
  virtual void Put(const double *x);

  /// Get deformable surface parameters
  ///
  /// This function can be used to store a backup of the current deformable
  /// surface parameters before an update such that these can be restored using
  /// the Put member function if the update did not result in the desired change
  /// of the overall energy.
  ///
  /// \param[in] x Current values of deformable surface parameters (DoFs).
  virtual void Get(double *x) const;

  /// Get function parameter value
  ///
  /// \returns Value of specified function parameter (DoF).
  virtual double Get(int) const;

  /// Add change (i.e., scaled gradient) to each deformable surface parameter
  ///
  /// This function updates each parameter of the deformable surface model
  /// given a vector of desired changes, i.e., the computed gradient of the
  /// energy function.
  ///
  /// \param[in] dx Change of each function parameter (DoF) as computed by the
  ///               Gradient member function and scaled by a chosen step length.
  ///
  /// \returns Maximum change of transformation parameter.
  virtual double Step(const double *dx);

  /// Update internal state after change of parameters
  virtual void Update(bool = true);

  /// Update energy function after convergence
  ///
  /// For example, fiducial registration error (FRE) terms may update the
  /// point correspondences before another gradient-based optimization of
  /// the new FRE term.
  ///
  /// \returns Whether the energy function has changed.
  virtual bool Upgrade();

  /// Perform local adaptive remeshing
  virtual bool Remesh();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Query/evaluate initial value of energy function
  double InitialValue();

  /// Get initial value of n-th energy term
  double InitialValue(int);

  /// Evaluate energy function
  double Value();

  /// Get value of n-th energy term computed upon last evaluation
  double Value(int);

  /// Evaluate gradient of energy function
  ///
  /// This gradient corresponds to the weighted sum of external and internal
  /// forces of the deformable surface model.
  ///
  /// \param[in]  dx      Gradient of energy function.
  /// \param[in]  step    Step length for finite differences.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  void Gradient(double *dx, double step = .0, bool *sgn_chg = NULL);

  /// Compute norm of gradient of energy function
  ///
  /// This norm can, for example, be the maximum absolute parameter change,
  /// the maximum control point displacement if a FFD transformation is used
  /// to deform the initial surface mesh, or the maximum vertex displacement.
  double GradientNorm(const double *) const;

  /// Adjust step length range
  ///
  /// \param[in]    dx  Gradient of objective function.
  /// \param[inout] min Minimum step length.
  /// \param[inout] max Maximum step length.
  void GradientStep(const double *dx, double &min, double &max) const;

  /// Evaluate energy function
  ///
  /// This function first updates the internal state of the function object
  /// if required due to a previous change of the function parameters and then
  /// evaluates the current energy function value.
  ///
  /// \param[in]  step    Step length for finite differences.
  /// \param[out] dx      Gradient of objective function.
  ///                     If \c NULL, only the function value is computed.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  ///                     Ignord if \p dx is \c NULL.
  virtual double Evaluate(double *dx = NULL, double step = .0, bool *sgn_chg = NULL);

  /// Enforce hard constraints on surface model deformation
  ///
  /// This function clamps a nodes' displacement vector (velocity times \delta t),
  /// if otherwise the hard constraints of the deformable surface model would
  /// be violated. Common hard constraints are non-self-intersection and a
  /// maximum total node displacement. If the surface model is deformed by
  /// a parametric transformation, this function does nothing as hard constraints
  /// can only be enforced during the optimization when the parameters of the
  /// deformable surface model are the positions of the individual surface nodes.
  ///
  /// \param[inout] dx (Scaled) gradient of objective function.
  virtual void EnforceHardConstraints(double *dx) const;

  /// Laplacian smooth surface displacements such that neighboring points move coherently
  virtual void SmoothGradient(double *dx) const;

  // ---------------------------------------------------------------------------
  // Debugging

  /// Get unweighted and unnormalized value of n-th energy term
  /// \remarks Use for progress reporting only.
  double RawValue(int);

  /// Write input of data force terms
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of force terms
  virtual void WriteGradient(const char *, const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkPointSet> irtkDeformableSurfaceModel::Output() const
{
  return _Surface.PointSet();
}

// -----------------------------------------------------------------------------
inline int irtkDeformableSurfaceModel::NumberOfPoints() const
{
  return _Surface.NumberOfPoints();
}


#endif

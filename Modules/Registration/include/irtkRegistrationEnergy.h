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

#ifndef _IRTKREGISTRATIONENERGY_H
#define _IRTKREGISTRATIONENERGY_H

#include <irtkObjectiveFunction.h>
#include <irtkTransformation.h>

class irtkEnergyTerm;


/**
 * Registration energy term which sums up the individual composite terms
 */
class irtkRegistrationEnergy : public irtkObjectiveFunction
{
  irtkObjectMacro(irtkRegistrationEnergy);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of pre-update function delegate
  typedef fastdelegate::FastDelegate1<bool> PreUpdateFunctionType;

private:

  // ---------------------------------------------------------------------------
  // Attributes

  /// Transformation with free parameters of energy function
  irtkPublicAggregateMacro(irtkTransformation, Transformation);

  /// Individual terms of registration energy function
  vector<irtkEnergyTerm *> _Term;

  /// Cached values of each individual energy term
  ///
  /// The irtkRegistrationEnergy::Evaluate method which is called by
  /// irtkEnergyTerm::Value, updates this vector of registration energy values.
  /// It can be used, for example, to report the value of individual terms.
  vector<double> _Value;

  /// Delegate of pre-update function, e.g., to update input to energy terms
  irtkPublicAttributeMacro(PreUpdateFunctionType, PreUpdateFunction);

  /// Whether to normalize each individually energy term gradient
  ///
  /// If this option is enabled, also the weights which determine the influence
  /// of each energy term gradient are normalized by the sum of the weights.
  /// Note that the weights used to weigh the contribution of the energy term
  /// values themselves are not normalized, however, to allow more freedom in
  /// adjusting the relative magnitude of each term.
  ///
  /// \attention The use of this option is experimental!
  irtkPublicAttributeMacro(bool, NormalizeGradients);

  /// Energy gradient preconditioning sigma used to supress noise.
  /// A non-positive value disables the preconditioning all together.
  ///
  /// Zikic, D., Baust, M., Kamen, A., & Navab, N. A General Preconditioning
  /// Scheme for Difference Measures in Deformable Registration. In ICCV 2011.
  irtkPublicAttributeMacro(double, Preconditioning);

  /// Forward events of energy terms to observers of the energy function
  irtkEventDelegate _EventDelegate;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy constructor
  /// \note Intentionally not implemented!
  irtkRegistrationEnergy(const irtkRegistrationEnergy &);

  /// Assignment operator
  /// \note Intentionally not implemented!
  irtkRegistrationEnergy &operator =(const irtkRegistrationEnergy &);

public:

  /// Constructor
  irtkRegistrationEnergy();

  /// Destructor
  virtual ~irtkRegistrationEnergy();

  // ---------------------------------------------------------------------------
  // Energy terms

  /// Initialize energy terms once input and parameters have been set
  virtual void Initialize();

  /// Delete previously added energy terms
  void Clear();

  /// Whether energy function has no terms
  bool Empty() const;

  /// Number of energy terms
  int NumberOfTerms() const;

  /// Add energy term and take over ownership of the object
  void Add(irtkEnergyTerm *);

  /// Remove energy term and revoke ownership of the object
  void Sub(irtkEnergyTerm *);

  /// Get the n-th energy term
  irtkEnergyTerm *Term(int);

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

  /// Get number of transformation parameters
  ///
  /// \returns Number of transformation parameters (DoFs).
  virtual int NumberOfDOFs() const;
  
  /// Set transformation parameters
  ///
  /// This function can be used to restore the transformation parameters after
  /// a failed update which did not result in the desired improvement.
  ///
  /// \param[in] x Value of transformation parameters (DoFs).
  virtual void Put(const double *x);

  /// Get transformation parameters
  ///
  /// This function can be used to store a backup of the current transformation
  /// parameters before an update such that these can be restored using the Put
  /// member function if the update did not result in the desired change of the
  /// overall registration energy.
  ///
  /// \param[in] x Current values of transformation parameters (DoFs).
  virtual void Get(double *x) const;

  /// Get function parameter value
  ///
  /// \returns Value of specified function parameter (DoF).
  virtual double Get(int) const;

  /// Add change (i.e., scaled gradient) to each transformation parameter
  ///
  /// This function updates each parameter of the registration energy function
  /// given a vector of desired changes, i.e., the computed gradient of the
  /// energy function w.r.t. these parameters or a desired change computed
  /// otherwise such as through matching patches (in a very loose sense an
  /// approximate "gradient") or a discrete assignment (c.f. MRF registration).
  /// It moreover triggers an update of the internal state of each energy term
  /// which depends on the current transformation. A data similarity term,
  /// for example, will update its moving (transformed) input data.
  ///
  /// \note How this change is applied to the transformation parameters depends
  ///       on the particular underlying transformation model. Often \p dx is
  ///       simply added to \c x, the current transformation parameters. Some
  ///       transformations may compute a slightly different update from \p dx,
  ///       however, such as converting displacements into velocities first.
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

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Query/evaluate initial value of registration energy
  double InitialValue();

  /// Get initial value of n-th energy term
  double InitialValue(int);

  /// Evaluate registration energy
  double Value();

  /// Get value of n-th energy term computed upon last evaluation
  double Value(int);

  /// Normalize energy gradient
  ///
  /// Zikic, D., Baust, M., Kamen, A., & Navab, N. A General Preconditioning
  /// Scheme for Difference Measures in Deformable Registration. In ICCV 2011.
  void NormalizeGradient(double *dx);

  /// Evaluate gradient of registration energy
  ///
  /// Note that this gradient need not necessarily correspond to the analytic
  /// gradient of the energy function. It can also be a desired change of each
  /// function parameter computed otherwise such as through a local patch
  /// match. The only requirement is that the change will increase the energy
  /// function value such that successively walking in the direction of the
  /// computed "gradient", a (local) maximum of the registration energy will be
  /// be reached.
  ///
  /// \param[in]  dx      Gradient of registration energy function.
  /// \param[in]  step    Step length for finite differences.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  void Gradient(double *dx, double step = .0, bool *sgn_chg = NULL);

  /// Compute norm of gradient of registration energy
  ///
  /// This norm can, for example, be the maximum absolute parameter change or
  /// the maximum control point displacement in case of a FFD transformation model.
  double GradientNorm(const double *) const;

  /// Adjust step length range
  ///
  /// \param[in]    dx  Gradient of objective function.
  /// \param[inout] min Minimum step length.
  /// \param[inout] max Maximum step length.
  void GradientStep(const double *dx, double &min, double &max) const;

  /// Evaluate registration energy
  ///
  /// This function first updates the internal state of the function object
  /// if required due to a previous change of the function parameters and then
  /// evaluates the current registration energy.
  ///
  /// \param[in]  step    Step length for finite differences.
  /// \param[out] dx      Gradient of objective function.
  ///                     If \c NULL, only the function value is computed.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  ///                     Ignord if \p dx is \c NULL.
  virtual double Evaluate(double *dx = NULL, double step = .0, bool *sgn_chg = NULL);

  // ---------------------------------------------------------------------------
  // Debugging

  /// Get unweighted and unnormalized value of n-th energy term
  /// \remarks Use for progress reporting only.
  double RawValue(int);

  /// Write input of data fidelity terms
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity terms w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

};


#include <irtkEnergyTerm.h>


#endif

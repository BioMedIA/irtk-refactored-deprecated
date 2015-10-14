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

#ifndef _IRTKOBJECTIVEFUNCTION_H
#define _IRTKOBJECTIVEFUNCTION_H

#include <irtkCommon.h>
#include <irtkObject.h>


/**
 * Interface for objective function of an optimization problem.
 *
 * An objective function can be minimized/maximized using a specialization
 * of the irtkLocalOptimizer base class.
 */
class irtkObjectiveFunction : public irtkObservable
{
  irtkAbstractMacro(irtkObjectiveFunction);

protected:

  /// Default step length used for gradient approximation
  irtkAttributeMacro(double, StepLength);

public:

  /// Destructor
  virtual ~irtkObjectiveFunction() = 0;

  // ---------------------------------------------------------------------------
  // Function parameters (DoFs)

  /// Get number of DoFs
  ///
  /// \returns Number of free function parameters.
  virtual int NumberOfDOFs() const = 0;

  /// Set function parameter values
  ///
  /// This is function can be used to set the parameters of the objective function
  /// to particular values. In particular, it can be used to restore the function
  /// parameters after a failed incremental update which did not result in the
  /// desired improvement.
  ///
  /// \param[in] x Function parameter (DoF) values.
  virtual void Put(const double *x) = 0;

  /// Get function parameter value
  ///
  /// \param[in] i Function parameter (DoF) index.
  ///
  /// \returns Value of specified function parameter (DoF).
  virtual double Get(int i) const = 0;

  /// Get function parameter values
  ///
  /// This function can be used to store a backup of the current funtion parameter
  /// values before an update such that these can be restored using the Put
  /// member function if the update did not result in the desired change of the
  /// overall objective function value.
  ///
  /// \param[in] x Function parameter (DoF) values.
  virtual void Get(double *x) const = 0;

  /// Add change (i.e., scaled gradient) to each parameter value
  ///
  /// This function updates each DoF of the objective function given a vector
  /// of corresponding changes, i.e., the computed gradient of the objective
  /// function w.r.t. these parameters or a desired change computed otherwise.
  ///
  /// \param[in] dx Change of each function parameter (DoF) as computed by the
  ///               Gradient member function and scaled by a chosen step length.
  ///
  /// \returns Maximum change of function parameter.
  virtual double Step(const double *dx) = 0;

  /// Update internal state after change of parameters
  ///
  /// \param[in] gradient Update also internal state required for evaluation of
  ///                     gradient of objective function.
  virtual void Update(bool gradient = true);

  /// Update objective function after convergence
  ///
  /// This function may be called by the optimizer after the optimization of
  /// the current objective function has converged or the improvement is
  /// insignificant. It allows the objective function to modify itself
  /// considering the current estimate of the parameters (DoFs). For example,
  /// the fiducial registration error term of a registration energy function
  /// could update the point correspondences.
  ///
  /// \returns Whether the objective function has changed.
  virtual bool Upgrade();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate objective function value
  virtual double Value() = 0;

  /// Evaluate gradient of objective function w.r.t its DoFs
  ///
  /// \param[int] step    Step length for finite differences.
  /// \param[out] dx      Gradient of objective function.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  virtual void Gradient(double *dx, double step = .0, bool *sgn_chg = NULL) = 0;

  /// Compute norm of gradient of objective function
  ///
  /// This norm is used to define a unit for the step length used by gradient
  /// descent methods. It is, for example, the maximum absolute value norm for
  /// linear transformations and the maximum control point displacement for FFDs.
  /// The computation of the norm may be done after conjugating the gradient
  /// vector obtained using the Gradient member function.
  ///
  /// \param[in] dx Gradient of objective function.
  virtual double GradientNorm(const double *dx) const = 0;

  /// Adjust step length range
  ///
  /// \param[in]    dx  Gradient of objective function.
  /// \param[inout] min Minimum step length.
  /// \param[inout] max Maximum step length.
  virtual void GradientStep(const double *dx, double &min, double &max) const;

  /// Evaluate objective function
  ///
  /// This function first updates the internal state of the function object
  /// if required due to a previous change of the function parameters and then
  /// evaluates the current function value. If \p dx is not \c NULL, the
  /// gradient of the objective function is computed as well.
  ///
  /// \param[int] step    Step length for finite differences.
  /// \param[out] dx      Gradient of objective function.
  ///                     If \c NULL, only the function value is computed.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  ///                     Ignord if \p dx is \c NULL.
  ///
  /// \returns Value of the objective function.
  virtual double Evaluate(double *dx = NULL, double step = .0, bool *sgn_chg = NULL);

  /// Enforce hard constraints
  ///
  /// This function is in particular called by optimizers of a deformable surface
  /// model to enforce that the surface does not self-intersect after the
  /// gradient step by setting the displacement of points to zero if necessary.
  ///
  /// \param[inout] dx (Scaled) gradient of objective function.
  virtual void EnforceHardConstraints(double *dx) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline irtkObjectiveFunction::~irtkObjectiveFunction()
{
}

// -----------------------------------------------------------------------------
inline void irtkObjectiveFunction::Update(bool)
{
}

// -----------------------------------------------------------------------------
inline bool irtkObjectiveFunction::Upgrade()
{
  return false;
}

// -----------------------------------------------------------------------------
inline void irtkObjectiveFunction::GradientStep(const double *, double &, double &) const
{
  // By default, step length range chosen by user/optimizer
}

// -----------------------------------------------------------------------------
inline double irtkObjectiveFunction::Evaluate(double *dx, double step, bool *sgn_chg)
{
  this->Update(dx != NULL);
  if (dx) {
    if (step <= .0) step = _StepLength;
    this->Gradient(dx, step, sgn_chg);
  }
  return this->Value();
}

// -----------------------------------------------------------------------------
inline void irtkObjectiveFunction::EnforceHardConstraints(double *) const
{
}


#endif

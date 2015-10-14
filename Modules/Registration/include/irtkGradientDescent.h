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

#ifndef _IRTKGRADIENTDESCENT_H
#define _IRTKGRADIENTDESCENT_H

#include <irtkLocalOptimizer.h>
#include <irtkLineSearch.h>
#include <irtkEventDelegate.h>


/**
 * Minimizes objective function using gradient descent
 */
class irtkGradientDescent : public irtkLocalOptimizer
{
  irtkOptimizerMacro(irtkGradientDescent, GradientDescent);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of gradient descent steps
  irtkPublicAttributeMacro(int, NumberOfSteps);

  /// Maximum number of restarts after upgrade of energy function
  irtkPublicAttributeMacro(int, NumberOfRestarts);

  /// Maximum streak of unsuccessful restarts without improvement
  irtkPublicAttributeMacro(int, NumberOfFailedRestarts);

  /// Line search strategy
  irtkPublicAttributeMacro(irtkLineSearchStrategy, LineSearchStrategy);

  /// (Possible) Line search parameter
  irtkAttributeMacro(irtkParameterList, LineSearchParameter);

  /// Line search optimization method
  irtkReadOnlyAggregateMacro(irtkLineSearch, LineSearch);

  /// Whether line search object is owned by this optimizer
  irtkAttributeMacro(bool, LineSearchOwner);

  /// Forwards line search event messages to observers of optimization
  irtkEventDelegate _EventDelegate;

protected:

  /// Allocated memory for line search direction
  irtkComponentMacro(double, Gradient);

  /// Whether to allow function parameter sign to change
  ///
  /// If \c false for a particular function parameter, the line search sets
  /// the function parameter to zero whenever the sign of the parameter would
  /// change when taking a full step along the scaled gradient direction.
  irtkComponentMacro(bool, AllowSignChange);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkGradientDescent(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkGradientDescent(const irtkGradientDescent &);

  /// Assignment operator
  irtkGradientDescent &operator =(const irtkGradientDescent &);

  /// Destructor
  virtual ~irtkGradientDescent();

  // Import overloads from base class
  using irtkLocalOptimizer::Function;

  /// Set objective function
  virtual void Function(irtkObjectiveFunction *);

  /// Set line search object
  virtual void LineSearch(irtkLineSearch *, bool = false);

  // ---------------------------------------------------------------------------
  // Parameters
  using irtkLocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize gradient descent
  ///
  /// This member funtion is implicitly called by Run. It can, however, be
  /// called prior to Run explicitly in order to be able to set up the line
  /// search instance. Otherwise, use the generic Set member function to
  /// change the line search parameters and simply have Run call Initialize.
  virtual void Initialize();

  /// Optimize objective function using gradient descent
  virtual double Run();

protected:

  /// Finalize gradient descent
  virtual void Finalize();

  /// Compute descent direction and corresponding step length unit
  virtual double Gradient(double *, double = .0, bool * = NULL);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline double irtkGradientDescent::Gradient(double *gradient, double step, bool *sgn_chg)
{
         Function()->Gradient    (gradient, step, sgn_chg);
  return Function()->GradientNorm(gradient);
}


#endif

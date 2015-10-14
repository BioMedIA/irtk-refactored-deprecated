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

#ifndef _IRTKINEXACTLINESEARCH_H
#define _IRTKINEXACTLINESEARCH_H

#include <irtkLineSearch.h>


/**
 * Searches sufficiently optimal step length along search direction
 *
 * This local optimizer implements an inexact line search with adaptive step
 * size control, increasing the step size while steps are accepted, and
 * decreasing it when a step did not yield a sufficient improvement.
 */
class irtkInexactLineSearch : public irtkLineSearch
{
  irtkAbstractMacro(irtkInexactLineSearch);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number consecutive rejected steps
  irtkPublicAttributeMacro(int, MaxRejectedStreak);

  /// Whether to start new search using step length of previous search
  irtkPublicAttributeMacro(bool, ReusePreviousStepLength);

  /// Whether (reused) step length is strictly limited by [min, max]
  ///
  /// 0: Step length is allowed to exceed max which may be the case when the
  ///    previously accumulated total step length is being reused.
  /// 1: Incremental steps are strictly limited to [min, max].
  /// 2: Accumulated total step length is strictly limited to [min, max].
  irtkPublicAttributeMacro(int, StrictStepLengthRange);

  /// Whether to refuse any function parameter sign changes
  ///
  /// If \c false for a particular function parameter, the line search sets
  /// the function parameter to zero whenever the sign of the parameter would
  /// change when taking a full step along the scaled gradient direction.
  irtkPublicAggregateMacro(bool, AllowSignChange);

  /// Previous function parameters
  irtkComponentMacro(double, CurrentDoFValues);

  /// Line search direction scaled by current step length
  irtkComponentMacro(double, ScaledDirection);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:
  using irtkLineSearch::Function;

  /// Constructor
  irtkInexactLineSearch(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkInexactLineSearch(const irtkInexactLineSearch &);

  /// Assignment operator
  irtkInexactLineSearch &operator =(const irtkInexactLineSearch &);

  /// Destructor
  virtual ~irtkInexactLineSearch();

  /// Set objective function
  virtual void Function(irtkObjectiveFunction *);

  // ---------------------------------------------------------------------------
  // Parameters
  using irtkLineSearch::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization
protected:

  /// Take step in search direction
  ///
  /// \returns Maximum change of DoF
  double Advance(double);

  /// Revert step in search direction
  void Retreat(double);

};


#endif

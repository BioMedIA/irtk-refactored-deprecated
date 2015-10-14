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

#ifndef _IRTKMINJACOBIANCONSTRAINT_H
#define _IRTKMINJACOBIANCONSTRAINT_H

#include <irtkJacobianConstraint.h>


/**
 * Penalizes non-diffeomorphic transformations by enforcing a minimum Jacobian determinant
 *
 * Rueckert et al., Diffeomorphic Registration Using B-Splines, MICCAI 2006.
 */
class irtkMinJacobianConstraint : public irtkJacobianConstraint
{
  irtkObjectMacro(irtkMinJacobianConstraint);

  /// Lower Jacobian determinant threshold, only determinant values less or
  /// equal this threshold are penalized
  irtkPublicAttributeMacro(double, Gamma);

public:

  /// Constructor
  irtkMinJacobianConstraint(const char * = "");

  /// Destructor
  virtual ~irtkMinJacobianConstraint();

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using irtkJacobianConstraint::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


#endif

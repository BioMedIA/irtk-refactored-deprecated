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

#ifndef _IRTKSPARSITYCONSTRAINT_H
#define _IRTKSPARSITYCONSTRAINT_H

#include <irtkTransformationConstraint.h>


/**
 * Sparsity contraint based on L1 norm of transformation parameters
 *
 * Wenzhe Shi et al., Registration using sparse free-form deformations,
 * MICCAI 2012 (http://www.ncbi.nlm.nih.gov/pubmed/23286105)
 */
class irtkSparsityConstraint : public irtkTransformationConstraint
{
  irtkObjectMacro(irtkSparsityConstraint);

public:

  /// Constructor
  irtkSparsityConstraint(const char * = "");

  /// Evaluate gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the evaluated gradient of this
  ///                         energy term is added to with its resp. weight.
  /// \param[in]     step     Step length for finite differences (unused).
  /// \param[out]    sgn_chg  Whether sign change of parameter should be allowed.
  void Gradient(double *gradient, double step, bool *sgn_chg = NULL);

  /// Adjust step length range
  ///
  /// \param[in]     gradient Gradient of objective function.
  /// \param[in,out] min      Minimum step length.
  /// \param[in,out] max      Maximum step length.
  virtual void GradientStep(const double *gradient, double &min, double &max) const;

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double, bool *);

};


#endif

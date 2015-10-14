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

#include <irtkSparsityConstraint.h>


// -----------------------------------------------------------------------------
irtkSparsityConstraint::irtkSparsityConstraint(const char *name)
:
  irtkTransformationConstraint(name)
{
  _ConstrainPassiveDoFs = false;
}

// -----------------------------------------------------------------------------
void irtkSparsityConstraint::Gradient(double *gradient, double step, bool *sgn_chg)
{
  this->EvaluateGradient(gradient, step, this->_Weight, sgn_chg);
}

// -----------------------------------------------------------------------------
void irtkSparsityConstraint::GradientStep(const double *, double &min, double &max) const
{
  const int ndofs = Transformation()->NumberOfDOFs();

  double norm    = .0;
  int    nactive = 0;

  for (int dof = 0; dof < ndofs; ++dof) {
    if (_ConstrainPassiveDoFs || (Transformation()->GetStatus(dof) == Active)) {
      norm += fabs(Transformation()->Get(dof));
      ++nactive;
    }
  }

  if (norm == .0 || nactive == 0) return;
  norm /= nactive;

  if (min < norm) {
    min = norm / 128;
    max = norm;
  }
}

// -----------------------------------------------------------------------------
double irtkSparsityConstraint::Evaluate()
{
  const int ndofs = Transformation()->NumberOfDOFs();

  double penalty = .0;
  int    nactive = 0;

  for (int dof = 0; dof < ndofs; ++dof) {
    if (_ConstrainPassiveDoFs || (Transformation()->GetStatus(dof) == Active)) {
      penalty += fabs(Transformation()->Get(dof));
      ++nactive;
    }
  }

  return (nactive > 0 ? (penalty / nactive) : .0);
}

// -----------------------------------------------------------------------------
void irtkSparsityConstraint::EvaluateGradient(double *gradient, double, double weight, bool *sgn_chg)
{
  const int ndofs   = Transformation()->NumberOfDOFs();
  const int nactive = Transformation()->NumberOfActiveDOFs();
  if (nactive == 0) return;

  const double norm = weight / nactive;

  double value;
  for (int dof = 0; dof < ndofs; ++dof) {
    if (_ConstrainPassiveDoFs || (Transformation()->GetStatus(dof) == Active)) {
      value = Transformation()->Get(dof);
      if      (value < .0) value = - norm; // = w * x / |x|
      else if (value > .0) value = + norm;
      gradient[dof] -= value;
      // Disallow sign change of parameter value if caused by sparsity gradient
      if (sgn_chg && (gradient[dof] * value) < 0) sgn_chg[dof] = false;
    }
  }
}

// -----------------------------------------------------------------------------
void irtkSparsityConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
  this->EvaluateGradient(gradient, step, weight, NULL);
}

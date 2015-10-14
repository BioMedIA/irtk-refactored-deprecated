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

#include <irtkEulerMethodWithMomentum.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkEulerMethodWithMomentum::irtkEulerMethodWithMomentum(irtkObjectiveFunction *f)
:
  irtkEulerMethod(f),
  _BodyMass(1e-3)
{
  _DampingFactor = 1e-3;
}

// -----------------------------------------------------------------------------
irtkEulerMethodWithMomentum::irtkEulerMethodWithMomentum(const irtkEulerMethodWithMomentum &other)
:
  irtkEulerMethod(other),
  _BodyMass(other._BodyMass)
{
}

// -----------------------------------------------------------------------------
irtkEulerMethodWithMomentum &irtkEulerMethodWithMomentum::operator =(const irtkEulerMethodWithMomentum &other)
{
  irtkEulerMethod::operator =(other);
  _BodyMass = other._BodyMass;
  return *this;
}

// -----------------------------------------------------------------------------
irtkEulerMethodWithMomentum::~irtkEulerMethodWithMomentum()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkEulerMethodWithMomentum::Set(const char *name, const char *value)
{
  if (strcmp(name, "Deformable surface mass") == 0) {
    return FromString(value, _BodyMass);
  }
  return irtkEulerMethod::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkEulerMethodWithMomentum::Parameter() const
{
  irtkParameterList params = irtkEulerMethod::Parameter();
  Insert(params, "Deformable surface mass", _BodyMass);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEulerMethodWithMomentum::ComputeVelocity(double *v, const double *f) const
{
  const double s = _StepLength / _BodyMass;
  v[0] += s * (f[0] - _DampingFactor * v[0]);
  v[1] += s * (f[1] - _DampingFactor * v[1]);
  v[2] += s * (f[2] - _DampingFactor * v[2]);
}

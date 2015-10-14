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


#include <irtkSurfaceConstraint.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkSurfaceConstraint::irtkSurfaceConstraint(const char *name, double weight)
:
  irtkPointSetConstraint(name, weight)
{
  _SurfaceForce = true;
}

// -----------------------------------------------------------------------------
irtkSurfaceConstraint::irtkSurfaceConstraint(const irtkSurfaceConstraint &other)
:
  irtkPointSetConstraint(other)
{
}

// -----------------------------------------------------------------------------
irtkSurfaceConstraint &irtkSurfaceConstraint::operator =(const irtkSurfaceConstraint &other)
{
  irtkPointSetConstraint::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkSurfaceConstraint::~irtkSurfaceConstraint()
{
}

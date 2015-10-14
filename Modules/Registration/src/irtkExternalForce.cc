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

#include <irtkExternalForce.h>


// =============================================================================
// Factory
// =============================================================================

// -----------------------------------------------------------------------------
irtkExternalForce *irtkExternalForce::New(irtkPointSetForceMeasure force)
{
  switch (force) {
    default:
      cerr << "irtkExternalForce::New: Unknown external point set force: " << force << endl;
      exit(1);
  }
  return NULL;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkExternalForce::Copy(const irtkExternalForce &other)
{
  _Image = other._Image;
}

// -----------------------------------------------------------------------------
irtkExternalForce::irtkExternalForce(const char *name, double weight)
:
  irtkPointSetForce(name, weight),
  _Image(NULL)
{
}

// -----------------------------------------------------------------------------
irtkExternalForce::irtkExternalForce(const irtkExternalForce &other)
:
  irtkPointSetForce(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkExternalForce &irtkExternalForce::operator =(const irtkExternalForce &other)
{
  irtkPointSetForce::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkExternalForce::~irtkExternalForce()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkExternalForce::Initialize()
{
  // Initialize base class
  irtkPointSetForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Check input
  if (_Image == NULL) {
    cerr << "irtkExternalForce::Initialize: Reference image is NULL" << endl;
    exit(1);
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkExternalForce::Update(bool gradient)
{
  // Before irtkPointSetForce::Update which sets _InitialUpdate = false
  if (_InitialUpdate || _Image->Transformation()) {
    _Image->Update(_InitialUpdate && _Image->SelfUpdate());
  }
  irtkPointSetForce::Update(gradient);
}

// -----------------------------------------------------------------------------
double irtkExternalForce::Evaluate()
{
  // Not all external forces may have a corresponding energy term expression
  // if the deformable surface model is rather based on reaching an equilibrium
  // of internal and external forces. In this case, return infinity such that
  // the optimizer only checks the equilibrium condition
  return numeric_limits<double>::infinity();
}
